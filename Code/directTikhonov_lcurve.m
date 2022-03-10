function [SusceptibilityMap, Alpha, Corr] = directTikhonov_lcurve(fieldMap, mask, parameters)
% DESCRIPTION: SusceptibilityMap = directTikhonov(FieldMap, Mask, Parameters)
% Calculates the susceptibility map by solving a least squares
% problem with closed-form Tikhonov regularisation:
%   x = argmin_x || B - B0(x*d) ||_2 + alpha||x||_2 where * denotes
%   convolution
%
% INPUTS:
% fieldMap: field map [ppm]
% mask: binary tissue mask
% parameters: Parameters.Alpha - Enter a single value of alpha or a range 
%               of alpha for l-curve optimisation. Default scans alpha from 10E-5 to 1
%             Parameters.Resolution(double vector) - Image resolution vector [dx dy dz] in mm
%             Parameters.Orientation(double vector) - 3-element vector of the B0 direction (default: [0 0 1])
%             Parameters.PSFCorr(string) - PSF corresciton: 'Yes' (default) or 'No'
%
% OUTPUT: 
% SusceptibilityMap(double matrix) - Susceptibility map in ppm
%
% AUTHORS: 
% Emma Biondetti, Emma Dixon, Anita Karsa, (2017) University College London, UK 

%% Initialisation of parameters
%Resolution
Resolution = parameters.Resolution;
%Orientation
if isfield(parameters,'Orientation')
    z_prime = parameters.Orientation;
    z_prime = z_prime/norm(z_prime);
else
    z_prime = [0 0 1];
end
%Matrix Size
MatrixSize = size(fieldMap);
%Regularisation parameter
if isfield(parameters,'Alpha')
    if length(parameters.Alpha)==1   
        Alpha = parameters.Alpha;
    else
        Alpha_range = parameters.Alpha;
        Alpha = l_curve(fieldMap, mask, Resolution, MatrixSize, z_prime, Alpha_range);
    end
else
    Alpha_range = logspace(-5,0,100);
    Alpha = l_curve(fieldMap, mask, Resolution, MatrixSize, z_prime, Alpha_range);
end

%% Fourier transform image

Fk = fftshift(fftn(ifftshift(fieldMap)));

%% Create tikhonov k-space kernel

Dk = dipole(Resolution, MatrixSize,z_prime);
Kernel = Dk./(Dk.^2 + Alpha);

%% Calculate PSF correction factor

Corr = ifftn(ifftshift(Dk.*Kernel));
Corr = Corr(1);

%% 'Deconvolution'

switch isfield(parameters,'PSFCorr')
    case 0
        Chik = Fk.*Kernel/Corr;
    case 1
        switch parameters.PSFCorr
            case 'Yes'
                Chik = Fk.*Kernel/Corr;
            case 'No'
                Chik = Fk.*Kernel;
        end
end

%% Transform back to image space
SusceptibilityMap = fftshift(ifftn(ifftshift(Chik)));
SusceptibilityMap = real(SusceptibilityMap).*mask;

end

%% Functions
function Alpha_opt = l_curve(FieldMap, Mask,Resolution, MatrixSize,z_prime,Alpha)
%Calculate l_curve and find point of maximum curvature 

%L-curve: Hansen PC. The L-curve and its use in the numerical treatment of
%inverse problems. IMM, Department of Mathematical Modelling, Technical
%Universityof Denmark; 1999.


    Dk = dipole(Resolution, MatrixSize,z_prime);
    
    %pre-allocate
    Residual = zeros(length(Alpha),1);
    Solution = zeros(length(Alpha),1);
    
    for i = 1:length(Alpha)
        Kernel = Dk./(Dk.^2 + Alpha(i)); 
        
        %calculate residuals and solution norms
        b = FieldMap;
        x_alpha_ft = Kernel.*fftshift(fftn(ifftshift(FieldMap))); 
        x_alpha = fftshift(ifftn(ifftshift(x_alpha_ft)));
        Ax = fftshift(ifftn(ifftshift(Dk.*x_alpha_ft)));
        
        Residual(i) = norm(Ax(Mask==1) - b(Mask==1));
        Solution(i) = norm(x_alpha(Mask==1));
    end
    

    %find maximum curvature of l_curve
    M = [0 3 0 0; 0 0 2 0; 0 0 0 1; 0 0 0 0];

    X_alpha = log(Residual.^2);
    Y_alpha = log(Solution.^2);

    % Fitting a spline to X_alpha
    pp_X_alpha = spline(Alpha, X_alpha);
    pp_X_alpha.coefs = pp_X_alpha.coefs * M;
    X_alpha_del = ppval(pp_X_alpha, Alpha);
    pp_X_alpha.coefs = pp_X_alpha.coefs * M;
    X_alpha_del2 = ppval(pp_X_alpha, Alpha);

    % Fitting a spline to Y_alpha
    pp_Y_alpha = spline(Alpha, Y_alpha);
    pp_Y_alpha.coefs = pp_Y_alpha.coefs * M;
    Y_alpha_del = ppval(pp_Y_alpha, Alpha);
    pp_Y_alpha.coefs = pp_Y_alpha.coefs * M;
    Y_alpha_del2 = ppval(pp_Y_alpha, Alpha);

    % Calculate curvature
    Kappa = 2 * (X_alpha_del2 .* Y_alpha_del - X_alpha_del .* Y_alpha_del2) ./ ...
        ((X_alpha_del.^2 + Y_alpha_del.^2).^1.5);

    idx = Kappa == max(Kappa);
    Alpha_opt = Alpha(idx);    
    
    % Comment/uncomment to produce figure of L-curve
    figure
    subplot(121)
    loglog(Residual, Solution, '-db')
    hold on
    loglog(Residual(idx), Solution(idx), '-or');
    hold off
    xlabel('log(Consistency)')
    ylabel('log(Regularisation)')    
    subplot(122)
    semilogx(Alpha, Kappa, '-xb')
    hold on
    plot(Alpha(idx), Kappa(idx), '-or')
    hold off
    xlabel('log(Regularisation parameter)')
    ylabel('Curvature')
    
end

function Dk = dipole(Resolution, MatrixSize,z_prime)
    dkx = 1/Resolution(1)/MatrixSize(1);
    dky = 1/Resolution(2)/MatrixSize(2);
    dkz = 1/Resolution(3)/MatrixSize(3);

    Center=ceil(MatrixSize/2);

    Range = cell(3,1);

    for Index = 1:3
        if mod(MatrixSize(Index),2)==0
            Range{Index} = -Center(Index):Center(Index)-1;
        else
            Range{Index} = -Center(Index)+1:Center(Index)-1;
        end
    end

    [Y,X,Z] = meshgrid(Range{2}*dky,...
                       Range{1}*dkx,...
                       Range{3}*dkz);

    K2 = X.^2+Y.^2+Z.^2;

    Dk = 1/3 - (z_prime(1)*X + z_prime(2)*Y + z_prime(3)*Z).^2./K2;
    Dk(isnan(Dk)) = 0;

end
