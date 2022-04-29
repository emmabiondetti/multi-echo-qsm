% AUTHOR: Emma Biondetti, University of Chieti-Pescara, Italy
% DATE: 10/03/2022
% DESCRIPTION: this code generates a multi-echo complex data set based on
% the 1-mm isotropic numerical head phantom with variable M0, T2* and
% susceptibility values in each region

clear; close all; clc

%% Loading phantom data

nii = load_nii('zubalPhantom_M0_wholeHead_ellipsoid_1iso.nii');
M0_1iso = double(nii.img);

nii = load_nii('zubalPhantom_susceptibility_wholeHead_ellipsoid_1iso.nii');
susc_1iso = double(nii.img);

nii = load_nii('zubalPhantom_T2_star_wholeHead_ellipsoid_1iso.nii');
T2Star_1iso = double(nii.img);

nii = load_nii('zubalPhantom_brainMask_1iso.nii');
brain_mask_1iso = double(nii.img);

% Loading 4th order polynomial fit to the in vivo data
% Available in Data\ZubalPhantom
load('par_phi0_4th_order_poly.mat')

% Imaging parameters
matrix_size = size(M0_1iso); % [voxels]
orig = matrix_size/2+1;
voxel_size = [1 1 1]; % [mm]
TE_first = 3 * 1e-3; % [s]
Delta_TE = 5.4 * 1e-3; % [s]
n_echoes = 5;
TE = TE_first : Delta_TE : Delta_TE*n_echoes; % [s]

% Calculating head mask
head_mask = ones(matrix_size);
head_mask(M0_1iso == 0) = 0;

%% Total field perturbation calculated with a forward model [T]
% Convolution of susceptibility with the unit dipole kernel

[D, ~] = unitDipole(matrix_size, voxel_size, 0, [0 0 1]');
D_shifted = ifftshift(D);

% 1e-6: correction factor for the susceptibility values, which are in ppm.
% B0: correction factor for field strength, as the convolution d*x gives
% the field map normalised to the field strength [adimensional]
B0 = 3; % [T]
totalField_T = ifftn(fftn(susc_1iso) .* D_shifted) * 1e-6 * B0; % [T]

gamma = 42.6 * 2*pi * 1e6; % hydrogen proton gyromagnetic ratio [rad/s/T]

totalField_Hz = totalField_T * gamma /2/pi; % [Hz]
totalField_rad_s = totalField_T * gamma; % [rad/s]

%% Reference scan method

susc_refScan = zeros(matrix_size);
susc_refScan(susc_1iso == 9.4) = 9.4;

bgField_T = ...
    ifftn(fftn(susc_refScan) .* D_shifted) * 1e-6 * B0; % [T]

bgField_Hz = bgField_T * gamma /2/pi; % [Hz]
bgField_rad_s = bgField_T * gamma; % [rad/s]

% True local field
localField_Hz = totalField_Hz - bgField_Hz;

%% Simulating phase offset
% Only in the voxels in the head mask
hp = find(head_mask == 1);

[X,Y,Z] = ndgrid(-matrix_size(1)/2:matrix_size(1)/2-1,...
    -matrix_size(2)/2:matrix_size(2)/2-1,...
    -matrix_size(3)/2:matrix_size(3)/2-1);

X = X * (voxel_size(1));
Y = Y * (voxel_size(2));
Z = Z * (voxel_size(3));

% Building the model
[modelterms, varlist] = buildcompletemodel(4, 3);
A4 = zeros(length(hp), size(modelterms, 1));
for i = 1:size(modelterms, 1)
    i_exp = modelterms(i,:);
    A4(:,i) = (X(hp).^i_exp(1)).*(Y(hp).^i_exp(2)).*(Z(hp).^i_exp(3));
end
yhat = A4*parameters_phase_4th_order_fit(:);

y_est_tot = zeros(prod(matrix_size),1);
y_est_tot(hp) = yhat;
phase_offset_est = reshape(y_est_tot, matrix_size);

%% Saving phase offset

nii.img = phase_offset_est;
save_nii(nii, 'poly_phi0_1iso.nii.gz')

%% Complex signal simulation
% With and without realistic noise levels and a phase offset based on a 2nd
% order polynomial fit of true MRI data

% Parameters for noise in the real and imaginary channels
noise_mean = 0;
noise_sd = 0.07;

% complex_sim = zeros([matrixSize nEchoes]);
phase_phi0_sd = zeros([matrix_size n_echoes]);
magnitude_sd = zeros([matrix_size n_echoes]);

msk = zeros([matrix_size n_echoes]);
for echo = 1:n_echoes
    msk(:,:,:,echo) = head_mask; 
end %for 

for echo = 1:n_echoes
    magnitude_sd(:,:,:, echo) = ...
        M0_1iso .* exp(- TE(echo) / T2Star_1iso ) .* head_mask;
    % minus sign to get increasing phase over TEs
    phase_phi0_sd(:,:,:, echo) = ...
        (-totalField_rad_s * TE(echo) + phase_offset_est) .* head_mask; 
end %for

complexSignal = magnitude_sd .* exp(1i * phase_phi0_sd);

% Real and imaginary part of the complex signal
complexSignal_re = real(complexSignal);
complexSignal_im = imag(complexSignal);

% Adding gaussian noise to the real and imaginary parts of the complex
% signal
complexSignal_re_withNoise = complexSignal_re + ...
    random('Normal', noise_mean, noise_sd, [matrix_size n_echoes]);
complexSignal_im_withNoise = complexSignal_im + ...
    random('Normal', noise_mean, noise_sd, [matrix_size n_echoes]);

% Complex signal with noise
complexSignal_withNoise = ...
    msk .* (complexSignal_re_withNoise + 1i * complexSignal_im_withNoise);

% Magnitude and phase with noise
magnitude_sd_withNoise = abs(complexSignal_withNoise);
phase_phi0_sd_withNoise = angle(complexSignal_withNoise);

% Wrap phase without noise
phase_phi0_sd_wr = wrapToPi(phase_phi0_sd);

clear echo complexSignal_re complexSignal_im complexSignal_re_withNoise ...
    complexSignal_im_withNoise

%% Save 

save('magnitude_wholeHead_1iso_sd.mat', 'magnitude_sd_withNoise','magnitude_sd')
save('phase_wholeHead_1iso_sd.mat', 'phase_phi0_sd_withNoise','phase_phi0_sd_wr')
