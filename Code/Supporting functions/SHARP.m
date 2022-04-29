function y = SHARP(phase_data, brain_mask, sigma, voxel_size)
% SHARP Unwraps phase and removes background fields using a Laplacian-based
% approach
%
% INPUTS:
%   phase_data: 4D matrix with the multi-echo phase data
%   brain_mask: binary brain mask if SHARP is used for background 
%       field removal; a 3D matrix of ones of the same size as phase_data 
%       if SHARP is used for phase unerapping
%   sigma: a threshold to regularise laplacian kernel inversion
%   voxel_size: 3-element vector with the voxel resolution [mm]
%   mode: 2D for slice-by-slice applications or 3D for volume applications
%
% OUTPUT:
%   y: SHARP-processed phase image 
%
% AUTHOR: Emma Biondetti, University of Chieti-Pescara, Italy
% DATE: 05/07/2016

matrix_size = squeeze(size(phase_data));

% Calculating the Laplacian Kernel
kernel_size = matrix_size;
kernel = laplacianKernel_3D(kernel_size, voxel_size);

% Regularised inversion of Laplacian kernel
inv_kernel = 1./kernel;
inv_kernel(abs(kernel) <= sigma) = 0;

% Applying the forward Laplacian kernel
first_term = ifftn(kernel .* fftn(sin(phase_data), kernel_size));
second_term = ifftn(kernel .* fftn(cos(phase_data), kernel_size));

first_term = cos(phase_data) .* first_term;
second_term = sin(phase_data) .* second_term;

yy = brain_mask .* (first_term - second_term);

% Applying the inverse Laplacian kernel
y = brain_mask .* (ifftn(inv_kernel .* fftn(yy, kernel_size)));

end 
