function [susc_map_ppm, local_field_SHARP_Hz, ...
    total_field_unw_rad, N_std, alpha] = ...
    NLFit_pipeline(complex_data, brain_mask, voxel_size, delta_TE, ...
    B0, B0_dir, bm_ero, varargin)
% NLFIT_PIPELINE Applies the NLFit pipeline for susceptibility mapping
%
% INPUTS:
%   complex_data: 4D matrix with the multi-echo complex data 
%       complex = (magnitude*exp(1i*phase))
%   brain_mask: binary brain mask
%   voxel_size: 3-element vector with the voxel resolution [mm]
%   delta_TE: echo time difference [s]
%   B0: static field strength [T]
%   B0_dir: 3-element vector with the direction of the B0 field 
%       (e.g. [0 0 1]' for the z-axis)
%   bm_ero: 3-element vector with the number of voxels by which the mask 
%       must be eroded in each dimension [ero_x ero_y ero_z]
%   varargin: if provided, regularisation parameter value for Tikhonov
%       inversion
%
% OUTPUTS:
%   susc_map_ppm: 3D map of magnetic susceptibility values [ppm]
%   local_field_SHARP_Hz: 3D map of local field values [Hz]
%   total_field_unw_rad: 3D map of unwrapped field values [Hz]
%   N_std: 3D map of total field noise values
%   alpha: regularisation parameter value used for Tikhonov inversion
%
% AUTHOR: Emma Biondetti, University of Chieti-Pescara, Italy
% DATE: 21/03/2022

if length(varargin) == 1
    par.Alpha = varargin{1};
end

% 1) Non-linearly fit of the complex data over time
[total_field_rad, N_std, ~] = Fit_ppm_complex(complex_data);
s = size(complex_data);
matrix_size = s(1:3);

% 2) Spatial Laplacian phase unwrapping
total_field_unw_rad = ...
    unwrapLaplacian(total_field_rad, matrix_size, voxel_size);
brain_mask_ero3vox = brainMask_erosion(brain_mask, bm_ero);

% 3) Local field map estimation (SHARP + V-SHARP)
local_field_SHARP_rad = SHARP(total_field_unw_rad, brain_mask_ero3vox,...
    0.05, voxel_size);
local_field_SHARP_Hz = local_field_SHARP_rad/delta_TE/2/pi;

% 4) Susceptibility map estimation using Tikhonov regularisation
local_field_ppm = local_field_SHARP_rad/(42.6*2*pi*delta_TE*B0);

par.Resolution = voxel_size;
par.Orientation = B0_dir;
[susc_map_ppm, alpha, ~] = ...
    directTikhonov_lcurve(local_field_ppm, brain_mask_ero3vox, par);

end %function

