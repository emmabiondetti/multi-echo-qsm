function [susc_map_ppm, local_field_SHARP_Hz, ...
    total_field_unw_rad_s, phase_data_unw, alpha] = ...
    TEwAvg_pipeline(phase_data, TE, brain_mask, voxel_size, B0, B0_dir, ...
    bm_ero, varargin)
% TEWAVG_PIPELINE Applies the TE_wAvg pipeline for susceptibility mapping
%
% INPUTS:
%   phase_data: 4D matrix with the multi-echo phase data
%   TE: vector with the echo time values [s]
%   brain_mask: binary brain mask
%   voxel_size: 3-element vector with the voxel resolution [mm]
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
%   phase_data_unw: 3D map of unwrapped phase data values [rad]
%   alpha: regularisation parameter value used for Tikhonov inversion
%
% AUTHOR: Emma Biondetti, University of Chieti-Pescara, Italy
% DATE: 21/03/2022

if length(varargin) == 1
    par.Alpha = varargin{1};
end

% 1) Laplacian unwrapping at each TE
s = size(phase_data);
matrix_size = s(1:3);
n_echoes = s(4);

phase_data_unw = zeros([matrix_size n_echoes]);
for i = 1:n_echoes
    phase_data_unw(:,:,:,i) = ...
        SHARP(phase_data(:,:,:,i), ones(matrix_size), ...
        1e-10, voxel_size);
end %for

% 2) Total field estimation
total_field_unw_rad_s = sum(phase_data_unw, 4)/sum(TE);

brain_mask_ero3vox = brainMask_erosion(brain_mask, bm_ero, '3D');

% 3) Local field map estimation (SHARP + V-SHARP)
delta_TE = TE(2) - TE(1);
local_field_SHARP_rad = SHARP(total_field_unw_rad_s*delta_TE, brain_mask_ero3vox,...
    0.05, voxel_size);
local_field_SHARP_Hz = local_field_SHARP_rad/delta_TE/2/pi;

% 4) Susceptibility map estimation using Tikhonov regularisation
local_field_ppm = local_field_SHARP_rad/(delta_TE*42.6*2*pi*B0);

par.Resolution = voxel_size;
par.Orientation = B0_dir;
[susc_map_ppm, alpha, ~] = ...
    directTikhonov_lcurve(local_field_ppm, brain_mask_ero3vox, par);

end %function

