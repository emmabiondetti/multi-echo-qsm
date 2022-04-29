function [susc_map_ppm, local_field_SHARP_Hz, ...
    local_fields_SHARP_rad, alpha] = ...
    SNRwAvg_pipeline(phase_data, R2star_map, TE, brain_mask, ...
    voxel_size, B0, B0_dir, bm_ero, varargin)
% SNRWAVG_PIPELINE Applies the SNR_wAvg pipeline for susceptibility mapping
%
% INPUTS:
%   phase_data: 4D matrix with the multi-echo phase data
%   R2star_map: 3D matrix of R2* values [Hz]
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
%   local_field_SHARP_Hz: 3D map of local field values after multi-echo 
%       combination [Hz]
%   local_fields_SHARP_rad: 4D map of local field values [Hz]
%   alpha: regularisation parameter value used for Tikhonov inversion
%
% AUTHOR: Emma Biondetti, University of Chieti-Pescara, Italy
% DATE: 22/03/2022

if length(varargin) == 1
    par.Alpha = varargin{1};
end

% 1) Laplacian unwrapping at each TE
s = size(phase_data);
matrix_size = s(1:3);
n_echoes = s(4);
gamma = 42.6*2*pi*1e6; % [rad]/[s*T]

brain_mask_ero3vox = brainMask_erosion(brain_mask, bm_ero, '3D');

% 2) Local fields (SHARP)
local_fields_SHARP_rad = zeros([matrix_size n_echoes]);
for i = 1:n_echoes
    local_fields_SHARP_rad(:,:,:,i) = SHARP(phase_data(:,:,:,i), ...
        brain_mask_ero3vox, 0.05, voxel_size);
end

% 3) Local field combination
w_j = zeros([matrix_size n_echoes]);
for j=1:n_echoes
    w_j(:,:,:,j) = TE(j)*exp(-TE(j)*R2star_map);
end
w_j_sum = sum(w_j,4);

local_fields_w_T = zeros([matrix_size n_echoes]);
for i=1:n_echoes
    w_i = TE(i)*exp(-TE(i)*R2star_map)./w_j_sum;
    local_fields_w_T(:,:,:,i) = ...
        w_i.*local_fields_SHARP_rad(:,:,:,i)/(gamma*TE(i));
end
local_field_SHARP_T = sum(local_fields_w_T, 4);
local_field_SHARP_Hz = local_field_SHARP_T*gamma/(2*pi);

% 4) Susceptibility map estimation using Tikhonov regularisation
local_field_ppm = local_field_SHARP_T/(1e-6*B0);
local_field_ppm(isnan(local_field_ppm)) = 0;

par.Resolution = voxel_size;
par.Orientation = B0_dir;
[susc_map_ppm, alpha, ~] = ...
    directTikhonov_lcurve(local_field_ppm, brain_mask_ero3vox, par);

end %function

