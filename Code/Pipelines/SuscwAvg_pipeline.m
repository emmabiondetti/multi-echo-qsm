function [susc_map_ppm, susc_maps_ppm, local_fields_SHARP_Hz, alpha] = ...
    SuscwAvg_pipeline(phase_data, magnitude_data, TE, brain_mask, ...
    voxel_size, B0, B0_dir, bm_ero, varargin)
% TEWAVG_PIPELINE Applies the TE_wAvg pipeline for susceptibility mapping
%
% INPUTS:
%   phase_data: 4D matrix with the multi-echo phase data
%   magnitude_data: 4D matrix with the multi-echo magnitude data
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
%   susc_maps_ppm: 4D map of magnetic susceptibility values before 
%       multi-echo combination [ppm]
%   local_fields_SHARP_Hz: 4D map of local field values [Hz]
%   alpha: regularisation parameter value used for Tikhonov inversion
%
% AUTHOR: Emma Biondetti, University of Chieti-Pescara, Italy
% DATE: 31/03/2022

if length(varargin) == 1
    alpha_v = varargin{1};    
end

% 1) Laplacian unwrapping + background field removal at each TE
s = size(phase_data);
matrix_size = s(1:3);
n_echoes = s(4);

brain_mask_ero3vox = brainMask_erosion(brain_mask, bm_ero, '3D');
local_fields_SHARP_rad = zeros([matrix_size n_echoes]);
for i = 1:n_echoes
    local_fields_SHARP_rad(:,:,:,i) = SHARP(phase_data(:,:,:,i), ...
        brain_mask_ero3vox, 0.05, voxel_size);
    local_fields_SHARP_Hz = local_fields_SHARP_rad/(TE(i)*2*pi);
end

% 2) Susceptibility map estimation using Tikhonov regularisation
par.Resolution = voxel_size;
par.Orientation = B0_dir;
susc_maps_ppm = zeros([matrix_size n_echoes]);
alpha = zeros(n_echoes,1);
for i = 1:n_echoes
    local_field_ppm = local_fields_SHARP_rad(:,:,:,i)/(42.6*2*pi*B0*TE(i));
    if length(varargin) == 1
        par.Alpha = alpha_v(i);
        [susc_maps_ppm(:,:,:,i), alpha(i), ~] = ...
            directTikhonov_lcurve(local_field_ppm, brain_mask_ero3vox, par);
    else
        [susc_maps_ppm(:,:,:,i), alpha(i), ~] = ...
            directTikhonov_lcurve(local_field_ppm, brain_mask_ero3vox, par);
    end
end

% 3) Susceptibility map combination
w_j = zeros([matrix_size n_echoes]);
for j=1:n_echoes
    w_j(:,:,:,j) = magnitude_data(:,:,:,j).^2*TE(j)^2;
end
w_j_sum = sum(w_j,4);

susc_maps_w_ppm = zeros([matrix_size n_echoes]);
for i=1:n_echoes
    w_i = magnitude_data(:,:,:,i).^2*TE(i)^2./w_j_sum;
    susc_maps_w_ppm(:,:,:,i) = w_i.*susc_maps_ppm(:,:,:,i);
end
susc_map_ppm = brain_mask_ero3vox.*sum(susc_maps_w_ppm,4);
susc_map_ppm(isnan(susc_map_ppm)) = 0;

end %function

