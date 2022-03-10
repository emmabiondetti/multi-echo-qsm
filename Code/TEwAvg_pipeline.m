function [susc_map_ppm, local_field_Hz, total_field_unw_rad_s, ...
    phase_data_unw, NewMask, alpha] = ...
    TEwAvg_pipeline(phase_data, TE, brain_mask, voxel_size, B0, B0_dir)
% TEWAVG_PIPELINE Applies the TE_wAvg pipeline for susceptibility mapping
%
% @AUTHOR: Emma Biondetti, University of Chieti-Pescara, Italy
% @DATE: 07/03/2022

% 1) Laplacian unwrapping at each TE
s = size(phase_data);
matrix_size = s(1:3);
n_echoes = s(4);

phase_data_unw = zeros([matrix_size n_echoes]);
for i = 1:n_echoes
    phase_data_unw(:,:,:,i) = ...
        unwrapLaplacian(phase_data(:,:,:,i), matrix_size, voxel_size);
end %for

% 2) Total field estimation
total_field_unw_rad_s = sum(phase_data_unw,4)/sum(TE);

% 3) Local field map estimation using V-SHARP
[local_field_rad_s, NewMask] = ...
    V_SHARP(total_field_unw_rad_s, brain_mask, 'voxelsize', voxel_size, ...
    'smvsize', 40);
local_field_Hz = local_field_rad_s/2/pi;

% 4) Susceptibility map estimation using Tikhonov regularisation
local_field_ppm = local_field_rad_s/(42.6*2*pi*B0);

par.Resolution = voxel_size;
par.Orientation = B0_dir;
[susc_map_ppm, alpha, ~] = ...
    directTikhonov_lcurve(local_field_ppm, NewMask, par);

end %function

