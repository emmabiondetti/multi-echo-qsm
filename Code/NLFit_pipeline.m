function [susc_map_ppm, local_field_Hz, total_field_unw_rad, N_std, NewMask, alpha] = ...
    NLFit_pipeline(complex_data, brain_mask, voxel_size, delta_TE, ...
    B0, B0_dir)
% NLFIT_PIPELINE Applies the NLFit pipeline for susceptibility mapping
% AUTHOR: Emma Biondetti, University of Chieti-Pescara, Italy
% DATE: 07/03/2022

% 1) Non-linearly fit of the complex data over time
[total_field_rad, N_std, ~] = Fit_ppm_complex(complex_data);
s = size(complex_data);
matrix_size = s(1:3);

% 2) Spatial Laplacian phase unwrapping
total_field_unw_rad = ...
    unwrapLaplacian(total_field_rad, matrix_size, voxel_size);

% 3) Local field map estimation using V-SHARP
[local_field_rad, NewMask] = ...
    V_SHARP(total_field_unw_rad, brain_mask, 'voxelsize', voxel_size, ...
    'smvsize', 40);
local_field_Hz = local_field_rad/delta_TE/2/pi;

% 4) Susceptibility map estimation using Tikhonov regularisation
local_field_ppm = local_field_rad/(42.6*2*pi*delta_TE*B0);

par.Resolution = voxel_size;
par.Orientation = B0_dir;
[susc_map_ppm, alpha, ~] = directTikhonov_lcurve(local_field_ppm, NewMask, par);

end %function

