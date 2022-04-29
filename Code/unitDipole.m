function [D, D_centre] = ...
    unitDipole(dipole_size, voxel_size, dipole_centre, B0_dir)
% UNITDIPOLE Creates a unit magnetic dipole in k-space
%
% INPUTS:
%   dipole_size: 3-element vector with the size of the unit magnetic dipole
%   voxel_size: 3-element vector with the voxel resolution [mm]
%   dipoleCentre: value at the centre of the dipole (k_x = k_y = k_z = 0)
%   B0_dir: 3-element vector with the direction of the B0 field 
%       (e.g. [0 0 1]' for the z-axis)
%
% OUTPUTS:
%   D: unit magnetic dipole in the Fourier domain
%   D_centre: dipole centre coordinates
%
% AUTHOR: Emma Biondetti, University of Chieti-Pescara, Italy
% DATE: 10/03/2022

% Evaluating unit magnetic dipole over gridded domain
% x: columns of k-space
% y: rows of k-space
% z: across slices of k-space
[kx, ky, kz] = ndgrid(-dipole_size(1)/2 : dipole_size(1)/2-1, ...
    -dipole_size(2)/2 : dipole_size(2)/2-1, ...
    -dipole_size(3)/2 : dipole_size(3)/2-1);

% Scaling dipole dimensions to voxel size
kx = kx / (dipole_size(1) * voxel_size(1));
ky = ky / (dipole_size(2) * voxel_size(2));
kz = kz / (dipole_size(3) * voxel_size(3));

% Evaluating dipole kernel over gridded k-space, rotating according to the
% B0 direction
k_squared = kx.^2 + ky.^2 + kz.^2;
D = 1/3 - (kx .* B0_dir(1) + ky .* B0_dir(2) + kz .* B0_dir(3)).^2 ./ k_squared;
[D_centre(1), D_centre(2), D_centre(3)] = ind2sub(size(D), find(isnan(D)));
D(D_centre(1), D_centre(2), D_centre(3)) = dipole_centre;

end % function