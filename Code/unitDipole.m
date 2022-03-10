function [D, D_centre] = ...
    unitDipole(dipoleSize, voxelSize, dipole_centre, B0_dir)
% UNITDIPOLE Creates a unit magnetic dipole in k-space
%
% INPUT
% - dipoleSize: size of unit magnetic dipole
% - voxelSize: voxel size [mm mm mm]
% - dipoleCentre: value at the centre of the dipole (k_x = k_y = k_z = 0)
% - B0_dir: normalised B0 direction calculated as B0 = R \ [0 0 1]' then 
% B0 = B0 / norm(B0), with R rotation matrix from the image to the scanner 
% coordinates
%
% OUTPUT
% - D: unit magnetic dipole in the Fourier domain
% - D_centre: dipole centre coordinates
%
% AUTHOR: Emma Biondetti, University of Chieti-Pescara, Italy
% DATE: 10/03/2022

% Evaluating unit magnetic dipole over gridded domain
% x: columns of k-space
% y: rows of k-space
% z: across slices of k-space
[kx, ky, kz] = ndgrid(-dipoleSize(1)/2 : dipoleSize(1)/2-1, ...
    -dipoleSize(2)/2 : dipoleSize(2)/2-1, ...
    -dipoleSize(3)/2 : dipoleSize(3)/2-1);

% Scaling dipole dimensions to voxel size
kx = kx / (dipoleSize(1) * voxelSize(1));
ky = ky / (dipoleSize(2) * voxelSize(2));
kz = kz / (dipoleSize(3) * voxelSize(3));

% Evaluating dipole kernel over gridded k-space, rotating according to the
% B0 direction
k_squared = kx.^2 + ky.^2 + kz.^2;
D = 1/3 - (kx .* B0_dir(1) + ky .* B0_dir(2) + kz .* B0_dir(3)).^2 ./ k_squared;
[D_centre(1), D_centre(2), D_centre(3)] = ind2sub(size(D), find(isnan(D)));
D(D_centre(1), D_centre(2), D_centre(3)) = dipole_centre; % k_x = k_y = k_z = 0

end % function