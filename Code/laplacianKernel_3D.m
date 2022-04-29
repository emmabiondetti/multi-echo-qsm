function K = laplacianKernel_3D(kernel_size, voxel_size)
% LAPLACIANKERNEL Creates a Laplacian kernel in k-space
%
% INPUTS:
%   kernel_size: 3-element vector with the size of the kernel
%   voxel_size: 3-element vector with the voxel resolution [mm]
%
% OUTPUT:
%   K: 3D Laplacian kernel in k-space
%
% AUTHOR: Emma Biondetti, University of Chieti-Pescara, Italy
% DATE: 17/12/2015

k = zeros(kernel_size);

c1 = kernel_size(1)/2;
c2 = kernel_size(2)/2;
c3 = kernel_size(3)/2;

k(c1+1, c2+1, c3+1) = ...
    - 2 / voxel_size(1)^2 - 2 / voxel_size(2)^2 - 2 / voxel_size(3)^2;
k(c1,   c2+1, c3+1) = 1 / voxel_size(1)^2;
k(c1+2, c2+1, c3+1) = 1 / voxel_size(1)^2;
k(c1+1, c2,   c3+1) = 1 / voxel_size(2)^2;
k(c1+1, c2+2, c3+1) = 1 / voxel_size(2)^2;
k(c1+1, c2+1, c3)   = 1 / voxel_size(3)^2;
k(c1+1, c2+1, c3+2) = 1 / voxel_size(3)^2;

K = fftn(ifftshift(k));

end %function