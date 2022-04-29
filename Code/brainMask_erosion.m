function y = brainMask_erosion(brain_mask, erode_size, mode)
%BRAINMASK_EROSION Erodes a binary brain mask by the specified size
%
% INPUTS:
%   brain_mask: the original 3D brain mask (binary)
%   erode_size: 3-element vector of the number of voxels by which the mask 
%       must be eroded in each dimension
%   mode: 2D for slice-by-slice erosion or 3D for volume erosion
%
% OUTPUT:
%   y: the eroded brain mask
%
% AUTHOR: Emma Biondetti, University of Chieti-Pescara, Italy
% DATE: 30/03/2016

% Erosion across image columns
brainMask_eroded_tmp = ...
    imerode(brain_mask, strel('line', erode_size(2), 0));

% Erosion across image rows
brainMask_eroded_tmp = ...
    imerode(brainMask_eroded_tmp, strel('line', erode_size(1), 90));

if strcmp(mode, '2D')
    
    y = brainMask_eroded_tmp;
    
elseif strcmp(mode, '3D')
    
    % 2nd (columns) and 3rd (slices) dimensions are switched
    brainMask_eroded_tmp = permute(brainMask_eroded_tmp, [1,3,2]);
    
    % Erosion across image columns, i.e. slices
    brainMask_eroded_tmp = imerode(brainMask_eroded_tmp, strel('line', erode_size(3), 0));
    
    % 2nd (slices) and 3rd (columns) dimensions are switched
    brainMask_eroded = permute(brainMask_eroded_tmp, [1,3,2]);
    
    y = brainMask_eroded;
    
end 

end 

