function [R2_star,S0,TE_max] = R2star_lin(magnitude_data,phase_data,mask,TE)
%R2STAR_LIN Calculates R2* maps using a log linearisation of magnitude data
% In each voxel, the number of TEs to fit is determined based on a TE
% quality map derived from phase data (Storey et al. Artifact-free T2* 
% mapping without post hoc corrections, ISMRM 2015, abstract #442)
%
% INPUTS:
%   magnitude_data - 4D matrix of multi-echo magnitude data
%   phase_data - 4D matrix of multi-echo phase data [rad]
%   mask - binary mask of the region of interest (e.g. the brain)
%   TE - vector of TEs [s]
%
% OUTPUTS:
%   R2_star - R2star map [Hz]
%   S0 - S0 map
%   TE_max - map of maximum TE used for fitting
%
% AUTHOR: Emma Biondetti, University of Chieti-Pescara, Italy
% DATE: 01/03/2022

s = size(magnitude_data);
matrix_size = s(1:3);
n_echoes = s(4);

R2_star = zeros(matrix_size);
S0 = zeros(matrix_size);

% Matrix representation of the linear model
G = ones(n_echoes, 2);
G(:,1) = -2*TE(1:n_echoes);

%%% Calculating a map of TE_max for each voxel
brain_mask_p = padarray(mask, [1 1 1]);
phase_data_unw_p = padarray(phase_data, [1 1 1 0]);
alpha = 0.4;

% Calculating TE quality map
TE_max = zeros(matrix_size); % number of echoes to fit
for r = 1:matrix_size(1) %row
    for c = 1:matrix_size(2) %col
        for s = 1:matrix_size(3) %slice
            if mask(r, c, s) == 1 %only brain voxels
                
                active_voxel = [r+1 c+1 s+1]; %corresponding voxel in padded matrix
                
                % brain mask values in the neighbourhood of the current
                % voxel
                brain_mask_nd = brain_mask_p(active_voxel(1)-1:active_voxel(1)+1, ...
                    active_voxel(2)-1:active_voxel(2)+1, ...
                    active_voxel(3)-1:active_voxel(3)+1);
                brain_mask_nd = brain_mask_nd(:);
                brain_mask_nd(14) = [];
                
                echo = 1;
                phase_n_diff = 0;
                
                while ~any(phase_n_diff > alpha*pi) && echo < n_echoes
                    
                    % phase values in the neighbourhood of the current
                    % voxel
                    phase_neighbourhood = ...
                        phase_data_unw_p(active_voxel(1)-1:active_voxel(1)+1, ...
                        active_voxel(2)-1:active_voxel(2)+1, ...
                        active_voxel(3)-1:active_voxel(3)+1);
                    phase_neighbourhood = phase_neighbourhood(:);
                    phase_neighbourhood(14) = [];
                    phase_neighbourhood = ...
                        phase_neighbourhood(brain_mask_nd == 1);
                    
                    phase_active_voxel = ...
                        ones(size(phase_neighbourhood, 1), 1) * ...
                        phase_data_unw_p(active_voxel(1), ...
                        active_voxel(2), ...
                        active_voxel(3), echo);
                    phase_n_diff = ...
                        phase_active_voxel - phase_neighbourhood;
                    
                    echo = echo + 1;
                end % while
                TE_max(r,c,s)=echo;
                
            end %if
        end %for
    end %for
end %for

% Calculating R2* and S0 according to TE quality map
min_echoes_used = min(TE_max(:));
max_echoes_used = max(TE_max(:));
magnitude_data_log = log(magnitude_data);

for j=min_echoes_used:max_echoes_used
    
    % Get the voxels where this is the case
    ind_e = find(TE_max==j);
    
    % Get the log(magn_data) of those voxels
    S_e = zeros(j, length(ind_e));
    for t=1:j
        tmp = magnitude_data_log(:,:,:,t);
        S_e(t,:) = tmp(ind_e);
    end
    clear tmp
    
    % Fit data
    G_t = G(1:j,:);
    Q_t = (G_t'*G_t)\G_t';
    fit_e = Q_t*S_e;
    S0(ind_e) = exp(fit_e(2,:));
    R2_star(ind_e) = (fit_e(1,:));
    clear ind_e t S_e fit_e
    
end

R2_star = R2_star .* mask;
S0 = S0 .* mask;

R2_star(R2_star < 0) = 0;
R2_star(isinf(R2_star)) = nan; % can be excluded from average R2* calculations in ROIs

end

