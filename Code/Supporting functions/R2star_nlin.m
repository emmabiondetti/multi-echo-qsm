function [R2star,S0] = R2star_nlin(magnitude_data,mask,TE,R2_star_map_lin,S0_map_lin)
%R2STAR_NLIN Calculates R2* maps based on a nonlinear monoexponential fit
%
% INPUTS:
%   magnitude_data - 4D matrix of multi-echo magnitude data
%   mask - binary mask of the region of interest (e.g. the brain)
%   TE - vector of TEs [s]
%   R2_star_map_lin - map of initial R2star estimates calculated using 
%       the R2star_lin function 
%   S0_map_lin - map of initial S0 estimates calculated using 
%       the R2star_lin function 
%
% OUTPUTS:
%   R2_star - R2star map [Hz]
%   S0 - S0 map
%
% AUTHOR: Emma Biondetti, University of Chieti-Pescara, Italy
% DATE: 01/03/2022

sz = size(magnitude_data);
matrix_size = sz(1:3);
n_echoes = sz(4);

S0_map_lin(isnan(S0_map_lin))=0;
R2_star_map_lin(isnan(R2_star_map_lin))=0;

modelfun = @(b,x)(b(1)*exp(-x*b(2)));
opts = statset('nlinfit');
opts.FunValCheck = 'off'; % do not check for nan or inf values

jj = find(mask);
mag = reshape(magnitude_data,[prod(matrix_size) n_echoes]);
s0 = S0_map_lin(:);
r2 = R2_star_map_lin(:);
s0_values = s0(jj);
r2_values = r2(jj);
mag_values = mag(jj,:);
s0_nlin = zeros(numel(jj),1);
r2_nlin = zeros(numel(jj),1);

tic
parfor m = 1:numel(jj)
    s = mag_values(m,:);
    p0 = [s0_values(m), r2_values(m)];
    [cf_] = nlinfit(TE,s,modelfun,p0,opts);
    s0_nlin(m) = cf_(1);
    r2_nlin(m) = cf_(2);
end
toc
clear pool

tmp1 = zeros(matrix_size);
tmp1 = tmp1(:);
tmp1(jj) = r2_nlin;
R2star = reshape(tmp1, [matrix_size(1) matrix_size(2) matrix_size(3)]);

tmp2 = zeros(matrix_size);
tmp2 = tmp2(:);
tmp2(jj) = s0_nlin;
S0 = reshape(tmp2, [matrix_size(1) matrix_size(2) matrix_size(3)]);

end

