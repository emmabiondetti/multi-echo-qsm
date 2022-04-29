% AUTHOR: Emma Biondetti, University of Chieti-Pescara, Italy
% DATE: 10/03/2022
% DESCRIPTION: this code calculates the numerical T2* phantom with both
% fixed or variables T2* distributions (_sd suffix)

clear; close all; clc

%% Loading original phantom data

matrix_size = [256 256 128]; % [mm]
voxel_size = [1.5 1.5 1.5]; % [mm]

fid = fopen('det_head_u2med.dat', 'r');
zubal_phantom = flip(fread(fid)); % flipped to match standard orientation
fclose(fid);
zubal_phantom = reshape(zubal_phantom, matrix_size);

clear ans fid

%% Upsampling
% Upsampling the original Zubal phantom using the nearest interpolation
% method to preserve the labels

zubal_phantom_1iso = imresize3(zubal_phantom, 1.5, "Method", "nearest");
matrix_size_1iso = matrix_size*1.5;

%% Label list

outside_phantom = 0;    
skin = 1;               % non-brain
csf = 2;                % brain - csf
spinal_cord = 3;        % brain - other
skull = 4;              % non-brain
spine = 5;              % non-brain
skeletal_muscle = 9;    % non-brain
pharynx = 15;           % non-brain
fat1 = 22;              % non-brain
blood_pool = 23;        % brain - venous blood
fornix = 26;            % brain - white matter
cartilage = 30;         % non-brain
dens_of_axis = 70;      % non-brain
jaw_bone = 71;          % non-brain
parotid_gland = 72;     % non-brain
lacrimal_glands = 74;   % non-brain
spinal_canal = 75;      % brain - csf
hard_palate = 76;       % non-brain
cerebellum = 77;        % brain - cerebellum grey matter
tongue = 78;            % non-brain
horn_of_mandible = 81;  % non-brain
nasal_septum = 82;      % non-brain
white_matter = 83;      % brain - white matter
superior_sagittal_sinus = 84; % brain - venous blood
medulla_oblungata = 85; % brain - other
frontal_lobes = 89;     % brain - grey matter
pons = 91;              % brain - other
third_ventricle = 92;   % spinal cord/brain - csf
occipital_lobes = 95;   % brain - grey matter
hippocampus = 96;       % brain
pituitary_gland = 97;   % brain
fat2 = 98;              % non-brain fat around the skull
uncus = 99;             % brain - grey matter
turbinates = 100;       % non-brain
caudate_nucleus = 101;  % brain
zygoma = 102;           % non-brain
insula_cortex = 103;    % brain - grey matter
sinuses_mouthCavity = 104;  % part brain (venous blood), part air (mouth and frontal sinuses)
putamen = 105;          % brain
optic_nerve = 106;      % non-brain
internal_capsule = 107; % brain - white matter
septum_pellucidum = 108;   % brain - both grey and white matter (intermediate value?)
thalamus = 109;         % brain 
eyeball = 110;          % non-brain
corpus_callosum = 111;  % brain - white matter
specialRegion_frontal_lobes = 112;  % brain - grey matter
cerebral_falx = 113;    % brain - classify it as CSF
temporal_lobes = 114;   % brain - grey matter
fourth_ventricle = 115; % brain - csf
frontal_portion_eyes = 116; % non-brain
parietal_lobes = 117;   % brain - grey matter
amygdala = 118;         % brain
eye = 119;              % non-brain
globus_pallidus = 120;  % brain
lens = 121;             % non-brain
cerebral_aquaduct = 122;    % brain - ventricular system - csf
lateral_ventricles = 123;   % brain - ventricular system - csf
prefrontal_lobes = 124; % brain - grey matter
teeth = 125;            % non-brain

%% Assigning labels to different groups

non_brain = [skin skull spine skeletal_muscle pharynx fat1 cartilage ...
    dens_of_axis jaw_bone parotid_gland lacrimal_glands hard_palate ...
    tongue horn_of_mandible nasal_septum fat2 turbinates zygoma ...
    optic_nerve eyeball frontal_portion_eyes eye lens teeth];

white_matter = [internal_capsule fornix white_matter corpus_callosum];

grey_matter = [cerebellum frontal_lobes occipital_lobes uncus ...
    insula_cortex specialRegion_frontal_lobes temporal_lobes ...
    parietal_lobes prefrontal_lobes];

cerebrospinal_fluid = [csf spinal_canal third_ventricle cerebral_falx ...
    fourth_ventricle cerebral_aquaduct lateral_ventricles];

venous_blood = superior_sagittal_sinus;

otherBrainRegions = [spinal_cord pons medulla_oblungata septum_pellucidum blood_pool];

%% T2 star values [s]

coeff_variation = 0.04;

T2_star_grey_matter_avg = 60 * 1e-3;
T2_star_grey_matter_sd = T2_star_grey_matter_avg*coeff_variation;

T2_star_white_matter_avg = 55 * 1e-3;
T2_star_white_matter_sd = T2_star_white_matter_avg*coeff_variation;

T2_star_venous_blood_avg = 21.2 * 1e-3;
T2_star_venous_blood_sd = T2_star_venous_blood_avg*coeff_variation;

T2_star_caudate_nucleus_avg = 54.8 * 1e-3;
T2_star_caudate_nucleus_sd = T2_star_caudate_nucleus_avg*coeff_variation;

T2_star_putamen_avg = 53.6 * 1e-3;
T2_star_putamen_sd = T2_star_putamen_avg*coeff_variation;

T2_star_thalamus_avg = 76.8 * 1e-3;
T2_star_thalamus_sd = T2_star_thalamus_avg*coeff_variation;

T2_star_globus_pallidus_avg = 28.8 * 1e-3;
T2_star_globus_pallidus_sd = T2_star_globus_pallidus_avg*coeff_variation;

T2_star_cerebrospinal_fluid_avg = 200 * 1e-3;
T2_star_cerebrospinal_fluid_sd = T2_star_cerebrospinal_fluid_avg*coeff_variation;

%% T2* value assignment

default_avg = 70 * 1e-3;
default_sd = default_avg * coeff_variation;

zubal_phantom_1iso_T2_star = ones(matrix_size_1iso)*default_avg; 
zubal_phantom_1iso_sd_T2_star = default_sd * randn(matrix_size_1iso) ...
    + default_avg; 

% Non-brain regions
n_non_brain_voxels = ...
    length(zubal_phantom_1iso(ismember(zubal_phantom_1iso, non_brain) == 1));
zubal_phantom_1iso_T2_star(ismember(zubal_phantom_1iso, non_brain)) = ...
    default_avg;
zubal_phantom_1iso_sd_T2_star(ismember(zubal_phantom_1iso, non_brain)) = ...
    default_sd*randn(n_non_brain_voxels,1) + default_avg;

zubal_phantom_1iso_T2_star(zubal_phantom_1iso == outside_phantom) = 0;
zubal_phantom_1iso_sd_T2_star(zubal_phantom_1iso == outside_phantom) = 0;

zubal_phantom_1iso_frontalSinuses_mouthCavity = zeros(matrix_size_1iso);
zubal_phantom_1iso_frontalSinuses_mouthCavity(zubal_phantom_1iso == ...
    sinuses_mouthCavity) = 1;
zubal_phantom_1iso_frontalSinuses_mouthCavity(:,1:round(1.5*125),:) = 0;
zubal_phantom_1iso_frontalSinuses_mouthCavity(:,:,1.5*62:end) = 0;

zubal_phantom_1iso_T2_star(...
    zubal_phantom_1iso_frontalSinuses_mouthCavity == 1) = 0;
zubal_phantom_1iso_sd_T2_star(...
    zubal_phantom_1iso_frontalSinuses_mouthCavity == 1) = 0;

zubal_phantom_1iso_frontalSinuses_mouthCavity(zubal_phantom_1iso == ...
    sinuses_mouthCavity) = 1;
zubal_phantom_1iso_frontalSinuses_mouthCavity...
    (:,1.5*126:matrix_size_1iso(2),:) = 0;

n_sinuses_mouthCavity_voxels = ...
    length(zubal_phantom_1iso_frontalSinuses_mouthCavity...
    (zubal_phantom_1iso_frontalSinuses_mouthCavity == 1));

zubal_phantom_1iso_T2_star(zubal_phantom_1iso_frontalSinuses_mouthCavity == 1) = ...
    default_avg; % default
zubal_phantom_1iso_sd_T2_star(zubal_phantom_1iso_frontalSinuses_mouthCavity == 1) = ...
    default_sd * randn(n_sinuses_mouthCavity_voxels,1) + default_avg; % default

% Grey matter
grey_matter_voxels = ismember(zubal_phantom_1iso(:), grey_matter);
n_grey_matter_voxels = length(grey_matter_voxels(grey_matter_voxels == 1));

zubal_phantom_1iso_T2_star(ismember(zubal_phantom_1iso, grey_matter)) = ...
    T2_star_grey_matter_avg;
zubal_phantom_1iso_sd_T2_star(ismember(zubal_phantom_1iso, grey_matter)) = ...
    T2_star_grey_matter_sd*randn(n_grey_matter_voxels,1)+T2_star_grey_matter_avg;

% White matter
white_matter_voxels = ismember(zubal_phantom_1iso(:), white_matter);
n_white_matter_voxels = length(white_matter_voxels(white_matter_voxels == 1));

zubal_phantom_1iso_T2_star(ismember(zubal_phantom_1iso, white_matter)) = ...
    T2_star_white_matter_avg;
zubal_phantom_1iso_sd_T2_star(ismember(zubal_phantom_1iso, white_matter)) = ...
    T2_star_white_matter_sd*randn(n_white_matter_voxels,1)+T2_star_white_matter_avg;

% Cerebrospinal fluid
cerebrospinal_fluid_voxels = ismember(zubal_phantom_1iso(:), cerebrospinal_fluid);
n_cerebrospinal_fluid_voxels = length(cerebrospinal_fluid_voxels(cerebrospinal_fluid_voxels == 1));

zubal_phantom_1iso_T2_star(ismember(zubal_phantom_1iso, cerebrospinal_fluid)) = ...
    T2_star_cerebrospinal_fluid_avg;
zubal_phantom_1iso_sd_T2_star(ismember(zubal_phantom_1iso, cerebrospinal_fluid)) = ...
    T2_star_cerebrospinal_fluid_sd*randn(n_cerebrospinal_fluid_voxels,1)+T2_star_cerebrospinal_fluid_avg;

% Venous blood
venous_blood_voxels = ismember(zubal_phantom_1iso(:), venous_blood);
n_venous_blood_voxels = length(venous_blood_voxels(venous_blood_voxels == 1));

zubal_phantom_1iso_T2_star(ismember(zubal_phantom_1iso, venous_blood)) = ...
    T2_star_venous_blood_avg;
zubal_phantom_1iso_sd_T2_star(ismember(zubal_phantom_1iso, venous_blood)) = ...
    T2_star_venous_blood_sd*randn(n_venous_blood_voxels,1)+T2_star_venous_blood_avg;

% Caudate nucleus
caudate_nucleus_voxels = ismember(zubal_phantom_1iso(:), caudate_nucleus);
n_caudate_nucleus_voxels = length(caudate_nucleus_voxels(caudate_nucleus_voxels == 1));

zubal_phantom_1iso_T2_star(ismember(zubal_phantom_1iso, caudate_nucleus)) = ...
    T2_star_caudate_nucleus_avg;
zubal_phantom_1iso_sd_T2_star(ismember(zubal_phantom_1iso, caudate_nucleus)) = ...
    T2_star_caudate_nucleus_sd*randn(n_caudate_nucleus_voxels,1)+T2_star_caudate_nucleus_avg;

% Putamen
putamen_voxels = ismember(zubal_phantom_1iso(:), putamen);
n_putamen_voxels = length(putamen_voxels(putamen_voxels == 1));

zubal_phantom_1iso_T2_star(ismember(zubal_phantom_1iso, putamen)) = ...
    T2_star_putamen_avg;
zubal_phantom_1iso_sd_T2_star(ismember(zubal_phantom_1iso, putamen)) = ...
    T2_star_putamen_sd*randn(n_putamen_voxels,1)+T2_star_putamen_avg;

% Thalamus
thalamus_voxels = ismember(zubal_phantom_1iso(:), thalamus);
n_thalamus_voxels = length(thalamus_voxels(thalamus_voxels == 1));

zubal_phantom_1iso_T2_star(ismember(zubal_phantom_1iso, thalamus)) = ...
    T2_star_thalamus_avg;
zubal_phantom_1iso_sd_T2_star(ismember(zubal_phantom_1iso, thalamus)) = ...
    T2_star_thalamus_sd*randn(n_thalamus_voxels,1)+T2_star_thalamus_avg;

% Globus pallidus
globus_pallidus_voxels = ismember(zubal_phantom_1iso(:), globus_pallidus);
n_globus_pallidus_voxels = length(globus_pallidus_voxels(globus_pallidus_voxels == 1));

zubal_phantom_1iso_T2_star(ismember(zubal_phantom_1iso, globus_pallidus)) = ...
    T2_star_globus_pallidus_avg;
zubal_phantom_1iso_sd_T2_star(ismember(zubal_phantom_1iso, globus_pallidus)) = ...
    T2_star_globus_pallidus_sd*randn(n_globus_pallidus_voxels,1)+T2_star_globus_pallidus_avg;

%% Reshaping mouth cavity with ellipsoid

[x,y,z] = meshgrid(1:matrix_size_1iso(1),1:matrix_size_1iso(2),1:matrix_size_1iso(3));
x0 = 210*1.5;
y0 = 130*1.5;
z0 = 32*1.5;
rx = round(55*1.5);
ry = 32*1.5;
rz = 22*1.5;

% Inside the ellipsoid
ell = ((x - x0)/rx).^2 + ((y - y0)/ry).^2 + ((z - z0)/rz).^2 <= 1;

zubal_phantom_1iso_T2_star_2 = zubal_phantom_1iso_T2_star + ell;
zubal_phantom_1iso_sd_T2_star_2 = zubal_phantom_1iso_sd_T2_star + ell;

n_v = length(zubal_phantom_1iso_T2_star_2(zubal_phantom_1iso_T2_star_2 == 0));
zubal_phantom_1iso_T2_star_2(zubal_phantom_1iso_T2_star_2 == 0) = ...
   default_avg;
zubal_phantom_1iso_sd_T2_star_2(zubal_phantom_1iso_sd_T2_star_2 == 0) = ...
    default_sd*randn(n_v,1)+default_avg;

zubal_phantom_1iso_sd_T2_star_2(zubal_phantom_1iso == outside_phantom) = 0; % background
zubal_phantom_1iso_sd_T2_star_2(zubal_phantom_1iso_sd_T2_star_2 == 1) = 0;

zubal_phantom_1iso_T2_star_2(zubal_phantom_1iso == outside_phantom) = 0; % background
zubal_phantom_1iso_T2_star_2(zubal_phantom_1iso_T2_star_2 == 1) = 0;


zubal_phantom_1iso_T2_star_2(and(zubal_phantom_1iso_T2_star_2 > 1.06,...
    zubal_phantom_1iso_T2_star_2 < 1.08)) = default_avg;

n_v2 = length(zubal_phantom_1iso_sd_T2_star_2(and(zubal_phantom_1iso_sd_T2_star_2 > 1.06,...
    zubal_phantom_1iso_sd_T2_star_2 < 1.08)));
zubal_phantom_1iso_sd_T2_star_2(and(zubal_phantom_1iso_sd_T2_star_2 > 1.06,...
    zubal_phantom_1iso_sd_T2_star_2 < 1.08)) = ...
    default_sd*randn(n_v2,1)+default_avg;

%% Plot phantom

x_coord = round(127*1.5);
y_coord = round(137*1.5);
z_coord = 72*1.5;

figure
[ha, ~] = tight_subplot(2, 2, [0 0], [0 0], [0.01 0.01]);

axes(ha(1))
imagesc(rot90(squeeze(zubal_phantom_1iso_T2_star_2(:, y_coord, :))), [0 .1])
axis off, daspect(voxel_size), colormap('gray')

axes(ha(2))
imagesc(rot90(squeeze(zubal_phantom_1iso_T2_star_2(x_coord, :, :))), [0 .1])
axis off, daspect(voxel_size), colormap('gray')

axes(ha(3))
imagesc(rot90(squeeze(zubal_phantom_1iso_T2_star_2(:, :, z_coord))), [0 .1])
axis off, daspect(voxel_size), colormap('gray')

axes(ha(4))
c = colorbar;
c.Location = 'South';
c.Ticks = [0, 0.5, 1];
c.TickLabels = ({'0', '[a.u.]', '1'});
text(.3, .7, ['X = ' num2str(x_coord) '\newline' 'Y = ' ...
    num2str(y_coord) '\newline' 'Z = ' num2str(z_coord)], 'FontSize', 12)
axis off
title('1-mm iso with constant T2*')
set(gca, 'FontSize', 12)

figure
[ha, pos] = tight_subplot(2, 2, [0 0], [0 0], [0.01 0.01]);

axes(ha(1))
imagesc(rot90(squeeze(zubal_phantom_1iso_sd_T2_star_2(:, y_coord, :))), [0 .1])
axis off, daspect(voxel_size), colormap('gray')

axes(ha(2))
imagesc(rot90(squeeze(zubal_phantom_1iso_sd_T2_star_2(x_coord, :, :))), [0 .1])
axis off, daspect(voxel_size), colormap('gray')

axes(ha(3))
imagesc(rot90(squeeze(zubal_phantom_1iso_sd_T2_star_2(:, :, z_coord))), [0 .1])
axis off, daspect(voxel_size), colormap('gray')

axes(ha(4))
c = colorbar;
c.Location = 'South';
c.Ticks = [0, 0.5, 1];
c.TickLabels = ({'0', '[a.u.]', '1'});
text(.3, .7, ['X = ' num2str(x_coord) '\newline' 'Y = ' ...
    num2str(y_coord) '\newline' 'Z = ' num2str(z_coord)], 'FontSize', 12)
axis off
title('1-mm iso with variable T2*')
set(gca, 'FontSize', 12)

%% Save phantom

nii1 = make_nii(zubal_phantom_1iso_T2_star_2, [1 1 1], matrix_size_1iso/2+1);
save_nii(nii1, 'zubalPhantom_T2_star_wholeHead_ellipsoid_1iso.nii')

nii2 = make_nii(zubal_phantom_1iso_sd_T2_star_2, [1 1 1], matrix_size_1iso/2+1);
save_nii(nii2, 'zubalPhantom_T2_star_wholeHead_ellipsoid_1iso_sd.nii')
