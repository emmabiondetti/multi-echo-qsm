% AUTHOR: Emma Biondetti, University of Chieti-Pescara, Italy
% DATE: 10/03/2022
% DESCRIPTION: this code calculates the  numerical susceptibility phantom 
% with both fixed or variables susceptibility distributions (_sd suffix)

clear
close all
clc

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

%% Magnetic susceptibility values [ppm]

coeff_variation = 0.2;

susc_grey_matter_avg = 0;
susc_grey_matter_sd = 0.02;

susc_white_matter_avg = -0.05;
susc_white_matter_sd = abs(susc_white_matter_avg)*coeff_variation;

susc_venous_blood_avg = 0.3;
susc_venous_blood_sd = susc_venous_blood_avg*coeff_variation;

susc_caudate_nucleus_avg = 0.09;
susc_caudate_nucleus_sd = susc_caudate_nucleus_avg*coeff_variation;

susc_putamen_avg = 0.09;
susc_putamen_sd = susc_putamen_avg*coeff_variation;

susc_thalamus_avg = 0.07;
susc_thalamus_sd = susc_thalamus_avg*coeff_variation;

susc_globus_pallidus_avg = 0.19;
susc_globus_pallidus_sd = susc_globus_pallidus_avg*coeff_variation;

susc_cerebrospinal_fluid_avg = 0;
susc_cerebrospinal_fluid_sd = 0.02;

susc_hippocampus_avg = 0;
susc_hippocampus_sd = 0.02;

susc_amygdala_avg = 0;
susc_amygdala_sd = 0.02;

susc_pituitary_gland_avg = 0;
susc_pituitary_gland_sd = 0.02;

% Background fields aren't well removed with 9.4 ppm
susc_background_avg = 9.4;
susc_background_sd = susc_background_avg*coeff_variation;

%% Susceptibility phantom: 1-mm iso

zubal_phantom_1iso_susceptibility = zeros(matrix_size_1iso);
zubal_phantom_1iso_sd_susceptibility = 0.02*randn(matrix_size_1iso);

zubal_phantom_1iso_susceptibility(zubal_phantom_1iso == outside_phantom) = ...
    susc_background_avg; % background [ppm]
zubal_phantom_1iso_sd_susceptibility(zubal_phantom_1iso == outside_phantom) = ...
    susc_background_avg; % background [ppm]

zubal_phantom_frontalSinuses_mouthCavity = zeros(matrix_size_1iso);
zubal_phantom_frontalSinuses_mouthCavity(zubal_phantom_1iso == sinuses_mouthCavity) = 1;

x = zubal_phantom_frontalSinuses_mouthCavity(:,1:round(1.5*125),:);
x_n = (x == 1);
t = 0.02*randn(matrix_size_1iso);
zubal_phantom_frontalSinuses_mouthCavity(:,1:round(1.5*125),:) = ...
    t(:,1:round(1.5*125),:) .* x_n;

x2 = zubal_phantom_frontalSinuses_mouthCavity(:,:,1.5*62:end);
x2_n = (x2 == 1);
zubal_phantom_frontalSinuses_mouthCavity(:,:,1.5*62:end) = ...
    t(:,:,1.5*62:end) .* x2_n;

zubal_phantom_1iso_sd_susceptibility(zubal_phantom_frontalSinuses_mouthCavity == 1) = ...
    susc_background_avg; % background [ppm]
zubal_phantom_1iso_susceptibility(zubal_phantom_frontalSinuses_mouthCavity == 1) = ...
    susc_background_avg; % background [ppm]

n_grey_matter_voxels = ...
    length(zubal_phantom_1iso_sd_susceptibility(ismember(zubal_phantom_1iso, grey_matter)));
zubal_phantom_1iso_susceptibility(ismember(zubal_phantom_1iso, grey_matter)) = ...
    susc_grey_matter_avg;
zubal_phantom_1iso_sd_susceptibility(ismember(zubal_phantom_1iso, grey_matter)) = ...
    susc_grey_matter_sd*randn(n_grey_matter_voxels,1) + susc_grey_matter_avg;

n_white_matter_voxels = ...
    length(zubal_phantom_1iso_sd_susceptibility(ismember(zubal_phantom_1iso, white_matter)));
zubal_phantom_1iso_susceptibility(ismember(zubal_phantom_1iso, white_matter)) = ...
    susc_white_matter_avg;
zubal_phantom_1iso_sd_susceptibility(ismember(zubal_phantom_1iso, white_matter)) = ...
    susc_white_matter_sd*randn(n_white_matter_voxels,1) + susc_white_matter_avg;

n_cerebrospinal_fluid_voxels = ...
    length(zubal_phantom_1iso_sd_susceptibility(ismember(zubal_phantom_1iso, cerebrospinal_fluid)));
zubal_phantom_1iso_susceptibility(ismember(zubal_phantom_1iso, cerebrospinal_fluid)) = ...
    susc_cerebrospinal_fluid_avg;
zubal_phantom_1iso_sd_susceptibility(ismember(zubal_phantom_1iso, cerebrospinal_fluid)) = ...
    susc_cerebrospinal_fluid_sd*randn(n_cerebrospinal_fluid_voxels,1) + susc_cerebrospinal_fluid_avg;

n_venous_blood_voxels = ...
    length(zubal_phantom_1iso_sd_susceptibility(ismember(zubal_phantom_1iso, venous_blood)));
zubal_phantom_1iso_susceptibility(ismember(zubal_phantom_1iso, venous_blood)) = ...
    susc_venous_blood_avg;
zubal_phantom_1iso_sd_susceptibility(ismember(zubal_phantom_1iso, venous_blood)) = ...
    susc_venous_blood_sd*randn(n_venous_blood_voxels,1)+susc_venous_blood_avg;

n_caudate_nucleus_voxels = ...
    length(zubal_phantom_1iso_sd_susceptibility(zubal_phantom_1iso == caudate_nucleus));
zubal_phantom_1iso_susceptibility(zubal_phantom_1iso == caudate_nucleus) = ...
    susc_caudate_nucleus_avg;
zubal_phantom_1iso_sd_susceptibility(zubal_phantom_1iso == caudate_nucleus) = ...
    susc_caudate_nucleus_sd*randn(n_caudate_nucleus_voxels,1)+susc_caudate_nucleus_avg;

n_putamen_voxels = ...
    length(zubal_phantom_1iso_sd_susceptibility(zubal_phantom_1iso == putamen));
zubal_phantom_1iso_susceptibility(zubal_phantom_1iso == putamen) = ...
    susc_putamen_avg;
zubal_phantom_1iso_sd_susceptibility(zubal_phantom_1iso == putamen) = ...
    susc_putamen_sd*randn(n_putamen_voxels,1)+susc_putamen_avg;

n_thalamus_voxels = ...
    length(zubal_phantom_1iso_sd_susceptibility(zubal_phantom_1iso == thalamus));
zubal_phantom_1iso_susceptibility(zubal_phantom_1iso == thalamus) = ...
    susc_thalamus_avg;
zubal_phantom_1iso_sd_susceptibility(zubal_phantom_1iso == thalamus) = ...
    susc_thalamus_sd*randn(n_thalamus_voxels,1)+susc_thalamus_avg;

n_globus_pallidus_voxels = ...
    length(zubal_phantom_1iso_sd_susceptibility(zubal_phantom_1iso == globus_pallidus));
zubal_phantom_1iso_susceptibility(zubal_phantom_1iso == globus_pallidus) = ...
    susc_globus_pallidus_avg;
zubal_phantom_1iso_sd_susceptibility(zubal_phantom_1iso == globus_pallidus) = ...
    susc_globus_pallidus_sd*randn(n_globus_pallidus_voxels,1)+susc_globus_pallidus_avg;

n_hippocampus_voxels = ...
    length(zubal_phantom_1iso_sd_susceptibility(zubal_phantom_1iso == hippocampus));
zubal_phantom_1iso_susceptibility(zubal_phantom_1iso == hippocampus) = ...
    susc_hippocampus_avg;
zubal_phantom_1iso_sd_susceptibility(zubal_phantom_1iso == hippocampus) = ...
    susc_hippocampus_sd*randn(n_hippocampus_voxels,1)+susc_hippocampus_avg;

n_amygdala_voxels = ...
    length(zubal_phantom_1iso_sd_susceptibility(zubal_phantom_1iso == amygdala));
zubal_phantom_1iso_susceptibility(zubal_phantom_1iso == amygdala) = ...
    susc_amygdala_avg;
zubal_phantom_1iso_sd_susceptibility(zubal_phantom_1iso == amygdala) = ...
    susc_amygdala_sd*randn(n_amygdala_voxels,1)+susc_amygdala_avg;

n_pituitary_gland_voxels = ...
    length(zubal_phantom_1iso_sd_susceptibility(zubal_phantom_1iso == pituitary_gland));
zubal_phantom_1iso_susceptibility(zubal_phantom_1iso == pituitary_gland) = ...
    susc_pituitary_gland_avg;
zubal_phantom_1iso_sd_susceptibility(zubal_phantom_1iso == pituitary_gland) = ...
    susc_pituitary_gland_sd*randn(n_pituitary_gland_voxels,1)+susc_pituitary_gland_avg;

zubal_phantom_1iso_susceptibility(ismember(zubal_phantom_1iso, non_brain)) = ...
    0;
zubal_phantom_1iso_sd_susceptibility(ismember(zubal_phantom_1iso, non_brain)) = ...
    0.02*randn(length(zubal_phantom_1iso_sd_susceptibility(ismember(zubal_phantom_1iso, non_brain))),1);

zubal_phantom_1iso_susceptibility(ismember(zubal_phantom_1iso, otherBrainRegions)) = ...
    0;
zubal_phantom_1iso_sd_susceptibility(ismember(zubal_phantom_1iso, otherBrainRegions)) = ...
    0.02*randn(length(zubal_phantom_1iso_sd_susceptibility(ismember(zubal_phantom_1iso, otherBrainRegions))),1);

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

zubal_phantom_1iso_susceptibility_2 = ...
    zubal_phantom_1iso_susceptibility - ell;
zubal_phantom_1iso_sd_susceptibility_2 = ...
    zubal_phantom_1iso_sd_susceptibility - ell;

zubal_phantom_1iso_susceptibility_2(zubal_phantom_1iso_susceptibility_2 < -0.08) = ...
    0;
f = (zubal_phantom_1iso_sd_susceptibility_2 < -0.08);
zubal_phantom_1iso_sd_susceptibility_2(zubal_phantom_1iso_sd_susceptibility_2 < -0.08) = ...
    0.02*randn(length(f(f==1)),1);

zubal_phantom_1iso_susceptibility_2(zubal_phantom_1iso_susceptibility_2 == 9.4) = ...
    0;
g = zubal_phantom_1iso_sd_susceptibility_2 == 9.4;
zubal_phantom_1iso_sd_susceptibility_2(zubal_phantom_1iso_sd_susceptibility_2 == 9.4) = ...
    0.02*randn(length(g(g == 1)),1);

zubal_phantom_1iso_susceptibility_2(zubal_phantom_1iso == outside_phantom) = ...
    susc_background_avg; % background [ppm]
zubal_phantom_1iso_sd_susceptibility_2(zubal_phantom_1iso == outside_phantom) = ...
    susc_background_avg; % background [ppm]

zubal_phantom_1iso_susceptibility_2(zubal_phantom_1iso_susceptibility_2 == 8.4) = 9.4;
zubal_phantom_1iso_sd_susceptibility_2(zubal_phantom_1iso_sd_susceptibility_2 == 8.4) = 9.4;

%% Plot phantom

x_coord = round(127*1.5);
y_coord = round(137*1.5);
z_coord = 72*1.5;

figure
[ha, ~] = tight_subplot(2, 2, [0 0], [0 0], [0.01 0.01]);

axes(ha(1))
imagesc(rot90(squeeze(zubal_phantom_1iso_sd_susceptibility_2(:, y_coord, :))), [-0.3 0.3])
axis off, daspect(voxel_size), colormap('gray')

axes(ha(2))
imagesc(rot90(squeeze(zubal_phantom_1iso_sd_susceptibility_2(x_coord, :, :))), [-0.3 0.3])
axis off, daspect(voxel_size), colormap('gray')

axes(ha(3))
imagesc(rot90(squeeze(zubal_phantom_1iso_sd_susceptibility_2(:, :, z_coord))), [-0.3 0.3])
axis off, daspect(voxel_size), colormap('gray')

axes(ha(4))
c = colorbar;
c.Location = 'South';
c.Ticks = [0, 0.5, 1];
c.TickLabels = ({'0', '[a.u.]', '1'});
text(.3, .7, ['X = ' num2str(x_coord) '\newline' 'Y = ' ...
    num2str(y_coord) '\newline' 'Z = ' num2str(z_coord)], 'FontSize', 12)
axis off
title('1-mm iso with constant \chi')
set(gca, 'FontSize', 12)

figure
[ha, pos] = tight_subplot(2, 2, [0 0], [0 0], [0.01 0.01]);

axes(ha(1))
imagesc(rot90(squeeze(zubal_phantom_1iso_susceptibility_2(:, y_coord, :))), [-0.3 0.3])
axis off, daspect(voxel_size), colormap('gray')

axes(ha(2))
imagesc(rot90(squeeze(zubal_phantom_1iso_susceptibility_2(x_coord, :, :))), [-0.3 0.3])
axis off, daspect(voxel_size), colormap('gray')

axes(ha(3))
imagesc(rot90(squeeze(zubal_phantom_1iso_susceptibility_2(:, :, z_coord))), [-0.3 0.3])
axis off, daspect(voxel_size), colormap('gray')

axes(ha(4))
c = colorbar;
c.Location = 'South';
c.Ticks = [0, 0.5, 1];
c.TickLabels = ({'0', '[a.u.]', '1'});
text(.3, .7, ['X = ' num2str(x_coord) '\newline' 'Y = ' ...
    num2str(y_coord) '\newline' 'Z = ' num2str(z_coord)], 'FontSize', 12)
axis off
title('1-mm iso with variable \chi')
set(gca, 'FontSize', 12)

%% Save phantom

nii1 = make_nii(zubal_phantom_1iso_susceptibility_2, [1 1 1], matrix_size_1iso/2+1);
save_nii(nii1, 'zubalPhantom_susceptibility_wholeHead_ellipsoid_1iso.nii')

nii2 = make_nii(zubal_phantom_1iso_sd_susceptibility_2, [1 1 1], matrix_size_1iso/2+1);
save_nii(nii2, 'zubalPhantom_susceptibility_wholeHead_ellipsoid_1iso_sd.nii')
