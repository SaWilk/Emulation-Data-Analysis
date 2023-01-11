% Data Analysis Pursuit-Tracking and Pursuit-Occlusion Paradigm
% Emulation pilot study 2021

% This script contains: 
% MRI plots for DBSCAN clusters
% this scipt: for task A, B

% plots for manuscript!

% Adriana Böttcher
% 07.11.2022

%%
clear
clc

%% path definition

% load fieldtrip toolbox
addpath('R:\AG-Beste-Orga\Skripts\Toolbox\fieldtrip-20210212\');
ft_defaults;

% load custom ft toolbox
addpath('R:\AG-Beste-Orga\Skripts\Fieldtrip\ft_mmcustom');

% load pauls source subplot function
addpath('R:\AG-Beste-Orga\Skripts\PW code\misc');

% initialize input and output folder
inputpath_DICS_GAV  = 'R:\AG-Beste-Studien\Emulation\06_analysis\Emulation-Data-Output\11_DICS\ratio_GAV';
inputpath_DBSCAN    = 'R:\AG-Beste-Studien\Emulation\06_analysis\Emulation-Data-Output\12_DBSCAN';
outputpath_plots    = 'R:\AG-Beste-Studien\Emulation\06_analysis\Emulation-Data-Output\13_figures\dbscan';

% electrode coordinates & layout
load fieldtrip_EEG_KJP_elec_61.mat;                                       
elec        = ft_convert_units(elec,'cm');
elec.label  = upper(elec.label);                                            

% load neighbours
load fieldtrip_EEG_KJP_neighbours_61.mat;                                  

% load headmodel
load standard_bem;
headmodel   = ft_convert_units(vol,'cm');

% load mri
load standard_mri.mat
mri         = ft_convert_units(mri,'cm');

% load DICS output (ratio GAV)
source = load([inputpath_DICS_GAV filesep 'source_ratio_all_avg_split_interval.mat']);
source = source.source_ratio_all_avg;

% same for B
% load DICS output (ratio GAV)
source_B = load([inputpath_DICS_GAV filesep 'source_ratio_B_all_avg.mat']);
%% define clusters (anatomical regions)

conditions  = {'until_500_alpha', 'until_500_beta', 'until_500_theta', 'after_500_theta'};
conditions_B = {'theta_occl_nonoccl'};

files       = {'alpha_const_rand1_until_500_neg', ...
    'beta_const_rand1_until_500_neg', ...
    'theta_const_rand1_until_500_neg', ...
    'theta_const_rand1_after_500_pos'};
files_B = {'theta_occl_nonoccl_neg'};

cluster     = {};
cluster_B   = {};

source_data     = {};
source_data_B   = {};

% alpha
source_data{1} = source.alpha.const_rand1.until_500;
cluster{1} = {{1, 'Precuneus_L', 'Parietal_Sup_L', 'Cingulum_Mid_L', 'Occipital_Mid_L', 'Precuneus_R'}};

% beta
source_data{2} = source.beta.const_rand1.until_500;
cluster{2} = {{1, 'Postcentral_L', 'Parietal_Sup_L'}, ...
    {2, 'Frontal_Sup_L', 'Supp_Motor_Area_L', 'Paracentral_Lobule_L', 'Precentral_L'}};

% theta, until 500
source_data{3} = source.theta.const_rand1.until_500;
cluster{3} = {{1, 'Precuneus_L', 'Precuneus_R', 'Cingulum_Mid_L'}};

% theta, after 500
% todo adjust, clusters are wrong
source_data{4} = source.theta.const_rand1.after_500;
cluster{4} = {{1, 'Precentral_R', 'Postcentral_R'}, ...
    {2, 'Precentral_L', 'Postcentral_L'}, ...
    {3, 'Supp_Motor_Area_R'}};

% task B, theta, whole interval
source_data_B{1} = source_B.source_ratio_all_avg.theta.occl_nonoccl;
cluster_B{1} = {{1, 'Fusiform_R', 'Occipital_Inf_R'}, ...
    {2, 'Occipital_Inf_L', 'Temporal_Inf_L'}, ...
    {3, 'Temporal_Mid_R', 'Occipital_Mid_R', 'Occipital_Inf_R'}};

%%  interpolate source data on MRI and save

source_int = {};
source_int_B = {};

cfg             = [];
cfg.parameter   = 'pow';

for cond = 1:size(conditions, 2)
    source_int{cond} = ft_sourceinterpolate(cfg, source_data{cond}, mri);
end

source_int_B{1} = ft_sourceinterpolate(cfg, source_data_B{1}, mri);

% extract maximum power for plotting with appropriate colorbar
maxPow  = max(cellfun(@(x) max(max(abs(x.pow))), source_int));
% maxPow for B is smaller than maxPow, so this value is used also for B

%% task A
for cond = 1:size(conditions, 2)
    % load dbscan results for this condition
    sourcemodel = load([inputpath_DBSCAN filesep files{cond}]);
    sourcemodel = sourcemodel.data;

    for clust = 1:size(cluster{cond}, 2)

        idx = find(ismember(sourcemodel.tissuelabel,cluster{cond}{clust}(2:end)));

        % empty mask in size of dim
        sourcemodel.mask = [];
        sourcemodel.mask = false(sourcemodel.dim);       

        for region = 1:size(idx,2) % if cluster & label match, set mask to true
           if isa(cluster{cond}{clust}{1},'cell')
               for ind = 1:size(cluster{cond}{clust}{1},2)
                   sourcemodel.mask(sourcemodel.tissue == idx(region) & sourcemodel.clusterslabelmat == cluster{cond}{clust}{1}{ind}) = true;
               end
           else
               sourcemodel.mask(sourcemodel.tissue == idx(region) & sourcemodel.clusterslabelmat == cluster{cond}{clust}{1}) = true;
           end
        end % region loop

        % interpolate on mri
        cfg                 = [];
        cfg.parameter       = 'all';
        cfg.interpmethod    = 'nearest';
        sourcemodel_int     = ft_sourceinterpolate(cfg, sourcemodel, mri);
        
        
        temp        = source_int{cond};
        temp.mask   = sourcemodel_int.mask > 0;
        temp.pow(~temp.mask) = 0;
        temp.inside = [];
        temp.inside = false(mri.dim);
        temp.inside(temp.mask) = true;

        % plot

        cfg                     = [];
        cfg.method              = 'ortho';
        cfg.funparameter        = 'pow';
        cfg.funcolormap         = 'turbo'; %'viridis';
        cfg.funcolorlim         = [-maxPow-0.1 maxPow+0.1];
        cfg.maskparameter       = 'inside';
        % cfg.maskstyle           = 'colormix';
        % cfg.opacitystyle        = 'vdown'; % 'rampup'
        % cfg.opacitylim          = [-1 1];% [0 0.1]
        cfg.camlight            = 'off';
        cfg.crosshair           = 'no';
        % cfg.atlas               = atlas;
        % cfg.renderer            = 'opengl';
        ft_sourceplot(cfg, temp, mri);

        set(gcf,'color','k')
        
        hc                  = findall(gcf,'type','ColorBar');     
        hc.Color            = 'w';
        hc.FontSize         = 10;
        hc.Label.FontSize   = 12;

    end % cluster loop
end % condition loop

%% task b
for cond = 1:size(conditions_B, 2)
    % load dbscan results for this condition
    sourcemodel = load([inputpath_DBSCAN filesep files_B{cond}]);
    sourcemodel = sourcemodel.data;

    for clust = 1:size(cluster_B{cond}, 2)

        idx = find(ismember(sourcemodel.tissuelabel,cluster_B{cond}{clust}(2:end)));

        % empty mask in size of dim
        sourcemodel.mask = [];
        sourcemodel.mask = false(sourcemodel.dim);       

        for region = 1:size(idx,2) % if cluster & label match, set mask to true
           if isa(cluster_B{cond}{clust}{1},'cell')
               for ind = 1:size(cluster_B{cond}{clust}{1},2)
                   sourcemodel.mask(sourcemodel.tissue == idx(region) & sourcemodel.clusterslabelmat == cluster_B{cond}{clust}{1}{ind}) = true;
               end
           else
               sourcemodel.mask(sourcemodel.tissue == idx(region) & sourcemodel.clusterslabelmat == cluster_B{cond}{clust}{1}) = true;
           end
        end % region loop

        % interpolate on mri
        cfg                 = [];
        cfg.parameter       = 'all';
        cfg.interpmethod    = 'nearest';
        sourcemodel_int     = ft_sourceinterpolate(cfg, sourcemodel, mri);
        
        
        temp        = source_int{cond};
        temp.mask   = sourcemodel_int.mask > 0;
        temp.pow(~temp.mask) = 0;
        temp.inside = [];
        temp.inside = false(mri.dim);
        temp.inside(temp.mask) = true;

        % plot

        cfg                     = [];
        cfg.method              = 'ortho';
        cfg.funparameter        = 'pow';
        cfg.funcolormap         = 'turbo'; %'viridis';
        cfg.funcolorlim         = [-maxPow-0.1 maxPow+0.1];
        cfg.maskparameter       = 'inside';
        % cfg.maskstyle           = 'colormix';
        % cfg.opacitystyle        = 'vdown'; % 'rampup'
        % cfg.opacitylim          = [-1 1];% [0 0.1]
        cfg.camlight            = 'off';
        cfg.crosshair           = 'no';
        % cfg.atlas               = atlas;
        % cfg.renderer            = 'opengl';
        ft_sourceplot(cfg, temp, mri);

        set(gcf,'color','k')
        
        hc                  = findall(gcf,'type','ColorBar');     
        hc.Color            = 'w';
        hc.FontSize         = 10;
        hc.Label.FontSize   = 12;

    end % cluster loop
end % condition loop
