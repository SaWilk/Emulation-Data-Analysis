% Data Analysis Pursuit-Tracking and Pursuit-Occlusion Paradigm
% Emulation pilot study 2021

% This script contains: 
% * thresholding of source data (ratio)
% * interpolation of atlas on data
% * exclusion of cerebellum
% * run DBSCAN (Ester et al., 1996)

% Adriana Böttcher
% 28.09.2022

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
% input: segmented and preprocessed time domain data 
inputpath       = 'R:\AG-Beste-Studien\Emulation\06_analysis\Emulation-Data-Output\11_DICS\ratio_GAV';
outputpath      = 'R:\AG-Beste-Studien\Emulation\06_analysis\Emulation-Data-Output\12_DBSCAN';
plots_path      = 'R:\AG-Beste-Studien\Emulation\06_analysis\Emulation-Data-Output\12_DBSCAN\plots';

%% configuration

threshold       = 0.01; % theshold: at which % of the power distribution should be thresholded?
min_vox         = 2; % minimum voxels to form a cluster
gridsize        = 0.5; % gridsize used for DICS in cm
epsilon         = 1.5*gridsize; % adjusted 25/10/22: contains neighbouring voxels that are connected edge to edge but not corner to corner
plot            = {'glass', 'MRI'}; % which plots should be generated in the end

freq_labels     = {'alpha', 'beta', 'theta'}; 

identifier = {
    {'const_rand1', 'until_500', -1}; ...
    {'const_rand1', 'until_500', -1}; ...
    {'const_rand1', 'until_500', -1; 'const_rand1', 'after_500', 1}; ...
    }; % ratios and their direction to enter to the clustering/thresholding algorithm; each cell = one freq band

%read atlas
atlas       = ft_read_atlas([fileparts(which('ft_defaults')) '\template\atlas\aal\ROI_MNI_V4.nii']);

%set parameter to use: pow or nai
param       = 'pow';

%% loop: read data, interpolate atlas, threshold, apply dbscan

% read source data
sourcedata = load([inputpath filesep 'source_ratio_all_avg_split_interval']);
sourcedata = sourcedata.source_ratio_all_avg;


for frq = 1:size(freq_labels, 2)
    for contr = 1:size(identifier{frq}, 1)
        data = sourcedata.(freq_labels{frq}).(identifier{frq}{contr, 1}).(identifier{frq}{contr, 2});

        %% interpolate atlas on data

        cfg                 = [];
        cfg.parameter       = 'tissue';
        cfg.interpmethod    = 'nearest';
        atlas_small         = ft_sourceinterpolate(cfg, atlas, data);               % interpolate atlas on data
        
        % MOVE FIELDS FROM INTERPOLATED ATLAS TO DATA
        data.tissue         = atlas_small.tissue;
        data.tissuelabel    = atlas_small.tissuelabel;

        %% threshold data

        %create atlas mask: exclude unlabeled voxels and cerebellum
        atlas_mask  = false(data.dim);
        atlas_mask(data.tissue ~= 0 & data.tissue < 90) = true;                     % excludes non-label and cerebellar structures
        atlas_mask  = reshape(atlas_mask,size(data.(param))); 

        % get contrast direction from identifier
        dir = identifier{frq}{contr, 3};

        switch dir    
            case 1 % positive cluster
                dirtext         = 'pos';                                            % set up string for output files
                vox_pow         = data.(param)(data.(param) > 0 & atlas_mask);      % power in all eligible voxels (that match direction & atlas)
                vox_pow         = sort(vox_pow,1,'descend');                        % sort power values in eligible voxels
                thres_ind       = round(threshold*size(vox_pow,1));            % how much is the threshold of all eligible voxels?
                check_pow(vox_pow,thres_ind)
                thres_pow       = vox_pow(thres_ind);                               % what is the power value in the threshold voxel?
                VOI             = data.(param) > thres_pow;                         % logical mask of thresholded voxels of interest (VOI)
                
            case -1 % NEGATIVE CLUSTERS
                dirtext         = 'neg';                                            % see above
                vox_pow         = data.(param)(data.(param) < 0 & atlas_mask);
                vox_pow         = sort(vox_pow,1,'ascend');
                thres_ind       = round(threshold*size(vox_pow,1));
                check_pow(vox_pow,thres_ind)
                thres_pow       = vox_pow(thres_ind);
                VOI             = data.(param) < thres_pow;
        end % dir switch

        VOI        = VOI & atlas_mask;
        fprintf('Selected %g voxels (direction: %g)\n',sum(VOI),dir);
        
        %% RUN DBSCAN
        % dbscan(data, epsilon, minpnts)
        % Epsilon is the neighbourhood search radius. This should be set to the grid
        % size of sourcemodel/leadfield to make sure that voxels are neighbouring on
        % an edge, not just on a corner, i.e.   |v||v|  , and not |v|
        %                                          |v|               |v|
        
        clusters        = dbscan(data.pos(VOI,:), epsilon, min_vox);
        
        % WRITE TO ORIGINAL DATA FILE
        data.clusterslabelmat       = zeros(data.dim);
        data.clusterslabelmat(VOI)  = clusters;
        
        num_clusters    = max(unique(data.clusterslabelmat(~isnan(data.clusterslabelmat)))); % number of clusters is maximum in labelmat

        %% PRINT CLUSTERS TO .TXT FILE (NUMBER OF VOXELS PER LABEL)
        
        % DEFINE FILE NAMING
        filename    = [freq_labels{frq} '_' identifier{frq}{contr, 1} '_' identifier{frq}{contr, 2} '_' dirtext];
        
%         for clust = 1:num_clusters     % loop through all clusters
%             
%             data.clustlabels{clust} = data.tissuelabel(data.tissue(data.clusterslabelmat == clust));
%             
%             write2txt([outputpath filesep filename '.txt'],['Cluster ',num2str(clust),':']);
%             
%             uniques = unique(data.clustlabels{clust});
%             to_sort = [];
%             
%             for i = 1:numel(uniques)
%                 to_sort = [to_sort; i numel(data.clustlabels{clust}(strcmp(data.clustlabels{clust},uniques{i})))];
%             end
%             
%             sorted = sortrows(to_sort, 2, 'descend');
%             
%             for i = 1:size(to_sort,1)
%                 write2txt([outputpath filesep filename '.txt'],[num2str(sorted(i,2),'%03.f'),': ', num2str(uniques{sorted(i,1)})]);
%             end
%             
%         end % cluster loop

        save([outputpath filesep filename '.mat'], 'data');

       %% PLOTTING
        
        % PLOT CLUSTERS ON ORTHOGONAL MRI SLICES
%         if sum(contains(plot,'MRI','IgnoreCase',true)) > 0
%             
%             % LOAD MRI AND ATLAS
%             load standard_mri.mat;
%             
%             mri     = ft_convert_units(mri,'cm');
%             mri     = ft_volumereslice([],mri);
%             
%             atlas   = ft_read_atlas([fileparts(which('ft_defaults')) '\template\atlas\aal\ROI_MNI_V4.nii']);
%             atlas   = ft_convert_units(atlas,'cm');
%             
%             % INTERPOLATE DATA ON MRI
%             cfg                     = [];
%             cfg.parameter           = 'all';
%             cfg.interpmethod        = 'nearest';
%             data_int                = ft_sourceinterpolate(cfg,data,mri);
%             
%             % CREATE MASK (WHICH DOES NOT WORK AS INTENDED)
%             data_int.mask           = false(data_int.dim);
%             data_int.mask(data_int.clusterslabelmat > 0) = true;
%             
%             % CREATE ORTHO-PLOT
%             cfg                     = [];
%             cfg.method              = 'ortho';
%             cfg.funcolorlim         = 'maxabs';
%             cfg.funparameter        = 'clusterslabelmat';
%             cfg.maskparameter       = 'mask';
%             cfg.maskstyle           = 'opacity';
%             cfg.funcolormap         = 'bipolar';
%             cfg.camlight            = 'off';
%             cfg.atlas               = atlas;
%             cfg.title               = filename;   % not working
%             ft_sourceplot(cfg, data_int, mri);
%             
%             title(gca, filename,'Interpreter','none')
%             
%             % SAVE AS INTERACTIVE FIGURE
%             savefig([plots_path filesep filename])
%         end %mri plot
%         
%         % GLASS BRAIN PLOTTING ON DBSCAN OUTPUT
%         if sum(contains(plot,'glass','IgnoreCase',true)) > 0
%             
%             % LOAD MRI AND ATLAS
%             load standard_mri.mat;
%             
%             atlas   = ft_read_atlas([fileparts(which('ft_defaults')) '\template\atlas\aal\ROI_MNI_V4.nii']);
%             atlas   = ft_convert_units(atlas,'cm');
%             
%             % INTERPOLATE DATA ON MRI
%             cfg                     = [];
%             cfg.parameter           = 'all';
%             cfg.interpmethod        = 'nearest';
%             data_int                = ft_sourceinterpolate(cfg,data,mri);
% 
%             % load discrete colormap
%             load('cb_discrete.mat');
%             
%             if num_clusters > 12
%                 cb = cat(1,cb,cb);
%             end
%             
%             resFig = figure;
%             
%             figure;
%             fieldtrip_mesh = load('surface_pial_both.mat');
%             ft_plot_mesh(fieldtrip_mesh.mesh, 'facecolor', 'cortex', 'facealpha', 0.2);
%             % PW standard facecolor: [0.78, 0.76, 0.66]
%             hold on;
%             
%             for c = 1:num_clusters
%                 pos{c}  = data_int.pos(data_int.clusterslabelmat == c,:);
%                 xs      = pos{c}(:, 1);
%                 ys      = pos{c}(:, 2);
%                 zs      = pos{c}(:, 3);
%                 
%                 scatter3(xs, ys, zs, 5,cb(c,:),'filled');
%                 hold on;
%             end
%             
%             label = strcat('Cluster ',cellstr(num2str((1:c)'))');
%             
%             % CREATE DIFFERENT VIEWS AND SAVE
%             view_angle = {[-180,0],[90,0],[0,0],[-90,0],[0,90]};
%             view_angle_labels = {'front', 'left', 'back', 'right', 'dorsal'};
%             
%             for ind = 1:numel(view_angle)              
%                 pw_sourceSubplot(resFig,2,3,ind)
%                 view(view_angle{ind});
%                 title(view_angle_labels{ind});
%             end
%             
%             figure(resFig)
%             
%             leg = legend('brain',label{:});
%             set(leg,'Position',[0.75,0.2,0.1,0.2]);
%             set(gcf, 'color', 'white');
%             set(gca, 'color', 'white');
%             
%             % SAVE AS GRAPHIC
%             exportgraphics(gcf,[plots_path filesep filename '.png'],'Resolution',1000,'BackgroundColor','white');           
%         end %glass brain

    end % contrast loop
end % freq loop
    

%% sub functions

function check_pow(data,thres) % check whether variable has values (to indicate potentially wrong dir)
    if isempty(data)
        error('Less eligble voxels than values above/below threshold! Selected direction may be wrong!')
    elseif numel(data) < thres
        error('Less eligble voxels than values above/below threshold! Selected direction may be wrong!')
    end
end
%%
function write2txt(file,text) % check whether file exists or create a new one, then add text in a new line
    if exist(file,'file')
        fid = fopen(file,'at');
        fprintf(fid,'%s\n%s\n',text);
        fclose(fid);
    else
        dlmwrite(file,text,'delimiter','','newline','pc');
    end
end