% The neurophysiology of continuous action monitoring
% Saskia Wilken, Adriana Böttcher, Nico Adelhöfer, Markus Raab, Sven
% Hoffmann, Christian Beste

% thresholding data, interpolation on atlas, apply DBSCAN (Ester et al., 1996)
% for experiment 1, 1 contrast & 3 frequency bands

% created by:
% Adriana Boettcher, Cognitive Neurophysiology TU Dresden
% 2022

% based on code by Paul Wendiggensen

%%
clear
clc

%% path definition

% load fieldtrip toolbox

% initialize input and output folder

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
        
        % move fields from interpolated atlas to data
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
                thres_ind       = round(threshold*size(vox_pow,1));           
                check_pow(vox_pow,thres_ind)
                thres_pow       = vox_pow(thres_ind);                               
                VOI             = data.(param) > thres_pow;                         % logical mask of thresholded voxels of interest (VOI)
                
            case -1 % negative clusters
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
        
        %% run dbscan
        
        clusters        = dbscan(data.pos(VOI,:), epsilon, min_vox);
        
        % write to original data file
        data.clusterslabelmat       = zeros(data.dim);
        data.clusterslabelmat(VOI)  = clusters;
        
        num_clusters    = max(unique(data.clusterslabelmat(~isnan(data.clusterslabelmat)))); % number of clusters is maximum in labelmat

        %% print clusters to txt file
        
        % define file naming
        filename    = [freq_labels{frq} '_' identifier{frq}{contr, 1} '_' identifier{frq}{contr, 2} '_' dirtext];
        
        for clust = 1:num_clusters     % loop through all clusters
            
            data.clustlabels{clust} = data.tissuelabel(data.tissue(data.clusterslabelmat == clust));
            
            write2txt([outputpath filesep filename '.txt'],['Cluster ',num2str(clust),':']);
            
            uniques = unique(data.clustlabels{clust});
            to_sort = [];
            
            for i = 1:numel(uniques)
                to_sort = [to_sort; i numel(data.clustlabels{clust}(strcmp(data.clustlabels{clust},uniques{i})))];
            end
            
            sorted = sortrows(to_sort, 2, 'descend');
            
            for i = 1:size(to_sort,1)
                write2txt([outputpath filesep filename '.txt'],[num2str(sorted(i,2),'%03.f'),': ', num2str(uniques{sorted(i,1)})]);
            end
            
        end % cluster loop

        save([outputpath filesep filename '.mat'], 'data');
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