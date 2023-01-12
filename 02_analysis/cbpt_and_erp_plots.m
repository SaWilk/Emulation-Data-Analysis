% ERP analyses

% creates plots and performs further analyses to investigate
% ERP

% author: Saskia Wilken
% creation date: 21.09.2022 

%% empty everything

format compact
format long G
clear
clc
close all

%% Define Paths

% file location
file_path = matlab.desktop.editor.getActiveFilename;
file_path = fileparts(file_path);
addpath(file_path);

% get data paths for parent dirs
filepath_parts = strsplit(file_path, filesep);
parent_dir_2 = filepath_parts(1:end-2);
parent_dir_2 = strjoin(parent_dir_2, filesep);
output_dir = strjoin([parent_dir_2, "Emulation-Data-Output"], filesep);
mean_matrices_path = strjoin([output_dir, 'mean_matrices'], filesep);
mkdir(mean_matrices_path);
mean_matrices_peaks_epoched_path = strjoin([mean_matrices_path, 'peaks'], filesep);
mkdir(mean_matrices_peaks_epoched_path);
% output
csd_dir = strjoin([output_dir, "csd_transform"], filesep);
% cluster based permutation test
neighbors_dir = strjoin([parent_dir_2, "Emulation-Data-Input", "EEG_files"], filesep);


%% Load CDS data

eeglab;

% long epochs
cd(csd_dir);
%list all *.set files in inputpath
file_names = dir('*.set');
%concatenate into one cell array
files2read = {file_names.name};
files2read = files2read(~contains(files2read, 'KMY6K')); % remove subject who 
% did not properly track during task B
% load eeg data epoched around peaks
for idx = 1:length(files2read)
    ALLEEG2(idx) = pop_loadset('filename',files2read{idx});
end

ALLEEG = ALLEEG2;
clear ALLEEG2


%% Prepare for Plotting

cd(mean_matrices_peaks_epoched_path)

load('z_mean_struct_csd.mat')
load('mean_struct_csd.mat')
load("condition_indices_csd.mat")
mean_struct = renameStructField(mean_struct, 'diff_occ', 'diff_occ_vis');
z_mean_struct = renameStructField(z_mean_struct, 'diff_occ', 'diff_occ_vis');
% remove everything reminding of random2 since we are not analyzing it.
rm_fields = {'rand2', 'diff_rand1_rand2', 'diff_const_rand2', 'num_epochs_rand2'};

for field = 1:length(rm_fields)
    if isfield(z_mean_struct, rm_fields{field})
        z_mean_struct = rmfield(z_mean_struct, rm_fields{field});
    end
    if isfield(mean_struct, rm_fields{field})
        mean_struct = rmfield(mean_struct, rm_fields{field});
    end
end

cond_labs = {'all', 'occ', 'vis', 'rand1', 'const', 'diff_occ_vis', ...
    'diff_const_rand1'};

avg_peak_dist = 300;
% average across subjects
for lab = 1:length(cond_labs)
    epoch_mean.(cond_labs{lab}) = mean(mean_struct.(cond_labs{lab}),3);
end

% specify parameters for (topo)plotting
base_dur = 500;
num_topos = 16;
topo_rows = 4;
topo_cols = 4;
topo_dur = 30;
epoch_end = 480 - topo_dur;
latencies_onset = linspace(0, epoch_end, num_topos);
latencies_offset = latencies_onset + topo_dur;


%% Prepare ERP Plotting 

cd(mean_matrices_peaks_epoched_path)

base_dur = 500;

title_string_z = strcat(...
    ['z-normalized ERP, CSD transformed of all subjects averaged across epochs around peaks']);
subtitle_string = strcat(['all channels']);


%% Cluster Based Permutation Test Attempt with Fieldtrip

% Script by Adriana Böttcher
% Adapted by Saskia Wilken
% which of these clusters are significant? 


%% add fieldtrip path

% fieldtrip path
% Add fieldtrip path here

ft_defaults;


%% Prepare Cluster Based Permutaiton Testing

cd(mean_matrices_peaks_epoched_path)
load("condition_indices_csd.mat")

clear all_eeg_field all_eeg_cell TMPEEG
% All data, compare with zero
% convert to fieldtrip:
for lab = find(~contains(cond_labs, 'diff'))
    for s = 1:length(ALLEEG)
        TMPEEG = ALLEEG(s);
        TMPEEG.data = ALLEEG(s).CSD_data;
        % create eeglab temporary datasets that only contain condition
        % epochs
        TMPEEG.data = TMPEEG.data(:,:,cond_ind(s).(cond_labs{lab}));
        all_eeg_field(s).(cond_labs{lab}) = eeglab2fieldtrip(TMPEEG, 'timelockanalysis', 'none');
        all_eeg_field(s).(cond_labs{lab}).elec.label = upper(all_eeg_field(s).(cond_labs{lab}).elec.label);
        all_eeg_field(s).(cond_labs{lab}).label = upper(all_eeg_field(s).(cond_labs{lab}).label);
        all_eeg_cell.(cond_labs{lab}){s} = all_eeg_field(s).(cond_labs{lab});
    end
end

% creating a null condition dataset
for s = 1:length(all_eeg_field)
    all_eeg_field(s).('zero') = all_eeg_field.all;
    all_eeg_field(s).zero.avg = zeros(size(all_eeg_field(s).zero.avg));
    all_eeg_field(s).zero.var = zeros(size(all_eeg_field(s).zero.avg));
    all_eeg_field(s).zero.elec.label = upper(all_eeg_field(s).zero.elec.label);
    all_eeg_field(s).zero.label = upper(all_eeg_field(s).zero.label);
    all_eeg_cell.zero{s} = all_eeg_field(s).zero;
end


%% specify design matrix

design = zeros(2, size(all_eeg_cell.all, 2)*2);
for i = 1:size(all_eeg_cell.all, 2)
    design(1,i) = i;
end
for i = 1:size(all_eeg_cell.all, 2)
    design(1,size(all_eeg_cell.all, 2)+i) = i;
end
design(2,1:size(all_eeg_cell.all, 2))        = 1;
design(2,size(all_eeg_cell.all, 2)+1:2*size(all_eeg_cell.all, 2)) = 2;


%% Get Neighbours

% generates a file called 'neighbours'
load(strjoin([neighbors_dir, 'fieldtrip_EEG_KJP_neighbours_61.mat'], filesep))

% adjust neighbours file
idx = find(strcmp({neighbours.label}, 'P11'));
neighbours(idx).neighblabel % even though it looks strange on the topography, 
% % the neighbors are likely correct

cfg_neighb        = [];
cfg_neighb.method = 'distance';
neighbours        = neighbours; % ft_prepare_neighbours(cfg_neighb, dataFC_LP);

% CLUSTER TEST
calpha{1}  = 0.001;
alpha  = 0.05;

% cfg is the configuraiton structure of fieldtrip
cfg                     = [];
cfg.design              = design;
cfg.uvar                = 1;
cfg.ivar                = 2;
cfg.channel             = {'all'};
cfg.avgovertime         = 'yes';
cfg.method              = 'montecarlo';
cfg.statistic           = 'depsamplesT'; % indepsamples With comparison against 0
cfg.correctm            = 'cluster';
cfg.clusteralpha        = calpha{1};  
cfg.clusterstatistic    = 'maxsum';
cfg.minnbchan           = 2;
cfg.neighbours          = neighbours;
cfg.tail                = 0;
cfg.clustertail         = 0;
cfg.alpha               = alpha;        
cfg.numrandomization    = 3000;
cfg.latency             = [0, 0.75];% time range to perform the CBP-Test on 

clear stats


%% For Reporting CBPT Results, do CBPT test without averaging for 
% the whole epochs

cfg.avgovertime         = 'no';
cfg.clusteralpha        = 0.01;
cfg.latency             = [0, 0.75];

cond_labs_short = {"occ", "vis", "const", "rand1"};
clear no_avg_cbpt

for cond = 1:2:length(cond_labs_short)-1
    no_avg_cbpt.(strjoin(["diff_", cond_labs_short{cond}, "_", ...
        cond_labs_short{cond+1}], "")) = ft_timelockstatistics(cfg, ...
        all_eeg_cell.(cond_labs_short{cond}){:}, ...
        all_eeg_cell.(cond_labs_short{cond+1}){:});
end


%% Prepare for Topoplotting

clear topo_struct

epoch_mean = renameStructField(epoch_mean, 'diff_all_zero', 'all');
zlims_all = [-0.0225, 0.0225]; % for microvolts/mm²


%% Prepare ERP plot of electrodes that are different in between conditions, 
% averaged together for each cluster

FONTSIZE = 9;
LINEWIDTH = 1;
PLOTSIZE = [0 , 0, 10, 7];

% Plot the ERPs and the time-points of significant differences
% Parameters: 
chan_locs = ALLEEG(1).chanlocs;
chan_lab = {chan_locs.labels};
plot_area =  find(z_mean_struct.time_vec == 0):length(z_mean_struct.time_vec);

% constant - random, occluded - visible

% crop the baseline from the data to make indexing consistent within this
% bit
cond_field_names = fieldnames(epoch_mean);

for i = 1:length(cond_field_names)
    clear cur_data
    cur_data = epoch_mean.(cond_field_names{i});
    epoch_mean_no_base.(cond_field_names{i}) = cur_data(:, plot_area);
    clear cur_data
    cur_data = epoch_mean.(cond_field_names{i});
    epoch_mean_no_base.(cond_field_names{i}) = cur_data(:, plot_area);
end

% Define Parameters
cond_names = {'occluded','visible'; 'constant', 'random1'};
cond_labs_mat = reshape(cond_labs_short, 2, numel(cond_labs_short)/2);
contrasts = fieldnames(no_avg_cbpt);
cluster_orient = {"pos", "neg"};
SAMP_DIST = 4 ; % distance between samples in ms


%% Create Plot of both ERP as well as topography of current clusters 

for cont = 1:length(fieldnames(no_avg_cbpt))
    cur_contrast = no_avg_cbpt.(contrasts{cont});
    field_names = fieldnames(cur_contrast);
    clust_no_idx = find(endsWith(field_names, 'clusters'));
    clust_mat_idx = find(endsWith(field_names, 'clusterslabelmat'));

    for orient = 1:length(clust_no_idx)

        for clust_no = 1:3

            cluster_summary = cur_contrast.(field_names{clust_no_idx(orient)})(clust_no); % negitive cluster with largest t-value sum descriptive
            cluster_elecxsamp = cur_contrast.(field_names{clust_mat_idx(orient)}) == clust_no;
            cluster_start_idx = min(find(any(cluster_elecxsamp,1))); % find start of anything being part of the cluster
            cluster_start_ms = cluster_start_idx*SAMP_DIST;
            cluster_end_idx = max(find(any(cluster_elecxsamp,1))); % ... and the end
            cluster_end_ms = cluster_end_idx*SAMP_DIST;

            % From when to when to plot?
            topo_bins_idx = [cluster_start_idx:cluster_end_idx];
            topo_bins_ms = [cluster_start_ms, cluster_end_ms];

            % Prepare data for topoplot
            clear clust_plot_struct clust_chans
            clear tmp current_samps
            % Put the mean of the time inteval into
            clust_plot_struct = mean(epoch_mean_no_base.(contrasts{cont})(:, topo_bins_idx),2);
            % get idx of electrodes that are part of the cluster in the bin
            [tmp, ~] = find(cluster_elecxsamp(:, topo_bins_idx));
            clust_chans= unique(tmp);

            % average the electrodes of one cluster
            avg_chans_plot_struct.(cond_labs_mat{1, cont}) = mean( ...
                epoch_mean.(cond_labs_mat{1, cont})(clust_chans, :),1);
            avg_chans_plot_struct.(cond_labs_mat{2, cont}) = mean( ...
                epoch_mean.(cond_labs_mat{2, cont})(clust_chans, :),1);
            x_axis_ticks = unique(sort([min(mean_struct.time_vec(plot_area)):250:0, 0, ...
                ceil(max(mean_struct.time_vec(plot_area))/10)*10:-250:0]));

            % actual plot code
            figure()
            set(gcf,'color','w') % set background color to white
            title_string_z = strjoin([...
                "Averaged ERP of Electrodes ", strjoin(chan_lab(clust_chans), ", "), ...
                " in time window ", num2str(topo_bins_ms(1)), " ms to ", num2str(topo_bins_ms(2)), " ms, " ...
                clust_no, "th largest ", cluster_orient{orient}, " cluster"],"");
            min_idx = 1;
            plot_area = min_idx:length(mean_struct.time_vec);
            hold on
            title(strjoin(["t-value sum = ", num2str(cluster_summary.clusterstat)],""))
            plot(mean_struct.time_vec(plot_area), ...
                avg_chans_plot_struct.(cond_labs_mat{1, cont})(plot_area), 'LineWidth', LINEWIDTH*1.5)
            plot(mean_struct.time_vec(plot_area), ...
                avg_chans_plot_struct.(cond_labs_mat{2, cont})(plot_area), 'LineWidth', LINEWIDTH*1.5)
            xline(0, 'LineWidth', LINEWIDTH*1.5)
            xline(topo_bins_ms(1), 'LineWidth', 1, 'Color', 'k', 'Linestyle', '--','LineWidth', LINEWIDTH)
            xline(topo_bins_ms(2), 'LineWidth', 1, 'Color', 'k', 'Linestyle', '--','LineWidth', LINEWIDTH)
            yline(0, 'LineWidth', LINEWIDTH);
            hold off
            legend({cond_names{cont, :}}, 'location', 'northwest')
            set(gca,'fontsize', FONTSIZE*1.3);
            xlabel("time (ms)")
            xlim([min(x_axis_ticks)-5, max(x_axis_ticks)+5]);
            xticks(x_axis_ticks)
            box off
            ylabel("amplitude (µV/mm²)")

            set(gca,'linewidth',LINEWIDTH)
            f = gcf;
            % copygraphics(gcf, 'BackgroundColor', 'none', 'ContentType', 'vector');
            % copies the graph to the clipboard
            sgtitle(title_string_z, 'Interpreter', 'None')
            set(gcf,'PaperUnits','centimeters','PaperPosition',PLOTSIZE)
            filename = strjoin(["Event-Related_Potential_of", cond_names{cont, 1}, ...
                "_-_", cond_names{cont, 2}, "_", cluster_orient{orient}, "_Clust_", ...
                clust_no],"")
            print('-djpeg', strjoin([filename, ".jpg"],""), '-r600');
            saveas(gcf,strjoin([filename, ".svg"],""),'svg')
               
            % topoplot
            figure()
            topoplot(clust_plot_struct, chan_locs, ...
                'emarker2', {clust_chans, 'x', 'k', 16, 1}, ...
                'maplimits', zlims_all, ...
                'style', 'map', ...
                'electrodes', 'off');
            a = colorbar;
            position = get(a,'Position');
            ylabel(a,"amplitude (µV/mm²)",'FontSize',26,'Rotation',270, 'Position', [position(1)+3 position(2)]);
            colormap 'parula'
            set(gca,'fontsize', FONTSIZE*5);
            % This copies the plot to the clipboard 
            copygraphics(gcf, 'BackgroundColor', 'none', 'ContentType', 'vector');
            pause % to copy the graphics into inkscape
            % can then insert into inkscape using Crtl + V
            % figure
            sgtitle(title_string_z, 'Interpreter', 'None')
            set(gcf,'PaperUnits','centimeters','PaperPosition',PLOTSIZE*0.5)
            filename = strjoin(["Average_of_Difference_of_", cond_names{cont, 1}, ...
                "_-_", cond_names{cont, 2}, "_", cluster_orient{orient}, "_Clust_", ...
                clust_no],"")
            print('-djpeg', strjoin([filename, ".jpg"],""), '-r600');
            saveas(gcf,strjoin([filename, ".svg"],""),'svg')
        end % cluster number
    end % cluster orientation (negative, positive)
end % contrast


%% Make Scree Plot of Cluster t-values

FONTSIZE = 26;
LINEWIDTH = 2;
PLOTSIZE = [0, 0, 5, 5];
CRITERION = [3, 3];
% Scree plot of cluster t-values
for cont = 1:length(fieldnames(no_avg_cbpt))
    for orient = 1:length(clust_no_idx)
        cur_data = [no_avg_cbpt.(contrasts{cont}).(strjoin([cluster_orient{orient}, "clusters"],"")).clusterstat];
        figure()
        plot(abs(cur_data), 'LineWidth', LINEWIDTH*1.5)
        title(strjoin(["Scree Plot of T-Value Sums of All ", cluster_orient{orient}, " Clusters"],""))
        subtitle(strjoin(["Contrast:", cond_names{cont, 1} "-", cond_names{cont, 2}]))
        xlabel("cluster number")
        ylabel("cluster magnitude (t-value sum of cluster")
        xticks(unique([1:3, round(linspace(4, length(cur_data), 5))]))
        xline(CRITERION(cont), 'Linestyle', "--", "LineWidth", LINEWIDTH)
        box off
        set(gca,'fontsize', FONTSIZE*0.75, 'linewidth',LINEWIDTH);
        set(gcf,'color','w','PaperUnits','inches','PaperPosition',PLOTSIZE)
        filename = strjoin(["Scree plot of T-Value Sums of all ", cluster_orient{orient}, ...
            ", ", cond_names{cont, 1} "-", cond_names{cont, 2}]);
        print('-djpeg', strjoin([filename, ".jpg"]), '-r600');
        print('-dsvg', strjoin([filename, ".svg"]), '-vector');
    end
end


%% Export Clusters to report

% Extract the cluster locations and labels from the clusters we want to
% investigate

% const rand
clust_no = 3;
cluster_elecxsamp = no_avg_cbpt.diff_const_rand1.negclusterslabelmat == clust_no;
tmp = find(any(cluster_elecxsamp, 2))'; % find the electrodes inside the cluster
clust{1} = chan_lab(tmp);

clust_no = 1;
cluster_elecxsamp = no_avg_cbpt.diff_const_rand1.negclusterslabelmat == clust_no;
tmp = find(any(cluster_elecxsamp, 2))'; % find the electrodes inside the cluster
clust{2} = chan_lab(tmp);

% occ vis
clust_no = 1;
cluster_elecxsamp = no_avg_cbpt.diff_occ_vis.negclusterslabelmat == clust_no;
tmp = find(any(cluster_elecxsamp, 2))'; % find the electrodes inside the cluster
clust{3} = chan_lab(tmp);

clust_no = 3;
cluster_elecxsamp = no_avg_cbpt.diff_occ_vis.negclusterslabelmat == clust_no;
tmp = find(any(cluster_elecxsamp, 2))'; % find the electrodes inside the cluster
clust{4} = chan_lab(tmp);

clust_no = 1;
cluster_elecxsamp = no_avg_cbpt.diff_occ_vis.posclusterslabelmat == clust_no;
tmp = find(any(cluster_elecxsamp, 2))'; % find the electrodes inside the cluster
clust{5} = chan_lab(tmp);

clust_no = 2;
cluster_elecxsamp = no_avg_cbpt.diff_occ_vis.posclusterslabelmat == clust_no;
tmp = find(any(cluster_elecxsamp, 2))'; % find the electrodes inside the cluster
clust{6} = chan_lab(tmp);

clust_no = 3;
cluster_elecxsamp = no_avg_cbpt.diff_occ_vis.posclusterslabelmat == clust_no;
tmp = find(any(cluster_elecxsamp, 2))'; % find the electrodes inside the cluster
clust{7} = chan_lab(tmp);

save clust_chans.mat clust
