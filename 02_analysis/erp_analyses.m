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
parent_dir = filepath_parts(1:end-1);
parent_dir = strjoin(parent_dir, filesep);
parent_dir_2 = filepath_parts(1:end-2);
parent_dir_2 = strjoin(parent_dir_2, filesep);
output_dir = strjoin([parent_dir_2, "Emulation-Data-Output"], filesep);
mean_matrices_path = strjoin([output_dir, 'mean_matrices'], filesep);
mkdir(mean_matrices_path);
mean_matrices_peaks_epoched_path = strjoin([mean_matrices_path, 'peaks'], filesep);
mkdir(mean_matrices_peaks_epoched_path);
peak_erp_plots_dir = strjoin([parent_dir, "plots", "erp_plots", "peaks"], filesep);
mkdir(peak_erp_plots_dir);
% output
csd_dir = strjoin([output_dir, "csd_transform"], filesep);
% input
epochs_plus_error_dir = strjoin([output_dir, "07_epochs_with_extra_fields"], filesep);
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
mean_struct = renameStructField(mean_struct, 'diff_occ', 'diff_vis_occ');
z_mean_struct = renameStructField(z_mean_struct, 'diff_occ', 'diff_vis_occ');

cond_labs = {'all', 'occ', 'vis', 'rand1', 'const', 'rand2', 'diff_vis_occ', ...
    'diff_const_rand1', 'diff_const_rand2', 'diff_rand1_rand2'};

avg_peak_dist = 300;
% average across subjects
for lab = 1:length(cond_labs)
    z_epoch_mean.(cond_labs{lab}) = mean(z_mean_struct.(cond_labs{lab}),3);
    epoch_mean.(cond_labs{lab}) = mean(mean_struct.(cond_labs{lab}),3);

end

% specify parameters for (topo)plotting
base_dur = 500;
num_topos = 20;
topo_rows = 5;
topo_cols = 4;
topo_dur = 30;
epoch_end = 600 - topo_dur;
latencies_onset = linspace(0, epoch_end, num_topos);
latencies_offset = latencies_onset + topo_dur;


%% Prepare for Topoplotting

clear topo_struct

% extract means across topo_dur ms periods for topoplot
for lab = 1:length(cond_labs)
    topo_struct.(cond_labs{lab}) = nan(size(z_epoch_mean.(cond_labs{lab}),1), length(latencies_offset));
    for i = 1:length(latencies_onset)
        % get sample that is closest to on- and offset
        [~, lat_idx_onset] = min(abs(z_mean_struct.time_vec - latencies_onset(i)));
        [~, lat_idx_offset] = min(abs(z_mean_struct.time_vec - latencies_offset(i)));
        % get the mean across that time window
        topo_struct.(cond_labs{lab})(:, i) = mean(...
            z_epoch_mean.(cond_labs{lab})(:, lat_idx_onset:lat_idx_offset), 2);
    end
end

chan_locs = ALLEEG(1).chanlocs;

% make sure all plots have the same color scale
for lab = 1:length(cond_labs)
        zlims_cond(:, lab) = [min(topo_struct.(cond_labs{lab}), [], 'all'), max(topo_struct.(cond_labs{lab}), [], 'all')]';
end
zlims_all =  [min(zlims_cond(1,:)), max(zlims_cond(2,:))];


%% Create Topoplot of the time points of interest

cd(peak_erp_plots_dir)

% doc topoplot: https://rdrr.io/cran/erpR/man/topoplot.html

% all conditions

for lab = 1:length(cond_labs)
    figure(lab)
        for i = 1:length(latencies_onset)

        subplot(topo_rows, topo_cols, i)
        topoplot(topo_struct.(cond_labs{lab})(:, i), chan_locs, 'maplimits',  zlims_all)
        title(strcat(['average of ', num2str(latencies_onset(i)), ' to ', ...
            num2str(latencies_offset(i)), ' ms']), 'Interpreter', 'none')
        colorbar()
        end
        sgtitle(strjoin(["epochs around peaks, ", cond_labs{lab} ...
            ", CSD transformed, normalized and weighted per subj"]), 'Interpreter', 'none')
        fig_name = strjoin([cond_labs{lab}, "topo_csd_z"], "_")
        savefig(fig_name)
end

%% Prepare ERP Plotting 

cd(mean_matrices_peaks_epoched_path)

base_dur = 500;

load('count_peaks_csd.mat')

title_string = strcat(...
    ['ERP CSD transformed of all subjects averaged across epochs around peaks, n = ', ...
    num2str(count_peaks.all), ' epochs']);
title_string_z = strcat(...
    ['z-normalized ERP, CSD transformed of all subjects averaged across epochs around peaks, n = ', ...
    num2str(count_peaks.all), ' epochs']);
subtitle_string = strcat(['all channels']);


%% Cluster Based Permutation Test Attempt with Fieldtrip

% Script by Adriana Böttcher
% Adapted by Saskia Wilken
% which of these clusters are significant? 


%% add fieldtrip path
% eeglab;
% close;

% fieldtrip path
addpath 'C:\Users\swilk\AppData\Roaming\MathWorks\MATLAB Add-Ons\Collections\FieldTrip'
% file_path = matlab.desktop.editor.getActiveFilename;
% file_path = fileparts(file_path);
% addpath(file_path);
% filepath_parts = strsplit(file_path, filesep);
% parent_dir = filepath_parts(1:end-1);
% parent_dir = strjoin(parent_dir, filesep);
% parent_dir_2 = filepath_parts(1:end-2);
% parent_dir_2 = strjoin(parent_dir_2, filesep);
% output_dir_AR = strjoin([parent_dir_2, "Emulation-Data-Output\06_artifact_rejection"], filesep);
% output_dir_AR_peaks = strjoin([output_dir_AR, 'peaks_-500_750'], filesep);
% input_path = strjoin([output_dir_AR, 'peaks_-500_750'], filesep);
% addpath([char(parent_dir) filesep 'functions']);

ft_defaults;

% load('fieldtrip_EEG_KJP_elec_61');
% load('fieldtrip_EEG_KJP_neighbours_61.mat');
% elec.label = upper(elec.label);


%% Prepare Cluster Based Permutaiton Testing

cd(mean_matrices_peaks_epoched_path)
load("condition_indices_csd.mat")

clear all_eeg_field all_eeg_cell TMPEEG
% All data, compare with zero
% convert to fieldtrip:
for lab = find(~contains(cond_labs, 'diff'))
    for s = 1:length(ALLEEG)
        TMPEEG = ALLEEG(s);
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

% I NEED THIS: ALTERNATIVELY, I CAN GENERATE IT MYSELF
% load('fieldtrip_EEG_KJP_elec_61');
% load('fieldtrip_EEG_KJP_neighbours_61.mat');
% make sure the plotting of the CBP tests works properly. 

% load(['C:\wilken\Emulation-Data-Input' filesep 'freq_all_theta']);
%
% % was already here:
% constant = cellfun(@(x) x.constant, freq_all, 'UniformOutput', false);
% rand1 = cellfun(@(x) x.rand1, freq_all, 'UniformOutput', false);
% rand2 = cellfun(@(x) x.rand2, freq_all, 'UniformOutput', false);
% % freq all is a cell 1x32, each cell contains subject data and it is split
% % in two cells, one for each condition.
% erp_cond = mean_struct.all; % ??
% zero_cond = zeros(size(mean_struct.all)); % ??



%% init
% conds = {'constant', 'rand1', 'rand2'};
% conds: zero, erp signal at ~200
% cond_erp = erp_signal(:,'200ms') % pseudo code
% cond_zero = zeros(size(erp_signal, 1), size(erp_signal, 2)) % pseudo code

% weirdly complicated way to specify design matrix
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

cfg_neighb        = [];
cfg_neighb.method = 'distance';
neighbours        = neighbours; % ft_prepare_neighbours(cfg_neighb, dataFC_LP);

% CLUSTER TEST
% erp vs zero
calpha  = 0.001;
alpha   = 0.001;

% cfg is the configuraiton structure of fieldtrip
cfg                     = [];
cfg.design              = design;
cfg.uvar                = 1;
cfg.ivar                = 2;
cfg.channel             = {'all'};
cfg.avgovertime         = 'yes';
cfg.method              = 'montecarlo';
cfg.statistic           = 'depsamplesT'; % really indepsamples??? With comparison against 0
cfg.correctm            = 'cluster';
cfg.clusteralpha        = calpha;               % 0.05;
cfg.clusterstatistic    = 'maxsum';
cfg.minnbchan           = 0;
cfg.neighbours          = neighbours;
cfg.tail                = 0;
cfg.clustertail         = 0;
cfg.alpha               = alpha;               % 0.025;
cfg.numrandomization    = 1000;

clear stats

cbpt_cond_labs = fieldnames(all_eeg_cell); % condition labels
diff_labels = {'diff_all_zero', 'diff_vis_occ', 'diff_const_rand1', ...
    'diff_const_rand2', 'diff_rand1_rand2'}; % contrast labels
cond_comp_order = [1, 2, 5, 5, 4;
    7, 3, 4, 6, 6]; % conditions to compare
latency = [latencies_onset; latencies_offset]/1000;

% actual cluster tests
for lab = 1:length(diff_labels)
    for lat = 1:length(latency)
        cfg.latency  = latency(:,lat)'; % time interval over which the experimental
        % conditions must be compared (in seconds)
        tmp = ft_timelockstatistics(cfg, ...
            all_eeg_cell.(cbpt_cond_labs{cond_comp_order(1,lab)}){:}, ...
            all_eeg_cell.(cbpt_cond_labs{cond_comp_order(2,lab)}){:});
        % if the test yields no significant clusters, add the fields to
        % avoid an error. 
        if ~isfield(tmp, 'posclusters')
            tmp.posclusters = [];
            tmp.posclusterslabelmat = [];
            tmp.posdistribution = [];
            tmp.negclusters = [];
            tmp.negclusterslabelmat = [];
            tmp.negdistribution = [];
        end
        stats.(diff_labels{lab})(lat) = tmp; % store results in structure
        clear tmp
    end
end

% The field prob contains the proportion of draws from the permutation
% distribution with a maximum cluster-level statistic that is larger than
% the cluster-level test statistic
% stats.prob
% stats.stat %  cluster-level test statistic (here with maxsum: the sum of
% the T-values in this cluster).

%% Plot the CBPT Results

% % make sure all plots have the same color scale
for lab = 1:2 %length(diff_labels)
    for lat = 1:2%length(latency)
        zlims_cond_cbpt(:, lab, lat) = [min( ...
            stats.(diff_labels{lab})(lat).stat, [], 'all'), max( ...
            stats.(diff_labels{lab})(lat).stat, [], 'all')]';
    end
end

zlims_all_cbpt =  [min(min(zlims_cond_cbpt(1,:, :))), max(max(zlims_cond_cbpt(2,:, :)))];
marker_style = 'labelpoint';
head_style = 'both';
topo_rows = 5;
topo_cols = 4;
cd(peak_erp_plots_dir);
chan_locs = ALLEEG(1).chanlocs;

for i = 1:length({chan_locs.labels})
    chan_locs(i).labels = upper(chan_locs(i).labels);
end
isequal({chan_locs.labels}', stats.diff_all_zero.label); % check that chan_

close all

cd(peak_erp_plots_dir)

% actual plot code
for lab = 1:length(diff_labels)
    figure(lab)
    for lat = 1:length(latencies_onset)
        % prepare plot
        subplot(topo_rows, topo_cols, lat);
        if ~isempty(stats.(diff_labels{lab})(lat).posclusters)
        % Make a vector of all p-values associated with the clusters from ft_timelockstatistics.
        pos_cluster_pvals = [stats.(diff_labels{lab})(lat).posclusters(:).prob];

        % Then, find which clusters are deemed interesting to visualize, here we use a cutoff criterion based on the
        % cluster-associated p-value, and take a 5% two-sided cutoff (i.e. 0.025 for the positive and negative clusters,
        % respectively
        pos_clust = find(pos_cluster_pvals < calpha);
        pos       = ismember(stats.(diff_labels{lab})(lat).posclusterslabelmat, pos_clust);

        % and now for the negative clusters...
        neg_cluster_pvals = [stats.(diff_labels{lab})(lat).negclusters(:).prob];
        neg_clust         = find(neg_cluster_pvals < calpha);
        neg               = ismember(stats.(diff_labels{lab})(lat).negclusterslabelmat, neg_clust);

        clust_chans = unique([find(neg); find(pos)]);
%         chans_to_plot = find(stats.(diff_labels{lab})(i).mask);
        end
        topoplot(stats.(diff_labels{lab})(:, lat).stat, chan_locs, 'maplimits',  ...
            zlims_all_cbpt, 'electrodes','on', 'ecolor', 'black',...
            'emarker2', {clust_chans, '*', 'white'}, 'style', head_style);
        title(strcat(['average of ', num2str(latencies_onset(lat)), ' to ', ...
            num2str(latencies_offset(lat)), ' ms']), 'Interpreter', 'none')
        a = colorbar;
        position = get(a,'Position');
        ylabel(a,'T-Values','FontSize',7,'Rotation',270, 'Position', [position(1)+3 position(2)]);

        head_style = 'both';
        marker_style = 'labelpoint';
        clear neg neg_clust pos pos_clust neg_cluster_pvals pos_cluster_pvals

    end
    sgtitle(strjoin(["CBPT, epochs ar. peaks,", ...
        "contr: ", diff_labels{lab}, ...
        ", CSD, normalized & weighted per subj, calpha = ", ...
        calpha, ", alpha = ", alpha]), 'Interpreter', 'none')
    fig_name = strjoin([diff_labels{lab}, "topo_CBPT_csd_z"], "_")
    savefig(fig_name);

end

% https://github.com/fieldtrip/fieldtrip/blob/master/ft_clusterplot.m
% https://www.fieldtriptoolbox.org/tutorial/plotting/
% https://github.com/fieldtrip/fieldtrip/blob/release/ft_topoplotER.m

%% Plot ERPs as grand average and save plots

cd(peak_erp_plots_dir)

% Average ERP Plots
for lab = 1:length(cond_labs)
    figure
    title_string_z = strjoin(...
        ["z-normalized ERP, CSD transformed of all subjects averaged across epochs around peaks, n = ", ...
        num2str(count_peaks.(cond_labs{lab})), " epochs"]);
        subtitle_string = cond_labs{lab};

    [~, min_ind] = plot_ERP(z_epoch_mean.(cond_labs{lab}), mean_struct.time_vec, base_dur, title_string_z, subtitle_string);
    fig_name = strjoin([cond_labs{lab}, "butterfly_csd_z"], "_")
    savefig(fig_name)
end


%% Get electrodes of largest deflection
% prepare single-channel plotting

clear index_struct
interesting_lats = [185, 245, 275, 575]; % lower limit of latencies that 
lat_labs = {"lat_180_210", "lat_240_270", "lat_270_300", "lat_570_600"};
num_chans = 1; % number of minimal/maximal channels to extract
% might be interesting
diff_labels = {'diff_all_zero', 'diff_vis_occ', 'diff_const_rand1', ...
    'diff_const_rand2', 'diff_rand1_rand2'}; % contrast labels
% index_struct.diff_all_zero. = struct();
for lab = 2:length(diff_labels)
%     index_struct.(diff_labels{lab}) = struct();
    for lat = 1:length(interesting_lats)
        % find minimal deflection
        interesting_window = max(find(latencies_onset < interesting_lats(lat))); % find the topo window
        % that is the most interesting in previous ERP plots.
        [~, min_chan_idx] = mink(topo_struct.(diff_labels{lab})(:, ...
            interesting_window),num_chans); % find the min_chan channels
        % with the lowest deflection
        [~, max_chan_idx] = maxk(topo_struct.(diff_labels{lab})(:, ...
            interesting_window),num_chans); % find the min_chan channels
        % with the highest deflection
        chan_lab_min = {ALLEEG(1).chanlocs(min_chan_idx).labels};
        chan_lab_max = {ALLEEG(1).chanlocs(max_chan_idx).labels};

        index_struct.(diff_labels{lab}).(num2str(lat_labs{lat}))(2).idx = lat_labs(lat);
        index_struct.(diff_labels{lab}).(num2str(lat_labs{lat}))(2).idx = interesting_window;
        index_struct.(diff_labels{lab}).(num2str(lat_labs{lat}))(1).min = chan_lab_min;
        % min channels are the ones where the first condition is more
        % negative than the second
        index_struct.(diff_labels{lab}).(num2str(lat_labs{lat}))(1).max = chan_lab_max;
        % max channels are the ones where the first condition is more
        % positive than the second
        index_struct.(diff_labels{lab}).(num2str(lat_labs{lat}))(2).min = min_chan_idx;
        index_struct.(diff_labels{lab}).(num2str(lat_labs{lat}))(2).max = max_chan_idx;
        clear chan_lab_max chan_lab_min max_chan_idx min_chan_idx interesting_window
    end
end
sign([[-0.5,5, 0.01, -54 ];[1.5, -5, 11, 50 ]] )


%% Single Channel ERP Plot + Topo
% Plots the electrode that shows the biggest amplitude in the ERP
minmax = {'min', 'max'};
num_figs = 1:[length(lat_labs) * length(diff_labels)-1]*length(minmax);
count = 1;
z_epoch_mean = renameStructField(z_epoch_mean, 'diff_occ', 'diff_vis_occ');
cond_comp_order = [1, 3, 5, 5, 4;
                   7, 2, 4, 6, 6];
close all

for lab = 2:length(diff_labels)
    for lat = 1:length(lat_labs)
        for m = 1:2
            chan_lab = index_struct.(diff_labels{lab}).(num2str(lat_labs{lat}))(1).(minmax{m});
            chan_idx = index_struct.(diff_labels{lab}).(num2str(lat_labs{lat}))(2).(minmax{m});
            lat_idx = index_struct.(diff_labels{lab}).(lat_labs{lat})(2).idx;
            cond{1} = cond_labs{cond_comp_order(1,lab)};
            cond{2} = cond_labs{cond_comp_order(2,lab)};
            % ERP plot
            figure(num_figs(count))
            count = 1 + count;
            subplot(1, 2, 1)
            title_string_z = strjoin([...
                "channel ", chan_lab, ", all subjects, epoch around peak"],"")
            subtitle_string = strjoin([cond{1}, " - ", cond{2}], "");
            hold on
            plot(mean_struct.time_vec, z_epoch_mean.(cond{1})(chan_idx, :), 'LineWidth', 1)
            plot(mean_struct.time_vec, z_epoch_mean.(cond{2})(chan_idx, :), 'LineWidth', 1)
            xline(0, 'LineWidth', 1)
            xline(latencies_onset(lat_idx), 'LineWidth', 1, 'Color', 'b', 'Linestyle', ':')
            xline(latencies_offset(lat_idx), 'LineWidth', 1, 'Color', 'b', 'Linestyle', ':')
            hold off
            title(title_string_z)
            subtitle(subtitle_string)
            legend({cond{:}, "peak", "topo on and offset"})

            % topoplot
            subplot(1, 2, 2)
            topoplot(topo_struct.(diff_labels{lab})(:, lat_idx), chan_locs, 'maplimits',  ...
                zlims_all, 'emarker2', {chan_idx,'*','m'}, 'electrodes','labelpoint')
            title(strcat(['average of ', num2str(latencies_onset(lat_idx)), ' to ', ...
                num2str(latencies_offset(lat_idx)), ' ms, ', diff_labels{lab}]), ...
                'Interpreter', 'None');
            a = colorbar;
            position = get(a,'Position');
            ylabel(a,'z-Values','FontSize',12,'Rotation',270, 'Position', [position(1)+1 position(2)]);

            % figure
            sgtitle(strjoin([diff_labels{lab}, ...
                ', ep ar. peaks, all subjects, CSD & normalized & weighted, all trials, ', ...
                chan_lab, ', ' minmax{m}, ' amp in time window ', num2str(latencies_onset(lat_idx)), ' to ', ...
                num2str(latencies_offset(lat_idx)), ' ms'],""), 'Interpreter', 'None')
            file_name = strjoin([diff_labels{lab}, "ERP_topo_single_csd_z"], "_")

%             savefig('file_name')
        end
    end
end


%% ERP image

clear z_ALLEEG TMPEEG z_mean_struct mean_struct z_tmpEEG
sum(arrayfun(@(s) length(s.epoch), ALLEEG)) % how many epochs we need. 
s = 1;
while s <= length(ALLEEG)
    EEG = ALLEEG(s);
    del = 0;
    event_idx = zeros(length(EEG.epoch),1);
    for ep = 1:length(EEG.epoch)
        event_idx(ep) = get_epoch_center(EEG, ep);
%         if ~(strcmp(EEG.event(ev-del).type, 'S 40') | strcmp(EEG.event(ev-del).type, 'S 50'))
%             EEG.event(ev-del) = [];
%             del = del +1;
%         end
    end
    log_idx = unfind(event_idx, length(EEG.event));
    EEG.event(~log_idx) = [];
    EEG = eeg_checkset(EEG);
    ALLEEG(s) = EEG;
    s = s + 1;
end
% WATCH OUT: FROM NOW ON get_epoch_center won't work anymore since the
% EEG.epoch structure is no longer aligned with the EEG.event structure. 

for s = 1:length(ALLEEG)
    EEG = ALLEEG(s);
    EEG = pop_selectevent( EEG, 'artifact_error',{'VALID'},'deleteevents',...
        'off','deleteepochs','on','invertepochs','off');
    ALLEEG(s) = EEG;
end
CONDEEG = ALLEEG;
% remove occluded VISIBLE
for s = 1:length(CONDEEG)
    EEG = CONDEEG(s);
    EEG = pop_selectevent( EEG, 'OCCL',{'OFF'},'deleteevents',...
        'off','deleteepochs','on','invertepochs','off');
    CONDEEG(s) = EEG;
end
CONDEEG2 = ALLEEG;
% remove visible OCCLUDED
for s = 1:length(CONDEEG2)
    EEG = CONDEEG2(s);
    EEG = pop_selectevent( EEG, 'OCCL',{'ON'},'deleteevents',...
        'off','deleteepochs','on','invertepochs','off');
    CONDEEG2(s) = EEG;
end
% remove random CONSTANT
for s = 1:length(MODEEG)
    EEG = CONDEEG(s);
    EEG = pop_selectevent( EEG, 'TRAJ',{'RANDOM1', 'RANDOM2'},'deleteevents',...
        'off','deleteepochs','on','invertepochs','off');
    CONDEEG(s) = EEG;
end
% remove constant RANDOM
for s = 1:length(MODEEG)
    EEG = CONDEEG(s);
    EEG = pop_selectevent( EEG, 'TRAJ',{'RANDOM1', 'RANDOM2'},'deleteevents',...
        'off','deleteepochs','on','invertepochs','off');
    CONDEEG(s) = EEG;
end


% clear MODEEG
fieldnames(cond_ind)

interesting_lats = [185, 245, 275, 575]; % lower limit of latencies that 
lat = 1;
interesting_window = max(find(latencies_onset < interesting_lats(lat))); % find the topo window
% that is the most interesting in previous ERP plots.
[~, min_chan_idx] = mink(topo_struct.('all')(:, ...
    interesting_window),num_chans); % find the min_chan channels
% with the lowest deflection
[~, max_chan_idx] = maxk(topo_struct.('all')(:, ...
    interesting_window),num_chans); % find the min_chan channels
chan_lab_min = {ALLEEG(1).chanlocs(min_chan_idx).labels};
chan_lab_max = {ALLEEG(1).chanlocs(max_chan_idx).labels};

% TMPEEG.data = TMPEEG.data(:,:,cond_ind(s).(cond_labs{lab})); % this is for selecting conditions. 
% initialize dataset
TMPEEG = CONDEEG(1);
TMPEEG.data = normalize(CONDEEG(1).data,2, 'zscore');
% put all epochs of all subjects in one dataset
% get all data and all event structures
lab = find(strcmp(cond_labs, 'all'));
for s = 2:length(CONDEEG)
    % normalize each subjects' data
    z_tmpEEG = CONDEEG(s);
    % normalize data
    z_tmpEEG.data = [];%MODEEG(s).data(:,:,cond_ind(s).(cond_labs{lab}));
    z_tmpEEG.data = normalize(CONDEEG(s).data,2, 'zscore');
    % get location where new data should go
%     data_idx = size(TMPEEG.data, 3)+1:size(TMPEEG.data, 3)+size(z_tmpEEG.data, 3);
    % concatenate subject data in one dataset
    TMPEEG.data = cat(3, z_tmpEEG.data, TMPEEG.data);
%      TMPEEG.data(:,:,data_idx) =
    % get epoch vectors
    epoch_vec = TMPEEG.event(end).epoch + [z_tmpEEG.event.epoch];
    % get event indices
    event_idx = size(TMPEEG.event,2)+1:size(TMPEEG.event,2)+size(z_tmpEEG.event, 2);
    % add events from previous subject to TMPEEG
    TMPEEG.event(:, event_idx) = z_tmpEEG.event;
    for ev = 1:length(event_idx)
        TMPEEG.event(:, event_idx(ev)).epoch = epoch_vec(ev);
    end
    TMPEEG.epoch(:, size(TMPEEG.epoch,2)+1:...
        size(TMPEEG.epoch,2)+size(z_tmpEEG.epoch, 2)) = z_tmpEEG.epoch;
    epoch_vec = [];
end
TMPEEG.trials = length(TMPEEG.epoch);

TMPEEG2 = CONDEEG2(1);
TMPEEG2.data = normalize(CONDEEG2(1).data,2, 'zscore');
% put all epochs of all subjects in one dataset
% get all data and all event structures
lab = find(strcmp(cond_labs, 'all'));
for s = 2:length(CONDEEG2)
    % normalize each subjects' data
    z_tmpEEG = CONDEEG2(s);
    % normalize data
    z_tmpEEG.data = [];%MODEEG(s).data(:,:,cond_ind(s).(cond_labs{lab}));
    z_tmpEEG.data = normalize(CONDEEG2(s).data,2, 'zscore');
    % get location where new data should go
%     data_idx = size(TMPEEG.data, 3)+1:size(TMPEEG.data, 3)+size(z_tmpEEG.data, 3);
    % concatenate subject data in one dataset
    TMPEEG2.data = cat(3, z_tmpEEG.data, TMPEEG2.data);
%      TMPEEG.data(:,:,data_idx) =
    % get epoch vectors
    epoch_vec = TMPEEG2.event(end).epoch + [z_tmpEEG.event.epoch];
    % get event indices
    event_idx = size(TMPEEG2.event,2)+1:size(TMPEEG2.event,2)+size(z_tmpEEG.event, 2);
    % add events from previous subject to TMPEEG
    TMPEEG2.event(:, event_idx) = z_tmpEEG.event;
    for ev = 1:length(event_idx)
        TMPEEG2.event(:, event_idx(ev)).epoch = epoch_vec(ev);
    end
    TMPEEG2.epoch(:, size(TMPEEG2.epoch,2)+1:...
        size(TMPEEG2.epoch,2)+size(z_tmpEEG.epoch, 2)) = z_tmpEEG.epoch;
    epoch_vec = [];
end

TMPEEG2.trials = length(TMPEEG2.epoch);
count = 1;
clear chan_lab_cell
clear chan_lab_vec
diff_labels = {'diff_all_zero', 'diff_vis_occ', 'diff_const_rand1', ...
    'diff_const_rand2', 'diff_rand1_rand2'};
% get vector of all interesting channels
% for lab = 2:length(diff_labels)
%     for lat = 1:length(lat_labs)
%         for m = 1:2
%             chan_lab_cell{count} = index_struct.(diff_labels{lab}).(num2str(lat_labs{lat}))(1).(minmax{m}){:};
%             chan_idx_vec(count) = index_struct.(diff_labels{lab}).(num2str(lat_labs{lat}))(2).(minmax{m})(:);
%             count = count + 1;
%         end
%     end
% end
% chan_lab_cell_plus = {chan_lab_cell{:},chan_lab_min{:},chan_lab_max{:}}
% chan_idx_vec_plus = [chan_idx_vec,min_chan_idx,max_chan_idx]

% actual erp image parameters
% sorting_cond = {'epoch_error', 'pursuit_lat', 'epoch_error_z', 'pursuit_lat_z'};
% cond = 1;
smooth = round(length(TMPEEG.epoch)/15,0); % number of trials
dec_fac = smooth/2; % ratio of trials in to plot out. 
% if the smoothing width is larger than twice the decimation factor, no
% info is lost. 
plot_type = 1; % 1 = channel, 0 = component
project_channel = [[]];
sort_type = {'S 40' 'S 50'};
% sort_win = [-0.1 0.1];
% align = 0;
unit = 'z-values';
% chan = 32;
% chan_lab = chan_lab_cell{chan};
chan_lab = 'P8'
chan_no = find(strcmp({ALLEEG(s).chanlocs.labels}, chan_lab));
alpha = 0.001;
% chan_no = chan_idx_vec(chan);
% eeg_c = lab;
% eeg_condition = cond_labs{eeg_c};
eeg_condition = 'visible';
% sort_event_field = sorting_cond{cond};
sort_event_field = 'pursuit_lat';
% [min_val, max_val] = bounds(TMPEEG.data, 'all')
% prctile(TMPEEG.data(:,:,1:10000), [1, 10, 25, 50, 75, 90, 100], 'all')
% prctile(ALLEEG(s).data(:,:,1:1000), [1, 10, 25, 50, 75, 90, 100], 'all')
% prctile(CONDEEG(s).data(:,:,1:1000), [1, 10, 25, 50, 75, 90, 100], 'all')
% prctile(z_tmpEEG.data(:,:,1:1000), [1, 10, 25, 50, 75, 90, 100], 'all')

% for cond = 1:length(sorting_cond)
%     for chan = 1:length(chan_idx_vec)
%         for eeg_c = 1:length(cond_labs)
        % dynamically changing image parameters
% 
% [ALLEEG TMPEEG index] = eeg_store(ALLEEG, TMPEEG);
% eeglab redraw

% plotting of ERP image
figure; 
title = strjoin(["All subjects, normalized & weighted, csd, Smoothing factor "...
    smooth,", decimate " dec_fac ", chan " chan_lab],"");
pop_erpimage(TMPEEG,plot_type, chan_no,[[]], title,smooth, dec_fac, ...
    sort_type,[],sort_event_field ,'yerplabel', unit,'erp', 'on', ...
    'cbar','on','cbar_title',unit, 'topo', ...
    { chan_no TMPEEG.chanlocs TMPEEG.chaninfo });
sgtitle(strjoin(["contrast: " eeg_condition, ", sorted by: ", ...
    sort_event_field " from top (higher) to bottom (lower)"],""), 'Interpreter', 'None');

eeg_condition = 'occluded';
smooth = round(length(TMPEEG2.epoch)/15,0); % number of trials
dec_fac = smooth/2; % ratio of trials in to plot out. 
% plotting of ERP image
figure; 
title = strjoin(["All subjects, normalized & weighted, csd, Smoothing factor "...
    smooth,", decimate " dec_fac ", chan " chan_lab],"");
pop_erpimage(TMPEEG2,plot_type, chan_no,[[]], title,smooth, dec_fac, ...
    sort_type,[],sort_event_field ,'yerplabel', unit,'erp', 'on', ...
    'cbar','on','cbar_title',unit, 'topo', ...
    { chan_no TMPEEG2.chanlocs TMPEEG2.chaninfo });
sgtitle(strjoin(["contrast: " eeg_condition, ", sorted by: ", ...
    sort_event_field " from top (higher) to bottom (lower)"],""), 'Interpreter', 'None');



