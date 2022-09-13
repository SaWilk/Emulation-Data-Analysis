% CSD transformation

% this script calculated the surface laplacians for the eeg data and saves
% the transformed ERP datasets

% author: Saskia Wilken
% creation date: 05.09.2022 

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


%% Load data for CSD transformation

cd(epochs_plus_error_dir);
%list all *.set files in inputpath
file_names = dir('*.set');
%concatenate into one cell array
files2read = {file_names.name};
eeglab;
close all

for idx = 1:length(files2read)

    EEG = pop_loadset('filename',files2read{idx});
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG);

end


%% Surface Laplacians CSD toolbox (aka current source density)

% keep default settings
addpath(genpath('C:\wilken\CSDtoolbox')); % CSD toolbox needs to be installed for this to work

% tutorial: https://psychophysiology.cpmc.columbia.edu/Software/CSDtoolbox/tutorial.html
% common errors: https://psychophysiology.cpmc.columbia.edu/Software/CSDtoolbox/errors.html
% FAQ: https://psychophysiology.cpmc.columbia.edu/Software/CSDtoolbox/FAQ.html

% exploit the EEG montage information included in any EEGlab data file to 
% generate the montage-dependent "Channels × Channels" transformation 
% matrices G and H (EEGlab_Make_G_H.m), and how to use these 
% transformation matrices with individual EEGlab data files

% Get usable list of electrodes from EEGlab data structure
for site = 1:length({ALLEEG(1).chanlocs.labels})
    trodes{site}=(ALLEEG(1).chanlocs(site).labels);
end
trodes=trodes';

% attempt to get more precise electrode locations, but ended poorly. 
% cd('C:\wilken\Emulation-Data-Input\EEG_montage')
% chan_locs = loadbvef('BC-64.bvef') % chan locs we use in Hagen
% locs = ALLEEG(1).chanlocs
% montage_64.("lab") = {locs.labels}.';
% montage_64.("theta") = [locs.sph_theta].';
% montage_64.("phi") = [locs.sph_phi].';
% montage_64.("xy") = [[locs.X].',[locs.Y].'];
% montage_64.xy(:,1) = montage_64.xy(:,1)*(-0.5)+0.5;
% montage_64.xy(:,2) = montage_64.xy(:,2)*0.5+0.5;
% MapMontage(montage_64); % looks even worse ^^
figure; 
topoplot([],ALLEEG(1).chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo);

% Get Montage for use with CSD Toolbox
% '10-5-System_Mastoids_EGI129.csd' is a data file with electrode locations
% according to the 10-5 system. I manually added the strange electrodes 
% P11, P12, O9 and O10. However, the final electrode locations still look a
% bit fishy. A problem is that EEGLAB apparently uses a different scaling
% for the X and Y coordinates, which makes the approach above not very
% helpful
Montage_64=ExtractMontage('10-5-System_Mastoids_EGI129_SW.csd',trodes); % _SW is my custom version
% results in a montage of the following shape:
%   struct with fields:
%       lab: {60×1 cell}
%     theta: [60×1 double]
%       phi: [60×1 double]
%        xy: [60×2 double]

MapMontage(Montage_64); % a wonky eeg location plot. 


%% Derive G and H!
[G,H] = GetGH(Montage_64); % default m = 4

%% Save G and H to later import when doing the CSD transform on files
% save('G:\PhysioData\MN_Fear\G.mat', 'G');
% save('G:\PhysioData\MN_Fear\H.mat', 'H');

% revised method to store G and H matrices with CSD montage for later import
Montage = Montage_64;                             % use generic variable name
save(strjoin([csd_dir,'CSD_montage_64.mat'], filesep), 'G', 'H', 'Montage'); % save variables to Matlab file
clear G H Montage;                                % remove variables from workspace
load(strjoin([csd_dir,'CSD_montage_64.mat'], filesep)); % restore variables to the workspace

%% Calculate the CSD

for s = 1:length(ALLEEG)
    tic
    for ne = 1:length(ALLEEG(s).epoch)               % loop through all epochs
        myEEG = single(ALLEEG(s).data(:,:,ne));      % reduce data precision to reduce memory demand
        MyResults = CSD(myEEG,G,H,1.0e-5, 9);            % compute CSD for <channels-by-samples> 2-D epoch
        % default lambda but 9 as head radius (slightly larger than 56
        % cm circumference)
        data(:,:,ne) = MyResults;              % assign data output
    end
    ALLEEG(s).("CSD_data") = double(data);          % final CSD data
    set_name = strsplit(ALLEEG(s).setname, '_');
    new_set_name = strjoin([set_name{1}, "CSD"], '_');
    ALLEEG(s).setname = new_set_name;
    clear data myEEG MyResults new_set_name
    ALLEEG(s) = pop_saveset(ALLEEG(s), 'filename', char(ALLEEG(s).setname), 'filepath', ...
        char(csd_dir));
    toc
end

figure()
topoplot(ALLEEG(s).CSD_data(:, 200, 1), EEG.chanlocs)
figure()
topoplot(ALLEEG(s).data(:, 200, 1), EEG.chanlocs)


%% Load CDS data

eeglab;

% long epochs
cd(csd_dir);
%list all *.set files in inputpath
file_names = dir('*.set');
%concatenate into one cell array
files2read = {file_names.name};
% load eeg data epoched around peaks
for idx = 1:length(files2read)
    ALLEEG2(idx) = pop_loadset('filename',files2read{idx});
end

ALLEEG = ALLEEG2;
clear ALLEEG2


%% Calculate means and normalize data

count_peaks.all = 0;
count_peaks.occ = 0;
count_peaks.vis = 0;
count_peaks.rand1 = 0;
count_peaks.const = 0;
count_peaks.rand2 = 0;

for s = 1:size(ALLEEG, 2)

    EEG = ALLEEG(s);
    % all conditions
    mean_struct.all(:,:,s) = mean(EEG.CSD_data,3);
    % get number of epochs per subject for weighting
    z_mean_struct.num_epochs_all(s) = size(EEG.CSD_data, 3);
    % z-normalized data
    z_mean_struct.all(:, :, s) = normalize(mean_struct.all(:, :, s),2, 'zscore');
    % get indices of peaks in event structure that also belong to a certain
    % condition
    occl_off_log = false(size(EEG.event,2),1);
    occl_on_log = false(size(EEG.event,2),1);
    traj_rand1_log = false(size(EEG.event,2),1);
    traj_const_log = false(size(EEG.event,2),1);
    traj_rand2_log = false(size(EEG.event,2),1);

    for ev = 1:size(EEG.event,2)
        % occlusion
        occl_off_log(ev) = strcmp(EEG.event(ev).type, 'S 40') && strcmp(EEG.event(ev).OCCL, 'OFF');
        occl_on_log(ev) = strcmp(EEG.event(ev).type, 'S 40') && strcmp(EEG.event(ev).OCCL, 'ON');
        % constant
        traj_rand1_log(ev) = strcmp(EEG.event(ev).type, 'S 40') && strcmp(EEG.event(ev).TRAJ, 'RANDOM1');
        traj_const_log(ev) = strcmp(EEG.event(ev).type, 'S 40') && strcmp(EEG.event(ev).TRAJ, 'CONST');
        traj_rand2_log(ev) = strcmp(EEG.event(ev).type, 'S 40') && strcmp(EEG.event(ev).TRAJ, 'RANDOM2');
    end

    % get epoch numbers for respective conditions
    cond_ind(s).epoch_vis = unique([EEG.event(occl_off_log).epoch]);
    cond_ind(s).epoch_occ = unique([EEG.event(occl_on_log).epoch]);
    cond_ind(s).epoch_rand1 = unique([EEG.event(traj_rand1_log).epoch]);
    cond_ind(s).epoch_const = unique([EEG.event(traj_const_log).epoch]);
    cond_ind(s).epoch_rand2 = unique([EEG.event(traj_rand2_log).epoch]);

    % calculate mean across condition
    mean_struct.occ_off(:,:,s) = mean(EEG.CSD_data(:,:,cond_ind(s).epoch_vis),3);
    mean_struct.occ_on(:,:,s) = mean(EEG.CSD_data(:,:,cond_ind(s).epoch_occ),3);
    mean_struct.rand1(:,:,s) = mean(EEG.CSD_data(:,:,cond_ind(s).epoch_rand1),3);
    mean_struct.const(:,:,s) = mean(EEG.CSD_data(:,:,cond_ind(s).epoch_const),3);
    mean_struct.rand2(:,:,s) = mean(EEG.CSD_data(:,:,cond_ind(s).epoch_rand2),3);

    % get z-normalized data
    z_mean_struct.occ_off(:, :, s) = normalize(mean_struct.occ_off(:, :, s),2, 'zscore');
    z_mean_struct.occ_on(:, :, s) = normalize(mean_struct.occ_on(:, :, s),2, 'zscore');
    z_mean_struct.rand1(:, :, s) = normalize(mean_struct.rand1(:, :, s),2, 'zscore');
    z_mean_struct.const(:, :, s) = normalize(mean_struct.const(:, :, s),2, 'zscore');
    z_mean_struct.rand2(:, :, s) = normalize(mean_struct.rand2(:, :, s),2, 'zscore');

    % get weights
    z_mean_struct.num_epochs_occ_off(s) = size(EEG.data(:,:,cond_ind(s).epoch_vis), 3);
    z_mean_struct.num_epochs_occ_on(s) = size(EEG.data(:,:,cond_ind(s).epoch_occ), 3);
    z_mean_struct.num_epochs_rand1(s) = size(EEG.data(:,:,cond_ind(s).epoch_rand1), 3);
    z_mean_struct.num_epochs_const(s) = size(EEG.data(:,:,cond_ind(s).epoch_const), 3);
    z_mean_struct.num_epochs_rand2(s) = size(EEG.data(:,:,cond_ind(s).epoch_rand2), 3);

    % calculate contrasts by simple subtraction
    mean_struct.diff_occ = mean_struct.occ_off(:,:,s) - mean_struct.occ_on(:,:,s);
    mean_struct.diff_const_rand1 = mean_struct.const(:,:,s) - mean_struct.rand1(:,:,s);
    mean_struct.diff_const_rand2 = mean_struct.const(:,:,s) - mean_struct.rand2(:,:,s);

    % get time vec
    mean_struct.time_vec = EEG.times;
    z_mean_struct.time_vec = EEG.times;


    % add peaks to know how much power each condition has
    count_peaks.all = size(EEG.CSD_data, 3) + count_peaks.all;
    count_peaks.occ = size(EEG.CSD_data(:,:,cond_ind(s).epoch_occ), 3) + count_peaks.occ;
    count_peaks.vis = size(EEG.CSD_data(:,:,cond_ind(s).epoch_vis), 3) + count_peaks.vis;
    count_peaks.rand1 = size(EEG.CSD_data(:,:,cond_ind(s).epoch_rand1), 3) + count_peaks.rand1;
    count_peaks.const = size(EEG.CSD_data(:,:,cond_ind(s).epoch_const), 3) + count_peaks.const;
    count_peaks.rand2 = size(EEG.CSD_data(:,:,cond_ind(s).epoch_rand2), 3) + count_peaks.rand2;

end

epoch_weights = z_mean_struct.num_epochs_all/mean(z_mean_struct.num_epochs_all);
epoch_weights_occ_on = z_mean_struct.num_epochs_all/mean(z_mean_struct.num_epochs_occ_on);
epoch_weights_occ_off = z_mean_struct.num_epochs_all/mean(z_mean_struct.num_epochs_occ_off);
epoch_weights_rand1 = z_mean_struct.num_epochs_all/mean(z_mean_struct.num_epochs_rand1);
epoch_weights_const = z_mean_struct.num_epochs_all/mean(z_mean_struct.num_epochs_const);
epoch_weights_rand2 = z_mean_struct.num_epochs_all/mean(z_mean_struct.num_epochs_rand2);


for ep = 1:size(z_mean_struct.all, 3)

    z_mean_struct.all(:,:,ep) = z_mean_struct.all(:,:,ep) * epoch_weights(ep);
    z_mean_struct.occ_on(:,:,ep) = z_mean_struct.occ_on(:,:,ep) * epoch_weights_occ_on(ep);
    z_mean_struct.occ_off(:,:,ep) = z_mean_struct.occ_off(:,:,ep) * epoch_weights_occ_off(ep);
    z_mean_struct.rand1(:,:,ep) = z_mean_struct.rand1(:,:,ep) * epoch_weights_rand1(ep);
    z_mean_struct.const(:,:,ep) = z_mean_struct.const(:,:,ep) * epoch_weights_const(ep);
    z_mean_struct.rand2(:,:,ep) = z_mean_struct.rand2(:,:,ep) * epoch_weights_rand2(ep);

end


%% Save the Mean Matrices for Convenient Loading

cd(mean_matrices_peaks_epoched_path)
save("mean_struct_csd.mat", 'mean_struct')
save("condition_indices_csd.mat", 'cond_ind')
save("count_peaks_csd.mat", 'count_peaks')
save("z_mean_struct_csd.mat", 'z_mean_struct')


%% Prepare for Topoplotting

cd(mean_matrices_peaks_epoched_path)

load('z_mean_struct_csd.mat')

avg_peak_dist = 300;
% average across subjects
z_epoch_mean_all = mean(z_mean_struct.all,3);
z_epoch_mean_occ = mean(z_mean_struct.occ_on,3, 'omitnan'); % omitnan is
% % necessary because subject 18 has no task B
z_epoch_mean_vis = mean(z_mean_struct.occ_off,3);
z_epoch_mean_rand1 = mean(z_mean_struct.rand1,3);
z_epoch_mean_const = mean(z_mean_struct.const,3);
z_epoch_mean_rand2 = mean(z_mean_struct.rand2,3);

% specify parameters for topoplotting
base_dur = 500;
num_topos = 20;
topo_rows = 5;
topo_cols = 4;
topo_dur = 30;
epoch_end = 600 - topo_dur;
latencies_onset = linspace(0, epoch_end, num_topos);
latencies_offset = latencies_onset + topo_dur;

% extract means across topo_dur ms periods for topoplot
for i = 1:length(latencies_onset)

    [~, lat_idx_onset] = min(abs(z_mean_struct.time_vec - latencies_onset(i)));
    [~, lat_idx_offset] = min(abs(z_mean_struct.time_vec - latencies_offset(i)));
    topo_struct.all(:, i) = mean(z_epoch_mean_all(:,lat_idx_onset:lat_idx_offset), 2);
    topo_struct.vis(:, i) = mean(z_epoch_mean_vis(:,lat_idx_onset:lat_idx_offset), 2);
    topo_struct.occ(:, i) = mean(z_epoch_mean_occ(:,lat_idx_onset:lat_idx_offset), 2);
    topo_struct.const(:, i) = mean(z_epoch_mean_const(:,lat_idx_onset:lat_idx_offset), 2);
    topo_struct.rand1(:, i) = mean(z_epoch_mean_rand1(:,lat_idx_onset:lat_idx_offset), 2);
    topo_struct.rand2(:, i) = mean(z_epoch_mean_rand2(:,lat_idx_onset:lat_idx_offset), 2);
%     topo_struct.diff_occ(:, i) = mean(z_mean_struct.diff_occ(:,lat_idx_onset:lat_idx_offset), 2);
%     topo_struct.diff_const_rand1(:, i) = mean(z_mean_struct.diff_const_rand1(:,lat_idx_onset:lat_idx_offset), 2);
%     topo_struct.diff_const_rand2(:, i) = mean(z_mean_struct.diff_const_rand2(:,lat_idx_onset:lat_idx_offset), 2);

end

chan_locs = ALLEEG(1).chanlocs


%% Create Topoplot of the time points of interest

cd(peak_erp_plots_dir)

% doc topoplot: https://rdrr.io/cran/erpR/man/topoplot.html


%% all conditions

figure()
zlims = [min(topo_struct.all, [], 'all'), max(topo_struct.all, [], 'all')]';
for i = 1:length(latencies_onset)

    subplot(topo_rows, topo_cols, i)
    topoplot(topo_struct.all(:, i), chan_locs, 'maplimits',  zlims)
    title(strcat(['average of ', num2str(latencies_onset(i)), ' to ', ...
        num2str(latencies_offset(i)), ' ms']))
    colorbar()
end
sgtitle(strcat(['epochs around peaks all subjects all trials, CSD transformed, normalized and weighted per subj']))

savefig('all_condition_peaks_topo_csd_z')


%% Prepare ERP Plotting 

cd(mean_matrices_peaks_epoched_path)

load('mean_struct_csd.mat')

avg_peak_dist = 300;
% average across subjects
epoch_mean_all = mean(mean_struct.all,3);
epoch_mean_occ = mean(mean_struct.occ_on,3, 'omitnan'); % omitnan is
% % necessary because subject 18 has no task B
epoch_mean_vis = mean(mean_struct.occ_off,3);
epoch_mean_rand1 = mean(mean_struct.rand1,3);
epoch_mean_const = mean(mean_struct.const,3);
epoch_mean_rand2 = mean(mean_struct.rand2,3);

load('z_mean_struct_csd.mat')

% get z-weighted mean for all subjects
z_epoch_mean_all = mean(z_mean_struct.all,3);
z_epoch_mean_occ_on = mean(z_mean_struct.occ_on,3, 'omitnan');
z_epoch_mean_occ_off = mean(z_mean_struct.occ_off,3, 'omitnan');
z_epoch_mean_rand1 = mean(z_mean_struct.rand1,3);
z_epoch_mean_const = mean(z_mean_struct.const,3);
z_epoch_mean_rand2 = mean(z_mean_struct.rand2,3);

base_dur = 500;

load('count_peaks.mat')

title_string = strcat(...
    ['ERP CSD transformed of all subjects averaged across epochs around peaks, n = ', ...
    num2str(count_peaks.all), ' epochs']);
title_string_z = strcat(...
    ['z-normalized ERP, CSD transformed of all subjects averaged across epochs around peaks, n = ', ...
    num2str(count_peaks.all), ' epochs']);
subtitle_string = strcat(['all channels']);


%% Plot ERPs as grand average and save plots


cd(peak_erp_plots_dir)

% Average ERP Plots
figure()
[~, min_ind] = plot_ERP(epoch_mean_all, mean_struct.time_vec, base_dur, title_string, subtitle_string);
savefig('all_condition_peaks_erp_csd')

cd(peak_erp_plots_dir)

% Average ERP Plots of z-score

figure()
[~, min_ind] = plot_ERP(z_epoch_mean_all, mean_struct.time_vec, base_dur, title_string_z, subtitle_string);
savefig('all_condition_peaks_erp_z_csd')


%% Get electrodes of largest deflection
% prepare single-channel plotting

% find minimal deflection
interesting_window = max(find(latencies_onset < 190 )); % find the topo window
% that is the most interesting in previous ERP plots. 
[~, min_chan_idx] = mink(topo_struct.all(:, interesting_window),3); % find the 
% three channels with the lowest defleciton
chan_lab_min = {ALLEEG(1).chanlocs(min_chan_idx).labels};

% find maximal defleciton
interesting_window = max(find(latencies_onset < 190 )); % find the topo window
% that is the most interesting in previous ERP plots. 
[~, max_chan_idx] = maxk(topo_struct.all(:, interesting_window),3); % find the 
% three channels with the lowest defleciton
chan_lab_max = {ALLEEG(1).chanlocs(max_chan_idx).labels};

% find maximal defleciton at 240-270
interesting_window = max(find(latencies_onset < 240 )); % find the topo window
% that is the most interesting in previous ERP plots. 
[~, max_chan_idx_240_270] = maxk(topo_struct.all(:, interesting_window),3); % find the 
% three channels with the lowest defleciton
chan_lab_max_240_270 = {ALLEEG(1).chanlocs(max_chan_idx_240_270).labels};

% Single Channel ERP Plot + Topo
% Plots the electrode that shows the biggest amplitude in the ERP

figure()
subplot(1, 2, 1)
title_string = strjoin(['electrode ', chan_lab_min, ', all subjects, all conditions around peaks.'])
subtitle_string = strjoin(['channel ', chan_lab_min]);
plot(mean_struct.time_vec, z_epoch_mean_all(min_chan_idx, :), 'LineWidth', 1.1)
title(title_string_z)
subtitle(subtitle_string)
hold on 
xline(0)
hold off
legend({chan_lab_min{:}, "peak"})
subplot(1, 2, 2)
topoplot(topo_struct.all(:, interesting_window), chan_locs, 'maplimits',  ...
    zlims, 'emarker2', {min_chan_idx,'s','m'}, 'electrodes','labelpoint')
title(strcat(['average of ', num2str(latencies_onset(interesting_window)), ' to ', ...
    num2str(latencies_offset(interesting_window)), ' ms']))
sgtitle(strjoin(['epochs around peaks all subjects all trials, elec ', chan_lab_min, ' are highlighted (3 min amp at 200ms win)']))

savefig('all_condition_peaks_min_180-210_topo_single_csd_z')


figure()
subplot(1, 2, 1)
title_string = strjoin(['electrode ', chan_lab_max, ', all subjects, all conditions around peaks.'])
subtitle_string = strjoin(['channel ', chan_lab_max]);
plot(mean_struct.time_vec, z_epoch_mean_all(max_chan_idx, :), 'LineWidth', 1.1)
title(title_string_z)
subtitle(subtitle_string)
hold on 
xline(0)
hold off
legend({chan_lab_max{:}, "peak"})
subplot(1, 2, 2)
topoplot(topo_struct.all(:, interesting_window), chan_locs, 'maplimits',  ...
    zlims, 'emarker2', {max_chan_idx,'s','m'}, 'electrodes','labelpoint')
title(strcat(['average of ', num2str(latencies_onset(interesting_window)), ' to ', ...
    num2str(latencies_offset(interesting_window)), ' ms']))
sgtitle(strjoin(['epochs around peaks all subjects all trials, elec ', chan_lab_max, ' are highlighted (3 min amp at 200ms win)']))

savefig('all_condition_peaks_max_180-210_topo_single_csd_z')


figure()
subplot(1, 2, 1)
title_string = strjoin(['electrode ', chan_lab_max_240_270, ', all subjects, all conditions around peaks.'])
subtitle_string = strjoin(['channel ', chan_lab_max_240_270]);
plot(mean_struct.time_vec, z_epoch_mean_all(max_chan_idx_240_270, :), 'LineWidth', 1.1)
title(title_string_z)
subtitle(subtitle_string)
hold on 
xline(0)
hold off
legend({chan_lab_max_240_270{:}, "peak"})
subplot(1, 2, 2)
topoplot(topo_struct.all(:, interesting_window), chan_locs, 'maplimits',  ...
    zlims, 'emarker2', {max_chan_idx_240_270,'s','m'}, 'electrodes','labelpoint')
colorbar()
title(strcat(['average of ', num2str(latencies_onset(interesting_window)), ' to ', ...
    num2str(latencies_offset(interesting_window)), ' ms']))
sgtitle(strjoin(['epochs around peaks all subjects all trials, elec ', chan_lab_max_240_270, ' are highlighted (3 min amp at 200ms win)']))

savefig('all_condition_peaks_max_240-270_topo_single_csd_z')

%% ERP image

% put all epochs of all subjects in one dataset
% initialize dataset
% TODO: Current Bugs: Subject 1 does not yet have a pursuit latency field
% TODO: Last Subject has no negative peaks (only S40) and no pursuit peaks
clear z_ALLEEG TMPEEG z_mean_struct mean_struct
sum(arrayfun(@(s) length(s.epoch), ALLEEG)) % how many epochs we need. 
TMPEEG = ALLEEG(2);

% get all data and all event structures
for s = 1:length(ALLEEG)-2
    % normalize each subjects' data
    z_tmpEEG = ALLEEG(s+1);
    % normalize data
    z_tmpEEG.data = normalize(ALLEEG(s+1).data,2, 'zscore');
    % get location where new data should go
    data_idx = size(TMPEEG.data, 3)+1:size(TMPEEG.data, 3)+size(z_tmpEEG.data, 3);
    % concatenate subject data in one dataset
    TMPEEG.data(:,:,data_idx) = z_tmpEEG.data;
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

TMPEEG.trials = length(TMPEEG.epoch)

% actual erp image parameters
chan_lab = 'F1';
chan_no = find(strcmp({TMPEEG.chanlocs.labels}, chan_lab));
plot_type = 1; % 1 = channel, 0 = component
project_channel = [[]];
smooth = round(length(TMPEEG.epoch)/15,0); % number of trials
dec_fac = smooth/2; % ratio of trials in to plot out. 
% if the smoothing width is larger than twice the decimation factor, no
% info is lost. 
title = strjoin(["Smoothing factor " smooth,", decimate " dec_fac ", chan " chan_lab]);
sort_type = {'S 40' 'S 50'};
sort_win = [-0.1 0.1];
sort_event_field = 'pursuit_lat';
align = Inf;
% 
[ALLEEG TMPEEG index] = eeg_store(ALLEEG, TMPEEG);
% eeglab redraw

% plotting of ERP image
figure; 
thingy = pop_erpimage(TMPEEG,plot_type, chan_no,[[]], title,smooth, dec_fac, ...
    sort_type,[],sort_event_field ,'yerplabel','\muV','erp','on', ...
    'cbar','on','align',Inf,'topo', { chan_no TMPEEG.chanlocs TMPEEG.chaninfo } );

