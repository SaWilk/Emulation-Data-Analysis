%% Make ERP Plot

% Author: Saskia Wilken
% creation Date: 28.06.2022


%% Empty

format compact
format long G
clear
clc


%% Folders and Dependencies

% add path and start EEGlab toolbox
% addpath('R:\AG-Beste-Orga\Skripts\Toolbox\eeglab2021.0');
eeglab;
close;

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

% add functions in parent dir
addpath([char(parent_dir) filesep 'functions']);

% set input & output directory
input_dir = strjoin([parent_dir_2, "Emulation-Data-Output\03_parallelize_with_traj"], filesep);
output_dir = strjoin([parent_dir_2, "Emulation-Data-Output"], filesep);
output_dir_epoched = strjoin([parent_dir_2, "Emulation-Data-Output\04_epoched"], filesep);
output_dir_baseline = strjoin([parent_dir_2, "Emulation-Data-Output\05_baseline"], filesep);
output_dir_AR = strjoin([parent_dir_2, "Emulation-Data-Output\06_artifact_rejection"], filesep);

subdir_const_rand = strjoin([output_dir_epoched, "const_rand"], filesep);
subdir_occl = strjoin([output_dir_epoched filesep "occl_nonoccl"], filesep);
% out dirs for peaks for multiple epoch lengths
subdir_peaks = strjoin([output_dir_epoched filesep 'peaks_-500_750'], filesep);

study_dir = strjoin([output_dir, "study"], filesep);

savepath_baseline_const_rand = strjoin([output_dir_baseline, "const_rand"], filesep);
savepath_baseline_occl = strjoin([output_dir_baseline, 'occl_nonoccl'], filesep);
% out dirs for peaks for multiple epoch lengths
savepath_baseline_peaks = strjoin([output_dir_baseline, 'peaks_-500_750'], filesep);

output_dir_AR_const = strjoin([output_dir_AR, 'const'], filesep);
output_dir_AR_occl = strjoin([output_dir_AR, 'occl'], filesep);
% out dirs for peaks for multiple epoch lengths
output_dir_AR_peaks = strjoin([output_dir_AR, 'peaks_-500_750'], filesep);

mean_matrices_path = strjoin([output_dir, 'mean_matrices'], filesep);
mkdir(mean_matrices_path);
mean_matrices_peaks_epoched_path = strjoin([mean_matrices_path, 'peaks'], filesep);
mkdir(mean_matrices_peaks_epoched_path);

peak_erp_plots_dir = strjoin([parent_dir, "plots", "erp_plots", "peaks"], filesep);
mkdir(peak_erp_plots_dir);

% neighbors dir
neighbors_dir = strjoin([parent_dir_2, "Emulation-Data-Input", "EEG_files"], filesep);


%% Load Epoched data

eeglab;

% long epochs
cd(output_dir_AR_peaks);
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


%% Calculate Means Across Subjects

<<<<<<< Updated upstream

count_peaks = 0;
count_peaks_occ = 0;
count_peaks_vis = 0;
count_peaks_rand1 = 0;
count_peaks_const = 0;
count_peaks_rand2 = 0;
=======
count_peaks.all = 0;
count_peaks.occ = 0;
count_peaks.vis = 0;
count_peaks.rand1 = 0;
count_peaks.const = 0;
count_peaks.rand2 = 0;
>>>>>>> Stashed changes

for s = 1:size(ALLEEG, 2)

    EEG = ALLEEG(s);
    % all conditions
    mean_struct.all(:,:,s) = mean(EEG.data,3);
    % get number of epochs per subject for weighting
    z_mean_struct.num_epochs_all(s) = size(EEG.data, 3);
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
    mean_struct.occ_off(:,:,s) = mean(EEG.data(:,:,cond_ind(s).epoch_vis),3);
    mean_struct.occ_on(:,:,s) = mean(EEG.data(:,:,cond_ind(s).epoch_occ),3);
    mean_struct.rand1(:,:,s) = mean(EEG.data(:,:,cond_ind(s).epoch_rand1),3);
    mean_struct.const(:,:,s) = mean(EEG.data(:,:,cond_ind(s).epoch_const),3);
    mean_struct.rand2(:,:,s) = mean(EEG.data(:,:,cond_ind(s).epoch_rand2),3);

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

    % add peaks to know how much power each condition has
    count_peaks = size(EEG.data, 3) + count_peaks;
    count_peaks_occ = size(EEG.data(:,:,cond_ind(s).epoch_occ), 3) + count_peaks_occ;
    count_peaks_vis = size(EEG.data(:,:,cond_ind(s).epoch_vis), 3) + count_peaks_vis;
    count_peaks_rand1 = size(EEG.data(:,:,cond_ind(s).epoch_rand1), 3) + count_peaks_rand1;
    count_peaks_const = size(EEG.data(:,:,cond_ind(s).epoch_const), 3) + count_peaks_const;
    count_peaks_rand2 = size(EEG.data(:,:,cond_ind(s).epoch_rand2), 3) + count_peaks_rand2;

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
save("mean_struct.mat", 'mean_struct')
save("condition_indices.mat", 'cond_ind')
<<<<<<< Updated upstream
=======
save("count_peaks.mat", 'count_peaks')
save("z_mean_struct.mat", 'z_mean_struct')

>>>>>>> Stashed changes


%% Prepare for ERP Plotting

cd(mean_matrices_peaks_epoched_path)

load('mean_struct.mat')

avg_peak_dist = 300;
% average across subjects
epoch_mean_all = mean(mean_struct.all,3);
epoch_mean_occ = mean(mean_struct.occ_on,3, 'omitnan'); % omitnan is
% % necessary because subject 18 has no task B
epoch_mean_vis = mean(mean_struct.occ_off,3);
epoch_mean_rand1 = mean(mean_struct.rand1,3);
epoch_mean_const = mean(mean_struct.const,3);
epoch_mean_rand2 = mean(mean_struct.rand2,3);

load('z_mean_struct.mat')

% get z-weighted mean for all subjects
z_epoch_mean_all = mean(z_mean_struct.all,3);
z_epoch_mean_occ_on = mean(z_mean_struct.occ_on,3, 'omitnan');
z_epoch_mean_occ_off = mean(z_mean_struct.occ_off,3, 'omitnan');
z_epoch_mean_rand1 = mean(z_mean_struct.rand1,3);
z_epoch_mean_const = mean(z_mean_struct.const,3);
z_epoch_mean_rand2 = mean(z_mean_struct.rand2,3);


base_dur = 500;

title_string = strcat(...
    ['ERP of all subjects averaged across epochs around peaks, n = ', ...
    num2str(count_peaks.all), ' epochs']);
title_string_z = strcat(...
    ['z-normalized ERP of all subjects averaged across epochs around peaks, n = ', ...
    num2str(count_peaks.all), ' epochs']);
subtitle_string = strcat(['all channels']);


%% Plot ERPs as grand average and save plots


cd(peak_erp_plots_dir)

% Average ERP Plots
<<<<<<< Updated upstream
title_string = strcat(...
    ['ERP of all subjects averaged across epochs around peaks, n = ', ...
    num2str(count_peaks), ' epochs']);
subtitle_string = strcat(['all channels']);
=======

>>>>>>> Stashed changes
figure()
[~, min_ind] = plot_ERP(epoch_mean_all, mean_struct.time_vec, base_dur, title_string, subtitle_string);
savefig('all_condition_peaks_erp')

cd(peak_erp_plots_dir)

% Average ERP Plots of z-score

figure()
[~, min_ind] = plot_ERP(z_epoch_mean_all, mean_struct.time_vec, base_dur, title_string_z, subtitle_string);
savefig('all_condition_peaks_erp_z')

% Single Channel ERP
% Plots the electrode that shows the biggest amplitude in the ERP

min_chan_idx = min_ind(1);
chan_lab = EEG.chanlocs(min_chan_idx).labels;
subtitle_string = strcat(['channel ', chan_lab]);
figure()
plot_ERP(epoch_mean_all(min_chan_idx, :), mean_struct.time_vec, base_dur, title_string, subtitle_string)
% TODO: Instead plot the electrode with the largest t-value in the CBP test

% Average ERP Plots occluded vs non-occluded

figure()
subplot(2, 1, 1)
<<<<<<< Updated upstream
title_string = strcat(...
    ['ERP of all subjects averaged across epochs around peaks, n = ', ...
    num2str(count_peaks_vis), ' epochs']);
=======
>>>>>>> Stashed changes
subtitle_string = strcat(['all channels, only non-occluded']);
[~, min_ind_vis] = plot_ERP(epoch_mean_vis, mean_struct.time_vec, base_dur, title_string, subtitle_string);
ax1 = gca;

subplot(2, 1, 2)
<<<<<<< Updated upstream
title_string = strcat(...
    ['ERP of all subjects averaged across epochs around peaks, n = ', ...
    num2str(count_peaks_occ), ' epochs']);
=======
>>>>>>> Stashed changes
subtitle_string = strcat(['all channels, only occluded']);
[~, min_ind_occ] = plot_ERP(epoch_mean_occ, mean_struct.time_vec, base_dur, title_string, subtitle_string);
ax2 = gca;
linkaxes([ax1 ax2],'xy')

savefig('occlusion_vs_visible_peaks_erp')


% Average ERP Plots of z-score occluded vs non-occluded

figure()
subplot(2, 1, 1)
title_string_z = strcat(...
    ['z-normalized ERP of all subjects averaged across epochs around peaks, n = ', ...
    num2str(count_peaks.vis), ' epochs']);
subtitle_string = strcat(['all channels, only non-occluded']);
[~, min_ind_vis] = plot_ERP(z_epoch_mean_occ_off, mean_struct.time_vec, base_dur, title_string_z, subtitle_string);
ax1 = gca;

subplot(2, 1, 2)
title_string_z = strcat(...
    ['z-normalized ERP of all subjects averaged across epochs around peaks, n = ', ...
    num2str(count_peaks.occ), ' epochs']);
subtitle_string = strcat(['all channels, only occluded']);
[~, min_ind_occ] = plot_ERP(z_epoch_mean_occ_on, mean_struct.time_vec, base_dur, title_string_z, subtitle_string);
ax2 = gca;
linkaxes([ax1 ax2],'xy')

savefig('occlusion_vs_visible_peaks_erp_z')


% Average ERP Plots constant vs random

figure()
subplot(3, 1, 1)
title_string = strcat(...
    ['ERP of all subjects averaged across epochs around peaks, n = ', ...
    num2str(count_peaks_const), ' epochs']);
subtitle_string = strcat(['all channels, only constant']);
[~, min_ind_const] = plot_ERP(epoch_mean_const, mean_struct.time_vec, base_dur, title_string, subtitle_string);
ax1 = gca;
subplot(3, 1, 2)
title_string = strcat(...
    ['ERP of all subjects averaged across epochs around peaks, n = ', ...
    num2str(count_peaks_rand1), ' epochs']);
subtitle_string = strcat(['all channels, only random1']);
[~, min_ind_rand1] = plot_ERP(epoch_mean_rand1, mean_struct.time_vec, base_dur, title_string, subtitle_string);
ax2 = gca;

subplot(3, 1, 3)
title_string = strcat(...
    ['ERP of all subjects averaged across epochs around peaks, n = ', ...
    num2str(count_peaks_rand2), ' epochs']);
subtitle_string = strcat(['all channels, only random2']);
[~, min_ind_rand2] = plot_ERP(epoch_mean_rand2, mean_struct.time_vec, base_dur, title_string, subtitle_string);
ax3 = gca;

linkaxes([ax1 ax2 ax3],'xy')
savefig('constant_vs_random_peaks_erp')


%% Plot Contrasts of ERPS as simple subtractions and save plots

% seems pretty useless. 

% Average ERP Plots
figure()
title_string = strcat(...
    ['ERP-difference of all subjects averaged across epochs around peaks']);
subtitle_string = strcat(['all channels, non-occluded - occluded']);
figure()
[~, min_ind] = plot_ERP(diff_occ, mean_struct.time_vec, base_dur, title_string, subtitle_string);
savefig('diff_occ_vis')

% Average ERP Plots
figure()
subplot(2, 1, 1)
title_string = strcat(...
    ['ERP-difference of all subjects averaged across epochs around peaks']);
subtitle_string = strcat(['all channels, constat - random1']);
[~, min_ind] = plot_ERP(diff_const_rand1, mean_struct.time_vec, base_dur, title_string, subtitle_string);
ax1 = gca;
linkaxes([ax1 ax2],'xy')

subplot(2, 1, 2)
title_string = strcat(...
    ['ERP-difference of all subjects averaged across epochs around peaks']);
subtitle_string = strcat(['all channels, constat - random2']);
[~, min_ind] = plot_ERP(diff_const_rand2, mean_struct.time_vec, base_dur, title_string, subtitle_string);
ax2 = gca;
linkaxes([ax1 ax2],'xy')
savefig('diff_const_rand')


%% Prepare for Topoplotting

cd(mean_matrices_peaks_epoched_path)

load('mean_struct.mat')

avg_peak_dist = 300;
% average across subjects
epoch_mean_all = mean(mean_struct.all,3);
epoch_mean_occ = mean(mean_struct.occ_on,3, 'omitnan'); % omitnan is
% % necessary because subject 18 has no task B
epoch_mean_vis = mean(mean_struct.occ_off,3);
epoch_mean_rand1 = mean(mean_struct.rand1,3);
epoch_mean_const = mean(mean_struct.const,3);
epoch_mean_rand2 = mean(mean_struct.rand2,3);

base_dur = 500;
num_topos = 15;
latencies_onset = linspace(0, 700, num_topos);
latencies_offset = latencies_onset + 50;

% extract means across 50 ms periods for topoplot
for i = 1:length(latencies_onset)

    [~, lat_idx_onset] = min(abs(mean_struct.time_vec - latencies_onset(i)));
    [~, lat_idx_offset] = min(abs(mean_struct.time_vec - latencies_offset(i)));
    topo_struct.all(:, i) = mean(epoch_mean_all(:,lat_idx_onset:lat_idx_offset), 2);
    topo_struct.vis(:, i) = mean(epoch_mean_vis(:,lat_idx_onset:lat_idx_offset), 2);
    topo_struct.occ(:, i) = mean(epoch_mean_occ(:,lat_idx_onset:lat_idx_offset), 2);
    topo_struct.const(:, i) = mean(epoch_mean_const(:,lat_idx_onset:lat_idx_offset), 2);
    topo_struct.rand1(:, i) = mean(epoch_mean_rand1(:,lat_idx_onset:lat_idx_offset), 2);
    topo_struct.rand2(:, i) = mean(epoch_mean_rand2(:,lat_idx_onset:lat_idx_offset), 2);
    topo_struct.diff_occ(:, i) = mean(mean_struct.diff_occ(:,lat_idx_onset:lat_idx_offset), 2);
    topo_struct.diff_const_rand1(:, i) = mean(mean_struct.diff_const_rand1(:,lat_idx_onset:lat_idx_offset), 2);
    topo_struct.diff_const_rand2(:, i) = mean(mean_struct.diff_const_rand2(:,lat_idx_onset:lat_idx_offset), 2);

end

chan_locs = ALLEEG(1).chanlocs


%% Create Topoplot of the time points of interest

cd(peak_erp_plots_dir)

% doc topoplot: https://rdrr.io/cran/erpR/man/topoplot.html

% all conditions
figure()
zlims = [min(topo_struct.all, [], 'all'), max(topo_struct.all, [], 'all')]';
for i = 1:length(latencies_onset)

    subplot(num_topos/5, num_topos/3, i)
    topoplot(topo_struct.all(:, i), chan_locs, 'maplimits',  zlims)
    title(strcat(['average of ', num2str(latencies_onset(i)), ' to ', ...
        num2str(latencies_offset(i)), ' ms']))
    colorbar()
end
sgtitle('epochs around peaks all subjects all trials')

savefig('all_condition_peaks_topo')


% occlusion vs non-occlusion
figure();
zlims = [min([topo_struct.occ,topo_struct.vis] , [], 'all'), max([topo_struct.occ,topo_struct.vis], [], 'all')]';
rows = 5;
cols = 7;
sub = [1:3, 8:10, 15:17, 22:24, 29:31];
sub2 = [5:7, 12:14, 19:21, 26:28, 33:35];
subvec = [sub; sub2];
subvec = subvec(:);
pos = 0;
for i = 1:length(latencies_onset)
    pos =  1 + pos;
    subplot(rows, cols, subvec(pos));
    topoplot(topo_struct.occ(:, i), chan_locs, 'maplimits',  zlims);
    title('occluded');
    subtitle(strcat(['average of ', num2str(latencies_onset(i)), ' to ', ...
        num2str(latencies_offset(i)), ' ms']));
    colorbar();
    pos =  1 + pos;
    subplot(rows, cols, subvec(pos));
    topoplot(topo_struct.vis(:, i), chan_locs, 'maplimits',  zlims);
    title('visible');
    subtitle(strcat(['average of ', num2str(latencies_onset(i)), ' to ', ...
        num2str(latencies_offset(i)), ' ms']));
    colorbar();
end
sgtitle("peak epoch ERPs: occlusion (left) vs non-occlusion (right)");
savefig('occclusion_vs_visible_peaks_topo');


% constant vs random
figure();
zlims = [min([topo_struct.const,topo_struct.rand1,topo_struct.rand2] , [], 'all'), ...
    max([topo_struct.const,topo_struct.rand1,topo_struct.rand2], [], 'all')]';
rows = 5;
cols = 11;
sub = [1:3, 12:14, 23:25, 34:36, 45:47];
sub2 = [5:7, 16:18, 27:29, 38:40, 49:51];
sub3 = [9:11, 20:22, 31:33, 42:44, 53:55];
subvec = [sub; sub2; sub3];
subvec = subvec(:);
pos = 0;
for i = 1:length(latencies_onset)
    pos =  1 + pos;
    subplot(rows, cols, subvec(pos));
    topoplot(topo_struct.const(:, i), chan_locs, 'maplimits',  zlims);
    title('constant');
    subtitle(strcat(['average of ', num2str(latencies_onset(i)), ' to ', ...
        num2str(latencies_offset(i)), ' ms']));
    colorbar();
    pos =  1 + pos;
    subplot(rows, cols, subvec(pos));
    topoplot(topo_struct.rand1(:, i), chan_locs, 'maplimits',  zlims);
    title('random1');
    subtitle(strcat(['average of ', num2str(latencies_onset(i)), ' to ', ...
        num2str(latencies_offset(i)), ' ms']));
    colorbar();
    pos =  1 + pos;
    subplot(rows, cols, subvec(pos));
    topoplot(topo_struct.rand2(:, i), chan_locs, 'maplimits',  zlims);
    title('random2');
    subtitle(strcat(['average of ', num2str(latencies_onset(i)), ' to ', ...
        num2str(latencies_offset(i)), ' ms']));
    colorbar();
end
sgtitle("peak epochs: constant (left) vs random1 (middle) vs random2 (right)");
savefig('constant_vs_random_peaks_topo');


% difference of occ vs vis
figure()
zlims = [min(topo_struct.diff_occ, [], 'all'), max(topo_struct.diff_occ, [], 'all')]';
for i = 1:length(latencies_onset)

    subplot(num_topos/5, num_topos/3, i)
    topoplot(topo_struct.diff_occ(:, i), chan_locs, 'maplimits',  zlims);
    title(strcat(['average of ', num2str(latencies_onset(i)), ' to ', ...
        num2str(latencies_offset(i)), ' ms']))
    colorbar()
end
sgtitle('peak epoch ERP differences: visible - occluded')

savefig('diff_occ_vis_peaks_topo')


%% Load Study

[ STUDY ALLEEG ] = pop_loadstudy('filename', 'study_peaks_epoched.study', 'filepath', study_dir)


%% Perform cluster-Based Permutation Test TODO Attempt with EEGLAB

% https://www.fieldtriptoolbox.org/tutorial/cluster_permutation_timelock/
% not sure whether I am using the fieldtrip parameters correctly or whether
% I am supposed to specify tzhings like alpha in eeglab or fieltrip syntax.
% try it out.
% [STUDY neighbors] = std_prepare_neighbors( STUDY, ALLEEG)
load(strjoin([neighbors_dir, 'fieldtrip_EEG_KJP_neighbours_61.mat'], filesep))
% stat_cond
% TODO: Add a condition in which all data is zero by copying the data and

eeglab redraw

% Load Previously Generated Mean Structures

cd(mean_matrices_peaks_epoched_path)

load('mean_struct.mat')

minimal_latency = mean_struct.time_vec(min_ind(2)); % the point in which
% % the signal is as negative as possible

peak_eeg = mean_struct.all;
zero_eeg = zeros(size(peak_eeg));

peak_across_subj = mean(peak_eeg, 3);
peak_across_subj_min_lat = peak_across_subj(:,minimal_latency);


% create data for a zero condition

% for subj = 1:length(ALLEEG)
%     zero_eeg(subj).("data") = zeros(size(ALLEEG(subj).data));
% end


For example, to compute mean ERPs statistics from a
STUDY for epochs of 800 frames in two conditions from three
groups of 12 subjects:

>> erp_data_for_all_subjects_possibly_from_STUDY = { [800x12] [800x12] [800x12];... % 3 groups, cond 1
    [800x12] [800x12] [800x12] };  % 3 groups, cond 2
[pcond, pgroup, pinter, statscond, statsgroup, statsinter] = std_stat(...
    {peak_eeg; zero_eeg}, 'groupstats', 'off', ...
    'condstats', 'off', 'method', 'permutation', 'naccu', '1000', 'alpha', '0.9', ...
    'mcorrect', 'none', 'mode', 'fieldtrip', 'fieldtripnaccu', 'numrandomization', ...
    'fieldtripalpha', '0.9', 'fieldtripmethod', 'montecarlo', ...
    'fieldtripmcorrect','cluster', 'fieldtripclusterparam', {'clusterstatistic', ...
    'maxsum', 'statistic', 'indepsamplesT'}, 'fieldtripchannelneighbor', 'neighbours' )
% 'fieldtripclusterparam' string or cell array for optional parameters
%                              for cluster correction method, see function
%                              ft_statistics_montecarlo for more information.
%    'fieldtripchannelneighbor' - Fieldtrip channel neighbour structure for
%                                 cluster correction method, see function
%                                 std_prepare_neighbors for more information.
% use this to find the electrode cluster
% In EEGLAB 12, press the statistics button in the channel or component STUDY
% plotting interface. Then you can select cluster statistics (click on
% Fieldtrip, then select monte-carlo statistics, then select cluster correction).
% It works as the other STUDY statistics.

%% Cluster Based Permutation Test Attempt with Fieldtrip

% Script by Adriana Böttcher
% Adapted by Saskia Wilken


%% add fieldtrip path
% eeglab;
% close;

% file location
addpath 'C:\Users\swilk\AppData\Roaming\MathWorks\MATLAB Add-Ons\Collections\FieldTrip'
file_path = matlab.desktop.editor.getActiveFilename;
file_path = fileparts(file_path);
addpath(file_path);
filepath_parts = strsplit(file_path, filesep);
parent_dir = filepath_parts(1:end-1);
parent_dir = strjoin(parent_dir, filesep);
parent_dir_2 = filepath_parts(1:end-2);
parent_dir_2 = strjoin(parent_dir_2, filesep);
output_dir_AR = strjoin([parent_dir_2, "Emulation-Data-Output\06_artifact_rejection"], filesep);
output_dir_AR_peaks = strjoin([output_dir_AR, 'peaks_-500_750'], filesep);
input_path = strjoin([output_dir_AR, 'peaks_-500_750'], filesep);
addpath([char(parent_dir) filesep 'functions']);

ft_defaults;

% load('fieldtrip_EEG_KJP_elec_61');
% load('fieldtrip_EEG_KJP_neighbours_61.mat');
% elec.label = upper(elec.label);


%% Prepare Cluster Based Permutaiton Testing


cd(mean_matrices_peaks_epoched_path)
load("condition_indices.mat")

% convert to fieldtrip:
% all data, compare with zero
for s = 1:length(ALLEEG)
    all_eeg_field(s) = eeglab2fieldtrip(ALLEEG(s), 'timelockanalysis', 'none');
end
zero_field = all_eeg_field;
for i = 1:length(all_eeg_field)
    zero_field(i).avg = zeros(size(zero_field(i).avg));
    zero_field(i).var = zeros(size(zero_field(i).avg));
end
for i = 1:length(all_eeg_field)
    eeg_cell{i} = all_eeg_field(i);
    zero_cell{i} = zero_field(i);
end

% Comparing Occlusion and Visible
for s = 1:length(ALLEEG)
    OCCEEG = ALLEEG(s);
    OCCEEG.data = OCCEEG.data(:,:,cond_ind(s).epoch_occ);
    % TODO: FOR SOME STRANGE REASON STOPS AT s = 2
    occ_eeg_field(i) = eeglab2fieldtrip(OCCEEG, 'timelockanalysis', 'none');
end

% I NEED THIS: ALTERNATIVELY, I CAN GENERATE IT MYSELF
% load('fieldtrip_EEG_KJP_elec_61');
% load('fieldtrip_EEG_KJP_neighbours_61.mat');
% this data is missing, Adriana will look it up late.r
% elec.label = upper(elec.label);

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

design = zeros(2, size(eeg_cell, 2)*2);
for i = 1:size(eeg_cell, 2)
    design(1,i) = i;
end
for i = 1:size(eeg_cell, 2)
    design(1,size(eeg_cell, 2)+i) = i;
end
design(2,1:size(eeg_cell, 2))        = 1;
design(2,size(eeg_cell, 2)+1:2*size(eeg_cell, 2)) = 2;


%% Get Neighbours

% generates a file called 'neighbours'
load(strjoin([neighbors_dir, 'fieldtrip_EEG_KJP_neighbours_61.mat'], filesep))

cfg_neighb        = [];
cfg_neighb.method = 'distance';
neighbours        = neighbours % ft_prepare_neighbours(cfg_neighb, dataFC_LP);


%%
% CLUSTER TEST
% erp vs zero
calpha  = 0.001;
alpha   = 0.001;
latency = [0.18, 0.22];

% cfg is the configuraiton structure of fieldtrip
cfg                     = [];
cfg.design              = design;
cfg.uvar                = 1;
cfg.ivar                = 2;
cfg.channel             = {'all'};
cfg.avgovertime         = 'no';
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
cfg.latency             = latency; % time interval over which the experimental
% conditions must be compared (in seconds)

stats = ft_timelockstatistics(cfg, eeg_cell{:}, zero_cell{:});
% The field prob contains the proportion of draws from the permutation
% distribution with a maximum cluster-level statistic that is larger than
% the cluster-level test statistic
stats.prob
stats.stat %  cluster-level test statistic (here with maxsum: the sum of
% the T-values in this cluster).

%% Plot the CBP results Nicos Approach
% https://github.com/fieldtrip/fieldtrip/blob/master/ft_clusterplot.m
% https://www.fieldtriptoolbox.org/tutorial/plotting/
% https://github.com/fieldtrip/fieldtrip/blob/release/ft_topoplotER.m


% ft_prepare_layout
load(strjoin([neighbors_dir, 'fieldtrip_EEG_KJP_layout_61.mat'], filesep))

cfg.layout = lay;
% cfg.style = 'blank';
% cfg.contournum = 0;
cfg.highlightcolorpos         = [0 0 0.75];
cfg.highlightcolorneg         = [0.75 0 0];
% cfg.highlightsymbolseries     = ['x', 'x', 'x', 'x', 'x'];
cfg.subplotsize    = [4, 4];
% cfg.saveaspng = "cluster_GO_19_pre.png";

figure(2)
ft_clusterplot(cfg, stats);
sgtitle(strjoin(["Significant clusters, calpha = ", calpha, " alpha = ", alpha], ""));
colorbar()

% topoplot verwenden um t-werte der channels zu plotten
% 'mask' um nur überschwellige Werte anzuzeigen
% Bedingungsvergleiche testen


%% Plot the CBP Results website approach

% To plot the results of the permutation test, we use the plotting function
% ft_topoplotER. In doing so, we will plot a topography of the difference
% between the two experimental conditions (FIC and FC). On top of that, and
% for each timestep of interest, we will highlight the sensors which are
% members of significant clusters. First, however, we must calculate the
% difference between conditions using ft_math.

cfg    = [];
avgFIC = ft_timelockanalysis(cfg, dataFIC_LP);
avgFC  = ft_timelockanalysis(cfg, dataFC_LP);

% Then take the difference of the averages using ft_math
cfg           = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
raweffectFICvsFC = ft_math(cfg, avgFIC, avgFC);

% We then construct a boolean matrix indicating whether a channel/time point
% belongs to a cluster that we deem interesting to inspect. This matrix has
% size [Number_of_MEG_channels x Number_of_time_samples], like
% stat.posclusterslabelmat. We’ll make two such matrices: one for positive
% clusters (named pos), and one for negative (neg). All (channel,time)-pairs
% belonging to the large clusters whose probability of occurrence is
% sufficiently low in relation to the associated randomization distribution
% of clusterstats will be coded in the new boolean matrix as 1, and all those
% that don’t will be coded as 0.

% Make a vector of all p-values associated with the clusters from ft_timelockstatistics.
pos_cluster_pvals = [stat.posclusters(:).prob];

% Then, find which clusters are deemed interesting to visualize, here we use a cutoff criterion based on the
% cluster-associated p-value, and take a 5% two-sided cutoff (i.e. 0.025 for the positive and negative clusters,
% respectively
pos_clust = find(pos_cluster_pvals < 0.025);
pos       = ismember(stat.posclusterslabelmat, pos_clust);

% and now for the negative clusters...
neg_cluster_pvals = [stat.negclusters(:).prob];
neg_clust         = find(neg_cluster_pvals < 0.025);
neg               = ismember(stat.negclusterslabelmat, neg_clust);

% Alternatively, we can manually select which clusters we want to plot. If
% we only want to see the extent of the first (i.e. most significant)
% positive and negative clusters, for instance, we can do so as follows:

pos = stat.posclusterslabelmat == 1; % or == 2, or 3, etc.
neg = stat.negclusterslabelmat == 1;

% To plot a sequence of twenty topographic plots equally spaced between 0 and
% 1 second, we define the vector j of time steps. These time intervals
% correspond to the samples m in stat and in the variables pos. and neg. m
% and j must, therefore, have the same length.

% To be sure that your sample-based time windows align with your time windows
% in seconds, check the following:

timestep      = 0.05; % timestep between time windows for each subplot (in seconds)
sampling_rate = dataFC_LP.fsample; % Data has a temporal resolution of 300 Hz
sample_count  = length(stat.time);
% number of temporal samples in the statistics object
j = [0:timestep:1]; % Temporal endpoints (in seconds) of the ERP average computed in each subplot
m = [1:timestep*sampling_rate:sample_count]; % temporal endpoints in M/EEG samples
To plot the data use the following for-loop:

% First ensure the channels to have the same order in the average and in the statistical output.
% This might not be the case, because ft_math might shuffle the order
[i1,i2] = match_str(raweffectFICvsFC.label, stat.label);

for k = 1:20
    subplot(4,5,k);
    cfg = [];
    cfg.xlim = [j(k) j(k+1)];   % time interval of the subplot
    cfg.zlim = [-2.5e-13 2.5e-13];
    % If a channel is in a to-be-plotted cluster, then
    % the element of pos_int with an index equal to that channel
    % number will be set to 1 (otherwise 0).

    % Next, check which channels are in the clusters over the
    % entire time interval of interest.
    pos_int = zeros(numel(raweffectFICvsFC.label),1);
    neg_int = zeros(numel(raweffectFICvsFC.label),1);
    pos_int(i1) = all(pos(i2, m(k):m(k+1)), 2);
    neg_int(i1) = all(neg(i2, m(k):m(k+1)), 2);

    cfg.highlight   = 'on';
    % Get the index of the to-be-highlighted channel
    cfg.highlightchannel = find(pos_int | neg_int);
    cfg.comment     = 'xlim';
    cfg.commentpos  = 'title';
    cfg.layout      = 'CTF151_helmet.mat';
    cfg.interactive = 'no';
    cfg.figure      = 'gca'; % plots in the current axes, here in a subplot
    ft_topoplotER(cfg, raweffectFICvsFC);
end
% In this for-loop, cfg.xlim defines the time interval of each subplot. The
% variables pos_int and neg_int boolean vectors indicating which channels of
% pos and neg are significant in the time interval of interest. This is
% defined in cfg.highlightchannel. The for-loop plots 20 subplots covering a
% time interval of 50 ms each. Running this for-loop creates the following figure:

%%
% CLUSTER TEST
%constant vs. rand2
% calpha  = 0.05;
% alpha   = 0.05;
%
% cfg                     = [];
% cfg.design              = design;
% cfg.uvar                = 1;
% cfg.ivar                = 2;
% cfg.channel             = {'all'};
% cfg.frequency           = [4 7];
% cfg.avgoverfreq         = 'yes';
% cfg.latency             = [min(constant{1}.time) max(constant{1}.time)]; %??
% cfg.avgovertime         = 'yes';
% cfg.method              = 'montecarlo';
% cfg.statistic           = 'depsamplesT';
% cfg.correctm            = 'cluster';
% cfg.clusteralpha        = calpha;               % 0.1;
% cfg.clusterstatistic    = 'maxsum';
% cfg.minnbchan           = 0;
% cfg.neighbours          = neighbours;
% cfg.tail                = -1;
% cfg.clustertail         = -1;
% cfg.alpha               = alpha;               % 0.025;
% cfg.numrandomization    = 1000;
%
% stats      = ft_freqstatistics(cfg, constant{:}, rand2{:} );


%% Create ERP Plot of the data TODO


% TODO: sort by error size in epoch
% sort by
'[STUDY ALLEEG] = std_editset( STUDY, ALLEEG, 'name','emulation_study_old','updatedat','on','rmclust','off' );
[STUDY ALLEEG] = std_checkset(STUDY, ALLEEG);'
'[STUDY ALLEEG] = std_precomp(STUDY, ALLEEG, {},'savetrials','on','interp','on','recompute','on','erp','on','erpparams',{'rmbase' [-500 0] },'erpim','on','erpimparams',{'nlines' 10 'smoothing' 10});['' ...
    ''][STUDY EEG] = pop_savestudy( STUDY, EEG, 'filename','study_peaks_epoched.study','filepath','C:\\wilken\\Emulation-Data-Output\\study\\')
std_erspplot(STUDY, ALLEEG)


%% Do TF Decomposition on data


%% Surface Laplacians CSD toolbox (aka current source density)

% keep default settings
addpath(genpath('C:\wilken\eeglab2019_0\CSDtoolbox'))

% tutorial: https://psychophysiology.cpmc.columbia.edu/Software/CSDtoolbox/tutorial.html
% common errors: https://psychophysiology.cpmc.columbia.edu/Software/CSDtoolbox/errors.html
% FAQ: https://psychophysiology.cpmc.columbia.edu/Software/CSDtoolbox/FAQ.html

% Get usable list of electrodes from EEGlab data structure
for site = 1:60
    trodes{site}=(ALLEEG(1).chanlocs(site).labels);
end;
trodes=trodes';

% Get Montage for use with CSD Toolbox
Montage_64=ExtractMontage('10-5-System_Mastoids_EGI129.csd',trodes);
MapMontage(Montage_64);

%% Derive G and H!
[G,H] = GetGH(Montage_64);

%% Save G and H to later import when doing the CSD transform on files
% save('G:\PhysioData\MN_Fear\G.mat', 'G');
% save('G:\PhysioData\MN_Fear\H.mat', 'H');

% revised method to store G and H matrices with CSD montage for later import
Montage = Montage_64;                             % use generic variable name
save G:\PhysioData\CSDmontage_64.mat G H Montage; % save variables to Matlab file
clear G H Montage;                                % remove variables from workspace
load G:\PhysioData\CSDmontage_64.mat;             % restore variables to the workspace

%% SAMPLE TO OPEN BIOSEMI -- can use instead of neuroscan version above
% [FileName,PathName,FilterIndex] = uigetfile('L:\Physiodata\*.bdf','Choose File to Open')
% FullName = [PathName FileName];
% [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
% EEG = pop_biosig(FullName, 'ref',[65 66] ,'blockepoch','off'); % Choose ref sites per your montage
% [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'gui','off');
% eeglab redraw;

tic                                        % stopwatch on
for ne = 1:length(ALLEEG(1).epoch)               % loop through all epochs
    myEEG = single(ALLEEG(1).data(:,:,ne));      % reduce data precision to reduce memory demand
    MyResults = CSD(myEEG,G,H);            % compute CSD for <channels-by-samples> 2-D epoch
    data(:,:,ne) = MyResults;              % assign data output
end
looping_CSD_final = double(data);          % final CSD data
looping_time = toc                         % stopwatch off

data = NaN;                         % re-initialize data output


%% Stuff I tried out.

STUDY = std_erpplot(STUDY,ALLEEG)

% figure;
% pop_plottopo(EEG, [1:60] , ''5STJS'', 0, ''ydir'',1);
%
% pop_topoplot(EEG, 1, [-500 748] ,'5STJS',[1 2] ,0,'electrodes','on')
%
% pop_erpplot()
%
%
%
% https://www.fieldtriptoolbox.org/tutorial/cluster_permutation_timelock/


