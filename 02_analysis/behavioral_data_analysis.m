%% Behavioral Data Analyses

% Author: Saskia Wilken
% Creation Date: 18.10.2022

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
track_data_dir = strjoin([output_dir, "03_parallelize_with_traj"], filesep);
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
track_data_dir = strjoin([output_dir, "03_parallelize_with_traj"], filesep);


%% Load Tracking Data

cd(track_data_dir)
load('all_tracking_data.mat', 'track_data') 

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

%% Plot Heatmap of Error

frame_rate = 250;
error_heatmap(track_data, frame_rate)



%% Prepare Scatterplot: Create one dataset with all the errors


clear behav_data behav_mean
cnt = 1;
for s = 1:length(ALLEEG)
    EEG = ALLEEG(s);
    for ep = 1:max([EEG.event.epoch])
        if isempty([EEG.event(round(median(find([EEG.event.epoch] == ep))) ...
                ).epoch_error]) | isempty([EEG.event(round(median(find([EEG.event.epoch] == ep))) ...
                ).pursuit_lat])
            behav_data(cnt,:) = nan(1,4);
        else
            behav_data(cnt,1) = [EEG.event(round(median(find([EEG.event.epoch] == ep)))).epoch_error] * 1080;
            behav_data(cnt,2) = [EEG.event(round(median(find([EEG.event.epoch] == ep)))).pursuit_lat];
            behav_data(cnt,3) = [EEG.event(round(median(find([EEG.event.epoch] == ep)))).epoch_error_z];
            behav_data(cnt,4) = [EEG.event(round(median(find([EEG.event.epoch] == ep)))).pursuit_lat_z];
        end
        cnt = cnt + 1;
    end
    one_subject_idx(s) = cnt;
    clear area_of_subj
    if s ~= 1
        area_of_subj = one_subject_idx(s-1):one_subject_idx(s)-1;
    else
        area_of_subj = 1:one_subject_idx(s)-1;
    end
    behav_mean(s, :) = mean(behav_data(area_of_subj,:), 1, 'omitnan');
end

cd(output_dir)
save('behav_data.mat', 'behav_data')
save('behav_mean.mat', 'behav_mean')
% DOn't run again, takes ages. Rather, load the saved dataset:


%% Plot scatterplot of pursuit latency and epoch error

cd(output_dir)
load('behav_data.mat', 'behav_data') % this is just the error and pursuit latency
behav_data = rmmissing(behav_data);

clear ALLEEG D EEG
% make scatterplot of the raw values
[corrs, pvals] = corrcoef(behav_data(:,1:2));

x = behav_data(:,1);
y = behav_data(:,2);
% R = 25;                        % circle radius
rand_idx = randi(numel(D), [1,100000]);
R = prctile(D(rand_idx), 1)
D = pdist2([x y],[x y]);        % create combinations of distances (eucledian)
D = pdist2([x y],[x y]);        % create combinations of distances
C = sum(D<R);                   % how mayn points inside circle
figure()
scatter(x,y,25,C,'fill')
colorbar
% plot(behav_data(:,1), behav_data(:,2), 'o', size);
title('scatter plot of pursuit latency and epoch error')
subtitle(['no significant correlation, rho = ',num2str(round(corrs(2),4)), ...
    ', p = ', num2str(round(pvals(2), 4)), ' ,' , num2str(size(behav_data,1)), ' observation pairs'])
xlabel('epoch error in px')
ylabel('pursuit latency in ms')
set(gca,'fontsize', 14);
% mean for each subject, scatterplot

clear ALLEEG D EEG

% make scatterplot of the z-transformed values
corrs = corrcoef(behav_data(:, 3:4));

x = behav_data(:,3);
y = behav_data(:,4);
% R = 25;                        % circle radius
rand_idx = randi(numel(D), [1,100000]);
R = prctile(D(rand_idx), 1)
D = pdist2([x y],[x y]);        % create combinations of distances (eucledian)
C = sum(D<R);                   % how mayn points inside circle
figure()
scatter(x,y,25,C,'fill')
colorbar
% plot(behav_data(:,1), behav_data(:,2), 'o', size);
title('scatter plot of pursuit latency and epoch error, z-transformed')
subtitle(['no significant correlation, ', num2str(size(behav_data,1)), ' observation pairs'])
xlabel('epoch error in px')
ylabel('pursuit latency in ms')
set(gca,'fontsize', 14);


%% Check what the pursuit latencies look like

% TODO: Put this in a separate script
% criteria:
% - pursuit peak with prominence = 0.05
% - a pursuit peak must follow the trajectory peak in the same epoch
% - pursuit and trajectory peaks must go in the same direction (up or down)
% - pursuit peak latency is between 80 and 300 ms. 
% - central peak of the epoch is not in the first half a second of a trial

disp(strjoin([round(no_purs/((purs_count+no_purs)/100),0), ...
    "% (", no_purs, " out of ", purs_count+no_purs, ...
    ") trajectory peaks are not followed by a pursuit peak."],""))
% at this time, 31472 out of 65627 trajectory peaks are not followed by a 
% pursuit peak (47.96%)
pursuit_lats = [];
all_pursuit_lats = [];
pursuit_lats_z = [];
all_pursuit_lats_z = [];
epoch_errors = [];
all_epoch_errors = [];
epoch_errors_z = [];
all_epoch_errors_z = [];

for s = 1:length(TMPEEG)
    % pursuit latencies raw
    pursuit_lats = [TMPEEG(s).event.pursuit_lat];
    % concatenate pursuit latencies of current subj with previous subj
    all_pursuit_lats = [all_pursuit_lats, pursuit_lats];
    % normalized pursuit latencies
    pursuit_lats_z = [TMPEEG(s).event.pursuit_lat_z];
    all_pursuit_lats_z = [all_pursuit_lats_z, pursuit_lats_z];
    % epoch error raw
    epoch_errors = [TMPEEG(s).event.epoch_error];
    all_epoch_errors = [all_epoch_errors, epoch_errors];
    % normalized epoch error
    epoch_errors_z = [TMPEEG(s).event.epoch_error_z];
    all_epoch_errors_z = [all_epoch_errors, epoch_errors_z];
end

% keep only one of the copies of the pursuit latency field values each. 
unique_lats_idx = find(diff(all_pursuit_lats));
% this method is not 100% failsafe - there is a chance that two following
% pursuit latencies are exactly the same and they would be missed by this
% algorithm. However, most pursuit lats should be in there. 
all_pursuit_lats_unique = all_pursuit_lats(unique_lats_idx);
unique_lats_idx_z = find(diff(all_pursuit_lats_z));
all_pursuit_lats_unique_z = all_pursuit_lats_z(unique_lats_idx_z);

BINWIDTH = 12;
range(all_pursuit_lats_unique)
PROMINENCE_THRESH = 0.05; 
ACCEPTABLE_LATS = [80, 300];
histogram(all_pursuit_lats_unique, ...
    round(range(all_pursuit_lats_unique)/BINWIDTH));
title(['Pursuit Latency'])
subtitle(strjoin(["Prominence thresh ", PROMINENCE_THRESH, ...
    ", peaks in same direction, latency lims between ", ...
    ACCEPTABLE_LATS, " ms", "binwidth = ", BINWIDTH, " ms"]))
mink(all_pursuit_lats_unique, 100)
% purs_count and length(all_pursuit_lats_unique) do not have the same size.
% which is odd. 

PROMINENCE_THRESH = 0.05; 
ACCEPTABLE_LATS = [80, 300];
histogram(all_pursuit_lats_unique_z, ...
    round(range(all_pursuit_lats_unique)/BINWIDTH));
title(['Pursuit Latency z-transformed per subject'])
subtitle(strjoin(["Prominence thresh ", PROMINENCE_THRESH, ...
    ", peaks in same direction, latency lims between ", ...
    ACCEPTABLE_LATS, " ms", "binwidth = ", BINWIDTH, " ms"]))
mink(all_pursuit_lats_unique, 100)

% check the epoch error distribution

figure()
BINWIDTH = 0.02;
histogram(all_epoch_errors, round(range(all_epoch_errors)/BINWIDTH));
title(['Epoch Error'])
subtitle(strjoin(["All epochs, all errors over time, binwidth = ", ...
    BINWIDTH*100, " % of screen"]))

figure()
histogram(all_epoch_errors_z, round(range(all_epoch_errors)/BINWIDTH));
title(['Epoch Error z-transformed'])
subtitle(strjoin(["All epochs, all errors over time, binwidth = ", ...
    BINWIDTH*100, " % of screen"]))

