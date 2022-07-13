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
output_dir = strjoin([parent_dir_2, "Emulation-Data-Output\"], filesep);
output_dir_epoched = strjoin([parent_dir_2, "Emulation-Data-Output\04_epoched"], filesep);
output_dir_baseline = strjoin([parent_dir_2, "Emulation-Data-Output\05_baseline"], filesep);
output_dir_AR = strjoin([parent_dir_2, "Emulation-Data-Output\06_artifact_rejection"], filesep);

subdir_const_rand = strjoin([output_dir_epoched, "const_rand"], filesep);
subdir_occl = strjoin([output_dir_epoched filesep "occl_nonoccl"], filesep);
% out dirs for peaks for multiple epoch lengths
subdir_peaks_medium = strjoin([output_dir_epoched filesep 'peaks_-750_500'], filesep);
subdir_peaks_long = strjoin([output_dir_epoched filesep 'peaks_-1000_750'], filesep);
subdir_peaks_short = strjoin([output_dir_epoched filesep 'peaks_-500_250'], filesep);


savepath_baseline_const_rand = strjoin([output_dir_baseline, "const_rand"], filesep);
savepath_baseline_occl = strjoin([output_dir_baseline, 'occl_nonoccl'], filesep);
% out dirs for peaks for multiple epoch lengths
savepath_baseline_peaks_medium = strjoin([output_dir_baseline, 'peaks_-750_500'], filesep);
savepath_baseline_peaks_short = strjoin([output_dir_baseline, 'peaks_-500_250'], filesep);
savepath_baseline_peaks_long = strjoin([output_dir_baseline, 'peaks_-1000_750'], filesep);

output_dir_AR_const = strjoin([output_dir_AR, 'const'], filesep);
output_dir_AR_occl = strjoin([output_dir_AR, 'occl'], filesep);
% out dirs for peaks for multiple epoch lengths
output_dir_AR_peaks_medium = strjoin([output_dir_AR, 'peaks_-750_500'], filesep);
output_dir_AR_peaks_short = strjoin([output_dir_AR, 'peaks_-500_250'], filesep);
output_dir_AR_peaks_long = strjoin([output_dir_AR, 'peaks_-1000_750'], filesep);


%% Load Epoched data 

eeglab;

% short epochs
cd(output_dir_AR_peaks_short);
%list all *.set files in inputpath
file_names = dir('*.set');
%concatenate into one cell array
files2read = {file_names.name};
% load all eeg files
ALLEEG = pop_loadset('filename',files2read);
epochs_short = ALLEEG;

% medium epochs
cd(output_dir_AR_peaks_medium);
%list all *.set files in inputpath
file_names = dir('*.set');
%concatenate into one cell array
files2read = {file_names.name};
% load all eeg files
ALLEEG = pop_loadset('filename',files2read);
epochs_medium = ALLEEG;

% long epochs
cd(output_dir_AR_peaks_long);
%list all *.set files in inputpath
file_names = dir('*.set');
%concatenate into one cell array
files2read = {file_names.name};
% load all eeg files
ALLEEG = pop_loadset('filename',files2read);
epochs_long = ALLEEG;


%% Calculate Means Across Subjects 

% TODO: Put the means in a structure. 
% One mean each epoch across all subjects and conditions
% One mean each epoch across all subjects and separate for constant/random
% One mean each epoch across all subjects and separate for occluded/visible

for i = 1:size(ALLEEG, 2) 

    EEG = epochs_short(i);
    mean_struct_short.all(:,:,i) = mean(EEG.data,3);
    mean_struct_short.all_time_vec = EEG.times;

    EEG = epochs_medium(i);
    mean_struct_medium.all(:,:,i) = mean(EEG.data,3);
    mean_struct_medium.all_time_vec = EEG.times;

    EEG = epochs_long(i);
    mean_struct_long.all(:,:,i) = mean(EEG.data,3);
    mean_struct_long.all_time_vec = EEG.times;
%     EEG = OCC_ALLEEG(i);
%     % TODO: assign here the data from the different occlusion event structs
%     mean_struct.occ(:,:,i) = mean(EEG.data,3);
%     mean_struct.vis(:,:,i) = mean(EEG.data,3);
%     
%     mean_struct.occ_time_vec = EEG.times;

end


%% Save the Mean Matrices for Convenient Loading TODO

save()

%% Plot ERPs as grand average

neg_time = 60; 
avg_peak_dist = 300;

% Average ERP Plots

figure()
% short epochs
subplot(3, 1, 1)

short_mean = mean(mean_struct_short.all(:,:,:),3);
[epoch_lims(1), epoch_lims(2)] = bounds(epochs_short(8).times)
plot(mean_struct_short.all_time_vec, short_mean);
title(strcat(['ERP of all subjects averaged across epochs around peaks, n ~= ', num2str( size(epochs_short(1).data,3)), ' per subject']));
subtitle(strcat(['all channels, epochs from ', num2str(epoch_lims(1)), ' to ' num2str(epoch_lims(2))]))
hold on
h(1) = xline(0);
h(2) = xline(neg_time, 'r');
h(4:5) = xline([-avg_peak_dist,avg_peak_dist], 'b');
h(3) = xline(epoch_lims(1)+500, 'g')
legend(h, {'peak', 'error negativity peak',  'baseline period', 'average distance of next peak'});
hold off

%medium epochs
subplot(3, 1, 2)

medium_mean = mean(mean_struct_medium.all(:,:,:),3);
[epoch_lims(1), epoch_lims(2)] = bounds(epochs_medium(8).times)
plot(mean_struct_medium.all_time_vec, medium_mean);
title(strcat(['ERP of all subjects averaged across epochs around peaks, n ~= ', num2str( size(epochs_medium(1).data,3)), ' per subject']));
subtitle(strcat(['all channels, epochs from ', num2str(epoch_lims(1)), ' to ' num2str(epoch_lims(2))]))
hold on
h(1) = xline(0);
h(2) = xline(neg_time, 'r');
h(4:5) = xline([-avg_peak_dist,avg_peak_dist], 'b');
h(3) = xline(epoch_lims(1)+500, 'g')
legend(h, {'peak', 'error negativity peak',  'baseline period', 'average distance of next peak'});
hold off

% long epocjs
subplot(3, 1, 3)

long_mean = mean(mean_struct_long.all(:,:,:),3);
[epoch_lims(1), epoch_lims(2)] = bounds(epochs_long(8).times)
plot(mean_struct_long.all_time_vec, long_mean);
title(strcat(['ERP of all subjects averaged across epochs around peaks, n ~= ', num2str( size(epochs_long(1).data,3)), ' per subject']));
subtitle(strcat(['all channels, epochs from ', num2str(epoch_lims(1)), ' to ' num2str(epoch_lims(2))]))
hold on
h(1) = xline(0);
h(2) = xline(neg_time, 'r');
h(4:5) = xline([-avg_peak_dist,avg_peak_dist], 'b');
h(3) = xline(epoch_lims(1)+500, 'g')
legend(h, {'peak', 'error negativity peak',  'baseline period', 'average distance of next peak'});
hold off

% Single Channel ERP
% Identify the electrode that shows the biggest amplitude in the ERP

[~, neg_time_idx] = min(abs(mean_struct_short.all_time_vec - neg_time));
[~, min_chan_idx_short] = min(short_mean(:, neg_time_idx));

[~, neg_time_idx] = min(abs(mean_struct_medium.all_time_vec - neg_time));
[~, min_chan_idx_medium] = min(medium_mean(:, neg_time_idx));

[~, neg_time_idx] = min(abs(mean_struct_long.all_time_vec - neg_time));
[~, min_chan_idx_long] = min(long_mean(:, neg_time_idx));



figure()

subplot(3, 1, 1)

chan_lab = EEG.chanlocs(min_chan_idx_short).labels;
[epoch_lims(1), epoch_lims(2)] = bounds(epochs_short(8).times)
plot(mean_struct_short.all_time_vec, mean(mean_struct_short.all(min_chan_idx_short,:,:),3));
title(strcat(['ERP of all subjects averaged across epochs around peaks, channel ', chan_lab]));
subtitle(strcat(['epochs from ', num2str(epoch_lims(1)), ' to ' num2str(epoch_lims(2))]))
hold on
h(1) = xline(0);
h(2) = xline(neg_time, 'r');
h(4:5) = xline([-avg_peak_dist,avg_peak_dist], 'b');
h(3) = xline(epoch_lims(1)+500, 'g')
legend(h, {'peak', 'error negativity peak',  'baseline period', 'average distance of next peak'});
hold off

subplot(3, 1, 2)
[epoch_lims(1), epoch_lims(2)] = bounds(epochs_medium(8).times)
chan_lab = EEG.chanlocs(min_chan_idx_medium).labels;
plot(mean_struct_medium.all_time_vec, mean(mean_struct_medium.all(min_chan_idx_medium,:,:),3));
title(strcat(['ERP of all subjects averaged across epochs around peaks, channel ', chan_lab]));
subtitle(strcat(['epochs from ', num2str(epoch_lims(1)), ' to ' num2str(epoch_lims(2))]))
hold on
h(1) = xline(0);
h(2) = xline(neg_time, 'r');
h(4:5) = xline([-avg_peak_dist,avg_peak_dist], 'b');
h(3) = xline(epoch_lims(1)+500, 'g')
legend(h, {'peak', 'error negativity peak',  'baseline period', 'average distance of next peak'});
hold off

subplot(3, 1, 3)
chan_lab = EEG.chanlocs(min_chan_idx_long).labels;
[epoch_lims(1), epoch_lims(2)] = bounds(epochs_long(8).times)
plot(mean_struct_long.all_time_vec, mean(mean_struct_long.all(min_chan_idx_long,:,:),3));
title(strcat(['ERP of all subjects averaged across epochs around peaks, channel ', chan_lab]));
subtitle(strcat(['epochs from ', num2str(epoch_lims(1)), ' to ' num2str(epoch_lims(2))]))
hold on
h(1) = xline(0);
h(2) = xline(neg_time, 'r');
h(4:5) = xline([-avg_peak_dist,avg_peak_dist], 'b');
h(3) = xline(epoch_lims(1)+500, 'g')
legend(h, {'peak', 'error negativity peak',  'baseline period', 'average distance of next peak'});
hold off


%% Create Topoplot of the time points of interest TODO

MEAN = mean(mean_struct.all,3);
topoplot(mean(MEAN(:,neg_time)), EEG.chanlocs)
neg_time  = 60/4


%% Create ERP Plot of the data TODO

% TODO: sort by error size in epoch

pop_erpimage()


%% Stuff I tried out. 


% STUDY = std_erpplot(STUDY,ALLEEG)

% figure; 
% pop_plottopo(EEG, [1:60] , ''5STJS'', 0, ''ydir'',1);
% 
% pop_topoplot(EEG, 1, [-500 748] ,'5STJS',[1 2] ,0,'electrodes','on')
% 
% pop_erpplot()
% 
% 
% 
