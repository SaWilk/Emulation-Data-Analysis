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
subdir_peaks = strjoin([output_dir_epoched filesep 'peaks_-500_750'], filesep);


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


%% Load Epoched data 

eeglab;

% long epochs
cd(output_dir_AR_peaks);
%list all *.set files in inputpath
file_names = dir('*.set');
%concatenate into one cell array
files2read = {file_names.name};
% load all eeg files
ALLEEG = pop_loadset('filename',files2read);
epochs = ALLEEG;


%% Calculate Means Across Subjects 

% TODO: Put the means in a structure. 
% One mean each epoch across all subjects and conditions
% One mean each epoch across all subjects and separate for constant/random
% One mean each epoch across all subjects and separate for occluded/visible

for i = 1:size(ALLEEG, 2) 

    EEG = epochs(i);
    mean_struct.all(:,:,i) = mean(EEG.data,3);
    mean_struct.all_time_vec = EEG.times;
%     EEG = OCC_ALLEEG(i);
%     % TODO: assign here the data from the different occlusion event structs
%     mean_struct.occ(:,:,i) = mean(EEG.data,3);
%     mean_struct.vis(:,:,i) = mean(EEG.data,3);
%     
%     mean_struct.occ_time_vec = EEG.times;

end


%% Save the Mean Matrices for Convenient Loading

cd(mean_matrices_peaks_epoched_path)
save("mean_struct.mat", 'mean_struct')


%% Load Previously Generated Mean Structures

load('mean_struct.mat')


%% Plot ERPs as grand average

neg_time = 200; 
avg_peak_dist = 300;

epoch_mean = mean(mean_struct.all(:,:,:),3);

% Average ERP Plots

figure()
% long epocjs

[epoch_lims(1), epoch_lims(2)] = bounds(epochs(8).times)
plot(mean_struct.all_time_vec, epoch_mean);
title(strcat(['ERP of all subjects averaged across epochs around peaks, n ~= ', num2str( size(epochs(1).data,3)), ' per subject']));
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
neg_time = 60; 

[~, neg_time_idx] = min(abs(mean_struct.all_time_vec - neg_time));
[~, min_chan_idx] = min(epoch_mean(:, neg_time_idx));

figure()

chan_lab = EEG.chanlocs(min_chan_idx).labels;
[epoch_lims(1), epoch_lims(2)] = bounds(epochs(8).times)
plot(mean_struct.all_time_vec, mean(mean_struct.all(min_chan_idx,:,:),3));
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

latencies_onset = linspace(0, 550, 12);
latencies_offset = latencies_onset + 50;
% TODO: Display all topoplots in the same color range
% doc topoplot: https://rdrr.io/cran/erpR/man/topoplot.html

% extract means across 50 ms periods for topoplot
for i = 1:length(latencies_onset)

    [~, lat_idx_onset] = min(abs(mean_struct.all_time_vec - latencies_onset(i)));
    [~, lat_idx_offset] = min(abs(mean_struct.all_time_vec - latencies_offset(i)));
    topo_mean(:, i) = mean(epoch_mean(:,lat_idx_onset:lat_idx_offset), 2);

end

zlims = [min(topo_mean), max(topo_mean)];

figure()
for i = 1:length(latencies)

    subplot(3, 4, i)
    topoplot(topo_mean(:, i), EEG.chanlocs, 'zlim', zlims)
    title('epochs around peaks all subjects all trials, long epochs')
    subtitle(strcat(['latency range: ', num2str(latencies_onset(i)), ' to ', ...
        num2str(latencies_onset(i)), ' ms relative to peak']))
    colorbar()

end


%% Perform cluster-Based Permutation Test TODO

% https://www.fieldtriptoolbox.org/tutorial/cluster_permutation_timelock/

% use this to find the electrode cluster 


%% Create ERP Plot of the data TODO

% TODO: sort by error size in epoch
% sort by

pop_erpimage()


%% Do TF Decomposition on data

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
% https://www.fieldtriptoolbox.org/tutorial/cluster_permutation_timelock/


