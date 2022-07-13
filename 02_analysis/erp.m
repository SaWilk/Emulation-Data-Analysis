%% Make ERP Plot

% Author: Saskia Wilken
% creation Date: 28.06.2022


%% Empty

format compact
format long G
clear
clc


%% folders and dependencies

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
subdir_peaks = strjoin([output_dir_epoched filesep "peaks"], filesep);

savepath_baseline_const_rand = strjoin([output_dir_baseline, "const_rand"], filesep);
savepath_baseline_occl = strjoin([output_dir_baseline, 'occl_nonoccl'], filesep);
savepath_baseline_peaks = strjoin([output_dir_baseline, 'peaks'], filesep);

output_dir_AR_const = strjoin([output_dir_AR, 'const'], filesep);
output_dir_AR_occl = strjoin([output_dir_AR, 'occl'], filesep);
output_dir_AR_peaks = strjoin([output_dir_AR, 'peaks'], filesep);




%% Load all the data 

cd(output_dir_AR_peaks);

%list all *.set files in inputpath
file_names = dir('*.set');

%concatenate into one cell array
files2read = {file_names.name};

eeglab;
% load all eeg files
ALLEEG = pop_loadset('filename',files2read);
% load(dir('*.mat').name);


%% Calculate Means Across Subjects 

% TODO: Put the means in a structure. 
% One mean each epoch across all subjects and conditions
% One mean each epoch across all subjects and separate for constant/random
% One mean each epoch across all subjects and separate for occluded/visible

for i = 1:size(ALLEEG, 2)-1 % because last subject has no latencies yet. 

    EEG = ALLEEG(i);
    mean_struct.all(:,:,i) = mean(EEG.data,3);
    mean_struct.all_time_vec = EEG.times;
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

chan_no = 2; % some channel 
chan_lab = EEG.chanlocs(chan_no).labels;
figure()
plot(mean_struct.all_time_vec, mean(mean_struct.all(chan_no,:,:),3));
title(strcat(['ERP of all subjects averaged across epochs around peaks, channel ', chan_lab]));


figure()
plot(mean_struct.all_time_vec, mean(mean_struct.all(:,:,:),3));
title(strcat(['ERP of all subjects averaged across epochs around peaks, all channels']));
hold on
h(1) = xline(0);
h(2) = xline(60, 'r');
h(3:4) = xline([-300,300], 'b');
legend(h, {'peak', 'error negativity peak', 'average distance of next peak'});
hold off


% TODO: Experiment with shorter epoch durations

% TODO: Identify the electrode that shows the biggest amplitude in the ERP
mean(mean_struct.all(chan_no,:,:),3)
% max amplitude 


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
