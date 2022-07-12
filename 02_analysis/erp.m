%% Make ERP Plot

% Author: Saskia Wilken
% creation Date: 28.06.2022


%% Empty

format compact
format long G
clear
clc


%% Set Paths

% file location
file_path = matlab.desktop.editor.getActiveFilename;
file_path = fileparts(file_path);
addpath(file_path);

% data_path = strcat(file_path, "\01_merged_data"); % get root path to data folder
epoch_peaks_path = strcat(file_path, "\02_epoched_peaks_data"); 
if ~mkdir(epoch_peaks_path)
    mkdir(epoch_peaks_path);
end
cd(epoch_peaks_path);


%% Load all the data 

% get all set files in the folder
files2read = {dir('*.set').name};

eeglab;
% load all eeg files
ALLEEG = pop_loadset('filename',files2read);
load(dir('*.mat').name);


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