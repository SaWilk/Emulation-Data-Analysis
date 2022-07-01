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

data_path = strcat(file_path, "\01_merged_data"); % get root path to data folder

cd(data_path);

%% Load all the data 


% get all set files in the folder
files2read = {dir('*.set').name};

eeglab;
% load all eeg files
ALLEEG = pop_loadset('filename',files2read);
load(dir('*.mat').name);
track_data = copy_data;


%% Calculate ERPs


for i = 1:size(ALLEEG, 2)-1 % because last subject has no latencies yet. 

    EEG = ALLEEG(i);
%     max(diff([EEG.event(strcmp({EEG.event.type}, 'S 40')).trial_latency_ms]) )
    EEG = pop_epoch( EEG, {  'S 40'  }, [-1  0.75], 'newname', EEG.subject, 'epochinfo', 'yes');
    EEG = eeg_checkset( EEG );
    EEG = pop_rmbase(EEG, [-1000 -500]);
    [ALLEEG, EEG, i] = eeg_store(ALLEEG, EEG, i);
    mean_erps(:,:,i)  = mean(EEG.data,3);

end
mean_across_subj_erps= mean(mean_erps,3);
time_vec = EEG.times;

figure; 
pop_plottopo(EEG, [1:60] , ''5STJS'', 0, ''ydir'',1);

pop_topoplot(EEG, 1, [-500 748] ,'5STJS',[1 2] ,0,'electrodes','on')

pop_erpplot()

STUDY = std_erpplot(STUDY,ALLEEG)

topoplot(mean_across_subj_erps(2,:), EEG.chanlocs)

% TMPEEG = pop_saveset(TMPEEG,'filename',[files2read{ind}(1:end-7) '_A_epoched'], 'filepath', char(savepath));

chan_no = 2; % some channel 
chan_lab = EEG.chanlocs(chan_no).labels;
figure()
plot(time_vec, mean_across_subj_erps(chan_no,:));
title(strcat(['ERP of all subjects averaged across epochs around peaks, channel ', chan_lab]));


figure()
plot(time_vec, mean_across_subj_erps);
title(strcat(['ERP of all subjects averaged across epochs around peaks, all channels']));







