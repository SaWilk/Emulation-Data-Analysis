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
epoch_peaks_path = strcat(file_path, "\02_epoched_peaks_data"); 
if ~mkdir(epoch_peaks_path)
    mkdir(epoch_peaks_path);
end
cd(data_path);


%% Load all the data 


% get all set files in the folder
files2read = {dir('*.set').name};

eeglab;
% load all eeg files
ALLEEG = pop_loadset('filename',files2read);
load(dir('*.mat').name);
track_data = copy_data;


%% Set Parameters

EPOCH_LIMS = [-.8, 0.6]; % if epochs are shorter I get an out of memory error
BASE_LIMS = [-800 -300];
EVENT = 'S 40';


%% Epoch data around peaks and save

% this creates as many epochs as possible that do not overlap. Each epoch,
% however, will have multiple peaks

epoched_data_suffix = '_epoched_peaks';
% COPY_ALLEEG = ALLEEG;

for i = 1:size(ALLEEG, 2)-1 % because last subject has no latencies yet. 

    EEG = ALLEEG(i);
    EEG = pop_epoch( EEG, {  EVENT  }, EPOCH_LIMS, 'newname', EEG.subject, 'epochinfo', 'yes');
    EEG = eeg_checkset( EEG );
    EEG = pop_rmbase(EEG, BASE_LIMS);
    event_idx = find(strcmp({EEG.event.type}, EVENT));
    % plot that shows that peak distances are not normally distributed
%     plot(sort(diff([EEG.event(event_idx).latency])*1000/EEG.srate))
    peak_distances(i) = median(diff([EEG.event(event_idx).latency])*1000/EEG.srate);
    % save epoched data
    file_name = strcat([EEG.subject, epoched_data_suffix]);
    EEG = pop_saveset(EEG, 'filename', file_name, 'filepath', char(epoch_peaks_path));

    [ALLEEG, EEG, i] = eeg_store(ALLEEG, EEG, i);

end

mean(peak_distances)
% peaks are usually around 300 ms apart...
PEAK_ALLEEG = ALLEEG;


%% Set Parameters for occlusion

EPOCH_LIMS = [-750, 0.5];
BASE_LIMS = [-750 -250];
EVENT = 'S 20';
EVENT2 = 'S 21';


%% Epoch data around Occlusion events


epoched_data_suffix = '_epoched_occlusion';
ALLEEG = COPY_ALLEEG;

for i = 1:size(ALLEEG, 2)-1 % because last subject has no latencies yet. 

    EEG = ALLEEG(i);
    % TODO: Create two different sets of epochs probably using different
    % epoching calls. 
    EEG = pop_epoch( EEG, {  EVENT  }, EPOCH_LIMS, 'newname', EEG.subject, 'epochinfo', 'yes');
    EEG = eeg_checkset( EEG );
    EEG = pop_rmbase(EEG, BASE_LIMS);
    file_name = strcat([EEG.subject, epoched_data_suffix]);
    EEG = pop_saveset(EEG, 'filename', file_name, 'filepath', char(epoch_peaks_path));

    [ALLEEG, EEG, i] = eeg_store(ALLEEG, EEG, i);

end


OCC_ALLEEG = ALLEEG;


%% Calculate Means Across Subjects 

% TODO: Put the means in a structure. 
% One mean each epoch across all subjects and conditions
% One mean each epoch across all subjects and separate for constant/random
% One mean each epoch across all subjects and separate for occluded/visible

for i = 1:size(ALLEEG, 2)-1 % because last subject has no latencies yet. 

    EEG = PEAK_ALLEEG(i);
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
