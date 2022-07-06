%% Epoching Data

% Author: Saskia Wilken (and soon Adriana Böttcher)
% Creation Date: 06.07.2022

% all the epoching options in one script including preprocessing steps that
% require epoched data like baseline correction and rejection of noisy epochs

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

% TODO: Clean Up this Script and bring everything in order. Each epoching
% step should save datasets in a separate folder which can be loaded
% conveniently for analyses (done in a separate script)


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


   %% separate datasets & epochs
    
    %epoch continuous data for constant and random traj and save as
    %TMPEEG_A
    TMPEEG_A = pop_epoch( TMPEEG, {'S 23' 'S 24' 'S 27'}, [-1.5 3], 'newname', [TMPEEG.setname '_A_epoched'], 'epochinfo', 'yes');
    TMPEEG_A.comment = TMPEEG_A.comment + "  *** epoch for traj parts";
    TMPEEG_A.setname = [filename '_icaclean_epoched_A'];
    TMPEEG_A = pop_saveset(TMPEEG_A, 'filename', TMPEEG_A.setname, 'filepath', char(savepath_epoched));

    %epoch continuous data for occlusion/ non-occlusion and save as
    %TMPEEG_B
    TMPEEG_B = pop_epoch( TMPEEG, {'S 20' 'S 21'}, [-1.5 2], 'newname', [TMPEEG.setname '_B_epoched'], 'epochinfo', 'yes');
    TMPEEG_B.comment = TMPEEG_B.comment + "  *** epoch for occl parts";
    TMPEEG_B.setname = [filename '_icaclean_epoched_B'];
    TMPEEG_B = pop_saveset(TMPEEG_B, 'filename', TMPEEG_B.setname, 'filepath', char(savepath_epoched));

    %% apply baseline correction and save

    TMPEEG_A = pop_rmbase(TMPEEG_A, [-1000 0] ,[]);
    TMPEEG_A.comment = TMPEEG_A.comment + "  *** apply baseline corr";
    TMPEEG_A.setname = [filename '_complete_preprocessing_A'];
    TMPEEG_A = pop_saveset(TMPEEG_A, 'filename', TMPEEG_A.setname, 'filepath', char(savepath_baseline));

    TMPEEG_B = pop_rmbase(TMPEEG_B, [-1000 0] ,[]);
    TMPEEG_B.comment = TMPEEG_B.comment + "  *** apply baseline corr";
    TMPEEG_B.setname = [filename '_complete_preprocessing_B'];
    TMPEEG_B = pop_saveset(TMPEEG_B, 'filename', TMPEEG_B.setname, 'filepath', char(savepath_baseline));

    %% artifact rejection and save
    
    %parameters: input, ICs (0) or EEG data (1), electrodes, thresholds in
    %stdv for electrodes, global threshold for all electrodes, superpose
    %prelabelling with previous labelling, keep (0) or reject (1) trials

    TMPEEG_A = pop_jointprob(TMPEEG_A, 1, 1:TMPEEG_A.nbchan, 5, 5, 0, 1);
    TMPEEG_A.comment = TMPEEG_A.comment + "  *** apply artifact rejection";
    TMPEEG_A.setname = [filename '_complete_preprocessing_withAR_A'];
    TMPEEG_A = pop_saveset(TMPEEG_A, 'filename', TMPEEG_A.setname, 'filepath', char(savepath_artifactrej));

    TMPEEG_B = pop_jointprob(TMPEEG_B, 1, 1:TMPEEG_B.nbchan, 5, 5, 0, 1);
    TMPEEG_B.comment = TMPEEG_B.comment + "  *** apply artifact rejection";
    TMPEEG_B.setname = [filename '_complete_preprocessing_withAR_B'];
    TMPEEG_B = pop_saveset(TMPEEG_B, 'filename', TMPEEG_B.setname, 'filepath', char(savepath_artifactrej));
