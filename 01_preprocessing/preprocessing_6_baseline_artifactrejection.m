% Data Analysis Pursuit-Tracking and Pursuit-Occlusion Paradigm
% Emulation pilot study 2021

% This script contains: 
% baseline correction
% remove bad epochs
% exclude epochs before practice trials

% export data at different steps

%% clear workspace
clear;
clc;
format compact
format long G


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





%% add this to loop


    %% delete practice trials

    exp_start_marker = 11;

    % first, task A
    % extract the marker_index from urevent S11
    marker_index = 1;
    found = false;
    while found == false
        if str2double(TMPEEG_A.urevent(marker_index).type(2:end)) == exp_start_marker
            found = true;
        else
            marker_index = marker_index + 1;
        end
    end
    % extract event index for this marker in event structure
    event_index = 1;
    found = false;
    while found == false
        if TMPEEG_A.event(event_index).bvmknum > marker_index
            found = true;
        else
            event_index = event_index + 1;
        end
    end
    % reduce events to those after practice, save events with practice to
    % extra cell
    TMPEEG_A.event_with_practice = TMPEEG_A.event;
    TMPEEG_A.event = TMPEEG_A.event_with_practice(event_index:end);
    %save without practice
    TMPEEG_A.setname = [filename '_epoched_nopractice_A'];
    TMPEEG_A = pop_saveset(TMPEEG_A, 'filename', TMPEEG_A.setname, 'filepath', char(savepath_no_practice));

    % second, task B
    % extract the marker_index from urevent S11
    marker_index = 1;
    found = false;
    while found == false
        if str2double(TMPEEG_B.urevent(marker_index).type(2:end)) == exp_start_marker
            found = true;
        else
            marker_index = marker_index + 1;
        end
    end
    % extract event index for this marker in event structure
    event_index = 1;
    found = false;
    while found == false
        if TMPEEG_B.event(event_index).bvmknum > marker_index
            found = true;
        else
            event_index = event_index + 1;
        end
    end
    % reduce events to those after practice, save events with practice to
    % extra cell
    TMPEEG_B.event_with_practice = TMPEEG_B.event;
    TMPEEG_B.event = TMPEEG_B.event_with_practice(event_index:end);
    %save without practicew
    TMPEEG_B.setname = [filename '_epoched_nopractice_B'];
    TMPEEG_B = pop_saveset(TMPEEG_B, 'filename', TMPEEG_B.setname, 'filepath', char(savepath_no_practice));


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


