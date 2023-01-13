% The neurophysiology of continuous action monitoring

% preprocessing task A&B, pt. 5

% created by:
% Adriana Boettcher, Cognitive Neurophysiology TU Dresden
% 2022

%% clear workspace
clear;
clc;


%% folders and dependencies

% add path and start EEGlab toolbox
eeglab;
close;

%list all *.set files in inputpath
file_names = dir('*_EEG.set');

%concatenate into one cell array
files2read = {file_names.name};

%% set parameters

epoch_lims_trial = [-1 13];
epoch_lims_traj = [-1 4]; % big time interval for wavelet analysis, reduce later to [0 3]
epoch_lims_occl = [-1 4]; % big time interval for wavelet analysis, reduce later to [0 2]

event_fixcross = 'S 26';
event_traj = {'S 27' 'S 23' 'S 24' };
event_occl = {'S 20' 'S 21'};

base_lims_trial = [-750 0];

%% read files, 

for ind = 1:length(file_names)
    
    file_name = files2read{ind};
    TMPEEG = pop_loadset(file_name);
    set_name = TMPEEG.subject;

    %% select only events of task A/B

    task_b = false;
    start_b = 1;

    while task_b == false
        if strcmp(TMPEEG.event(start_b).task, 'task_b')
            task_b = true;
        else
            start_b = start_b + 1;
        end
    end

    TMPEEG_A = TMPEEG;
    TMPEEG_A.event = TMPEEG_A.event(1:start_b);
    TMPEEG_A.setname = [TMPEEG_A.setname  '_without_B_events'];

    TMPEEG_B = TMPEEG;
    TMPEEG_B.event = TMPEEG_B.event(start_b:end);
    TMPEEG_B.setname = [TMPEEG_B.setname  '_without_A_events'];

    %% epoch data for trials in task A

    %epoch continuous data for whole trial
    TMPEEG_A = pop_epoch( TMPEEG_A, {event_fixcross}, epoch_lims_trial, 'newname', ...
        [TMPEEG.setname '_epoched_trials'], 'epochinfo', 'no');
    TMPEEG_A.comment = TMPEEG_A.comment + "  *** epoch for trials";
    TMPEEG_A.setname = [set_name '_epoched_trials_A.set'];

     %% baseline correction for fixed interval around fixation cross task A

    TMPEEG_A = pop_rmbase(TMPEEG_A, base_lims_trial ,[]);
    TMPEEG_A.comment = TMPEEG_A.comment + "  *** apply baseline corr around fixcross";
    TMPEEG_A.setname = [set_name '_baseline_corr_whole_trial_A.set'];

     %% epoch data for trials in task B

    %epoch continuous data for whole trial
    if ~strcmp(TMPEEG_A.subject, 'KMY6K') % skip subject 18 Task B because too many trials are labelled as erroneous
        TMPEEG_B = pop_epoch( TMPEEG_B, {event_fixcross}, epoch_lims_trial, 'newname', ...
            [TMPEEG.setname '_epoched_trials'], 'epochinfo', 'no');
        TMPEEG_B.comment = TMPEEG_B.comment + "  *** epoch for trials";
        TMPEEG_B.setname = [set_name '_epoched_trials_B.set'];
    end

     %% baseline correction for fixed interval around fixcross for task B
    if ~strcmp(TMPEEG_A.subject, 'KMY6K') 
        TMPEEG_B = pop_rmbase(TMPEEG_B, base_lims_trial ,[]);
        TMPEEG_B.comment = TMPEEG_B.comment + "  *** apply baseline corr around fixcross";
        TMPEEG_B.setname = [set_name '_baseline_corr_whole_trial_B.set'];
    end
     %% epoch data again, for traj parts now

    TMPEEG_A = pop_epoch( TMPEEG_A, event_traj , epoch_lims_traj, 'newname', ...
        [TMPEEG_A.setname '_A_epoched_traj'], 'epochinfo', 'yes');
    TMPEEG_A.comment = TMPEEG_A.comment + "  *** epoch for traj parts";
    TMPEEG_A.setname = [set_name '_epoched_traj_parts.set'];

    %% epoch data again, for occlusion parts now
    if ~strcmp(TMPEEG_A.subject, 'KMY6K') 
        TMPEEG_B = pop_epoch(TMPEEG_B, event_occl , epoch_lims_occl, 'newname',...
            [TMPEEG_B.setname '_B_epoched'], 'epochinfo', 'yes');
        TMPEEG_B.comment = TMPEEG_B.comment + "  *** epoch for occl parts";
        TMPEEG_B.setname = [set_name '_epoched_occl_parts.set'];
    end
     %% artifact rejection task A with new baseline 

    TMPEEG_A = pop_jointprob(TMPEEG_A, 1, 1:TMPEEG_A.nbchan, 5, 5, 0, 1);
    TMPEEG_A.comment = TMPEEG_A.comment + "  *** apply artifact rejection";
    TMPEEG_A.setname = [set_name '_complete_preprocessing_withAR_A'];
    TMPEEG_A = pop_saveset(TMPEEG_A, 'filename', TMPEEG_A.setname, 'filepath', ...
        char(output_dir_AR));

    %% artifact rejection task B with new baseline 
    if ~strcmp(TMPEEG_A.subject, 'KMY6K') 
        TMPEEG_B = pop_jointprob(TMPEEG_B, 1, 1:TMPEEG_B.nbchan, 5, 5, 0, 1);
        TMPEEG_B.comment = TMPEEG_B.comment + "  *** apply artifact rejection";
        TMPEEG_B.setname = [set_name '_complete_preprocessing_withAR_B'];
        TMPEEG_B = pop_saveset(TMPEEG_B, 'filename', TMPEEG_B.setname, 'filepath', ...
            char(output_dir_AR));
    end
end