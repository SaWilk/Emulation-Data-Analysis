% Data Analysis Pursuit-Tracking and Pursuit-Occlusion Paradigm
% Emulation pilot study 2021

% This script contains: 
% * epoching data for whole trial (S26 + 9s)
% * removing baseline of 750 ms before fixation cross
% * epoching again for smaller intervals
% * artifact rejection 

% Adriana Böttcher
% 13.07.22

%% clear workspace
clear;
clc;


%% folders and dependencies

% add path and start EEGlab toolbox
addpath('R:\AG-Beste-Orga\Skripts\Toolbox\eeglab2021.0');
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
output_dir_epoched = strjoin([parent_dir_2, "Emulation-Data-Output\04_epoched_b"], filesep);
output_dir_AR = strjoin([parent_dir_2, "Emulation-Data-Output\06_artifact_rejection_b"], filesep);

cd(input_dir);

%list all *.set files in inputpath
file_names = dir('*_EEG.set');

%concatenate into one cell array
files2read = {file_names.name};

%% set parameters

epoch_lims_trial = [-0.75 13];
epoch_lims_traj = [0 3];
epoch_lims_occl = [0 2];

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
    %TMPEEG_A = add_traj_events(TMPEEG_A);
    TMPEEG_A.setname = [TMPEEG_A.setname  '_without_B_events'];

    TMPEEG_B = TMPEEG;
    TMPEEG_B.event = TMPEEG_B.event(start_b:end);
    %TMPEEG_B = add_occl_events(TMPEEG_B);
    TMPEEG_B.setname = [TMPEEG_B.setname  '_without_A_events'];

    %% epoch data for trials in task A

    %epoch continuous data for whole trial
    TMPEEG_A = pop_epoch( TMPEEG_A, {event_fixcross}, epoch_lims_trial, 'newname', ...
        [TMPEEG.setname '_epoched_trials'], 'epochinfo', 'no');
    TMPEEG_A.comment = TMPEEG_A.comment + "  *** epoch for trials";
    TMPEEG_A.setname = [set_name '_epoched_trials_A.set'];
    %TMPEEG_A = pop_saveset(TMPEEG_A, 'filename', TMPEEG_A.setname, 'filepath', ...
     %  char(output_dir_epoched));

     %% baseline correction for fixed interval around fixation cross task A

    TMPEEG_A = pop_rmbase(TMPEEG_A, base_lims_trial ,[]);
    TMPEEG_A.comment = TMPEEG_A.comment + "  *** apply baseline corr around fixcross";
    TMPEEG_A.setname = [set_name '_baseline_corr_whole_trial_A.set'];
    %TMPEEG_A = pop_saveset(TMPEEG_A, 'filename', TMPEEG_A.setname, 'filepath', ...
     %   char(savepath_baseline_const_rand));

     %% epoch data for trials in task B

    %epoch continuous data for whole trial
    if ~strcmp(TMPEEG_A.subject, 'KMY6K') % skip subject 18 Task B cause it is only bad trials
        TMPEEG_B = pop_epoch( TMPEEG_B, {event_fixcross}, epoch_lims_trial, 'newname', ...
            [TMPEEG.setname '_epoched_trials'], 'epochinfo', 'no');
        TMPEEG_B.comment = TMPEEG_B.comment + "  *** epoch for trials";
        TMPEEG_B.setname = [set_name '_epoched_trials_B.set'];
        %TMPEEG_B = pop_saveset(TMPEEG_B, 'filename', TMPEEG_B.setname, 'filepath', ...
         %  char(output_dir_epoched));
    end

     %% baseline correction for fixed interval around fixcross for task B
if ~strcmp(TMPEEG_A.subject, 'KMY6K') % skip subject 18 Task B cause it is only bad trials
    TMPEEG_B = pop_rmbase(TMPEEG_B, base_lims_trial ,[]);
    TMPEEG_B.comment = TMPEEG_B.comment + "  *** apply baseline corr around fixcross";
    TMPEEG_B.setname = [set_name '_baseline_corr_whole_trial_B.set'];
    %TMPEEG_B = pop_saveset(TMPEEG_B, 'filename', TMPEEG_B.setname, 'filepath', ...
     %   char(savepath_baseline_const_rand));
end
     %% epoch data again, for traj parts now

    TMPEEG_A = pop_epoch( TMPEEG_A, event_traj , epoch_lims_traj, 'newname', ...
        [TMPEEG_A.setname '_A_epoched_traj'], 'epochinfo', 'yes');
    TMPEEG_A.comment = TMPEEG_A.comment + "  *** epoch for traj parts";
    TMPEEG_A.setname = [set_name '_epoched_traj_parts.set'];
    %TMPEEG_A = pop_saveset(TMPEEG_A, 'filename', TMPEEG_A.setname, 'filepath', ...
    %    char(subdir_const_rand));

    %% epoch data again, for occlusion parts now
if ~strcmp(TMPEEG_A.subject, 'KMY6K') % skip subject 18 Task B cause it is only bad trials
    TMPEEG_B = pop_epoch(TMPEEG_B, event_occl , epoch_lims_occl, 'newname',...
        [TMPEEG_B.setname '_B_epoched'], 'epochinfo', 'yes');
    TMPEEG_B.comment = TMPEEG_B.comment + "  *** epoch for occl parts";
    TMPEEG_B.setname = [set_name '_epoched_occl_parts.set'];
%     TMPEEG_B = pop_saveset(TMPEEG_B, 'filename', TMPEEG_B.setname, 'filepath', ...
%         char(subdir_occl));
end
     %% artifact rejection task A with new baseline 

    TMPEEG_A = pop_jointprob(TMPEEG_A, 1, 1:TMPEEG_A.nbchan, 5, 5, 0, 1);
    TMPEEG_A.comment = TMPEEG_A.comment + "  *** apply artifact rejection";
    TMPEEG_A.setname = [set_name '_complete_preprocessing_withAR_A'];
    TMPEEG_A = pop_saveset(TMPEEG_A, 'filename', TMPEEG_A.setname, 'filepath', ...
        char(output_dir_AR));

    %% artifact rejection task B with new baseline 
if ~strcmp(TMPEEG_A.subject, 'KMY6K') % skip subject 18 Task B cause it is only bad trials
    TMPEEG_B = pop_jointprob(TMPEEG_B, 1, 1:TMPEEG_B.nbchan, 5, 5, 0, 1);
    TMPEEG_B.comment = TMPEEG_B.comment + "  *** apply artifact rejection";
    TMPEEG_B.setname = [set_name '_complete_preprocessing_withAR_B'];
    TMPEEG_B = pop_saveset(TMPEEG_B, 'filename', TMPEEG_B.setname, 'filepath', ...
        char(output_dir_AR));
end
end


%% check mean fixation cross duration
% diff = [];
% index = 1;
% for i = 1:size(TMPEEG.event, 2)
%     if TMPEEG.event(i).type == "S 26"
%         diff(index) = TMPEEG.event(i).trial_latency - TMPEEG.event(i-1).trial_latency;
%         index = index + 1;
%     end
% end
% mean(diff)