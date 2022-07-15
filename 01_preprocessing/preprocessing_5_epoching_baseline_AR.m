% Data Analysis Pursuit-Tracking and Pursuit-Occlusion Paradigm
% Emulation pilot study 2021

% This script contains: 
% * separating the data for both tasks by:
% epoching based on respective eventsdeleted plots and problems
% * baseline correction
% * remove bad epochs
% * export data at different steps

% Adriana Böttcher, edits by Saskia Wilken
% 08.07.22


%% clear workspace
clear;
clc;


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
% out dirs for peaks for multiple epoch lengths
subdir_peaks_medium = strjoin([output_dir_epoched filesep 'peaks_-750_500'], filesep);
mkdir(subdir_peaks_medium)
subdir_peaks_long = strjoin([output_dir_epoched filesep 'peaks_-1000_750'], filesep);
mkdir(subdir_peaks_long)
subdir_peaks_short = strjoin([output_dir_epoched filesep 'peaks_-500_250'], filesep);
mkdir(subdir_peaks_short)


savepath_baseline_const_rand = strjoin([output_dir_baseline, "const_rand"], filesep);
savepath_baseline_occl = strjoin([output_dir_baseline, 'occl_nonoccl'], filesep);
% out dirs for peaks for multiple epoch lengths
savepath_baseline_peaks_medium = strjoin([output_dir_baseline, 'peaks_-750_500'], filesep);
mkdir(savepath_baseline_peaks_medium)
savepath_baseline_peaks_short = strjoin([output_dir_baseline, 'peaks_-500_250'], filesep);
mkdir(savepath_baseline_peaks_short)
savepath_baseline_peaks_long = strjoin([output_dir_baseline, 'peaks_-1000_750'], filesep);
mkdir(savepath_baseline_peaks_long)


output_dir_AR_const = strjoin([output_dir_AR, 'const'], filesep);
output_dir_AR_occl = strjoin([output_dir_AR, 'occl'], filesep);
% out dirs for peaks for multiple epoch lengths
output_dir_AR_peaks_medium = strjoin([output_dir_AR, 'peaks_-750_500'], filesep);
mkdir(output_dir_AR_peaks_medium)
output_dir_AR_peaks_short = strjoin([output_dir_AR, 'peaks_-500_250'], filesep);
mkdir(output_dir_AR_peaks_short)
output_dir_AR_peaks_long = strjoin([output_dir_AR, 'peaks_-1000_750'], filesep);
mkdir(output_dir_AR_peaks_long)


cd(input_dir);

%list all *.set files in inputpath
file_names = dir('*_EEG.set');

%concatenate into one cell array
files2read = {file_names.name};


%% set parameters for epoching

epoch_lims_A = [-1.5 3];
epoch_lims_B =  [-1.5 2];
epoch_lims_peaks_long = [-1, 0.75]; 
epoch_lims_peaks_medium = [-.75, 0.5]; 
epoch_lims_peaks_short = [-.5, 0.25]; 


event_A = {'S 23' 'S 24' 'S 27'};
event_B = {'S 20' 'S 21'};
event_peaks = 'S 40';

base_lims_A = [-1000 0];
base_lims_B = [-1000 0];
base_lims_peaks_medium = [-750 -250];
base_lims_peaks_long = [-1000 -500];
base_lims_peaks_short = [-500 0];


%% loop through files and save epoched data

for ind = 1:length(file_names)
    
    file_name = files2read{ind};
    TMPEEG = pop_loadset(file_name);
    set_name = TMPEEG.subject;
    TMPEEG.old_urevent = TMPEEG.urevent;
    TMPEEG.urevent = TMPEEG.event;
    
    %% epoching & remove baseline for A, B, peaks

    %epoch continuous data for constant and random traj and save as
    %TMPEEG_A
    TMPEEG_A = pop_epoch( TMPEEG, event_A , epoch_lims_A, 'newname', ...
        [TMPEEG.setname '_A_epoched'], 'epochinfo', 'yes');
    TMPEEG_A.comment = TMPEEG_A.comment + "  *** epoch for traj parts";
    TMPEEG_A.setname = [set_name '_epoched_A.set'];
    TMPEEG_A = pop_saveset(TMPEEG_A, 'filename', TMPEEG_A.setname, 'filepath', ...
        char(subdir_const_rand));
    % remove baseline
    TMPEEG_A = pop_rmbase(TMPEEG_A, base_lims_A ,[]);
    TMPEEG_A.comment = TMPEEG_A.comment + "  *** apply baseline corr";
    TMPEEG_A.setname = [set_name '_complete_preprocessing_A.set'];
    TMPEEG_A = pop_saveset(TMPEEG_A, 'filename', TMPEEG_A.setname, 'filepath', ...
        char(savepath_baseline_const_rand));

if ~strcmp(TMPEEG_A.subject, 'KMY6K') % skip subject 18 Task B cause it is only bad trials
    %epoch continuous data for occlusion/ non-occlusion and save as
    %TMPEEG_B
    TMPEEG_B = pop_epoch(TMPEEG, event_B , epoch_lims_B, 'newname',...
        [TMPEEG.setname '_B_epoched'], 'epochinfo', 'yes');
    TMPEEG_B.comment = TMPEEG_B.comment + "  *** epoch for occl parts";
    TMPEEG_B.setname = [set_name '_epoched_B.set'];
    TMPEEG_B = pop_saveset(TMPEEG_B, 'filename', TMPEEG_B.setname, 'filepath', ...
        char(subdir_occl));
    % remove baseline
    TMPEEG_B = pop_rmbase(TMPEEG_B, base_lims_B ,[]);
    TMPEEG_B.comment = TMPEEG_B.comment + "  *** apply baseline corr";
    TMPEEG_B.setname = [set_name '_complete_preprocessing_B.set'];
    TMPEEG_B = pop_saveset(TMPEEG_B, 'filename', TMPEEG_B.setname, 'filepath', ...
        char(savepath_baseline_occl));
end

    % epoch for peaks & remove baseline, short duration
    TMPEEG_peaks_short = pop_epoch(TMPEEG, { event_peaks }, epoch_lims_peaks_short,...
        'newname', TMPEEG.subject, 'epochinfo', 'yes');
    TMPEEG_peaks_short = eeg_checkset( TMPEEG_peaks_short );
    TMPEEG_peaks_short = pop_rmbase(TMPEEG_peaks_short, base_lims_peaks_short);
    event_idx = find(strcmp({TMPEEG_peaks_short.event.type}, event_peaks));
    % save epoched, baseline-corrected data
    TMPEEG_peaks_short.setname = strjoin([set_name, "_epoched_peaks.set"],"");
    TMPEEG_peaks_short.setname = pop_saveset(TMPEEG_peaks_short, 'filename', ...
        char(TMPEEG_peaks_short.setname), 'filepath', char(savepath_baseline_peaks_short));

    % epoch for peaks & remove baseline, medium duration
    TMPEEG_peaks_medium = pop_epoch(TMPEEG, { event_peaks }, epoch_lims_peaks_medium,...
        'newname', TMPEEG.subject, 'epochinfo', 'yes');
    TMPEEG_peaks_medium = eeg_checkset( TMPEEG_peaks_medium );
    TMPEEG_peaks_medium = pop_rmbase(TMPEEG_peaks_medium, base_lims_peaks_medium);
    event_idx = find(strcmp({TMPEEG_peaks_medium.event.type}, event_peaks));
    % save epoched, baseline-corrected data
    TMPEEG_peaks_medium.setname = strjoin([set_name, "_epoched_peaks.set"],"");
    TMPEEG_peaks_medium.setname = pop_saveset(TMPEEG_peaks_medium, 'filename', ...
        char(TMPEEG_peaks_medium.setname), 'filepath', char(savepath_baseline_peaks_medium));

        % epoch for peaks & remove baseline, long duration
    TMPEEG_peaks_long = pop_epoch(TMPEEG, { event_peaks }, epoch_lims_peaks_long,...
        'newname', TMPEEG.subject, 'epochinfo', 'yes');
    TMPEEG_peaks_long = eeg_checkset( TMPEEG_peaks_long );
    TMPEEG_peaks_long = pop_rmbase(TMPEEG_peaks_long, base_lims_peaks_long);
    event_idx = find(strcmp({TMPEEG_peaks_long.event.type}, event_peaks));
    % save epoched, baseline-corrected data
    TMPEEG_peaks_long.setname = strjoin([set_name, "_epoched_peaks.set"],"");
    TMPEEG_peaks_long.setname = pop_saveset(TMPEEG_peaks_long, 'filename', ...
        char(TMPEEG_peaks_long.setname), 'filepath', char(savepath_baseline_peaks_long));



    %% artifact rejection
    %parameters: input, ICs (0) or EEG data (1), electrodes, thresholds in
    %stdv for electrodes, global threshold for all electrodes, superpose
    %prelabelling with previous labelling, keep (0) or reject (1) trials

    TMPEEG_A = pop_jointprob(TMPEEG_A, 1, 1:TMPEEG_A.nbchan, 5, 5, 0, 1);
    TMPEEG_A.comment = TMPEEG_A.comment + "  *** apply artifact rejection";
    TMPEEG_A.setname = [set_name '_complete_preprocessing_withAR_A'];
    TMPEEG_A = pop_saveset(TMPEEG_A, 'filename', TMPEEG_A.setname, 'filepath', ...
        char(output_dir_AR_const));

    TMPEEG_B = pop_jointprob(TMPEEG_B, 1, 1:TMPEEG_B.nbchan, 5, 5, 0, 1);
    TMPEEG_B.comment = TMPEEG_B.comment + "  *** apply artifact rejection";
    TMPEEG_B.setname = [set_name '_complete_preprocessing_withAR_B'];
    TMPEEG_B = pop_saveset(TMPEEG_B, 'filename', TMPEEG_B.setname, 'filepath', char(output_dir_AR_occl));

    % artifact rejection short duration
    TMPEEG_peaks_short = pop_jointprob(TMPEEG_peaks_short, 1, 1:TMPEEG_peaks_short.nbchan, 5, 5, 0 ,1);
    TMPEEG_peaks_short.comment = TMPEEG_peaks_short.comment + "  *** apply artifact rejection";
    TMPEEG_peaks_short.setname = [set_name '_complete_preprocessing_withAR_peaks'];
    TMPEEG_peaks_short = pop_saveset(TMPEEG_peaks_short, 'filename', ...
        TMPEEG_peaks_short.setname, 'filepath', char(output_dir_AR_peaks_short));

    % artifact rejection medium duration
    TMPEEG_peaks_medium = pop_jointprob(TMPEEG_peaks_medium, 1, 1:TMPEEG_peaks_medium.nbchan, 5, 5, 0 ,1);
    TMPEEG_peaks_medium.comment = TMPEEG_peaks_medium.comment + "  *** apply artifact rejection";
    TMPEEG_peaks_medium.setname = [set_name '_complete_preprocessing_withAR_peaks'];
    TMPEEG_peaks_medium = pop_saveset(TMPEEG_peaks_medium, 'filename', ...
        TMPEEG_peaks_medium.setname, 'filepath', char(output_dir_AR_peaks_medium));

    % artifact rejection long duration
    TMPEEG_peaks_long = pop_jointprob(TMPEEG_peaks_long, 1, 1:TMPEEG_peaks_long.nbchan, 5, 5, 0 ,1);
    TMPEEG_peaks_long.comment = TMPEEG_peaks_long.comment + "  *** apply artifact rejection";
    TMPEEG_peaks_long.setname = [set_name '_complete_preprocessing_withAR_peaks'];
    TMPEEG_peaks_long = pop_saveset(TMPEEG_peaks_long, 'filename', ...
        TMPEEG_peaks_long.setname, 'filepath', char(output_dir_AR_peaks_long));

end