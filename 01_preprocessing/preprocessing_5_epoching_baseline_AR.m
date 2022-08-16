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
subdir_peaks = strjoin([output_dir_epoched filesep 'peaks_-500_750'], filesep);
mkdir(subdir_peaks)


savepath_baseline_const_rand = strjoin([output_dir_baseline, "const_rand"], filesep);
savepath_baseline_occl = strjoin([output_dir_baseline, 'occl_nonoccl'], filesep);
% out dirs for peaks for multiple epoch lengths
savepath_baseline_peaks = strjoin([output_dir_baseline, 'peaks_-500_750'], filesep);
mkdir(savepath_baseline_peaks)


output_dir_AR_const = strjoin([output_dir_AR, 'const'], filesep);
output_dir_AR_occl = strjoin([output_dir_AR, 'occl'], filesep);
% out dirs for peaks for multiple epoch lengths
output_dir_AR_peaks = strjoin([output_dir_AR, 'peaks_-500_750'], filesep);
mkdir(output_dir_AR_peaks)


cd(input_dir);

% list all *.set files in inputpath
file_names = dir('*parallelized.set');

% concatenate into one cell array
files2read = {file_names.name};


%% set parameters for epoching

epoch_lims_A = [-1.5 3];
epoch_lims_B =  [-1.5 2];
epoch_lims_peaks = [-0.5, 0.75]; 


event_A = {'S 23' 'S 24' 'S 27'};
event_B = {'S 20' 'S 21'};
event_peaks = 'S 40';

base_lims_A = [-1000 0];
base_lims_B = [-1000 0];
base_lims_peaks = [-500 0];


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

        % epoch for peaks & remove baseline, long duration
    TMPEEG_peaks = pop_epoch(TMPEEG, { event_peaks }, epoch_lims_peaks,...
        'newname', TMPEEG.subject, 'epochinfo', 'yes');
    TMPEEG_peaks = eeg_checkset( TMPEEG_peaks );
    TMPEEG_peaks = pop_rmbase(TMPEEG_peaks, base_lims_peaks);
    event_idx = find(strcmp({TMPEEG_peaks.event.type}, event_peaks));
    % save epoched, baseline-corrected data
    TMPEEG_peaks.setname = strjoin([set_name, "_epoched_peaks.set"],"");
    TMPEEG_peaks.setname = pop_saveset(TMPEEG_peaks, 'filename', ...
        char(TMPEEG_peaks.setname), 'filepath', char(savepath_baseline_peaks));


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

    % artifact rejection long duration
    TMPEEG_peaks = pop_jointprob(TMPEEG_peaks, 1, 1:TMPEEG_peaks.nbchan, 5, 5, 0 ,1);
    TMPEEG_peaks.comment = TMPEEG_peaks.comment + "  *** apply artifact rejection";
    TMPEEG_peaks.setname = [set_name '_complete_preprocessing_withAR_peaks'];
    TMPEEG_peaks = pop_saveset(TMPEEG_peaks, 'filename', ...
        TMPEEG_peaks.setname, 'filepath', char(output_dir_AR_peaks));

end

% TODO: Add study, for group analyses