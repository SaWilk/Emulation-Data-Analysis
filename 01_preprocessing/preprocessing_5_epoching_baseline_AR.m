% Data Analysis Pursuit-Tracking and Pursuit-Occlusion Paradigm
% Emulation pilot study 2021

% This script contains: 
% * separating the data for both tasks by:
% epoching based on respective events
% * baseline correction
% * remove bad epochs
% * export data at different steps

% Adriana Böttcher
% 08.07.22

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
filepath_parts = strsplit(file_path, "\");
parent_dir = filepath_parts(1:end-1);
parent_dir = strjoin(parent_dir, "\");
parent_dir_2 = filepath_parts(1:end-2);
parent_dir_2 = strjoin(parent_dir_2, "\");

% add functions in parent dir
addpath([char(parent_dir) filesep char("functions")]);

% set input & output directory
input_dir = [parent_dir_2 filesep "Emulation-Data-Output\03_parallelize_with_traj"];
output_dir = [parent_dir_2 filesep "Emulation-Data-Output\"];
output_dir_epoched = [parent_dir_2 filesep "Emulation-Data-Output\04_epoched"];
output_dir_baseline = [parent_dir_2 filesep "Emulation-Data-Output\05_baseline"];
output_dir_AR = [parent_dir_2 filesep "Emulation-Data-Output\06_artifact_rejection"];

subdir_const_rand = [output_dir_epoched filesep "const_rand"];
subdir_occl = [output_dir_epoched filesep "occl_nonoccl"];
subdir_peaks = [output_dir_epoched filesep "peaks"];

cd(input_dir);

%list all *.set files in inputpath
% @saskia: add file ending produced by your script
%filenames = dir('**.set');

%concatenate into one cell array
files2read = {filenames.name};

%% set parameters for epoching

epoch_lims_A = [-1.5 3];
epoch_lims_B =  [-1.5 2];
epoch_lims_peaks = [-.8, 0.6]; 

event_A = {'S 23' 'S 24' 'S 27'};
event_B = {'S 20' 'S 21'};
event_peaks = 'S 40';

base_lims_A = [-1000 0];
base_lims_B = [-1000 0];
base_lims_peaks = [-800 -300];

%% loop through files and save epoched data
for ind = 1:length(filenames)
    
    % todo: read data and save as tmpeeg
    
    %% epoching & remove baseline for A, B, peaks

    %epoch continuous data for constant and random traj and save as
    %TMPEEG_A
    TMPEEG_A = pop_epoch( TMPEEG, event_A , epoch_lims_A, 'newname', [TMPEEG.setname '_A_epoched'], 'epochinfo', 'yes');
    TMPEEG_A.comment = TMPEEG_A.comment + "  *** epoch for traj parts";
    TMPEEG_A.setname = [filename '_epoched_A'];
    TMPEEG_A = pop_saveset(TMPEEG_A, 'filename', TMPEEG_A.setname, 'filepath', char(subdir_const_rand));
    % remove baseline
    TMPEEG_A = pop_rmbase(TMPEEG_A, base_lims_A ,[]);
    TMPEEG_A.comment = TMPEEG_A.comment + "  *** apply baseline corr";
    TMPEEG_A.setname = [filename '_complete_preprocessing_A'];
    TMPEEG_A = pop_saveset(TMPEEG_A, 'filename', TMPEEG_A.setname, 'filepath', char(savepath_baseline));

    %epoch continuous data for occlusion/ non-occlusion and save as
    %TMPEEG_B
    TMPEEG_B = pop_epoch(TMPEEG, event_B , epoch_lims_B, 'newname', [TMPEEG.setname '_B_epoched'], 'epochinfo', 'yes');
    TMPEEG_B.comment = TMPEEG_B.comment + "  *** epoch for occl parts";
    TMPEEG_B.setname = [filename '_epoched_B'];
    TMPEEG_B = pop_saveset(TMPEEG_B, 'filename', TMPEEG_B.setname, 'filepath', char(subdir_occl));
    % remove baseline
    TMPEEG_B = pop_rmbase(TMPEEG_B, base_lims_B ,[]);
    TMPEEG_B.comment = TMPEEG_B.comment + "  *** apply baseline corr";
    TMPEEG_B.setname = [filename '_complete_preprocessing_B'];
    TMPEEG_B = pop_saveset(TMPEEG_B, 'filename', TMPEEG_B.setname, 'filepath', char(savepath_baseline));

    % epoch for peaks & remove baseline 
    TMPEEG_peaks = pop_epoch(TMPEEG, { event_peaks }, epoch_lims_peaks, 'newname', TMPEEG.subject, 'epochinfo', 'yes');
    TMPEEG_peaks = eeg_checkset( TMPEEG_peaks );
    TMPEEG_peaks = pop_rmbase(TMPEEG_peaks, base_lims_peaks);
    event_idx = find(strcmp({TMPEEG_peaks.event.type}, event_peaks));
    % save epoched, baseline-corrected data
    file_name = strcat([EEG.subject, "_epoched_peaks"]);
    EEG = pop_saveset(EEG, 'filename', file_name, 'filepath', char(subdir_peaks));

    %% artifact rejection
    %parameters: input, ICs (0) or EEG data (1), electrodes, thresholds in
    %stdv for electrodes, global threshold for all electrodes, superpose
    %prelabelling with previous labelling, keep (0) or reject (1) trials

    TMPEEG_A = pop_jointprob(TMPEEG_A, 1, 1:TMPEEG_A.nbchan, 5, 5, 0, 1);
    TMPEEG_A.comment = TMPEEG_A.comment + "  *** apply artifact rejection";
    TMPEEG_A.setname = [filename '_complete_preprocessing_withAR_A'];
%     TMPEEG_A = pop_saveset(TMPEEG_A, 'filename', TMPEEG_A.setname, 'filepath', char(output_dir_AR));

    TMPEEG_B = pop_jointprob(TMPEEG_B, 1, 1:TMPEEG_B.nbchan, 5, 5, 0, 1);
    TMPEEG_B.comment = TMPEEG_B.comment + "  *** apply artifact rejection";
    TMPEEG_B.setname = [filename '_complete_preprocessing_withAR_B'];
%     TMPEEG_B = pop_saveset(TMPEEG_B, 'filename', TMPEEG_B.setname, 'filepath', char(output_dir_AR));

    TMPEEG_peaks = pop_jointprob(TMPEEG_peaks, 1, 1:TMPEEG_peaks.nbchan, 5, 5, 0 ,1);
    TMPEEG_peaks.comment = TMPEEG_peaks.comment + "  *** apply artifact rejection";
    TMPEEG_peaks.setname = [filename '_complete_preprocessing_withAR_peaks'];
%     TMPEEG_peaks = pop_saveset(TMPEEG_peaks, 'filename', TMPEEG_B.setname, 'filepath', char(output_dir_AR));

end