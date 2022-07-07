% Data Analysis Pursuit-Tracking and Pursuit-Occlusion Paradigm
% Emulation pilot study 2021

% This script contains: 
% separating the data for both tasks by:
% epoching based on respective events
% export data at different steps

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

% get data paths
filepath_parts = strsplit(file_path, "\");
filepath_parts = filepath_parts(1:end-2);
parent_dir_2 = strjoin(filepath_parts, "\");

% set input & output directory
input_dir = [parent_dir_2 filesep "Emulation-Data-Output\03_parallelize_with_traj"];
output_dir = [parent_dir_2 filesep "Emulation-Data-Output\04_epoched"];

subdir_const_rand = [output_dir filesep "const_rand"];
subdir_occl = [output_dir filesep "occl_nonoccl"];
subdir_peaks = [output_dir filesep "peaks"];

cd(input_dir);

%list all *.set files in inputpath
% @saskia: add file ending produced by your script
filenames = dir('*continuous*.set');

%concatenate into one cell array
files2read = {filenames.name};

%% set parameters for peaks epoching

EPOCH_LIMS = [-.8, 0.6]; 
BASE_LIMS = [-800 -300];
EVENT = 'S 40';

%% loop through files and save epoched data
for ind = 1:length(filenames)
    
    % todo: read data and save as tmpeeg
    
    %epoch continuous data for constant and random traj and save as
    %TMPEEG_A
    TMPEEG_A = pop_epoch( TMPEEG, {'S 23' 'S 24' 'S 27'}, [-1.5 3], 'newname', [TMPEEG.setname '_A_epoched'], 'epochinfo', 'yes');
    TMPEEG_A.comment = TMPEEG_A.comment + "  *** epoch for traj parts";
    TMPEEG_A.setname = [filename '_epoched_A'];
    TMPEEG_A = pop_saveset(TMPEEG_A, 'filename', TMPEEG_A.setname, 'filepath', char(savepath_epoched));


    %epoch continuous data for occlusion/ non-occlusion and save as
    %TMPEEG_B
    TMPEEG_B = pop_epoch( TMPEEG, {'S 20' 'S 21'}, [-1.5 2], 'newname', [TMPEEG.setname '_B_epoched'], 'epochinfo', 'yes');
    TMPEEG_B.comment = TMPEEG_B.comment + "  *** epoch for occl parts";
    TMPEEG_B.setname = [filename '_epoched_B'];
    TMPEEG_B = pop_saveset(TMPEEG_B, 'filename', TMPEEG_B.setname, 'filepath', char(savepath_epoched));

    % epoch for peaks & remove baseline for these epochs (other are
    % baseline corrected in the next script)
    EEG_peaks = pop_epoch( EEG, {  EVENT  }, EPOCH_LIMS, 'newname', EEG.subject, 'epochinfo', 'yes');
    EEG_peaks = eeg_checkset( EEG_peaks );
    EEG = pop_rmbase(EEG, BASE_LIMS);
    event_idx = find(strcmp({EEG.event.type}, EVENT));
    % save epoched data
    file_name = strcat([EEG.subject, "_epoched_peaks"]);
    EEG = pop_saveset(EEG, 'filename', file_name, 'filepath', char(subdir_peaks));
end