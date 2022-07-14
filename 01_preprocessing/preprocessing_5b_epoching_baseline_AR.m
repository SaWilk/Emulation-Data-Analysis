% Data Analysis Pursuit-Tracking and Pursuit-Occlusion Paradigm
% Emulation pilot study 2021

% This script contains: 
% * epoching data for whole trial (S26 + 9s)
% * removing baseline of 1s (?) before fixation cross
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

cd(input_dir);

%list all *.set files in inputpath
file_names = dir('*_EEG.set');

%concatenate into one cell array
files2read = {file_names.name};


%% read files, 

for ind = 1%:length(file_names)
    
    file_name = files2read{ind};
    TMPEEG = pop_loadset(file_name);
    set_name = TMPEEG.subject;
    
    %epoch continuous data for constant and random traj and save as
    %TMPEEG_A
    TMPEEG_A = pop_epoch( TMPEEG, {"S 26"}, [-750 9000], 'newname', ...
        [TMPEEG.setname '_epoched_trials'], 'epochinfo', 'no');
    TMPEEG_A.comment = TMPEEG_A.comment + "  *** epoch for trials";
    TMPEEG_A.setname = [set_name '_epoched_trials.set'];
    %TMPEEG_A = pop_saveset(TMPEEG_A, 'filename', TMPEEG_A.setname, 'filepath', ...
       % char(output_dir_epoched));

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
% mean(diff);