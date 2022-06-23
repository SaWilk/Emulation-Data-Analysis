%% ERP Analysis

% Creation date: 22.06.2022
% Author: Saskia Wilken

% This script performs ERP analysis on the data output by
% merge_eeg_and_tracking_data.m which should be located in the same file. 

%% Empty

format compact
format long G
clear
clc


%% Set Paths

% file location
file_path = matlab.desktop.editor.getActiveFilename;
file_path = fileparts(file_path);
% tracking data
track_data_path = strcat(file_path, "\00_npz_files"); % get root path to data folder
subj_paths = genpath(track_data_path); % generate subfolder paths
subj_paths = strsplit(subj_paths, ";"); % split the string vector at ;

subj_paths(1) = []; % remove top level entry
subj_paths(end) = []; % remove final entry

% eeg data
eeg_data_path = strcat(file_path, "\eeg_files"); % get root path to data folder
% contains raw and preprocessed folders

% get subject ids from folder names
[~, tmp] = fileparts(subj_paths);
subj_ids = split(tmp, '_');
subj_ids = subj_ids(:, :, 1);
subj_ids = unique(subj_ids);
num_subj = length(subj_ids); % get num of subfolders

% initialize data structure
track_data = struct();
clear tmp
% eeglab



%% Epoch around Peaks

s = 1
TMPEEG = all_data_struct(s)
TMPEEG.event = all_data_struct(s).track_event_peaks
TMPEEG = pop_epoch( TMPEEG, {'S 40'}, [-0.25 0.75], 'newname', [TMPEEG.setname '_A_peak_epoched'], 'epochinfo', 'yes');

% this needs lots of improvement still....
pop_erpimage(TMPEEG, 1)
figure;
pop_erpimage(EEG,1, [1],[[]],'Cz',10,1,{},[],'' ,'yerplabel',...
    '\muV','erp','on','cbar','on','topo', { [1] EEG.chanlocs EEG.chaninfo } );
plottopo

% 60 channels because...
% ground wird gleich reingferechnet
% drei EOG-Elektroden, die aber nicht verwendet werden, werden nicht
% aufgezeichnet
% reference ist Fpz
