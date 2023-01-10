%% Transform EEGLAB to LORETA-suitable data
% use the eeglab plugin LORETA to save eeglab data suitable for the LORETA
% standalone software that needs to be installed. Follow instrucitons in
% this video: 
% https://www.youtube.com/watch?v=amttvN_Sb6A
% keep in mind that you need non-CSD-transformed data for that

% Author: Saskia Wilken
% Creation date: 09.10.2022


%% empty everything

format compact
format long G
clear
clc
close all


%% Define Paths

% file location
file_path = matlab.desktop.editor.getActiveFilename;
file_path = fileparts(file_path);
addpath(file_path);

% get data paths for parent dirs
filepath_parts = strsplit(file_path, filesep);
parent_dir_2 = filepath_parts(1:end-2);
parent_dir_2 = strjoin(parent_dir_2, filesep);
output_dir = strjoin([parent_dir_2, "Emulation-Data-Output"], filesep);
input_dir = strjoin([parent_dir_2, "Emulation-Data-Output", "07_epochs_with_extra_fields"], filesep);
loreta_file = strjoin([parent_dir_2, "Emulation-Data-Output", "loreta-ready_files"], filesep);
mean_matrices_path = strjoin([output_dir, 'mean_matrices'], filesep);
mean_matrices_peaks_epoched_path = strjoin([mean_matrices_path, 'peaks'], filesep);


%% Load non-transformed data for creating one big dataset out of all the data

eeglab;

% long epochs
cd(input_dir);
%list all *.set files in inputpath
file_names = dir('*.set');
%concatenate into one cell array
files2read = {file_names.name};
files2read = files2read(~contains(files2read, 'KMY6K')); % remove subject who 
% did not properly track during task B
% load eeg data epoched around peaks
for idx = 1:length(files2read)
    ALLEEG2(idx) = pop_loadset('filename',files2read{idx});
end

ALLEEG = ALLEEG2;
clear ALLEEG2

%% Delete all events except for the Peak events so there is only one per epoch

clear z_ALLEEG TMPEEG z_mean_struct mean_struct z_tmpEEG CONDEEG CONDEEG2 TMPEEG2
sum(arrayfun(@(s) length(s.epoch), ALLEEG)); % how many epochs we need. 
s = 1;
while s <= length(ALLEEG)
    EEG = ALLEEG(s);
    del = 0;
    event_idx = zeros(length(EEG.epoch),1);
    for ep = 1:length(EEG.epoch)
        event_idx(ep) = get_epoch_center(EEG, ep);
    end
    log_idx = unfind(event_idx, length(EEG.event));
    EEG.event(~log_idx) = [];
    EEG = eeg_checkset(EEG);
    ALLEEG(s) = EEG;
    s = s + 1;
end
% WATCH OUT: FROM NOW ON get_epoch_center won't work anymore since the
% EEG.epoch structure is no longer aligned with the EEG.event structure. 


%% Data Export

cd(mean_matrices_peaks_epoched_path)

load('mean_struct.mat')

% file and folder names
cd(loreta_file)
file_names_out = {"complete_occvis"};

% define time intervals of clusters
SR = 250;
time_windows = [0.004, 0.748]*SR;
condition_cell = {"occ", "vis"};

% make sure files are saved as ASCII
feature('DefaultCharacterSet', 'ASCII');

cond_labs = {'vis', 'occ', 'rand1', 'const'};

% delete channels that are non-standard
del_chans = {'O9', 'O10', 'P11', 'P12'};

% get index of useless channels. 
for chan = 1:length(del_chans)
    del_idx(chan) = find(strcmp({ALLEEG(1).chanlocs.labels}, del_chans{chan}));
end
% deletion from mean structure
for cond = 1:length(cond_labs)
    mean_struct.(cond_labs{cond})(del_idx,:,:) = [];
end

% Deletion from ALLEEG structure
ALLEEG(1).chanlocs(del_idx) = [];
chanlocs = ALLEEG(1).chanlocs;

% export files adapted for sLORETA
for cont = 1:size(condition_cell, 1)
    cur_dir = strjoin([loreta_file, file_names_out{cont}], "\");
    mkdir(cur_dir);
    cd(cur_dir);
    for cond = 1:size(condition_cell(cont, :), 2)
        for s = 1:size(mean_struct.(condition_cell{cont, cond}),3)
            eeglab2loreta(chanlocs, mean_struct.(condition_cell{cont, cond})(:, ...
                time_windows(cont, 1):time_windows(cont, 2), s), 'exporterp', 'on');
            % rename the new file because I cannot define
            % a name in the eeglab2loreta function.
            new_filename = strjoin([ALLEEG(s).subject, "Cluster", file_names_out{cont}, ...
                "Condition", condition_cell{cont, cond}, "Loreta-ready", ...
                "ERP", "Time-window", time_windows(cont, 1)*1000/SR,"-", ...
                time_windows(cont, 2)*1000/SR, ".txt"], '_');
            movefile("erp.txt", new_filename);
        end
    end
end

% save the time vector as well.
writematrix(EEG.times', 'time_vector.txt');


