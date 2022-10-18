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
parent_dir = filepath_parts(1:end-1);
parent_dir = strjoin(parent_dir, filesep);
parent_dir_2 = filepath_parts(1:end-2);
parent_dir_2 = strjoin(parent_dir_2, filesep);
output_dir = strjoin([parent_dir_2, "Emulation-Data-Output"], filesep);
input_dir = strjoin([parent_dir_2, "Emulation-Data-Output", "07_epochs_with_extra_fields"], filesep);
loreta_file = strjoin([parent_dir_2, "Emulation-Data-Output", "loreta-ready_files"], filesep);
mean_matrices_path = strjoin([output_dir, 'mean_matrices'], filesep);
mean_matrices_peaks_epoched_path = strjoin([mean_matrices_path, 'peaks'], filesep);
condition_paths = strjoin([parent_dir_2, "Emulation-Data-Output", "Condition-Separated-EEGLAB-sets"], filesep)


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
%         if ~(strcmp(EEG.event(ev-del).type, 'S 40') | strcmp(EEG.event(ev-del).type, 'S 50'))
%             EEG.event(ev-del) = [];
%             del = del +1;
%         end
    end
    log_idx = unfind(event_idx, length(EEG.event));
    EEG.event(~log_idx) = [];
    EEG = eeg_checkset(EEG);
    ALLEEG(s) = EEG;
    s = s + 1;
end
% WATCH OUT: FROM NOW ON get_epoch_center won't work anymore since the
% EEG.epoch structure is no longer aligned with the EEG.event structure. 


%% Export Data for Brain irgendwas. 

for s = 1:length(ALLEEG)
    EEG = ALLEEG(s);
    EEG = pop_selectevent( EEG, 'artifact_error',{'VALID'},'deleteevents',...
        'off','deleteepochs','on','invertepochs','off');
    ALLEEG(s) = EEG;
end

% remove occluded VISIBLE
mkdir(strjoin([condition_paths, "vis"], filesep))

CONDEEG = ALLEEG;
for s = 1:length(CONDEEG)
    EEG = CONDEEG(s);
    EEG = pop_selectevent( EEG, 'OCCL',{'OFF'},'deleteevents',...
        'off','deleteepochs','on','invertepochs','off');
    file_name = strjoin([EEG.subject, "vis_only.set"], '_');
    file_path = strjoin([condition_paths, "vis"], filesep);
    pop_saveset(EEG, 'filename', char(file_name), 'filepath', char(file_path));
    CONDEEG(s) = EEG;
end

mkdir(strjoin([condition_paths, "occl"], filesep))

CONDEEG2 = ALLEEG;
% remove visible OCCLUDED
for s = 1:length(CONDEEG2)
    EEG = CONDEEG2(s);
    EEG = pop_selectevent( EEG, 'OCCL',{'ON'},'deleteevents',...
        'off','deleteepochs','on','invertepochs','off');
        file_name = strjoin([EEG.subject, "occ_only.set"], '_');
    file_path = strjoin([condition_paths, "occl"], filesep);
    pop_saveset(EEG, 'filename', char(file_name), 'filepath', char(file_path));
    CONDEEG2(s) = EEG;
end

mkdir(strjoin([condition_paths, "const"], filesep))

clear CONDEEG CONDEEG2
% remove random CONSTANT
CONDEEG = ALLEEG;
for s = 1:length(CONDEEG)
    EEG = CONDEEG(s);
    EEG = pop_selectevent( EEG, 'TRAJ',{'RANDOM1', 'RANDOM2'},'deleteevents',...
        'off','deleteepochs','on','invertepochs','off');
        file_name = strjoin([EEG.subject, "const_only.set"], '_');
    file_path = strjoin([condition_paths, "const"], filesep);
    pop_saveset(EEG, 'filename', char(file_name), 'filepath', char(file_path));
    CONDEEG(s) = EEG;
end

mkdir(strjoin([condition_paths, "rand1"], filesep))

CONDEEG2 = ALLEEG;
% remove constant RANDOM1
for s = 1:length(CONDEEG2)
    EEG = CONDEEG2(s);
    EEG = pop_selectevent( EEG, 'TRAJ',{'RANDOM2', 'CONST'},'deleteevents',...
        'off','deleteepochs','on','invertepochs','off');
        file_name = strjoin([EEG.subject, "rand1_only.set"], '_');
    file_path = strjoin([condition_paths, "rand1"], filesep);
    pop_saveset(EEG, 'filename', char(file_name), 'filepath', char(file_path));
    CONDEEG2(s) = EEG;
end


%% Actual Data Export

cd(mean_matrices_peaks_epoched_path)

load('mean_struct.mat')

cd(loreta_file)
% make sure files are saved as ASCII
feature('DefaultCharacterSet', 'ASCII');

cond_labs = {'vis', 'occ', 'rand1', 'const'};
% delete channels that are non-standard

% delete the following channels: 
del_chans = {'O9', 'O10', 'P11', 'P12'};

% get index of useless channels. 
for chan = 1:length(del_chans)
    del_idx(chan) = find(strcmp({ALLEEG(1).chanlocs.labels}, del_chans{chan}));
end
% deletion from mean structure
for cond = 1:length(cond_labs)
    mean_struct.(cond_labs{cond})(del_idx,:,:) = [];
end

% Delection from ALLEEG structure
ALLEEG(1).chanlocs(del_idx) = [];
chanlocs = ALLEEG(1).chanlocs;

% export all data
for s = 1:size(mean_struct.all,3)
    for cond = 1:length(cond_labs)
        EEG = ALLEEG(s);
        eeglab2loreta(chanlocs, mean_struct.(cond_labs{cond})(:,:,s), 'exporterp', 'on');
        % rename the new file because the output is stupid and I cannot define
        % a name in the eeglab2loreta function.
        new_filename = strjoin([EEG.subject, cond_labs{cond}, "loreta", "ERP.txt"], '_');
        movefile("erp.txt", new_filename);
    end
end

% export only a time window
SR = 250;
time_window = round([0.175, 0.225]*SR);

for s = 1:size(mean_struct.all,3)
    for cond = 1:length(cond_labs)
        EEG = ALLEEG(s);
        eeglab2loreta(chanlocs, mean_struct.(cond_labs{cond})(:,time_window,s), 'exporterp', 'on');
        % rename the new file because the output is stupid and I cannot define
        % a name in the eeglab2loreta function.
        new_filename = strjoin([EEG.subject, cond_labs{cond}, "loreta", "ERP_only_P2.txt"], '_');
        movefile("erp.txt", new_filename);
    end
end

% save the time vector as well. 
writematrix(EEG.times', 'time_vector.txt');


