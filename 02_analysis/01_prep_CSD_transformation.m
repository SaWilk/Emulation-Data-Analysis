% CSD transformation

% this script calculated the surface laplacians for the eeg data and saves
% the transformed ERP datasets

% author: Saskia Wilken
% creation date: 05.09.2022 


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
mean_matrices_path = strjoin([output_dir, 'mean_matrices'], filesep);
mkdir(mean_matrices_path);
mean_matrices_peaks_epoched_path = strjoin([mean_matrices_path, 'peaks'], filesep);
mkdir(mean_matrices_peaks_epoched_path);

% output
csd_dir = strjoin([output_dir, "csd_transform"], filesep);
% input
epochs_plus_error_dir = strjoin([output_dir, "07_epochs_with_extra_fields"], filesep);


%% Load data for CSD transformation

cd(epochs_plus_error_dir);
%list all *.set files in inputpath
file_names = dir('*.set');
%concatenate into one cell array
files2read = {file_names.name};
eeglab;
close all

for idx = 1:length(files2read)

    EEG = pop_loadset('filename',files2read{idx});
    [ALLEEG, EEG] = eeg_store(ALLEEG, EEG);

end


%% Surface Laplacians CSD toolbox (resp. current source density)

% keep default settings
addpath(genpath('C:\wilken\CSDtoolbox')); % CSD toolbox needs to be installed for this to work

% tutorial: https://psychophysiology.cpmc.columbia.edu/Software/CSDtoolbox/tutorial.html
% common errors: https://psychophysiology.cpmc.columbia.edu/Software/CSDtoolbox/errors.html
% FAQ: https://psychophysiology.cpmc.columbia.edu/Software/CSDtoolbox/FAQ.html

% exploit the EEG montage information included in any EEGlab data file to 
% generate the montage-dependent "Channels × Channels" transformation 
% matrices G and H (EEGlab_Make_G_H.m), and how to use these 
% transformation matrices with individual EEGlab data files

% Get usable list of electrodes from EEGlab data structure

% instead of using the custom electrode position motage, I am renaming the
% electrodes in ALLEEG.chanlocs. 
for site = 1:length({ALLEEG(1).chanlocs.labels})
    trodes{site}=(ALLEEG(1).chanlocs(site).labels);
end
trodes=trodes';

figure; 
topoplot([],ALLEEG(1).chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo);

% Get Montage for use with CSD Toolbox
% '10-5-System_Mastoids_EGI129.csd' is a data file with electrode locations
% according to the 10-5 system. I manually added the strange electrodes 
% P11, P12, O9 and O10. However, the final electrode locations still look a
% bit fishy. A problem is that EEGLAB apparently uses a different scaling
% for the X and Y coordinates, which makes the approach above not very
% helpful
Montage_64=ExtractMontage('10-5-System_Mastoids_EGI129_SW.csd',trodes); % _SW is my custom version

MapMontage(Montage_64); % a wonky eeg location plot. 


%% Derive G and H!
[G,H] = GetGH(Montage_64); % default m = 4


%% Save G and H to later import when doing the CSD transform on files

% revised method to store G and H matrices with CSD montage for later import
Montage = Montage_64;                             % use generic variable name
save(strjoin([csd_dir,'CSD_montage_64.mat'], filesep), 'G', 'H', 'Montage'); % save variables to Matlab file
clear G H Montage;                                % remove variables from workspace
load(strjoin([csd_dir,'CSD_montage_64.mat'], filesep)); % restore variables to the workspace


%% Calculate the CSD

[~, base_end] = min(abs(ALLEEG(1).times)); % find end of baseline period
for s = 1:length(ALLEEG)
    tic
    for ne = 1:length(ALLEEG(s).epoch)               % loop through all epochs
        myEEG = single(ALLEEG(s).data(:,:,ne));      % reduce data precision to reduce memory demand
        MyResults = CSD(myEEG,G,H,1.0e-5, 9);        % compute CSD for <channels-by-samples> 2-D epoch
        % default lambda but 9 as head radius (slightly larger than 56
        % cm circumference)
        data(:,:,ne) = MyResults;              % assign data output
    end
    ALLEEG(s).("CSD_data") = double(data);          % final CSD data
    set_name = strsplit(ALLEEG(s).setname, '_');
    new_set_name = strjoin([set_name{1}, "CSD"], '_');
    ALLEEG(s).setname = new_set_name;
    clear data myEEG MyResults new_set_name
    ALLEEG(s) = pop_saveset(ALLEEG(s), 'filename', char(ALLEEG(s).setname), 'filepath', ...
        char(csd_dir));
    toc
end


%% Calculate means and normalize data

clear z_mean_struct mean_struct z_epoch_mean cond_ind cond_labs epoch_weights count peaks

% remove the useless subject
filenames = {ALLEEG.filename};
count = 1;
for i = 1:length(filenames)
    if ~contains(filenames{i}, 'KMY6K')
    TMPEEG(count) = ALLEEG(i);
    count = count +1;
    end
end

clear ALLEEG
ALLEEG = TMPEEG;
clear TMPEEG

% conditions labels for looping
cond_labs = {'all', 'occ', 'vis', 'rand1', 'const', 'rand2'};
% count the peaks of the data
for lab = 1:length(cond_labs)
    count_peaks.(cond_labs{lab}) = 0;
end

count_peaks.diff_occ = [];
count_peaks.diff_const_rand1 = [];
count_peaks.diff_const_rand2 = [];
count_peaks.diff_rand1_rand2 = [];

for s = 1:size(ALLEEG, 2)

    EEG = ALLEEG(s);
    % remove epochs in first 500 ms area of each trial

    % get indices of peaks in event structure that also belong to a certain
    % condition
    for lab = 1:length(cond_labs)
        log_vecs.(cond_labs{lab}) = false(size(EEG.event,2),1);
    end

    for ev = 1:size(EEG.event,2)
        log_vecs.(cond_labs{1})(ev) = (strcmp(EEG.event(ev).type, 'S 40') | ...
            strcmp(EEG.event(ev).type, 'S 50'));
    end

    for ev = 1:size(EEG.event,2)
        % occlusion
        log_vecs.(cond_labs{2})(ev) = log_vecs.(cond_labs{1})(ev) && ...
            strcmp(EEG.event(ev).task, 'task_b') && ...
            strcmp(EEG.event(ev).OCCL, 'ON');
        log_vecs.(cond_labs{3})(ev) = log_vecs.(cond_labs{1})(ev) && ...
            strcmp(EEG.event(ev).task, 'task_b') && ...
            strcmp(EEG.event(ev).OCCL, 'OFF');
        % constant
        log_vecs.(cond_labs{4})(ev) = log_vecs.(cond_labs{1})(ev) && ...
            strcmp(EEG.event(ev).task, 'task_a') && ...
            strcmp(EEG.event(ev).TRAJ, 'RANDOM1');
        log_vecs.(cond_labs{5})(ev) = log_vecs.(cond_labs{1})(ev) && ...
            strcmp(EEG.event(ev).task, 'task_a') && ...
            strcmp(EEG.event(ev).TRAJ, 'CONST');
        log_vecs.(cond_labs{6})(ev) = log_vecs.(cond_labs{1})(ev) && ...
            strcmp(EEG.event(ev).task, 'task_a') && ...
            strcmp(EEG.event(ev).TRAJ, 'RANDOM2');
    end

    % get epoch numbers for respective conditions
    for lab = 1:length(cond_labs)
        % get the epochs that satisfy the condition conditions
        relevant_epochs = unique([ ...
            EEG.event(log_vecs.(cond_labs{lab})).epoch]);
        count = 1;
        % initialize structure field of unknown size
        cond_ind(s).(cond_labs{lab}) = [];
        for ep = relevant_epochs
            % get epochs of the zero latency events that match the
            % conditions
            event_idx = get_epoch_center(EEG, ep);
            % if this is not in the artifactual range of the epochs...
            if  strcmp(EEG.event(event_idx).artifact_error, 'VALID')
                % store the epoch number in cond_inds
                cond_ind(s).(cond_labs{lab})(count) = EEG.event(event_idx).epoch;
                count = count + 1;
            end
        end
        % calculate mean of EEG data across condition
        mean_struct.(cond_labs{lab})(:,:,s) = mean(...
            EEG.CSD_data(:, :, cond_ind(s).(cond_labs{lab})), 3);
        % get z-normalized data
        z_mean_struct.(cond_labs{lab})(:, :, s) = normalize( ...
            mean_struct.(cond_labs{lab})(:, :, s),2, 'zscore');
        % get weights
        z_mean_struct.(strjoin(["num_epochs_", cond_labs{lab}],""))(s) = size( ...
            EEG.CSD_data(:,:,cond_ind(s).(cond_labs{lab})), 3);
    end

    % calculate contrasts by subtraction
    mean_struct.diff_occ_vis(:,:,s) = mean_struct.occ(:,:,s) - mean_struct.vis(:,:,s);
    mean_struct.diff_const_rand1(:,:,s) = mean_struct.const(:,:,s) - mean_struct.rand1(:,:,s);
    mean_struct.diff_const_rand2(:,:,s) = mean_struct.const(:,:,s) - mean_struct.rand2(:,:,s);
    mean_struct.diff_rand1_rand2(:,:,s) = mean_struct.rand1(:,:,s) - mean_struct.rand2(:,:,s);

    z_mean_struct.diff_occ_vis(:,:,s) = z_mean_struct.occ(:,:,s) - z_mean_struct.vis(:,:,s);
    z_mean_struct.diff_const_rand1(:,:,s) = z_mean_struct.const(:,:,s) - z_mean_struct.rand1(:,:,s);
    z_mean_struct.diff_const_rand2(:,:,s) = z_mean_struct.const(:,:,s) - z_mean_struct.rand2(:,:,s);
    z_mean_struct.diff_rand1_rand2(:,:,s) = z_mean_struct.rand1(:,:,s) - z_mean_struct.rand2(:,:,s);

    % get time vec
    mean_struct.time_vec = EEG.times;
    z_mean_struct.time_vec = EEG.times;
end


%% Create the weights
for lab = 1:length(cond_labs)
    % divide the respective conditions for each subject by the total number
    % of trials that subject had. The resulting porportion is your weight. 
    epoch_weights.(cond_labs{lab}) = median(z_mean_struct.(...
        strjoin(["num_epochs_", cond_labs{lab}],""))) ./ z_mean_struct.(...
        strjoin(["num_epochs_", cond_labs{lab}],""));
end

%% Weigh the data according to trials per subject
for lab = 1:length(cond_labs)
    for s = 1:size(z_mean_struct.all, 3)
        z_mean_struct.(cond_labs{lab})(:,:,s) = z_mean_struct.(...
            cond_labs{lab})(:,:,s) * epoch_weights.(cond_labs{lab})(s);
        mean_struct.(cond_labs{lab})(:,:,s) = mean_struct.(...
            cond_labs{lab})(:,:,s) * epoch_weights.(cond_labs{lab})(s);
    end
end


%% Save the Mean Matrices for Convenient Loading

cd(mean_matrices_peaks_epoched_path)
save("mean_struct_csd.mat", 'mean_struct')
save("condition_indices_csd.mat", 'cond_ind')
save("z_mean_struct_csd.mat", 'z_mean_struct')
