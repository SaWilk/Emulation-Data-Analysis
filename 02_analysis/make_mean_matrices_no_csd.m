%% Create Mean Matrices of Data without CSD Transformation before. 

% Author: Saskia Wilken
% creation Date: 28.06.2022


%% Empty

format compact
format long G
clear
clc


%% Folders and Dependencies

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

output_dir = strjoin([parent_dir_2, "Emulation-Data-Output"], filesep);

epochs_plus_error_dir = strjoin([output_dir, "07_epochs_with_extra_fields"], filesep);

mean_matrices_path = strjoin([output_dir, 'mean_matrices'], filesep);
mkdir(mean_matrices_path);
mean_matrices_peaks_epoched_path = strjoin([mean_matrices_path, 'peaks'], filesep);
mkdir(mean_matrices_peaks_epoched_path);

peak_erp_plots_dir = strjoin([parent_dir, "plots", "erp_plots", "peaks"], filesep);
mkdir(peak_erp_plots_dir);

% neighbors dir
neighbors_dir = strjoin([parent_dir_2, "Emulation-Data-Input", "EEG_files"], filesep);


%% Load Epoched data

eeglab;

% long epochs
cd(epochs_plus_error_dir);
%list all *.set files in inputpath
file_names = dir('*.set');
%concatenate into one cell array
files2read = {file_names.name};
% load eeg data epoched around peaks
for idx = 1:length(files2read)
    ALLEEG2(idx) = pop_loadset('filename',files2read{idx});
end

ALLEEG = ALLEEG2;
clear ALLEEG2


%% Calculate Means Across Subjects
 
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
%
for s = 1:size(ALLEEG, 2)

    EEG = ALLEEG(s);
    % remove epochs in artifactual area of each trial

    % get indices of peaks in event structure that also belong to a certain
    % condition
    for lab = 1:length(cond_labs)
        log_vecs.(cond_labs{lab}) = false(size(EEG.event,2),1);
    end

    for ev = 1:size(EEG.event,2)
        log_vecs.(cond_labs{1})(ev) = (strcmp(EEG.event(ev).type, 'S 40') | ...
            strcmp(EEG.event(ev).type, 'S 50'));
    end
%     test = cond_ind.(cond_labs{5});
%     for i = 1:length([EEG.event])
%         EEG.event(i).('log') = test(i);
%     end
    for ev = 1:size(EEG.event,2)
        % occlusion
        log_vecs.(cond_labs{2})(ev) = log_vecs.(cond_labs{1})(ev) && strcmp(EEG.event(ev).OCCL, 'ON');
        log_vecs.(cond_labs{3})(ev) = log_vecs.(cond_labs{1})(ev) && strcmp(EEG.event(ev).OCCL, 'OFF');
        % constant
        log_vecs.(cond_labs{4})(ev) = log_vecs.(cond_labs{1})(ev) && strcmp(EEG.event(ev).TRAJ, 'RANDOM1');
        log_vecs.(cond_labs{5})(ev) = log_vecs.(cond_labs{1})(ev) && strcmp(EEG.event(ev).TRAJ, 'CONST');
        log_vecs.(cond_labs{6})(ev) = log_vecs.(cond_labs{1})(ev) && strcmp(EEG.event(ev).TRAJ, 'RANDOM2');
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
%         cond_ind(s).(cond_labs{lab}) = unique([ ...
%             EEG.event(log_vecs.(cond_labs{lab})).epoch]);
        % calculate mean of EEG data across condition
        mean_struct.(cond_labs{lab})(:,:,s) = mean(...
            EEG.data(:, :, cond_ind(s).(cond_labs{lab})), 3);
        % get z-normalized data
        z_mean_struct.(cond_labs{lab})(:, :, s) = normalize( ...
            mean_struct.(cond_labs{lab})(:, :, s),2, 'zscore');
        % get weights
        z_mean_struct.(strjoin(["num_epochs_", cond_labs{lab}],""))(s) = size( ...
            EEG.data(:,:,cond_ind(s).(cond_labs{lab})), 3);
    end

    % calculate contrasts by simple subtraction
    mean_struct.diff_vis_occ(:,:,s) = mean_struct.vis(:,:,s) - mean_struct.occ(:,:,s);
    mean_struct.diff_const_rand1(:,:,s) = mean_struct.const(:,:,s) - mean_struct.rand1(:,:,s);
    mean_struct.diff_const_rand2(:,:,s) = mean_struct.const(:,:,s) - mean_struct.rand2(:,:,s);
    mean_struct.diff_rand1_rand2(:,:,s) = mean_struct.rand1(:,:,s) - mean_struct.rand2(:,:,s);

    % calculate contrasts by simple subtraction
    z_mean_struct.diff_vis_occ(:,:,s) = z_mean_struct.vis(:,:,s) - z_mean_struct.occ(:,:,s);
    z_mean_struct.diff_const_rand1(:,:,s) = z_mean_struct.const(:,:,s) - z_mean_struct.rand1(:,:,s);
    z_mean_struct.diff_const_rand2(:,:,s) = z_mean_struct.const(:,:,s) - z_mean_struct.rand2(:,:,s);
    z_mean_struct.diff_rand1_rand2(:,:,s) = z_mean_struct.rand1(:,:,s) - z_mean_struct.rand2(:,:,s);

    % get time vec
    mean_struct.time_vec = EEG.times;
    z_mean_struct.time_vec = EEG.times;

    % add peaks to know how much power each condition has
    for lab = 1:length(cond_labs)
        count_peaks.(cond_labs{lab}) = length(cond_ind(s).(cond_labs{lab})) + ...
            count_peaks.(cond_labs{lab});
    end
    % add vague estimates of the number of peacs for the peak counts.
    count_peaks.diff_occ = (length(cond_ind(s).occ) + length(cond_ind(s).vis))/2 +...
        count_peaks.diff_occ;
    count_peaks.diff_const_rand1 = (length(cond_ind(s).const) + length(cond_ind(s).rand1))/2 +...
        count_peaks.diff_const_rand1;
    count_peaks.diff_const_rand2 = (length(cond_ind(s).const) + length(cond_ind(s).rand2))/2 +...
        count_peaks.diff_const_rand2;
    count_peaks.diff_rand1_rand2 = (length(cond_ind(s).rand1) + length(cond_ind(s).rand2))/2 +...
        count_peaks.diff_rand1_rand2;

end

% create the weights
for lab = 1:length(cond_labs)
    % divide the respective conditions for each subject by the total number
    % of trials that subject had. The resulting porportion is your weight. 
    epoch_weights.(cond_labs{lab}) = median(z_mean_struct.(...
        strjoin(["num_epochs_", cond_labs{lab}],""))) ./ z_mean_struct.(...
        strjoin(["num_epochs_", cond_labs{lab}],""));
end

% weigh the data according to trials per subject
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
save("mean_struct.mat", 'mean_struct')
save("condition_indices.mat", 'cond_ind')
save("count_peaks.mat", 'count_peaks')
save("z_mean_struct.mat", 'z_mean_struct')


