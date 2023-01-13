% The neurophysiology of continuous action monitoring
% Saskia Wilken, Adriana Böttcher, Nico Adelhöfer, Markus Raab, Sven
% Hoffmann & Christian Beste

% Add additional Fields to Event structure
% Adds event fuelds:
% - error per epoch
% - latency between pursuit peak and trajectory peak

% Created by: 
% Saskia Wilken, General Psychology: Judgement, Decision Making & Action, 
% University of Hagen
% 05.09.2022


%% Empty

format compact
format long G
clear
clc


%% Define Paths

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
input_dir = strjoin([parent_dir_2, "Emulation-Data-Output"], filesep);
subdir_peaks = strjoin([input_dir '06_artifact_rejection' 'peaks_-500_750'], filesep);
track_dir = strjoin([input_dir '03_parallelize_with_traj'], filesep);
% output dir
epochs_plus_error_dir = strjoin([input_dir, "07_epochs_with_extra_fields"], filesep);
mkdir(epochs_plus_error_dir);


%% Load Data

% tracking data
load(strjoin([track_dir, 'all_tracking_data.mat'], '\'));
% load data into eeglab
cd(subdir_peaks);
%list all *.set files in inputpath
file_names = dir('*.set');
%concatenate into one cell array
files2read = {file_names.name};
% load eeg data
for idx = 1:length(files2read)
    EEG = pop_loadset('filename',files2read{idx});
    [ALLEEG, EEG] = eeg_store(ALLEEG, EEG);
end
eeglab redraw
TMPEEG = ALLEEG;


%% Calculate Error per Epoch and put in event field

for s = 1:length(TMPEEG)

    EEG = TMPEEG(s);
    all_epoch_errors = cell(1, size(EEG.data, 3));
    % get plus / minus latencies around peaks
    epoch_lims = [EEG.xmin, EEG.xmax];
    % transform ms time points in samples
    epoch_lims = epoch_lims * EEG.srate;
    % get the latencies of the peaks around which epoching was done

    for ep = 1:size(TMPEEG(s).data, 3)

        % get central peak of epoch
        peak_idx = get_epoch_center(EEG, ep);
        % get task we are in
        epoch_task = EEG.event(peak_idx).task;
        % get number of trial we are in
        epoch_trial = EEG.event(peak_idx).trial_number;
        % do all this only if we are dealing with a valid trial
        if epoch_trial ~= 9999
            % and if the epoch is not within the first half a second of a
            % trial
            if strcmp(EEG.event(peak_idx).artifact_error, 'VALID')
                % get latency of peak inside trial
                epoch_trial_latency = EEG.event(peak_idx).trial_latency;

                % calculate error in epoch
                %  get current error
                cur_error = abs(track_data(s).upsamp_data.(epoch_task)(epoch_trial).error);
                cur_purs = track_data(s).upsamp_data.(epoch_task)(epoch_trial).purs_y;
                cur_track = track_data(s).upsamp_data.(epoch_task)(epoch_trial).traj_y;
                epoch_start = epoch_trial_latency+epoch_lims(1);
                epoch_end = epoch_trial_latency+epoch_lims(2);
                % account for first and last epochs
                if epoch_start < 1
                    epoch_start = 1;
                end
                if epoch_end > length(cur_error)
                    epoch_end = length(cur_error);
                end
                error_of_epoch = mean(cur_error(epoch_start:epoch_end), 'omitnan');

                % work on event field
                % get event fields that belong to epoch
                epoch_idx = find([EEG.event.epoch] == ep);
                % store epoch errors for normalization
                all_epoch_errors{ep} = error_of_epoch;
                % put error in epoch in event field
                for srow = 1:length(epoch_idx)
                    EEG.event(epoch_idx(srow)).('epoch_error') = error_of_epoch;
                end
            end
        end
    end
    % calculate the normalized epoch errors
    % get the indices of the valid pursuit latencies
    all_epoch_errors_idx = find(~cellfun(@isempty, all_epoch_errors));
    % initialize normalization cell
    all_epoch_errors_z = cell(size(all_epoch_errors));
    % normalize and put in cell at according locations
    all_epoch_errors_z(all_epoch_errors_idx) = num2cell(normalize([all_epoch_errors{:}], 'zscore'));
    % loop through all epochs
    for ep_err = 1:length(all_epoch_errors_idx)
        % initialize vector to hold the indices of the current epoch event
        % fields
        current_ep_idx = [];
        % get these indices
        current_ep_idx = find([EEG.event.epoch] == all_epoch_errors_idx(ep_err));
        % distribute the pursuit latencies to the corresponding epochs
        [EEG.event(current_ep_idx).('epoch_error_z')] = deal( ...
            all_epoch_errors_z{all_epoch_errors_idx(ep_err)});
    end
    TMPEEG(s) = EEG;
end

eeglab redraw


%% Calculate latency of pursuit peak following trajectory peak and put in event field

% get event peak
% check if next event is pursuit peak in same direction
% if yes, get latency of pursuit peak and write in epoch-linked event field
% if no, count one plus and skip to next event
% count how many trajectory peaks had no following pursuit peak.
ACCEPTABLE_LATS = [80, 300]; % 80 ms is the reflex threshold, 300 is the 
% mean distance between peaks. 

no_purs = 0;
purs_count = 0;

for s = 1:length(TMPEEG)
    all_purs_lat = cell(1, size(EEG.data, 3));
    % get plus / minus latencies around peaks
    epoch_lims = [TMPEEG(s).xmin, TMPEEG(s).xmax];
    % transform ms time points in samples
    epoch_lims = epoch_lims * TMPEEG(s).srate;
    % get the latencies of the peaks around which epoching was done

    EEG = TMPEEG(s);

    for ep = 1:size(TMPEEG(s).data, 3)

        % get central peak of epoch
        peak_idx = EEG.epoch(ep).event([EEG.epoch(ep).eventlatency{:}] == 0);
        % in the rare case that there are multiple events with a latency of
        % 0, determine which one is the peak
        if length(peak_idx) > 1
            which_peak_idx = find(strcmp({EEG.event(peak_idx).type}, {'S 40'}) );
            which_peak_idx_2 = find(strcmp({EEG.event(peak_idx).type}, {'S 50'}) );
            if ~isempty(which_peak_idx)
                peak_idx = peak_idx(which_peak_idx);
            elseif ~isempty(which_peak_idx_2)
                peak_idx = peak_idx(which_peak_idx_2);
            end
        end

        % get task we are in
        epoch_task = EEG.event(peak_idx).task;
        % get number of trial we are in
        epoch_trial = EEG.event(peak_idx).trial_number;
        % do all this only if we are dealing with a valid trial
        idx_in_trial = find(EEG.epoch(ep).event == peak_idx);

        if epoch_trial ~= 9999
            % and if the epoch is not within the first half a second of a
            % trial
            if strcmp(EEG.event(peak_idx).artifact_error, 'VALID')
                % get latency inside trial
                epoch_trial_latency = EEG.event(peak_idx).trial_latency;
                peak_type = EEG.epoch(ep).eventtype(idx_in_trial);
                % calculate latency to next pursuit peak
                % is there a pursuit peak in the same direction in the epoch?
                all_events = EEG.epoch(ep).eventtype;
                % positive peaks...
                if strcmp(peak_type, 'S 40')
                    purs_peak_idx = find(strcmp(all_events, 'S 41'));
                    % if there is a pursuit peak in the same direction and it
                    % is after the trajectory peak...
                    if ~isempty(purs_peak_idx) && any(purs_peak_idx > idx_in_trial)
                        % if there are multiple pursuit peaks, keep only the first one
                        % after the peak
                        if length(purs_peak_idx) > 1
                            purs_peak_idx = purs_peak_idx(min(find(purs_peak_idx > idx_in_trial)));
                        end
                        % purs latency in ms
                        purs_lat = EEG.epoch(ep).eventlatency{purs_peak_idx};
                        purs_count = purs_count + 1;
                    else
                        no_purs = no_purs + 1;
                        purs_lat = [];
                    end
                    % negative peaks...
                elseif strcmp(peak_type, 'S 50')
                    purs_peak_idx = find(strcmp(all_events, 'S 51'));
                    % if there is a pursuit peak in the same direction and it
                    % is after the trajectory peak...
                    if ~isempty(purs_peak_idx) && any(purs_peak_idx > idx_in_trial)
                        % if there are multiple pursuit peaks, keep only the first one
                        % after the peak
                        if length(purs_peak_idx) > 1
                            purs_peak_idx = purs_peak_idx(min(find(purs_peak_idx > idx_in_trial)));
                        end
                        % purs latency in ms
                        purs_lat = EEG.epoch(ep).eventlatency{purs_peak_idx};
                        purs_count = purs_count + 1;
                    else
                        no_purs = no_purs + 1;
                        purs_lat = [];
                    end
                else
                    error('There is no peak at latency 0 of the epoch')
                end

                % work on event field
                % get event fields that belong to epoch
                epoch_idx = find([EEG.event.epoch] == ep);
                % keep only latencies in a sensible range (80 ms to 500 ms)
                if ~isempty(purs_lat) && purs_lat > ACCEPTABLE_LATS(1) && purs_lat < ACCEPTABLE_LATS(2)
                    % put latency of pursuit in epoch in event field
                    % store latency for normalization
                    all_purs_lat{ep} = purs_lat;
                    for srow = 1:length(epoch_idx)
                        EEG.event(epoch_idx(srow)).('pursuit_lat') = purs_lat;
                    end
                else
                    for srow = 1:length(epoch_idx)
                        EEG.event(epoch_idx(srow)).('pursuit_lat') = [];
                        % also add normalized latency
                        EEG.event(epoch_idx(srow)).('pursuit_lat_z') = [];
                    end
                    no_purs = no_purs + 1;
                    purs_count = purs_count - 1;
                end
            end
        end

    end
    % calculate the normalized latencies
    % get the indices of the valid pursuit latencies
    all_purs_lat_idx = find(~cellfun(@isempty, all_purs_lat));
    % initialize normalization cell
    all_purs_lat_z = cell(size(all_purs_lat));
    % normalize and put in cell at according locations
    all_purs_lat_z(all_purs_lat_idx) = num2cell(normalize([all_purs_lat{:}], 'zscore'));
    % loop through all epochs
    for ep_purs = 1:length(all_purs_lat_idx)
        % initialize vector to hold the indices of the current epoch event
        % fields
        current_ep_idx = [];
        % get these indices
        current_ep_idx = find([EEG.event.epoch] == all_purs_lat_idx(ep_purs));
        % distribute the pursuit latencies to the corresponding epochs
        [EEG.event(current_ep_idx).('pursuit_lat_z')] = deal( ...
            all_purs_lat_z{all_purs_lat_idx(ep_purs)});
    end

    pop_saveset(EEG, 'filename', EEG.setname, 'filepath', ...
        char(epochs_plus_error_dir));
    TMPEEG(s) = EEG;
end

eeglab redraw

