% Add additional Fields to Event structure

% Adds error per epoch, latency between pursuit peak and trajectory peak

% Author: Saskia Wilken
% Creation Date: 05.09.2022

%% Empty

format compact
format long G
clear
clc


%% Define Paths

% TODO: Subject 03 something s =1 does neither have a pursuit latency nor
% an artifact_error field. WHY?!

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
for idx = 1:length(files2read)

    EEG = pop_loadset('filename',files2read{idx});
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG);

end
eeglab redraw


%% Calculate Error per Epoch and put in event field

for s = 1:length(ALLEEG)
    % get plus / minus latencies around peaks
    epoch_lims = [ALLEEG(s).xmin, ALLEEG(s).xmax];
    % transform ms time points in samples
    epoch_lims = epoch_lims * ALLEEG(s).srate;
    % get the latencies of the peaks around which epoching was done
    %     count = 0;
    %     for peak = 1:length(ALLEEG(s).urevent)
    %         count = strcmp(ALLEEG(s).urevent(peak).type, 'S 40') + count;
    %     end % ok, so peak latencies contains all the peak event latencies.

    EEG = ALLEEG(s);

    for ep = 1:size(ALLEEG(s).data, 3)

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
        if epoch_trial ~= 9999
            % and if the epoch is not within the first half a second of a
            % trial
            if strcmp(EEG.event(peak_idx).artifact_error, 'VALID')
                % get latency inside trial
                epoch_trial_latency = EEG.event(peak_idx).trial_latency;

                % calculate error in epoch
                %  get current error
                cur_error = abs(track_data(s).upsamp_data.(epoch_task)(epoch_trial).error);
                epoch_start = epoch_trial_latency+epoch_lims(1);
                epoch_end = epoch_trial_latency+epoch_lims(2);
                % account for first and last epochs
                if epoch_start < 1
                    epoch_start = 1;
                end
                if epoch_end > length(cur_error)
                    epoch_end = length(cur_error);
                end
                error_of_epoch = cur_error(epoch_start:epoch_end);
                % work on event field
                % get event fields that belong to epoch
                epoch_idx = find([EEG.event.epoch] == ep);
                % put error in epoch in event field
                for srow = 1:length(epoch_idx)
                    EEG.event(epoch_idx(srow)).('epoch_error') = mean(error_of_epoch, 'omitnan');
                end
            end
        end

    end
    ALLEEG(s) = EEG;
end

eeglab redraw

% eeg_retrieve() % Retrieve an EEG dataset from the variable
%                    containing all datasets (standard: ALLEEG).


%% Calculate latency of pursuit peak following trajectory peak and put in event field

% get event peak
% check if next event is pursuit peak in same direction
% if yes, get latency of pursuit peak and write in epoch-linked event field
% if no, count one plus and skip to next event
% count how many trajectory peaks had no following pursuit peak.
no_purs = 0;
purs_count = 0;

for s = 1:length(ALLEEG)
    % get plus / minus latencies around peaks
    epoch_lims = [ALLEEG(s).xmin, ALLEEG(s).xmax];
    % transform ms time points in samples
    epoch_lims = epoch_lims * ALLEEG(s).srate;
    % get the latencies of the peaks around which epoching was done
    %     count = 0;
    %     for peak = 1:length(ALLEEG(s).urevent)
    %         count = strcmp(ALLEEG(s).urevent(peak).type, 'S 40') + count;
    %     end % ok, so peak latencies contains all the peak event latencies.

    EEG = ALLEEG(s);

    for ep = 1:size(ALLEEG(s).data, 3)

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
                if ~isempty(purs_lat) && purs_lat > 80 && purs_lat < 500
                    % put latency of pursuit in epoch in event field
                    for srow = 1:length(epoch_idx)
                        EEG.event(epoch_idx(srow)).('pursuit_lat') = purs_lat;
                    end
                else
                    for srow = 1:length(epoch_idx)
                        EEG.event(epoch_idx(srow)).('pursuit_lat') = [];
                    end
                    no_purs = no_purs + 1;
                    purs_count = purs_count - 1;
                end
            end
        end

    end

    pop_saveset(EEG, 'filename', EEG.setname, 'filepath', ...
        char(epochs_plus_error_dir));
    ALLEEG(s) = EEG;
end

eeglab redraw


%% Check what the pursuit latencies look like

% criteria:
% - pursuit peak with prominence = 0.05
% - a pursuit peak must follow the trajectory peak in the same epoch
% - pursuit and trajectory peaks must go in the same direction (up or down)
% - pursuit peak latency is between 80 and 500 ms. 
% - central peak of the epoch is not in the first half a second of a trial

no_purs/((purs_count+no_purs)/100)
% at this time, 31472 out of 65627 trajectory peaks are not followed by a 
% pursuit peak (47.96%)
pursuit_lats = [];
all_pursuit_lats = [];
for s = 1:length(ALLEEG)
    pursuit_lats = [ALLEEG(s).event.pursuit_lat];
    all_pursuit_lats = [all_pursuit_lats, pursuit_lats];
end
% keep only one of the copies of the pursuit latency field values each. 
unique_lats_idx = find(diff(all_pursuit_lats));
all_pursuit_lats_unique = all_pursuit_lats(unique_lats_idx);

plot(sort(all_pursuit_lats_unique))
histogram(all_pursuit_lats_unique)
mink(all_pursuit_lats_unique, 100)
% purs_count and length(all_pursuit_lats_unique) do not have the same size.
% which is odd. 


