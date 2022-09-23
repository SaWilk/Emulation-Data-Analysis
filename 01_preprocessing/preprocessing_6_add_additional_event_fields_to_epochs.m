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
    %     count = 0;
    %     for peak = 1:length(TMPEEG(s).urevent)
    %         count = strcmp(TMPEEG(s).urevent(peak).type, 'S 40') + count;
    %     end % ok, so peak latencies contains all the peak event latencies.

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

% eeg_retrieve() % Retrieve an EEG dataset from the variable
%                    containing all datasets (standard: TMPEEG).


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
    %     count = 0;
    %     for peak = 1:length(TMPEEG(s).urevent)
    %         count = strcmp(TMPEEG(s).urevent(peak).type, 'S 40') + count;
    %     end % ok, so peak latencies contains all the peak event latencies.

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
                        % also add normalized latency
                        % EEG.event(epoch_idx(srow)).('pursuit_lat_z') = normalize(purs_lat, 'zscore');

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


%% Check what the pursuit latencies look like

% criteria:
% - pursuit peak with prominence = 0.05
% - a pursuit peak must follow the trajectory peak in the same epoch
% - pursuit and trajectory peaks must go in the same direction (up or down)
% - pursuit peak latency is between 80 and 300 ms. 
% - central peak of the epoch is not in the first half a second of a trial

disp(strjoin([round(no_purs/((purs_count+no_purs)/100),0), ...
    "% (", no_purs, " out of ", purs_count+no_purs, ...
    ") trajectory peaks are not followed by a pursuit peak."],""))
% at this time, 31472 out of 65627 trajectory peaks are not followed by a 
% pursuit peak (47.96%)
pursuit_lats = [];
all_pursuit_lats = [];
pursuit_lats_z = [];
all_pursuit_lats_z = [];
epoch_errors = [];
all_epoch_errors = [];
epoch_errors_z = [];
all_epoch_errors_z = [];

for s = 1:length(TMPEEG)
    % pursuit latencies raw
    pursuit_lats = [TMPEEG(s).event.pursuit_lat];
    % concatenate pursuit latencies of current subj with previous subj
    all_pursuit_lats = [all_pursuit_lats, pursuit_lats];
    % normalized pursuit latencies
    pursuit_lats_z = [TMPEEG(s).event.pursuit_lat_z];
    all_pursuit_lats_z = [all_pursuit_lats_z, pursuit_lats_z];
    % epoch error raw
    epoch_errors = [TMPEEG(s).event.epoch_error];
    all_epoch_errors = [all_epoch_errors, epoch_errors];
    % normalized epoch error
    epoch_errors_z = [TMPEEG(s).event.epoch_error_z];
    all_epoch_errors_z = [all_epoch_errors, epoch_errors_z];
end

% keep only one of the copies of the pursuit latency field values each. 
unique_lats_idx = find(diff(all_pursuit_lats));
% this method is not 100% failsafe - there is a chance that two following
% pursuit latencies are exactly the same and they would be missed by this
% algorithm. However, most pursuit lats should be in there. 
all_pursuit_lats_unique = all_pursuit_lats(unique_lats_idx);
unique_lats_idx_z = find(diff(all_pursuit_lats_z));
all_pursuit_lats_unique_z = all_pursuit_lats_z(unique_lats_idx_z);

BINWIDTH = 12;
range(all_pursuit_lats_unique)
PROMINENCE_THRESH = 0.05; 
ACCEPTABLE_LATS = [80, 300];
histogram(all_pursuit_lats_unique, ...
    round(range(all_pursuit_lats_unique)/BINWIDTH));
title(['Pursuit Latency'])
subtitle(strjoin(["Prominence thresh ", PROMINENCE_THRESH, ...
    ", peaks in same direction, latency lims between ", ...
    ACCEPTABLE_LATS, " ms", "binwidth = ", BINWIDTH, " ms"]))
mink(all_pursuit_lats_unique, 100)
% purs_count and length(all_pursuit_lats_unique) do not have the same size.
% which is odd. 

PROMINENCE_THRESH = 0.05; 
ACCEPTABLE_LATS = [80, 300];
histogram(all_pursuit_lats_unique_z, ...
    round(range(all_pursuit_lats_unique)/BINWIDTH));
title(['Pursuit Latency z-transformed per subject'])
subtitle(strjoin(["Prominence thresh ", PROMINENCE_THRESH, ...
    ", peaks in same direction, latency lims between ", ...
    ACCEPTABLE_LATS, " ms", "binwidth = ", BINWIDTH, " ms"]))
mink(all_pursuit_lats_unique, 100)

% check the epoch error distribution

figure()
BINWIDTH = 0.02;
histogram(all_epoch_errors, round(range(all_epoch_errors)/BINWIDTH));
title(['Epoch Error'])
subtitle(strjoin(["All epochs, all errors over time, binwidth = ", ...
    BINWIDTH*100, " % of screen"]))

figure()
histogram(all_epoch_errors_z, round(range(all_epoch_errors)/BINWIDTH));
title(['Epoch Error z-transformed'])
subtitle(strjoin(["All epochs, all errors over time, binwidth = ", ...
    BINWIDTH*100, " % of screen"]))

