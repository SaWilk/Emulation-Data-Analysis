% Merge eeg and tracking data

% creation date: 02.06.2022 
% Author: Saskia Wilken

% this script reads in tracking data, upsamples it, calculates tracking
% error (TBD), loads in preprocessed continuous (TBD) EEG files, imputes 
% triggers in the event structure field and creates a latency vector that allows
% for synchronizing the tracking and the eeg data. It saves eeglab .sets
% that contain the tracking data (TBD)

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


%% Load Tracking data

warning("off","all");
% for now only for a
% TODO: Add in functionality for b
a_paths = subj_paths(endsWith(subj_paths, "_A"));
% what is contained in the different columns
col_names = {'traj_x', 'traj_y', 'purs_x', 'purs_y', 'error', 'error_abs'};

for p = 1:length(a_paths)
    if ~isempty(a_paths{p}) % check if there is data
        cd(a_paths{p});
        file_names_curr_subj = {dir('*traj_purs*.csv').name};
        track_data(p).("subject") = subj_ids{p};
        track_data(p).("group") = 'control'; % TODO: Adjust if there is ever another group
        track_data(p).("condition") = {'task_A'}; % TODO: adjust for other tasks
        track_data(p).("path") = a_paths{p};
        for t = 1:length(file_names_curr_subj)
            tmp_mat = table2array(readtable(file_names_curr_subj{t}));
            for c = 1:size(col_names,2)
                track_data(p).("trials")(t).(col_names{c}) = tmp_mat(:,c);
            end
            clear tmp_mat
        end
    else fprintf('%s contains no data', a_paths{p});
    end
end
warning("on","all");


%% Upsample Tracking Data to Match EEG Data SR

% decided on spline interpolation because of this convincing comparison:
% https://de.mathworks.com/help/matlab/ref/pchip.html
% spline works best to display oscillatory movements.

FR = 60;
SR = 250;

% keep only subject ids that are also a field in track_data.

for s = 1:length(track_data)

        % this is such smart code...
        % it will loop through subj ids but give an index where the
        % elements of subj ids are in track_data
        tmp_idx(s) = find(cellfun(@isequal, subj_ids , ...
            repmat({track_data(s).subject}, size(subj_ids))));

end
subj_ids = subj_ids(tmp_idx);
clear tmp_idx

% upsample tracking data and transfer into time domain
for s = 1:length(subj_ids)
    for t = 1:size(track_data(s).trials, 2)
        tmp_traj_y = track_data(s).trials(t).traj_y;
        time_vec = linspace(0, length(tmp_traj_y)/FR*1000, round(length(tmp_traj_y))).';
        inter_time_vec = linspace(time_vec(1), time_vec(end), round(length(time_vec)/FR*SR));
        inter_traj_y = spline(time_vec, tmp_traj_y, inter_time_vec);

        % proof that it works well:
        %         figure
        %         plot(no_points, inter_traj_y, 'ob')
        %         hold on
        %         plot(time_vec, tmp_traj_y, '.r')
        %         hold off
        % same for pursuit
        tmp_purs_y = track_data(s).trials(t).purs_y;
        inter_purs_y = spline(time_vec, tmp_purs_y, inter_time_vec);
        % note that x axis points are milliseconds now.

        % proof that it works well:
        %         figure
        %         plot(tmp_purs_x, tmp_purs_y, 'ob')
        %         hold on
        %         plot(no_points, inter_purs_y, '.r')
        %         hold off

        % for now, I am just deleting the error columns
        % TODO: adjust or re-calculate error in Matlab
        tmp_mat = [inter_time_vec; inter_traj_y; inter_purs_y].';
        for c = 1:size(tmp_mat, 2)
            track_data(s).("upsamp_trials")(t).(col_names{c}) = tmp_mat(:,c);
        end
    end
    clear tmp_traj_y tmp_purs_y tmp_mat

end


%% Calculate difference between pursuit and trajectory


%% Load EEG data

% now we need eeglab
eeglab;
cd(strcat(eeg_data_path, '\preprocessed\'));
% get all set files in the folder
files2read = {dir('*.set').name};
% only get files for a for now
% TODO: adjust for b
files2read_a = files2read(~cellfun(@isempty, regexp(files2read,'_A_')));

eeg_struct = pop_loadset('filename',files2read);


%% Put Tracking and EEG Data in the same structure

all_data_struct = eeg_struct;

% remove all subjects that have no correcponding tracking data

for s = 1:length(all_data_struct)

        % this is such smart code...
        % it will loop through all data struct but give an index where the 
        % elements of all_data_struct are in track_data
        tmp_idx(s) = find(cellfun(@isequal, {track_data.subject}, ...
            repmat({all_data_struct(s).subject}, size({track_data.subject}))));

end
copy_data = track_data(tmp_idx);

field_names = fieldnames(copy_data);
% nope, it does not seem to be possible to concatenate structures
% without a loop
for s = 1:size(copy_data,2)
    if strcmp(copy_data(s).subject, all_data_struct(s).subject)
        for f = 1:length(field_names)
            all_data_struct(s).(field_names{f}) = copy_data(s).(field_names{f});
        end
    else
        disp(strcat(all_data_struct(s).subject, " has no matching tracking data"));
    end
end


%% Get Latency Vector that Shows Beginnings of Trials
% (Necessary cause Peaks are only given in Trial Latency)
%   "instruction":  10,
%     "exp_start":    11,
%     "fix":          12,
%     "trial_start_L":13,
%     "trial_start_R":14,
%     "trial_end":    15, # superfluous with fix trigger (except for break)
%     "pause_start":  16,
%     "pause_end":    17,
%     "exp_end":      18,
%     "button":       19,
%     "occlusion":    20,
%     "reappear":     21,
%     "fourth_trial": 22,
%     "start_constant": 23,
%     "end_constant": 24,
%     "fourth trial": 25,
%     "start_startvec": 26,
%     "end_startvec": 27,
%     "C_too_early":  28,
%     "C_just_right": 29,
%     "C_too_late":   30,
% custom: "S 15": End Trial
% custom: "S 40": Peak


for s = 1:size(all_data_struct,2)

    % copy urevent field to other field
    all_data_struct(s).track_event = all_data_struct(s).urevent;
    % copy event field to convenient format
    event = all_data_struct(s).track_event;
    % sanity check of triggers - are we DEALING WITH TASK A EVEN???
    % TODO: Alter check once pipeline is adjusted for Task B as well.
    event_cat = categorical({event.type});
    categories(event_cat);
    % if you should ever feel the need to count cats, there you go: countcats(event_cat)

    if any(strcmp(event_cat, "S 30")) | any(strcmp(event_cat, "S 28"))
        warning(strcat(['Caution! Subject', all_data_struct(s).subject, 'seems to be a task C trigger file. Skipping']));
        continue
    elseif any(strcmp(event_cat, "S 20")) | any(strcmp(event_cat, "S 21"))
        warning(strcat(['Caution! Subject', all_data_struct(s).subject, ' seems to be a task B trigger file. Skipping']));
        continue
    end

    % delete all triggers before S11
    event = event(find(strcmp({event.type}, "S 11")):end);
    % copy event types column to convenient format

    % insert end trial triggers

    first_12 = true;
    last_trigger = "none";
    % loop through elements of event_types.
    i = 0;
    while i < length({event.type})
        i = i + 1;
        % if an event type is either 16 or 12...
        if (strcmp(event(i).type, "S 16") | strcmp(event(i).type, "S 12"))
            % make sure the first entry is skipped (since the first fix cross
            % is not a trial end)

            if first_12
                first_12 = false;
                continue
            end

            % and exclude the first 12 after a pause, cause it's also not a
            % trial end.
            if strcmp(last_trigger, "S 16")
                last_trigger = "none";
                continue
            end
            %             disp(strcat("is true for ", num2str(i)))
            %             disp(strcat("following event", event(i+1).type))
            % save current trigger for next loop iteration
            last_trigger = event(i).type;
            % copy the row in question
            event = [event(1:i-1), event(i), event(i:end)];
            % replace the copied rows' type with S 15
            event(i).type = 'S 15';
            event(i).code = 'inserted';
            % increment i additionally because we just inserted a trigger
            i = i + 1;
        end

    end
    % add an end trial trigger to very end of trigger list
    event = [event, event(end)];
    event(end).type = 'S 15';
    event(end).code = 'inserted';
    % copy the new event_types struct into track event
    all_data_struct(s).track_event = event;

    % start trial indices
    % get the indices of the "S 27" in event.type
    trial_start_ind = strcmp({event.type}, "S 27");
    trial_starts = find(trial_start_ind);

    % end trial indices
    % get the indices of the "S 15" in event.type
    trial_end_ind = strcmp({event.type}, "S 15");
    trial_ends = find(trial_end_ind);

    latencies = [event.latency];

    % add new field: trial_latency
    tmp = num2cell(latencies);
    % todo here: latencies set to 0 at beginnings of trials so I can use
    % the findpeaks function to get the right timing
    % probably requires a loop.

    % this syntax is confusing as hell. left hand side of the assignment needs
    % to be a vector, right hand side single outputs and matlab is happy...
    [all_data_struct(s).track_event.trial_latency] = tmp{:};

    % update event
    event = all_data_struct(s).track_event;

    if length(trial_starts) == length(trial_ends)
        for idx = 1:length(trial_starts)
            % set latencies
            subtracted =  num2cell([event(trial_starts(idx):trial_ends(idx)).trial_latency] - ...
                event(trial_starts(idx)).trial_latency);
            [event(trial_starts(idx):trial_ends(idx)).trial_latency] = subtracted{:};
        end
    else
        warning(strcat(['Subject ', all_data_struct(s).subject, ...
            ' has an unequal number of trial start and trial end triggers!', ...
            ' No trial latencies will be computed. ']));
    end

    all_data_struct(s).track_event = event;

end


% Notes: 91L3HA is excluded because recording started too late.
% WM87B is excluded cause there is no behavioral data. 


%% Find Peaks of Tracking Data and put them in Events

% "S 40": Peak

for s = 1:size(all_data_struct,2)

    event = all_data_struct(s).track_event;

    for t = 1:size(all_data_struct(s).upsamp_trials,2)

        % get peak indices within trial trajectory data
        [~, index_max_traj] = findpeaks(all_data_struct(s).upsamp_trials(t).traj_y);
        [~, index_min_traj] = findpeaks(-all_data_struct(s).upsamp_trials(t).traj_y);
        % locate trial starts in eeg event
        trial_starts = find(strcmp({event.type}, "S 27"));
        % locate trial ends in eeg event
        trial_ends = find(strcmp({event.type}, "S 15"));

        if length(trial_starts) ~= length(trial_ends)
            disp(strcat(['Subject ', all_data_struct(s).subject, ...
                ' does not have an equal number of start and end trial triggers.']))
            continue
        end

        % get the latencies of the current trial in the eeg data
        current_start_latency = event(trial_starts(t)).latency;
        current_end_latency = event(trial_ends(t)).latency;
        % add current start latency to the indicies of the peaks (which are in the same sampling rate
        % as the eeg signal due to upsampling)
        current_trial_peak_latencies = sort([index_max_traj;index_min_traj]);
        current_peak_latencies = current_start_latency+current_trial_peak_latencies;

        % for each of the current peak latencies...
        for idx = 1:length(current_peak_latencies)

            % find the event that is one event before the peak
            % (which is the largest negative latency)
            tmp = [[event.latency]- current_peak_latencies(idx)];
            current_event_idx = max(find(tmp <= 0 ));

            % insert the peak event markers at position in question
            % (copy event before it and adjust its values)
            event = [event(1:current_event_idx-1), event(current_event_idx), event(current_event_idx:end)];
            % replace the copied rows' type with S 40, making sure the
            % second event is overwritten
            event(current_event_idx+1).type = 'S 40';
            event(current_event_idx+1).code = 'inserted';
            event(current_event_idx+1).latency = current_peak_latencies(idx);
            event(current_event_idx+1).trial_latency = current_trial_peak_latencies(idx);

        end
    end
    % copy the event structure into a new field in all_data_struct
    all_data_struct(s).track_event_peaks =  event;
end
%plot(copy_data.task_a.(field_names{s}).(trial_name{t})(:,1), copy_data.task_a.(field_names{s}).(trial_name{t})(:,2))

%get pixel coordinates of extracted maxima to calculate euclidean distance
% max_traj_coords = [trialdata.traj_x_pix(index_max_traj), trialdata.traj_y_pix(index_max_traj)];

%create vector of 1 to append to matrix (marking maxima)
%id = zeros(size(index_max_traj,1),1)+1;

%create matrix for analyses
%max_traj = [max_traj_coords, index_max_traj, id];

%get pixel coordinates of extracted minima to calculate euclidean distance
% min_traj_coords = [trialdata.traj_x_pix(index_min_traj), trialdata.traj_y_pix(index_min_traj)];


%% Replace Constant Traj Trigger from .vmrk with calculated from .npz
%align constant trajectories of several trials
%compare constant trajectory and trial trajectory at every point
%calculate variance of difference
%should be minimal at constant trajectory

% [shift, start_ind, end_ind] = align_const_traj(const_traj, trial_traj);


%% save sets

    TMPEEG = pop_saveset(TMPEEG,'filename',[files2read{ind}(1:end-7) '_A_epoched'], 'filepath', char(savepath));
