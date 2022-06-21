% Read in Tracking Datafs


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


%% playground
%
% vec_x = [1, 1.5, 1.75]%, 2.1, 2.4, 3, 3.1]; % 7 points
% diff(vec_x)
% vec_y = [0, 7, 8]% 5.5, 2, -4, 0];
% diff(vec_y)
% % upsample to 14 points
% interp = spline(vec_x, vec_y, linspace(1, vec_x(end), length(vec_x)*2));


%% Upsample Tracking Data to Match EEG Data SR

% decided on spline interpolation because of this convincing comparison:
% https://de.mathworks.com/help/matlab/ref/pchip.html
% spline works best to display oscillatory movements.
% keep only subject ids that are also a field in track_data.
% this is quite experimental code cause it is late and I wanna get done.
% prone to cause errors.
FR = 60;
SR = 250;

tmp_idx = zeros(1,length(subj_ids));
for s = 1:length(subj_ids)
    try
        tmp_idx(s) = cell2mat(strfind({track_data.subject}, subj_ids(s)));
    end
end
subj_ids = subj_ids(logical(tmp_idx));

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

end


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

field_names = fieldnames(track_data);
% nope, it does not seem to be possible to concatenate structures
% without a loop
for s = 1:size(track_data,2)
    if strcmp(track_data(s).subject, all_data_struct(s).subject)
        for f = 1:length(field_names)
            all_data_struct(s).(field_names{f}) = track_data(s).(field_names{f});
        end
    else
        disp(strcat(all_data_struct(s).subject, " has no matching tracking data"));
    end
end


%% Get Latency Vector that Shows Beginnings of Trials
% (Necessary cause Peaks are only given in Trial Latency)
% # List of Trigger Meanings
% # sort(unique(trig_list[[subj_ids[i]]][["trig_data"]][["Description"]]))
% # "S 11": Start of Experiment
% # "S 12": Fixation cross appears
% # "S 13": Start of Trial (left)
% # "S 14": Start of Trial (right)
% # "S 16": Pause start
% # "S 17": Pause end
% # "S 20": Cursor Occlusion
% # "S 21": Cursor Reappearance
% # "S 23": Start Constant Trajectory
% # "S 24": End Constant Trajectory
% # "S 26": Start Start Vector
% # "S 27": End Start Vector
% custom: "S 15": End Trial
% custom: "S 30": Peak


for s = 1:size(all_data_struct,2)

    % copy urevent field to other field
    all_data_struct(s).track_event = all_data_struct(s).urevent;
    % copy event field to convenient format
    event = all_data_struct(s).track_event;
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


    for idx = 1:length(trial_starts)
        % set latencies
        subtracted =  num2cell([event(trial_starts(idx):trial_ends(idx)).trial_latency] - ...
            event(trial_starts(idx)).trial_latency);
        [event(trial_starts(idx):trial_ends(idx)).trial_latency] = subtracted{:};
    end

    all_data_struct(s).track_event = event;

end

% Notes: 91L3HA is excluded because recording started too late.

% add latency of trials

% # end of start vector is beginning of trial
% StartPattern = list(c("S 26","S 27"), c("S 26", "", "S 27")) # possible sequence
% # of triggers that indicates that the trial starts
%
% # fixation cross of the next trial or pause is the end of the trial
% EndPattern = list(c("S 16","S 17"), c("S 16", "", "S 17"), c("S 12", "S 13"),
%                   c("S 12", "", "S 13"), c("S 12", "S 14"), c("S 12", "", "S 14"))
% # possible sequence of triggers that indicates that the trial ends
%
% # there is a weird anomaly in the data: subject 91L3H has a stray "S 16" without
% # a matching "S 17" following at position 187. Because it messes with the
% # algorithm, I remove it manually
% trig_list[["91L3H"]][["trig_data"]] = trig_list[["91L3H"]][["trig_data"]][-187,]

% ISSUE: currently there is an error: when s = 16, there are only 54 end trial
% triggers but 72 start trial triggers. why?

%% Find Peaks of Tracking Data and put them in Events

% "S 30": Peak

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
            % replace the copied rows' type with S 30, making sure the
            % second event is overwritten
            event(current_event_idx+1).type = 'S 30';
            event(current_event_idx+1).code = 'inserted';
            event(current_event_idx+1).latency = current_peak_latencies(idx);
            event(current_event_idx+1).trial_latency = current_trial_peak_latencies(idx);

        end
    end
    % copy the event structure into a new field in all_data_struct
    all_data_struct(s).track_event_peaks =  event;
end
%plot(track_data.task_a.(field_names{s}).(trial_name{t})(:,1), track_data.task_a.(field_names{s}).(trial_name{t})(:,2))

%get pixel coordinates of extracted maxima to calculate euclidean distance
% max_traj_coords = [trialdata.traj_x_pix(index_max_traj), trialdata.traj_y_pix(index_max_traj)];

%create vector of 1 to append to matrix (marking maxima)
%id = zeros(size(index_max_traj,1),1)+1;

%create matrix for analyses
%max_traj = [max_traj_coords, index_max_traj, id];

%get pixel coordinates of extracted minima to calculate euclidean distance
% min_traj_coords = [trialdata.traj_x_pix(index_min_traj), trialdata.traj_y_pix(index_min_traj)];




%% Plot Tracking Data

s = 2; % pick subject
trials = 1:10;

figure()
for t = trials
    trial_starts = find(strcmp({all_data_struct(s).track_event_peaks.type}, "S 27"));
    trial_ends = find(strcmp({all_data_struct(s).track_event_peaks.type}, "S 15"));
    trial_x_only = all_data_struct(s).track_event_peaks(trial_starts(t):trial_ends(t));
    peaks_only = find(strcmp({trial_x_only.type}, "S 30"));
    peak_latencies = [trial_x_only(peaks_only).trial_latency];
    const_only = find(strcmp({trial_x_only.type}, "S 23"));
    const_latencies = [trial_x_only(const_only).trial_latency];
    const_only_end = find(strcmp({trial_x_only.type}, "S 24"));
    const_latencies_end = [trial_x_only(const_only_end).trial_latency];
    end_trial_only = find(strcmp({trial_x_only.type}, "S 15"));
    trial_latencies_end = [trial_x_only(end_trial_only).trial_latency];
    start_trial_only = find(strcmp({trial_x_only.type}, "S 27"));
    trial_latencies_start = [trial_x_only(start_trial_only).trial_latency];
    traj_y = all_data_struct(s).upsamp_trials(t).traj_y;
    traj_x_temp = 1:length(traj_y);

    % sanity check that trigger plotting works


    subplot(5, 2, t)
    plot(traj_x_temp, traj_y, "linewidth", 1.5)
    hold on
    peak_handle = xline(peak_latencies,"linewidth", 1.5);
    const_handle = xline(const_latencies, "r", "linewidth", 1.5);
    xline(const_latencies_end, "r", "linewidth", 1.5)
    start_handle = xline(trial_latencies_start, "g", "linewidth", 1.5);
    hold off
    legend([peak_handle(1), const_handle, start_handle], {"peaks", "const", "trial"});

end



%% Epoch around Peaks




