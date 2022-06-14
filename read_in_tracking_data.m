% Read in Tracking Data


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
SR = 250; % sampling rate of EEG data
FR = 60; % frame rate of tracking data
copy_struct = track_data;

for s = 1:length(subj_ids)
    if strcmp(track_data(s).subject, subj_ids{s}) % check if there is a field with the current subject id
        for t = 1:size(track_data(s).trials, 2)
            tmp_traj_x = track_data(s).trials(t).traj_x;
            tmp_traj_y = track_data(s).trials(t).traj_y;
            inter_traj_x = linspace(tmp_traj_x(1), tmp_traj_x(end), round(length(tmp_traj_x)/FR*SR));
            inter_traj_y = spline(tmp_traj_x, tmp_traj_y, inter_traj_x);

            % % proof that it works well:
            % figure
            % plot(tmp_traj_x, tmp_traj_y, 'ob')
            % hold on
            % plot(inter_traj_x, inter_traj_y, '.r')
            % hold off

            tmp_purs_x = track_data(s).trials(t).purs_x;
            tmp_purs_y = track_data(s).trials(t).purs_y;
            inter_purs_x = linspace(tmp_purs_x(1), tmp_purs_x(end), round(length(tmp_purs_x)/FR*SR));
            inter_purs_y = spline(tmp_purs_x, tmp_purs_y, inter_purs_x);

            % % proof that it works well:
            % figure
            % plot(tmp_purs_x, tmp_purs_y, 'ob')
            % hold on
            % plot(inter_purs_x, inter_purs_y, '.r')
            % hold off

            % for now, I am just deleting the error columns
            % TODO: adjust or re-calculate error in Matlab
            tmp_mat = [inter_traj_x; inter_traj_y; inter_purs_x; inter_purs_y].';
            for c = 1:size(tmp_mat, 2)
                track_data(s).("upsamp_trials")(t).(col_names{c}) = tmp_mat(:,c);
            end
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

% find trial start

start_pattern = {["S 26","S 27"], ["S 26", "", "S 27"]}; % possible sequence of
% start trials

strfind(start_pattern{1}(1), char({all_data_struct(s).urevent.type}))
% what is boundary? only sitmuli 
[all_data_struct(s).event.latency]
char({all_data_struct(s).event.type})


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



%% Find Peaks of Tracking Data and put them in Events


[~, index_max_traj] = findpeaks(track_data(s).upsamp_trials(.traj
% 1.1.2. target trajectory minima
[~, index_min_traj] = findpeaks(-track_data.task_a.(field_names{s}).(trial_name{t})(:,2));

plot(track_data.task_a.(field_names{s}).(trial_name{t})(:,1), track_data.task_a.(field_names{s}).(trial_name{t})(:,2))

%get pixel coordinates of extracted maxima to calculate euclidean distance
% max_traj_coords = [trialdata.traj_x_pix(index_max_traj), trialdata.traj_y_pix(index_max_traj)];

%create vector of 1 to append to matrix (marking maxima)
%id = zeros(size(index_max_traj,1),1)+1;

%create matrix for analyses
%max_traj = [max_traj_coords, index_max_traj, id];

%get pixel coordinates of extracted minima to calculate euclidean distance
% min_traj_coords = [trialdata.traj_x_pix(index_min_traj), trialdata.traj_y_pix(index_min_traj)];




%% Plot Tracking Data


%% Epoch around Peaks




