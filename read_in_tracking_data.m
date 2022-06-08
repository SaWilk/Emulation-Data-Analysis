% Read in Tracking Data


%% Empty

format compact
format long G
clear
clc


%% Set Paths

% Load in trajectory data

file_path = matlab.desktop.editor.getActiveFilename;
file_path = fileparts(file_path);
data_path = strcat(file_path, "\00_npz_files"); % get root path to data folder
subj_paths = genpath(data_path); % generate subfolder paths
subj_paths = strsplit(subj_paths, ";"); % split the string vector at ;

subj_paths(1) = []; % remove top level entry
subj_paths(end) = []; % remove final entry

% get subject ids from folder names
[~, tmp] = fileparts(subj_paths);
subj_ids = split(tmp, '_');
subj_ids = subj_ids(:, :, 1);
subj_ids = unique(subj_ids);
num_subj = length(subj_ids); % get num of subfolders

% initialize data structure
track_data = struct(); % create empty structure
clear tmp
% eeglab


%% Load Tracking data

warning("off","all");
% for now only for a
% TODO: Add in functionality for b
a_paths = subj_paths(endsWith(subj_paths, "_A"));
field_names = strcat("s_",subj_ids); % necessary to make subj ids valid
% field names
for p = 1:length(a_paths)
    if ~isempty(a_paths{p}) % check if there is data
        cd(a_paths{p});
        file_names_curr_subj = {dir('*traj_purs*.csv').name};
        for t = 1:length(file_names_curr_subj)
            trial_name = strcat("t_", num2str(t));
            track_data.task_a.(field_names{p}).(trial_name) = table2array(readtable(file_names_curr_subj{t}));
        end
    else fprintf('%s contains no data', a_paths{p});
    end
end
warning("on","all");

% what is contained in the different columns
legend = {'traj_x', 'traj_y', 'purs_x', 'purs_y', 'error', 'error_abs'}


%% Upsample Tracking Data to Match EEG Data SR

% decided on spline interpolation because of this convincing comparison:
% https://de.mathworks.com/help/matlab/ref/pchip.html
% spline works best to display oscillatory movements.
SR = 250; % sampling rate of EEG data
FR = 60; % frame rate of tracking data

for s = 1:length(field_names)
    if isfield(track_data.task_a, field_names{s}) % check if there is a field with the current subject id
        for t = 1:length(trial_name)
            tmp_traj_x = track_data.task_a.(field_names{s}).(trial_name{t})(:, 1);
            tmp_traj_y = track_data.task_a.(field_names{s}).(trial_name{t})(:, 2);
            inter_traj_x = linspace(tmp_traj_x(1), tmp_traj_x(end), round(length(tmp_traj_x)/FR*SR));
            inter_traj_y = spline(tmp_traj_x, tmp_traj_y, inter_traj_x);

            % % proof that it works well:
            % figure
            % plot(tmp_traj_x, tmp_traj_y, 'ob')
            % hold on
            % plot(inter_traj_x, inter_traj_y, '.r')
            % hold off

            tmp_purs_x = track_data.task_a.(field_names{s}).(trial_name{t})(:, 3);
            tmp_purs_y = track_data.task_a.(field_names{s}).(trial_name{t})(:, 4);
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
            track_data.task_a.(field_names{s}).(trial_name{t}) = [inter_traj_x; inter_traj_y; inter_purs_x; inter_purs_y].';
        end
    end

end


%% Find Peaks of Tracking Data and put them in a latency vector

% 1.1. target trajectory
% 1.1.1. target trajectory maxima/peaks
[~, index_max_traj] = findpeaks(track_data.task_a.(field_names{s}).(trial_name{t})(:,2));

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


%% Plot Tracking Data



%% Load EEG data

TMPEEG = pop_loadset(TMPEEG,'filename',[files2read{ind}(1:end-7) '_B_epoched'], 'filepath', char(savepath));


%% Put Tracking and EEG Data in the same structure


%% Epoch around Peaks

%%

all_data.tracking{s} = 1



