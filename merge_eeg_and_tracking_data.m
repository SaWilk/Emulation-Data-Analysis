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
paths{1} = subj_paths(endsWith(subj_paths, "_A"));
paths{2} = subj_paths(endsWith(subj_paths, "_B"));

% what is contained in the different columns
col_names = {'traj_x', 'traj_y', 'purs_x', 'purs_y', 'error', 'error_abs'};
task_names = {"task_a", "task_b"};
for task = 1:2
    for p = 1:length(paths{task})
        if ~isempty(paths{task}{p}) % check if there is data
            cd(paths{task}{p});
            file_names_curr_subj = {dir('*traj_purs*.csv').name};
            track_data(p).("subject") = extractBefore([file_names_curr_subj{1}],'_');
            track_data(p).("group") = 'control'; % TODO: Adjust if there is ever another group
            %             track_data(p).("condition") = {'task_A'}; % TODO: adjust for other tasks
            track_data(p).("path").(task_names{task}) = paths{task}{p};
            for t = 1:length(file_names_curr_subj)
                tmp_mat = table2array(readtable(file_names_curr_subj{t}));
                for c = 1:size(col_names,2)
                    track_data(p).("trials").(task_names{task})(t).(col_names{c}) = tmp_mat(:,c);
                end
                clear tmp_mat
            end
        else fprintf('%s contains no data', a_paths{p});
        end
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
    for task = 1:2
        for t = 1:size(track_data(s).trials.(task_names{task}), 2)

            tmp_traj_y = track_data(s).trials.(task_names{task})(t).traj_y;
            % create time vec
            time_vec = linspace(0, length(tmp_traj_y)/FR*1000, round(length(tmp_traj_y))).';
            % intepolated time vec
            inter_time_vec = linspace(time_vec(1), time_vec(end), round(length(time_vec)/FR*SR));
            % interpolate using spline interpolation
            inter_traj_y = spline(time_vec, tmp_traj_y, inter_time_vec);

            % proof that it works well:
            %         figure
            %         plot(no_points, inter_traj_y, 'ob')
            %         hold on
            %         plot(time_vec, tmp_traj_y, '.r')
            %         hold off
            % same for pursuit
            tmp_purs_y = track_data(s).trials.(task_names{task})(t).purs_y;
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
                tmp_mat
                track_data(s).("upsamp_data").(task_names{task})(t).(col_names{c}) = tmp_mat(:,c);
            end
        end
        clear tmp_traj_y tmp_purs_y tmp_mat
    end
end


%% Calculate difference between pursuit and trajectory


for s = 1:length(subj_ids)
    for task = 1:2
        for t = 1:size(track_data(s).trials.(task_names{task}), 2)
            cur_traj = track_data(s).upsamp_data.(task_names{task})(t).traj_y;
            cur_purs = track_data(s).upsamp_data.(task_names{task})(t).purs_y;
            cur_err = abs(cur_traj - cur_purs);
            track_data(s).upsamp_data.(task_names{task})(t).("error") = 
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
% there is only merged a and b available, so this line is no longer needed.
% Was for separating a and b task files.
% files2read_a = files2read(~cellfun(@isempty, regexp(files2read,'_A_')));

EEG = pop_loadset('filename',files2read);
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
eeg_struct = ALLEEG;


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
% keep only subjects with eeg and tracking data
copy_data = track_data(tmp_idx);


% Deprecated: Put tracking and eeg data in same strucure
% TODO: write them directly in a 3D Matrix instead of in fields inside a
% structure
% field_names = fieldnames(copy_data);
% % nope, it does not seem to be possible to concatenate structures
% % without a loop
% for s = 1:size(copy_data,2)
%     if strcmp(copy_data(s).subject, all_data_struct(s).subject)
%         for f = 1:length(field_names)
%             % put tracking data in all data struct
%             all_data_struct(s).(field_names{f}) = copy_data(s).(field_names{f});
%         end
%     else
%         disp(strcat(all_data_struct(s).subject, " has no matching tracking data"));
%     end
% end


%% Reject Trials Behaviorally




%% Run The Same Script that Adriana used to reject Trials TODO


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

for s = 1:size(eeg_struct,2)

    % copy urevent field to other field
    eeg_struct(s).track_event = eeg_struct(s).urevent;
    % copy event field to convenient format
    event = eeg_struct(s).track_event;
    % sanity check of triggers - are we DEALING WITH TASK A EVEN???
    % TODO: Alter check once pipeline is adjusted for Task B as well.
    event_cat = categorical({event.type});
    categories(event_cat);
    % if you should ever feel the need to count cats, there you go:
    countcats(event_cat); % this gives you how many triggers of each type there are.

    % is this task c?
    if any(strcmp({event.type}, repmat("S 30",size({event.type},1), size({event.type},2) ))) | ...
            any(strcmp({event.type}, repmat("S 28",size({event.type},1), size({event.type},2) )))
        warning(strcat(['Caution! Subject', eeg_struct(s).subject, 'seems to be a task C trigger file. Skipping']));
        continue
        % if it is not, is this task b?
        %     elseif any(strcmp({event.type}, repmat("S 20",size({event.type},1), size({event.type},2) ))) | ...
        %            any(strcmp({event.type}, repmat("S 21",size({event.type},1), size({event.type},2) )))
        %         warning(strcat(['Caution! Subject', eeg_struct(s).subject, ' seems to be a task B trigger file. Skipping']));
        %         continue
    end

    % delete all triggers before S11
    start_pract = find(strcmp({event.code}, "New Segment"))+1;
    end_pract = find(strcmp({event.type}, "S 11"))-1;
    event([start_pract(1):end_pract(1), start_pract(2):end_pract(2)]) = [];

    % insert end trial triggers

    last_trigger = "none";

    % loop through elements of event_types.
    i = 0;
    while i < length({event.type})
        i = i + 1;
        start_exp = find(strcmp({event.code}, "New Segment"));

        % if an event type is either 16 or 12...
        if (strcmp(event(i).type, "S 16") | strcmp(event(i).type, "S 12"))
            % make sure the first entry after S 11 is skipped (since the first fix cross
            % is not a trial end) via checking whether i is shortly after S
            % 11
            if any(start_exp(1):start_exp(1)+4 == i) | any(start_exp(2):start_exp(2)+4 == i)
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
    % and before S 11
    event = [event(1:start_exp(2)-1), event(start_exp(2)-1), event(start_exp(2):end)];
    event(start_exp(2)-1).type = 'S 15';
    event(start_exp(2)-1).code = 'inserted';
    % copy the new event_types struct into track event
    eeg_struct(s).track_event = event;

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
    [eeg_struct(s).track_event.trial_latency] = tmp{:};

    % update event
    event = eeg_struct(s).track_event;
    % initialize trial number field
    tmp = num2cell(repmat(9999, size(eeg_struct(s).track_event)));
    [event.("trial_number")] = tmp{:};

    if length(trial_starts) == length(trial_ends)
        % loop through trial start indices
        for idx = 1:length(trial_starts)
            % set per trial latencies of events in trial
            trial_event_latencies = [event(trial_starts(idx):trial_ends(idx)).trial_latency] - ...
                event(trial_starts(idx)).trial_latency;
            trial_event_latencies_cell =  num2cell(trial_event_latencies);
            [event(trial_starts(idx):trial_ends(idx)).trial_latency] = trial_event_latencies_cell{:};
            % add field: trial number
            if idx <= 72
            trial_num = idx;
            elseif 72 < idx <= 72*2
                trial_num = idx - 72;
            else
                trial_num = idx - 72*2;
            end
            trial_idx = num2cell(repmat(trial_num, size(trial_event_latencies_cell)));
            [event(trial_starts(idx):trial_ends(idx)).("trial_number")] = trial_idx{:};
            % add field: trial latency ms
            trial_event_latencies_ms = num2cell(trial_event_latencies/250*1000);
            [event(trial_starts(idx):trial_ends(idx)).("trial_latency_ms")] = trial_event_latencies_ms{:};
        end
    else
        warning(strcat(['Subject ', eeg_struct(s).subject, ...
            ' has an unequal number of trial start and trial end triggers!', ...
            ' No trial latencies will be computed. ']));
        NOLATENCIES = eeg_struct(s).subject;
    end

    eeg_struct(s).track_event = event;

end


% Notes: 91L3HA is excluded because recording started too late.
% WM87B is excluded cause there is no behavioral data.
% TODO: Warning: Subject ZV583 has an unequal number of trial start and trial end triggers! No trial latencies will be computed.


%% Find Peaks of Tracking Data and put them in Events

% "S 40": Peak

for s = 1:size(eeg_struct,2)

    event = eeg_struct(s).track_event;
    if ~strcmp(NOLATENCIES, eeg_struct(s).subject)
        for task = 1:2
            % get eeg triggers of current task
            exp_ind = [find(strcmp({event.code}, "New Segment"))-1, length({event.type})];
            exp_ind(1) = 1;
            % add field to events: task
            task_labels = cellstr(repmat(task_names{task}, [1, length(exp_ind(task):exp_ind(task+1))]));
            [event(exp_ind(task):exp_ind(task+1)).("task")] = task_labels{:};
            event_cur_task = event(find(strcmp({event.task}, task_names{task})));

            for t = 1:size(copy_data(s).upsamp_data.(task_names{task}),2)

                current_trial_traj = copy_data(s).upsamp_data.(task_names{task})(t).traj_y;

                % get peak indices within trial trajectory data
                [~, index_max_traj] = findpeaks(current_trial_traj);
                [~, index_min_traj] = findpeaks(-current_trial_traj);

                % locate current trial in eeg event
                cur_trial_idx = find([event_cur_task.trial_number] == t);

                %
                %             if length(trial_starts) ~= length(trial_ends)
                %                 warning(strcat(['Subject ', eeg_struct(s).subject, ...
                %                     ' does not have an equal number of start and end trial ' ...
                %                     'triggers in ', char(task_names{task}), '. Skipping.']))
                %                 break
                %             end

                % get the latencies of the current trial in the eeg data
                current_start_latency = event_cur_task(cur_trial_idx(1)).latency;
                current_end_latency = event_cur_task(cur_trial_idx(end)).latency;
                % add current start latency to the indicies of the peaks (which are in the same sampling rate
                % as the eeg signal due to upsampling)
                current_trial_peak_latencies = sort([index_max_traj; index_min_traj]);
                current_peak_latencies = current_start_latency + current_trial_peak_latencies;

                % for each of the current peak latencies...
                for idx = 1:length(current_peak_latencies)

                    % find the event that is one event before the peak
                    % (which is the largest negative latency)
                    tmp = [[event_cur_task.latency] - current_peak_latencies(idx)];
                    current_event_idx = max(find(tmp <= 0 ));

                    % insert the peak event markers at position in question
                    % (copy event before it and adjust its values)
                    event_cur_task = [event_cur_task(1:current_event_idx-1), ...
                        event_cur_task(current_event_idx), event_cur_task(current_event_idx:end)];
                    % replace the copied rows' type with S 40, making sure the
                    % second event is overwritten
                    event_cur_task(current_event_idx+1).type = 'S 40';
                    event_cur_task(current_event_idx+1).code = 'inserted';
                    event_cur_task(current_event_idx+1).latency = current_peak_latencies(idx);
                    event_cur_task(current_event_idx+1).trial_latency = current_trial_peak_latencies(idx);
                    event_cur_task(current_event_idx+1).trial_latency_ms = event_cur_task(current_event_idx+1).trial_latency/250*1000;
                    event_cur_task(current_event_idx+1).trial_number = t;

                end
            end
            % copy the event structure into a temporary field in eeg_struct
            eeg_struct(s).(task_names{task}) =  event_cur_task;
        end
        % concatenate the two fields so eeglab can work with them
        eeg_struct(s).track_event_peaks = [eeg_struct(s).(task_names{1}), eeg_struct(s).(task_names{2})];
    else
        warning(strcat(['Skipping subject ', eeg_struct(s).subject, ...
            ' due to unequal number of trial start and trial end triggers!', ...
            ' And therefore missing trial latencies. ']));
    end
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


%% Remove historical fields; can be commented out here easily for re-tracing changes

% tmp_struct = eeg_struct;

tmp_struct = all_data_struct;

for s = 1:size(tmp_struct, 2)

    tmp_struct(s).event = tmp_struct(s).track_event_peaks;

end
tmp_struct = rmfield(tmp_struct, ["track_event", "track_event_peaks"]);

% isequal(fieldnames(ALLEEG), fieldnames(tmp_struct))


%% Replace Constant Traj Trigger from .vmrk with calculated from .npz

%align constant trajectories of several trials
%compare constant trajectory and trial trajectory at every point
%calculate variance of difference
%should be minimal at constant trajectory

% [shift, start_ind, end_ind] = align_const_traj(const_traj, trial_traj);

% TOOD: once you are done, tell Adriana


%% Save sets

% This does not work... :(
% go over this
% https://eeglab.org/tutorials/ConceptsGuide/Data_Structures.html

out_path = strcat([file_path, '\01_merged_data\']);
EEG_file_name_suffix = '_EEG';
track_file_name_suffix = 'all_tracking_data.mat';

for s = 1:size(copy_data, 2)

    file_name = strcat([tmp_struct(s).subject, EEG_file_name_suffix]);
    pop_saveset(tmp_struct(s), 'filename', file_name, 'filepath', out_path);

end
file_name = strcat([out_path, track_file_name_suffix]);
save(file_name, 'copy_data');

