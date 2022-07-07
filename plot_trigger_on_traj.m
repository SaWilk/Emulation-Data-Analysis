function [plot_of_trigger_and_traj] = plot_trigger_on_traj(EEG_sets, traj_struct, subject, trial, task, plot_purs)
%plot_trigger_on_traj plots both all trigger and trajectory in same plot
%   requires the output of merge_eeg_and_tracking_data script

% epoch_lim, epoch_cen  % moved from list of input arguments
% EEG_sets = eeg_struct;
% traj_struct = track_data;
% subject = 4;
% trial = 50;
% task = 'task_a';
% plot_purs = true;
% epoch_lim = [-1000,750]; % epoch limits in ms
% epoch_cen = "S 40"; % epoch center

clear handle trig_idx trigger_labels trigger_names trigger_dict triggers handle_legend_idx

trigger_dict = {["exp_start", "S 11"], ["fix", "S 12"], ["trial_start_L", "S 13"], ...
    [ "trial_start_R", "S 14"], ["trial_end", "S 15"], ["pause_start", "S 16"], ...
    ["pause_end", "S 17"], ["exp_end", "S 18"], ["button", "S 19"], ...
    ["occlusion", "S 20"], ["reappear", "S 21"], ["start_constant", "S 23"], ...
    ["end_constant", "S 24"], ["start_startvec", "S 26"], ["end_startvec", "S 27"], ...
    ["C_too_early", "S 28"], ["C_just_right", "S 29"], ["C_too_late", "S 30"], ...
    ["peak", "S 40"]};

if isnumeric(task)
    task_names = {'task_a', 'task_b', 'task_c'};
    task = task_names{task};
elseif ischar(task)
else
    error("task must be either numeric or character")
end

% select events
% find the indicies of the task in the events structure
cor_task_idx = find(strcmp({EEG_sets(subject).event.task}, task));
[min_idx, max_idx] = bounds(cor_task_idx);
% and within that, the indicies of the events for the trial
cor_trial_idx = find([EEG_sets(subject).event(cor_task_idx).trial_number] == trial);
% then get all the trigger codes within that trial
triggers = {EEG_sets(subject).event(cor_task_idx(cor_trial_idx)).type};
% get the unique labels
trigger_labels = unique(triggers, 'stable');
% and then the trial latency in ms of the events for plotting
time_points = [EEG_sets(subject).event(cor_task_idx(cor_trial_idx)).trial_latency_ms]; 
% get idx of central event
% epoch_cen_times = time_points(find(strcmp(triggers, epoch_cen)));

% get epoch boundaries
% tmp = repmat(epoch_lim, size(epoch_cen_times,2), 1);
% epoch_bounds_times = bsxfun(@plus, epoch_cen_times', tmp);

% get upsampled trajectory
x = traj_struct(subject).upsamp_data.(task)(trial).traj_x;
y =  traj_struct(subject).upsamp_data.(task)(trial).traj_y;
y_purs =  traj_struct(subject).upsamp_data.(task)(trial).purs_y;


for trig_types = 1:length(trigger_labels)
    % get the indices of triggers of a certain type and put in handle_idx cell
    handle_idx{trig_types} = find(strcmp(triggers, trigger_labels{trig_types}));
    % translate the trigger codes into trigger names
    trigger_names = cellfun(@(v)v(1), trigger_dict);
    %overwrite trigger labels with names
    trigger_labels{trig_types} = trigger_names(strcmp(cellfun(@(v)v(2), trigger_dict), trigger_labels{trig_types}));
end

% prepare color palette
colors = {'m', 'b', 'r', 'g', 'k', 'c', 'y'};

% v = [0 0; 1 0; 1 1; 0 1];
% f = [1 2 3 4];

% plot traj
plot(x, y, 'linewidth', 1.25)
% for each of the different trigger types..
for trig = 1:length(trigger_labels)
    hold on
    % draw an xline of the corresponding time point and use the index
    handle{trig} = xline(time_points(handle_idx{trig}), 'color', colors{trig}, 'linewidth', 1.25);
end
if plot_purs
    plot(x, y_purs, 'linewidth', 1.25, 'color', 'r')
end

% [min_y, max_y] = bounds(y);
% y_bounds = repmat([min_y, max_y], size(epoch_bounds_times,1), 1)
% rectangle('Position',[epoch_bounds_times(1,1), y_bounds(1,1), epoch_bounds_times(1,2), y_bounds(1,2)])
% only keep first elements since we only want one legend per plot type. 
handle_legend_idx = cellfun(@(v)v(1),handle);

% TODO: Add in functionality that the epoch starts and ends are visible. 
% h=fill([0,1000,1,0],[0,0,2,2],'red');
% h.FaceAlpha=0.3;
legend(handle_legend_idx, trigger_labels);
hold off

