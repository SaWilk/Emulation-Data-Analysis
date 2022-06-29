function [plot_of_trigger_and_traj] = plot_trigger_on_traj(EEG_sets, traj_struct, subject, trial, task)
%plot_trigger_on_traj plots both all trigger and trajectory in same plot
%   requires the output of merge_eeg_and_tracking_data script

% EEG_sets = ALLEEG;
% traj_struct = track_data;
% subject = 1;
% trial = 1;
% task = 'task_a';

% select events

cor_task_idx = find(strcmp({EEG_sets(subject).event.task}, task));
cor_trial_idx = find([EEG_sets(subject).event(cor_task_idx).trial_number] == trial);
triggers = {EEG_sets(subject).event(cor_trial_idx).type};
trigger_labels = unique(triggers, 'stable');
time_points = [EEG_sets(subject).event(cor_trial_idx).trial_latency_ms]; 

x = traj_struct(subject).upsamp_data.(task).traj_x;
y =  traj_struct(subject).upsamp_data.(task).traj_y;

for trig_types = 1:length(trigger_labels)
    handle_idx{trig_types} = find(strcmp(triggers, trigger_labels{trig_types}));
end

colors = {'m','b','r','g','k', 'c'};

plot(x, y, 'linewidth', 1.5)
for trig = 1:length(trigger_labels)
    hold on
    handle{trig} = xline(time_points(handle_idx{trig}), 'color', colors{trig}, 'linewidth', 1.5);
end
% only keep first elements since we only want one legend per plot type. 
handle_legend_idx = cellfun(@(v)v(1),handle);
legend(handle_legend_idx, trigger_labels);
hold off
