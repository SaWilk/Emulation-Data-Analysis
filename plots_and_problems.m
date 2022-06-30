%% Plots and Issues

% author: Saskia Wilken
% creation date: 21.06.2022

% this script only works when the data was read in by
% read_in_tracking_data.m and is loaded in the workspace

s = 4; % picking some random trials and subjects
t = 5;

tmp_traj = track_data(s).trials(t).traj_x;
tmp_traj(:,2) = track_data(s).trials(t).traj_y;
for i = 1:size(tmp_traj,1)-1
    distance(i) = sqrt((tmp_traj(i, 1) - tmp_traj(i+1, 1))^2 + (tmp_traj(i, 2) - tmp_traj(i+1, 2))^2);
end
% what is shit: distance along traj should be the same, but is not. 
figure()
subplot(2,1,1)
plot(1:length(distance), distance, "linewidth", 1.5)
title("Distance travelled per frame")
hold on
plot(1:length(distance), [diff(distance), diff(distance(end-1:end))])
legend({"absolute distance", "difference to previous frame"});
hold off
subplot(2,1,2)
plot(tmp_traj(:,1), tmp_traj(:,2), "linewidth", 1.5)



%% Test Plot Peaks Triggers

s = 5
t = 1
task = 'task_b'

figure()
plot_trigger_on_traj(ALLEEG, track_data, s, t, task)

% gap array is 30 events times upsampling factor before constant and after
% task =
%     'task_b'
% Index exceeds the number of array elements. Index must not exceed 26.
% Error in plot_trigger_on_traj (line 32)
%     handle{trig} = xline(time_points(handle_idx{trig}), 'color', colors{trig}, 'linewidth', 1.5); 

