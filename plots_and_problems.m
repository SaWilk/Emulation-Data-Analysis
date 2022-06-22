%% Plots and Issues

% author: Saskia Wilken
% creation date: 21.06.2022

% this script only works when the data was read in by
% merge_eeg_and_tracking_data.m and is loaded in the workspace

%% target speed alongt traj


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


% Note: This is NOT an issue, because the target speed was calculated using
% radius stuff and the eucledian distance always connects two points with a
% straight line, which is not as accurate as the way it was calculated for
% the experiment. So the problem is in the plot, not in the data!

%% Plot Tracking Data

s = 2; % pick subject
trials = 1:10; % pick some trials

figure()
for t = trials
    trial_starts = find(strcmp({all_data_struct(s).track_event_peaks.type}, "S 27"));
    trial_ends = find(strcmp({all_data_struct(s).track_event_peaks.type}, "S 15"));
    trial_x_only = all_data_struct(s).track_event_peaks(trial_starts(t):trial_ends(t));
    peaks_only = find(strcmp({trial_x_only.type}, "S 40"));
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

    % sanity check that trigger and eeg are aligned

    subplot(5, 2, t)
    plot(traj_x_temp, traj_y, "linewidth", 1.5)
    hold on
    peak_handle = xline(peak_latencies,"linewidth", 1.5);
    const_handle = xline(const_latencies, "r", "linewidth", 1.5);
    xline(const_latencies_end, "r", "linewidth", 1.5)
    start_handle = xline(trial_latencies_start, "g", "linewidth", 1.5);
    xline(trial_latencies_end, "g", "linewidth", 1.5);
    hold off
    legend([peak_handle(1), const_handle, start_handle], {"peaks", "const", "trial"});

end
 
% PROBLEM: this shows that the constant traj triggers are not really around 
% the constant traj, which explains the need for calculating the onset of the
% constant traj from the data

