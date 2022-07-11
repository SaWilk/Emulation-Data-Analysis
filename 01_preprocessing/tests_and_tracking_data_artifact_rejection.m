%% Tests and Tracking Data Artifact Rejection
%
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

% TODO: Test if it looks good when the constant triggers are moved by the 
% duration of the gap array. It is 30 points times upsampling factor before 
% % constant and after. Put some trials in a subplot to check

%% Trial Rejection 

% Requires that parallelize tracking and eeg data is loaded. 

% Visual artifact rejection 

tmp = cellfun(@(v)v(1),error_cell);
tmp_vec = [tmp{:}];
% [val, idx] = maxk(tmp_vec,10);

% check sope of error distribitopn
% plot error distribution
figure()
plot(sort(tmp_vec))
cutoff = mean(tmp_vec) + std(tmp_vec)*3.5;
hold on
idx = find(tmp_vec > cutoff)
suspicious_traj = cellfun(@(v)v(2),{error_cell{idx}});
yline(cutoff)
hold off

tmp2 = cellfun(@(v)v(1),error_max_cell);
tmp2_vec = [tmp2{:}];
% plot error max distribution
figure()
plot(sort(tmp2_vec))
cutoff = mean(tmp2_vec) + std(tmp2_vec)*3.5;
hold on
idx2 = find(tmp2_vec > cutoff)
tmp_cell = cellfun(@(v)v(2),{error_max_cell{idx2}});
suspicious_traj(end+1:end+length(tmp_cell)) = tmp_cell;

yline(cutoff)
hold off
suspicious_traj{1}
suspicious_traj{2}
while i <= length(suspicious_traj)

    if sum(cellfun(@isequal,suspicious_traj, repmat(suspicious_traj(i),1, length(suspicious_traj)))) > 1
        suspicious_traj(i) = [];
    end
    i = i +1 ;
end
% for some odd reason, they are exactly the same.. 
% plot suspicious trials
for i = 1:length(suspicious_traj)
    figure()
plot_trigger_on_traj(eeg_struct, track_data, suspicious_traj{i}(1), suspicious_traj{i}(3), suspicious_traj{i}(2), true)
title(strcat(["Traj and Pursuit of subject", suspicious_traj{i}(1), ...
    ", trial ", suspicious_traj{i}(2), ", task ", suspicious_traj{i}(3)]));
end

% Joystick moved very far from zero point. 
figure()
plot_trigger_on_traj(eeg_struct, track_data, 31, 19, 2, true)
title(strcat(["Traj and Pursuit of subject", 31, ...
    ", trial ", 40, ", task ", 2]));

% subject 18 seems to have lost motivation on task b. I am not sure what
% he is doing, but it doesn't look like tracking. 
for i = [1:5, 21:25, 65:70]
    figure()
    plot_trigger_on_traj(eeg_struct, track_data, 18, i, 2, true)
    title(strcat(["Traj and Pursuit of subject", 18, ...
    ", trial ", i, ", task ", 2]));
end

% subject 17 does not really like to move the joystick too much during task
% a, but seems to have improved during task b. I would keep him. 
for i = [1:5, 21:25, 65:70]
    figure()
    plot_trigger_on_traj(eeg_struct, track_data, 17, i, 2, true)
    title(strcat(["Traj and Pursuit of subject", 17, ...
    ", trial ", i, ", task ", 2]));
end

% % looks like it starts to get really weird at 1.5

% % trials to exclude: 
% % subject, task, trial
[24, 2, 46];
[18, 2, 1];
[16, 2, 70];
[30, 2, 2];
[19, 1, 22];
[18, 2, 21];
[18, 2, 2];
[3, 2, 38];
[2, 1, 53];
[20, 2, 40];
[20, 2, 10];
[18, 2, 1];
[17, 1, 20];
[17, 1, 19];
[16, 2, 70];
[18, 2, 20];
% all of subject 18 task B






