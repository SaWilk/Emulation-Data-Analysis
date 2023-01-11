%% Sanity Check for Sample Entropy Calculation

% Author: Saskia Wilken

% Only works if main script entropy.m has been run.


%% Check out random samples of data with entropy values to get a feeling for
% the measure. 

s = randi(30)
task = randi(2)
t = randi(72)
[val t] = min(tmp(1:2:end));
figure()
plot(track_data(s).upsamp_data.(task_names{task})(t).error(1:2000/SAMP_DIST), 'LineWidth', 2)
title(["First 500: ", track_data(s).upsamp_data.(task_names{task})(t).entropy(1), ...
    "Second segment: ", track_data(s).upsamp_data.(task_names{task})(t).entropy(2)])
xline(125)
