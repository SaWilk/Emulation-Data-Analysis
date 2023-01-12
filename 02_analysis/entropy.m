%% Entropy Calculation 
%
% - Calculates sample entropy for each pursuit of each
% subject using track_data.mat
% - Tests whether the first 500 ms of each trial have significantly different 
% sample entropy values than the rest
% - Plots the result
% - Verison that calculates the Sample Etrnopy per Subject included

% Author: Saskia Wilken
% Creation date: 10.11.2022

%% empty everything

format compact
format long G
clear
clc
close all

%% Define Paths

% file location
file_path = matlab.desktop.editor.getActiveFilename;
file_path = fileparts(file_path);
addpath(file_path);

% get data paths for parent dirs
filepath_parts = strsplit(file_path, filesep);
parent_dir_2 = filepath_parts(1:end-2);
parent_dir_2 = strjoin(parent_dir_2, filesep);
output_dir = strjoin([parent_dir_2, "Emulation-Data-Output"], filesep);
track_data_dir = strjoin([output_dir, "03_parallelize_with_traj"], filesep);


%% Read in tracking data

cd(track_data_dir);

load all_tracking_data.mat


%% Calculate Sample Entropy for each segment and put in matrix

% Set Parameters
M = 2; % embedded dimension, 2 is the default
TAU = 1; % downsampling factor, can be used for multi-scale entropy 
% % calculations (then loop through tau)
task_names = fieldnames(track_data(1).upsamp_data);
SAMP_DIST = 4;
count = 1;
entropy_vals = [];
R = 0.2; % tolerance threshold
% "to make sure we know the variance of the data and it is the same for
% every sample, we z-transform the data." Source: DOI: 10.3390/e21060541

% Entropy Calculation per Trial
for s = 1:size(track_data,2)
    for task = 1:2
        for t = 1:size(track_data(s).upsamp_data.(task_names{task}),2)
            cur_data = track_data(s).upsamp_data.(task_names{task})(t).error;
            for seg = 1:2
                switch seg
                    case 1
                        segment = cur_data(1:round((500/SAMP_DIST)));
                    case 2
                        segment = cur_data(round(501/SAMP_DIST):round(2000/SAMP_DIST));
                end
                segment_z = normalize(segment, 1, 'zscore');
                entropy_vals(count, seg) = SampEn(M, R, segment_z, TAU);
            end % segment
            track_data(s).upsamp_data.(task_names{task})(t).entropy = entropy_vals(count, :);
            count = count + 1;
        end % trial
    end % task
end % subject

% Uncomment this if you want to calculate one entropy value per subject
% only (for statistical tests)
% count = 1;
% clear subject_entropy
% % version with one value per subject
% for s = 1:size(track_data,2)
%     for task = 1:2
%         for t = 1:size(track_data(s).upsamp_data.(task_names{task}),2)
%             tmp = [track_data(s).upsamp_data.(task_names{task})(t).entropy];
%             subject_entropy(count, :, s) = tmp;
%             count = count+1;
%         end % trial
%     end % task
%     count = 1;
% end % subject

% segment1 = squeeze(median(subject_entropy(:, 1, :), 1, "omitnan"));
% segment2 = squeeze(median(subject_entropy(:, 2, :), 1, "omitnan"));


%% Plot Resulting Sample Entropy distributions

data1 = entropy_vals(:,1);
data2 = entropy_vals(:,2);
% data1 = segment1;
% data2 = segment2;
% data1 = subject_entropy(:, 1, 3);
% data2 = subject_entropy(:, 2, 3);
N_BINS = 80;
ALPHA = 0.05;
EDGES = linspace(min([data1, data2],[], 'all'), max([data1, data2], [], 'all'), N_BINS+1);
[p,h,stats] = signrank(data1,data2, 'Alpha', ALPHA)

% Plot Code
figure()
histogram(data1, 'FaceAlpha', 0.3, 'FaceColor', 'b', 'BinEdges', EDGES)
hold on
histogram(data2, 'FaceAlpha', 0.3, 'FaceColor', 'r', 'BinEdges', EDGES)
xline(median(data1, 'omitnan'), 'LineWidth', 2, 'Color', 'b')
xline(median(data2, 'omitnan'), 'LineWidth', 2, 'Color', 'r')
% xline(mean(data1, 'omitnan'), 'LineWidth', 2, 'Color', 'g') %
% demonstrates that mean operations would underestimate the difference
% severely
% xline(mean(data2, 'omitnan'), 'LineWidth', 2, 'Color', 'y')
hold off
box off
% text
legend({"1 - 500 ms", "501 - 2000 ms", "median 1 - 500 ms", "median 501 - 2000 ms"})
title(strjoin("Comparison of Sample Entropy (m = 2, r = 0.2*std) of Tracking Error of First 500 ms to the Following 1500 ms"))
subtitle(strjoin([length(data1) - sum(isnan(data1)), ...
    "observations, signed rank test result: Z = ", round(stats.zval, 2), ", p = ", p ...
    , " r = Z/sqrt(N) = ", round(stats.zval/sqrt(numel(data1)), 4)]))
xlabel("sample entropy")
ylabel("frequency")
set(gca, 'FontSize', 26, 'LineWidth', 2)
% save figure
set(gcf,'color', 'w', 'PaperUnits','inches','PaperPosition',[0 0 13 8])
filename = strjoin(["Comparison of Sample Entropy first 500 to following 1500 ms.jpg"], "");
print('-djpeg', filename, '-r600');
filename = strjoin(["Comparison of Sample Entropy first 500 to following 1500 ms.svg"], "");
print('-dsvg', filename, '-vector');

