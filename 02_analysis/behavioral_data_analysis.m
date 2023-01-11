%% Behavioral Data Analyses

% Analyses the tracking data, error measures, their association, how error
% changes over the course of the trials, etc...
% Author: Saskia Wilken
% Creation Date: 18.10.2022

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
% input
track_data_dir = strjoin([output_dir, "03_parallelize_with_traj"], filesep);


%% Load Tracking Data

cd(track_data_dir)
load('all_tracking_data.mat', 'track_data') 

%% Load CDS data

eeglab;

% long epochs
cd("C:\wilken\Emulation-Data-Output\03_parallelize_with_traj"); % load non-epoched data
%list all *.set files in inputpath
file_names = dir('*.set');
%concatenate into one cell array
files2read = {file_names.name};
files2read = files2read(~contains(files2read, 'KMY6K')); % remove subject who 
% did not properly track during task B
% load eeg data epoched around peaks
for idx = 1:length(files2read)
    ALLEEG2(idx) = pop_loadset('filename',files2read{idx});
end

ALLEEG = ALLEEG2;
clear ALLEEG2


%% Plot Tracking Data

SUBJ = randi(30); % 4 % 24 % 23
TRIAL = randi(72); % 66 % 70 % 29
SR = 250;
NUM_SEC = 2;
SCREEN_HEIGHT = 1920;
FONTSIZE = 26;
LINEWIDTH = 2;
PLOTSIZE = [0, 0, 8, 5];

time_vec = linspace(0, (length(track_data(SUBJ).upsamp_data.task_a(TRIAL).traj_x)/250)*1000, ...
    length(track_data(SUBJ).upsamp_data.task_a(TRIAL).traj_x));

PLOT_TIME = 1:SR*NUM_SEC; % plot only the first two seconds

figure()
plot(time_vec(PLOT_TIME), ...
    track_data(SUBJ).upsamp_data.task_a(TRIAL).traj_y(PLOT_TIME)*SCREEN_HEIGHT, 'LineWidth', LINEWIDTH*1.5)
hold on
plot(time_vec(PLOT_TIME), ...
    track_data(SUBJ).upsamp_data.task_a(TRIAL).purs_y(PLOT_TIME)*SCREEN_HEIGHT, 'LineWidth', LINEWIDTH*1.5)
% xline((SR/2)*(1000/SR), '--', 'LineWidth', LINEWIDTH*1.5)
xline(find(abs(track_data(SUBJ).upsamp_data.task_a(TRIAL).purs_y(PLOT_TIME)*SCREEN_HEIGHT - ...
    track_data(SUBJ).upsamp_data.task_a(TRIAL).purs_y(1)*SCREEN_HEIGHT) > 5, 1)*1000/SR ...
    , '--', 'LineWidth', LINEWIDTH*1.5) % mark position where cursor first deviates 5 pixels from its starting point

legend({"trajectory", "pursuit"})

% titles and axes labels
title(strjoin(["Trajectory and pursuit of subject ", SUBJ, ", trial ", TRIAL], ""))
xlabel("time (ms)")
ylabel("position on screen (pixels)")
box off
% saveas(gcf,strjoin(["Scree plot of T-Value Sums of all positive clusters"]),'pdf')
set(gca,'fontsize', FONTSIZE*0.75, 'linewidth',LINEWIDTH);
set(gcf,'color','w','PaperUnits','inches','PaperPosition',PLOTSIZE)
filename = strjoin(["Trajectory and pursuit of subject ", SUBJ, ", trial ", TRIAL], "");
print('-djpeg', strjoin([filename, ".jpg"]), '-r600');
print('-dsvg', strjoin([filename, ".svg"]), '-vector');


%% Plot Heatmap of Error

cd("C:\wilken\Emulation-Data-Analysis-Journal\Behavioral_Data")
frame_rate = 60; % using non-interpolated data
error_heatmap(track_data, frame_rate)


%% Retrieve Occluded Error 

SR = 250;
clear vis_data occ_data
for s = 1:length(track_data)
    for t = 1:length(track_data(s).upsamp_data.task_b)
        i = 1;
        count_even = 1;
        count_odd = 1;
        while (i+1)*2*SR < length(track_data(s).upsamp_data.task_b(t).traj_y)
            % taking advantage of the fact that occlusion is every 2
            % seconds. 
            frame_of_change_1 = i*2*SR+1; % skipping ifrst sgment due to different error
            frame_of_change_2 = (i+1)*2*SR+1;
            clear cur_error
            cur_error = abs(track_data(s).upsamp_data.task_b(t).error_pix(frame_of_change_1:frame_of_change_2));
            cur_purs = track_data(s).upsamp_data.task_b(t).purs_y;
            cur_track = track_data(s).upsamp_data.task_b(t).traj_y;
            if rem(i, 2) == 0 % if we have an even number (which means occlusion is OFF)
                vis_data(t, count_even, s) = mean(cur_error, 'omitnan');
                count_even = count_even+1;
            elseif rem(i, 2) ~= 0 % if we have an odd number (which means occlusion is ON)
                occ_data(t, count_odd, s) = mean(cur_error, 'omitnan');
                count_odd = count_odd+1;
            else
                error("WTF")
            end % odd even condition
            i = i + 1;
        end % i loop
    end % trial loop
end % subject loop

occ_data(find(occ_data == 0)) = nan; % replace all tirlas with exactly 0 
% % error with NaN (to remove superfluous zeroes)
vis_data(find(vis_data == 0)) = nan;

vis_data(find(vis_data == 0)) = nan;


%% Retrieve Constant Error 
%     "start_constant": 23,
%     "end_constant": 24,
% custom: "S 15": End Trial

SR = 250;
clear const_data rand_data % easier code cause I am taking advantage of the 
% % fact that there is exactly one random and once constant segment per trial
for s = 1:length(track_data)
    EEG = ALLEEG(s);
    for t = 1:length(track_data(s).upsamp_data.task_a)
        cur_trial_events = EEG.event(find([EEG.event.trial_number] == t)); % weirdly enough, I have two triggers per peak. (S 41 and S 40 )
        % however, this shouldn't affect these analyses...
        start_rand1_lat = 0; % trial latency of start of rand1 is always 0
        for i = 1:size(cur_trial_events, 1)
            if strcmp(cur_trial_events(i).type, "S 23")
                end_rand1_lat = cur_trial_events(i).trial_latency-1;
                start_const_lat = cur_trial_events(i).trial_latency;
            elseif strcmp(cur_trial_events(i).type, "S 24")
                end_const_lat = cur_trial_events(i).trial_latency;
            end % we do not analyze random2
            cur_rand_err = track_data(s).upsamp_data.task_a(t).error_pix(start_rand1_lat:end_rand1_lat);
            cur_const_err = track_data(s).upsamp_data.task_a(t).error_pix(start_const_lat:end_const_lat);
            rand_data(t, s) = cur_rand_err;
            const_data(t, s) = cur_const_err;
        end % event
    end % trial
end % subjects

%% mean comparisons

BOOT_REP = 1000;
ALPHA = 0.001;

mean_per_subj{1} = squeeze(mean(mean(occ_data, 2, 'omitnan'), 1, 'omitnan'))';
mean_per_subj{2} = squeeze(mean(mean(vis_data, 2, 'omitnan'), 1, 'omitnan'))';
[stats, df, pvals, surrog_behav] = statcond(mean_per_subj, 'paired', 'on', 'method', 'bootstrap', ...
                'naccu', BOOT_REP, 'alpha', ALPHA, 'structoutput', 'on');
[h,p,ci,stats] = ttest(mean_per_subj{1}, mean_per_subj{2})
%          p: 1.95354304348172e-09 
%     tstat: 8.56599143061268
%        df: 29
mean(mean_per_subj{1})
%         0.106824454701965
mean(mean_per_subj{2})
%         0.0928398813482515

meanEffectSize(double(squeeze(surrog_behav)), Effect="cohen", Alpha=ALPHA)

%                     Effect                   ConfidenceIntervals          
%                _________________    ______________________________________
%     CohensD    -2.21088964143018    -2.40472442050528    -2.01852376515657


%% Prepare Scatterplot: Create one dataset with all the errors

clear behav_data behav_mean
cnt = 1;
for s = 1:length(ALLEEG)
    EEG = ALLEEG(s);
    for ep = 1:max([EEG.event.epoch])
        if isempty([EEG.event(round(median(find([EEG.event.epoch] == ep))) ...
                ).epoch_error]) | isempty([EEG.event(round(median(find([EEG.event.epoch] == ep))) ...
                ).pursuit_lat])
            behav_data(cnt,:) = nan(1,4);
        else
            behav_data(cnt,1) = [EEG.event(round(median(find([EEG.event.epoch] == ep)))).epoch_error] * 1080;
            behav_data(cnt,2) = [EEG.event(round(median(find([EEG.event.epoch] == ep)))).pursuit_lat];
            behav_data(cnt,3) = [EEG.event(round(median(find([EEG.event.epoch] == ep)))).epoch_error_z];
            behav_data(cnt,4) = [EEG.event(round(median(find([EEG.event.epoch] == ep)))).pursuit_lat_z];
        end
        cnt = cnt + 1;
    end
    one_subject_idx(s) = cnt;
    clear area_of_subj
    if s ~= 1
        area_of_subj = one_subject_idx(s-1):one_subject_idx(s)-1;
    else
        area_of_subj = 1:one_subject_idx(s)-1;
    end
    behav_mean(s, :) = mean(behav_data(area_of_subj,:), 1, 'omitnan');
end

cd(output_dir)
save('behav_data.mat', 'behav_data')
save('behav_mean.mat', 'behav_mean')
% DOn't run again, takes ages. Rather, load the saved dataset:


%% Plot scatterplot of pursuit latency and epoch error - Parameters

cd(output_dir)
load('behav_data.mat', 'behav_data') % this is just the error and pursuit latency
behav_data = rmmissing(behav_data);

clear ALLEEG D EEG
% make scatterplot of the raw values
[corrs, pvals] = corrcoef(behav_data(:,1:2));

x = behav_data(:,1);
y = behav_data(:,2);
% R = 25;                        % circle radius
rand_idx = randi(numel(D), [1,100000]);
R = prctile(D(rand_idx), 1)
D = pdist2([x y],[x y]);        % create combinations of distances (eucledian)
D = pdist2([x y],[x y]);        % create combinations of distances
C = sum(D<R);                   % how mayn points inside circle
figure()
scatter(x,y,25,C,'fill')
colorbar
% plot(behav_data(:,1), behav_data(:,2), 'o', size);
title('scatter plot of pursuit latency and epoch error')
subtitle(['no significant correlation, rho = ',num2str(round(corrs(2),4)), ...
    ', p = ', num2str(round(pvals(2), 4)), ' ,' , num2str(size(behav_data,1)), ' observation pairs'])
xlabel('epoch error in px')
ylabel('pursuit latency in ms')
set(gca,'fontsize', 14);
% mean for each subject, scatterplot

clear ALLEEG D EEG

% make scatterplot of the z-transformed values
corrs = corrcoef(behav_data(:, 3:4));

x = behav_data(:,3);
y = behav_data(:,4);
% R = 25;                        % circle radius
rand_idx = randi(numel(D), [1,100000]);
R = prctile(D(rand_idx), 1)
D = pdist2([x y],[x y]);        % create combinations of distances (eucledian)
C = sum(D<R);                   % how mayn points inside circle


%% Plot scatterplot of pursuit latency and epoch error - Plot Code

figure()
scatter(x,y,25,C,'fill')
colorbar
% plot(behav_data(:,1), behav_data(:,2), 'o', size);
title('scatter plot of pursuit latency and epoch error, z-transformed')
subtitle(['no significant correlation, ', num2str(size(behav_data,1)), ' observation pairs'])
xlabel('epoch error in px')
ylabel('pursuit latency in ms')
set(gca,'fontsize', 14);


%% Check what the pursuit latencies look like

% criteria:
% - pursuit peak with prominence = 0.05
% - a pursuit peak must follow the trajectory peak in the same epoch
% - pursuit and trajectory peaks must go in the same direction (up or down)
% - pursuit peak latency is between 80 and 300 ms. 
% - central peak of the epoch is not in the first half a second of a trial

disp(strjoin([round(no_purs/((purs_count+no_purs)/100),0), ...
    "% (", no_purs, " out of ", purs_count+no_purs, ...
    ") trajectory peaks are not followed by a pursuit peak."],""))
% at this time, 31472 out of 65627 trajectory peaks are not followed by a 
% pursuit peak (47.96%)
pursuit_lats = [];
all_pursuit_lats = [];
pursuit_lats_z = [];
all_pursuit_lats_z = [];
epoch_errors = [];
all_epoch_errors = [];
epoch_errors_z = [];
all_epoch_errors_z = [];

for s = 1:length(TMPEEG)
    % pursuit latencies raw
    pursuit_lats = [TMPEEG(s).event.pursuit_lat];
    % concatenate pursuit latencies of current subj with previous subj
    all_pursuit_lats = [all_pursuit_lats, pursuit_lats];
    % normalized pursuit latencies
    pursuit_lats_z = [TMPEEG(s).event.pursuit_lat_z];
    all_pursuit_lats_z = [all_pursuit_lats_z, pursuit_lats_z];
    % epoch error raw
    epoch_errors = [TMPEEG(s).event.epoch_error];
    all_epoch_errors = [all_epoch_errors, epoch_errors];
    % normalized epoch error
    epoch_errors_z = [TMPEEG(s).event.epoch_error_z];
    all_epoch_errors_z = [all_epoch_errors, epoch_errors_z];
end

% keep only one of the copies of the pursuit latency field values each. 
unique_lats_idx = find(diff(all_pursuit_lats));
% this method is not 100% failsafe - there is a chance that two following
% pursuit latencies are exactly the same and they would be missed by this
% algorithm. However, most pursuit lats should be in there. 
all_pursuit_lats_unique = all_pursuit_lats(unique_lats_idx);
unique_lats_idx_z = find(diff(all_pursuit_lats_z));
all_pursuit_lats_unique_z = all_pursuit_lats_z(unique_lats_idx_z);


%% Histograms of Pursuit Latency

BINWIDTH = 12;
range(all_pursuit_lats_unique)
PROMINENCE_THRESH = 0.05; 
ACCEPTABLE_LATS = [80, 300];
histogram(all_pursuit_lats_unique, ...
    round(range(all_pursuit_lats_unique)/BINWIDTH));
title(['Pursuit Latency'])
subtitle(strjoin(["Prominence thresh ", PROMINENCE_THRESH, ...
    ", peaks in same direction, latency lims between ", ...
    ACCEPTABLE_LATS, " ms", "binwidth = ", BINWIDTH, " ms"]))
mink(all_pursuit_lats_unique, 100)
% purs_count and length(all_pursuit_lats_unique) do not have the same size.
% which is odd. 

PROMINENCE_THRESH = 0.05; 
ACCEPTABLE_LATS = [80, 300];
histogram(all_pursuit_lats_unique_z, ...
    round(range(all_pursuit_lats_unique)/BINWIDTH));
title(['Pursuit Latency z-transformed per subject'])
subtitle(strjoin(["Prominence thresh ", PROMINENCE_THRESH, ...
    ", peaks in same direction, latency lims between ", ...
    ACCEPTABLE_LATS, " ms", "binwidth = ", BINWIDTH, " ms"]))
mink(all_pursuit_lats_unique, 100)


%% Plot epoch error distribution

figure()
BINWIDTH = 0.02;
histogram(all_epoch_errors, round(range(all_epoch_errors)/BINWIDTH));
title(['Epoch Error'])
subtitle(strjoin(["All epochs, all errors over time, binwidth = ", ...
    BINWIDTH*100, " % of screen"]))

figure()
histogram(all_epoch_errors_z, round(range(all_epoch_errors)/BINWIDTH));
title(['Epoch Error z-transformed'])
subtitle(strjoin(["All epochs, all errors over time, binwidth = ", ...
    BINWIDTH*100, " % of screen"]))

