%% Latency Analysis for Cluster Occluded Visible Negative 1

% Calculates the latency between occluded and visible condition in one
% cluster
% Author: Saskia Wilken
% Creation date: 09.12.2022


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
parent_dir = filepath_parts(1:end-1);
parent_dir = strjoin(parent_dir, filesep);
parent_dir_2 = filepath_parts(1:end-2);
parent_dir_2 = strjoin(parent_dir_2, filesep);
output_dir = strjoin([parent_dir_2, "Emulation-Data-Output"], filesep);
track_data_dir = strjoin([output_dir, "03_parallelize_with_traj"], filesep);
% output
csd_dir = strjoin([output_dir, "csd_transform"], filesep);


%% Load CDS data

eeglab;
close all

% long epochs
cd(csd_dir);
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


%% Set Parameters

load clust_chans.mat

contrasts = {"all", "constrand", "occvis"};
conditionxclust = [{[1, 2]}, ...
    {3:length(clust)}];
cluster_labels = {"constrand_neg_3", "constrand_neg_1", "occvis_neg_1", ...
    "occvis_neg_3", "occvis_pos_1", "occvis_pos_2", "occvis_pos_3"};
clusters = clust;
INTER_SAMPLE_DIST = 4;
BASE_INT_MS = 500;
TIME_RANGE_MS = [BASE_INT_MS+200:INTER_SAMPLE_DIST:500+BASE_INT_MS];
TIME_RANGE_SAMP = round(TIME_RANGE_MS/INTER_SAMPLE_DIST);
TIME_RANGE_MS = TIME_RANGE_MS - BASE_INT_MS;
BOOT_REP = 5000;
CLUST = clusters{3}
for lab = 1:length(CLUST)
    chan_idx(lab) = find(strcmp({ALLEEG(1).chanlocs.labels}, CLUST{lab}));
end


%% Remove all unimportant event field entries

% copy csd data to data to avoid working with the worng data
for s = 1:length(ALLEEG)
    ALLEEG(s).data = ALLEEG(s).CSD_data;
    tmpEEG(s) = rmfield(ALLEEG(s), "CSD_data");
end

ALLEEG = tmpEEG;
clear tmpEEG

clear z_ALLEEG TMPEEG z_mean_struct mean_struct CONDEEG CONDEEG2 TMPEEG2
sum(arrayfun(@(s) length(s.epoch), ALLEEG)); % how many epochs we need.

% actual removal of unimportant event fields
s = 1;
while s <= length(ALLEEG)
    EEG = ALLEEG(s);
    del = 0;
    event_idx = zeros(length(EEG.epoch),1);
    for ep = 1:length(EEG.epoch)
        event_idx(ep) = get_epoch_center(EEG, ep);
    end
    log_idx = unfind(event_idx, length(EEG.event));
    EEG.event(~log_idx) = [];
    EEG = eeg_checkset(EEG);
    ALLEEG(s) = EEG;
    s = s + 1;
end
% WATCH OUT: FROM NOW ON get_epoch_center won't work anymore since the
% EEG.epoch structure is no longer aligned with the EEG.event structure.

% Remove epochs within first 500 ms

for s = 1:length(ALLEEG)
    EEG = ALLEEG(s);
    EEG = pop_selectevent( EEG, 'artifact_error',{'VALID'},'deleteevents',...
        'off','deleteepochs','on','invertepochs','off'); % remove the first
    % 500 ms peaks from data because we want to associate eeg with behavioral data
    ALLEEG(s) = EEG;
end

CONTRAST = contrasts{3};

CONDEEG = ALLEEG;
% remove occluded VISIBLE
for s = 1:length(CONDEEG)
    EEG = CONDEEG(s);
    EEG = pop_selectevent( EEG, 'task', {'task_b'},'deleteevents',...
        'off','deleteepochs','on','invertepochs','off');
    EEG = pop_selectevent( EEG, 'OCCL',{'OFF'},'deleteevents',...
        'off','deleteepochs','on','invertepochs','off');
    CONDEEG(s) = EEG;
    cond_lab = "visible";
    [CONDEEG(s).event.("subject")] = deal(s);
end

CONDEEG2 = ALLEEG;
% remove visible OCCLUDED
for s = 1:length(CONDEEG2)
    EEG = CONDEEG2(s);
    EEG = pop_selectevent( EEG, 'OCCL',{'ON'},'deleteevents',...
        'off','deleteepochs','on','invertepochs','off');
    CONDEEG2(s) = EEG;
    cond_lab2 = "occluded";
    [CONDEEG2(s).event.("subject")] = deal(s);
end


%% Find the peaks in the data

close all
clear cur_data_2 cur_data_1 all_subj_data_1 all_subj_data_2

BOOT_REP = 5000;
ALPHA = 0.05;

for s = 1:size(CONDEEG, 2)
    cur_data_1 = mean(mean(CONDEEG2(s).data(chan_idx, :, :), 1), 3);
    cur_data_2 = mean(mean(CONDEEG(s).data(chan_idx, :, :), 1), 3);
    all_subj_data_1(s, :) = cur_data_1;
    all_subj_data_2(s, :) = cur_data_2;
end

% use the minimum and maximum of the data per subject as thresholds for the
% prominence
lims(1) = min(all_subj_data_1(s, :));
lims(2) = max(all_subj_data_1(s, :));
PROM_THRESH = (lims(2) - lims(1))*0.33;

whole_epoch_time = -500:INTER_SAMPLE_DIST:750;
all_subj_locs_1 = nan(size(CONDEEG, 2), 3);
all_subj_locs_2 = nan(size(CONDEEG, 2), 3);


%% Plot all subjects' data with peaks indicated

figure()
for s = 1:size(CONDEEG, 2)
    cur_data_1 = all_subj_data_1(s, :);
    cur_data_2 = all_subj_data_2(s, :);
    [pks1, locs1] = findpeaks(cur_data_1, whole_epoch_time, 'MinPeakProminence', PROM_THRESH);
    val_idx = locs1 > min(TIME_RANGE_MS) & locs1 < max(TIME_RANGE_MS);
    pks1_i = pks1(val_idx);
    locs1_i = locs1(val_idx);
    clear val_idx
    [pks2, locs2] = findpeaks(cur_data_2, whole_epoch_time, 'MinPeakProminence', PROM_THRESH);
    val_idx = locs2 > min(TIME_RANGE_MS) & locs2 < max(TIME_RANGE_MS);
    locs2_i = locs2(val_idx);
    pks2_i = pks2(val_idx);
    clear val_idx
    subplot(5, 6, s)
    plot(TIME_RANGE_MS, cur_data_1(TIME_RANGE_SAMP))
    hold on
    plot(TIME_RANGE_MS, cur_data_2(TIME_RANGE_SAMP))
    plot(locs1_i, pks1_i, '*')
    plot(locs2_i, pks2_i, '*')
    hold off
    all_subj_locs_1(s, 1:length(locs1_i)) = locs1_i;
    all_subj_locs_2(s, 1:length(locs2_i)) = locs2_i;
end
legend({"occluded", "visible"})
sgtitle("Part of Time of Visible Occluded Negative 1 Cluster, ERP peak per subject, detected by findpeaks function")


%% Only keep subjects who have only one clear peak per condition

for s = 1:size(CONDEEG, 2)
    if sum(~isnan(all_subj_locs_1(s,:))) ~= 1 || sum(~isnan(all_subj_locs_2(s,:))) ~= 1
        del_idx(s) = logical(1);
    else
        del_idx(s) = logical(0);
    end
end
all_subj_locs_1(del_idx, :) = [];
all_subj_locs_2(del_idx, :) = [];
data{1}= all_subj_locs_1(:,1)';
data{2}= all_subj_locs_2(:,1)';


%% Calculate statistics

[stats, df, pvals, surrog] = statcond(data, 'paired', 'off', 'method', 'bootstrap', ...
    'naccu', BOOT_REP, 'alpha', ALPHA, 'structoutput', 'on');
round(mean(all_subj_locs_1(:, 1)))
round(mean(all_subj_locs_2(:, 1)))

tmp = meanEffectSize(data{1}, data{2}, Effect="cohen", Alpha=ALPHA);
