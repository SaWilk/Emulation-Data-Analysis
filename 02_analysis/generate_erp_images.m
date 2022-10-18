% ERP analyses

% creates plots and performs further analyses to investigate
% ERP

% author: Saskia Wilken
% creation date: 21.09.2022 

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
mean_matrices_path = strjoin([output_dir, 'mean_matrices'], filesep);
mkdir(mean_matrices_path);
mean_matrices_peaks_epoched_path = strjoin([mean_matrices_path, 'peaks'], filesep);
mkdir(mean_matrices_peaks_epoched_path);
peak_erp_plots_dir = strjoin([parent_dir, "plots", "erp_plots", "peaks"], filesep);
mkdir(peak_erp_plots_dir);
% output
csd_dir = strjoin([output_dir, "csd_transform"], filesep);
% input
epochs_plus_error_dir = strjoin([output_dir, "07_epochs_with_extra_fields"], filesep);
% cluster based permutation test
neighbors_dir = strjoin([parent_dir_2, "Emulation-Data-Input", "EEG_files"], filesep);


%% Load CDS data

eeglab;

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


%% Get zlims

clear zlims

all_zlims = nan(2, length(ALLEEG));
% % make sure all plots have the same color scale
for s = 1:length(ALLEEG) %length(diff_labels)
    all_zlims(:, s) = [prctile( ...
        ALLEEG(s).CSD_data, 35, 'all'), prctile( ...
        ALLEEG(s).CSD_data, 65, 'all')]';
    disp(s)
end

zlims =  [min(min(all_zlims(1,:))), max(max(all_zlims(2,:)))];


%% Remove all unimportant event field entries

for s = 1:length(ALLEEG)
    ALLEEG(s).data = ALLEEG(s).CSD_data;
end

clear z_ALLEEG TMPEEG z_mean_struct mean_struct z_tmpEEG CONDEEG CONDEEG2 TMPEEG2
sum(arrayfun(@(s) length(s.epoch), ALLEEG)); % how many epochs we need. 
s = 1;
while s <= length(ALLEEG)
    EEG = ALLEEG(s);
    del = 0;
    event_idx = zeros(length(EEG.epoch),1);
    for ep = 1:length(EEG.epoch)
        event_idx(ep) = get_epoch_center(EEG, ep);
%         if ~(strcmp(EEG.event(ev-del).type, 'S 40') | strcmp(EEG.event(ev-del).type, 'S 50'))
%             EEG.event(ev-del) = [];
%             del = del +1;
%         end
    end
    log_idx = unfind(event_idx, length(EEG.event));
    EEG.event(~log_idx) = [];
    EEG = eeg_checkset(EEG);
    ALLEEG(s) = EEG;
    s = s + 1;
end
% WATCH OUT: FROM NOW ON get_epoch_center won't work anymore since the
% EEG.epoch structure is no longer aligned with the EEG.event structure. 


%% Create Condition-Specific ALLEEG Datasets

for s = 1:length(ALLEEG)
    EEG = ALLEEG(s);
    EEG = pop_selectevent( EEG, 'artifact_error',{'VALID'},'deleteevents',...
        'off','deleteepochs','on','invertepochs','off');
    ALLEEG(s) = EEG;
end
CONDEEG = ALLEEG;
% remove occluded VISIBLE
for s = 1:length(CONDEEG)
    EEG = CONDEEG(s);
    EEG = pop_selectevent( EEG, 'OCCL',{'OFF'},'deleteevents',...
        'off','deleteepochs','on','invertepochs','off');
    CONDEEG(s) = EEG;
end
CONDEEG2 = ALLEEG;
% remove visible OCCLUDED
for s = 1:length(CONDEEG2)
    EEG = CONDEEG2(s);
    EEG = pop_selectevent( EEG, 'OCCL',{'ON'},'deleteevents',...
        'off','deleteepochs','on','invertepochs','off');
    CONDEEG2(s) = EEG;
end

% remove random CONSTANT
% CONDEEG = ALLEEG;
% 
% for s = 1:length(CONDEEG)
%     EEG = CONDEEG(s);
%     EEG = pop_selectevent( EEG, 'TRAJ',{'RANDOM1', 'RANDOM2'},'deleteevents',...
%         'off','deleteepochs','on','invertepochs','off');
%     CONDEEG(s) = EEG;
% end
% CONDEEG2 = ALLEEG;
% % remove constant RANDOM1
% for s = 1:length(CONDEEG2)
%     EEG = CONDEEG2(s);
%     EEG = pop_selectevent( EEG, 'TRAJ',{'RANDOM2', 'CONST'},'deleteevents',...
%         'off','deleteepochs','on','invertepochs','off');
%     CONDEEG2(s) = EEG;
% end


%% Put all subjects' data in the same dataset

% initialize dataset
TMPEEG = CONDEEG(1);
TMPEEG.data = normalize(CONDEEG(1).data,2, 'zscore');
% put all epochs of all subjects in one dataset
% get all data and all event structures
for s = 2:length(CONDEEG)
    % normalize each subjects' data
    z_tmpEEG = CONDEEG(s);
    % normalize data
    z_tmpEEG.data = [];%MODEEG(s).data(:,:,cond_ind(s).(cond_labs{lab}));
    z_tmpEEG.data = normalize(CONDEEG(s).data,2, 'zscore');
    % get location where new data should go
%     data_idx = size(TMPEEG.data, 3)+1:size(TMPEEG.data, 3)+size(z_tmpEEG.data, 3);
    % concatenate subject data in one dataset
    TMPEEG.data = cat(3, z_tmpEEG.data, TMPEEG.data);
%      TMPEEG.data(:,:,data_idx) =
    % get epoch vectors
    epoch_vec = TMPEEG.event(end).epoch + [z_tmpEEG.event.epoch];
    % get event indices
    event_idx = size(TMPEEG.event,2)+1:size(TMPEEG.event,2)+size(z_tmpEEG.event, 2);
    % add events from previous subject to TMPEEG
    TMPEEG.event(:, event_idx) = z_tmpEEG.event;
    for ev = 1:length(event_idx)
        TMPEEG.event(:, event_idx(ev)).epoch = epoch_vec(ev);
    end
    TMPEEG.epoch(:, size(TMPEEG.epoch,2)+1:...
        size(TMPEEG.epoch,2)+size(z_tmpEEG.epoch, 2)) = z_tmpEEG.epoch;
    epoch_vec = [];
end
TMPEEG.trials = length(TMPEEG.epoch);

clear vec
vec = [TMPEEG.event.epoch_error]*1080;
for i = 1:length(vec)
    TMPEEG.event(i).('epoch_error_pix') = vec(i);
end


TMPEEG2 = CONDEEG2(1);
TMPEEG2.data = normalize(CONDEEG2(1).data,2, 'zscore');
% put all epochs of all subjects in one dataset
% get all data and all event structures
for s = 2:length(CONDEEG2)
    % normalize each subjects' data
    z_tmpEEG = CONDEEG2(s);
    % normalize data
    z_tmpEEG.data = [];%MODEEG(s).data(:,:,cond_ind(s).(cond_labs{lab}));
    z_tmpEEG.data = normalize(CONDEEG2(s).data,2, 'zscore');
    % get location where new data should go
%     data_idx = size(TMPEEG.data, 3)+1:size(TMPEEG.data, 3)+size(z_tmpEEG.data, 3);
    % concatenate subject data in one dataset
    TMPEEG2.data = cat(3, z_tmpEEG.data, TMPEEG2.data);
%      TMPEEG.data(:,:,data_idx) =
    % get epoch vectors
    epoch_vec = TMPEEG2.event(end).epoch + [z_tmpEEG.event.epoch];
    % get event indices
    event_idx = size(TMPEEG2.event,2)+1:size(TMPEEG2.event,2)+size(z_tmpEEG.event, 2);
    % add events from previous subject to TMPEEG
    TMPEEG2.event(:, event_idx) = z_tmpEEG.event;
    for ev = 1:length(event_idx)
        TMPEEG2.event(:, event_idx(ev)).epoch = epoch_vec(ev);
    end
    TMPEEG2.epoch(:, size(TMPEEG2.epoch,2)+1:...
        size(TMPEEG2.epoch,2)+size(z_tmpEEG.epoch, 2)) = z_tmpEEG.epoch;
    epoch_vec = [];
end
TMPEEG2.trials = length(TMPEEG2.epoch);
clear vec
vec = [TMPEEG2.event.epoch_error]*1080;
for i = 1:length(vec)
    TMPEEG2.event(i).('epoch_error_pix') = vec(i);
end
count = 1;
clear chan_lab_cell
clear chan_lab_vec
% diff_labels = {'diff_all_zero', 'diff_vis_occ', 'diff_const_rand1', ...
%     'diff_const_rand2', 'diff_rand1_rand2'};
diff_labels = {'diff_all_zero', 'diff_vis_occ', 'diff_const_rand1'};
% without rand2
% get vector of all interesting channels
% for lab = 2:length(diff_labels)
%     for lat = 1:length(lat_labs)
%         for m = 1:2
%             chan_lab_cell{count} = index_struct.(diff_labels{lab}).(num2str(lat_labs{lat}))(1).(minmax{m}){:};
%             chan_idx_vec(count) = index_struct.(diff_labels{lab}).(num2str(lat_labs{lat}))(2).(minmax{m})(:);
%             count = count + 1;
%         end
%     end
% end
% chan_lab_cell_plus = {chan_lab_cell{:},chan_lab_min{:},chan_lab_max{:}}
% chan_idx_vec_plus = [chan_idx_vec,min_chan_idx,max_chan_idx]


%% Actual erp image plotting and parameters
% sorting_cond = {'epoch_error', 'pursuit_lat', 'epoch_error_z', 'pursuit_lat_z'};
% cond = 1;
smooth = round(length(TMPEEG.epoch)/15,0); % number of trials
dec_fac = smooth/2; % ratio of trials in to plot out. 
% if the smoothing width is larger than twice the decimation factor, no
% info is lost. 
plot_type = 1; % 1 = channel, 0 = component
project_channel = [[]];
sort_type = {'S 40' 'S 50'};
% sort_win = [-0.1 0.1];
% align = 0;
unit = 'z-transformed \muV';
% chan = 32;
% chan_lab = chan_lab_cell{chan};
chan_lab = 'Cz';
chan_no = find(strcmp({ALLEEG(s).chanlocs.labels}, chan_lab));
alpha = 0.001;
% chan_no = chan_idx_vec(chan);
% eeg_c = lab;

% sort_event_field = sorting_cond{cond};
sort_event_field = 'epoch_error';
% [min_val, max_val] = bounds(TMPEEG.data, 'all')
% prctile(TMPEEG.data(:,:,1:10000), [1, 10, 25, 50, 75, 90, 100], 'all')
% prctile(ALLEEG(s).data(:,:,1:1000), [1, 10, 25, 50, 75, 90, 100], 'all')
% prctile(CONDEEG(s).data(:,:,1:1000), [1, 10, 25, 50, 75, 90, 100], 'all')
% prctile(z_tmpEEG.data(:,:,1:1000), [1, 10, 25, 50, 75, 90, 100], 'all')

% for cond = 1:length(sorting_cond)
%     for chan = 1:length(chan_idx_vec)
%         for eeg_c = 1:length(cond_labs)
        % dynamically changing image parameters
% 
% [ALLEEG TMPEEG index] = eeg_store(ALLEEG, TMPEEG);
% eeglab redraw
time_range = [-200, 750];
if strcmp(sort_event_field, 'epoch_error')
    plot_sortvar  = 'off';
else
    plot_sortvar  = 'on';
end

if strcmp(sort_event_field, 'epoch_error')
    trial_ax_lab = 'epoch error in fraction of screen';
    trial_ax_ticks = 0:0.01:0.2;
elseif strcmp(sort_event_field, 'pursuit_lat')
    trial_ax_lab = 'pursuit latency in ms';
    trial_ax_ticks = 0:25:325;
end
% plotting of ERP image
% eeg_condition = cond_labs{eeg_c};
eeg_condition = 'all';
figure; 
title = strjoin(["All subjects, normalized & weighted, csd, Smoothing factor "...
    smooth,", decimate " dec_fac ", chan " chan_lab],"");
pop_erpimage(TMPEEG,plot_type, chan_no,[[]], title,smooth, dec_fac, ...
    sort_type,[],sort_event_field ,'yerplabel', unit,'erp', 'on', 'erpalpha', 0.001,...
    'cbar','on','cbar_title',unit, 'limits', [time_range, zlims], 'caxis', zlims, ...
    'topo', { chan_no TMPEEG.chanlocs TMPEEG.chaninfo }, ...
    'img_trialax_label', trial_ax_lab, 'erp_grid')..., ... % sets the x axis ticks in units of the sorting variable
%     'noplot', plot_sortvar); 
colormap 'parula'
sgtitle(strjoin(["contrast: " eeg_condition, ", sorted by: ", ...
    sort_event_field],""), 'Interpreter', 'None');

eeg_condition = 'occluded';
smooth = round(length(TMPEEG2.epoch)/15,0); % number of trials
dec_fac = smooth/2; % ratio of trials in to plot out. 
% plotting of ERP image
figure; 
title = strjoin(["All subjects, normalized & weighted, csd, Smoothing factor "...
    smooth,", decimate " dec_fac ", chan " chan_lab],"");
pop_erpimage(TMPEEG2,plot_type, chan_no,[[]], title,smooth, dec_fac, ...
    sort_type,[],sort_event_field ,'yerplabel', unit,'erp', 'on', 'erpalpha', 0.001,...
    'cbar','on','cbar_title',unit, 'limits',  [time_range, zlims], 'caxis', zlims, ...
    'topo', { chan_no TMPEEG2.chanlocs TMPEEG2.chaninfo }, ...
    'img_trialax_label', trial_ax_lab)%, ... % sets the x axis ticks in units of the sorting variable
    %'noplot', plot_sortvar); 
colormap 'parula'
sgtitle(strjoin(["contrast: " eeg_condition, ", sorted by: ", ...
    sort_event_field],""), 'Interpreter', 'None');

