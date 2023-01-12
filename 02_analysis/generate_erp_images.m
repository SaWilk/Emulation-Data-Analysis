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
parent_dir_2 = filepath_parts(1:end-2);
parent_dir_2 = strjoin(parent_dir_2, filesep);
output_dir = strjoin([parent_dir_2, "Emulation-Data-Output"], filesep);
% output
csd_dir = strjoin([output_dir, "csd_transform"], filesep);


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
for s = 1:length(ALLEEG) 
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

% Remove the behavioral events of the first 500 ms of each trial
for s = 1:length(ALLEEG)
    EEG = ALLEEG(s);
    EEG = pop_selectevent( EEG, 'artifact_error',{'VALID'},'deleteevents',...
        'off','deleteepochs','on','invertepochs','off');
    ALLEEG(s) = EEG;
end

% Uncomment this if you want to investigate the contrast OCCLUDED - VISIBLE
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

% Uncomment this if you want to investigate the contrast CONSTANT - RANDOM
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
    z_tmpEEG.data = [];
    z_tmpEEG.data = normalize(CONDEEG(s).data,2, 'zscore');
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

% create pixel unit epoch error
clear vec
vec = [TMPEEG.event.epoch_error]*1080;
for i = 1:length(vec)
    TMPEEG.event(i).('epoch_error_pix') = vec(i);
end


%% same for the second condiiton

TMPEEG2 = CONDEEG2(1);
TMPEEG2.data = normalize(CONDEEG2(1).data,2, 'zscore');
% put all epochs of all subjects in one dataset
% get all data and all event structures
for s = 2:length(CONDEEG2)
    % normalize each subjects' data
    z_tmpEEG = CONDEEG2(s);
    % normalize data
    z_tmpEEG.data = [];
    z_tmpEEG.data = normalize(CONDEEG2(s).data,2, 'zscore');
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

% create pixel unit epoch error
clear vec
vec = [TMPEEG2.event.epoch_error]*1080;
for i = 1:length(vec)
    TMPEEG2.event(i).('epoch_error_pix') = vec(i);
end


%% Set ERPimage parameters

smooth = round(length(TMPEEG.epoch)/15,0); % number of trials
dec_fac = smooth/2; % ratio of trials in to plot out. 
% if the smoothing width is larger than twice the decimation factor, no
% info is lost. 
plot_type = 1; % 1 = channel, 0 = component
project_channel = [[]];
sort_type = {'S 40' 'S 50'};
unit = 'z-transformed \muV';
chan_lab = 'P3';
chan_no = find(strcmp({ALLEEG(s).chanlocs.labels}, chan_lab));
alpha = 0.001;

sort_event_field = 'pursuit_lat'; % change depending on which sorting variable you want
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


%% plotting of ERP image

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


%% erpimage of only one condition

eeg_condition = 'occluded'; % switch name here depending on which condition to plot
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

