%% Make ERP Plot

% Author: Saskia Wilken
% creation Date: 28.06.2022


%% Empty

format compact
format long G
clear
clc


%% Folders and Dependencies

% add path and start EEGlab toolbox
% addpath('R:\AG-Beste-Orga\Skripts\Toolbox\eeglab2021.0');
eeglab;
close;

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

% add functions in parent dir
addpath([char(parent_dir) filesep 'functions']);

% set input & output directory
input_dir = strjoin([parent_dir_2, "Emulation-Data-Output\03_parallelize_with_traj"], filesep);
output_dir = strjoin([parent_dir_2, "Emulation-Data-Output"], filesep);
output_dir_epoched = strjoin([parent_dir_2, "Emulation-Data-Output\04_epoched"], filesep);
output_dir_baseline = strjoin([parent_dir_2, "Emulation-Data-Output\05_baseline"], filesep);
output_dir_AR = strjoin([parent_dir_2, "Emulation-Data-Output\06_artifact_rejection"], filesep);

subdir_const_rand = strjoin([output_dir_epoched, "const_rand"], filesep);
subdir_occl = strjoin([output_dir_epoched filesep "occl_nonoccl"], filesep);
% out dirs for peaks for multiple epoch lengths
subdir_peaks = strjoin([output_dir_epoched filesep 'peaks_-500_750'], filesep);

study_dir = strjoin([output_dir, "study"], filesep);

savepath_baseline_const_rand = strjoin([output_dir_baseline, "const_rand"], filesep);
savepath_baseline_occl = strjoin([output_dir_baseline, 'occl_nonoccl'], filesep);
% out dirs for peaks for multiple epoch lengths
savepath_baseline_peaks = strjoin([output_dir_baseline, 'peaks_-500_750'], filesep);

output_dir_AR_const = strjoin([output_dir_AR, 'const'], filesep);
output_dir_AR_occl = strjoin([output_dir_AR, 'occl'], filesep);
% out dirs for peaks for multiple epoch lengths
output_dir_AR_peaks = strjoin([output_dir_AR, 'peaks_-500_750'], filesep);

mean_matrices_path = strjoin([output_dir, 'mean_matrices'], filesep);
mkdir(mean_matrices_path);
mean_matrices_peaks_epoched_path = strjoin([mean_matrices_path, 'peaks'], filesep);
mkdir(mean_matrices_peaks_epoched_path);

% neighbors dir
neighbors_dir = strjoin([parent_dir_2, "Emulation-Data-Input", "EEG_files"], filesep);


%% Load Epoched data 

eeglab;

% long epochs
cd(output_dir_AR_peaks);
%list all *.set files in inputpath
file_names = dir('*.set');
%concatenate into one cell array
files2read = {file_names.name};
% load all eeg files
% ALLEEG2 = struct();
for idx = 1:length(files2read)
     ALLEEG2(idx) = pop_loadset('filename',files2read{idx});
%      ALLEEG2(idx) = EEG;
end

% ALLEEG
% 
% pop_study


%% Calculate Means Across Subjects 

% TODO: Put the means in a structure. 
% One mean each epoch across all subjects and conditions
% One mean each epoch across all subjects and separate for constant/random
% One mean each epoch across all subjects and separate for occluded/visible

for i = 1:size(ALLEEG, 2) 

    EEG = epochs(i);
    mean_struct.all(:,:,i) = mean(EEG.data,3);
    mean_struct.all_time_vec = EEG.times;
%     EEG = OCC_ALLEEG(i);
%     % TODO: assign here the data from the different occlusion event structs
%     mean_struct.occ(:,:,i) = mean(EEG.data,3);
%     mean_struct.vis(:,:,i) = mean(EEG.data,3);
%     
%     mean_struct.occ_time_vec = EEG.times;

end


%% Save the Mean Matrices for Convenient Loading

cd(mean_matrices_peaks_epoched_path)
save("mean_struct.mat", 'mean_struct')


%% Load Previously Generated Mean Structures

cd(mean_matrices_peaks_epoched_path)

load('mean_struct.mat')


%% Plot ERPs as grand average


avg_peak_dist = 300;

epoch_mean = mean(mean_struct.all(:,:,:),3);

[~, min_ind] = min(epoch_mean,[],'all');
[min_ind(1), min_ind(2)] = ind2sub(size(epoch_mean), min_ind)
% Average ERP Plots
figure()
% long epochs

[epoch_lims(1), epoch_lims(2)] = bounds(epochs(8).times)
plot(mean_struct.all_time_vec, epoch_mean);
title(strcat(['ERP of all subjects averaged across epochs around peaks, n ~= ', num2str( size(epochs(1).data,3)), ' per subject']));
subtitle(strcat(['all channels, epochs from ', num2str(epoch_lims(1)), ' to ' num2str(epoch_lims(2))]))
hold on
h(1) = xline(0);
h(2) = xline(mean_struct.all_time_vec(min_ind(2)), 'r');
h(4:5) = xline([-avg_peak_dist,avg_peak_dist], 'b');
h(3) = xline(epoch_lims(1)+500, 'g')
legend(h, {'peak', 'error negativity peak',  'baseline period', 'average distance of next peak'});
hold off

% Single Channel ERP
% Identify the electrode that shows the biggest amplitude in the ERP
neg_time = 60; 

[~, neg_time_idx] = min(abs(mean_struct.all_time_vec - neg_time));
[~, min_chan_idx] = min(epoch_mean(:, neg_time_idx));

figure()

chan_lab = EEG.chanlocs(min_chan_idx).labels;
[epoch_lims(1), epoch_lims(2)] = bounds(epochs(8).times)
plot(mean_struct.all_time_vec, mean(mean_struct.all(min_chan_idx,:,:),3));
title(strcat(['ERP of all subjects averaged across epochs around peaks, channel ', chan_lab]));
subtitle(strcat(['epochs from ', num2str(epoch_lims(1)), ' to ' num2str(epoch_lims(2))]))
hold on
h(1) = xline(0);
h(2) = xline(mean_struct.all_time_vec(min_ind(2)), 'r');
h(4:5) = xline([-avg_peak_dist,avg_peak_dist], 'b');
h(3) = xline(epoch_lims(1)+500, 'g')
legend(h, {'peak', 'error negativity peak',  'baseline period', 'average distance of next peak'});
hold off


%% Create Topoplot of the time points of interest 

epoch_mean = mean(mean_struct.all(:,:,:),3);

latencies_onset = linspace(0, 725, 30);
latencies_offset = latencies_onset + 25;
% doc topoplot: https://rdrr.io/cran/erpR/man/topoplot.html

% extract means across 50 ms periods for topoplot
for i = 1:length(latencies_onset)

    [~, lat_idx_onset] = min(abs(mean_struct.all_time_vec - latencies_onset(i)));
    [~, lat_idx_offset] = min(abs(mean_struct.all_time_vec - latencies_offset(i)));
    topo_mean(:, i) = mean(epoch_mean(:,lat_idx_onset:lat_idx_offset), 2);

end

zlims = [min(topo_mean, [], 'all'), max(topo_mean, [], 'all')]';
EEG = ALLEEG(1)

figure()
for i = 1:length(latencies_onset)

    subplot(6, 5, i)
    topoplot(topo_mean(:, i), EEG.chanlocs, 'maplimits',  zlims)
    title('epochs around peaks all subjects all trials, long epochs')
    subtitle(strcat(['latency range: average of ', num2str(latencies_onset(i)), ' to ', ...
        num2str(latencies_offset(i)), ' ms']))
    colorbar()

end


%% Load Study

[ STUDY ALLEEG ] = pop_loadstudy('filename', 'study_peaks_epoched.study', 'filepath', study_dir)


%% Perform cluster-Based Permutation Test TODO Attempt with EEGLAB

% https://www.fieldtriptoolbox.org/tutorial/cluster_permutation_timelock/
% not sure whether I am using the fieldtrip parameters correctly or whether
% I am supposed to specify tzhings like alpha in eeglab or fieltrip syntax.
% try it out.
% [STUDY neighbors] = std_prepare_neighbors( STUDY, ALLEEG)
load(strjoin([neighbors_dir, 'fieldtrip_EEG_KJP_neighbours_61.mat'], filesep))
% stat_cond
% TODO: Add a condition in which all data is zero by copying the data and 

eeglab redraw

% Load Previously Generated Mean Structures

cd(mean_matrices_peaks_epoched_path)

load('mean_struct.mat')

minimal_latency = mean_struct.all_time_vec(min_ind(2)); % the point in which 
% % the signal is as negative as possible

peak_eeg = mean_struct.all;
zero_eeg = zeros(size(peak_eeg));

peak_across_subj = mean(peak_eeg, 3);
peak_across_subj_min_lat = peak_across_subj(:,minimal_latency);


% create data for a zero condition

% for subj = 1:length(ALLEEG)
%     zero_eeg(subj).("data") = zeros(size(ALLEEG(subj).data));
% end


For example, to compute mean ERPs statistics from a
            STUDY for epochs of 800 frames in two conditions from three
            groups of 12 subjects:
 
            >> erp_data_for_all_subjects_possibly_from_STUDY = { [800x12] [800x12] [800x12];... % 3 groups, cond 1
                    [800x12] [800x12] [800x12] };  % 3 groups, cond 2
[pcond, pgroup, pinter, statscond, statsgroup, statsinter] = std_stat(...
    {peak_eeg; zero_eeg}, 'groupstats', 'off', ...
    'condstats', 'off', 'method', 'permutation', 'naccu', '1000', 'alpha', '0.9', ...
    'mcorrect', 'none', 'mode', 'fieldtrip', 'fieldtripnaccu', 'numrandomization', ...
    'fieldtripalpha', '0.9', 'fieldtripmethod', 'montecarlo', ...
    'fieldtripmcorrect','cluster', 'fieldtripclusterparam', {'clusterstatistic', ...
    'maxsum', 'statistic', 'indepsamplesT'}, 'fieldtripchannelneighbor', 'neighbours' )
% 'fieldtripclusterparam' string or cell array for optional parameters
%                              for cluster correction method, see function
%                              ft_statistics_montecarlo for more information.
%    'fieldtripchannelneighbor' - Fieldtrip channel neighbour structure for 
%                                 cluster correction method, see function
%                                 std_prepare_neighbors for more information.
% use this to find the electrode cluster 
% In EEGLAB 12, press the statistics button in the channel or component STUDY 
% plotting interface. Then you can select cluster statistics (click on 
% Fieldtrip, then select monte-carlo statistics, then select cluster correction).
% It works as the other STUDY statistics.

%% Cluster Based Permutation Test Attempt with Fieldtrip

% Script by Adriana Böttcher
% Adapted by Saskia Wilken


%% add fieldtrip path
eeglab;
close;

% file location
addpath 'C:\Users\swilk\AppData\Roaming\MathWorks\MATLAB Add-Ons\Collections\FieldTrip'
file_path = matlab.desktop.editor.getActiveFilename;
file_path = fileparts(file_path);
addpath(file_path);
filepath_parts = strsplit(file_path, filesep);
parent_dir = filepath_parts(1:end-1);
parent_dir = strjoin(parent_dir, filesep);
parent_dir_2 = filepath_parts(1:end-2);
parent_dir_2 = strjoin(parent_dir_2, filesep);
output_dir_AR = strjoin([parent_dir_2, "Emulation-Data-Output\06_artifact_rejection"], filesep);
output_dir_AR_peaks = strjoin([output_dir_AR, 'peaks_-500_750'], filesep);
input_path = strjoin([output_dir_AR, 'peaks_-500_750'], filesep);
addpath([char(parent_dir) filesep 'functions']);

ft_defaults;

% load('fieldtrip_EEG_KJP_elec_61');
% load('fieldtrip_EEG_KJP_neighbours_61.mat');
% elec.label = upper(elec.label);


%% read data and change structure

% load([input_path filesep 'freq_all_theta']);
% long epochs
cd(output_dir_AR_peaks);
load(strjoin([neighbors_dir, 'fieldtrip_EEG_KJP_neighbours_61.mat'], filesep))



% convert to fieldtrip
for i = 1:length(ALLEEG2)
    eeg_field(i) = eeglab2fieldtrip(ALLEEG2(i), 'timelockanalysis', 'none');
end
% plot(eeg_field(26).time, eeg_field(26).avg')

zero_field = eeg_field; 
for i = 1:length(eeg_field)
    zero_field(i).avg = zeros(size(zero_field(i).avg));
end

for i = 1:length(eeg_field)
    eeg_cell{i} = eeg_field(i);
    zero_cell{i} = zero_field(i);
end

% I NEED THIS: ALTERNATIVELY, I CAN GENERATE IT MYSELF
% load('fieldtrip_EEG_KJP_elec_61');
% load('fieldtrip_EEG_KJP_neighbours_61.mat');
% this data is missing, Adriana will look it up late.r 
% elec.label = upper(elec.label);

% load(['C:\wilken\Emulation-Data-Input' filesep 'freq_all_theta']);
% 
% % was already here: 
% constant = cellfun(@(x) x.constant, freq_all, 'UniformOutput', false);
% rand1 = cellfun(@(x) x.rand1, freq_all, 'UniformOutput', false);
% rand2 = cellfun(@(x) x.rand2, freq_all, 'UniformOutput', false);
% % freq all is a cell 1x32, each cell contains subject data and it is split
% % in two cells, one for each condition. 
% erp_cond = mean_struct.all; % ??
% zero_cond = zeros(size(mean_struct.all)); % ??



%% init
% conds = {'constant', 'rand1', 'rand2'};
% conds: zero, erp signal at ~200
% cond_erp = erp_signal(:,'200ms') % pseudo code
% cond_zero = zeros(size(erp_signal, 1), size(erp_signal, 2)) % pseudo code

design = zeros(2, size(eeg_cell, 2)*2);
for i = 1:size(eeg_cell, 2)
    design(1,i) = i;
end
for i = 1:size(eeg_cell, 2)
    design(1,size(eeg_cell, 2)+i) = i;
end
design(2,1:size(eeg_cell, 2))        = 1;
design(2,size(eeg_cell, 2)+1:2*size(eeg_cell, 2)) = 2;


%% Get Neighbours

% generates a file called 'neighbours'
load(strjoin([neighbors_dir, 'fieldtrip_EEG_KJP_neighbours_61.mat'], filesep))

cfg_neighb        = [];
cfg_neighb.method = 'distance';
neighbours        = neighbours % ft_prepare_neighbours(cfg_neighb, dataFC_LP);


%%
% CLUSTER TEST
% erp vs zero
calpha  = 0.0009;
alpha   = 0.0009;
latency = [0.172, 0.228];

% cfg is the configuraiton structure of fieldtrip
cfg                     = [];
cfg.design              = design;
cfg.uvar                = 1;
cfg.ivar                = 2;
cfg.channel             = {'all'};
cfg.avgovertime         = 'no';
cfg.method              = 'montecarlo';
cfg.statistic           = 'depsamplesT'; % really indepsamples??? With comparison against 0 
cfg.correctm            = 'cluster';
cfg.clusteralpha        = calpha;               % 0.05;
cfg.clusterstatistic    = 'maxsum';
cfg.minnbchan           = 0;
cfg.neighbours          = neighbours;
cfg.tail                = 0;
cfg.clustertail         = 0;
cfg.alpha               = alpha;               % 0.025;
cfg.numrandomization    = 1000;
cfg.latency             = latency; % time interval over which the experimental
                                 % conditions must be compared (in seconds)


stats      = ft_timelockstatistics(cfg, eeg_cell{:}, zero_cell{:});


% Plot the CBP results

% ft_prepare_layout
load(strjoin([neighbors_dir, 'fieldtrip_EEG_KJP_layout_61.mat'], filesep))

cfg.layout = lay;
cfg.style = 'blank';
% cfg.contournum = 0;
cfg.highlightcolorpos         = [0 0 0.75];
cfg.highlightcolorneg         = [0.75 0 0];
cfg.highlightsymbolseries     = ['x', 'x', 'x', 'x', 'x'];
% cfg.saveaspng = "cluster_GO_19_pre.png";

figure(1)
ft_clusterplot(cfg, stats);
sgtitle(strjoin(["Significant clusters, calpha = ", calpha, " alpha = ", alpha], ""));

%%
% CLUSTER TEST
%constant vs. rand2
% calpha  = 0.05;
% alpha   = 0.05;
% 
% cfg                     = [];
% cfg.design              = design;
% cfg.uvar                = 1;
% cfg.ivar                = 2;
% cfg.channel             = {'all'};
% cfg.frequency           = [4 7];
% cfg.avgoverfreq         = 'yes';
% cfg.latency             = [min(constant{1}.time) max(constant{1}.time)]; %??
% cfg.avgovertime         = 'yes';
% cfg.method              = 'montecarlo';
% cfg.statistic           = 'depsamplesT';
% cfg.correctm            = 'cluster';
% cfg.clusteralpha        = calpha;               % 0.1;
% cfg.clusterstatistic    = 'maxsum';
% cfg.minnbchan           = 0;
% cfg.neighbours          = neighbours;
% cfg.tail                = -1;
% cfg.clustertail         = -1;
% cfg.alpha               = alpha;               % 0.025;
% cfg.numrandomization    = 1000;
% 
% stats      = ft_freqstatistics(cfg, constant{:}, rand2{:} );


%% Create ERP Plot of the data TODO


% TODO: sort by error size in epoch
% sort by
'[STUDY ALLEEG] = std_editset( STUDY, ALLEEG, 'name','emulation_study_old','updatedat','on','rmclust','off' );
     [STUDY ALLEEG] = std_checkset(STUDY, ALLEEG);'
'[STUDY ALLEEG] = std_precomp(STUDY, ALLEEG, {},'savetrials','on','interp','on','recompute','on','erp','on','erpparams',{'rmbase' [-500 0] },'erpim','on','erpimparams',{'nlines' 10 'smoothing' 10});['' ...
    ''][STUDY EEG] = pop_savestudy( STUDY, EEG, 'filename','study_peaks_epoched.study','filepath','C:\\wilken\\Emulation-Data-Output\\study\\')
std_erspplot(STUDY, ALLEEG)


%% Do TF Decomposition on data


%% Surface Laplacians CSD toolbox (aka current source density)

% keep default settings
addpath(genpath('C:\wilken\eeglab2019_0\CSDtoolbox'))

% tutorial: https://psychophysiology.cpmc.columbia.edu/Software/CSDtoolbox/tutorial.html
% common errors: https://psychophysiology.cpmc.columbia.edu/Software/CSDtoolbox/errors.html
% FAQ: https://psychophysiology.cpmc.columbia.edu/Software/CSDtoolbox/FAQ.html

% Get usable list of electrodes from EEGlab data structure
for site = 1:60 
    trodes{site}=(ALLEEG(1).chanlocs(site).labels);
end;
trodes=trodes';

% Get Montage for use with CSD Toolbox
Montage_64=ExtractMontage('10-5-System_Mastoids_EGI129.csd',trodes);
MapMontage(Montage_64);

%% Derive G and H!
[G,H] = GetGH(Montage_64);

%% Save G and H to later import when doing the CSD transform on files
% save('G:\PhysioData\MN_Fear\G.mat', 'G');
% save('G:\PhysioData\MN_Fear\H.mat', 'H');

% revised method to store G and H matrices with CSD montage for later import
Montage = Montage_64;                             % use generic variable name
save G:\PhysioData\CSDmontage_64.mat G H Montage; % save variables to Matlab file
clear G H Montage;                                % remove variables from workspace
load G:\PhysioData\CSDmontage_64.mat;             % restore variables to the workspace

%% SAMPLE TO OPEN BIOSEMI -- can use instead of neuroscan version above
% [FileName,PathName,FilterIndex] = uigetfile('L:\Physiodata\*.bdf','Choose File to Open')
% FullName = [PathName FileName];
% [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
% EEG = pop_biosig(FullName, 'ref',[65 66] ,'blockepoch','off'); % Choose ref sites per your montage
% [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'gui','off'); 
% eeglab redraw;

tic                                        % stopwatch on
for ne = 1:length(ALLEEG(1).epoch)               % loop through all epochs
    myEEG = single(ALLEEG(1).data(:,:,ne));      % reduce data precision to reduce memory demand
    MyResults = CSD(myEEG,G,H);            % compute CSD for <channels-by-samples> 2-D epoch
    data(:,:,ne) = MyResults;              % assign data output
end
looping_CSD_final = double(data);          % final CSD data
looping_time = toc                         % stopwatch off

data(:,:,:) = NaN;                         % re-initialize data output


%% Stuff I tried out. 

STUDY = std_erpplot(STUDY,ALLEEG)

% figure; 
% pop_plottopo(EEG, [1:60] , ''5STJS'', 0, ''ydir'',1);
% 
% pop_topoplot(EEG, 1, [-500 748] ,'5STJS',[1 2] ,0,'electrodes','on')
% 
% pop_erpplot()
% 
% 
% 
% https://www.fieldtriptoolbox.org/tutorial/cluster_permutation_timelock/


