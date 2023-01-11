% Data Analysis Pursuit-Tracking and Pursuit-Occlusion Paradigm
% Emulation pilot study 2021

% This script contains: 
% TF analysis on single trial level
% for task A data: constant vs. random1/random2
% for task B data: occl vs. non-occl

% load segmented data, perform wavelet analysis and save average data for
% each subject to calculate grand average in next script

% adjusted for plots 26/10/22: changed padding to 5 for accurate frequency bins; changed foi to
% 1:0.5_35 for higher resolution

% Adriana Böttcher
% 17.08.22

%% 
clc;
clearvars;

% load fieldtrip toolbox
addpath('R:\AG-Beste-Orga\Skripts\Toolbox\fieldtrip-20210212\');
ft_defaults;

% load custom ft toolbox
addpath('R:\AG-Beste-Orga\Skripts\Fieldtrip\ft_mmcustom');

% load local EEG configuration (electrodes, neighbours, layout)
load fieldtrip_EEG_KJP_elec_61;
load fieldtrip_EEG_KJP_layout_61;
load fieldtrip_EEG_KJP_neighbours_61;

% initialize input and output folder
inputpath_A     = 'R:\AG-Beste-Studien\Emulation\06_analysis\Emulation-Data-Output\08_ft_segments_A';
inputpath_B     = 'R:\AG-Beste-Studien\Emulation\06_analysis\Emulation-Data-Output\08_ft_segments_B';

outputpath_A     = 'R:\AG-Beste-Studien\Emulation\06_analysis\Emulation-Data-Output\09_B_segments_A_TFT_adjusted_for_figures';
outputpath_B     = 'R:\AG-Beste-Studien\Emulation\06_analysis\Emulation-Data-Output\09_B_segments_B_TFT_adjusted_for_figures';

% load IDs of included subjects
load subjects;

% conditions
conds_A = {'const', 'rand1', 'rand2'};
conds_B = {'occl', 'nonoccl'};

% variables for TFT
gwidth  = 3;
width   = 5;
foi     = 1:1:35;

%% task A
% loop through conditions, then loop through subjects

% change directory to load data
cd(inputpath_A);

for cond = 1:size(conds_A,2)
    cond = char(conds_A(cond));
    
    for sbj = 1:size(subjects,2)
        sbj = char(subjects(sbj));
        filename = [char(sbj) '_' char(cond)];
        data = load(filename);
        
        % get to data struct
        data = data.(cond);
        
        % configuration for TFT
        cfg             = [];
        cfg.method      = 'wavelet';
        cfg.width       = width;
        cfg.gwidth      = gwidth; % length of the used wavelets in standard deviations of the implicit Gaussian kernel
        cfg.foi         = foi; % frequencies of interest
        cfg.toi         = data.time{1}; % use time field of data struct
        cfg.trials      = 'all'; % use all trials
        cfg.pad         = 5; %adjusted 26/10/22 due to incoherent freq bins %'nextpow2'; %padding; rounds the maximum trial length up to the next power of 2
        cfg.output      = 'pow';
        cfg.keeptrials  = 'no';
        
        % perform TFT with configuration
        [freq] = ft_freqanalysis(cfg, data);
        
        % save results of TFT in outputpath
        outputname = [outputpath_A filesep char(sbj) '_' char(cond) '_freq'];
        save(outputname, 'freq')
    end 
end
    
%% task B
% loop through conditions, then loop through subjects

% todo: remove:
% just for now: remove sbj 18 (KMY6K) since data could not be transferred
% to ft
subjects(18) = [];

% change directory to load data
cd(inputpath_B);

for cond = 1:size(conds_B,2)
    cond = char(conds_B(cond));
    
    for sbj = 1:size(subjects,2)
        sbj = char(subjects(sbj));
        filename = [char(sbj) '_' char(cond)];
        data = load(filename);
        
        % get to data struct
        data = data.(cond);
        
        % configuration for TFT
        cfg             = [];
        cfg.method      = 'wavelet';
        cfg.width       = width;
        cfg.gwidth      = gwidth; % length of the used wavelets in standard deviations of the implicit Gaussian kernel
        cfg.foi         = foi; % frequencies of interest
        cfg.toi         = data.time{1}; % use time field of data struct
        cfg.trials      = 'all'; % use all trials
        cfg.pad         = 5; %adjusted 26/10/22 due to incoherent freq bins 'nextpow2'; %padding; rounds the maximum trial length up to the next power of 2
        cfg.output      = 'pow';
        cfg.keeptrials  = 'no';
        
        % perform TFT with configuration
        [freq] = ft_freqanalysis(cfg, data);
        
        % save results of TFT in outputpath
        outputname = [outputpath_B filesep char(sbj) '_' char(cond) '_freq'];
        save(outputname, 'freq')
    end 
end

