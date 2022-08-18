% Data Analysis Pursuit-Tracking and Pursuit-Occlusion Paradigm
% Emulation pilot study 2021

% This script contains: 
% TF analysis on single trial level
% for task A data: constant vs. random1/random2

% load segmented data, perform wavelet analysis and save average data for
% each subject to calculate grand average in next script

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
inputpath   = 'R:\AG-Beste-Studien\Emulation\06_analysis\Emulation-Data-Output\08_ft_segments_A';
outputpath  = 'R:\AG-Beste-Studien\Emulation\06_analysis\Emulation-Data-Output\09_segments_A_TFT';

% load IDs of included subjects
load subjects;

% conditions
conds = {'const', 'rand1', 'rand2'};

% change directory to load data
cd(inputpath);

%%
% loop through conditions, then loop through subjects

for cond = 1:size(conds,2)
    cond = char(conds(cond));
    
    for sbj = 1:size(subjects,2)
        sbj = char(subjects(sbj));
        filename = [char(sbj) '_' char(cond)];
        data = load(filename);
        
        % get to data struct
        data = data.(cond);
        
        % configuration for TFT
        cfg             = [];
        cfg.method      = 'wavelet';
        cfg.width       = 5;
        cfg.gwidth      = 3; % length of the used wavelets in standard deviations of the implicit Gaussian kernel
        cfg.foi         = 1:1:15; % frequencies of interest
        cfg.toi         = data.time{1}; % use time field of data struct
        cfg.trials      = 'all'; % use all trials
        cfg.pad         = 'nextpow2'; %padding; rounds the maximum trial length up to the next power of 2
        cfg.output      = 'pow';
        cfg.keeptrials  = 'no';
        
        % perform TFT with configuration
        [freq] = ft_freqanalysis(cfg, data);
        
        % save results of TFT in outputpath
        outputname = [outputpath filesep char(sbj) '_' char(cond) '_freq'];
        save(outputname, 'freq')
    end 
end
    


