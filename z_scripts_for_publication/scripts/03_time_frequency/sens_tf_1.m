% The neurophysiology of continuous action monitoring
% Saskia Wilken, Adriana Böttcher, Nico Adelhöfer, Markus Raab, Sven
% Hoffmann, Christian Beste

% wavelet analysis on single trial level

% created by:
% Adriana Boettcher, Cognitive Neurophysiology TU Dresden
% 2022

%% 
clc;
clearvars;

% load fieldtrip toolbox

% load local EEG configuration (electrodes, neighbours, layout)

% initialize input and output folder

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
        cfg.pad         = 'nextpow2'; %padding; rounds the maximum trial length up to the next power of 2
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
        cfg.pad         = 'nextpow2'; %padding; rounds the maximum trial length up to the next power of 2
        cfg.output      = 'pow';
        cfg.keeptrials  = 'no';
        
        % perform TFT with configuration
        [freq] = ft_freqanalysis(cfg, data);
        
        % save results of TFT in outputpath
        outputname = [outputpath_B filesep char(sbj) '_' char(cond) '_freq'];
        save(outputname, 'freq')
    end 
end

