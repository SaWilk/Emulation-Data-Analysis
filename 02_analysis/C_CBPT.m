%% Cluster Based Permutation Test 

% Script by Adriana Böttcher
% Adapted by Saskia Wilken

clear
clc


%% add fieldtrip path
eeglab;
close;

% file location
file_path = matlab.desktop.editor.getActiveFilename;
file_path = fileparts(file_path);
addpath(file_path);
% TODO: Add Fieldtrip Path 
input_path = strjoin([output_dir_AR, 'peaks_-500_750'], filesep);
ft_defaults;

load('fieldtrip_EEG_KJP_elec_61');
load('fieldtrip_EEG_KJP_neighbours_61.mat');
elec.label = upper(elec.label);


%% read data and change structure

load([input_folder filesep 'freq_all_theta']);
% convert to fieldtrip
FIELD_EEG = eeglab2fieldtrip(EEG)
constant = cellfun(@(x) x.constant, freq_all, 'UniformOutput', false);
rand1 = cellfun(@(x) x.rand1, freq_all, 'UniformOutput', false);
rand2 = cellfun(@(x) x.rand2, freq_all, 'UniformOutput', false);

%% init
conds = {'constant', 'rand1', 'rand2'};

%create design matrix for within test
design = zeros(2, size(constant, 2)*2);
for i = 1:size(constant, 2)
    design(1,i) = i;
end
for i = 1:size(constant, 2)
    design(1,size(constant, 2)+i) = i;
end
design(2,1:size(constant, 2))        = 1;
design(2,size(constant, 2)+1:2*size(constant, 2)) = 2;


%%
% CLUSTER TEST
%constant vs. rand1
calpha  = 0.9;
alpha   = 0.9;

cfg                     = [];
cfg.design              = design;
cfg.uvar                = 1;
cfg.ivar                = 2;
cfg.channel             = {'all'};
cfg.frequency           = [4 7];
cfg.avgoverfreq         = 'no';
cfg.avgovertime         = 'no';
%cfg.latency             = [max(constant{1}.time)]; %??
cfg.method              = 'montecarlo';
cfg.statistic           = 'depsamplesT';
cfg.correctm            = 'cluster';
cfg.clusteralpha        = calpha;               % 0.1;
cfg.clusterstatistic    = 'maxsum';
cfg.minnbchan           = 0;
cfg.neighbours          = neighbours;
cfg.tail                = 0;
cfg.clustertail         = 0;
cfg.alpha               = alpha;               % 0.025;
cfg.numrandomization    = 1000;

stats      = ft_freqstatistics(cfg, constant{:}, rand1{:} );

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