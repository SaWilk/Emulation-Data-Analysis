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
%list all *.set files in inputpath
file_names = dir('*.set');
%concatenate into one cell array
files2read = {file_names.name};


% convert to fieldtrip
hdr = ft_read_header( files2read{1} );
data = ft_read_data( files2read{1}, 'header', hdr );
events = ft_read_event( files2read{1}, 'header', hdr );


% I NEED THIS: ALTERNATIVELY, I CAN GENERATE IT MYSELF
% load('fieldtrip_EEG_KJP_elec_61');
% load('fieldtrip_EEG_KJP_neighbours_61.mat');
% this data is missing, Adriana will look it up late.r 
% elec.label = upper(elec.label);

% was already here: 
constant = cellfun(@(x) x.constant, freq_all, 'UniformOutput', false);
rand1 = cellfun(@(x) x.rand1, freq_all, 'UniformOutput', false);
rand2 = cellfun(@(x) x.rand2, freq_all, 'UniformOutput', false);
% freq all is a cell 1x32, each cell contains subject data and it is split
% in two cells, one for each condition. 

%% init
% conds = {'constant', 'rand1', 'rand2'};
% conds: zero, erp signal at ~200
cond_erp = erp_signal(:,'200ms') % pseudo code
cond_zero = zeros(size(erp_signal, 1), size(erp_signal, 2)) % pseudo code

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

constant = rand(1, 100)
random = rand(1, 100)


%% Get Neighbours

cfg_neighb        = [];
cfg_neighb.method = 'distance';
neighbours        = ft_prepare_neighbours(cfg_neighb, dataFC_LP);


%%
% CLUSTER TEST
%constant vs. rand1
calpha  = 0.9;
alpha   = 0.9;

% cfg is the configuraiton structure of fieldtrip
cfg                     = [];
cfg.design              = [ones(1,50), ones(1,50)*2]; % design matrix, pseudo code
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
cfg.statistic           = 'indepsamplesT';

stats      = ft_timelockstatistics(cfg, constant{:}, rand1{:});


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