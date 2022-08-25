% Data Analysis Pursuit-Tracking and Pursuit-Occlusion Paradigm
% Emulation pilot study 2021

% This script contains: 
% CBPT for constant vs. random trajectory, task A
% for freqs: alpha, beta, theta

% Adriana Böttcher
% 23.08.22


% parameters and design matrix
% design matrix: 1s and 2s for contrast
% dependent samples t test

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

% load IDs of included subjects
load subjects;

%remove subject with missing data for B
subjects(18) = [];

% conditions
conds_A = {'const', 'rand1', 'rand2'};
conds_B = {'occl', 'nonoccl'};

% initialize input and output folder
inputpath_B     = 'R:\AG-Beste-Studien\Emulation\06_analysis\Emulation-Data-Output\09_segments_B_TFT';

%% create design matrix
% within design: row 1 contains subject numbers 2 times
% row 2 contains 1s and 2s for all subject, marking the within contrast

design = zeros(2, size(subjects, 2)*2);
for i = 1:size(subjects, 2)
    design(1,i) = i;
end
for i = 1:size(subjects, 2)
    design(1,size(subjects, 2)+i) = i;
end
design(2,1:size(subjects, 2))        = 1;
design(2,size(subjects, 2)+1:2*size(subjects, 2)) = 2;

%% load data
load([inputpath_B filesep 'freq_all_sbj_B']);

%% configuration

% occlusion vs non-occlusion
% for time intervals, averaged over time
calpha  = 0.001;
alpha   = 0.001;

cfg                     = [];
cfg.design              = design;
cfg.uvar                = 1; % row of design matrix containing unit variables, i.e. number of subjects
cfg.ivar                = 2; % row of design matrix containing conditions
cfg.channel             = {'all'};
cfg.frequency           = [4 7];
cfg.avgoverfreq         = ['no'];
cfg.avgovertime         = 'no';
cfg.latency             = [0 2]; % begin end in s
cfg.method              = 'montecarlo';
cfg.statistic           = 'depsamplesT';
cfg.correctm            = 'cluster';
cfg.clusteralpha        = calpha;               
cfg.clusterstatistic    = 'maxsum';
cfg.minnbchan           = 0;
cfg.neighbours          = neighbours;
cfg.tail                = 0;
cfg.clustertail         = 0;
cfg.alpha               = alpha;               
cfg.numrandomization    = 1000;


%% perform CBPT

stats = ft_freqstatistics(cfg, freq_all_B.(conds_B{1}){:}, freq_all_B.(conds_B{2}){:});


% The field prob contains the proportion of draws from the permutation 
% distribution with a maximum cluster-level statistic that is larger than 
% the cluster-level test statistic
stats.prob
stats.stat %  cluster-level test statistic (here with maxsum: the sum of 
% the T-values in this cluster).

%% Plot the CBP results Nicos Approach
% https://github.com/fieldtrip/fieldtrip/blob/master/ft_clusterplot.m
% https://www.fieldtriptoolbox.org/tutorial/plotting/
% https://github.com/fieldtrip/fieldtrip/blob/release/ft_topoplotER.m

cfg.layout = lay;
% cfg.style = 'blank';
% cfg.contournum = 0;
cfg.highlightcolorpos         = [0 0 0.75];
cfg.highlightcolorneg         = [0.75 0 0];
% cfg.highlightsymbolseries     = ['x', 'x', 'x', 'x', 'x'];
cfg.subplotsize    = [4, 4];
% cfg.saveaspng = "cluster_GO_19_pre.png";

figure()
ft_clusterplot(cfg, stats);
sgtitle(strjoin(["Significant clusters, calpha = ", calpha, " alpha = ", alpha], ""));
colorbar()

% topoplot verwenden um t-werte der channels zu plotten 
% 'mask' um nur überschwellige Werte anzuzeigen
% Bedingungsvergleiche testen

%% explore the results a little bit further
% 
% avg_diff_B = freq_all_B;
% avg_diff.powspctrm = avg_constant.powspctrm - avg_rand1.powspctrm;
% 


%% paul plot code
% cfg.style                   = 'blank';
cfg.alpha                   = alpha;
cfg.parameter               = 'stat';
cfg.highlightcolorpos       = [0 0 0];
cfg.highlightcolorneg       = [1 0 0];
cfg.highlightsymbolseries   = ['x', 'x', 'x', 'x', 'x'];
cfg.highlightsizeseries     = [8,8,8,8,8];
cfg.subplotsize             = [1 1];

ft_clusterplot(cfg, stats);