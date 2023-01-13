% The neurophysiology of continuous action monitoring
% Saskia Wilken, Adriana Böttcher, Nico Adelhöfer, Markus Raab, Sven
% Hoffmann, Christian Beste

% cluster-based permutation tests for alpha, beta, theta frequency bands
% for experiment 2

% created by:
% Adriana Boettcher, Cognitive Neurophysiology TU Dresden
% 2022

%% 
clc;
clearvars;

% load fieldtrip toolbox

% add path for custom functions

% load local EEG configuration (electrodes, neighbours, layout)

% load IDs of included subjects
load subjects;

% conditions
conds_B = {'occl', 'nonoccl'};

% initialize input and output folder

%% cbpt parameters

avg_freq    = 'yes';
avg_time    = 'no';

theta_lims      = [4 7];
alpha_lims      = [8 12];
beta_lims       = [13 30];

timelim_B = [0 2];

calpha  = 0.001;
alpha   = 0.001;

%% data export struct

CBPT_avg    = {};
CBPT_avg.B  = {};

%% design matrix B

design_B = zeros(2, size(subjects_B, 2)*2);
for i = 1:size(subjects_B, 2)
    design_B(1,i) = i;
end
for i = 1:size(subjects_B, 2)
    design_B(1,size(subjects_B, 2)+i) = i;
end
design_B(2,1:size(subjects_B, 2))        = 1;
design_B(2,size(subjects_B, 2)+1:2*size(subjects_B, 2)) = 2;

%% load data
load([inputpath_B filesep 'freq_all_sbj_B']);

%% configuration task B

% occlusion vs non-occlusion
% for time intervals, averaged over time

cfg                     = [];
cfg.design              = design_B;
cfg.uvar                = 1; % row of design matrix containing unit variables, i.e. number of subjects
cfg.ivar                = 2; % row of design matrix containing conditions
cfg.channel             = {'all'};
cfg.frequency           = theta_lims;
cfg.avgoverfreq         = 'yes';
cfg.avgovertime         = 'yes';
cfg.latency             = timelim_B; % begin end in s
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

stats               = ft_freqstatistics(cfg, freq_all_B.(conds_B{1}){:}, freq_all_B.(conds_B{2}){:});
CBPT_avg.B.theta    = stats;

[sig_times, first, last]    = CBPT_sig_timewindow(CBPT_avg.B.theta, alpha);
CBPT_avg.B.theta.sig_times  = sig_times;
CBPT_avg.B.theta.first_sig  = first;
CBPT_avg.B.theta.last_sig   = last;

cfg.frequency       = alpha_lims;
stats               = ft_freqstatistics(cfg, freq_all_B.(conds_B{1}){:}, freq_all_B.(conds_B{2}){:});
CBPT_avg.B.alpha    = stats;

[sig_times, first, last]    = CBPT_sig_timewindow(CBPT_avg.B.alpha, alpha);
CBPT_avg.B.alpha.sig_times  = sig_times;
CBPT_avg.B.alpha.first_sig  = first;
CBPT_avg.B.alpha.last_sig   = last;

cfg.frequency       = beta_lims;
stats               = ft_freqstatistics(cfg, freq_all_B.(conds_B{1}){:}, freq_all_B.(conds_B{2}){:});
CBPT_avg.B.beta    = stats;

[sig_times, first, last]    = CBPT_sig_timewindow(CBPT_avg.B.beta, alpha);
CBPT_avg.B.beta.sig_times  = sig_times;
CBPT_avg.B.beta.first_sig  = first;
CBPT_avg.B.beta.last_sig   = last;
