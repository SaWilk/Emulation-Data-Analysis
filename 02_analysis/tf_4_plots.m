% Data Analysis Pursuit-Tracking and Pursuit-Occlusion Paradigm
% Emulation pilot study 2021

% This script contains: 
% TF and CBPT plots

% Adriana Böttcher
% 25.08.22

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
subjects_B      = subjects;
subjects_B(18)    = [];

% conditions
conds_A = {'const', 'rand1', 'rand2'};
conds_B = {'occl', 'nonoccl'};

% initialize input and output folder
inputpath_A     = 'R:\AG-Beste-Studien\Emulation\06_analysis\Emulation-Data-Output\09_segments_A_TFT';
inputpath_B     = 'R:\AG-Beste-Studien\Emulation\06_analysis\Emulation-Data-Output\09_segments_B_TFT';
inputpath_CBPT  = 'R:\AG-Beste-Studien\Emulation\06_analysis\Emulation-Data-Output\10_CBPT';

%% load data

% GAVs
load([inputpath_A filesep 'freq_GAV_A']);
load([inputpath_B filesep 'freq_GAV_B']);

% CBPT results
load([inputpath_CBPT filesep 'CBPT']);

%% lims
theta_lims   = [4 7];
alpha_lims   = [8 12];
all_freq     = [2 20];

timelim_A = [0 3];
timelim_B = [0 2];

%% calculate differences for plotting contrasts

freq_GAV_A.diff_const_rand1 = freq_GAV_A.const;
freq_GAV_A.diff_const_rand1.powspctrm = freq_GAV_A.const.powspctrm - freq_GAV_A.rand1.powspctrm;

freq_GAV_A.diff_const_rand2 = freq_GAV_A.const;
freq_GAV_A.diff_const_rand2.powspctrm = freq_GAV_A.const.powspctrm - freq_GAV_A.rand2.powspctrm;

freq_GAV_A.diff_rand1_rand2 = freq_GAV_A.rand1;
freq_GAV_A.diff_rand1_rand2.powspctrm = freq_GAV_A.rand1.powspctrm - freq_GAV_A.rand2.powspctrm;

freq_GAV_B.diff_occl_nonoccl = freq_GAV_B.nonoccl;
freq_GAV_B.diff_occl_nonoccl = freq_GAV_B.occl.powspctrm - freq_GAV_B.nonoccl.powspctrm;

%% plotting tf results task A

% multiplot constant

cfg = [];
cfg.layout = lay;
layout = ft_prepare_layout(cfg);

cfg.ylim = all_freq;
cfg.xlim = timelim_A;
cfg.zlim = 'maxmin';
cfg.layout = layout;
cfg.title = 'Task A: avg const';
ft_multiplotTFR(cfg, freq_GAV_A.const);
ft_topoplotTFR(cfg, freq_GAV_A.const);

% multiplot random 1

cfg.ylim = all_freq;
cfg.xlim = timelim_A;
cfg.zlim = 'maxmin';
cfg.layout = layout;
cfg.title = 'Task A: avg rand1';
ft_multiplotTFR(cfg, freq_GAV_A.rand1);
ft_topoplotTFR(cfg, freq_GAV_A.rand1);

% multiplot random 2

cfg.ylim = all_freq;
cfg.xlim = timelim_A;
cfg.zlim = 'maxmin';
cfg.layout = layout;
cfg.title = 'Task A: avg rand2';
ft_multiplotTFR(cfg, freq_GAV_A.rand2);

% multiplotcfg.ylim = all_freq;
cfg.xlim = timelim_A;
cfg.ylim = theta_lims;
cfg.zlim = 'maxmin';
cfg.layout = layout;
cfg.title = 'Task A: difference betweeen avg const and avg rand1';
ft_multiplotTFR(cfg, freq_GAV_A.diff_const_rand1); 
ft_topoplotTFR(cfg, freq_GAV_A.diff_const_rand1);
colorbar;

%% plotting stuff task B

% multiplot constant

cfg = [];
cfg.layout = lay;
layout = ft_prepare_layout(cfg);

cfg.ylim = all_freq;
cfg.xlim = timelim_B;
cfg.zlim = 'maxmin';
cfg.layout = layout;
cfg.title = 'Task B: avg occl';
ft_multiplotTFR(cfg, freq_GAV_B.occl);

% multiplot random 1

cfg.ylim = all_freq;
cfg.xlim = timelim_B;
cfg.zlim = 'maxmin';
cfg.layout = layout;
cfg.title = 'Task B: avg nonoccl';
ft_multiplotTFR(cfg, freq_GAV_B.nonoccl);

% constast/difference
cfg.title = 'Task B: avg occluded - non_occluded';
ft_multiplotTFR(cfg, freq_GAV_B.diff_occl_nonoccl);
ft_topoplotTFR(cfg, freq_GAV_B.diff_occl_nonoccl);
colorbar;

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