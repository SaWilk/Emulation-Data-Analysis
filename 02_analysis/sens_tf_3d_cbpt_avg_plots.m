% Data Analysis Pursuit-Tracking and Pursuit-Occlusion Paradigm
% Emulation pilot study 2021

% This script contains: 
% plots for CBPT ( A & B, A split for time intervals )

% 30.09.22

%% 
clc;
clearvars;

% add path for custom functions
addpath('R:\AG-Beste-Studien\Emulation\06_analysis\Emulation-Data-Analysis\functions');

% load fieldtrip toolbox
addpath('R:\AG-Beste-Orga\Skripts\Toolbox\fieldtrip-20210212\');
ft_defaults;

% load custom ft toolbox
addpath('R:\AG-Beste-Orga\Skripts\Fieldtrip\ft_mmcustom');

% load local EEG configuration (electrodes, neighbours, layout)
load fieldtrip_EEG_KJP_elec_61;
load fieldtrip_EEG_KJP_layout_61;
load fieldtrip_EEG_KJP_neighbours_61;

% load CBPT output
load('R:\AG-Beste-Studien\Emulation\06_analysis\Emulation-Data-Output\10_CBPT\CBPT_avg_exclude_start_all.mat');

outputpath = 'R:\AG-Beste-Studien\Emulation\06_analysis\Emulation-Data-Output\10_CBPT\plots';

%% 2 loops: task A and B
% loop and reduce data (downsampling) and plot downsampled data

CBPT_reduced = {};

for intrv = 2%:size(fields(CBPT_avg.A))
    intervals = fields(CBPT_avg.A);
    for frq = 1%:size(fields(CBPT_avg.A.(intervals{intrv})))

        freqs = fields(CBPT_avg.A.(intervals{intrv}));
        
        % set filename
        filename = ['CBPT_A_' intervals{intrv} '_' freqs{frq} '_const_rand1'];

        % extract time window for this test
        length = size(CBPT_avg.A.(intervals{intrv}).(freqs{frq}).const_rand1.time, 2);

        % reduce time to 1:5:length
        CBPT_reduced.A.(intervals{intrv}).(freqs{frq}).const_rand1 = CBPT_select_timepoints(CBPT_avg.A.(intervals{intrv}).(freqs{frq}).const_rand1, [1:10:length]);

        % and plot
        cfg.alpha                   = 0.001;
        cfg.parameter               = 'stat';
        cfg.highlightcolorpos       = [0 0 0];
        cfg.highlightcolorneg       = [1 0 0];
        cfg.highlightsymbolseries   = ['x', 'x', 'x', 'x', 'x'];
        cfg.highlightsizeseries     = [8,8,8,8,8];
        cfg.subplotsize             = [4 6];
         
        
        ft_clusterplot(cfg, CBPT_reduced.A.(intervals{intrv}).(freqs{frq}).const_rand1);

        %exportgraphics(gcf,[outputpath filesep filename '.png'],'Resolution',1000,'BackgroundColor','white');    
    end % freq
end % interval

%% task B

% extract time window for this test
length = size(CBPT_avg.B.theta.time, 2);

% reduce time to 1:5:length
CBPT_reduced.B.theta = CBPT_select_timepoints(CBPT_avg.B.theta, [1:10:length]);

% and plot
cfg.alpha                   = 0.001;
cfg.parameter               = 'stat';
cfg.highlightcolorpos       = [0 0 0];
cfg.highlightcolorneg       = [1 0 0];
cfg.highlightsymbolseries   = ['x', 'x', 'x', 'x', 'x'];
cfg.highlightsizeseries     = [8,8,8,8,8];
cfg.subplotsize             = [1 1];


ft_clusterplot(cfg, CBPT_reduced.B.theta);

% 
% %% paul plot code
% % cfg.style                   = 'blank';
% cfg.alpha                   = 0.001;
% cfg.parameter               = 'stat';
% cfg.highlightcolorpos       = [0 0 0];
% cfg.highlightcolorneg       = [1 0 0];
% cfg.highlightsymbolseries   = ['x', 'x', 'x', 'x', 'x'];
% cfg.highlightsizeseries     = [8,8,8,8,8];
% cfg.subplotsize             = [3 5];
% 
% ft_clusterplot(cfg, CBPT_avg.A.until_500.theta.const_rand1);
% 
% CBPT_theta_until500_constrand1_reduced = CBPT_select_timepoints(CBPT_avg.A.until_500.theta.const_rand1, [1:5:size(CBPT_avg.A.until_500.theta.const_rand1.time, 2)]);
% 
% cfg.alpha                   = 0.001;
% cfg.parameter               = 'stat';
% cfg.highlightcolorpos       = [0 0 0];
% cfg.highlightcolorneg       = [1 0 0];
% cfg.highlightsymbolseries   = ['x', 'x', 'x', 'x', 'x'];
% cfg.highlightsizeseries     = [8,8,8,8,8];
% cfg.subplotsize             = [4 5];
% 
% ft_clusterplot(cfg, CBPT_theta_until500_constrand1_reduced);
% 
% % %% Plot the CBP results Nicos Approach
% % % https://github.com/fieldtrip/fieldtrip/blob/master/ft_clusterplot.m
% % % https://www.fieldtriptoolbox.org/tutorial/plotting/
% % % https://github.com/fieldtrip/fieldtrip/blob/release/ft_topoplotER.m
% % 
% % cfg.layout = lay;
% % % cfg.style = 'blank';
% % % cfg.contournum = 0;
% % cfg.highlightcolorpos         = [0 0 0.75];
% % cfg.highlightcolorneg         = [0.75 0 0];
% % % cfg.highlightsymbolseries     = ['x', 'x', 'x', 'x', 'x'];
% % cfg.subplotsize    = [4, 4];
% % % cfg.saveaspng = "cluster_GO_19_pre.png";
% % 
% % figure()
% % ft_clusterplot(cfg, stats);
% % sgtitle(strjoin(["Significant clusters, calpha = ", calpha, " alpha = ", alpha], ""));
% % colorbar()
% % 
% % % topoplot verwenden um t-werte der channels zu plotten 
% % % 'mask' um nur überschwellige Werte anzuzeigen
% % % Bedingungsvergleiche testen
% % 
% % %% explore the results a little bit further
% % % 
% % % avg_diff_B = freq_all_B;
% % % avg_diff.powspctrm = avg_constant.powspctrm - avg_rand1.powspctrm;
% % % 
% % 
% % 
% % %% paul plot code
% % % cfg.style                   = 'blank';
% % cfg.alpha                   = alpha;
% % cfg.parameter               = 'stat';
% % cfg.highlightcolorpos       = [0 0 0];
% % cfg.highlightcolorneg       = [1 0 0];
% % cfg.highlightsymbolseries   = ['x', 'x', 'x', 'x', 'x'];
% % cfg.highlightsizeseries     = [8,8,8,8,8];
% % cfg.subplotsize             = [1 1];
% 
% ft_clusterplot(cfg, stats);