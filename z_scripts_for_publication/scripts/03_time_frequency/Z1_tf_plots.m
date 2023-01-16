% The neurophysiology of continuous action monitoring
% Saskia Wilken, Adriana Böttcher, Nico Adelhöfer, Markus Raab, Sven
% Hoffmann, Christian Beste

% TF and CBPT plots for sig time intervals
% for contrast: const/rand1 in experiment 1 (separated intervals) 
% & occl/nonoccl in experiment 2

% created by:
% Adriana Boettcher, Cognitive Neurophysiology TU Dresden
% 2022

%%
clc;
clearvars;

% load fieldtrip toolbox

% load local EEG configuration (electrodes, neighbours, layout)

% load IDs of included subjects
load subjects;

% conditions
conds_A = {'const', 'rand1'};
conds_B = {'occl', 'nonoccl'};

% initialize input and output folder

%% lims
theta_lims_plot   = [3 10];
alpha_lims_plot   = [7 14];
beta_lims_plot    = [12 24];

timelim_A = [0 3];
timelim_B = [0 2];

calpha  = 0.001;
alpha   = 0.001;

%% load data
% GAVs
load([inputpath_TF_A filesep 'freq_GAV_A']);
load([inputpath_TF_B filesep 'freq_GAV_B']);

% CBPT results
load([inputpath_CBPT filesep 'CBPT']);
load([inputpath_CBPT filesep 'CBPT_avg_exclude_start_all']);

% change electrode labels to match neighbour and layout file
freq_GAV_A.const.label = upper(freq_GAV_A.const.label);
freq_GAV_A.rand1.label = upper(freq_GAV_A.rand1.label);
freq_GAV_A.rand2.label = upper(freq_GAV_A.rand2.label);

freq_GAV_B.occl.label = upper(freq_GAV_B.occl.label);
freq_GAV_B.nonoccl.label = upper(freq_GAV_B.nonoccl.label);

%% reduce trial length according to condition
% rand1: 500-3000
% const: 0-2500
% then, perform CBPT on whole interval (cfg.latency = 'all')

freq_GAV_A_after_500 = {};

timelims_rand1 = [0.5 3];
timelims_other = [0 2.5];

for cond = 1:size(conds_A, 2)
    % set timelims according to condition
    if strcmp(conds_A{cond}, 'rand1')
        timelims = timelims_rand1;
    else
        timelims = timelims_other;
    end

    % copy original structure
    freq_GAV_A_after_500.(conds_A{cond}) = freq_GAV_A.(conds_A{cond});
   
    % reduce time window
    % get start and end index
    start_ind = find(freq_GAV_A_after_500.(conds_A{cond}).time > timelims(1), 1, 'first');
    end_ind = find(freq_GAV_A_after_500.(conds_A{cond}).time > timelims(2), 1, 'first');
    % apply time window to time info and power spectrum
    freq_GAV_A_after_500.(conds_A{cond}).time      = freq_GAV_A_after_500.(conds_A{cond}).time(:, start_ind:end_ind);
    freq_GAV_A_after_500.(conds_A{cond}).powspctrm = freq_GAV_A_after_500.(conds_A{cond}).powspctrm(:, :, start_ind:end_ind);
end

% new dataframe for: 0-500 ms (same for all conditions)

freq_GAV_A_until_500 = {};
timelims = [0 0.5];

for cond = 1:size(conds_A, 2)
    % copy original structure
    freq_GAV_A_until_500.(conds_A{cond}) = freq_GAV_A.(conds_A{cond});
   
    % reduce time window
    % get start and end index
    start_ind = find(freq_GAV_A_until_500.(conds_A{cond}).time > timelims(1), 1, 'first');
    end_ind = find(freq_GAV_A_until_500.(conds_A{cond}).time > timelims(2), 1, 'first');
    % apply time window to time info and power spectrum
    freq_GAV_A_until_500.(conds_A{cond}).time      = freq_GAV_A_until_500.(conds_A{cond}).time(:, start_ind:end_ind);
    freq_GAV_A_until_500.(conds_A{cond}).powspctrm = freq_GAV_A_until_500.(conds_A{cond}).powspctrm(:, :, start_ind:end_ind);
end

%% calculate differences for plotting contrasts

freq_GAV_A_after_500.diff_const_rand1 = freq_GAV_A_after_500.const;
freq_GAV_A_after_500.diff_const_rand1.powspctrm = freq_GAV_A_after_500.const.powspctrm - freq_GAV_A_after_500.rand1.powspctrm;

freq_GAV_A_until_500.diff_const_rand1 = freq_GAV_A_until_500.const;
freq_GAV_A_until_500.diff_const_rand1.powspctrm = freq_GAV_A_until_500.const.powspctrm - freq_GAV_A_until_500.rand1.powspctrm;

freq_GAV_B.diff_occl_nonoccl = freq_GAV_B.nonoccl;
freq_GAV_B.diff_occl_nonoccl.powspctrm = freq_GAV_B.occl.powspctrm - freq_GAV_B.nonoccl.powspctrm;

%% singleplot & topoplot for significant electrodes

% black markers for positive clusters, red for negative clusters

% 1.) task A, until 500
%       select data
freq_data = freq_GAV_A_until_500.diff_const_rand1;

% 1.1.1) task A, until 500, alpha, const vs. rand1

%       get significant electrodes and CBPT output
CBPT_output = CBPT_avg.A.until_500.alpha.const_rand1;
[~, labels] = CBPT_get_sig_elec_from_timewindow(CBPT_output);

%       singleplot for sig electrodes
cfg         = [];
cfg.layout  = lay;
cfg.ylim    = alpha_lims_plot;
cfg.xlim    = timelims;
cfg.zlim    = 'maxmin';
cfg.channel = labels;
ft_singleplotTFR(cfg, freq_data); 
colorbar;

%       topoplot for all electrodes, highlighted sig electrodes
cfg                     = [];
cfg.layout              = lay;
cfg.ylim                = alpha_lims_plot;
cfg.xlim                = timelims;
cfg.zlim                = 'maxmin';
cfg.channel             = 'all'; 
cfg.marker              = 'off';
cfg.highlight           = 'on'; % highlight electrodes yes/no
cfg.highlightchannel    = labels;
cfg.highlightsymbol     = 'x';
cfg.highlightcolor      = [1 0 0]; % negative cluster = red
ft_topoplotTFR(cfg, freq_data)
colorbar;


% 1.1.2) task A, until 500, beta, const vs. rand1

%       get significant electrodes and CBPT output
CBPT_output = CBPT_avg.A.until_500.beta.const_rand1;
[~, labels] = CBPT_get_sig_elec_from_timewindow(CBPT_output);

%       singleplot for sig electrodes
cfg         = [];
cfg.layout  = lay;
cfg.ylim    = beta_lims_plot;
cfg.xlim    = timelims;
cfg.zlim    = 'maxmin';
cfg.channel = labels;
ft_singleplotTFR(cfg, freq_data); 
colorbar;

%       topoplot for all electrodes, highlighted sig electrodes
cfg                     = [];
cfg.layout              = lay;
cfg.ylim                = beta_lims_plot;
cfg.xlim                = timelims;
cfg.zlim                = 'maxmin';
cfg.channel             = 'all'; 
cfg.marker              = 'off';
cfg.highlight           = 'on'; 
cfg.highlightchannel    = labels;
cfg.highlightsymbol     = 'x';
cfg.highlightcolor      = [1 0 0]; % negative cluster = red
ft_topoplotTFR(cfg, freq_data)
colorbar;


% 1.1.3) task A, until 500, theta, const vs. rand1

%       get significant electrodes and CBPT output
CBPT_output = CBPT_avg.A.until_500.theta.const_rand1;
[~, labels] = CBPT_get_sig_elec_from_timewindow(CBPT_output);

%       singleplot for sig electrodes
cfg         = [];
cfg.layout  = lay;
cfg.ylim    = theta_lims_plot;
cfg.xlim    = timelims;
cfg.zlim    = 'maxmin';
cfg.channel = labels;
ft_singleplotTFR(cfg, freq_data); 
colorbar;

%       topoplot for all electrodes, highlighted sig electrodes
cfg                     = [];
cfg.layout              = lay;
cfg.ylim                = theta_lims_plot;
cfg.xlim                = timelims;
cfg.zlim                = 'maxmin';
cfg.channel             = 'all'; 
cfg.marker              = 'off';
cfg.highlight           = 'on';
cfg.highlightchannel    = labels;
cfg.highlightsymbol     = 'x';
cfg.highlightcolor      = [1 0 0]; % negative cluster = red
ft_topoplotTFR(cfg, freq_data)
colorbar;


% 1.2.1) task A, after 500, theta, const vs. rand1

%       select data
freq_data = freq_GAV_A_after_500.diff_const_rand1;

%       get significant electrodes and CBPT output
CBPT_output = CBPT_avg.A.after_500.theta.const_rand1;
[~, labels] = CBPT_get_sig_elec_from_timewindow(CBPT_output);

%       singleplot for sig electrodes
cfg         = [];
cfg.layout  = lay;
cfg.ylim    = theta_lims_plot;
cfg.xlim    = timelims_other;
cfg.zlim    = 'maxmin';
cfg.channel = labels;
ft_singleplotTFR(cfg, freq_data); 
colorbar;

%       topoplot for all electrodes, highlighted sig electrodes
cfg                     = [];
cfg.layout              = lay;
cfg.ylim                = theta_lims_plot;
cfg.xlim                = timelims_other;
cfg.zlim                = 'maxmin';
cfg.channel             = 'all'; 
cfg.marker              = 'off';
cfg.highlight           = 'on'; % highlight electrodes yes/no
cfg.highlightchannel    = labels;
cfg.highlightsymbol     = 'x';
cfg.highlightcolor      = [0 0 0]; % positive cluster = black
ft_topoplotTFR(cfg, freq_data)
colorbar;

% 1.2.1) task B,theta, occl vs visible

%       select data
freq_data = freq_GAV_B.diff_occl_nonoccl;

%       get significant electrodes and CBPT output
CBPT_output = CBPT_avg.B.theta;
[~, labels] = CBPT_get_sig_elec_from_timewindow(CBPT_output);

%       singleplot for sig electrodes
cfg         = [];
cfg.layout  = lay;
cfg.ylim    = theta_lims_plot;
cfg.xlim    = timelim_B;
cfg.zlim    = 'maxmin';
cfg.channel = labels;
ft_singleplotTFR(cfg, freq_data); 
colorbar;

%       topoplot for all electrodes, highlighted sig electrodes
cfg                     = [];
cfg.layout              = lay;
cfg.ylim                = theta_lims_plot;
cfg.xlim                = timelim_B;
cfg.zlim                = 'maxmin';
cfg.channel             = 'all';
cfg.marker              = 'off';
cfg.highlight           = 'on'; 
cfg.highlightchannel    = labels;
cfg.highlightsymbol     = 'x';
cfg.highlightcolor      = [1 0 0]; % negative cluster = red
ft_topoplotTFR(cfg, freq_data)
colorbar;