% Data Analysis Pursuit-Tracking and Pursuit-Occlusion Paradigm
% Emulation pilot study 2021

% This script contains: 
% CBPT_avg for constant vs. random trajectory, task A
% for freqs: alpha, beta, theta

% new:
% * now with loop over freqs and contrasts
% * exclude subject 18 in task A and B (KMY6K, missing data for B)
% * for two time intervals separately: 
%   * 0-500 ms, averaged for time & freq
%   * 500 ms onwards, averaged only for freq

% Adriana Böttcher
% 14/09/22

% update 16/09/22: change time windows: 
%   1) compare 500-3000 ms rand1 to 0-2500 ms rand2 & const
%   2) compare 0-500 ms rand1 to 0-500 ms rand2 & const

% update 20/09/2022: save as _exclude_start

%% configuration

clc;
clearvars;

% contrasts with conditions
contrasts_A = {
    {'const', 'rand1'}, ...
    {'const', 'rand2'}, ...
    {'rand1', 'rand2'} ... 
};

% freq bands
freq_bands      = {[4 7] [8 12] [13 30]};
freq_label      = {'theta' 'alpha' 'beta'};

% CBP parameters
avg_freq        = 'yes'; % always avg over freq
avg_time_A      = 'no';
avg_time_B      = 'yes'; % avg over time in B, A not

calpha  = 0.001;
alpha   = 0.001;

% two time limits for A (0-500 ms, 500 ms onwards)
interval_labels_A   = {'until_500' 'after_500'};
interval_lims_A     = {[0 0.5] 'all'};

% time limits for B
timelim_B = [0 2];

%% path definition

% load fieldtrip toolbox
addpath('R:\AG-Beste-Orga\Skripts\Toolbox\fieldtrip-20210212\');
ft_defaults;

% load custom ft toolbox
addpath('R:\AG-Beste-Orga\Skripts\Fieldtrip\ft_mmcustom');

% add path for custom functions
addpath('R:\AG-Beste-Studien\Emulation\06_analysis\Emulation-Data-Analysis\functions');

% load local EEG configuration (electrodes, neighbours, layout)
load fieldtrip_EEG_KJP_elec_61;
load fieldtrip_EEG_KJP_layout_61;
load fieldtrip_EEG_KJP_neighbours_61;

% load IDs of included subjects
load subjects;
%remove subject with missing data for B
subjects(18)    = [];

% initialize input and output folder
inputpath_A     = 'R:\AG-Beste-Studien\Emulation\06_analysis\Emulation-Data-Output\09_segments_A_TFT';
inputpath_B     = 'R:\AG-Beste-Studien\Emulation\06_analysis\Emulation-Data-Output\09_segments_B_TFT';

outputpath      = 'R:\AG-Beste-Studien\Emulation\06_analysis\Emulation-Data-Output\10_CBPT';
%output_suffix   = '_without_0-500_ms'; % old suffix: exlude 0-500 for all conditions
output_suffix   = '_exclude_start'; % new suffix: exclude 0-500 only for rand1

%% data export structure

CBPT_avg    = {};
CBPT_avg.A  = {};
CBPT_avg.B  = {};

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

load([inputpath_A filesep 'freq_all_sbj_A']);
load([inputpath_B filesep 'freq_all_sbj_B']);

% remove 18th subject from freq data task A (each cell)
conds = fieldnames(freq_all_A);
for cond = 1:length(fieldnames(freq_all_A))
    freq_all_A.(conds{cond})(18) = [];
end

%% create empty export structure 

for intrv = 1:size(interval_labels_A, 2)
    for frq = 1:size(freq_label, 2)
        CBPT_avg.A.(interval_labels_A{intrv}).(freq_label{frq}) = {};
    end
end

%% loop over freqbands and contrasts and perform CBPT for A
% first interval: before reducing the trial duration
for intrv = 1%:size(interval_labels_A, 2)
    for frq = 1:size(freq_label, 2)
        for contr = 1:size(contrasts_A, 2)
            this_cond1 = contrasts_A{contr}{1};
            this_cond2 = contrasts_A{contr}{2};

            this_contr = [this_cond1 '_' this_cond2];
    
            cfg                     = [];
            cfg.design              = design;
            cfg.uvar                = 1; % row of design matrix containing unit variables, i.e. number of subjects
            cfg.ivar                = 2; % row of design matrix containing conditions
            cfg.channel             = {'all'};
            cfg.frequency           = [freq_bands{frq}];
            cfg.avgoverfreq         = avg_freq;
            cfg.avgovertime         = avg_time_A;
            cfg.latency             = interval_lims_A{intrv}; % time interval was set before
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
    
            stats                           = ft_freqstatistics(cfg, freq_all_A.(contrasts_A{contr}{1}){:}, freq_all_A.(contrasts_A{contr}{2}){:});
            CBPT_avg.A.(interval_labels_A{intrv}).(freq_label{frq}).([contrasts_A{contr}{1} '_' contrasts_A{contr}{2}]) = stats;
    
            [sig_times, first, last] = CBPT_sig_timewindow(CBPT_avg.A.(interval_labels_A{intrv}).(freq_label{frq}).([contrasts_A{contr}{1} '_' contrasts_A{contr}{2}]), alpha);
            CBPT_avg.A.(interval_labels_A{intrv}).(freq_label{frq}).([contrasts_A{contr}{1} '_' contrasts_A{contr}{2}]).sig_times  = sig_times;
            CBPT_avg.A.(interval_labels_A{intrv}).(freq_label{frq}).([contrasts_A{contr}{1} '_' contrasts_A{contr}{2}]).first_sig  = first;
            CBPT_avg.A.(interval_labels_A{intrv}).(freq_label{frq}).([contrasts_A{contr}{1} '_' contrasts_A{contr}{2}]).last_sig   = last;
    
        end % contrast loop
    end % freq loop
end % interval loop

%% reduce trial length according to condition
% rand1: 500-3000
% rand2/const: 0-2500
% then, perform CBPT on whole interval (cfg.latency = 'all')

freq_all_A_short_short = {};
conds_A = {'rand1' 'const' 'rand2'};

timelims_rand1 = [0.5 3];
timelims_other = [0 2.5];

for cond = 1:size(conds_A, 2)

    % set timelims according to condition
    if strcmp(conds_A{cond}, 'rand1')
        timelims = timelims_rand1;
    else
        timelims = timelims_other;
    end

    for sbj = 1:size(freq_all_A.(conds_A{cond}), 2)
        % copy original structure
        freq_all_A_short.(conds_A{cond}){sbj} = freq_all_A.(conds_A{cond}){sbj};
       
        % reduce time window
        % get start and end index
        start_ind = find(freq_all_A_short.(conds_A{cond}){sbj}.time > timelims(1), 1, 'first');
        end_ind = find(freq_all_A_short.(conds_A{cond}){sbj}.time > timelims(2), 1, 'first');
        % apply time window to time info and power spectrum
        freq_all_A_short.(conds_A{cond}){sbj}.time      = freq_all_A_short.(conds_A{cond}){sbj}.time(:, start_ind:end_ind);
        freq_all_A_short.(conds_A{cond}){sbj}.powspctrm = freq_all_A_short.(conds_A{cond}){sbj}.powspctrm(:, :, start_ind:end_ind);
    end
end

%% loop over freqbands and contrasts and perform CBPT for A
% second interval: after reducing the trial duration

for intrv = 2%1:size(interval_labels_A, 2)
    for frq = 1:size(freq_label, 2)
        for contr = 1:size(contrasts_A, 2)
            this_cond1 = contrasts_A{contr}{1};
            this_cond2 = contrasts_A{contr}{2};

            this_contr = [this_cond1 '_' this_cond2];
    
            cfg                     = [];
            cfg.design              = design;
            cfg.uvar                = 1; % row of design matrix containing unit variables, i.e. number of subjects
            cfg.ivar                = 2; % row of design matrix containing conditions
            cfg.channel             = {'all'};
            cfg.frequency           = [freq_bands{frq}];
            cfg.avgoverfreq         = avg_freq;
            cfg.avgovertime         = avg_time_A;
            cfg.latency             = interval_lims_A{intrv}; % time interval was set before
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
    
            stats                           = ft_freqstatistics(cfg, freq_all_A_short.(contrasts_A{contr}{1}){:}, freq_all_A_short.(contrasts_A{contr}{2}){:});
            CBPT_avg.A.(interval_labels_A{intrv}).(freq_label{frq}).([contrasts_A{contr}{1} '_' contrasts_A{contr}{2}]) = stats;
    
            [sig_times, first, last] = CBPT_sig_timewindow(CBPT_avg.A.(interval_labels_A{intrv}).(freq_label{frq}).([contrasts_A{contr}{1} '_' contrasts_A{contr}{2}]), alpha);
            CBPT_avg.A.(interval_labels_A{intrv}).(freq_label{frq}).([contrasts_A{contr}{1} '_' contrasts_A{contr}{2}]).sig_times  = sig_times;
            CBPT_avg.A.(interval_labels_A{intrv}).(freq_label{frq}).([contrasts_A{contr}{1} '_' contrasts_A{contr}{2}]).first_sig  = first;
            CBPT_avg.A.(interval_labels_A{intrv}).(freq_label{frq}).([contrasts_A{contr}{1} '_' contrasts_A{contr}{2}]).last_sig   = last;
    
        end % contrast loop
    end % freq loop
end % interval loop


%% save structure

save([outputpath filesep 'CBPT_avg' output_suffix], 'CBPT_avg');