% Data Analysis Pursuit-Tracking and Pursuit-Occlusion Paradigm
% Emulation pilot study 2021

% This script contains: 
% CBPT for significant time window derived in exploratory CBPT (sens_tf_3c)
% * read significant time window from CBPT 
% * perform CBPT with avg over time & freq for significant time window

% Adriana Böttcher
% 19/10/22

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

% load IDs of included subjects
load subjects;
% remove subject with missing data for B
subjects(18)    = [];

% set paths
inputpath_A     = 'R:\AG-Beste-Studien\Emulation\06_analysis\Emulation-Data-Output\09_segments_A_TFT';
inputpath_B     = 'R:\AG-Beste-Studien\Emulation\06_analysis\Emulation-Data-Output\09_segments_B_TFT';

outputpath = 'R:\AG-Beste-Studien\Emulation\06_analysis\Emulation-Data-Output\10_CBPT';

%% configuration

% contrasts with conditions
contrasts_A = {
    {'const', 'rand1'}, ...
};

% freq bands
freq_bands      = {{[4 7] [8 12] [13 30]}, ...
    {[4 7]}};
freq_label      = {{'theta' 'alpha' 'beta'}, ...
    {'theta'}};

% CBP parameters
avg_freq        = 'yes'; % always avg over freq
avg_time        = 'yes';

calpha  = 0.001;
alpha   = 0.001;

% time limits are set to significant time windows in exploratory CBPT
interval_labels_A   = {'until_500' 'after_500'};

%% prepare CBPT

% data export structure
CBPT_sig    = {};
CBPT_sig.A  = {};
CBPT_sig.B  = {};

% create design matrix
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

% remove 18th subject from freq data task A (each cell)
conds = fieldnames(freq_all_A);
for cond = 1:length(fieldnames(freq_all_A))
    freq_all_A.(conds{cond})(18) = [];
end

%% create empty export structure 

for intrv = 1:size(interval_labels_A, 2)
    for frq = 1:size(freq_label, 2)
        CBPT_sig.A.(interval_labels_A{intrv}).(freq_label{frq}) = {};
    end
end

%% create reduced tf data, reduced to significant time window for every condition

for intrv = 1%:size(interval_labels_A, 2)
    for frq = 1:size(freq_label{intrv}, 2)
        for contr = 1:size(contrasts_A, 2)

            this_cond1 = contrasts_A{contr}{1};
            this_cond2 = contrasts_A{contr}{2};
            this_contr = [this_cond1 '_' this_cond2];
    

            % get time interval from CBPT
            this_CBPT_result  	= CBPT_avg.A.(interval_labels_A{intrv}).(freq_label{intrv}{frq}).(this_contr);
            time_interval       = [this_CBPT_result.first_sig, this_CBPT_result.last_sig];

            % now calculate CBPT for significant time window and generate
            % plot
            cfg                     = [];
            cfg.design              = design;
            cfg.uvar                = 1; % row of design matrix containing unit variables, i.e. number of subjects
            cfg.ivar                = 2; % row of design matrix containing conditions
            cfg.channel             = {'all'};
            cfg.frequency           = [freq_bands{intrv}{frq}];
            cfg.avgoverfreq         = avg_freq;
            cfg.avgovertime         = avg_time;
            cfg.latency             = time_interval;
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
    
            stats = ft_freqstatistics(cfg, freq_all_A.(this_cond1){:}, freq_all_A.(this_cond2){:});
            CBPT_sig.A.(interval_labels_A{intrv}).(freq_label{intrv}{frq}).([this_cond1 '_' this_cond2]) = stats;
            CBPT_sig.A.(interval_labels_A{intrv}).(freq_label{intrv}{frq}).([this_cond1 '_' this_cond2]).interval = time_interval;
        end % contrast loop
    end % freq loop
end % interval loop


%% second interval: after_500

% since the two time intervals do not match (const: 0-2500, rand: 500-3000
% ms), the TF data are reduced to the significant time window (applied to
% the respective condition
% the significant time window in the CBPT output is given in the time
% dimension of rand1 (starting at 0.5) --> to get the relevant time window
% for const: minus 0.5 s

freq_all_A_short = {};

for intrv = 2%:size(interval_labels_A, 2)
    for frq = 1:size(freq_label{intrv}, 2)
        for contr = 1:size(contrasts_A, 2)

            this_cond1 = contrasts_A{contr}{1};
            this_cond2 = contrasts_A{contr}{2};
            this_contr = [this_cond1 '_' this_cond2];

            % get time interval from CBPT
            this_CBPT_result  	= CBPT_avg.A.(interval_labels_A{intrv}).(freq_label{intrv}{frq}).(this_contr);
            time_interval       = [this_CBPT_result.first_sig, this_CBPT_result.last_sig];
            
            % apply limits of significant time window to tf data and save
            % in extra structure (freq_all_short)
            for cond = 1:size(contrasts_A{contr},2)
                conds_A = contrasts_A{contr};

                % change significant time window depending on condition
                if strcmp(conds_A(cond), 'const')
                    time_interval = time_interval - 0.5;
                else
                    time_interval = [this_CBPT_result.first_sig, this_CBPT_result.last_sig];
                end % if statement

                for sbj = 1:size(freq_all_A.(conds_A{cond}), 2)
                    % copy original structure
                    freq_all_A_short.(conds_A{cond}){sbj} = freq_all_A.(conds_A{cond}){sbj};
                   
                    % reduce time window
                    % get start and end index
                    start_ind = find(freq_all_A_short.(conds_A{cond}){sbj}.time > time_interval(1), 1, 'first');
                    end_ind = find(freq_all_A_short.(conds_A{cond}){sbj}.time > time_interval(2), 1, 'first');
                    % apply time window to time info and power spectrum
                    freq_all_A_short.(conds_A{cond}){sbj}.time      = freq_all_A_short.(conds_A{cond}){sbj}.time(:, start_ind:end_ind);
                    freq_all_A_short.(conds_A{cond}){sbj}.powspctrm = freq_all_A_short.(conds_A{cond}){sbj}.powspctrm(:, :, start_ind:end_ind);
                end %subject loop
            end % condition loop

            % now calculate CBPT for significant time window and generate
            % plot
            cfg                     = [];
            cfg.design              = design;
            cfg.uvar                = 1; % row of design matrix containing unit variables, i.e. number of subjects
            cfg.ivar                = 2; % row of design matrix containing conditions
            cfg.channel             = {'all'};
            cfg.frequency           = [freq_bands{intrv}{frq}];
            cfg.avgoverfreq         = avg_freq;
            cfg.avgovertime         = avg_time;
            cfg.latency             = 'all';
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
    
            stats = ft_freqstatistics(cfg, freq_all_A_short.(this_cond1){:}, freq_all_A_short.(this_cond2){:});
            CBPT_sig.A.(interval_labels_A{intrv}).(freq_label{intrv}{frq}).([this_cond1 '_' this_cond2]) = stats;
            CBPT_sig.A.(interval_labels_A{intrv}).(freq_label{intrv}{frq}).([this_cond1 '_' this_cond2]).interval = time_interval;

        end % contrast loop
    end % freq loop
end % interval loop

%% save results

save([outputpath filesep 'CBPT_avg_in_sig_timewindow'], 'CBPT_sig');