% Data Analysis Pursuit-Tracking and Pursuit-Occlusion Paradigm
% Emulation pilot study 2021

% This script contains: 
% * load source data for task A
% * calculate ratio for all contrasts & freq bands per subject
% * calculate grand average for all contrasts
% * plot grand average for each contrast on standard brain

% Adriana Böttcher
% 12.09.2022

%% configuration
tic

clc;
clearvars;

% contrasts with conditions
contrasts_A = {
    {'const', 'rand1'}, ...
    {'const', 'rand2'}, ...
    %{'rand1', 'rand2'} ... % not significant 
    };

% freq bands
freq_bands      = {[4 7] [8 12] [13 30]};
freq_label      = {'theta' 'alpha' 'beta'};

%% path definition

% load fieldtrip toolbox
addpath('R:\AG-Beste-Orga\Skripts\Toolbox\fieldtrip-20210212\');
ft_defaults;

% load custom ft toolbox
addpath('R:\AG-Beste-Orga\Skripts\Fieldtrip\ft_mmcustom');

% load subject list and remove subject with missing data for task B 
load('subjects.mat');
subjects(18) = [];

% initialize input and output folder
% input: segmented and preprocessed time domain data 
datapath      = 'R:\AG-Beste-Studien\Emulation\06_analysis\Emulation-Data-Output\11_DICS';

%% load MRI

load standard_bem; % BEM headmodel
load standard_mri.mat; % structural MRI (colin27)
load fieldtrip_EEG_KJP_elec_61; % electrodes & positions
load fieldtrip_EEG_KJP_layout_61; % really needed?
load fieldtrip_EEG_KJP_neighbours_61; % really needed?

% prepare MRI
mri = ft_convert_units(mri,'cm');           % better make sure units match
mri = ft_volumereslice([],mri);             % reslice for good measure

% change electrode labels
elec.label = upper(elec.label);             % upper case letters

%%

source_ratio            = [];
source_ratio_all        = [];
source_ratio_all_avg    = [];

for frq = 1:size(freq_label, 2)

    % load data for respective frequency band
    data = load([datapath filesep 'source_data_all_' freq_label{frq}]);
    data = data.(['source_data_all_' freq_label{frq}]);
    source_ratio.(freq_label{frq}) = data;

    for contr = 1%:size(contrasts_A, 2)

        % save conds of this contrast for indexing
        cond1 = contrasts_A{contr}{1};
        cond2 = contrasts_A{contr}{2};

        source_ratio.(freq_label{frq}).([char(cond1) '_' char(cond2)]) = [];

        for subj = 1:size(subjects, 2)

            % calculate ratio for conds of this contrast and save to extra
            % df to calculate GAV in next step

            % copy data to save ratio
            source_ratio.(freq_label{frq}).([char(cond1) '_' char(cond2)]){subj} = data.([char(cond1) '_' char(cond2)]){subj};
            source_ratio.(freq_label{frq}).([char(cond1) '_' char(cond2)]){subj}.ratio = source_ratio.(freq_label{frq}).([char(cond1) '_' char(cond2)]){subj}.(cond1);
            
            % calculate ratio and save in copied structure ((x-y)./(x+y))
            source_ratio.(freq_label{frq}).([char(cond1) '_' char(cond2)]){subj}.ratio.avg.pow = ...
                (source_ratio.(freq_label{frq}).([char(cond1) '_' char(cond2)]){subj}.(cond1).avg.pow - ...
                source_ratio.(freq_label{frq}).([char(cond1) '_' char(cond2)]){subj}.(cond2).avg.pow) ./ ...
                (source_ratio.(freq_label{frq}).([char(cond1) '_' char(cond2)]){subj}.(cond1).avg.pow + ...
                source_ratio.(freq_label{frq}).([char(cond1) '_' char(cond2)]){subj}.(cond2).avg.pow);

        end % subject loop
    
        ratio_this_contr = cellfun(@(x) x.ratio, source_ratio.(freq_label{frq}).([char(cond1) '_' char(cond2)]), 'UniformOutput',false);

        % calculate GAV for this contrast
        cfg                     = [];
        cfg.parameter           = 'avg.pow';
        cfg.keepindividual      = 'no' ;
        
        ratio_this_contr_avg    = ft_sourcegrandaverage(cfg, ratio_this_contr{:});

        % create plot for this contrast

        % todo
        source_ratio_all.(freq_label{frq}).([char(cond1) '_' char(cond2)])      = ratio_this_contr;
        source_ratio_all_avg.(freq_label{frq}).([char(cond1) '_' char(cond2)])  = ratio_this_contr_avg;
    end % contrast loop
end % freq loop

save([datapath filesep 'source_ratio_all_avg'], "source_ratio_all_avg");
