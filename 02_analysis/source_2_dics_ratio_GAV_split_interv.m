% Data Analysis Pursuit-Tracking and Pursuit-Occlusion Paradigm
% Emulation pilot study 2021

% This script contains: 
% * load source data for task A
% * calculate ratio for all contrasts & freq bands per subject
% * calculate grand average for all contrasts
% * interpolate data
% * plot grand average for each contrast on standard brain

% split for the two intervals: 0-500 ms; 500 ms onwards

% Adriana Böttcher
% 22.09.2022

%% configuration
tic

clc;
clearvars;

% contrasts with conditions
contrasts_A = {
    {'const', 'rand1'}, ...
    %{'const', 'rand2'}, ...
    %{'rand1', 'rand2'} ... % not significant 
    };

% freq bands, now separetely for time intervals 
% (0-500 ms: sign. differences in all freqbands; 500 onwards only theta)
freq_bands      = {
    {[4 7] [8 12] [13 30]}, ... % first interval, all freqs
    {[4 7]} ... % second interval, only theta
    };
freq_label      = {
    {'theta' 'alpha' 'beta'}, ...
    {'theta'}, ...
    };

% two time limits for A (0-500 ms, 500 ms onwards)
interval_labels_A = {'until_500' 'after_500'};

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
plots_path    = 'R:\AG-Beste-Studien\Emulation\06_analysis\Emulation-Data-Output\11_DICS\plots';

output_suffix = '_split_interval';

%% load MRI for plotting

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

%% calculate ratio per subject and GAV for each contrast

source_ratio            = [];
source_ratio_all        = [];
source_ratio_all_avg    = [];

for intrv = 1:size(interval_labels_A, 2)
    for frq = 1:size(freq_label{intrv}, 2)
    
        % load data for respective frequency band
        data = load([datapath filesep 'source_data_all_' freq_label{intrv}{frq} '_split_interv']);
        data = data.(['source_data_all_' freq_label{intrv}{frq}]);
        source_ratio.((freq_label{intrv}{frq}))= data;
    
        for contr = 1:size(contrasts_A, 2)
    
            % save conds of this contrast for indexing
            cond1 = contrasts_A{contr}{1};
            cond2 = contrasts_A{contr}{2};
    
            source_ratio.(freq_label{intrv}{frq}).([char(cond1) '_' char(cond2)]).(interval_labels_A{intrv}) = [];
    
            for subj = 1:size(subjects, 2)
    
                % calculate ratio for conds of this contrast and save to extra
                % df to calculate GAV in next step
    
                % copy data to save ratio
                source_ratio.(freq_label{intrv}{frq}).([char(cond1) '_' char(cond2)]).(interval_labels_A{intrv}){subj} = data.([char(cond1) '_' char(cond2)]).(interval_labels_A{intrv}){subj};
                source_ratio.(freq_label{intrv}{frq}).([char(cond1) '_' char(cond2)]).(interval_labels_A{intrv}){subj}.ratio = source_ratio.(freq_label{intrv}{frq}).([char(cond1) '_' char(cond2)]).(interval_labels_A{intrv}){subj}.(cond1);
                
                % calculate ratio and save in copied structure ((x-y)./(x+y))
                source_ratio.(freq_label{intrv}{frq}).([char(cond1) '_' char(cond2)]).(interval_labels_A{intrv}){subj}.ratio.avg.pow = ...
                    (source_ratio.(freq_label{intrv}{frq}).([char(cond1) '_' char(cond2)]).(interval_labels_A{intrv}){subj}.(cond1).avg.pow - ...
                    source_ratio.(freq_label{intrv}{frq}).([char(cond1) '_' char(cond2)]).(interval_labels_A{intrv}){subj}.(cond2).avg.pow) ./ ...
                    (source_ratio.(freq_label{intrv}{frq}).([char(cond1) '_' char(cond2)]).(interval_labels_A{intrv}){subj}.(cond1).avg.pow + ...
                    source_ratio.(freq_label{intrv}{frq}).([char(cond1) '_' char(cond2)]).(interval_labels_A{intrv}){subj}.(cond2).avg.pow);
    
            end % subject loop
        
            ratio_this_contr = cellfun(@(x) x.ratio, source_ratio.(freq_label{intrv}{frq}).([char(cond1) '_' char(cond2)]).(interval_labels_A{intrv}), 'UniformOutput',false);
    
            % calculate GAV for this contrast
            cfg                     = [];
            cfg.parameter           = 'avg.pow';
            cfg.keepindividual      = 'no' ;
            
            ratio_this_contr_avg    = ft_sourcegrandaverage(cfg, ratio_this_contr{:});
    
            % create plot for this contrast
    
            % todo
            source_ratio_all.(freq_label{intrv}{frq}).([char(cond1) '_' char(cond2)]).(interval_labels_A{intrv})     = ratio_this_contr;
            source_ratio_all_avg.(freq_label{intrv}{frq}).([char(cond1) '_' char(cond2)]).(interval_labels_A{intrv}) = ratio_this_contr_avg;
        end % contrast loop
    end % freq loop
end % interval loop

save([datapath filesep 'source_ratio_all_avg' output_suffix], "source_ratio_all_avg", '-v7.3');

%% interpolate data on MRI and save

%load averaged source data
%load([datapath filesep 'source_ratio_all_avg' output_suffix]);

source_ratio_avg_int = [];

% create new data frame with interpolated data for each freq and contrast
for intrv = 1:size(interval_labels_A, 2)
    for frq = 1:size(freq_label{intrv}, 2)
        for contr = 1:size(contrasts_A, 2)
    
            this_contrast = [contrasts_A{contr}{1} '_' contrasts_A{contr}{2}];
    
            cfg                         = [];
            cfg.parameter               = 'pow';
            source_ratio_avg_int.(freq_label{intrv}{frq}).(this_contrast).(interval_labels_A{intrv}) = ft_sourceinterpolate(cfg, source_ratio_all_avg.(freq_label{intrv}{frq}).(this_contrast).(interval_labels_A{intrv}), mri);
        end % contrast loop
    end % freq loop
end % interval loop

% save interpolated data
save([datapath filesep 'source_ratio_avg_int' output_suffix], "source_ratio_avg_int", '-v7.3');

%% plotting

cfg                 = [];
cfg.funparameter    = 'pow';
cfg.maskparameter   = 'mask';
%cfg.maskparameter   = cfg.funparameter;
cfg.method          = 'slice';
cfg.camlight        = 'off';
cfg.slicedim        = 1;
cfg.nslices         = 20;
cfg.funcolormap     = 'bipolar';
cfg.funcolorlim     = 'maxabs';

for intrv = 1:size(interval_labels_A, 2)
    for frq = 1:size(freq_label{intrv}, 2)
        for contr = 1:size(contrasts_A, 2)
            ft_sourceplot(cfg, source_ratio_avg_int.(freq_label{intrv}{frq}).([contrasts_A{contr}{1} '_' contrasts_A{contr}{2}]).(interval_labels_A{intrv}));
            set(gcf, 'color', 'black') % black background                     
            set(gcf,'inverthardcopy','off'); %black background also for printed figure
            title([freq_label{intrv}{frq} ' interval ' interval_labels_A{intrv}  ' power for contrast: ' contrasts_A{contr}{1} ' vs. ' contrasts_A{contr}{2}], 'Color', 'w', 'Interpreter', 'None');
            c = colorbar;
            ylabel(c,['[(' contrasts_A{contr}{1} '-' contrasts_A{contr}{2} ')/(' contrasts_A{contr}{1} '+' contrasts_A{contr}{2} ')]']);
            c.Color = 'w';
            print(gcf,fullfile(plots_path, ['DICS_A_GAV_' freq_label{intrv}{frq} '_' interval_labels_A{intrv} '_' contrasts_A{contr}{1} '_' contrasts_A{contr}{2} '_5mm.png']),'-dpng','-r1000')
        end % contrast loop
    end % freqband loop
end % interval loop