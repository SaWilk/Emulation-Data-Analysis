% Data Analysis Pursuit-Tracking and Pursuit-Occlusion Paradigm
% Emulation pilot study 2021

% This script contains: 
% * load source data for task A
% * calculate ratio for all contrasts & freq bands per subject
% * calculate grand average for all contrasts
% * interpolate data
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

for frq = 1:size(freq_label, 2)

    % load data for respective frequency band
    data = load([datapath filesep 'source_data_all_' freq_label{frq}]);
    data = data.(['source_data_all_' freq_label{frq}]);
    source_ratio.(freq_label{frq}) = data;

    for contr = 1:size(contrasts_A, 2)

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

%% interpolate data on MRI and save
% 
%load averaged source data
load([datapath filesep 'source_ratio_all_avg']);

source_ratio_avg_int = [];

% create new data frame with interpolated data for each freq and contrast
for frq = 1:size(freq_label, 2)
    for contr = 1%:size(contrasts_A, 2)

        this_contrast = [contrasts_A{contr}{1} '_' contrasts_A{contr}{2}];

        cfg                         = [];
        cfg.parameter               = 'pow';
        source_ratio_avg_int.(freq_label{frq}).(this_contrast) = ft_sourceinterpolate(cfg, source_ratio_all_avg.(freq_label{frq}).(this_contrast), mri);
    end % contrast loop
end % freq loop
% 
%% plotting

% add main title with freqband and condition

cfg                 = [];
cfg.funparameter    = 'pow';
cfg.maskparameter   = 'mask';
%cfg.maskparameter   = cfg.funparameter;
cfg.method          = 'slice';
cfg.camlight        = 'off';
cfg.slicedim        = 1;
cfg.nslices         = 20;
cfg.funcolormap     = bipolar;
cfg.funcolorlim     = 'maxabs';

ft_sourceplot(cfg, source_ratio_avg_int.theta.const_rand1);
set(gcf, 'color', 'black')                              % make it look nice
set(gcf,'inverthardcopy','off');
%set(gcf, 'Position', get(0, 'Screensize'));             % make it big
c = colorbar;
ylabel(c,'RATIO [(const-rand1)/(const+rand1)]');
c.Color = 'w';
print(gcf,fullfile(PATHOUT,['DICS_GA_',conf.band{freqband},'_post_5mm','.png']),'-dpng','-r1000')



% 
%        cfg                         = [];
%         cfg.parameter               = 'pow';
%         sourcePost_RATIO_avg_INT    = ft_sourceinterpolate(cfg, sourcePost_RATIO_avg, mri);
%         
%         cfg                         = [];
%         cfg.parameter               = 'nai';
%         sourcePre_NAI_avg_INT       = ft_sourceinterpolate(cfg, sourcePre_NAI_avg, mri);
%         
%         
%         cfg                 = [];
%         cfg.funparameter    = 'pow';
%         cfg.maskparameter   = 'mask';
%         cfg.method          = 'slice';
%         cfg.camlight        = 'off';
%         cfg.slicedim        = 1;
%         cfg.nslices         = 20;
%         cfg.funcolormap     = bipolar;
%         % cfg.funcolorlim     = [-5 5];
%         cfg.funcolorlim     = 'maxabs';
%         
%         % WITHIN
%         ft_sourceplot(cfg, sourcePost_RATIO_avg_INT);
%         
%         set(gcf, 'color', 'black')                              % make it look nice
%         set(gcf,'inverthardcopy','off');
%         %set(gcf, 'Position', get(0, 'Screensize'));             % make it big
%         c = colorbar;
%         ylabel(c,'RATIO [(SAME-DIFF)/(SAME+DIFF)]');
%         c.Color = 'w';
%         print(gcf,fullfile(PATHOUT,['DICS_GA_',conf.band{freqband},'_post_5mm','.png']),'-dpng','-r1000')
% 
% 
% 
% 
% 
