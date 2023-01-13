% The neurophysiology of continuous action monitoring
% Saskia Wilken, Adriana Böttcher, Nico Adelhöfer, Markus Raab, Sven
% Hoffmann, Christian Beste

% calculate ratio for contrast, for theta band
% calculate grand average for contrast, interpolate data
% for experiment 2

% created by:
% Adriana Boettcher, Cognitive Neurophysiology TU Dresden
% 2022

%% configuration
tic

clc;
clearvars;

% contrasts with conditions
contrasts_B = {
    {'occl', 'nonoccl'} ...
    };

% freq bands
freq_bands      = {[4 7]};
freq_label      = {'theta'};

%% path definition

% load fieldtrip toolbox

% load subject list and remove subject with missing data for task B 
load('subjects.mat');

% initialize input and output folder
% input: segmented and preprocessed time domain data 


%% load MRI for plotting

% load bem headmodel, structural MRI, electrode positions

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
    data = load([datapath filesep 'source_data_all_B_theta']);
    data = data.source_data_all;

    for contr = 1:size(contrasts_B, 2)

        % save conds of this contrast for indexing
        cond1 = contrasts_B{contr}{1};
        cond2 = contrasts_B{contr}{2};

        source_ratio.([char(cond1) '_' char(cond2)]) = data;

        for subj = 1:size(subjects, 2)

            % calculate ratio for conds of this contrast and save to extra
            % df to calculate GAV in next step

            % copy data to save ratio
            source_ratio.([char(cond1) '_' char(cond2)]){subj} = data{subj};
            source_ratio.([char(cond1) '_' char(cond2)]){subj}.ratio = source_ratio.([char(cond1) '_' char(cond2)]){subj}.(cond1);
            
            % calculate ratio and save in copied structure ((x-y)./(x+y))
            source_ratio.([char(cond1) '_' char(cond2)]){subj}.ratio.avg.pow = ...
                (source_ratio.([char(cond1) '_' char(cond2)]){subj}.(cond1).avg.pow - ...
                source_ratio.([char(cond1) '_' char(cond2)]){subj}.(cond2).avg.pow) ./ ...
                (source_ratio.([char(cond1) '_' char(cond2)]){subj}.(cond1).avg.pow + ...
                source_ratio.([char(cond1) '_' char(cond2)]){subj}.(cond2).avg.pow);

        end % subject loop
    
        ratio_this_contr = cellfun(@(x) x.ratio, source_ratio.([char(cond1) '_' char(cond2)]), 'UniformOutput',false);

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

save([datapath filesep 'source_ratio_B_all_avg'], "source_ratio_all_avg", '-v7.3');

%% interpolate data on MRI and save

%load averaged source data
load([datapath filesep 'source_ratio_B_all_avg']);

source_ratio_avg_int = [];

% create new data frame with interpolated data for each freq and contrast
for frq = 1:size(freq_label, 2)
    for contr = 1:size(contrasts_B, 2)

        this_contrast = [contrasts_B{contr}{1} '_' contrasts_B{contr}{2}];

        cfg                         = [];
        cfg.parameter               = 'pow';
        source_ratio_avg_int.(freq_label{frq}).(this_contrast) = ft_sourceinterpolate(cfg, source_ratio_all_avg.(freq_label{frq}).(this_contrast), mri);
    end % contrast loop
end % freq loop

% save interpolated data
save([datapath filesep 'source_ratio_B_avg_int'], "source_ratio_avg_int", '-v7.3');
