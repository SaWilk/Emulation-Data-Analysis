% The neurophysiology of continuous action monitoring
% Saskia Wilken, Adriana Böttcher, Nico Adelhöfer, Markus Raab, Sven
% Hoffmann, Christian Beste

% calculate cross spectral density for each condition + appended for contrast
% separately for the two time intervals for experiment 1
% DICS beamformer applied to calculate common spatial filter for contrast:
%       const vs. rand1
% apply spatial filter for each condition
% for TBA, ABA, BBA

% created by:
% Adriana Boettcher, Cognitive Neurophysiology TU Dresden
% 2022

%% configuration
tic

clc;
clearvars;

% contrasts with conditions
contrasts_A = {
    {'const', 'rand1'}, ...
    };

% freq bands, separetely for time intervals 
% (0-500 ms: sign. differences in all freqbands; 500 onwards only theta)
freq_bands      = {
    {[4 7] [8 12] [13 30]}, ... % first interval, all freqs
    {[4 7]} ... % second interval, only theta
    };
freq_label      = {
    {'theta' 'alpha' 'beta'}, ...
    {'theta'}, ...
    };

% two time limits for the intervals (0-500 ms, 500 ms onwards)
interval_labels_A = {'until_500' 'after_500'};

%% path definition

% load fieldtrip toolbox

% load subject list and remove subject with missing data for task B 
load('subjects.mat');

% initialize input and output folder
% input: segemented and preprocessed time domain data

% for significant time window, load CBPT output
output_suffix   = '_split_interv';

%% Headmodel, MRI, leadfield

% load bem headmodel, structural MRI, electrode positions

% prepare MRI
mri = ft_convert_units(mri,'cm');           % better make sure units match
mri = ft_volumereslice([],mri);             % reslice for good measure

% change electrode labels
elec.label = upper(elec.label);             % upper case letters

% align electrode positions to headmodel
[headmodel, sens]       = ft_prepare_vol_sens(vol,elec);
sens                    = ft_convert_units(sens,'cm');          % change sens unit
headmodel               = ft_convert_units(headmodel,'cm');     % change headmodel unit

% create leadfield
cfg                     = [];
cfg.elec                = sens;
cfg.headmodel           = headmodel;
cfg.reducerank          = 3;                                    % default is 3 for EEG, 2 for MEG
cfg.resolution          = 0.5;                                  % use a 3D grid with a 0.5cm (/5mm) resolution
cfg.sourcemodel.unit    = 'cm';
cfg.tight               = 'yes';
leadfield               = ft_prepare_leadfield(cfg,[]);

% check: plot leadfield
figure; hold on;
ft_plot_headmodel(headmodel,'facecolor','cortex','edgecolor','none','facealpha',0.5);
ft_plot_sens(sens,'label','label','fontsize',8,'style','*r')    % this includes electrode labels in size 8, sensors are red *
ft_plot_mesh(leadfield.pos(leadfield.inside,:));                % plot inside positions of leadfield only
view([-90 0 0]);

%% 

% separate data frames for each contrast to be filled in loop
source_const_rand1_theta_until_500     = [];
source_const_rand1_theta_after_500     = [];
source_const_rand1_alpha_until_500     = [];
source_const_rand1_alpha_after_500     = [];
source_const_rand1_beta_until_500      = [];
source_const_rand1_beta_after_500      = [];

for intrv = 1:size(interval_labels_A, 2)
    for frq = 1:size(freq_label{intrv}, 2)
        for contr = 1:size(contrasts_A, 2)
        
            % save significant time window in CBPT for this contrast
            sig_window_start = CBPT_avg.A.(interval_labels_A{intrv}).(freq_label{intrv}{frq}).([contrasts_A{contr}{1} '_' contrasts_A{contr}{2}]).first_sig;
            sig_window_end = CBPT_avg.A.(interval_labels_A{intrv}).(freq_label{intrv}{frq}).([contrasts_A{contr}{1} '_' contrasts_A{contr}{2}]).last_sig;
        
            for subj = 1:size(subjects, 2)
                try
                % save subject ID
                this_subj   = subjects(subj);
        
                % create dataframes for raw data and freq data for this subject
                % will later be appended to overall dataframes (_all)
                raw_data    = [];
                freq_data   = [];
                source_data = [];
        
                for cond = 1:size(contrasts_A{contr}, 2)
        
                    this_cond = contrasts_A{contr}{cond};
        
                    % create file name from subject ID and condition and load
                    % respective file from the inputpath
                    this_filename   = [char(this_subj) '_' char(this_cond)];
                    temp_data       = load([inputpath_segmented_A filesep this_filename]);
                    temp_data       = temp_data.(this_cond);
                    temp_data.label = upper(temp_data.label);
        
                    % shorter trial length (was originally segmented for - 1 to 4
                    % for wavelet analysis)

                    % for the second interval, there are two different
                    % significant time windows: the window for rand1, given
                    % in CBPT_avg (--> .time starts at 0.5 because first
                    % 500 ms were excluded); and the window for const,
                    % which is 0.5 s before --> therefore, if the condition
                    % == const, the time window = sig window-0.5
                    if intrv == 2 && strcmp(contrasts_A{contr}{cond}, 'const')
                        temp_shorter         = reduce_trial_length(temp_data, sig_window_start-0.5, sig_window_end-0.5);
                        raw_data.(this_cond) = temp_shorter;
                    else
                        temp_shorter         = reduce_trial_length(temp_data, sig_window_start, sig_window_end);
                        raw_data.(this_cond) = temp_shorter;
                    end
        
                    % calculate cross-spectral density matrix
                    cfg                     = [];
                    cfg.method              = 'mtmfft';
                    cfg.foilim              = freq_bands{intrv}{frq};                           
                    cfg.tapsmofrq           = 1/(sig_window_end-sig_window_start); % time window = 3 s
                    cfg.taper               = 'hanning';
                    cfg.output              = 'powandcsd'; % power and cross-spectral density
                    cfg.pad                 = 'nextpow2';
                
                    freq = ft_freqanalysis(cfg, temp_shorter);
                    freq_data.(this_cond) = freq;
                end % condition loop
        
                % append data for both conditions and calculate crossspectral
                % density
                raw_data.([contrasts_A{contr}{1} '_' contrasts_A{contr}{2}]) = ft_appenddata([], ...
                    raw_data.(contrasts_A{contr}{1}), raw_data.(contrasts_A{contr}{2}));
        
                cfg                     = [];
                cfg.method              = 'mtmfft';
                cfg.foilim              = freq_bands{intrv}{frq};                           
                cfg.tapsmofrq           = 1/(sig_window_end-sig_window_start); % time window = 3 s
                cfg.taper               = 'hanning';
                cfg.output              = 'powandcsd'; % power and cross-spectral density
                cfg.pad                 = 'nextpow2';
            
                freq = ft_freqanalysis(cfg, raw_data.([contrasts_A{contr}{1} '_' contrasts_A{contr}{2}]));
                freq_data.([contrasts_A{contr}{1} '_' contrasts_A{contr}{2}]) = freq;
        
                % calculate common filter based on appended data (cross-spectral
                % density matrix)
        
                cfg                     = [];
                cfg.method              = 'dics';
                cfg.frequency           = freq_bands{intrv}{frq};      
                cfg.grid                = leadfield;
                cfg.vol                 = headmodel;
                cfg.elec                = sens;
                cfg.senstype            = 'EEG';
                cfg.dics.keepfilter     = 'yes'; % save common spatial filter to apply to conditions
                cfg.dics.realfilter     = 'yes';                                   
                cfg.dics.projectnoise   = 'yes'; % can be used for NAI
                cfg.dics.lambda         = '5%';
                cfg.dics.reducerank     = 3;
                cfg.channel             = {'all'};
                         
                source_common           = ft_sourceanalysis(cfg, freq_data.([contrasts_A{contr}{1} '_' contrasts_A{contr}{2}]));
                source_data.filter      = source_common;
        
                % set common filter in configuration
                cfg.sourcemodel.filter     = source_common.avg.filter;
                % apply common filter to both conditions
                source_data.(contrasts_A{contr}{1}) = ft_sourceanalysis(cfg, freq_data.(contrasts_A{contr}{1}));
                source_data.(contrasts_A{contr}{2}) = ft_sourceanalysis(cfg, freq_data.(contrasts_A{contr}{2}));
        
                % save data for this subject (source)
                if frq == 1
                    if intrv == 1
                        source_const_rand1_theta_until_500{subj} = source_data;
                    elseif intrv == 2
                        source_const_rand1_theta_after_500{subj} = source_data;
                    end
    
                elseif frq == 2
                    if intrv == 1
                        source_const_rand1_alpha_until_500{subj} = source_data;
                    elseif intrv == 2
                        source_const_rand1_alpha_after_500{subj} = source_data;
                    end
    
                elseif frq == 3
                    if intrv == 1
                        source_const_rand1_beta_until_500{subj} = source_data;
                    elseif intrv == 2
                        source_const_rand1_beta_after_500{subj} = source_data;
                    end
                end
        
                catch 
                    disp([char('error for sbj ') num2str(subj)]);
                    continue
                end % try catch statement
    
            end % subject loop
    
        end % contrast loop
    
    end % freq loop

end % interval loop

source_data_all_theta = [];
source_data_all_theta.const_rand1.until_500 = source_const_rand1_theta_until_500;
source_data_all_theta.const_rand1.after_500 = source_const_rand1_theta_after_500;
save([outputpath filesep 'source_data_all_theta' output_suffix], "source_data_all_theta", '-v7.3');

source_data_all_alpha = [];
source_data_all_alpha.const_rand1.until_500 = source_const_rand1_alpha_until_500;
source_data_all_alpha.const_rand1.after_500 = source_const_rand1_alpha_after_500;
save([outputpath filesep 'source_data_all_alpha' output_suffix], "source_data_all_alpha", '-v7.3');

source_data_all_beta = [];
source_data_all_beta.const_rand1.until_500 = source_const_rand1_beta_until_500;
source_data_all_beta.const_rand1.after_500 = source_const_rand1_beta_after_500;
save([outputpath filesep 'source_data_all_beta' output_suffix], "source_data_all_beta", '-v7.3');

toc

