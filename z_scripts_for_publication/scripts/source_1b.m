% The neurophysiology of continuous action monitoring
% Saskia Wilken, Adriana Böttcher, Nico Adelhöfer, Markus Raab, Sven
% Hoffmann, Christian Beste

% calculate cross spectral density for each condition + appended for contrast
% for experiment 2, conditions: occluded & non_occluded/visible
% DICS beamformer applied to calculate common spatial filter for contrast:
%       occluded vs. visible
% apply spatial filter for each condition
% only for TBA

% created by:
% Adriana Boettcher, Cognitive Neurophysiology TU Dresden
% 2022

%% configuration
tic

clc;
clearvars;

% contrasts with conditions
% just one contrast for B
contrasts_B = {
    {'occl', 'nonoccl'} ...
    };

% freq bands
% only theta band for B, others were not significant
freq_bands      = {[4 7]};
freq_label      = {'theta'};

% time window
sig_window_start = 0;
sig_window_end   = 2;

%% path definition

% load fieldtrip toolbox

% load subject list and remove subject with missing data for task B 
load('subjects.mat');

% initialize input and output folder
% input: segmented and preprocessed time domain data

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

raw_data_all    = [];
freq_data_all   = [];
source_data_all = [];

for frq = 1:size(freq_label, 2)
    for contr = 1:size(contrasts_B, 2)    
        for subj = 1:size(subjects, 2)
            try
            % save subject ID
            this_subj   = subjects(subj);
    
            % create dataframes for raw data and freq data for this subject
            % will later be appended to overall dataframes (_all)
            raw_data    = [];
            freq_data   = [];
            source_data = [];
    
            for cond = 1:size(contrasts_B{contr}, 2)
    
                this_cond = contrasts_B{contr}{cond};
    
                % create file name from subject ID and condition and load
                % respective file from the inputpath
                this_filename   = [char(this_subj) '_' char(this_cond)];
                temp_data       = load([inputpath_segmented_B filesep this_filename]);
                temp_data       = temp_data.(this_cond);
                temp_data.label = upper(temp_data.label);
    
                % shorter trial length (was originally segmented for - 1 to 4
                % for wavelet analysis)
                temp_shorter         = reduce_trial_length(temp_data, sig_window_start, sig_window_end);
                raw_data.(this_cond) = temp_shorter;
    
                % calculate cross-spectral density matrix
                cfg                     = [];
                cfg.method              = 'mtmfft';
                cfg.foilim              = freq_bands{frq};                           
                cfg.tapsmofrq           = 1/(sig_window_end-sig_window_start); % time window = 2 s
                cfg.taper               = 'hanning';
                cfg.output              = 'powandcsd'; % power and cross-spectral density
                cfg.pad                 = 'nextpow2';
            
                freq = ft_freqanalysis(cfg, temp_shorter);
                freq_data.(this_cond) = freq;
            end % condition loop
    
            % append data for both conditions and calculate crossspectral
            % density
            raw_data.([contrasts_B{contr}{1} '_' contrasts_B{contr}{2}]) = ft_appenddata([], ...
                raw_data.(contrasts_B{contr}{1}), raw_data.(contrasts_B{contr}{2}));
    
            cfg                     = [];
            cfg.method              = 'mtmfft';
            cfg.foilim              = freq_bands{frq};                           
            cfg.tapsmofrq           = 1/(sig_window_end-sig_window_start); % time window = 2 s
            cfg.taper               = 'hanning';
            cfg.output              = 'powandcsd'; % power and cross-spectral density
            cfg.pad                 = 'nextpow2';
        
            freq = ft_freqanalysis(cfg, raw_data.([contrasts_B{contr}{1} '_' contrasts_B{contr}{2}]));
            freq_data.([contrasts_B{contr}{1} '_' contrasts_B{contr}{2}]) = freq;
    
            % calculate common filter based on appended data (cross-spectral
            % density matrix)
    
            cfg                     = [];
            cfg.method              = 'dics';
            cfg.frequency           = freq_bands{frq};      
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
                     
            source_common           = ft_sourceanalysis(cfg, freq_data.([contrasts_B{contr}{1} '_' contrasts_B{contr}{2}]));
            source_data.filter      = source_common;
    
            % set common filter in configuration
            cfg.sourcemodel.filter     = source_common.avg.filter;
            % apply common filter to both conditions
            source_data.(contrasts_B{contr}{1}) = ft_sourceanalysis(cfg, freq_data.(contrasts_B{contr}{1}));
            source_data.(contrasts_B{contr}{2}) = ft_sourceanalysis(cfg, freq_data.(contrasts_B{contr}{2}));
    
            % save data for this subject (source)
            % not separated for freqs and contrasts bc there is only one
            % freq band and one contrast

            source_data_all{subj}   = source_data;
            freq_data_all{subj}     = freq_data;
            raw_data_all{subj}      = raw_data;

            catch 
                disp([char('error for sbj ') num2str(subj)]);
                continue
            end % try catch statement

        end % subject loop

    end % contrast loop

end % freq loop

save([outputpath filesep 'source_data_all_B_theta'], "source_data_all", '-v7.3');

toc
