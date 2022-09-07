% Data Analysis Pursuit-Tracking and Pursuit-Occlusion Paradigm
% Emulation pilot study 2021

% This script contains: 
% DICS beamformer applied for constant vs. random trajectory segment in
% task A
% for TBA, ABA, BBA

% Adriana Böttcher
% 06.09.2022

%%
clc;
clearvars;

% load fieldtrip toolbox
addpath('R:\AG-Beste-Orga\Skripts\Toolbox\fieldtrip-20210212\');
ft_defaults;

% load custom ft toolbox
addpath('R:\AG-Beste-Orga\Skripts\Fieldtrip\ft_mmcustom');

% conditions
conds_A = {'const', 'rand1', 'rand2'};
conds_B = {'occl', 'nonoccl'};

% initialize input and output folder
% input: segemented and preprocessed time domain data
inputpath_preprocessed    = 'R:\AG-Beste-Studien\Emulation\06_analysis\Emulation-Data-Output\07_data_ft';
inputpath_segmented_A     = 'R:\AG-Beste-Studien\Emulation\06_analysis\Emulation-Data-Output\08_ft_segments_A';
inputpath_segmented_B     = 'R:\AG-Beste-Studien\Emulation\06_analysis\Emulation-Data-Output\08_ft_segments_B';

outputpath      = 'R:\AG-Beste-Studien\Emulation\06_analysis\Emulation-Data-Output\11_DICS';

%% Headmodel, MRI, leadfield

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
% figure; hold on;
% ft_plot_headmodel(headmodel,'facecolor','cortex','edgecolor','none','facealpha',0.5);
% ft_plot_sens(sens,'label','label','fontsize',8,'style','*r')    % this includes electrode labels in size 8, sensors are red *
% ft_plot_mesh(leadfield.pos(leadfield.inside,:));                % plot inside positions of leadfield only
% view([-90 0 0]);

%% 
% tf analysis with output: cross-spectral density and power
% for all conds (07_data_ft) and separately for conds (08_segmented_blabla)

% todo load data (across conditions)
% load data for each condition and save to data struct --> data.const,
% data.rand1, data.rand2

% 
% for cond = 1:size(conds_A,2)
%     cond = char(conds_A(cond));
%     
%     for sbj = 1:size(subjects,2)
%         sbj = char(subjects(sbj));
%         filename = [char(sbj) '_' char(cond)];
%         data = load(filename);
%         
%         % get to data struct
%         data = data.(cond);
%         



% code paul:

% cfg                     = [];
% cfg.method              = 'mtmfft';
% cfg.foilim              = [4 7];                                % theta
% cfg.tapsmofrq           = 1/conf.twd;                           % where conf.twd = length of time window in seconds
% cfg.taper               = 'hanning';
% cfg.output              = 'powandcsd';                          % specify that you want the cross-spectral density matrix
% cfg.pad                 = 'nextpow2';
%                
% freq_all                = ft_freqanalysis(cfg, data_all);       % For common spatial filter
% freq.GO                 = ft_freqanalysis(cfg, data.GO);        % condition 1
% freq.NOGO               = ft_freqanalysis(cfg, data.NOGO);      % condition 2



%%
% create common spatial filter

% pauls code:

% 
% cfg                     = [];
% cfg.method              = 'dics';
% cfg.frequency           = [4 7];                                % theta frequency band
% cfg.grid                = leadfield;
% cfg.vol                 = headmodel;
% cfg.elec                = sens;
% cfg.senstype            = 'EEG';
% cfg.dics.keepfilter     = 'yes';                                % important, if you want to apply a common spatial filter to conditions
% cfg.dics.realfilter     = 'yes';                                   
% cfg.dics.projectnoise   = 'yes';                                % important, can be used to calculate Neural Activity Index (NAI)
% cfg.dics.lambda         = '5%';
% cfg.dics.reducerank     = 3;
% cfg.channel             = {'all'};
%          
% source_all              = ft_sourceanalysis(cfg, freq_all);

%% 
% apply common spatial filter for separate conditions

