% Data Analysis Pursuit-Tracking and Pursuit-Occlusion Paradigm
% Emulation pilot study 2021

% This script contains: 
% Grand average for task A conditions: const, rand1, rand2

% load TF decomposed data for each subject, and calculate GAV

% Adriana Böttcher
% 18.08.22

%% 
clc;
clearvars;

% load fieldtrip toolbox
addpath('R:\AG-Beste-Orga\Skripts\Toolbox\fieldtrip-20210212\');
ft_defaults;

% load custom ft toolbox
addpath('R:\AG-Beste-Orga\Skripts\Fieldtrip\ft_mmcustom');

% load local EEG configuration (electrodes, neighbours, layout)
load fieldtrip_EEG_KJP_elec_61;
load fieldtrip_EEG_KJP_layout_61;
load fieldtrip_EEG_KJP_neighbours_61;

% initialize input and output folder
datapath = 'R:\AG-Beste-Studien\Emulation\06_analysis\Emulation-Data-Output\09_segments_A_TFT';

% load IDs of included subjects
load subjects;

% conditions
conds = {'const', 'rand1', 'rand2'};

% change directory to load data
cd(datapath);

%% prepare data struct to calculate GAV

freq_all = {};

for cond = 1:size(conds,2)
    cond = char(conds(cond));

    % section in struct to fill with subject data
    freq_all.(cond) = {};
    
    for sbj_num = 1:size(subjects,2)
        sbj = char(subjects(sbj_num));
        filename = [char(sbj) '_' char(cond) '_freq'];
        freq = load(filename);
        freq = freq.freq;
        freq_all.(cond){sbj_num} = freq;
    end

end

% save all data in datapath
outputname = [datapath filesep 'freq_all_sbj'];
save(outputname, 'freq_all');

%%  calculate GAV for the three conditions

freq_GAV = {};

for cond = 1:size(conds,2)
    cond = char(conds(cond));
    freq_data_cond = freq_all.(cond);

    cfg = [];
    freq_GAV.(cond) = ft_freqgrandaverage(cfg, freq_data_cond{:});
end


%% multiplot constant

cfg = [];
cfg.layout = lay;
layout = ft_prepare_layout(cfg);

cfg = [];
cfg.ylim = [1 15];
cfg.zlim = 'maxmin';
cfg.layout = layout;
cfg.title = 'Task A: avg const';
ft_multiplotTFR(cfg, freq_GAV.const);

%% multiplot random 1

cfg = [];
cfg.layout = lay;
layout = ft_prepare_layout(cfg);

cfg = [];
cfg.ylim = [1 15];
cfg.zlim = 'maxmin';
cfg.layout = layout;
cfg.title = 'Task A: avg rand1';
ft_multiplotTFR(cfg, freq_GAV.rand1);

%% multiplot random 2

cfg = [];
cfg.layout = lay;
layout = ft_prepare_layout(cfg);

cfg = [];
cfg.ylim = [1 15];
cfg.zlim = 'maxmin';
cfg.layout = layout;
cfg.title = 'Task A: avg rand2';
ft_multiplotTFR(cfg, freq_GAV.rand2);