% Data Analysis Pursuit-Tracking and Pursuit-Occlusion Paradigm
% Emulation pilot study 2021

% This script contains: 
% Grand average for task A conditions: const, rand1, rand2
% Grand average for task B conditions: occl, non-occl
% save GAVs for A and B

% load TF decomposed data for each subject, and calculate GAV

% adjusted for plots 27/10/22: GAV for data with adjusted padding for
% figures

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
datapath_A = 'R:\AG-Beste-Studien\Emulation\06_analysis\Emulation-Data-Output\09_B_segments_A_TFT_adjusted_for_figures';
datapath_B = 'R:\AG-Beste-Studien\Emulation\06_analysis\Emulation-Data-Output\09_B_segments_B_TFT_adjusted_for_figures';

% load IDs of included subjects
load subjects;

% conditions
conds_A = {'const', 'rand1', 'rand2'};
conds_B = {'occl', 'nonoccl'};

%% prepare data struct to calculate GAV for Task A

% change directory to load data
cd(datapath_A);

freq_all_A = {};

for cond = 1:size(conds_A,2)
    cond = char(conds_A(cond));

    % section in struct to fill with subject data
    freq_all.(cond) = {};
    
    for sbj_num = 1:size(subjects,2)
        sbj = char(subjects(sbj_num));
        filename = [char(sbj) '_' char(cond) '_freq'];
        freq = load(filename);
        freq = freq.freq;
        freq_all_A.(cond){sbj_num} = freq;
    end

end

% save all data in datapath
outputname = [datapath_A filesep 'freq_all_sbj_A'];
save(outputname, 'freq_all_A');

%%  calculate GAV for the three conditions

freq_GAV_A = {};

for cond = 1:size(conds_A,2)
    cond = char(conds_A(cond));
    freq_data_cond = freq_all_A.(cond);

    cfg = [];
    freq_GAV_A.(cond) = ft_freqgrandaverage(cfg, freq_data_cond{:});
end

% save GAV in datapath
outputname = [datapath_A filesep 'freq_GAV_A'];
save(outputname, 'freq_GAV_A');

%% prepare data struct to calculate GAV for Task B

% change directory to load data
cd(datapath_B);

% todo: remove:
% just for now: remove sbj 18 (KMY6K) since data could not be transferred
% to ft
subjects(18) = [];

freq_all_B = {};

for cond = 1:size(conds_B,2)
    cond = char(conds_B(cond));

    % section in struct to fill with subject data
    freq_all.(cond) = {};
    
    for sbj_num = 1:size(subjects,2)
        sbj = char(subjects(sbj_num));
        filename = [char(sbj) '_' char(cond) '_freq'];
        freq = load(filename);
        freq = freq.freq;
        freq_all_B.(cond){sbj_num} = freq;
    end

end

% save all data in datapath
outputname = [datapath_B filesep 'freq_all_sbj_B'];
save(outputname, 'freq_all_B');

%%  calculate GAV for the two conditions

freq_GAV_B = {};

for cond = 1:size(conds_B,2)
    cond = char(conds_B(cond));
    freq_data_cond = freq_all_B.(cond);

    cfg = [];
    freq_GAV_B.(cond) = ft_freqgrandaverage(cfg, freq_data_cond{:});
end

% save GAV in datapath
outputname = [datapath_B filesep 'freq_GAV_B'];
save(outputname, 'freq_GAV_B');
