% The neurophysiology of continuous action monitoring
% Saskia Wilken, Adriana Böttcher, Nico Adelhöfer, Markus Raab, Sven
% Hoffmann, Christian Beste

% grand average for conditions per experiment
% Experiment 1: constant and random trajectory
% experiment 2: occluded and visible trajectory parts

% created by:
% Adriana Boettcher, Cognitive Neurophysiology TU Dresden
% 2022

%% 
clc;
clearvars;

% load fieldtrip toolbox

% load local EEG configuration (electrodes, neighbours, layout)

% initialize input and output folder

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

