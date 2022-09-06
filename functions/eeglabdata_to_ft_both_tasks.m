% Data Analysis Pursuit-Tracking and Pursuit-Occlusion Paradigm
% Emulation pilot study 2021

% This script contains: 
% convert EEGlab data format to fieldtrip 

% Adriana Böttcher
% 18.08.22

%% initialize

% add path and start EEGlab toolbox
addpath('R:\AG-Beste-Orga\Skripts\Toolbox\eeglab2021.0');
eeglab;
close;

datapath = 'R:\AG-Beste-Studien\Emulation\06_analysis\Emulation-Data-Output\06_artifact_rejection_b';
cd(datapath);

outputpath = 'R:\AG-Beste-Studien\Emulation\06_analysis\Emulation-Data-Output\07_data_ft';
outputpath_segments_A = 'R:\AG-Beste-Studien\Emulation\06_analysis\Emulation-Data-Output\08_ft_segments_A';
outputpath_segments_B = 'R:\AG-Beste-Studien\Emulation\06_analysis\Emulation-Data-Output\08_ft_segments_B';

%list all *.set files in inputpath
filenames_A = dir('*_complete_preprocessing_withAR_A.set');
filenames_B = dir('*_complete_preprocessing_withAR_B.set');

%concatenate into one cell array
files2read_A = {filenames_A.name};
files2read_B = {filenames_B.name};

%% task A

for sbj_ind = 1:length(files2read_A)
    % import the data file
    TMPEEG_A = pop_loadset('filename', files2read_A(sbj_ind), 'filepath', char(datapath));
    
    %get the file name for saving later
    filename = TMPEEG_A.filename(1:5);

    % transform to ft data format
    % not sure about the 'raw' parameter
    ft_data = eeglab2fieldtrip(TMPEEG_A, 'raw');

    % save to output path
    save([char(outputpath) filesep char(filename) char('_preprocessed_A.mat')], 'ft_data');

    %% select constant and random trials/epochs and save to extra structure
    % add subject loop here or put both loops together
    % for sbj end
    
    % save trial/epoch indices
    is_constant = [];
    is_rand1 = [];
    is_rand2 = [];
    for i = 1:size(ft_data.trialinfo, 1) % loop through trials
        if strcmp(char(ft_data.trialinfo.TRAJ(i)),'CONST')
            is_constant(end+1) = i;
        elseif strcmp(char(ft_data.trialinfo.TRAJ(i)),'RANDOM1')
            is_rand1(end+1) = i;
        elseif strcmp(char(ft_data.trialinfo.TRAJ(i)),'RANDOM2')
            is_rand2(end+1) = i;
        end
    end
    % and create structs for the respective trajectory segments
    const           = struct();
    const.elec      = ft_data.elec;
    const.trial     = ft_data.trial(is_constant);
    const.time      = ft_data.time(is_constant);
    const.label     = ft_data.label;
    const.trialinfo = ft_data.trialinfo(is_constant, :);
    const.cfg       = ft_data.cfg;

    % save to output path
    save([char(outputpath_segments_A) filesep char(filename) char('_const.mat')], 'const');
    
    rand1           = struct();
    rand1.elec      = ft_data.elec;
    rand1.trial     = ft_data.trial(is_rand1);
    rand1.time      = ft_data.time(is_rand1);
    rand1.label     = ft_data.label;
    rand1.trialinfo = ft_data.trialinfo(is_rand1, :);
    rand1.cfg       = ft_data.cfg;
    
    % save to output path
    save([char(outputpath_segments_A) filesep char(filename) char('_rand1.mat')], 'rand1');

    rand2           = struct();
    rand2.elec      = ft_data.elec;
    rand2.trial     = ft_data.trial(is_rand2);
    rand2.time      = ft_data.time(is_rand2);
    rand2.label     = ft_data.label;
    rand2.trialinfo = ft_data.trialinfo(is_rand2, :);
    rand2.cfg       = ft_data.cfg;   

    % save to output path
    save([char(outputpath_segments_A) filesep char(filename) char('_rand2.mat')], 'rand2');
end

%% task B

%until now one subject missing: KMY6K

for sbj_ind = 1:length(files2read_A)
    % import the data file
    TMPEEG_B = pop_loadset('filename', files2read_B(sbj_ind), 'filepath', char(datapath));
    
    %get the file name for saving later
    filename = TMPEEG_B.filename(1:5);

    % transform to ft data format
    % not sure about the 'raw' parameter
    ft_data = eeglab2fieldtrip(TMPEEG_B, 'raw');

    % save to output path
    save([char(outputpath) filesep char(filename) char('_preprocessed_B.mat')], 'ft_data');

    %% select occl and non-occl trials/epochs and save to extra structure
    % add subject loop here or put both loops together
    % for sbj end
    
    % save trial/epoch indices
    is_occl         = [];
    is_nonoccl      = [];

    for i = 1:size(ft_data.trialinfo, 1) % loop through trials
        if strcmp(char(ft_data.trialinfo.OCCL(i)),'ON')
            is_occl(end+1) = i;
        elseif strcmp(char(ft_data.trialinfo.OCCL(i)),'OFF')
            is_nonoccl(end+1) = i;
        end
    end
    % and create structs for the respective trajectory segments
    occl           = struct();
    occl.elec      = ft_data.elec;
    occl.trial     = ft_data.trial(is_occl);
    occl.time      = ft_data.time(is_occl);
    occl.label     = ft_data.label;
    occl.trialinfo = ft_data.trialinfo(is_occl, :);
    occl.cfg       = ft_data.cfg;

    % save to output path
    save([char(outputpath_segments_B) filesep char(filename) char('_occl.mat')], 'occl');
    
    nonoccl           = struct();
    nonoccl.elec      = ft_data.elec;
    nonoccl.trial     = ft_data.trial(is_nonoccl);
    nonoccl.time      = ft_data.time(is_nonoccl);
    nonoccl.label     = ft_data.label;
    nonoccl.trialinfo = ft_data.trialinfo(is_nonoccl, :);
    nonoccl.cfg       = ft_data.cfg;
    
    % save to output path
    save([char(outputpath_segments_B) filesep char(filename) char('_nonoccl.mat')], 'nonoccl');

end