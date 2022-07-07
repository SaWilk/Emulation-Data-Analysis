% Data Analysis Pursuit-Tracking and Pursuit-Occlusion Paradigm
% Emulation pilot study 2021

% Script contains:
% IC-label for extracting neural components based on combined data of task
% A and B

% Adriana Boettcher
% 13.06.2022

%% clear workspace
clear;
clc;

%% folders and dependencies

% add path and start EEGlab toolbox
addpath('R:\AG-Beste-Orga\Skripts\Toolbox\eeglab2021.0');
eeglab;
close;

% set input path
inputpath = "R:\AG-Beste-Studien\Emulation\06_analysis\output_ICA_combined\01_ICA";
cd(inputpath);

% set export directory
savepath = "R:\AG-Beste-Studien\Emulation\06_analysis\output_ICA_combined\02_IClabel";  

%list all *.set files in inputpath
filenames = dir('*ICA_merged*.set');

%concatenate into one cell array
files2read = {filenames.name};

%% loop through files with ICA components based on both data sets

for ind = 1:length(filenames)
    
    % import the data file
    TMPEEG = pop_loadset('filename', files2read(ind), 'filepath', char(inputpath));
    
    %get the file name for saving later
    filename = TMPEEG.filename(1:end-4);

    % save copy of old components extracted based on single task
    TMPEEG.icawinv_singletask = TMPEEG.icawinv;
    TMPEEG.icasphere_singletask = TMPEEG.icasphere;
    TMPEEG.icaweights_singletask = TMPEEG.icaweights;

    % save new components based on both tasks with standard variable names
    TMPEEG.icawinv = TMPEEG.icawinv_merged;
    TMPEEG.icasphere = TMPEEG.icasphere_merged;
    TMPEEG.icaweights = TMPEEG.icaweights_merged;

    %run IClabel to classify components
    TMPEEG = iclabel(TMPEEG);

    %extract brain components
    neurocomps  = find(TMPEEG.etc.ic_classification.ICLabel.classifications(:,1) > 0.5);

    %extract occular components
    eyes = find(TMPEEG.etc.ic_classification.ICLabel.classifications(:,3) > 0.6);

    %identify components to exclude (all non-brain components)
    ic2rem = unique([eyes' setdiff(1:size(TMPEEG.icawinv,2), neurocomps)]);

    % save inverse and weights of comps excluded by IClabel
    TMPEEG.IClabel_excl.icawinv = TMPEEG.icawinv(:, ic2rem);
    TMPEEG.IClabel_excl.icasphere = TMPEEG.icasphere(:, ic2rem);
    TMPEEG.IClabel_excl.icaweights = TMPEEG.icaweights(:, ic2rem);

    % keep only brain components
    TMPEEG = pop_subcomp(TMPEEG, ic2rem, 0);
    TMPEEG.setname = [filename '-icaclean'];

    %SHs script contains a check here, todo include this as well
    % visual checks:
    pop_topoplot(TMPEEG, 0, [1:5] ,TMPEEG.setname,[2 3] ,0,'electrodes','on');
    pause;
    ics = input('Number of IC to remove?:');

    %save weights and inverse of manually excluded comps
    TMPEEG.manual_excl.icawinv = TMPEEG.icawinv(:, ics);
    TMPEEG.manual_excl.icasphere = TMPEEG.icasphere(:, ics);
    TMPEEG.manual_excl.icaweights = TMPEEG.icaweights(:, ics);

    pause;
    if isempty(ics)
        TMPEEG = pop_saveset(TMPEEG,'filename',TMPEEG.setname,'filepath', char(savepath));
    else
        TMPEEG = pop_subcomp( TMPEEG, ics, 0);
        TMPEEG = pop_saveset(TMPEEG,'filename',TMPEEG.setname,'filepath', char(savepath));
    end
end
