% Data Analysis Pursuit-Tracking and Pursuit-Occlusion Paradigm
% Emulation pilot study 2021

% This script contains: 
% IC label for combined data

% Adriana Böttcher
% 23.06.22

% update AB 7.7.22: epoching etc. in separate script

% update AB 08/07/22: relative paths

%% clear workspace
clear;
clc;

%% paths & dependencies

% add path and start EEGlab toolbox
addpath('R:\AG-Beste-Orga\Skripts\Toolbox\eeglab2021.0');
eeglab;
close;

% file location
file_path = matlab.desktop.editor.getActiveFilename;
file_path = fileparts(file_path);
addpath(file_path);

% get data paths for parent dirs
filepath_parts = strsplit(file_path, "\");
parent_dir = filepath_parts(1:end-1);
parent_dir = strjoin(parent_dir, "\");
parent_dir_2 = filepath_parts(1:end-2);
parent_dir_2 = strjoin(parent_dir_2, "\");

% add functions in parent dir
addpath([char(parent_dir) filesep char("functions")]);

% set input path
inputpath = [parent_dir_2 filesep "Emulation-Data-Output\01_ICA"];
cd(inputpath);

% set export directory
savepath_IClabel = [parent_dir_2 filesep "Emulation-Data-Output\02_IClabel"];

%list all *.set files in inputpath
filenames = dir('*ICA_merged_new*.set');

%concatenate into one cell array
files2read = {filenames.name};

%% loop through files with ICA components based on both data sets

for ind = 1:length(filenames)
    
    % import the data file
    TMPEEG = pop_loadset('filename', files2read(ind), 'filepath', char(inputpath));
    
    %get the file name for saving later
    filename = TMPEEG.filename(1:end-4);

    %% iclabel

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
    TMPEEG.comment = TMPEEG.comment + "  *** remove ICs automatically";

    % visual checks:
    pop_topoplot(TMPEEG, 0, [1:5] ,TMPEEG.setname,[2 3] ,0,'electrodes','on');
    pause;
    ics = input('Number of IC to remove?:');

    %save weights and inverse of manually excluded comps
    TMPEEG.manual_excl.icawinv = TMPEEG.icawinv(:, ics);
    TMPEEG.manual_excl.icasphere = TMPEEG.icasphere(:, ics);
    TMPEEG.manual_excl.icaweights = TMPEEG.icaweights(:, ics);

    TMPEEG.setname = [filename '_icaclean_continuous'];

    pause;
    if isempty(ics)
        TMPEEG = pop_saveset(TMPEEG,'filename', TMPEEG.setname,'filepath', char(savepath_IClabel));
    else
        TMPEEG = pop_subcomp( TMPEEG, ics, 0);
        TMPEEG.comment = TMPEEG.comment + "  *** remove ICs manually";
        TMPEEG.setname = [filename '_icaclean_continuous'];
        TMPEEG = pop_saveset(TMPEEG,'filename',TMPEEG.setname,'filepath', char(savepath_IClabel));
    end



end
