% Data Analysis Pursuit-Occlusion Paradigm
% Emulation pilot study 2021

% Script B contains:
% IC-label for extracting neural components

% Adriana Boettcher
% 02.06.2022

%% clear workspace
clear;
clc;

%% folders and dependencies

% add path and start EEGlab toolbox
addpath('R:\AG-Beste-Orga\Skripts\Toolbox\eeglab2021.0');
eeglab;
close;

% set input path
inputpath = "R:\AG-Beste-Studien\Emulation\06_analysis\output_analysis_task_B\01_preprocessed";
cd(inputpath);

% set export directory
savepath = "R:\AG-Beste-Studien\Emulation\06_analysis\output_analysis_task_B\02_icaclean";  

%list all *.set files in inputpath
filenames = dir('*epoched*.set');

%concatenate into one cell array
files2read = {filenames.name};

%% loop through files generated with script A (preprocessed EEG data)

for ind = 9 %1:length(filenames)
    
    % import the data file
    TMPEEG = pop_loadset('filename', files2read(ind), 'filepath', char(inputpath));
    
    %get the file name for saving later
    filename = TMPEEG.filename(1:end-4);

    %run IClabel to classify components
    TMPEEG = iclabel(TMPEEG);

    %extract brain components
    neurocomps  = find(TMPEEG.etc.ic_classification.ICLabel.classifications(:,1) > 0.5);

    %extract occular components
    eyes = find(TMPEEG.etc.ic_classification.ICLabel.classifications(:,3) > 0.6);

    %identify components to exclude (all non-brain components)
    ic2rem = unique([eyes' setdiff(1:size(TMPEEG.icawinv,2), neurocomps)]);

    % keep only brain components
    TMPEEG = pop_subcomp(TMPEEG, ic2rem, 0);
    TMPEEG.setname = [filename '-icaclean'];

    %SHs script contains a check here, todo include this as well
    % visual checks:
    pop_topoplot(TMPEEG, 0, [1:5] ,TMPEEG.setname,[2 3] ,0,'electrodes','on');
    pause;
    ics = input('Number of IC to remove?:');
    pause;
    if isempty(ics)
        TMPEEG = pop_saveset(TMPEEG,'filename',TMPEEG.setname,'filepath', char(savepath));
    else
        TMPEEG = pop_subcomp( TMPEEG, ics, 0);
        TMPEEG = pop_saveset(TMPEEG,'filename',TMPEEG.setname,'filepath', char(savepath));
    end

    %for now, just save cleaned data:
    %TMPEEG = pop_saveset(TMPEEG,'filename',TMPEEG.setname,'filepath', char(savepath));
end
