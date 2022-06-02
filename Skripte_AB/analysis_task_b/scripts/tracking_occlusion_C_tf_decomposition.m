% Data Analysis Pursuit-Occlusion Paradigm
% Emulation pilot study 2021

% Script C contains:
% TF decomposition for occlusion and reappearance

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
inputpath = "R:\AG-Beste-Studien\Emulation\06_analysis\output_analysis_task_B\02_icaclean";
cd(inputpath);

% set export directory
savepath = "R:\AG-Beste-Studien\Emulation\06_analysis\output_analysis_task_B\03_tf";  

%list all *.set files in inputpath
filenames = dir('*icaclean*.set');

%concatenate into one cell array
files2read = {filenames.name};

%% loop through files generated with script B (IClabel)

for ind = 1:length(filenames)

    % import the data file
    TMPEEG = pop_loadset('filename', files2read(ind), 'filepath', char(inputpath));
    
    

end