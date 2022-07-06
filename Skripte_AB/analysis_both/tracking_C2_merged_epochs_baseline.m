% Data Analysis Pursuit-Tracking and Pursuit-Occlusion Paradigm
% Emulation pilot study 2021

% This script contains: extra script for new epochs and baseline correction
% epoching based on respective events
% baseline correction
% export data at different steps

% Adriana Böttcher
% 24.06.22

%% clear workspace
clear;
clc;

%% folders and dependencies

% add path and start EEGlab toolbox
addpath('R:\AG-Beste-Orga\Skripts\Toolbox\eeglab2021.0');
eeglab;
close;

% set input path
inputpath = "R:\AG-Beste-Studien\Emulation\06_analysis\output_ICA_combined_new\02_IClabel";
cd(inputpath);

% set export directory 
savepath_epoched = "R:\AG-Beste-Studien\Emulation\06_analysis\output_ICA_combined_new\03_epoched";  
savepath_baseline = "R:\AG-Beste-Studien\Emulation\06_analysis\output_ICA_combined_new\04_baseline";  

%list all *.set files in inputpath
filenames = dir('*new_icaclean_continuous*.set');

%concatenate into one cell array
files2read = {filenames.name};


for ind = 1:length(filenames)
    
    % import the data file
    TMPEEG = pop_loadset('filename', files2read(ind), 'filepath', char(inputpath));
    
    %get the file name for saving later
    filename = TMPEEG.filename(1:5);
    
    %% separate datasets & epochs
    
    %epoch continuous data for constant and random traj and save as
    %TMPEEG_A
    TMPEEG_A = pop_epoch( TMPEEG, {'S 23' 'S 24' 'S 27'}, [-1.5 3], 'newname', [TMPEEG.setname '_A_epoched'], 'epochinfo', 'yes');
    TMPEEG_A.comment = TMPEEG_A.comment + "  *** epoch for traj parts";
    TMPEEG_A.setname = [filename '_icaclean_epoched_A'];
    TMPEEG_A = pop_saveset(TMPEEG_A, 'filename', TMPEEG_A.setname, 'filepath', char(savepath_epoched));

    %epoch continuous data for occlusion/ non-occlusion and save as
    %TMPEEG_B
    TMPEEG_B = pop_epoch( TMPEEG, {'S 20' 'S 21'}, [-1.5 2], 'newname', [TMPEEG.setname '_B_epoched'], 'epochinfo', 'yes');
    TMPEEG_B.comment = TMPEEG_B.comment + "  *** epoch for occl parts";
    TMPEEG_B.setname = [filename '_icaclean_epoched_B'];
    TMPEEG_B = pop_saveset(TMPEEG_B, 'filename', TMPEEG_B.setname, 'filepath', char(savepath_epoched));

    %% apply baseline correction and save

    TMPEEG_A = pop_rmbase(TMPEEG_A, [-1000 0] ,[]);
    TMPEEG_A.comment = TMPEEG_A.comment + "  *** apply baseline corr";
    TMPEEG_A.setname = [filename '_complete_preprocessing_A'];
    TMPEEG_A = pop_saveset(TMPEEG_A, 'filename', TMPEEG_A.setname, 'filepath', char(savepath_baseline));

    TMPEEG_B = pop_rmbase(TMPEEG_B, [-1000 0] ,[]);
    TMPEEG_B.comment = TMPEEG_B.comment + "  *** apply baseline corr";
    TMPEEG_B.setname = [filename '_complete_preprocessing_B'];
    TMPEEG_B = pop_saveset(TMPEEG_B, 'filename', TMPEEG_B.setname, 'filepath', char(savepath_baseline));

end