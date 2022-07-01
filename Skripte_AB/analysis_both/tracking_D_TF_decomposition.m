% Data Analysis Pursuit-Tracking and Pursuit-Occlusion Paradigm
% Emulation pilot study 2021

% This script contains: 
% create condition matrix for comparing average
% calculate ERSP for different trajs, different occl states
% calculate average and plot average
% tbc: statistics

% Adriana Böttcher
% 23.06.22

%% clear workspace
clear;
clc;

%% folders and dependencies

% add path and start EEGlab toolbox
addpath('R:\AG-Beste-Orga\Skripts\Toolbox\eeglab2021.0');
eeglab;
close;

% set input path
inputpath = "R:\AG-Beste-Studien\Emulation\06_analysis\output_ICA_combined_new\03_epoched";
cd(inputpath);

% set export directory
savepath = "R:\AG-Beste-Studien\Emulation\06_analysis\output_ICA_combined_new\05_TF");

%list all *.set files in inputpath
filenames_A = dir('*_ICA_merged_new_icaclean_epoched_A.set');
filenames_B = dir('*_ICA_merged_new_icaclean_epoched_B.set');

%concatenate into one cell array
files2read_A = {filenames_A.name};
files2read_B = {filenames_B.name};

%% loop for A

for sbj_ind = 1:length(files2read_A)

    % import the data file
    TMPEEG_A = pop_loadset('filename', files2read(sbj_ind), 'filepath', char(inputpath));
    
    %get the file name for saving later
    filename = TMPEEG_A.filename(1:5);

    % create empty traj variable for contrast
    traj = [];
    for trial_ind = 1:TMPEEG_A.trials
        tr = TMPEEG_A.epoch(trial_ind).eventTRAJ{find(cell2mat(TMPEEG_A.epoch(trial_ind).eventlatency)==0)};
        if strmatch(tr,'CONSTANT')
                traj(i,1) = 1;
        elseif strmatch(tr,'RANDOM1')
            traj(i,1) = 0;
        elseif strmatch(tr,'RANDOM2')
            traj(i,1) = 2;
        else
            traj(i,1) = NaN;
        end
    end

    %% loop through channels

    for chan_ind = 1:TMPEEG_A:nbchan
        [erspconstant(:,:,channum), ~, ~, ersptimes, erspfreqs] = newtimef( TMPEEG.data(channum,:, traj == 1), ...
                TMPEEG.pnts,[TMPEEG.xmin * 1000  TMPEEG.xmax*1000], TMPEEG.srate,[3   0.5] ,  'elocs', ...
                TMPEEG.chanlocs, 'chaninfo', TMPEEG.chaninfo, 'caption', chan, 'basenorm','on', ...
                'freqs', [1 30],  'nfreqs', 30,'plotphase','off', 'plotersp','off','plotitc','off',...
                'padratio', 1,'trialbase','full'); close;
                     
        [ersprandom(:,:,channum), ~, ~, ersptimes, erspfreqs] = newtimef( TMPEEG.data(channum,:, traj == 0), ...
            TMPEEG.pnts,[TMPEEG.xmin * 1000  TMPEEG.xmax * 1000], TMPEEG.srate,[3   0.5] ,  'elocs', ...
            TMPEEG.chanlocs, 'chaninfo', TMPEEG.chaninfo, 'caption', chan, 'basenorm','on', ...
            'freqs', [1 30],  'nfreqs', 30,'plotphase','off', 'plotersp','off','plotitc','off',...
             'padratio', 1,'trialbase','full'); close;



    end



end

%% loop for B

for i = 1:length(files2read_B)


end


