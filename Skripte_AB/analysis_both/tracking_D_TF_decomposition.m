% Data Analysis Pursuit-Tracking and Pursuit-Occlusion Paradigm
% Emulation pilot study 2021

% This script contains: 
% create condition matrix for comparing average
% calculate ERSP for different trajs, different occl states
% calculate average and plot average
% tbc: statistics

% Adriana Böttcher
% 28.06.22

%% clear workspace
clear;
clc;

%% folders and dependencies

% add path and start EEGlab toolbox
addpath('R:\AG-Beste-Orga\Skripts\Toolbox\eeglab2021.0');
eeglab;
close;

% set input path
inputpath = "R:\AG-Beste-Studien\Emulation\06_analysis\output_ICA_combined_new\05_artifact_rejection";
cd(inputpath);

% set export directory
savepath = "R:\AG-Beste-Studien\Emulation\06_analysis\output_ICA_combined_new\06_TF";

%list all *.set files in inputpath
filenames_A = dir('*_complete_preprocessing_withAR_A.set');
filenames_B = dir('*_complete_preprocessing_withAR_B.set');

%concatenate into one cell array
files2read_A = {filenames_A.name};
files2read_B = {filenames_B.name};

%% loop for A

for sbj_ind = 1%:length(files2read_A)

    % import the data file
    TMPEEG_A = pop_loadset('filename', files2read_A(sbj_ind), 'filepath', char(inputpath));
    
    %get the file name for saving later
    filename = TMPEEG_A.filename(1:5);

    % create empty traj variable for contrast
    traj = [];
    for trial_ind = 1:TMPEEG_A.trials
        tr = TMPEEG_A.epoch(trial_ind).eventTRAJ{find(cell2mat(TMPEEG_A.epoch(trial_ind).eventlatency)==0)};
        if strmatch(tr,'CONSTANT')
                traj(trial_ind,1) = 1;
        elseif strmatch(tr,'RANDOM1')
            traj(trial_ind,1) = 0;
        elseif strmatch(tr,'RANDOM2')
            traj(trial_ind,1) = 2;
        else
            traj(trial_ind,1) = NaN;
        end
    end

    %% loop through channels

    for chan_ind = 1:TMPEEG_A.nbchan
        [erspconstant(:,:,chan_ind), ~, ~, ersptimes, erspfreqs] = newtimef(TMPEEG_A.data(chan_ind,:, traj == 1), ...
                TMPEEG_A.pnts,[TMPEEG_A.xmin * 1000  TMPEEG_A.xmax*1000], TMPEEG_A.srate,[3   0.5] ,  'elocs', ...
                TMPEEG_A.chanlocs, 'chaninfo', TMPEEG_A.chaninfo,'basenorm','on', ...
                'freqs', [1 30],  'nfreqs', 30,'plotphase','off', 'plotersp','off','plotitc','off',...
                'padratio', 1,'trialbase','full'); close;
                     
        [ersprandom(:,:,chan_ind), ~, ~, ersptimes, erspfreqs] = newtimef(TMPEEG_A.data(chan_ind,:, traj == 0), ...
            TMPEEG_A.pnts,[TMPEEG_A.xmin * 1000  TMPEEG_A.xmax * 1000], TMPEEG_A.srate,[3   0.5] ,  'elocs', ...
            TMPEEG_A.chanlocs, 'chaninfo', TMPEEG_A.chaninfo, 'basenorm','on', ...
            'freqs', [1 30],  'nfreqs', 30,'plotphase','off', 'plotersp','off','plotitc','off',...
             'padratio', 1,'trialbase','full'); close;

      [ersprandom2(:,:,chan_ind),~, ~, ersptimes, erspfreqs] = newtimef(TMPEEG_A.data(chan_ind,:, traj == 2), ...
            TMPEEG_A.pnts,[TMPEEG_A.xmin * 1000  TMPEEG_A.xmax * 1000], TMPEEG_A.srate,[3   0.5] ,  'elocs', ...
            TMPEEG_A.chanlocs, 'chaninfo', TMPEEG_A.chaninfo, 'basenorm','on', ...
            'freqs', [1 30],  'nfreqs', 30,'plotphase','off', 'plotersp','off','plotitc','off',...
             'padratio', 1,'trialbase','full'); close;
    end


    %restructure for tftopo
    ERSP_A.constant{sbj_ind} = permute(erspconstant,[2,1,3]);
    ERSP_A.random{sbj_ind} = permute(ersprandom,[2,1,3]);
    ERSP_A.random2{sbj_ind} = permute(ersprandom2,[2,1,3]);

    ERSP_A.subject{sbj_ind} = TMPEEG_A.subject;
    ERSP_A.chanlocs{sbj_ind} = TMPEEG_A.chanlocs;
    ERSP_A.times{sbj_ind} = ersptimes;
    ERSP_A.freqs{sbj_ind} = erspfreqs;
    ERSP_A.TRAJ{sbj_ind} = traj;
    clear ersprandom2 ersprandom erspconstant;
    
end

%save data
save([char(savepath) filesep char('ERSP_A_all_subjects')], 'ERSP_A');

%% loop for B

for sbj_ind = 1:length(files2read_B)

    % import the data file
    TMPEEG_B = pop_loadset('filename', files2read_B(sbj_ind), 'filepath', char(inputpath));
    
    %get the file name for saving later
    filename = TMPEEG_B.filename(1:5);

    % create empty occl variable for contrast
    occl = [];
    for trial_ind = 1:TMPEEG_B.trials
        tr = TMPEEG_B.epoch(trial_ind).eventOCCL{find(cell2mat(TMPEEG_B.epoch(trial_ind).eventlatency)==0)};
        if strmatch(tr,'ON')
            occl(trial_ind,1) = 1;
        elseif strmatch(tr,'OFF')
            occl(trial_ind,1) = 0;
        else
            occl(trial_ind,1) = NaN;
        end
    end

    %% loop through channels

    for chan_ind = 1:TMPEEG_B.nbchan
        [erspocclon(:,:,chan_ind), ~, ~, ersptimes, erspfreqs] = newtimef(TMPEEG_B.data(chan_ind,:, occl == 1), ...
                TMPEEG_B.pnts,[TMPEEG_B.xmin * 1000  TMPEEG_B.xmax*1000], TMPEEG_B.srate,[3   0.5] ,  'elocs', ...
                TMPEEG_B.chanlocs, 'chaninfo', TMPEEG_B.chaninfo,'basenorm','on', ...
                'freqs', [1 30],  'nfreqs', 30,'plotphase','off', 'plotersp','off','plotitc','off',...
                'padratio', 1,'trialbase','full'); close;
                     
        [erspoccloff(:,:,chan_ind), ~, ~, ersptimes, erspfreqs] = newtimef(TMPEEG_B.data(chan_ind,:, occl == 0), ...
            TMPEEG_B.pnts,[TMPEEG_B.xmin * 1000  TMPEEG_B.xmax * 1000], TMPEEG_B.srate,[3   0.5] ,  'elocs', ...
            TMPEEG_B.chanlocs, 'chaninfo', TMPEEG_B.chaninfo, 'basenorm','on', ...
            'freqs', [1 30],  'nfreqs', 30,'plotphase','off', 'plotersp','off','plotitc','off',...
             'padratio', 1,'trialbase','full'); close;
    end


    %restructure for tftopo
    ERSP_B.occl_on{sbj_ind} = permute(erspocclon,[2,1,3]);
    ERSP_B.occl_off{sbj_ind} = permute(erspoccloff,[2,1,3]);

    ERSP_B.subject{sbj_ind} = TMPEEG_B.subject;
    ERSP_B.chanlocs{sbj_ind} = TMPEEG_B.chanlocs;
    ERSP_B.times{sbj_ind} = ersptimes;
    ERSP_B.freqs{sbj_ind} = erspfreqs;
    ERSP_B.OCCL{sbj_ind} = occl;
    clear erspoccloff erspocclon;
    
end

%save data
save([char(savepath) filesep char('ERSP_B_all_subjects')], 'ERSP_B');