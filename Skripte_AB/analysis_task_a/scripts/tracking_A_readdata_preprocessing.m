% Data Analysis Pursuit-Tracking Paradigm
% Emulation pilot study 2021

% Script A contains:
% - read raw EEG data
% - preprocessing (filter, re-referencing, ICA)
% - create epochs for occluded and non-occluded intervals
% - export in EEGlab format

% Adriana Boettcher
% 06.06.2022

%% clear workspace
clear;
clc;

%% folders and dependencies

% add path and start EEGlab toolbox
addpath('R:\AG-Beste-Orga\Skripts\Toolbox\eeglab2021.0');
eeglab;
close;

% get file dir and parent dir to create output folder
file_path = matlab.desktop.editor.getActiveFilename;
file_dir = fileparts(file_path);
parent_dir = fileparts(file_dir);

% add path for functions
addpath([parent_dir '\functions']);

% datadir for raw EEG data
rawdata_dir = 'R:\AG-Beste-Studien\Emulation\04_data\EEG\raw';
cd(rawdata_dir);

% data path for data saving
savepath = "R:\AG-Beste-Studien\Emulation\06_analysis\output_analysis_task_A\01_preprocessed";  

%list all *.vhdr files in data directory
filenames = dir('*A*.vhdr');

%concatenate into one cell array
files2read = {filenames.name};

% define empty bad data cell for collecting bad data files
BADDATA = {};

%% Loop through files
% load data
% apply preprocessing steps
% save preprocessed EEGlab data

for ind = 1:length(filenames) 
    
    %import raw data from brainvision
    TMPEEG = pop_loadbv(rawdata_dir, files2read{ind});
    
    %apply add_traj_events --> create TRAJ var in EEG.event (coding
    %traj = rand1/const/rand2
    TMPEEG = add_traj_events(TMPEEG);
    
    %remove all events for TRAJ = "none" --> occlusion events are the only
    %that remain
    TMPEEG = rm_traj_none(TMPEEG);    

    %define setname (remove file ending)
    TMPEEG.setname = files2read{ind}(1:end-4);
    
    %store original filename in set
    TMPEEG.comment = files2read{ind};
    TMPEEG.group = 'Control';    

    %define condition
    TMPEEG.condition = 'TaskA';

    % define subject
    substring = files2read{ind}(1:end-7);
    TMPEEG.subject = substring;
    
    % Edit channel info, mark Fpz as ref electrode
    TMPEEG = pop_chanedit(TMPEEG, 'append',60,'changefield',{61,'labels','FPz'});
    TMPEEG = pop_chanedit(TMPEEG, 'setref',{'1:60','FPz'});

    % downsampling tp 250 Hz
    TMPEEG = pop_resample( TMPEEG, 250);
    TMPEEG.setname = files2read{ind}(1:end-7);

    %hp filter
    TMPEEG = pop_eegfiltnew(TMPEEG, 1, []);

          
    %remove line noise
    TMPEEG = pop_cleanline(TMPEEG, 'bandwidth',2,'chanlist', 1:TMPEEG.nbchan ,...
       'computepower',1,'linefreqs',[50 100],'normSpectrum',0,'p',0.01,'pad',2,...
       'plotfigures',0,'scanforlines',1,'sigtype','Channels','tau',100,'verb',1,...
       'winsize',4,'winstep',4);
          


    %remove bad channels 
    TMPEEG.oldchanslocs = TMPEEG.chanlocs;
    TMPEEG = clean_rawdata(TMPEEG, 5, -1, 0.85, -1, -1, -1);
            
    %lowpass filter
    TMPEEG = pop_eegfiltnew(TMPEEG, [], 40);

    %interpolate removed channels
    TMPEEG = pop_interp(TMPEEG, TMPEEG.oldchanslocs, 'spherical');

    %rereference to average reference
    TMPEEG = pop_reref(TMPEEG, [ ]);
    
    %% prepare data for ICA (to exclude occular artifacts)
    
    %prepare ICA via data subset
    ICA = eeg_regepochs(TMPEEG);
    
    %detrend eeg data
    ICA = eeg_detrend(ICA); %code copied from https://github.com/widmann/erptools/blob/master/eeg_detrend.m

    %artifact rejection
    ICA = pop_jointprob(ICA, 1, 1:ICA.nbchan, 5, 5, 0, 1);
    
    %select only some random trials
    trl = 1:ICA.trials;
    trl = shuffle(trl);
    ICA = pop_select(ICA, 'trial', trl(1:round(length(trl)/2))) ;
    
    %prepare data for ICA 
    x = double(ICA.data);
    
    %reshape for ICA
    x = reshape(x,size(x,1),size(x,2)*size(x,3));
    
    %get rank (use function modified by SH)
    rnk = getrank(x);
    
    % now run ICA
    % extended infomax due to higher efficiency
    ICA = pop_runica(ICA, 'icatype', 'runica', 'extended',1, 'pca',rnk);
    TMPEEG.icaweights = ICA.icaweights;
    TMPEEG.icasphere = ICA.icasphere;
    TMPEEG.icawinv = ICA.icawinv; 
    TMPEEG.icachansind = ICA.icachansind;
    TMPEEG = eeg_checkset(TMPEEG);
   

     %save continuous ica  data
     %TMPEEG = pop_saveset(TMPEEG,'filename',[files2read{ind}(1:end-7) '_B_runica'], 'filepath', char(savepath));
     

    % create epochs for constant and random traj
    % set epoch limits to 0-2 
    TMPEEG = pop_epoch( TMPEEG, {'S 23' 'S 24' 'S 27'}, [0 2], 'newname', [TMPEEG.setname '_A_epoched'], 'epochinfo', 'yes');
    
    %baseline correction was applied for const/rand, skip that here
    %TMPEEG = pop_rmbase( TMPEEG, [-1000 0] ,[]);
    
    %save epoched data
    TMPEEG = pop_saveset(TMPEEG,'filename',[files2read{ind}(1:end-7) '_A_epoched'], 'filepath', char(savepath));
     
end
