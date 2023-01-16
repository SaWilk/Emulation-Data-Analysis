% The neurophysiology of continuous action monitoring

% preprocessing task B, pt. 1

% created by:
% Adriana Boettcher, Cognitive Neurophysiology TU Dresden
% 2022

%% clear workspace
clc;
clear;

%% folders and dependencies

% add path and start EEGlab toolbox
eeglab;
close;

% define directory rawdata_dir for raw EEG data
cd(rawdata_dir);

% define savepath for saving data

%list all *.vhdr files that end with _B in data directory
filenames = dir('*_B.vhdr');

%concatenate into one cell array
files2read = {filenames.name};

%% exclude files with missing data

% define files to be excluded
exclude = {'91L3H_B.vhdr', '607SE_B.vhdr', 'WM87B_B.vhdr'};

% empty vector for indexing files2read
exclude_files = zeros(1, length(files2read));

% loop through files to be excluded and remove from files2read
for i = 1:length(exclude)
    exclude_files = exclude_files + strcmp(files2read, exclude{i});
end

% remove chosen files
files2read(exclude_files~=0) = [];

%% Loop through files
% load data
% apply preprocessing steps
% save preprocessed EEGlab data

for ind = 1:length(files2read) 

    %import raw data from brainvision
    TMPEEG = pop_loadbv(rawdata_dir, files2read{ind});
    
    %apply add_occl_events --> create OCCL var in EEG.event (coding
    %occl = occl/non_occl
    TMPEEG = add_occl_events(TMPEEG);
    
    %remove all events for OCCL = "none" --> occlusion events are the only
    %that remain
    TMPEEG = rm_occl_none(TMPEEG);    

    %define setname (remove file ending)
    TMPEEG.setname = files2read{ind}(1:end-4);
    
    %store original filename in set
    TMPEEG.comment = files2read{ind};
    TMPEEG.group = 'Control';    

    %define condition
    TMPEEG.condition = 'TaskB';

    % define subject
    substring = files2read{ind}(1:end-7);
    TMPEEG.subject = substring;
    
    % Edit channel info, mark Fpz as ref electrode
    TMPEEG = pop_chanedit(TMPEEG, 'append',60,'changefield',{61,'labels','FPz'});
    TMPEEG = pop_chanedit(TMPEEG, 'setref',{'1:60','FPz'});
    TMPEEG.comment = TMPEEG.comment + "  *** edited channel info";

    % downsampling tp 250 Hz
    TMPEEG = pop_resample( TMPEEG, 250);
    TMPEEG.setname = files2read{ind}(1:end-7);
    TMPEEG.comment = TMPEEG.comment + "  *** apply downsampling";

    %hp filter
    TMPEEG = pop_eegfiltnew(TMPEEG, 1, []);
    TMPEEG.comment = TMPEEG.comment + "  *** applied hp filter";

          
    %remove line noise
    TMPEEG = pop_cleanline(TMPEEG, 'bandwidth',2,'chanlist', 1:TMPEEG.nbchan ,...
       'computepower',1,'linefreqs',[50 100],'normSpectrum',0,'p',0.01,'pad',2,...
       'plotfigures',0,'scanforlines',1,'sigtype','Channels','tau',100,'verb',1,...
       'winsize',4,'winstep',4);
    TMPEEG.comment = TMPEEG.comment + "  *** remove line noise";

    %remove bad channels 
    TMPEEG.oldchanslocs = TMPEEG.chanlocs;
    % clean_rawdata parameters: disable highpass, line noise and ASR
    % remove flat channels and channels with minimum channel correlation
    TMPEEG = clean_rawdata(TMPEEG, 5, -1, 0.85, -1, -1, -1);
    TMPEEG.comment = TMPEEG.comment + "  *** remove bad channels";
            
    %lowpass filter
    TMPEEG = pop_eegfiltnew(TMPEEG, [], 40);
    TMPEEG.comment = TMPEEG.comment + "  *** apply lp filter";

    %interpolate removed channels
    TMPEEG = pop_interp(TMPEEG, TMPEEG.oldchanslocs, 'spherical');
    TMPEEG.comment = TMPEEG.comment + "  *** interpolate bad channels";

    %rereference to average reference
    TMPEEG = pop_reref(TMPEEG, [ ]);
    TMPEEG.comment = TMPEEG.comment + "  *** rereferencing";
    
    %save preprocessed data
    TMPEEG = pop_saveset(TMPEEG,'filename',[files2read{ind}(1:end-7) '_B_preprocessed'], 'filepath', char(savepath));
     
end