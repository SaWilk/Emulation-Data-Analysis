% The neurophysiology of continuous action monitoring

% preprocessing task A&B, pt. 2

% created by:
% Adriana Boettcher, Cognitive Neurophysiology TU Dresden
% 2022

%% clear workspace
clear;
clc;

%% paths & dependencies

% add path and start EEGlab toolbox
eeglab;
close;

% set input path
% set output path

%list all *.set files in inputpath
cd(inputpath_A);
filenames_A = dir('*preprocessed*.set');
cd(inputpath_B);
filenames_B = dir('*preprocessed*.set');

%concatenate into one cell array
files2read_A = {filenames_A.name};
files2read_B = {filenames_B.name};

%% Loop through files
% load data
% merge data
% apply ICA
% save new data

for ind = 1:length(filenames_A)

    % check if every position of subject IDs match for data file A and B
    if sum(char(extractBefore(files2read_A(ind), 6)) == char(extractBefore(files2read_B(ind), 6))) ~= 5
        disp("datasets not matching");
        return
    end

    % load preprocessed EEG data for A and B
    TMPEEG_A = pop_loadset('filename', files2read_A(ind), 'filepath', char(inputpath_A));
    TMPEEG_B = pop_loadset('filename', files2read_B(ind), 'filepath', char(inputpath_B));

    %merge data
    TMPEEG = pop_mergeset(TMPEEG_A, TMPEEG_B, 1);

    % save sbj ID for merged set
    filename = TMPEEG.filename(1:5);

    %% prepare data for ICA (to exclude occular artifacts)
    
    %prepare ICA via data subset

    %create epochs
    ICA = eeg_regepochs(TMPEEG);

    %detrend eeg data
    ICA = eeg_detrend(ICA); % see https://github.com/widmann/erptools/blob/master/eeg_detrend.m

    %artifact rejection
    ICA = pop_jointprob(ICA, 1, 1:ICA.nbchan, 5, 5, 0, 1);
    
    %select only some random trials
    trl = 1:ICA.trials;
    trl = shuffle(trl);
    ICA = pop_select(ICA, 'trial', trl(1:round(length(trl)/4))) ;
    
    %prepare data for ICA 
    x = double(ICA.data);
    
    %reshape for ICA
    x = reshape(x,size(x,1),size(x,2)*size(x,3));
    
    %get rank (use function modified by SH)
    rnk = getrank(x);
    
    % now run ICA
    % extended infomax due to higher efficiency
    ICA = pop_runica(ICA, 'icatype', 'runica', 'extended',1, 'pca',rnk);
    TMPEEG.icaweights   = ICA.icaweights;
    TMPEEG.icasphere    = ICA.icasphere;
    TMPEEG.icawinv      = ICA.icawinv; 
    TMPEEG.icachansind  = ICA.icachansind;

    TMPEEG              = eeg_checkset(TMPEEG);

    %save new data
    TMPEEG = pop_saveset(TMPEEG,'filename',[filename '_ICA_merged_new'], 'filepath', char(savepath));
end
