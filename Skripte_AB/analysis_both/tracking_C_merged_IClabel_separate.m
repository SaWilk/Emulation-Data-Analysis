% Data Analysis Pursuit-Tracking and Pursuit-Occlusion Paradigm
% Emulation pilot study 2021

% This script contains: 
% IC label for combined data
% separating the data for both tasks
% epoching based on respective events
% baseline correction
% remove bad epochs
% exclude epochs before practice trials

% export data at different steps

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
inputpath = "R:\AG-Beste-Studien\Emulation\06_analysis\output_ICA_combined_new\01_ICA";
cd(inputpath);

% set export directory
savepath_IClabel = "R:\AG-Beste-Studien\Emulation\06_analysis\output_ICA_combined_new\02_IClabel";  
savepath_epoched = "R:\AG-Beste-Studien\Emulation\06_analysis\output_ICA_combined_new\03_epoched"; 
savepath_no_practice = "R:\AG-Beste-Studien\Emulation\06_analysis\output_ICA_combined_new\03_b_epoched_no_practice";

savepath_baseline = "R:\AG-Beste-Studien\Emulation\06_analysis\output_ICA_combined_new\04_baseline"; 
savepath_artifactrej = "R:\AG-Beste-Studien\Emulation\06_analysis\output_ICA_combined_new\05_artifact_rejection"; 

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

    %% delete practice trials

    exp_start_marker = 11;

    % first, task A
    % extract the marker_index from urevent S11
    marker_index = 1;
    found = false;
    while found == false
        if str2double(TMPEEG_A.urevent(marker_index).type(2:end)) == exp_start_marker
            found = true;
        else
            marker_index = marker_index + 1;
        end
    end
    % extract event index for this marker in event structure
    event_index = 1;
    found = false;
    while found == false
        if TMPEEG_A.event(event_index).bvmknum > marker_index
            found = true;
        else
            event_index = event_index + 1;
        end
    end
    % reduce events to those after practice, save events with practice to
    % extra cell
    TMPEEG_A.event_with_practice = TMPEEG_A.event;
    TMPEEG_A.event = TMPEEG_A.event_with_practice(event_index:end);
    %save without practice
    TMPEEG_A.setname = [filename '_epoched_nopractice_A'];
    TMPEEG_A = pop_saveset(TMPEEG_A, 'filename', TMPEEG_A.setname, 'filepath', char(savepath_no_practice));

    % second, task B
    % extract the marker_index from urevent S11
    marker_index = 1;
    found = false;
    while found == false
        if str2double(TMPEEG_B.urevent(marker_index).type(2:end)) == exp_start_marker
            found = true;
        else
            marker_index = marker_index + 1;
        end
    end
    % extract event index for this marker in event structure
    event_index = 1;
    found = false;
    while found == false
        if TMPEEG_B.event(event_index).bvmknum > marker_index
            found = true;
        else
            event_index = event_index + 1;
        end
    end
    % reduce events to those after practice, save events with practice to
    % extra cell
    TMPEEG_B.event_with_practice = TMPEEG_B.event;
    TMPEEG_B.event = TMPEEG_B.event_with_practice(event_index:end);
    %save without practicew
    TMPEEG_B.setname = [filename '_epoched_nopractice_B'];
    TMPEEG_B = pop_saveset(TMPEEG_B, 'filename', TMPEEG_B.setname, 'filepath', char(savepath_no_practice));


    %% apply baseline correction and save

    TMPEEG_A = pop_rmbase(TMPEEG_A, [-1000 0] ,[]);
    TMPEEG_A.comment = TMPEEG_A.comment + "  *** apply baseline corr";
    TMPEEG_A.setname = [filename '_complete_preprocessing_A'];
    TMPEEG_A = pop_saveset(TMPEEG_A, 'filename', TMPEEG_A.setname, 'filepath', char(savepath_baseline));


    TMPEEG_B = pop_rmbase(TMPEEG_B, [-1000 0] ,[]);
    TMPEEG_B.comment = TMPEEG_B.comment + "  *** apply baseline corr";
    TMPEEG_B.setname = [filename '_complete_preprocessing_B'];
    TMPEEG_B = pop_saveset(TMPEEG_B, 'filename', TMPEEG_B.setname, 'filepath', char(savepath_baseline));


    %% artifact rejection and save
    
    %parameters: input, ICs (0) or EEG data (1), electrodes, thresholds in
    %stdv for electrodes, global threshold for all electrodes, superpose
    %prelabelling with previous labelling, keep (0) or reject (1) trials

    TMPEEG_A = pop_jointprob(TMPEEG_A, 1, 1:TMPEEG_A.nbchan, 5, 5, 0, 1);
    TMPEEG_A.comment = TMPEEG_A.comment + "  *** apply artifact rejection";
    TMPEEG_A.setname = [filename '_complete_preprocessing_withAR_A'];
    TMPEEG_A = pop_saveset(TMPEEG_A, 'filename', TMPEEG_A.setname, 'filepath', char(savepath_artifactrej));


    TMPEEG_B = pop_jointprob(TMPEEG_B, 1, 1:TMPEEG_B.nbchan, 5, 5, 0, 1);
    TMPEEG_B.comment = TMPEEG_B.comment + "  *** apply artifact rejection";
    TMPEEG_B.setname = [filename '_complete_preprocessing_withAR_B'];
    TMPEEG_B = pop_saveset(TMPEEG_B, 'filename', TMPEEG_B.setname, 'filepath', char(savepath_artifactrej));


end
