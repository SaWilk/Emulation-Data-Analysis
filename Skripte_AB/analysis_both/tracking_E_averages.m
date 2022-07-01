% Data Analysis Pursuit-Tracking and Pursuit-Occlusion Paradigm
% Emulation pilot study 2021

% This script contains: 
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
datapath = "R:\AG-Beste-Studien\Emulation\06_analysis\output_ICA_combined_new\06_TF";
cd(datapath);

%%

load('ERSP_A_all_subjects.mat');
load('ERSP_B_all_subjects.mat');

%% average for each subject (A)

for id = 1:length(ERSP_A.subject)
    average_ersp_constant(:,:,:,id) = ERSP_A.constant{id};
    average_ersp_random(:,:,:,id) = ERSP_A.random{id};
    average_ersp_random2(:,:,:,id) = ERSP_A.random2{id};
end

%grand average
 ga_ersps_constant = mean(average_ersp_constant,4);
 ga_ersps_random = mean(average_ersp_random,4);
 ga_ersps_random2 = mean(average_ersp_random2,4);


chan = 'Cz';
 channr = find(strcmp(chan,{ERSP_A.chanlocs{1}.labels}));
  figure; 
   subplot(1,3,1);imagesc(ERSP_A.times{1},ERSP_A.freqs{1},squeeze(ga_ersps_random(:,:,channr)));title([chan '-random']);set(gca,'YDir','normal')
   subplot(1,3,2);imagesc(ERSP_A.times{1},ERSP_A.freqs{1},squeeze(ga_ersps_random2(:,:,channr)));title([chan '-random2']);set(gca,'YDir','normal')
   subplot(1,3,3);imagesc(ERSP_A.times{1},ERSP_A.freqs{1},squeeze(ga_ersps_constant(:,:,channr)));title([chan '-constant']);set(gca,'YDir','normal')


%% average for each subject (B)

for id = 1:length(ERSP_B.subject)
    average_ersp_occlon(:,:,:,id) = ERSP_B.occl_on{id};
    average_ersp_occloff(:,:,:,id) = ERSP_B.occl_off{id};
end

%grand average
 ga_ersps_occlon = mean(average_ersp_occlon,4);
 ga_ersps_occloff = mean(average_ersp_occloff,4);


chan = 'FCz';
 channr = find(strcmp(chan,{ERSP_B.chanlocs{1}.labels}));
  figure; 
   subplot(1,2,1);imagesc(ERSP_B.times{1},ERSP_B.freqs{1},squeeze(ga_ersps_occlon(:,:,channr)));title([chan '-occlusion on']);set(gca,'YDir','normal')
   subplot(1,2,2);imagesc(ERSP_B.times{1},ERSP_B.freqs{1},squeeze(ga_ersps_occloff(:,:,channr)));title([chan '-occlusion off']);set(gca,'YDir','normal')


