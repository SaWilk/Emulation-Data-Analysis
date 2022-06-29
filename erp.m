%% Make ERP Plot

% Author: Saskia Wilken
% creation Date: 28.06.2022

%% Empty

format compact
format long G
clear
clc


%% Set Paths

% file location
file_path = matlab.desktop.editor.getActiveFilename;
file_path = fileparts(file_path);

data_path = strcat(file_path, "\01_merged_data"); % get root path to data folder

cd(data_path);
% get all set files in the folder
files2read = {dir('*.set').name};

% load all eeg files
ALLEEG = pop_loadset('filename',files2read);
load(dir('*.mat').name);
track_data = copy_data;

x = track_data(1).upsamp_data.task_a(1).traj_x;
y = track_data(1).upsamp_data.task_a(1).traj_y;
trial_start = find([strcmp({EEG.event.type}, 'S 27')]);
trial_end = find([strcmp({EEG.event.type}, 'S 15')]);
peaks = [EEG.event(find(strcmp({EEG.event(trial_start(1):trial_end(1)).type}, 'S 40'))).latency];

figure()
plot(x, y)
hold on 
xline(peaks)

for i = 1:size(ALLEEG, 2)

    EEG = ALLEEG(i);
    max(diff([EEG.event(strcmp({EEG.event.type}, 'S 40')).latency]) )
    EEG = pop_epoch( EEG, {  'S 40'  }, [-1  0.75], 'newname', EEG.subject, 'epochinfo', 'yes');
    EEG = eeg_checkset( EEG );
    EEG = pop_rmbase(EEG, [-1000 -500]);
    EEG.data

end



