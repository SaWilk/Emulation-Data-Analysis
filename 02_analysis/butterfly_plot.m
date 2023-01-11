%% Make Butterfly Plots

% Author Saskia Wilken

% Creation Date 24.11.2022


%% empty everything

format compact
format long G
clear
clc
close all


%% Define Paths

% file location
file_path = matlab.desktop.editor.getActiveFilename;
file_path = fileparts(file_path);
addpath(file_path);

% get data paths for parent dirs
filepath_parts = strsplit(file_path, filesep);
parent_dir = filepath_parts(1:end-1);
parent_dir = strjoin(parent_dir, filesep);
parent_dir_2 = filepath_parts(1:end-2);
parent_dir_2 = strjoin(parent_dir_2, filesep);
output_dir = strjoin([parent_dir_2, "Emulation-Data-Output"], filesep);
track_data_dir = strjoin([output_dir, "03_parallelize_with_traj"], filesep);
mean_matrices_path = strjoin([output_dir, 'mean_matrices'], filesep);
mkdir(mean_matrices_path);
mean_matrices_peaks_epoched_path = strjoin([mean_matrices_path, 'peaks'], filesep);
mkdir(mean_matrices_peaks_epoched_path);
peak_erp_plots_dir = strjoin([parent_dir, "plots", "erp_plots", "peaks"], filesep);
mkdir(peak_erp_plots_dir);
% output
csd_dir = strjoin([output_dir, "csd_transform"], filesep);
% input
epochs_plus_error_dir = strjoin([output_dir, "07_epochs_with_extra_fields"], filesep);
% cluster based permutation test
neighbors_dir = strjoin([parent_dir_2, "Emulation-Data-Input", "EEG_files"], filesep);


%% Load CDS data

eeglab;

% long epochs
cd(csd_dir);
%list all *.set files in inputpath
file_names = dir('*.set');
%concatenate into one cell array
files2read = {file_names.name};
files2read = files2read(~contains(files2read, 'KMY6K')); % remove subject who 
% did not properly track during task B
% load eeg data epoched around peaks
for idx = 1:length(files2read)
    ALLEEG2(idx) = pop_loadset('filename',files2read{idx});
end

ALLEEG = ALLEEG2;
clear ALLEEG2


%% Prepare for Plotting

cd(mean_matrices_peaks_epoched_path)

% set(0,'DefaultTextFontName','Arial'); % change default font, if wanted

load('z_mean_struct_csd.mat')
load('mean_struct_csd.mat')
load("condition_indices_csd.mat")
mean_struct = renameStructField(mean_struct, 'diff_occ', 'diff_occ_vis');
z_mean_struct = renameStructField(z_mean_struct, 'diff_occ', 'diff_occ_vis');
% remove everything reminding of random2
rm_fields = {'rand2', 'diff_rand1_rand2', 'diff_const_rand2', 'num_epochs_rand2'};

for field = 1:length(rm_fields)
    if isfield(z_mean_struct, rm_fields{field})
        z_mean_struct = rmfield(z_mean_struct, rm_fields{field});
    end
    if isfield(mean_struct, rm_fields{field})
        mean_struct = rmfield(mean_struct, rm_fields{field});
    end
end

% cond_labs = {'all', 'occ', 'vis', 'rand1', 'const', 'rand2', 'diff_vis_occ', ...
%     'diff_const_rand1', 'diff_const_rand2', 'diff_rand1_rand2'};
cond_labs = {'all', 'occ', 'vis', 'rand1', 'const', 'diff_occ_vis', ...
    'diff_const_rand1'};

avg_peak_dist = 300;
% average across subjects
for lab = 1:length(cond_labs)
    z_epoch_mean.(cond_labs{lab}) = mean(z_mean_struct.(cond_labs{lab}),3);
    epoch_mean.(cond_labs{lab}) = mean(mean_struct.(cond_labs{lab}),3);

end


%% Plot code

cd("C:\wilken\Emulation-Data-Analysis-Journal\Butterfly Plots")

figure()
subplot(2, 1, 1)
plot(mean_struct.time_vec, epoch_mean.occ.')
title("Occluded")
xline(0)
set(gca,'fontsize', 12, 'color', 'w', 'linewidth', 1);
xlabel("time (ms)")
ylabel("amplitude (µV/mm²)")

subplot(2, 1, 2)
plot(mean_struct.time_vec, epoch_mean.vis.')
title("Visible")
xline(0)

sgtitle("CSD transformed mean across epochs and subjects")
set(gca,'fontsize', 12, 'color', 'w', 'linewidth', 1);
set(gcf, 'PaperUnits','inches','PaperPosition',[0 0 15 5])
xlabel("time (ms)")
ylabel("amplitude (µV/mm²)")
filename = strjoin(["Butterfly_Plot_of_Occluded-Visible.jpg"],"")
print('-djpeg', filename, '-r600');


figure()
subplot(2, 1, 1)
plot(mean_struct.time_vec, epoch_mean.const.')
title("Constant")
xline(0)
set(gca,'fontsize', 12, 'color', 'w', 'linewidth', 1);
xlabel("time (ms)")
ylabel("amplitude (µV/mm²)")

subplot(2, 1, 2)
plot(mean_struct.time_vec, epoch_mean.rand1.')
title("Random1")
xline(0)

sgtitle("CSD transformed mean across epochs and subjects")
set(gca,'fontsize', 12, 'color', 'w', 'linewidth', 1);
set(gcf, 'PaperUnits','inches','PaperPosition',[0 0 15 5])
xlabel("time (ms)")
ylabel("amplitude (µV/mm²)")
filename = strjoin(["Butterfly_Plot_of_Constant-Random.jpg"],"")
print('-djpeg', filename, '-r600');





