% Load in trajectory data

file_path = matlab.desktop.editor.getActiveFilename;
file_path = fileparts(file_path);
data_path = strcat(file_path, "\00_npz_files"); % get root path to data folder
subject_paths = genpath(data_path); % generate subfolder paths
subject_paths = strsplit(subject_paths, ";"); % split the string vector at ;
length(subject_paths); % get num of subfolders
all_data = struct(); % create empty structure
subject_paths(1) = []; % remove top level entry


%% Load data

for p = 1:length(subject_paths)
    if ~isempty(subject_paths{p})
        cd(subject_paths{p});
        file_names_curr_subj = {dir('*.csv').name}
        for s = 1:length(file_names_curr_subj)
            all_data.tracking.subjects{s} = readtable(file_names_curr_subj{s})
        end
    end
end


 all_data.tracking{s} = 1





 