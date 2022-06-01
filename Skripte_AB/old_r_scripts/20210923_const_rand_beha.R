
behav_dir = ("R:/AG-Beste-Studien/Emulation/04_data/behav")
markers_dir = ("R:/AG-Beste-Studien/Emulation/05_automagic/Analyzer/export")

behav_subdirs = grep(list.files(behav_dir, pattern = "*_A"), pattern = "Tracking|Traj", invert = TRUE, value = TRUE)
marker_files = list.files(markers_dir, pattern = "pp_AM.*_A_removed_practice_with_last_trigger.Markers")

# get subject ids from marker files (not all behavioral datasets have markerfiles yet)
subject_ids = vector()
for (i in 1:length(marker_files)) {
  subject_ids[i] = strsplit(marker_files[i], split = "_")[[1]][3]
}

# 27, 23, 24, 12
# extract: markers, samples
#add new columns: new samples (trigger_sample-sample_27), perc (new_sample/new_sample_12)
# behav_index*percentage --> behav frame of trigger
#subset of behavioral data for each sbj: random and constant
#calculate mean error for subsets, compare statistically

#load files & prepare export directory

for (j in 1:length(subject_ids)) {
  this_subdir = paste(behav_dir, "/", subject_ids[j], "_A", sep= "")
  this_markerfile = marker_files[j]
  setwd(this_subdir)
  trialfiles = list.files(this_subdir, pattern = "*_traj_")
  setwd(markers_dir)
  trialsamples = vector()
  markerfile = as.data.frame(read.table(this_markerfile, skip = 1, sep = ",", header = T)[, 2:3])
  rows_markers = which(markerfile[1] == " S 27" | markerfile[1] == " S 23" | markerfile[1] == " S 24" | markerfile[1] == " S 12" | markerfile[1] == " S 16")
  trial_markers = markerfile[rows_markers, ]
  rownames(trial_markers) = NULL
  trial_markers[ ,3] = NA
  trial_markers[ ,4] = NA
  colnames(trial_markers) <- c("marker", "eeg_sample", "trial_sample", "perc")
  rows_27 = which(trial_markers[,1] == " S 27")
  for (k in 1: length(rows_27)) {
    trial_markers[rows_27[k] +1, 3] = trial_markers[rows_27[k] +1, 2] - trial_markers[rows_27[k], 2]
    trial_markers[rows_27[k] +2, 3] = trial_markers[rows_27[k] +2, 2] - trial_markers[rows_27[k], 2]
    trial_markers[rows_27[k] +3, 3] = trial_markers[rows_27[k] +3, 2] - trial_markers[rows_27[k], 2]
    trial_markers[rows_27[k], 3] = 0
  }
  trial_markers = trial_markers[complete.cases(trial_markers[ , 3]), ]
  rows_12_16 = which(trial_markers[,1] == " S 12"  | trial_markers[,1] == " S 16")
  rows_27 = which(trial_markers[,1] == " S 27")
  for (l in 1: length(rows_12_16)) {
    trial_markers[rows_12_16[l], 4] = trial_markers[rows_12_16[l], 3] / trial_markers[rows_12_16[l], 3]
    trial_markers[rows_12_16[l] -1, 4] = trial_markers[rows_12_16[l] -1, 3] / trial_markers[rows_12_16[l], 3]
    trial_markers[rows_12_16[l] -2, 4] = trial_markers[rows_12_16[l] -2, 3] / trial_markers[rows_12_16[l], 3]
    trial_markers[rows_12_16[l] -3, 4] = trial_markers[rows_12_16[l] -3, 3] / trial_markers[rows_12_16[l], 3]
  }
  setwd(this_subdir)
  for (z in 1: length(trialfiles)) {
    #list_all_reduced[[j]][[z]] = list()
    trial_df = read.csv(trialfiles[z], header = T, sep = ",")
    start_const = (ceiling(trial_markers[rows_27[z]+ 1, 4] * nrow(trial_df)))
    end_const = (ceiling(trial_markers[rows_27[z]+ 2, 4] * nrow(trial_df))) 
    rand1 = trial_df[1:start_const, ]
    rand2 = trial_df[end_const: nrow(trial_df), ]
    const = trial_df[start_const:end_const, ]
    
    samples_const = 875
    
    rand1_interp = approx(rand1$error_abs, method = "linear", n = samples_const)
    const_interp = approx(const$error_abs, method = "linear", n = samples_const)
    rand2_interp = approx(rand2$error_abs, method = "linear", n = samples_const)
    
    file_const = paste(subject_ids[j], "_const_", z, ".csv", sep = "")
    file_rand1 = paste(subject_ids[j], "_rand1_", z, ".csv", sep = "")
    file_rand2 = paste(subject_ids[j], "_rand2_", z, ".csv", sep = "")
    
    file_const_interp = paste(subject_ids[j], "_const_interp", z, ".csv", sep = "")
    file_rand1_interp = paste(subject_ids[j], "_rand1_interp", z, ".csv", sep = "")
    file_rand2_interp = paste(subject_ids[j], "_rand2_interp", z, ".csv", sep = "")
    
    write.csv(const, file = file_const, row.names = F)
    write.csv(rand1, file = file_rand1, row.names = F)
    write.csv(rand2, file = file_rand2, row.names = F)
    
    write.csv(const_interp, file = file_const_interp, row.names = F)
    write.csv(rand1_interp, file = file_rand1_interp, row.names = F)
    write.csv(rand2_interp, file = file_rand2_interp, row.names = F)
#     
#     append(list_all_reduced[[j]][[z]], list(trial_df))
#     append(list_all_reduced[[j]][[z]], list(rand1))
#     # list_all_reduced[[j]][[z]][1]= trial_df
#     # list_all_reduced[[j]][[z]][2] = rand1
#     # list_all_reduced[[j]][[z]][3] = const$error_abs
#     # list_all_reduced[[j]][[z]][4] = rand2$error_abs
#     #names(list_all_reduced[[j]][[z]]) = c("whole_trial", "rand1", "const", "rand2")
    }
}


plot(trial_df$traj.x, trial_df$traj.y)