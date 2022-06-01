
behav_dir = ("R:/AG-Beste-Studien/Emulation/04_data/behav")
markers_dir = ("R:/AG-Beste-Studien/Emulation/05_automagic/Analyzer/export")

behav_subdirs = grep(list.files(behav_dir, pattern = "*_A"), pattern = "Tracking|Traj", invert = TRUE, value = TRUE)
marker_files = list.files(markers_dir, pattern = "pp_AM.*_A_removed_practice_with_last_trigger.Markers")

# get subject ids from marker files (not all behavioral datasets have markerfiles)
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

list_all_reduced = list()
for (j in 1) {
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
  for (k in 1: length(which(trial_markers[,1] == " S 27"))) {
    trial_markers[which(trial_markers[,1] == " S 27")[k] +1, 3] = trial_markers[which(trial_markers[,1] == " S 27")[k] +1, 2] - trial_markers[which(trial_markers[,1] == " S 27")[k], 2]
    trial_markers[which(trial_markers[,1] == " S 27")[k] +2, 3] = trial_markers[which(trial_markers[,1] == " S 27")[k] +2, 2] - trial_markers[which(trial_markers[,1] == " S 27")[k], 2]
    trial_markers[which(trial_markers[,1] == " S 27")[k] +3, 3] = trial_markers[which(trial_markers[,1] == " S 27")[k] +3, 2] - trial_markers[which(trial_markers[,1] == " S 27")[k], 2]
    trial_markers[which(trial_markers[,1] == " S 27")[k], 3] = 0
  }
  trial_markers = trial_markers[complete.cases(trial_markers[ , 3]), ]
  
  for (l in 1: length(which(trial_markers[,1] == " S 12"  | trial_markers[,1] == " S 16"))) {
    trial_markers[which(trial_markers[,1] == " S 12"  | trial_markers[,1] == " S 16")[l], 4] = trial_markers[which(trial_markers[,1] == " S 12" | trial_markers[,1] == " S 16")[l], 3] / trial_markers[which(trial_markers[,1] == " S 12" | trial_markers[,1] == " S 16")[l], 3]
    trial_markers[which(trial_markers[,1] == " S 12"  | trial_markers[,1] == " S 16")[l] -1, 4] = trial_markers[which(trial_markers[,1] == " S 12" | trial_markers[,1] == " S 16")[l] -1, 3] / trial_markers[which(trial_markers[,1] == " S 12" | trial_markers[,1] == " S 16")[l], 3]
    trial_markers[which(trial_markers[,1] == " S 12"  | trial_markers[,1] == " S 16")[l] -2, 4] = trial_markers[which(trial_markers[,1] == " S 12" | trial_markers[,1] == " S 16")[l] -2, 3] / trial_markers[which(trial_markers[,1] == " S 12" | trial_markers[,1] == " S 16")[l], 3]
    trial_markers[which(trial_markers[,1] == " S 12"  | trial_markers[,1] == " S 16")[l] -3, 4] = trial_markers[which(trial_markers[,1] == " S 12" | trial_markers[,1] == " S 16")[l] -3, 3] / trial_markers[which(trial_markers[,1] == " S 12" | trial_markers[,1] == " S 16")[l], 3]
  }
  setwd(this_subdir)
  j = 1
  list_all_reduced[[j]] = list()
  names(list_all_reduced)[j] = behav_subdirs[j]
  for (z in 1: length(trialfiles)) {
    list_all_reduced[[j]][[z]] = list()
    trial_df = read.csv(trialfiles[z], header = T, sep = ",")
    start_const = (ceiling(trial_markers[which(trial_markers[, 1] == " S 27")[1]+ 1, 4] * nrow(trial_df)))
    end_const = (ceiling(trial_markers[which(trial_markers[, 1] == " S 27")[1]+ 2, 4] * nrow(trial_df))) 
    rand1 = trial_df[1:start_const, ]
    const = trial_df[start_const:end_const, ]
    rand2 = trial_df[end_const: nrow(trial_df), ]
    list_all_reduced[[j]][[z]][1]= as.data.frame(trial_df)
    list_all_reduced[[j]][[z]][2] = rand1
    list_all_reduced[[j]][[z]][3] = const
    list_all_reduced[[j]][[z]][4] = rand2
    names(list_all_reduced[[j]][[z]]) = c("whole_trial", "rand1", "const", "rand2")
    }
}

#problem: nicht alle variablen werden in die liste überführt --> warum?