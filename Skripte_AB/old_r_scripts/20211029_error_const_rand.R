library(stringr)
#-------------------------------------------------------------------------------

behav_dir = ("R:/AG-Beste-Studien/Emulation/04_data/behav")
markers_dir = ("R:/AG-Beste-Studien/Emulation/05_automagic/Analyzer/export")

behav_subdirs = grep(list.files(behav_dir, pattern = "*_A"), pattern = "Tracking|Traj|*.npz|*.csv", invert = TRUE, value = TRUE)
marker_files = list.files(markers_dir, pattern = "pp_AM.*_A_Segmentation_no_practice_all.vmrk")

#-------------------------------------------------------------------------------

subject_ids = vector()
for (i in 1:length(marker_files)) {
  subject_ids[i] = strsplit(marker_files[i], split = "_")[[1]][3]
}

#-------------------------------------------------------------------------------

for (i in 1){#:length(subject_ids)){
  tryCatch({
    # setwd(behav_dir, "/", subject_ids[i], "_A", sep = "")
    setwd(markers_dir)
    markerfile = as.data.frame(read.table(marker_files[i], skip = 12, sep = ",", col.names = c("type", "trig", "sample", "length","x", "xx"), header = F, fill = T)[, 2:3])
    
    # select only those triggers that indicate start/stop of trial or const traj
    rows_markers = which(markerfile[1] == "S 27" | markerfile[1] == "S 23" | markerfile[1] == "S 24" | markerfile[1] == "S 12" | markerfile[1] == "S 16")
    trial_markers = markerfile[rows_markers, ]
    
    # reset row numbers
    rownames(trial_markers) = NULL
    
    # add third and fourth column for new samples and percentage values
    trial_markers[ ,3] = NA
    trial_markers[ ,4] = NA
    colnames(trial_markers) <- c("marker", "eeg_sample", "trial_sample", "perc")
    
    # look for trigger S27 indicating trial start to loop for trials
    rows_27 = which(trial_markers[,1] == "S 27")
    
    # check for trialnum
    if (length(rows_27) != 72) {
      print(paste("warning: trialnumber =", length(rows_27), "for subject", subject_ids[i]))
    }
    
    # if S16, check if there is S12 following, if TRUE then delete line of S12
    delete_lines_S12 = vector()
    for (ii in 1: length(rows_27)){
      if (trial_markers[rows_27[ii]+3, 1] == "S 16") {
        if (trial_markers[rows_27[ii]+4, 1] == "S 12") {
          delete_lines_S12 = c(delete_lines_S12, rows_27[ii] + 4)
        }
      }
    }
    trial_markers = trial_markers[-delete_lines_S12, ]
    # reset rownames for indexing
    rownames(trial_markers) = NULL
    # new rows_27 due to changed lines
    rows_27 = which(trial_markers[,1] == "S 27")
    
    trial_markers[ ,2] = as.numeric(trial_markers[ ,2])
    # for every S27 (trialstart), calculate new samples for following 3 triggers
    for (ii in 1: length(rows_27)) {
      tryCatch( {
      trial_markers[rows_27[ii] +1, 3] = trial_markers[rows_27[ii] +1, 2] - trial_markers[rows_27[ii], 2]
      trial_markers[rows_27[ii] +2, 3] = trial_markers[rows_27[ii] +2, 2] - trial_markers[rows_27[ii], 2]
      trial_markers[rows_27[ii] +3, 3] = trial_markers[rows_27[ii] +3, 2] - trial_markers[rows_27[ii], 2]
      trial_markers[rows_27[ii], 3] = 0
      # check: if the fourth trigger is not S27, give error message
      if (ii < length(rows_27) & trial_markers[rows_27[ii] +4, 1] != "S 27")
        print(paste("error: no S27 following trial", ii, "for subject", subject_ids[i]))
      }, error = function(x)
        print(paste("error s 27 sbj", subject_ids[i], "trial", ii))
        )
    }
    # delete missings (mostly S 12 at first line)
    # # flag rows for trial end
    trial_markers = trial_markers[complete.cases(trial_markers[ , 3]), ]
    rows_12_16 = which(trial_markers[,1] == "S 12"  | trial_markers[,1] == "S 16")
    rows_27 = which(trial_markers[,1] == "S 27")
    
    # calculate percentage
    
    for (ii in 1: length(rows_12_16)) {
      tryCatch({
      trial_markers[rows_12_16[ii], 4] = trial_markers[rows_12_16[ii], 3] / trial_markers[rows_12_16[ii], 3]
      trial_markers[rows_12_16[ii] -1, 4] = trial_markers[rows_12_16[ii] -1, 3] / trial_markers[rows_12_16[ii], 3]
      trial_markers[rows_12_16[ii] -2, 4] = trial_markers[rows_12_16[ii] -2, 3] / trial_markers[rows_12_16[ii], 3]
      trial_markers[rows_12_16[ii] -3, 4] = trial_markers[rows_12_16[ii] -3, 3] / trial_markers[rows_12_16[ii], 3]
      }, error = function(x)
        print(paste("error  s 27 sbj", subject_ids[i], "trial", ii))
        )
    }
    
    setwd(paste(behav_dir, "/", subject_ids[i], "_A", sep = ""))
    trialfiles = list.files(paste(behav_dir, "/", subject_ids[i], "_A", sep = ""), pattern = "*_traj_")
    #sort trialfiles to ensure correct trial order
    trialfiles = str_sort(trialfiles, numeric = T)
    # 
    # for (ii in 1: length(trialfiles)) {
    #   trial_df = read.csv(trialfiles[ii], header = T, sep = ",")
    #   start_const = (ceiling(trial_markers[rows_27[ii]+ 1, 4] * nrow(trial_df)))
    #   end_const = (ceiling(trial_markers[rows_27[ii]+ 2, 4] * nrow(trial_df))) 
    #   rand1 = trial_df[1:start_const, ]
    #   rand2 = trial_df[end_const: nrow(trial_df), ]
    #   const = trial_df[start_const:end_const, ]
    #   file_const = paste(subject_ids[i], "_const_", ii, ".csv", sep = "")
    #   file_rand1 = paste(subject_ids[i], "_rand1_", ii, ".csv", sep = "")
    #   file_rand2 = paste(subject_ids[i], "_rand2_", ii, ".csv", sep = "")
    #   write.csv(const, file = file_const)
    #   write.csv(rand1, file = file_rand1)
    #   write.csv(rand2, file = file_rand2)
    # }
    eeg_dur = trial_markers[which(trial_markers$marker == "S 12" | trial_markers$marker == "S 16"), 6]
    samps_psychopy = vector()
    for (ii in 1:length(trialfiles)){
      tryCatch({
        trial_df = read.csv(trialfiles[ii], header = T, sep = ",")
        samps_psychopy[ii] = nrow(trial_df)
      }, error = function(x)
        print(paste("error calculating trialdur from behavior for subject", subject_ids[i], "trial", ii))
        )
    }
    psychopy_dur = samps_psychopy/60
  }, error = function(x)
    print(paste("error sbj", subject_ids[i]))
    )
}

#delete gap array 

trial_markers$sec = trial_markers$eeg_sample/256
trial_markers$sec_trial = trial_markers$trial_sample/256
