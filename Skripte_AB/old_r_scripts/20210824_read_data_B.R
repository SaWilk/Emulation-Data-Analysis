#############################
# read data emulation study #
#############################

#TASK B
#-----------------------------------------
#preparations

# install package for python extension
if (require(reticulate) == FALSE) {
  install.packages("reticulate")
  library(reticulate)
}

#set data directory
this_dir = dirname(rstudioapi::getSourceEditorContext()$path)

# select emulation task data analysis main directory
parent_dir = strsplit(this_dir,"/")
parent_dir = parent_dir[[1]][1:(length(parent_dir[[1]])-2)]
parent_dir = paste(parent_dir,collapse ="/" )

datadir = paste(parent_dir, "00_npz_files", sep="/")
setwd(datadir)

#import python 
np <- import("numpy")

#read files + number of subjects
n_files_all = length(list.files(datadir, pattern = "Tracking"))

n_files_B = length(list.files(datadir, pattern = "Tracking.*_B-Trials.npz"))
files_B =  list.files(datadir, pattern = "Tracking.*_B-Trials.npz")
# one file for B is missing, was not saved 

#-----------------------------------------
#extract and enumerate subject ids
#cave: subject numbers are not in real order
#todo: extract from excel file 
subj_ids = matrix(NA, length(files_B), 2) 
for (i in 1:n_files_B) {
  subj_ids[i ,1] = substr(files_B[i], 10, 14)
  subj_ids[i, 2] = i
}

#-----------------------------------------
#set trialnum
trialnum = 72

#-----------------------------------------
#load data from npz, save as list

list_b = list()

for (i in 1:n_files_B) {
  filename = paste("Tracking-", subj_ids[i, 1], "_B-Trials.npz", sep = "")
  data <- np$load(filename, allow_pickle = TRUE)
  data_list = data$f[["arr_0"]][[1]]
  list_b[[i]] = data_list
}
rm(data_list)


#-----------------------------------------
#calculate tracking error for every trial and append to list

for (i in 1: n_files_B) {
  vector_rmsd = vector()
  for (j in 1:trialnum) {
    vector_rmsd[j] = sqrt(sum(((list_b[[i]]$Trajectory[[j]][ ,2] - list_b[[i]]$Pursuit[[j]][ ,2])^2))/trialnum)
    list_b[[i]]$rmsd = vector_rmsd
  }
}

#-----------------------------------------
#calculate error curve and add to list
#error curve = difference of traj and purs

for (i in 1:n_files_B) {
  for (j in 1:trialnum) {
    vector_error = vector()
    for (k in 1: length(list_b[[i]]$Trajectory[[j]][,2])) {
      vector_error[k] = (list_b[[i]]$Trajectory[[j]][ ,2][k] - list_b[[i]]$Pursuit[[j]][ ,2][k])
      list_b[[i]]$error[[j]] = vector_error
    }
  }
}

#abs error curve

for (i in 1:n_files_B) {
  for (j in 1:trialnum) {
    vector_error_abs = vector()
    for (k in 1: length(list_b[[i]]$Trajectory[[j]][,2])) {
      vector_error_abs[k] = abs(list_b[[i]]$Trajectory[[j]][ ,2][k] - list_b[[i]]$Pursuit[[j]][ ,2][k])
      list_b[[i]]$error_abs[[j]] = vector_error_abs
    }
  }
}

#-----------------------------------------
#median split based on rmsd of tracking error

for (i in 1:n_files_B) {
  median_split = vector()
  median_rmsd = median(list_b[[i]]$rmsd)
  for (j in 1:trialnum) {
    median_split[j] = ifelse(list_b[[i]]$rmsd[j] >= median_rmsd, 1, 0)
    list_b[[i]]$median_split = median_split
  }
} 

#-----------------------------------------
#data export

if (getwd() != datadir) {
  setwd(datadir)
}

#date = Sys.Date()
#tasks = c("A", "B", "C")

for (i in 1:n_files_B) {
  folder_name = paste(subj_ids[i], "_B", sep = "")
  file_name = paste(subj_ids[i], "_B_trialdata.csv", sep = "")
  dir.create(folder_name)
  new_path = paste(datadir, "/", folder_name, sep = "")
  setwd(new_path)
  data_export = as.data.frame(cbind(list_b[[i]]$Trial, list_b[[i]]$TrialDuration, list_b[[i]]$Direction, list_b[[i]]$rmsd, list_b[[i]]$median_split))
  colnames(data_export) = c("Trial", "Dur", "Dir", "rmsd", "median_split")
  write.csv(data_export, file = file_name, row.names = FALSE)
  for (j in 1:trialnum) {
    file_name_2 = paste(subj_ids[i], "_traj_purs_", j, ".csv", sep = "")
    data_export_2 = as.data.frame(cbind(list_b[[i]]$Trajectory[j][[1]], list_b[[i]]$Pursuit[j][[1]], list_b[[i]]$error[j][[1]], list_b[[i]]$error_abs[j][[1]]))
    colnames(data_export_2) <- c("traj-x", "traj-y", "purs-x", "purs-y", "error", "error_abs")
    write.csv(data_export_2, file = file_name_2, row.names = FALSE)
  }
  setwd(datadir)
  rm(data_export)
  rm(data_export_2)
}
