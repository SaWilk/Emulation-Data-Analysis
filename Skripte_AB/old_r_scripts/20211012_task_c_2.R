# read files created by task_c.R

#-------------------------------------------------------------------------------

path_data = "R:/AG-Beste-Studien/Emulation/04_data/behav"
setwd(path_data)

files_c = list.files(path_data, pattern = "*_trials_c.csv")

sbj_ids = vector()

for (i in 1: length(files_c)) {
  sbj_ids[i] = strsplit(files_c[i], split = "_")[[1]][1]
}

trials = 72

#-------------------------------------------------------------------------------
# read data for every subject, save means for diff velocities

df_all = matrix(nrow = length(sbj_ids), ncol = 5)

for (i in 1: length(sbj_ids)) {
  filename = paste(sbj_ids[i], "_trials_c.csv", sep = "")
  data_sbj = read.csv(filename, header = T)
  df_all[i, 1] = sbj_ids[i]
  df_all[i, 2] = mean(data_sbj$accuracy[which(data_sbj$vel == 0.7)])
  df_all[i, 3] = mean(data_sbj$accuracy[which(data_sbj$vel == 1.0)])
  df_all[i, 4] = mean(data_sbj$accuracy[which(data_sbj$vel == 1.3)])
  df_all[i, 5] = mean(data_sbj$accuracy)
}

colnames(df_all) <- c("sbj_id", "mean_slow", "mean_normal", "mean_fast", "mean_all")

setwd("R:/AG-Beste-Studien/Emulation/04_data")
write.csv(df_all, file = "means_c.csv")

df_all = as.data.frame(df_all)

# compare accuracy for different velocities

library(ggplot2)

mean_slow = data.frame(mean = as.numeric(df_all$mean_slow), target_velocity = as.factor("slow"))
mean_normal = data.frame(mean = as.numeric(df_all$mean_normal), target_velocity = as.factor("medium"))
mean_fast = data.frame(mean = as.numeric(df_all$mean_fast), target_velocity = as.factor("fast"))
means = rbind(mean_slow, mean_normal, mean_fast)

ggplot(means, aes(x = target_velocity, y = mean, color = target_velocity)) +
  geom_boxplot(notch =F, outlier.shape = NA) +
  geom_jitter(size = 1.5, width = 0.1, height = 0)+#, aes(shape = target_velocity))
  labs(y = "Mean accuracy of button press (diff in sec)")

#-------------------------------------------------------------------------------
# read data for tracking performance in task A

data_a = read.csv("behav_means.csv", header = T)

data_merge = merge(data_a, df_all, by.x = "subject", by.y = "sbj_id")
data_merge$mean_all = as.numeric(data_merge$mean_all)
data_merge$mean_slow = as.numeric(data_merge$mean_slow)
data_merge$mean_fast = as.numeric(data_merge$mean_fast)
data_merge$mean_normal = as.numeric(data_merge$mean_normal)

#-------------------------------------------------------------------------------
# add other data for tracking performance

cor.test(data_merge$wholetrial3, data_merge$mean_all)

data_add = read.csv("20210903_dataframe_means.csv", header = T)
data_merge = merge(data_add, data_merge, by.x = "subject_id", by.y = "subject")

cor.test(data_merge$mean_const, data_merge$mean_all)
cor.test(data_merge$const1, data_merge$mean_all)
cor.test(data_merge$const2, data_merge$mean_all)
cor.test(data_merge$const3, data_merge$mean_all)

plot(data_merge$mean_all)

cor.test(data_merge$mean_all, data_merge$mean_slow)

#-------------------------------------------------------------------------------
# merge with information on tracking performance

setwd(path_data)
traj_choice = read.csv("traj_choice.csv", header = T)

#data_merge = merge(data_merge, traj_choice, by.x = "subject_id", by.y = "subject")

# only 3 subjects
#need data for tracking performance on A for all available subjects first 

#read overall data for a

setwd("R:/AG-Beste-Studien/Emulation/04_data")
write.csv(data_merge, "overall_behav_data.csv")

