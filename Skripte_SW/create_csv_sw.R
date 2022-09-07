#script to transform unipark and other data to .csv file
#
#
#When checking for duplicates Script needs to be run through str+shift+enter or by using "Source"
#If there are no duplicates script can be run regularly


# COMMENT: For protability, include these lines in your scripts. You cannot know 
# whether the other person has the packages you are using installed. 
if (require(readxl) == FALSE) {
  install.packages("readxl", dependencies = T)
  library(readxl)
}
if (require(rstudioapi) == FALSE) {
  install.packages("rstudioapi", dependencies = T)
  library(rstudioapi)
}
if (require(dplyr) == FALSE) {
  install.packages("dplyr", dependencies = T)
  library(dplyr)
}

#if duplicates are not filtered out, this needs to be 1. If data is without duplicates -> 0
#the duplicate rows that should be removed, need to be marked in the corresponding section
#is_data_filtered <- 1 # COMMENT: Much better to use an automatic detection of 
# duplicates in the data, e.g. using unique(). Also might be sensible to add a warning() message

#script requires one .csv from unipark and one xlsx for behavioral data

###########################

current_folder = dirname(rstudioapi::getSourceEditorContext()$path)
setwd(current_folder)

#----read data from unipark----
data_name = grep("^data_project_", dir(current_folder),value=T) 
# COMMENT: was too error-prone, changed it for now, may need adjustment in the 
# future because the default unipark output name sucks. 

data = read.csv2(data_name, sep= ";", encoding="UTF-8")
dim(data)



#----filter data----
#show duplicates
 #names(which(table(data$sub_id) > 1))
 #if(is_data_filtered == 1){
# insert rownumbers of duplicate data
#data <- data[-c(6,33,32,30,38),] # COMMENT: Style remark: Do not use row numbers 
# to reference data. Very error-prone. Instead, use something like grep to 
# identify the sub_id to make sure you always select the correct rows. 
# Also more readable. 
#COMMENT_2: grep does not differentiate between duplicate entries.
#If one entry should be removed while another one is supposed to be kept,
#indexes, multiple greps or user input is required
 #}

#check for duplicates and let user decide how duplicates are being handled
if(any(duplicated(data$sub_id))) {
  
  print("There are duplicates in the provided dataset.")
  print("Should all duplicate entries be discarded? ")
  discard <- readline(prompt="If an ID occurs twice, both entries will be removed. yes/no: ") 
  
if(tolower(discard) == "y" | tolower(discard) == "yes"){
    duplicates <- data[-duplicated(data$sub_id, fromLast = T) | -duplicated(data$sub_id),]
    data <- dplyr::setdiff(data, duplicates)
  }
  if(tolower(discard) == "n" | tolower(discard) == "no"){
    print("Would you like to investigate the duplicate entries or would you like to quit ")
    investigate <- readline(prompt="Select either investigate or quit. investigate/quit: ") 
    
    if(tolower(investigate) == "investigate"){
  
      for(i in names(which(table(data$sub_id) > 1))) {
      cat(paste("ID: \"", i, "\" exists ", length(grep(i, data$sub_id)), "times."))
        invest_indiv <- readline(prompt = "Would you like to investigate them or remove them all? investigate/remove: ")
        if(tolower(invest_indiv == "investigate")) {
          print(data[data$sub_id == i,])
          print("which entries would you like to remove?")
          print("Make sure that either a single or no entry remains.")
          remove_entries <- readline(prompt = "For example type \"1, 3, 5\" and press enter :")
          remove_entries <- as.integer(unlist(strsplit(remove_entries, ",")))
          rm_df <- data[data$sub_id == i,]
          rm_df <- rm_df[-remove_entries,]
          data <- dplyr::setdiff(data, rm_df)
        } else {
          data <- data[data$sub_id != i,]
        }

      }
    } else {
      stopifnot(investigate == "investigate")
    }
  }
 
   

}

#-------- recode inverted autism items -----------------

autism_items_idx = grep("aut", names(data))
no_is_one_point = c("aut_01", "aut_03", "aut_08", "aut_10", "aut_11", "aut_14",
                    "aut_15", "aut_17", "aut_24", "aut_25", "aut_27", "aut_28", 
                    "aut_29", "aut_30", "aut_31", "aut_32", "aut_34", "aut_36", 
                    "aut_37", "aut_38", "aut_40", "aut_44", "aut_47", "aut_48",
                    "aut_49", "aut_50")
data[,no_is_one_point] = ifelse(data[,no_is_one_point] > 2, 1, 0) # recode the "nos" as one point
data[,autism_items_idx][
  !(grep("aut", names(data), value=T) %in%  no_is_one_point )] = ifelse(data[,autism_items_idx][
    !(grep("aut", names(data), value=T) %in%  no_is_one_point )] < 3, 1, 0) # recode the "yes" as one point

# recode gaming items

gaming_items_idx = grep("gam", names(data))
data[,gaming_items_idx] = data[,gaming_items_idx] -1 # never is 0 now, 10+ is 5


#----define video game category----

total_df <- data.frame()

for(i in 1:dim(data)[1]){
  latest_subj_idx = i
  data_latest_subj = data[latest_subj_idx,]
  
  attach(data_latest_subj, warn.conflicts = F)
  names(data_latest_subj)
  
  
  # --- Defining video game player categories ---------------
  try({
    game_items_idx_pres = grep("gam_0[0-8]", names(data))
    game_items_idx_pas = grep("gam_09|gam_1[0-9]", names(data))
    
    # past year
    nvgp_crit = (all(c(gam_01, gam_02, gam_03, gam_04) < 2) & # action video games
                   all(c(gam_05, gam_06, gam_07, gam_08) < 3) & # other games
                   sum(c(gam_01, gam_02, gam_03, gam_04, gam_05, gam_06, gam_07, gam_08)) < 5 & # overall gaming
                   # before the past year
                   all(c(gam_09, gam_10, gam_11, gam_12) < 2) & # action video games
                   all(c(gam_13, gam_14, gam_15, gam_16) < 3) & # other games
                   sum(c(gam_09, gam_10, gam_11, gam_12, gam_13, gam_14, gam_15, gam_16)) < 5) # overall gaming
    
    
    # OPTION 1: current high action game players
    avgp_crit = (((gam_01 | gam_02) > 3 & # shooter or action rpg
                    all(c(gam_05, gam_06, gam_07, gam_08) < 3)) | # all except for RTS and Sport
                   # OPTION 2: current medium action game players, past heavy gamers
                   ((gam_09 | gam_10) > 3 & # shooter or action rpg
                      # before the past year
                      all(c(gam_05, gam_06, gam_07, gam_08) < 3)) | # all except for RTS and Sport
                   # OPTION 3: current medium action gamers but high sports gamers
                   ((gam_01 | gam_02) > 3 & 
                      gam_03 > 3 & # sport games
                      all(c(gam_05, gam_06, gam_07, gam_08) < 2)) | # all genres except for sport 
                   # OPTION 4: current medium action players and heavy RTS/MOBA players
                   ((gam_01 | gam_02) > 3 & 
                      gam_04 > 3 & # RTS/MOBA
                      all(c(gam_05, gam_06, gam_07, gam_08) < 2)))  # all genres except for sport 
    
    
    # hourly means of all present game types
    hourly_means_pres = (0.5*sum(data_latest_subj[,game_items_idx_pres] == 1) +
                           2*sum(data_latest_subj[,game_items_idx_pres] == 2) + 
                           4*sum(data_latest_subj[,game_items_idx_pres] == 3) + 
                           7.5*sum(data_latest_subj[,game_items_idx_pres] == 4) + 
                           15*sum(data_latest_subj[,game_items_idx_pres] == 5))
    
    hourly_means_pas = (0.5*sum(data_latest_subj[,game_items_idx_pas] == 1) +
                          2*sum(data_latest_subj[,game_items_idx_pas] == 2) + 
                          4*sum(data_latest_subj[,game_items_idx_pas] == 3) + 
                          7.5*sum(data_latest_subj[,game_items_idx_pas] == 4) + 
                          15*sum(data_latest_subj[,game_items_idx_pas] == 5))
    
    # individuals in between action game players and non video game players
    low_tweener_crit = ((gam_01 & gam_02) < 2 & # shooter or action rpg
                          ((gam_03 < 3 & gam_04 < 2) | (gam_03 < 2 & gam_04 < 3)) & # sports and rts
                          (gam_09 & gam_10) < 3 & # past action games
                          (gam_11 & gam_12) < 4 & # past sport and rts
                          !nvgp_crit & # not a nvgp
                          # hourly means of all game types present
                          ((1 > hourly_means_pres & hourly_means_pres <= 10 ) |
                             # hourly means of all game types present
                             (hourly_means_pres <= 1  & 
                                # hourly means of all game types past
                                hourly_means_pas >= 5 )) &
                          # non-expert experience in any game
                          ((gam_01 & gam_02 & gam_03 & gam_04 & gam_05 & gam_06 & gam_07 & gam_08) < 4)) # all present games
    
    # to label all the poor people that do not fall into any category...
    other_gamers_crit = !avgp_crit & !nvgp_crit & !low_tweener_crit
    
    
    #---- split motor imagery scales -----------------------------
    
    evi_items_idx = grep("moi_0[1-9]|moi_1[0-2]", names(data)) # external visuell
    ivi_items_idx = grep("moi_1[3-9]|moi_2[0-4]", names(data)) # internal visuell
    kin_items_idx = grep("moi_2[5-9]|moi_3[0-6]", names(data)) # kinästhetisch
    
    
    #------- merge subjects together into a df ------------------
    
    subj_i <- data.frame(sub_id = data_latest_subj$sub_id, nvgp_crit, avgp_crit, low_tweener_crit, other_gamers_crit, 
                         evi_items_idx = sum(data_latest_subj[,evi_items_idx]), 
                         ivi_items_idx = sum(data_latest_subj[,ivi_items_idx]),
                         kin_items_idx = sum(data_latest_subj[,kin_items_idx]),
                         handedness = sum(data_latest_subj[,grep("hand", names(data))]),
                         autism = sum(data_latest_subj[,grep("aut", names(data))]))
    ifelse(i == 1, total_df <- subj_i, total_df[i,] <- subj_i)
  })
}
#----merge video game categories and imagery scales with original df----

data_up_merged <- merge(total_df, data, by = "sub_id")


#----read behavioral data----

# COMMENT: note: xlsx file must be closed, otherwise throws an error
data_name_behav <- grep("Verhaltensdaten", dir(current_folder),value=T)

data__behav_s1 <- read_xlsx(data_name_behav, sheet = 1)
data__behav_s2 <- read_xlsx(data_name_behav, sheet = 2)
data__behav_s3 <- read_xlsx(data_name_behav, sheet = 3)

data_total <- merge(merge(data__behav_s1, data__behav_s2, by = "VP_ID"), data__behav_s3, by = "VP_ID")

#----read coding list ----

# COMMENT: I adjusted the Verhaltensdaten file to include a gorup variable, pleas
# update the code here. 
#data_name_list <- grep("list", dir(current_folder),value=T)
#data_list <- read_xlsx(data_name_list, sheet = 1)
#names(data_list)[names(data_list) == "ID"] <- "sub_id" 
#data_list <- data_list[, c("sub_id", "Gruppe")] # COMMENT: style note - we should 
# name variables consistently in snake_case and in English. Increases efficiency 
# when coding.


#----merge data to one df and save as .csv---- 
names(data_total)[names(data_total) == "VP_ID"] <- "sub_id" 
data_merged <- merge(data_total, data_up_merged, by = "sub_id")
#data_merged <- merge(data_total_1, data_up_merged, by = "sub_id")

data_merged$lang_ger[data_merged$lang_ger != 1] <- "not_german"
data_merged$lang_ger[data_merged$lang_ger == 1] <- "german"

col_names = c("sub_id", "age", "edu", "lang_ger", "sex", "group","job", "sport_01", "sport_02", "athlete_01",
              "athlete_02", "athlete_03", "athlete_04",
              "block_span_forward", "block_span_backward", "block_sum_forward", "block_sum_backward",
              "corsi_points_forward", "corsi_points_backward", "cat_strat", "cat_const", "APM_sumscore",
              "digit_span_forward", "digit_span_backward", "digit_sum_forward", "digit_sum_backward",
              "task_order", "corsi_forward_percentile",	"digit_forward_percentile",	"raven_sample_percentile",
              "digit_backward_percentile",
              "nvgp_crit", "avgp_crit", "low_tweener_crit", "other_gamers_crit", "evi_items_idx", "ivi_items_idx", 
              "kin_items_idx", "handedness", "autism")

# debugging: check if all col_names are inside data_merged colnames
#col_names[!c(col_names %in% colnames(data_merged))] 
# colnames(data_merged) 
# COMMENT: I renamed the variable names in the xlsx file because they were awful.
# and yeah, I know I put them in the sample excel file like that... ^^

recoded_data <- data_merged[,col_names]
recoded_data[recoded_data == -66] <- NA
recoded_data[recoded_data == -77] <- NA


data_csv <- recoded_data

write.csv(data_csv, "combined_data.csv", row.names = F) # COMMENT: This would have
# caused the script to fail if I had not changed the original definition of what
# is the input to this script. In general it is good practice to use a separate
# folder for the outputs and inputs to a script, see dir.create()
#Comment_2: it would only cause a crash when running it twice in the same directory.

# FINAL COMMENTS:
# TODOS: Assign NA values to all numbers that code for "missing response" 
# (i.e. all the weird ones, like -77 and -66)
# Really Weird: we have lots of missing values in the column 'sex'... No idea 
# why so far... probably something about Unipark itself. Am on it. 
