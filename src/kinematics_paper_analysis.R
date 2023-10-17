# Kinematics paper analysis script 
##### Instructions ####
# Please run this script prior to knitting the kinematics RMD file in order 
# to reproduce the analyses in the paper. 
#####

# Loading required libraries:
library(dplyr)
library(tidyr)
library(forcats)

# Load scripts with necessary functions:
source("./src/kinematic_fxns.R")
source("./src/ethogram_fxns.R")

# Load kinematic data:
frame_rate <- 29.97

track_list <- load_event_list()
all_dat <- assemble_kin_data(frame_rate)
all_dat <- remove_prompted(all_dat)

tracks <- pivot_kinematic_data(all_dat)

# Importing list of runs to analyze
task_list <- load_task_list()

# Loads and shapes task event data
event_data <- load_ethograms(task_list, "individual")

# Lines up ethogram and kinematic tracking time scales
points_2E1H_info <- correct_ethogram_times(event_data$points_2E1H, frame_rate)
points_2E1H_info <- find_behavior_times(points_2E1H_info, all_dat)
points_Amm_info <- correct_ethogram_times(event_data$points_Amm, frame_rate)
points_Amm_info <- find_behavior_times(points_Amm_info, all_dat)
#states_2E1H_info_starts <- correct_ethogram_times(event_data$states_2E1H[
#  event_data$states_2E1H$status == "start", ], frame_rate)
#states_2E1H_info_starts <- find_behavior_times(states_2E1H_info_starts, all_dat)
#states_2E1H_info_stops <- correct_ethogram_times(event_data$states_2E1H[
#  event_data$states_2E1H$status == "stop", ], frame_rate)
#states_2E1H_info_stops <- find_behavior_times(states_2E1H_info_stops, all_dat)
#states_2E1H_info <- rbind(states_2E1H_info_starts, states_2E1H_info_stops)
states_2E1H_info <- correct_ethogram_times(event_data$states_2E1H, frame_rate)
states_2E1H_info <- find_behavior_times(states_2E1H_info, all_dat)
states_Amm_info <- correct_ethogram_times(event_data$states_Amm, frame_rate)
states_Amm_info <- find_behavior_times(states_Amm_info, all_dat)

event_times <- rbind(event_data$time_budget_2E1H, event_data$time_budget_Amm)
event_times <- event_times[event_times$trained == T,]

point_counts <- rbind(event_data$point_counts_2E1H, event_data$point_counts_Amm)
point_counts <- point_counts[point_counts$trained==T,]
# Loads and shapes task event data
task_data <- load_ethograms(task_list, "tasks")

tasks_2E1H_info <- correct_ethogram_times(task_data$tasks_2E1H, frame_rate)
tasks_2E1H_info <- find_behavior_times(tasks_2E1H_info, all_dat)

tasks_Amm_info <- correct_ethogram_times(task_data$tasks_Amm, frame_rate)
tasks_Amm_info <- find_behavior_times(tasks_Amm_info, all_dat)

tasks_info <- rbind(tasks_2E1H_info, tasks_Amm_info)

event_data <- compare_task_point(task_list, task_data, event_data)

state_assoc_counts <- summarize_linked_behaviors(event_data, 
                                                 task_list$Dog, "states_")
point_assoc_counts <- summarize_linked_behaviors(event_data, 
                                                 task_list$Dog, "points_")


tracks <- assign_tasks_tracks(tracks, tasks_info)

tracks$task <- factor(tracks$task, levels=c("t", "c", "o", "a"))

# Restricting data set to only Trial 3 (trained)
tracks <- tracks[tracks$trained == T,]


# Saving final kinematics track data to results folder
dir.create("./results/csv-files", showWarnings = FALSE)
write.csv(x = tracks, 
          file = paste0("./results/csv-files/kinematics-final-tracks_", 
                        Sys.Date(), ".csv"),
          row.names = FALSE)

trial_holes <- all_dat[["T3D1R1"]]$holes
write.csv(x = trial_holes, file = paste0("./results/csv-files/trial_holes_",
                                         Sys.Date(),".csv"),
          row.names = FALSE)

# Calculating distances away from the hot hole for casting and on-odor/detailing
e1h_tracks_casting <- calc_hot_dist(na.omit(tracks[tracks$point == "nose" & tracks$run == 1 & tracks$task == "c",]), 
                                    all_dat[["T3D1R1"]]$holes[6,])
e1h_tracks_onodor <- calc_hot_dist(na.omit(tracks[tracks$point == "nose" & tracks$run == 1 & tracks$task == "o",]), 
                                   all_dat[["T3D1R1"]]$holes[6,])
amm_tracks_casting <- calc_hot_dist(na.omit(tracks[tracks$point == "nose" & tracks$run == 2 & tracks$task == "c",]),
                                    all_dat[[1]]$holes[9,])
amm_tracks_onodor <-  calc_hot_dist(na.omit(tracks[tracks$point == "nose" & tracks$run == 2 & tracks$task == "o",]),
                                    all_dat[[1]]$holes[9,])

mean_hot_dist <- data.frame("chemical" = rep("Ammonia", nrow(amm_tracks_casting)), 
                            "task" = rep("casting", nrow(amm_tracks_casting)), 
                            "dog" = amm_tracks_casting$dog,
                            "value" = amm_tracks_casting$hot_dist)
mean_hot_dist <- rbind(mean_hot_dist, 
                       data.frame("chemical" = rep("Ammonia", nrow(amm_tracks_onodor)), 
                                  "task" = rep("detailing", nrow(amm_tracks_onodor)), 
                                  "dog" = amm_tracks_onodor$dog,
                                  "value" = amm_tracks_onodor$hot_dist), 
                       data.frame("chemical" = rep("2E1H", nrow(e1h_tracks_casting)), 
                                  "task" = rep("casting", nrow(e1h_tracks_casting)), 
                                  "dog" = e1h_tracks_casting$dog,
                                  "value" = e1h_tracks_casting$hot_dist),
                       data.frame("chemical" = rep("2E1H", nrow(e1h_tracks_onodor)), 
                                  "task" = rep("detailing", nrow(e1h_tracks_onodor)), 
                                  "dog" = e1h_tracks_onodor$dog,
                                  "value" = e1h_tracks_onodor$hot_dist))

dog_mean_hot_dist <- mean_hot_dist %>% 
  group_by(dog, task, chemical) %>%
  summarise(value = mean(value, na.rm = T))

write.csv(x = dog_mean_hot_dist, file = paste0("./results/csv-files/dog-mean-hot-dist_", 
                                               Sys.Date(), ".csv"),
          row.names = FALSE)

# Calculate total times in each task for each dog: 
task_times <- rbind(task_data[["tasks_2E1H"]], task_data[["tasks_Amm"]])
task_times <- task_times[task_times$trial == "3",]

dog_mean_total_task_time <- task_times %>% 
  group_by(dog, behavior, target) %>%
  summarise(value = sum(elapsed_time, na.rm = T))

colnames(dog_mean_total_task_time)<- c("dog", "task", "chemical", "value")

write.csv(x = dog_mean_total_task_time, file = paste0("./results/csv-files/dog-mean-total-task-time_",
                                                      Sys.Date(), ".csv"),
          row.names = FALSE)


dog_mean_event_times <- event_times[event_times$trained == T,] %>%
  group_by(dog, behavior, target) %>%
  summarise(value = sum(total_times, na.rm = T))
colnames(dog_mean_event_times)<- c("dog", "task",  "chemical", "value")

write.csv(x = dog_mean_event_times, file = paste0("./results/csv-files/dog-mean-event-times_",
                                                      Sys.Date(), ".csv"),
          row.names = FALSE)

dog_point_counts <- point_counts %>%
  group_by(dog, behavior, target) %>%
  summarise(value = sum(counts, na.rm = T))
colnames(dog_point_counts)<- c("dog", "task", "chemical", "value")

write.csv(x = dog_point_counts, file = paste0("./results/csv-files/dog-point-counts_",
                                                  Sys.Date(), ".csv"),
          row.names = FALSE)

alerts <- rbind(task_data$counts_2E1H, task_data$counts_Amm)
alerts <- alerts[alerts$behavior == "a" & alerts$trained == T, ]
alerts$success <- ifelse(alerts$counts > 0, T, F)
alerts <- alerts %>% mutate(dog = forcats::fct_reorder(dog, as.numeric(as.character(dog))))

write.csv(x = alerts, file = paste0("./results/csv-files/dog-alerts_",
                                              Sys.Date(), ".csv"),
          row.names = FALSE)
