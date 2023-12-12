# Ethogram functions

# Load the ethogram task list
load_task_list <- function(opt = "default"){
  if(opt == "default"){
    task_list <- read.csv("./data/ethogram/task-list.csv")
    task_list <- task_list[task_list$Dog != 21,] # Removing dog 21
    task_list$Dog <- factor(task_list$Dog)
  }else if(opt == "kinematics-paper"){
    task_list <- read.csv("./data/ethogram/task-list-kinematics-paper.csv")
    task_list <- task_list[,-grep("Notes", colnames(task_list))]
    task_list <- na.omit(task_list) # Removing rows with NA (not completed)
    task_list <- task_list[task_list$Dog != 10,] # Removing Dog 10
    #task_list <- task_list[task_list$Dog != 21,] # Removing dog 21
  }else{
    stop("Please choose a valid option.")
  }
  task_list$Dog <- factor(task_list$Dog)
  return(task_list)
}


# Head function that loads and reformats the data
load_ethograms <- function(task_list, opt){
  require(stringr)
  if(opt=="tasks"){
    
    # Initializing list to store data
    event_data <- list()
    
    # Loops over each dog in the task list and calculates: 
    for(k in levels(task_list$Dog)){
      dog_tasks <- task_list[task_list$Dog == k, ]
      #trial, dog, run, trained
      print(paste("Dog:",k,", Trained:", dog_tasks$Trained))
      dog_before_2E1H <- load_ethogram_tasks(dog_tasks$Trial[1], k, 1, 
                                             dog_tasks$Trained[1])
      dog_before_Amm <- load_ethogram_tasks(dog_tasks$Trial[2], k, 2, 
                                            dog_tasks$Trained[2])
      dog_after_2E1H <- load_ethogram_tasks(dog_tasks$Trial[3], k, 1, 
                                            dog_tasks$Trained[3])
      dog_after_Amm <- load_ethogram_tasks(dog_tasks$Trial[4], k, 2, 
                                           dog_tasks$Trained[4])
      if(is.null(dog_before_2E1H) | is.null(dog_after_2E1H) | is.null(dog_before_Amm) | 
         is.null(dog_before_Amm)){
        time_changes_2E1H <- NULL
        time_changes_Amm <- NULL
        count_changes_2E1H <- NULL
        count_changes_Amm <- NULL
      }else {
        time_changes_2E1H <- calc_time_differences(dog_before_2E1H$time_budget,
                                                   dog_after_2E1H$time_budget)
        time_changes_Amm <- calc_time_differences(dog_before_Amm$time_budget,
                                                  dog_after_Amm$time_budget)
        count_changes_2E1H <- calc_count_differences(dog_before_2E1H$event_counts,
                                                     dog_after_2E1H$event_counts)
        count_changes_Amm <- calc_count_differences(dog_before_Amm$event_counts,
                                                    dog_after_Amm$event_counts)
      }
      if(k == "1"){
        event_data$tasks_2E1H <- rbind(dog_before_2E1H$tasks, 
                                       dog_after_2E1H$tasks)
        event_data$tasks_Amm <- rbind(dog_before_Amm$tasks, 
                                      dog_after_Amm$tasks)
        event_data$time_budget_2E1H <- rbind(dog_before_2E1H$time_budget, 
                                             dog_after_2E1H$time_budget)
        event_data$time_budget_Amm <- rbind(dog_before_Amm$time_budget, 
                                            dog_after_Amm$time_budget)
        event_data$counts_2E1H <- rbind(dog_before_2E1H$event_counts, 
                                        dog_after_2E1H$event_counts)
        event_data$counts_Amm <- rbind(dog_before_Amm$event_counts, 
                                       dog_after_Amm$event_counts)
        event_data$time_changes_2E1H <- time_changes_2E1H
        event_data$time_changes_Amm <- time_changes_Amm
        event_data$count_changes_2E1H <- count_changes_2E1H
        event_data$count_changes_Amm <- count_changes_Amm
      }else{
        event_data$tasks_2E1H <- rbind(event_data$tasks_2E1H, 
                                       dog_before_2E1H$tasks, 
                                       dog_after_2E1H$tasks)
        event_data$tasks_Amm <- rbind(event_data$tasks_Amm,
                                      dog_before_Amm$tasks, 
                                      dog_after_Amm$tasks)
        event_data$time_budget_2E1H <- rbind(event_data$time_budget_2E1H,
                                             dog_before_2E1H$time_budget, 
                                             dog_after_2E1H$time_budget)
        event_data$time_budget_Amm <- rbind(event_data$time_budget_Amm, 
                                            dog_before_Amm$time_budget, 
                                            dog_after_Amm$time_budget)
        event_data$counts_2E1H <- rbind(event_data$counts_2E1H, 
                                        dog_before_2E1H$event_counts, 
                                        dog_after_2E1H$event_counts)
        event_data$counts_Amm <- rbind(event_data$counts_Amm, 
                                       dog_before_Amm$event_counts,
                                       dog_after_Amm$event_counts)
        event_data$time_changes_2E1H <- rbind(event_data$time_changes_2E1H, 
                                              time_changes_2E1H)
        event_data$time_changes_Amm <- rbind(event_data$time_changes_Amm, 
                                             time_changes_Amm)
        event_data$count_changes_2E1H <- rbind(event_data$count_changes_2E1H, 
                                               count_changes_2E1H)
        event_data$count_changes_Amm <- rbind(event_data$count_changes_Amm, 
                                              count_changes_Amm)
      }

    }
    return(event_data)
    
  } else if(opt == "individual"){
    individual_data <- list()
    
    for(k in levels(task_list$Dog)){
      dog_tasks <- task_list[task_list$Dog == k, ]
      print(paste("Dog:",k,", Trained:", dog_tasks$Trained))
      if(dog_tasks$Trial[1] > 3){
        dog_before_2E1H <- NULL
        dog_before_Amm <- NULL
        dog_after_2E1H <- load_ethogram_individuals(dog_tasks$Trial[1], k, 1, T)
        dog_after_Amm <- load_ethogram_individuals(dog_tasks$Trial[1], k, 2, T)
        warning("Trial 4 and 5 runs don't correspond to the same chemicals!")
      }else{
        dog_before_2E1H <- load_ethogram_individuals(1, k, 1, F)
        dog_before_Amm <- load_ethogram_individuals(1, k, 2, F)
        dog_after_2E1H <- load_ethogram_individuals(3, k, 1, T)
        dog_after_Amm <- load_ethogram_individuals(3, k, 2, T)
      }
      if(is.null(dog_before_2E1H) | is.null(dog_after_2E1H) | is.null(dog_before_Amm) | 
         is.null(dog_before_Amm)){
        time_changes_2E1H <- NULL
        time_changes_Amm <- NULL
        state_count_changes_2E1H <- NULL
        state_count_changes_Amm <- NULL
        point_count_changes_2E1H <- NULL
        point_count_changes_Amm <- NULL
      }else{
        time_changes_2E1H <- calc_time_differences(
          dog_before_2E1H$states_time_budget,
          dog_after_2E1H$states_time_budget)
        time_changes_Amm <- calc_time_differences(
          dog_before_Amm$states_time_budget,
          dog_after_Amm$states_time_budget)
        state_count_changes_2E1H <- calc_count_differences(
          dog_before_2E1H$state_counts,
          dog_after_2E1H$state_counts)
        state_count_changes_Amm <- calc_count_differences(
          dog_before_Amm$state_counts,
          dog_after_Amm$state_counts)
        point_count_changes_2E1H <- calc_count_differences(
          dog_before_2E1H$point_counts,
          dog_after_2E1H$point_counts)
        point_count_changes_Amm <- calc_count_differences(
          dog_before_Amm$point_counts,
          dog_after_Amm$point_counts)
      }
      if(k == "1"){
        individual_data$states_2E1H <- rbind(dog_before_2E1H$states, 
                                             dog_after_2E1H$states)
        individual_data$states_Amm <- rbind(dog_before_Amm$states, 
                                            dog_after_Amm$states)
        
        individual_data$time_budget_2E1H <- rbind(
          dog_before_2E1H$states_time_budget, 
          dog_after_2E1H$states_time_budget)
        individual_data$time_budget_Amm <- rbind(
          dog_before_Amm$states_time_budget, 
          dog_after_Amm$states_time_budget)
        
        individual_data$time_changes_2E1H <- time_changes_2E1H
        individual_data$time_changes_Amm <- time_changes_Amm
        
        individual_data$state_counts_2E1H <- rbind(
          dog_before_2E1H$state_counts, dog_after_2E1H$state_counts
        )
        individual_data$state_counts_Amm <- rbind(
          dog_before_Amm$state_counts, dog_after_Amm$state_counts
        )
        
        individual_data$state_counts_change_2E1H <- state_count_changes_2E1H
        individual_data$state_counts_change_Amm <- state_count_changes_Amm
        
        individual_data$points_2E1H <- rbind(dog_before_2E1H$points,
                                                   dog_after_2E1H$points)
        individual_data$points_Amm <- rbind(dog_before_Amm$points,
                                                  dog_after_Amm$points)
        
        individual_data$point_counts_2E1H <- rbind(dog_before_2E1H$point_counts,
                                             dog_after_2E1H$point_counts)
        individual_data$point_counts_Amm <- rbind(dog_before_Amm$point_counts,
                                            dog_after_Amm$point_counts)
        
        individual_data$point_counts_change_2E1H <- point_count_changes_2E1H
        individual_data$point_counts_change_Amm <- point_count_changes_Amm
        
      } else{
        individual_data$states_2E1H <- rbind(individual_data$states_2E1H, 
                                             dog_before_2E1H$states, 
                                             dog_after_2E1H$states)
        individual_data$states_Amm <- rbind(individual_data$states_Amm, 
                                            dog_before_Amm$states, 
                                            dog_after_Amm$states)
        
        individual_data$time_budget_2E1H <- rbind(
          individual_data$time_budget_2E1H,
          dog_before_2E1H$states_time_budget, 
          dog_after_2E1H$states_time_budget)
        individual_data$time_budget_Amm <- rbind(
          individual_data$time_budget_Amm,
          dog_before_Amm$states_time_budget, 
          dog_after_Amm$states_time_budget)
        
        individual_data$time_changes_2E1H <- rbind(
          individual_data$time_changes_2E1H, 
          time_changes_2E1H)
        individual_data$time_changes_Amm <- rbind(
          individual_data$time_changes_Amm,
          time_changes_Amm)
        
        individual_data$state_counts_2E1H <- rbind(
          individual_data$state_counts_2E1H,
          dog_before_2E1H$state_counts,
          dog_after_2E1H$state_counts)
        individual_data$state_counts_Amm <- rbind(
          individual_data$state_counts_Amm, 
          dog_before_Amm$state_counts,
          dog_after_Amm$state_counts)
        
        individual_data$state_counts_change_2E1H <- rbind(
          individual_data$state_counts_change_2E1H, 
          state_count_changes_2E1H
        )
        individual_data$state_counts_change_Amm <- rbind(
          individual_data$state_counts_change_Amm, 
          state_count_changes_Amm
        )
        
        individual_data$points_2E1H <- rbind(
          individual_data$points_2E1H, 
          dog_before_2E1H$points, dog_after_2E1H$points
        )
        individual_data$points_Amm <- rbind(
          individual_data$points_Amm, 
          dog_before_Amm$points, dog_after_Amm$points
        )
        
        individual_data$point_counts_2E1H <- rbind(
          individual_data$point_counts_2E1H, 
          dog_before_2E1H$point_counts, 
          dog_after_2E1H$point_counts)
        
        individual_data$point_counts_Amm <- rbind(
          individual_data$point_counts_Amm, 
          dog_before_Amm$point_counts, 
          dog_after_Amm$point_counts)
        
        individual_data$state_counts_change_2E1H <- rbind(
          individual_data$state_counts_change_2E1H, 
          state_count_changes_2E1H)
        individual_data$state_counts_change_Amm <- rbind(
          individual_data$state_counts_change_Amm, 
          state_count_changes_Amm)
        
      }
      
    }
    return(individual_data)
  }
  
}

# Loads task ethogram data 
load_ethogram_tasks <- function(trial, dog, run, trained){
  # Initialize list object to store data
  events <- list()
  # Loads data file
  if(!file.exists(paste0("./data/ethogram/task_csv_files/T", 
            trial, "D", dog, "R", run, 
            "_sideview_events_1.csv"))){
    events <- NULL
  }else{
    test_dat <- read.table(paste0("./data/ethogram/task_csv_files/T", trial, 
                                  "D", dog, "R", run, "_sideview_events_1.csv"), 
                           nrow=1, sep = ",")
    trialcode <- paste0("T", trial, "D", dog, "R", run)
    split_test_code <- strsplit(test_dat$V2, "_")
    if(trialcode != split_test_code[[1]][1]){
      warning(paste(trialcode,"file had something wrong, please check!"))
      events <- NULL
    }else{
      events_dat <- read.csv(paste0("./data/ethogram/task_csv_files/T", 
                                    trial, "D", dog, "R", run, 
                                    "_sideview_events_1.csv"), 
                             skip = 15)
      
      # Overall trial information
      events$media_location <- events_dat$Media.file.path[1]
      events$video_total_length <- events_dat$Total.length[1]
      events$fps <- events_dat$FPS[1]
      events$event_code <- paste0("T",trial,"D",dog,"R",run)
      events$target <- ifelse(run==1, "2E1H", "Ammonia")
      events$trained <- trained
      # Construct data frame for task events:
      events$tasks <- data.frame("dog" = rep(dog, nrow(events_dat)), 
                                 "trial" = rep(trial, nrow(events_dat)), 
                                 "run" = rep(run, nrow(events_dat)), 
                                 "target" = rep(events$target, nrow(events_dat)),
                                 "trained" = rep(events$trained, nrow(events_dat)),
                                 "behavior" = events_dat$Behavior, 
                                 "event_time" = events_dat$Time, 
                                 "status" = events_dat$Status
                                 )
      events$tasks$behavior <- factor(events$tasks$behavior, 
                                      levels = c("t", "c" , "o", "a"))
      events$tasks$status <- factor(str_to_lower(events$tasks$status), 
                                    levels = c("start", "stop"))
      
      events$tasks$elapsed_time <- find_elapsed_times(events$tasks)
      time_dat <- calc_time_budget(events$tasks)
      events$time_budget <- data.frame("dog" = rep(dog, nlevels(events$tasks$behavior) + 1), 
                                       "trained" = rep(trained, nlevels(events$tasks$behavior) + 1),
                                       "target" = rep(events$target, nlevels(events$tasks$behavior) + 1),
                                       "behavior" = c(NA, levels(events$tasks$behavior)), 
                                       "total_times" = time_dat[, "totals"],
                                       "fractions" = time_dat[, "fractions"])
      events$event_counts <- count_events(events$tasks)
    }
  }
  return(events)
}

# Loads ethogram data for individual behaviors (states and points)
load_ethogram_individuals <- function(trial, dog, run, trained){
  ind_behaviors <- list()
  if(!file.exists(paste0("./data/ethogram/individual_behavior_csv_files/B_T",
                         trial,"D",dog,"R",run,"_events.csv"))){
    ind_behaviors <- NULL
  }else{
    ind_dat <- read.csv(paste0("./data/ethogram/individual_behavior_csv_files/B_T",
                               trial,"D",dog,"R",run,"_events.csv"), 
                        skip = 15)
    ind_behaviors$media_location <- ind_dat$Media.file.path[1]
    ind_behaviors$video_total_length <- ind_dat$Total.length[1]
    ind_behaviors$fps <- ind_dat$FPS[1]
    ind_behaviors$event_code <- paste0("T",trial,"D",dog,"R",run)
    ind_behaviors$target <- ifelse(run==1, "2E1H", "Ammonia")
    ind_behaviors$trained <- trained
    # Construct data frame for task events:
    df_point <- ind_dat[ind_dat$Status == "POINT", ] # point behaviors
    df_states <- ind_dat[ind_dat$Status != "POINT", ] # state behaviors
    ind_behaviors$states <- data.frame("dog" = rep(dog, nrow(df_states)), 
                                      "target" = rep(ind_behaviors$target, 
                                                     nrow(df_states)),
                                      "trained" = rep(ind_behaviors$trained, 
                                                      nrow(df_states)),
                                      "behavior" = df_states$Behavior, 
                                      "event_time" = df_states$Time, 
                                      "status" = df_states$Status
    )
    # Replacing replicate task codes
    ind_behaviors$states$behavior <- factor(ind_behaviors$states$behavior, 
                                    levels = c("p", "j", "s", "l"))
    ind_behaviors$states$status <- factor(str_to_lower(ind_behaviors$states$status), 
                                  levels = c("start", "stop"))
    ind_behaviors$state_counts <- count_events(ind_behaviors$states)
    ind_behaviors$states$elapsed_time <- find_elapsed_times(ind_behaviors$states)
    time_dat <- calc_time_budget(ind_behaviors$states)
    ind_behaviors$states_time_budget <- data.frame("dog" = rep(dog, nlevels(ind_behaviors$states$behavior) + 1), 
                                     "trained" = rep(trained, nlevels(ind_behaviors$states$behavior) + 1),
                                     "target" = rep(ind_behaviors$target, nlevels(ind_behaviors$states$behavior) + 1),
                                     "behavior" = c(NA, levels(ind_behaviors$states$behavior)), 
                                     "total_times" = time_dat[, "totals"],
                                     "fractions" = time_dat[, "fractions"])
    ind_behaviors$points <- data.frame("dog" = rep(dog, nrow(df_point)), 
                                       "target" = rep(ind_behaviors$target, 
                                                      nrow(df_point)),
                                       "trained" = rep(ind_behaviors$trained, 
                                                       nrow(df_point)),
                                       "behavior" = df_point$Behavior,
                                       "event_time" = df_point$Time
    )
    ind_behaviors$points$behavior[ind_behaviors$points$behavior=="t"] <- "w"
    ind_behaviors$points$behavior[ind_behaviors$points$behavior=="c"] <- "k"
    ind_behaviors$points$behavior <- factor(ind_behaviors$points$behavior, 
                                            levels = c("h", "e", "f", "b", "w", "v",
                                                       "n", "k"))
    ind_behaviors$point_counts <- count_events(ind_behaviors$points)
  }
  return(ind_behaviors)
}

# Counts number of behavior events in a data set
count_events <- function(dat){
  counts <- matrix(0, nrow=1, ncol=nlevels(dat$behavior))
  for(j in 1:nlevels(dat$behavior)){
    colnames(counts) <- levels(dat$behavior)
    dat_sub <- dat[dat$behavior==levels(dat$behavior)[j],]
    if(nrow(dat_sub) != 0){
      starts <- dat_sub$event_time[dat_sub$status == "start"]
      if(length(starts)>0){
        counts[j] <- length(starts)
      }else{
        counts[j] <- nrow(dat_sub)
      }
      
    }
  }
  df_counts <- tidyr::pivot_longer(data.frame(counts), cols = everything())
  df_counts <- data.frame("dog" = dat$dog[1], 
                          "target" = dat$target[1],
                          "trained" = dat$trained[1],
                          "behavior" = df_counts$name, 
                          "counts" = df_counts$value)
  return(df_counts)
}

# Link individual behaviors with tasks
link_tasks_behaviors <- function(points_sub, events_sub){
  asso_task <- rep(NA, nrow(points_sub))
  for (j in 1:nrow(points_sub)){
    test_time <- points_sub$event_time[j] 
    for(k in 1:nlevels(events_sub$behavior)){
      events_extra_sub <- events_sub[events_sub$behavior == 
                                       levels(events_sub$behavior)[k], ]
      if(nrow(events_extra_sub) == 0){
        
      }else{
        starts <- events_extra_sub$event_time[
          events_extra_sub$status == "start"] <= test_time
        stops <- events_extra_sub$event_time[
          events_extra_sub$status == "stop"] >= test_time
        event_temp <- ifelse(starts & stops, T, F)
        if(all(is.na(event_temp))){
          asso_task[j] <- NA
        }else if(any(event_temp) & is.na(asso_task[j])){
          asso_task[j] <- levels(events_sub$behavior)[k]
        } else if(any(event_temp) & !is.na(asso_task[j])){
          asso_task[j] <- paste(asso_task[j], 
                                levels(events_sub$behavior)[k], 
                                sep = ",")
        }
      }
    }
  }
  return(asso_task)
}

# Compares all individual behaviors to task event times to see if they associate
compare_task_point <- function(task_list, event_data, individual_data){
  target <- c("2E1H", "Amm")
  trained <- c(F, T)
  for(pip in 1:2){
    type <- c("points_", "states_")[pip]
  for(p in 1:length(target)){
    dat_points <- individual_data[[paste0(type, target[p])]]
    if(type == "states_") dat_points <- dat_points[dat_points$status == "start",]
    dat_points$asso_task <- rep(NA, nrow(dat_points))
    dat_events <- event_data[[paste0("tasks_", target[p])]]
    for(n in 1:nlevels(task_list$Dog)){
      #print(levels(task_list$Dog)[n])
      for(u in 1:length(trained)){
        points_sub <- dat_points[ dat_points$dog == levels(task_list$Dog)[n] &
                                    dat_points$trained == trained[u], ]
        events_sub <- dat_events[ dat_events$dog == levels(task_list$Dog)[n] &
                                    dat_events$trained == trained[u], ]
        dat_points$asso_task[ dat_points$dog == 
                                levels(task_list$Dog)[n] &  
                                dat_points$trained == trained[u]] <-
          link_tasks_behaviors(points_sub, events_sub)
      }
    }
    if(type == "states_") {
      individual_data[[paste0(type, target[p])]]$asso_task[ 
        individual_data[[paste0(type, target[p])]]$status == "start"] <- 
        dat_points$asso_task
    }else{
      individual_data[[paste0(type, target[p])]]$asso_task <- 
      dat_points$asso_task
    }
  }
  }
  return(individual_data)
}

# Calculates the elapsed times for state behaviors
find_elapsed_times <- function(dat){
  elapsed_time <- rep(NA, nrow(dat))
  for(j in 1:nlevels(dat$behavior)){
    dat_sub <- dat[dat$behavior==levels(dat$behavior)[j],]
    if(nrow(dat_sub) != 0){
      starts <- dat_sub$event_time[dat_sub$status == "start"]
      stops <- dat_sub$event_time[dat_sub$status == "stop"]
      elapsed_time[dat$behavior == levels(dat$behavior)[j] & 
                                  dat$status == "stop"] <- stops - starts
    }
  }
  return(elapsed_time)
}
  
# Calculates the time budget (percent of total time) for state behaviors
calc_time_budget <- function(dat){
  total_time <- diff(range(dat$event_time))
  time_dat <- matrix(NA, nrow = (nlevels(dat$behavior)+1), ncol = 2)
  row.names(time_dat) <- c("total", levels(dat$behavior))
  #row.names(time_dat) <- c("total", "off_task", "casting", "on_odor", "alert")
  colnames(time_dat) <- c("totals", "fractions")
  time_dat[1,] <- c(total_time, 1)
  for(k in 1:nlevels(dat$behavior)){
    time_dat[k+1,] <- c(sum(dat$elapsed_time[dat$behavior == 
                                               levels(dat$behavior)[k]], 
                            na.rm = T), 
                      sum(dat$elapsed_time[dat$behavior == 
                                             levels(dat$behavior)[k]], 
                          na.rm = T)/total_time)
  }
  return(time_dat)
}

# Calculates the difference in time budgets before and after training
calc_time_differences <- function(time_budget_before, time_budget_after){
  change_time <- time_budget_after$total_times -
    time_budget_before$total_times
  change_frac <- time_budget_after$fractions -
    time_budget_before$fractions
  
  time_changes <- data.frame(
    "dog" = rep(time_budget_after$dog[1], 5),
    "target" = rep(time_budget_after$target[1], 5),
    "behavior" = time_budget_after$behavior[1:5],
    "time" = change_time,
    "frac" = change_frac
  )
  return(time_changes)
}

# Calculates the count of behavior differences before and after training
calc_count_differences <- function(counts_before, counts_after){
  behaviors <- counts_before$behavior
  counts <- matrix(NA, nrow=1, ncol = length(behaviors))
  colnames(counts) <- behaviors
  for(j in 1:length(behaviors)){
    counts[j] <- counts_after$counts[counts_after$behavior == behaviors[j]] - 
      counts_before$counts[counts_before$behavior == behaviors[j]]
  }
  
  counts_df <- tidyr::pivot_longer(data.frame(counts), cols=everything())
  counts_df <- data.frame(counts_before$dog[1], 
                          counts_before$target[1],
                          counts_df)
  colnames(counts_df) <- c("dog", "target", "behavior", "count_change")
  return(counts_df)
}

# Find and summarize behaviors linked with tasks
find_task_counts <- function(dat, behavior){
  temp_vec <- dat$asso_task[dat$behavior == behavior]
  temp_vec <- temp_vec[- which(is.na(temp_vec))] 
  temp_vec2 <- strsplit(temp_vec, split = ",")
  temp_vec2 <- unlist(temp_vec2)
  levels_vec <- c("t", "c", "o", "a")
  count <- rep(NA, length(levels_vec))
  for(i in 1:length(count)){
    count[i] <- sum(temp_vec2 == levels_vec[i])
  }
  return(count)
}

# Summarizing linked behaviors and tasks
summarize_linked_behaviors <- function(individual_data, dogs, type){
  targets <- c("2E1H", "Amm")
  traineds <- c(F, T)
  df <- data.frame()
  for(target in targets){
    for(dog in levels(dogs)){
      for(trained in traineds){
        dat <- individual_data[[paste0(type, target)]]
        dat <- dat[dat$dog == dog & 
                     dat$trained == trained,]
        behavior_levels <- levels(dat$behavior)
        for(j in 1:length(behavior_levels)){
          count_temp <- find_task_counts(dat,  
                                         behavior_levels[j])
          times_temp <- find_task_times(dat, behavior_levels[j])
          df_temp <- data.frame("dog" = rep(dog, 4),
                                "target" = rep(target, 4),
                                "trained" = rep(trained, 4),
                                "behavior" =
                                  rep(behavior_levels[j], 4), 
                                "task" = c("t","c","o","a"), 
                                "count"  = count_temp,
                                "times" = times_temp)
          df <- rbind(df, df_temp)
          rm(df_temp, count_temp)
        }
      }
    }
  }
  return(df)
}

find_task_times <- function(dat, behavior_level){
  task_levels <- c("t","c","o","a")
  dat_temp<- dat[dat$behavior == behavior_level,]
  times <- rep(NA, length(task_levels))
  for(i in 1:length(task_levels)){
    times[i] <- sum(dat_temp$event_time[dat_temp$asso_task == task_levels[i]], na.rm = T)
  }
  return(times)
}

# Other important things: 
point_behavior_labels <- c("shoulder shrug", "eye convergence", "facial tensing",
                           "head turn",  "lip lick",  "ear prick",
                           "head snap", "tail wag")
state_behavior_labels <- c("jumping", "look at handler", "pawing", 
                           "sniffing")
