# Functions for data analysis 

# Loads parameter list
load_sniff_parameters <- function(){
  dat <- read.csv("./data/data_parameters.csv", 
                  colClasses = c("factor", "factor", "factor", 
                                 "factor", "integer", "numeric",
                                 "numeric", "numeric"))
  dat$target <- ifelse(dat$chemical == "2E1H" | dat$chemical == "Ammonia", T, F)
  return(dat)
}

load_sniff_results <- function(thedate){
  dat <- read.csv(paste0("./results/csv-files/sniffdata_", thedate, ".csv"), 
                  row.names = NULL,
                  colClasses = c("factor", "factor", "factor", "factor", 
                                 "integer", "numeric", "numeric", "numeric", 
                                 "factor", "numeric", "numeric", "numeric",
                                 "numeric","integer","numeric"))
}

find_dog_heads <- function(dog){
  filename <- paste0("./data/kinematics/head_measurements/Dog", dog, 
                     "_top_measurements.csv")
  if(!file.exists(filename)) {
    head_mean <- NA
    warning(paste0("Dog ",dog," file not found!"))
  }else{
    head_meas <- read.csv(filename)
    head_mean <- mean(head_meas$Length) 
  }
  return(head_mean)
}

calculate_snout_lengths <- function(sniff_parameters){
  head_means <- rep(NA, nlevels(sniff_parameters$dog))
  for(i in 1:length(head_means)) {
    dognum <- levels(sniff_parameters$dog)[i]
    head_means[i] <- find_dog_heads(dognum)
  }
  dog_heads <- data.frame("dog" = factor(levels(sniff_parameters$dog)), 
                          "head_mean" = head_means)
  return(dog_heads)
}

# Loads data from files
load_data <- function(trial, dog, run, chemical, sample_rate, min_stop = 7.5, plot_it = TRUE ){
  require(signal)
  filename <- paste0("data/flowsensor/Trial", trial,"Run",run,"/T",trial,"D",dog,"R",run," ",
                     chemical,".edf")
  #print(filename)
  dat <- read.table(filename, header = T)
  dat <- dat[,-1]
  colnames(dat) <- c("date_time", "sample")
  dat$time <- seq(from=0, by=1/sample_rate, length.out = nrow(dat))
  dat <- dat[,c(1,3,2)]
  min_val <- which.min(dat$sample[dat$time < min_stop])
  zero_time <- dat$time[min_val]
  dat$adj_time <- dat$time - zero_time
  if(plot_it) {
    plot(dat$adj_time, dat$sample, type="l", 
         xlab="Time (seconds)", ylab = "Flow sensor data")
    points(dat$adj_time[min_val], dat$sample[min_val], col="red")
  }
  return(dat)
}

# Loads calibration data from files
load_cal_data <- function(chemical, cal_flow, rep, sample_rate, plot_it = TRUE ){
  require(signal)
  if(cal_flow == 0.8) {
    filename <- paste0("data/flowsensor/calibration/", chemical," .8 l.m.edf")
  } else{
    filename <- paste0("data/flowsensor/calibration/", chemical," 4.1 Trial ", rep ,".edf")
  }
  dat <- read.table(filename, header = T)
  dat <- dat[,-1]
  colnames(dat) <- c("date_time", "sample")
  dat$time <- seq(from=0, by=1/sample_rate, length.out = nrow(dat))
  dat <- dat[,c(1,3,2)]
  dat$adj_time <- dat$time
  if(plot_it) {
    plot(dat$time, dat$sample, type="l", 
         xlab="Time (seconds)", ylab = "Flow sensor data")
  }
  return(dat)
}

# Plots data over an optional range, will plot filtered data in red if present
plot_range <- function(dat, starttime = NULL, endtime = NULL){
  if(starttime >= endtime) stop("starttime must be before endtime!")
  if(is.null(starttime)) starttime <- min(dat$adj_time)
  if(is.null(endtime)) endtime <- max(dat$adj_time)
  dat1<-dat[dat$adj_time >= starttime & dat$adj_time <= endtime, ]
  plot(dat1$adj_time, dat1$sample, type="l", 
       xlab="Time (seconds)", ylab = "Flow sensor data", ylim = range(c(dat1$sample, dat1$sample_filtered)))
  lines(x = range(dat$adj_time), y=c(0,0), col = "gray90")
  if(!is.null(dat1$sample_filtered)){
    lines(dat1$adj_time, dat1$sample_filtered, col="red")
    legend("topleft",legend=c("raw","filtered"),col=c("black","red"), 
           lty = c(1,1))
  }
}

# Filters data with the signal package
filter_data <- function(dat, rezero = T, plotit = T, N, filsig){
  require(signal)
  # Butterworth filtering
  bf <- butter(N, filsig)
  dat$sample_filtered <- filtfilt(bf, dat$sample)
  # Setting up to rezero for sensor drift
  dat_subset<-dat[dat$adj_time >= 1,]
  dat_lim <- std(dat_subset$sample_filtered)
  i <- 3.0
  repeat{
    dat_sub_try <-dat_subset[dat_subset$adj_time >= i & 
                               dat_subset$adj_time <= i + 3,]
    check_sum <- sum(dat_sub_try$sample_filtered > 1.5*dat_lim & 
                       dat_sub_try$sample_filtered < -1.5*dat_lim)
    lm.try <- lm(sample_filtered~adj_time, data = dat_sub_try)
    p.value <- summary(lm.try)$coefficients[2,4]
    if(check_sum == 0 & p.value > 0.05) {
      break
    } else {
      i <- i + 1
    }
    if(i+3 > max(dat$adj_time)) stop("There was no subset of data appropriate to 
                                     calculate a rezero. :( ")
  }
  # Calculating mean to rezero
  mean_filt <- mean(dat_sub_try$sample_filtered, na.rm = T)
  if(rezero) {
    print(paste("Rezeroing by:", mean_filt))
    dat$sample_filtered <- dat$sample_filtered - mean_filt
  }
  # Plotting if indicated
  if(plotit) plot_range(dat, min(dat$adj_time), max(dat$adj_time))
  return(dat)
}

# Calculates sniff frequency, sniff volume, and sniff volume flow rate
sniff_it <- function(dat, starttime, endtime, checkpeaks = T){
  dat_subset <- dat[dat$adj_time <= endtime & 
                      dat$adj_time >= starttime,]
  peaks_high <- findpeaks(dat_subset$sample_filtered, minpeakdistance = 10, 
                          minpeakheight = 0)
  peaks_low <- findpeaks(-dat_subset$sample_filtered, minpeakdistance = 10, 
                         minpeakheight = -10)
  if(checkpeaks){
    plot(dat_subset$adj_time, dat_subset$sample_filtered, type = "l")
    points(y = peaks_high[, 1], x = dat_subset$adj_time[peaks_high[, 2]], 
           col = "red")
    points(y = dat_subset$sample_filtered[peaks_low[, 2]], 
           x = dat_subset$adj_time[peaks_low[, 2]], 
           col="blue")
  }
  peaks_high <- cbind(peaks_high, dat_subset$adj_time[peaks_high[, 2]])
  peaks_low <- cbind(peaks_low, dat_subset$adj_time[peaks_low[, 2]])
  if(nrow(peaks_high) < 2) {
    sniff_freq <- NA
    sniff_duration <- NA
  } else {
    sniff_freq <- 1/diff(sort(peaks_high[, 5]))
    sniff_duration <- diff(sort(peaks_low[,5]))
  }
  lows_list <- sort(peaks_low[, 2])
  if(length(lows_list) < 2){
    sniff_volume <- NA
    sniff_vfr <- NA
  }else{
    sniff_volume <- rep(NA, length(lows_list)-1)
    sniff_vfr <- rep(NA, length(lows_list)-1)
    for(i in 2:length(lows_list)){
      low_seq <- seq(lows_list[i-1],lows_list[i])
      sniff_volume[i-1] <- sum(dat_subset$sample_filtered[low_seq] - 
                                 min(dat_subset$sample_filtered[low_seq[1]]))
      sniff_vfr[i-1] <- sniff_volume[i-1]/(dat_subset$adj_time[max(low_seq)] - 
                                             dat_subset$adj_time[min(low_seq)])
    }
  }
  return(list("sniff_freq" = sniff_freq, 
              "sniff_duration" = sniff_duration, 
              "sniff_volume" = sniff_volume, 
              "sniff_vfr" = sniff_vfr))
}

convert_final_df <- function(sniff_parameters, sniff_dat){
  vec <- nrow(sniff_parameters)
  final_df <- data.frame(
    "sniff_freq" = rep(NA, vec),
    "sniff_duration" = rep(NA, vec),
    "sniff_volume"= rep(NA, vec),
    "sniff_vfr"= rep(NA, vec),
    "sniff_counts"= rep(NA, vec),
    "total_volume"= rep(NA, vec)
  )
  for(k in 1:vec){
      final_df[k,1] <- mean(unlist(sniff_dat[[k]][[1]]))
      final_df[k,2] <- mean(unlist(sniff_dat[[k]][[2]]))
      final_df[k,3] <- (mean(unlist(sniff_dat[[k]][[3]]))/60)*1000 # Converting from standard liters to mL
      final_df[k,4] <- (mean(unlist(sniff_dat[[k]][[4]]))/60)*1000
      final_df[k,5] <- length(unlist(sniff_dat[[k]][[3]]))
      final_df[k,6] <- (sum(unlist(sniff_dat[[k]][[3]]))/60)*1000
  }
  final_df2 <- cbind(sniff_parameters,final_df)
  return(final_df2)
}

correct_flowsensor_times <- function(dat, frame_rate){
  dat$adj_start_time <- rep(NA, nrow(dat))
  dat$adj_end_time <- rep(NA, nrow(dat))
  for(j in 1:nrow(dat)){
    dog <- dat$dog[j]
    run <- dat$run[j]
    trial <- dat$trial[j]
    chemical <- ifelse(as.character(dat$chemical[j]) == "2E1H", "X2E1H", 
                       as.character(dat$chemical[j]))
    offsets <- read.table(paste0("./data/flowsensor/FrameOffsets/T",
                                 trial, "R", run, "_FrameOffsets_flow.txt"),
                          sep = ",", header = T, skip = 1)
    trialcode <- find_trial_code(dat, j)
    which_chem <- read.csv("./data/flowsensor/chem_locations.csv", skip = 1)
    which_hole <- which_chem[which_chem$trial == trial & which_chem$run == run, chemical]
    time_offset <- offsets[ offsets$trial_code == trialcode, which_hole] / frame_rate
    dat$adj_start_time[j] <- dat$start_time[j] + time_offset
    dat$adj_end_time[j] <- dat$end_time[j] + time_offset
  }
  return(dat)
}

link_tasks_sniffs <- function(sniff_results, event_data){
  sniff_results$asso_tasks <- rep(NA, nrow(sniff_results))
  for(bip in 1:nrow(sniff_results)){
    test_sniff <- sniff_results[bip,]
    trial <- as.numeric(test_sniff$trial)
    dog <- as.numeric(test_sniff$dog)
    run <- as.numeric(test_sniff$run)
    chemical <- ifelse(run == 1, "2E1H", "Amm")
    event_name <- paste0("tasks_",chemical)
    tasks <- event_data[[event_name]]
    tasks_sub <- tasks[tasks$trial == trial & tasks$dog == dog, ]
    if(nrow(tasks_sub) == 0 ) next
    num_tasks <- sum(tasks_sub$status == "start")
    start_list <- which(tasks_sub$status == "start")
    end_list <- which(tasks_sub$status == "stop")
    for(ts in 1:num_tasks){
      start_time <- tasks_sub$event_time[start_list[ts]]
      end_time <- tasks_sub$event_time[end_list[ts]]
      spot_task <- as.character(tasks_sub$behavior[start_list[ts]])
      if(test_sniff$adj_start_time > start_time & test_sniff$adj_start_time < end_time  | 
         test_sniff$adj_end_time > start_time & test_sniff$adj_end_time < end_time  | 
         start_time > test_sniff$adj_start_time & end_time < test_sniff$adj_end_time){
        if(is.na(test_sniff$asso_tasks)){
          test_sniff$asso_tasks <- spot_task
        }else if (grepl(spot_task, test_sniff$asso_tasks, fixed = TRUE)){
        }else{
          test_sniff$asso_tasks <- paste(test_sniff$asso_tasks,spot_task, sep = ",")
        }
      }
    }
    sniff_results$asso_tasks[bip]<- test_sniff$asso_tasks
  }
  return(sniff_results)
}

find_sniffs_per_task <- function(sniff_results, task){
  list_sniffs <- grepl(task, sniff_results$asso_tasks)
  sniffs <- sniff_results[list_sniffs,]
  sniffs$asso_tasks <- task
  return(sniffs)
}



