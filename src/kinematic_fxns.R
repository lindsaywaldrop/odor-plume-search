# Kinematic analysis functions


#### Functions ####

load_dataset <- function(trial, dog, run, frame_rate, holes_condition = T){
  dat <- read.csv(paste0("./data/kinematics/Trial", trial, "Run", run,
                         "/T",trial,"D",dog,"R",run,"_DLTdv8_data_xyzpts.csv"))
  nose_dat <- data.frame("frame" = seq(1,nrow(dat)), 
                         "time" = seq(from = 0,by = (1/frame_rate), 
                                      length.out = nrow(dat)),
                         "x" = dat$pt1_X, 
                         "y" = dat$pt1_Y, 
                         "z" = dat$pt1_Z, 
                         "point" = rep("nose",nrow(dat)))
  stop_dat <- data.frame("frame" = seq(1,nrow(dat)), 
                         "time" = seq(from = 0,by = (1/frame_rate), 
                                      length.out = nrow(dat)),
                         "x" = dat$pt2_X, 
                         "y" = dat$pt2_Y, 
                         "z" = dat$pt2_Z,
                         "point" = rep("stop", nrow(dat)))
  if(holes_condition){
    hole_pts <- data.frame("x" = dat$pt3_X, 
                           "y" = dat$pt3_Y, 
                           "z" = dat$pt3_Z)
    hole_pts <- na.omit(hole_pts)
    if(nrow(hole_pts)!=10){
      stop("The calibration points were not marked correctly for this event.")
    }
    row.names(hole_pts) <- c("top_outer_4", "bottom_outer_4",
                                 "right_outer_4", "left_outer_4", 
                                 "floor_mat", "bottom_inner_4", 
                                 "bottom_inner_3", "bottom_inner_2",
                                 "bottom_inner_1", "high_point")
  }
  
  tracks <- rbind(nose_dat,stop_dat)
  tracks[tracks == "NaN"] <- NA
  #
  if(holes_condition){
    return(list("tracks" = tracks, "holes" = hole_pts))
  }else{
    return(list("tracks" = tracks))
  }
  
}

load_event_list <- function(opt = "default"){
  if(opt == "default"){
    dat <- read.csv("./data/kinematics/kinematic_events.csv", skip = 0, header = T, 
                    colClasses = c("factor", "factor", "factor"))
  }else if(opt == "kinematics-paper"){
    dat <- read.csv("./data/kinematics/kinematic_events_kinematics-paper.csv", 
                    skip = 0, header = T, 
                    colClasses = c("factor", "factor", "factor"))
  }else {
    stop("Unknown opt, please choose valid opt for event list.")
  }
  return(dat)
}

assemble_kin_data <- function(frame_rate, opt = "default"){
  # Load event list to process
  if(opt == "default"){
    event_list <- load_event_list()
  }else if(opt != "default"){
    event_list <- load_event_list(opt = opt)
  } else{
    stop("Opt is unknown, please choose valid opt.")
  }
  all_dat <- list()
  # load tracks data and correct it in 3D and align with fixed points
  for(i in 1:nrow(event_list)){
        trial <- event_list$trial[i]
        run <- event_list$run[i]
        dog <- event_list$dog[i]
        trained <- ifelse(trial=="1" | trial=="2", F, T)
        trialcode <- paste0("T",trial,"D",dog,"R",run)
        print(trialcode)
        set_dat <- load_dataset(trial, dog, run, frame_rate, T)
        values <- create_affine_correction(trial, dog, run, F)
        set_dat$tracks[,3:5] <- correct_3D_pts(set_dat$tracks[,3:5], values, trial, run, F)
        set_dat$holes <- correct_3D_pts(set_dat$holes, values, trial, run, F)
        set_dat$tracks$trial <- rep(trial, nrow(set_dat$tracks))
        set_dat$tracks$dog <- rep(dog, nrow(set_dat$tracks))
        set_dat$tracks$run <- rep(run, nrow(set_dat$tracks))
        set_dat$tracks$trained <- rep(trained, nrow(set_dat$tracks))
        all_dat[[trialcode]] <- set_dat
  }
  return(all_dat)
}

adjust_holes <- function(tholes, fholes){
  # Realigning along Y-Z plane
  tube_4 <- fholes[1,]
  tube_3 <- fholes[2,]
  tube_2 <- fholes[3,]
  tube_1 <- fholes[4,]
  
  straight_vert_line <- data.frame("Y" = c(high_pt$y, tube_4$y,
                                           tube_3$y, tube_2$y, 
                                           tube_1$y, floor$y),
                                   "Z" = c(high_pt$z, tube_4$z,
                                           tube_3$z, tube_2$z, 
                                           tube_1$z, floor$z))
  vert_line_mod <- lm (Z~Y, data = straight_vert_line)
  summary(vert_line_mod)
  if(plotit){
    plot(straight_vert_line$Y,straight_vert_line$Z)
    abline(vert_line_mod)
  }
  gamma_yz <- 1.57079633 - atan(vert_line_mod$coefficients[["Y"]]) 
  center_pt_yz <- c(floor$y, floor$z)
  
  
    # Correcting the points:
    blep <- correct_2D_pts(straight_vert_line, center_pt_yz, gamma_yz)
    
    straight_hor_line <- data.frame("X" = c(tube_1$x, tube_2$x, 
                                            tube_3$x, tube_4$x),
                                    "Y" = c(tube_1$y, tube_2$y, 
                                            tube_3$y, tube_4$y))
    
    hor_line_mod <- lm (Y~X, data = straight_hor_line)
    summary(hor_line_mod)
}
  
pivot_kinematic_data <- function(all_dat){
  n <- length(all_dat)
  tracks <- data.frame()
  for(i in 1:n){
    tracks <- rbind(tracks, na.omit(all_dat[[i]]$tracks))
  }
  return(tracks)
}

remove_prompted <- function(all_dat){
  n <- length(all_dat)
  for(i in 1:n){
    trialcode <- names(all_dat)[i]
    if(file.exists(paste0("./data/ethogram/unprompted_search_files/",
                          trialcode,"_sideview_unprompted_1.csv"))){
      timeline <- read.csv(paste0("./data/ethogram/unprompted_search_files/",
                                  trialcode,"_sideview_unprompted_1.csv"),
                           skip = 15, header = T)
      trial <- all_dat[[i]]$tracks$trial[1]
      run <- all_dat[[i]]$tracks$run[1]
      offsets <- read.table(paste0("./data/kinematics/FrameOffsets/T",
                                   trial, "R", run, "_FrameOffsets.txt"),
                            sep = ",", header = T, skip = 1)
      timeline$Time <- timeline$Time - ((offsets$movie1_side_offset[offsets$trial_code == trialcode]-1)/
                                        timeline$FPS[1])
      all_dat[[i]]$tracks$x <- ifelse(all_dat[[i]]$tracks$time < timeline$Time[1] | 
                                        all_dat[[i]]$tracks$time > timeline$Time[2], 
                                      NA, all_dat[[i]]$tracks$x)
      all_dat[[i]]$tracks$y <- ifelse(all_dat[[i]]$tracks$time < timeline$Time[1] | 
                                        all_dat[[i]]$tracks$time > timeline$Time[2], 
                                      NA, all_dat[[i]]$tracks$y)
      all_dat[[i]]$tracks$z <- ifelse(all_dat[[i]]$tracks$time < timeline$Time[1] | 
                                        all_dat[[i]]$tracks$time > timeline$Time[2], 
                                      NA, all_dat[[i]]$tracks$z)
    }else {
      warning(paste(trialcode,"does not have unprompted timing file!"))
    }
  }
  return(all_dat)
}

clip_to_event <- function(dat, start_time, end_time){
  dat <- dat[dat$time >= start_time & dat$time <= end_time, ]
  return(dat)
}

calculate_angles <- function(dat){
  n <- unique(dat$frame)
  #dat$yaw <- rep(NA, nrow(dat))
  #dat$pitch <- rep(NA, nrow(dat))
  for(j in n){
    #subsets data
    sub_dat <- dat[dat$frame == j,]
    if (nrow(sub_dat) > 2) stop("There are replicate points - check the run and dog again.")
    # tests data for NA's, skips if found 
    if(any(is.na(sub_dat$x)) || any(is.na(sub_dat$y)) || 
       any(is.na(sub_dat$z)) || nrow(sub_dat) < 2){
    }else {
      # rezeros data on nose point
      sub_dat$x_new <- sub_dat$x - sub_dat$x[sub_dat$point == "nose"]
      sub_dat$y_new <- sub_dat$y - sub_dat$y[sub_dat$point == "nose"]
      sub_dat$z_new <- sub_dat$z - sub_dat$z[sub_dat$point == "nose"]
      # calculates yaw in degrees
      yaw <- atan(sub_dat$y_new[sub_dat$point == "stop"] / 
                    sub_dat$x_new[sub_dat$point == "stop"]) * 
        (360 / (2 * pi))
      # calculates pitch in degrees
      pitch <- atan(sub_dat$z_new[sub_dat$point == "stop"] / 
                      sub_dat$y_new[sub_dat$point == "stop"]) * 
        (360 / (2 * pi))
      # puts yaw on a full 360 circle
      if(sub_dat$x_new[sub_dat$point == "stop"] < 0 &&
         sub_dat$y_new[sub_dat$point == "stop"] > 0){
        yaw <- yaw + 180
      } else if(sub_dat$x_new[sub_dat$point == "stop"] < 0 &&
                sub_dat$y_new[sub_dat$point == "stop"] < 0 ){
        yaw <- yaw + 270
        #} else if(sub_dat$x_new[sub_dat$point == "stop"] > 0 &&
        #          sub_dat$y_new[sub_dat$point == "stop"] < 0 ){
        #  yaw <- yaw + 270
      } else {
        yaw <- yaw
      }
      # puts pitch on a full 360 circle
      if(sub_dat$y_new[sub_dat$point == "stop"] < 0 &&
         sub_dat$z_new[sub_dat$point == "stop"] > 0){
        pitch <- pitch + 90
      } else if(sub_dat$y_new[sub_dat$point == "stop"] < 0 &&
                sub_dat$z_new[sub_dat$point == "stop"] < 0 ){
        pitch <- pitch - 270
        #} else if(sub_dat$x_new[sub_dat$point == "stop"] > 0 &&
        #         sub_dat$y_new[sub_dat$point == "stop"] < 0 ){
        # yaw <- yaw + 270
      } else {
        pitch <- pitch-90
      }
      dat$yaw[dat$frame==j] <- rep(yaw,2)
      dat$pitch[dat$frame==j] <- rep(pitch,2)
    }
    
  }
  return(dat)
}


count_close_frames <- function(tracks, sphere_center, radius, frame_rate, plotit=F){
  #if(any(is.na(tracks))) tracks <- na.omit(tracks)
  test_pts <- data.frame("x" = tracks$x, "y" = tracks$y, "z" = tracks$z)
  test_pts <- na.omit(test_pts)
  test_dist <- sqrt((test_pts[, 1] - sphere_center[,1])^2 + 
                      (test_pts[, 2] - sphere_center[,2])^2 + 
                      (test_pts[, 3] - sphere_center[,3])^2)
  test_results <- ifelse(na.omit(test_dist) <= radius, T, F)
  time_near_hole <- sum(test_results)/frame_rate # seconds
  if(plotit){
    clean_nose_tracks <- test_pts
    plot(clean_nose_tracks$x, 
         clean_nose_tracks$y,
         xlim = c(sphere_center$x - radius * 1.5, 
                  sphere_center$x + radius * 1.5),
         ylim = c(sphere_center$y - radius * 1.5, 
                  sphere_center$y + radius * 1.5))
    points(sphere_center$x, sphere_center$y, 
           col = "orange", pch = 19, cex = 3)
    points(clean_nose_tracks$x[test_results], 
           clean_nose_tracks$y[test_results], col = "red")
    plot(clean_nose_tracks$y, 
         clean_nose_tracks$z,
         xlim = c(sphere_center$y - radius * 1.5, 
                  sphere_center$y + radius * 1.5),
         ylim = c(sphere_center$z - radius * 1.5, 
                  sphere_center$z + radius * 1.5))
    points(sphere_center$y, sphere_center$z, 
           col = "orange", pch = 19, cex = 3)
    points(clean_nose_tracks$y[test_results], 
           clean_nose_tracks$z[test_results], col = "red")
  }
  return(time_near_hole)
}


magnitude_3D <- function(x1, x2, y1, y2, z1, z2) {
  return(sqrt((x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2))
}

calc_hot_dist <- function(dat, hole){
  dat$hot_dist <- magnitude_3D(dat$x, hole$x, dat$y, hole$y, dat$z, hole$z)
  return(dat)
}

correct_2D_pts <- function(pts, center_pt, alpha){
  rotm <- matrix(c(cos(alpha), sin(alpha), -sin(alpha), cos(alpha)), ncol = 2)
  if(!is.matrix(pts)) pts <- as.matrix(pts)
  rot_pts <- t(rotm %*% ( t(pts) - c(center_pt[1], center_pt[2])) + 
                 c(center_pt[1], center_pt[2]))
  corr_pts <- rot_pts
  corr_pts[,1] <- rot_pts[,1] - (center_pt[1])
  corr_pts[,2] <- rot_pts[,2] - (center_pt[2])
  return(corr_pts)
}

create_affine_correction <- function(trial, dog, run, plotit = FALSE){
  # Load points for affine translation
  dat <- load_dataset(trial, dog, run, 0, T)
  holes <- dat$holes
  
  # Realigning along Y-Z plane
  tube_4 <- holes[6,]
  tube_3 <- holes[7,]
  tube_2 <- holes[8,]
  tube_1 <- holes[9,]
  floor <- holes[5,]
  high_pt <- holes[10,]
  tube4_noon <- holes[1,]
  tube4_6 <- holes[2,]
  tube4_3 <- holes[3,]
  tube4_9 <- holes[4,]
  
  straight_vert_line <- data.frame("Y" = c(high_pt$y, tube_4$y,
                                           tube_3$y, tube_2$y, 
                                           tube_1$y, floor$y),
                                   "Z" = c(high_pt$z, tube_4$z,
                                           tube_3$z, tube_2$z, 
                                           tube_1$z, floor$z))
  vert_line_mod <- lm (Z~Y, data = straight_vert_line)
  summary(vert_line_mod)
  if(plotit){
    plot(straight_vert_line$Y,straight_vert_line$Z)
    abline(vert_line_mod)
  }
  gamma_yz <- 1.57079633 - atan(vert_line_mod$coefficients[["Y"]]) 
  center_pt_yz <- c(floor$y, floor$z)
  
  if(plotit){
    # Correcting the points:
    blep <- correct_2D_pts(straight_vert_line, center_pt_yz, gamma_yz)
    
    plot(tube_4$y, tube_4$z, xlim = c(-0.1, 0.1), ylim = c(-0.6, 0.7), 
         asp=1)
    points(blep[,1], blep[,2],col = "blue")
    points(tube_4$y,tube_4$z, col = "black")
    points(floor$y,floor$z, col = "red")
    points(high_pt$y, high_pt$z, col = "orange")
  }
  # Realigning along X-Y plane
  
  straight_hor_line <- data.frame("X" = c(tube_1$x, tube_2$x, 
                                          tube_3$x, tube_4$x),
                                  "Y" = c(tube_1$y, tube_2$y, 
                                          tube_3$y, tube_4$y))
  
  hor_line_mod <- lm (Y~X, data = straight_hor_line)
  summary(hor_line_mod)
  
  if(plotit){
    plot(straight_hor_line$X,straight_hor_line$Y)
    abline(hor_line_mod)
  }
  
  alpha_xy <- -atan(hor_line_mod$coefficients[["X"]]) 
  center_pt_xy <- c((tube_3$x - tube_2$x)/2, (tube_3$y - tube_2$y)/2)
  
  if(plotit){
    blepx <- correct_2D_pts(straight_hor_line, center_pt_xy, alpha_xy)
    
    plot(tube_4$x,tube_4$y, xlim = c(-0.1, 0.1), ylim = c(-0.5, 0.5), asp=1)
    points(tube_3$x,tube_3$y, col = "black")
    points(tube_1$x,tube_1$y, col = "red")
    points(tube_2$x, tube_2$y, col = "orange")
    points(blepx, col="blue")
  }
  
  ###########
  # Realigning along X-Z plane

  straight_hor_line <- data.frame("X" = c(tube_1$x, tube_2$x,
                                          tube_3$x, tube_4$x),
                                  "Z" = c(tube_1$z, tube_2$z,
                                          tube_3$z, tube_4$z))

  hor_line_mod <- lm (Z~X, data = straight_hor_line)
  summary(hor_line_mod)

  if(plotit){
    plot(straight_hor_line$X, straight_hor_line$Z)
    abline(hor_line_mod)
  }


  omega_xz <- -atan(hor_line_mod$coefficients[["X"]])
  center_pt_xz <- c((tube_3$x - tube_2$x)/2, (tube_3$y - tube_2$y)/2)

  if(plotit){
    blepz <- correct_2D_pts(straight_hor_line, center_pt_xz, omega_xz)

    plot(tube_4$x,tube_4$z, xlim = c(-0.1, 0.1), ylim = c(-0.5, 0.5), asp=1)
    points(tube_3$x,tube_3$z, col = "black")
    points(tube_1$x,tube_1$z, col = "red")
    points(tube_2$x, tube_2$z, col = "orange")
    points(tube4_noon$x, tube4_noon$z, col = "darkgreen")
    points(tube4_6$x, tube4_6$z, col = "darkgreen")
    points(blepz, col="blue")
  }
  
  ###########
  
  # Put together final table of values
  values <- data.frame("center_pt_xy_x" = center_pt_xy[1], 
                       "center_pt_xy_y" = center_pt_xy[2], 
                       "center_pt_yz_y" = center_pt_yz[1], 
                       "center_pt_yz_z" = center_pt_yz[2], 
                       "center_pt_xz_x" = center_pt_xz[1], 
                       "center_pt_xz_z" = center_pt_xz[2], 
                       "alpha_xy" = alpha_xy,
                       "gamma_yz" = gamma_yz, 
                       "omega_xz" = omega_xz)
  return(values)
}

correct_3D_pts <- function(size_dat, values, trial, run, plotit = FALSE){
  require(scatterplot3d)
  #values <- read.csv(paste0("./data/kinematics/affine_measurements/correction_values_T", 
  #                          trial, "R", run, ".csv"))
  size_dat_corrected <- size_dat
  center_pt_xz <- c(values$center_pt_xz_x, values$center_pt_xz_z)
  center_pt_yz <- c(values$center_pt_yz_y, values$center_pt_yz_z)
  center_pt_xy <- c(values$center_pt_xy_x, values$center_pt_xy_y)
  
  # adjust x, z: 
  size_dat_corrected[,c(1,3)] <- correct_2D_pts(size_dat_corrected[,c(1,3)], center_pt_xz, 
                                             values$omega_xz)
  
  # adjust y, z: 
  size_dat_corrected[,2:3] <- correct_2D_pts(size_dat_corrected[,2:3], center_pt_yz, 
                                             values$gamma_yz)
  # adjust x, y:
  size_dat_corrected[,1:2] <- correct_2D_pts(size_dat_corrected[,1:2], center_pt_xy, 
                                             values$alpha_xy)
  if(plotit){
    sp3 <- scatterplot3d(size_dat, angle=10, asp = 1)
    sp3$points3d(size_dat_corrected, col= "red", pch=19)
  }
  return(size_dat_corrected)
}

calculate_kin_error <- function(filepath_name, frame_rate){
  dat <- read.csv(filepath_name)
  nose_dat <- data.frame("frame" = seq(1,nrow(dat)), 
                         "time" = seq(from=0,by=(1/frame_rate), 
                                      length.out = nrow(dat)),
                         "x" = dat$pt1_X, "y" = dat$pt1_Y, "z" = dat$pt1_Z, 
                         "point" = rep("nose",nrow(dat)))
  nose_dat <- na.omit(nose_dat)
  
  stop_dat <- data.frame("frame" = seq(1,nrow(dat)), 
                         "time" = seq(from=0,by=(1/frame_rate), 
                                      length.out = nrow(dat)),
                         "x" = dat$pt2_X, "y" = dat$pt2_Y, "z" = dat$pt2_Z,
                         "point" = rep("stop",nrow(dat)))
  stop_dat <- na.omit(stop_dat)
  
  temp_dat <- data.frame("frame" = nose_dat$frame, 
                         "length" = rep(NA, nrow(nose_dat)))
  
  for (j in 1:nrow(nose_dat)){
    frame_number <- nose_dat$frame[j]
    if(any(stop_dat$frame==frame_number)) {
      temp_dat$length[j] <- magnitude_3D(stop_dat$x[stop_dat$frame==frame_number],
                                      nose_dat$x[nose_dat$frame==frame_number],
                                      stop_dat$y[stop_dat$frame==frame_number], 
                                      nose_dat$y[nose_dat$frame==frame_number],
                                      stop_dat$z[stop_dat$frame==frame_number], 
                                      nose_dat$z[nose_dat$frame==frame_number])
    }
  }
  
  face_lengths <- na.omit(temp_dat)
  print(shapiro.test(face_lengths$length))
  return(face_lengths)
}

correct_ethogram_times <- function(dat, frame_rate){
    dat$adj_time <- rep(NA, nrow(dat))
  for(j in 1:nrow(dat)){
    dog <- dat$dog[j]
    run <- ifelse(dat$target[j]=="2E1H", 1, 2)
    trial <- ifelse(as.numeric(as.character(dog)) < 22 & dat$trained[j] == F, 1,
                    ifelse(as.numeric(as.character(dog)) >= 22 & dat$trained[j] == F, 2,
                           3))
    offsets <- read.table(paste0("./data/kinematics/FrameOffsets/T",
                                 trial, "R", run, "_FrameOffsets.txt"),
                          sep = ",", header = T, skip = 1)
    trialcode <- find_trial_code(dat, j)
    time_offset <- offsets$movie1_side_offset[ offsets$trial_code == trialcode] / 
      frame_rate
    dat$adj_time[j] <- dat$event_time[j] - time_offset
  }
  return(dat)
}

find_behavior_times <- function(trial_info, all_dat){
  trial_info$x <- rep(NA, nrow(trial_info))
  trial_info$y <- rep(NA, nrow(trial_info))
  trial_info$z <- rep(NA, nrow(trial_info))
  
  for(i in 1:nrow(trial_info)){
    trialcode <- find_trial_code(trial_info, i)
    tracks <- all_dat[[trialcode]]$tracks
    if(!is.null(tracks)){
      frames <- tracks$frame[which(abs(tracks$time - 
                                         trial_info$event_time[i]) == 
                                     min(abs(tracks$time - 
                                               trial_info$event_time[i])))]
      trial_info$x[i] <- tracks$x[tracks$frame == frames[1] & tracks$point == "nose"]
      trial_info$y[i] <- tracks$y[tracks$frame == frames[1] & tracks$point == "nose"]
      trial_info$z[i] <- tracks$z[tracks$frame == frames[1] & tracks$point == "nose"]
    } else{
      trial_info$x[i] <- NA
      trial_info$y[i] <- NA
      trial_info$z[i] <- NA
    }
  }
  return(trial_info)
}

find_trial_code <- function(dat, j){
  dog <- dat$dog[j]
  if(is.null(dat$trial)){
    run <- ifelse(dat$target[j]=="2E1H", 1, 2)
    trial <- ifelse(as.numeric(as.character(dog)) < 22 & dat$trained[j] == F, 1,
                    ifelse(as.numeric(as.character(dog)) >= 22 & dat$trained[j] == F, 2,
                           3))
    
  } else {
    run <- dat$run[j]
    trial <- dat$trial[j]
  }
  trialcode <- paste0("T",trial, "D", dog, "R", run)
  return(trialcode)
}

assign_task <- function(track_sub, tasks_sub, task){
  starts <- tasks_sub$adj_time[tasks_sub$behavior == task & 
                                 tasks_sub$status == "start"]
  stops <- tasks_sub$adj_time[tasks_sub$behavior == task & 
                                tasks_sub$status == "stop"]
  if(length(starts) == 0){
  }else{
    for(i in 1:nrow(track_sub)){
      for(k in 1:length(starts)){
        if(track_sub$time[i] >= starts[k] & track_sub$time[i] <= stops[k]){
          track_sub$task[i] <- task
        }
      }
    }
  }
  
  return(track_sub)
}

assign_tasks_tracks <- function(tracks, tasks_info){
  tracks$task <- rep(NA, nrow(tracks))
  for (j in 1:nrow(tracks)){
    track_sub <- tracks[j,]
    dog <- tracks$dog[j]
    run <- tracks$run[j]
    trial <- tracks$trial[j]
    trained <- ifelse(trial == 3, T, F)
    tasks_sub <- tasks_info[tasks_info$dog == dog & 
                              tasks_info$run == run &
                              tasks_info$trained == trained, ]
    track_sub$task <- NA
    track_sub <- assign_task(track_sub, tasks_sub, "t")
    track_sub <- assign_task(track_sub, tasks_sub, "c")
    track_sub <- assign_task(track_sub, tasks_sub, "o")
    track_sub <- assign_task(track_sub, tasks_sub, "a")
    tracks$task[j] <- track_sub$task
    
  }
  tracks$task <- factor(tracks$task, levels=c("t", "c", "o", "a"))
  tracks$point <- factor(tracks$point)
  return(tracks)
}

num_dec_place <- function(x) {
  if ((x %% 1) != 0) {
    nchar(strsplit(as.character(x), ".", fixed=TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}

calc_kmeans <- function(scaled_tracks, run, task, n){
  require(cluster)
  tracks_cluster<- scaled_tracks[scaled_tracks$run == run  & scaled_tracks$task == task, 
                                 c("x","y","z")]
  kmeans_cluster <- kmeans(as.matrix(tracks_cluster), centers = n)
  gaps_cluster <- clusGap(tracks_cluster, FUNcluster = kmeans, K.max = n)
  cluster_data <- data.frame("center_x" = kmeans_cluster$centers[,1], 
                             "center_y" = kmeans_cluster$centers[,2],
                             "center_z" = kmeans_cluster$centers[,3],
                             "gap" = gaps_cluster$Tab[,3])
  return(list("cluster_data" = cluster_data, "cluster_points" = kmeans_cluster$cluster))
}

t_test_sniff <- function(sniff_summaries, task, n){
  
  sniff_2e1h <- sniff_summaries$value[sniff_summaries$task == task & 
                                        sniff_summaries$chemical == "2E1H"]
  sniff_Ammonia <- sniff_summaries$value[sniff_summaries$task == task & 
                                           sniff_summaries$chemical == "Ammonia"]
  p_2e <- shapiro.test(na.omit(sniff_2e1h))
  p_Am <- shapiro.test(na.omit(sniff_Ammonia))
  
  if(p_2e$p.value < 0.01 | p_Am$p.value < 0.01) {
    #warning("These data might not be normal!")
    output <- wilcox.test(sniff_2e1h, sniff_Ammonia)
    output$p.value <- p.adjust(output$p.value, method = "bonferroni", n = n)
  }else {
    output <- t.test(sniff_2e1h, sniff_Ammonia, paired = F)
    output$p.value <- p.adjust(output$p.value, method = "bonferroni", n = n)
  }
  return(output)
}

plot_angle_tracks <- function(tracks, dog, run, col_name, min_height = NULL){
  require(signal)
  require(pracma)
  cleaned_tracks <- na.omit(tracks[tracks$point=="nose" & tracks$run == run & 
                                     tracks$dog == dog,])
  if(nrow(cleaned_tracks) < 2) {
    warning(paste("Looks like dog", dog, "has an issue."))
    return(NULL)
    break
  }
  bf <- butter(4, 0.15)
  cleaned_tracks[[col_name]] <- filtfilt(bf, cleaned_tracks[[col_name]])
  
  cleaned_tracks$partition <- rep(NA, nrow(cleaned_tracks))
  blep <- 1
  for(i in 1:(nrow(cleaned_tracks)-1)){
    if(cleaned_tracks$frame[i+1]-cleaned_tracks$frame[i] == 1 &
       cleaned_tracks$task[i+1] == cleaned_tracks$task[i]){
      cleaned_tracks$partition[i] <- blep
    }else{
      blep <- blep + 1
    }
  }
  cleaned_tracks <- na.omit(cleaned_tracks)
  cleaned_tracks$partition <- factor(cleaned_tracks$partition)
  rezero <- function(dat){
    m <- mean(dat, na.rm = T)
    dat2 <- dat - m
    return(dat2)
  }
  if(is.null(min_height)) {
    peaks_high <- findpeaks(rezero(cleaned_tracks[[col_name]]), 
                            nups = 4, ndowns = 4)
    peaks_low <- findpeaks(-rezero(cleaned_tracks[[col_name]]), 
                           nups = 4, ndowns = 4)
  }else{
    peaks_high <- findpeaks(rezero(cleaned_tracks[[col_name]]), 
                            nups = 5, ndowns = 5, minpeakheight = min_height)
    peaks_low <- findpeaks(-rezero(cleaned_tracks[[col_name]]), 
                           nups = 5, ndowns = 5, minpeakheight =  min_height)
  }
  
  peaks_points <- data.frame("dog" = rep(dog, length(c(peaks_high,peaks_low))),
                             "run" = rep(run, length(c(peaks_high,peaks_low))),
                             "frames" = c(cleaned_tracks$frame[peaks_high[,2]],
                                          cleaned_tracks$frame[peaks_low[,2]]),
                             "times" = c(cleaned_tracks$time[peaks_high[,2]],
                                         cleaned_tracks$time[peaks_low[,2]]),
                             "peak_value" = c(cleaned_tracks[[col_name]][peaks_high[,2]],
                                              cleaned_tracks[[col_name]][peaks_low[,2]]),
                             "task" = c(cleaned_tracks$task[peaks_high[,2]],
                                        cleaned_tracks$task[peaks_low[,2]]),
                             "partition" = c(cleaned_tracks$partition[peaks_high[,2]],
                                             cleaned_tracks$partition[peaks_low[,2]])
  )
  peaks_points <- peaks_points %>% arrange(peaks_points$times)
  peaks_points <- unique(peaks_points)
  all_amps <- data.frame()
  for(j in 1:nlevels(cleaned_tracks$partition)){
    dat_subset <- peaks_points[peaks_points$partition == levels(peaks_points$partition)[j],]
    if(nrow(dat_subset) < 2 ) next
    amps <- data.frame("amp" = rep(NA, nrow(dat_subset)-1), 
                       "task" = rep(dat_subset$task[1], nrow(dat_subset)-1),
                       "dog" = rep(dog, nrow(dat_subset)-1),
                       "run" = rep(run, nrow(dat_subset)-1))
    
    for(k in 1:(nrow(dat_subset)-1)){
      amps$amp[k] <- abs(dat_subset$peak_value[k] - dat_subset$peak_value[k+1])
    }
    all_amps <- rbind(all_amps, amps)
  }
  all_signs <- data.frame()
  for(j in 1:nlevels(cleaned_tracks$partition)){
    dat_subset <- peaks_points[peaks_points$partition == levels(peaks_points$partition)[j],]
    if(nrow(dat_subset) < 2 ) next
    for(k in 2:nrow(dat_subset)){
      if(sign(dat_subset$peak_value[k-1]) != sign(dat_subset$peak_value[k])){
        sign_frame <- floor((dat_subset$frames[k] - 
                               dat_subset$frames[k-1])/2) + 
          dat_subset$frames[k-1]
        signs <- data.frame("dog" = rep(dog, nrow(dat_subset)-1),
                            "run" = rep(run, nrow(dat_subset)-1),
                            "frame" = sign_frame,
                            "time" = cleaned_tracks$time[cleaned_tracks$frame == 
                                                           sign_frame],
                            "task" = rep(dat_subset$task[1], nrow(dat_subset)-1),
                            "partition" = levels(peaks_points$partition)[j])
        all_signs <- rbind(all_signs, signs)
      }
    }
  }
  p1 <- ggplot(cleaned_tracks, aes(x = time, y = !!sym(col_name), color = task, group = partition)) + 
    geom_line(aes(color=NULL), color = "gray50") + 
    geom_point() + 
    geom_point(data = peaks_points, mapping = aes(x = times, y = peak_value), 
               col = "black") +
    
    theme_minimal()
  
  #peaks_points$peak_value <- abs(peaks_points$peak_value)
  return(list("p1" = p1,"peaks_points" = peaks_points, "amps" = all_amps, 
              "signs" = all_signs))
}

#### Constants #### 

chemical_holes <- data.frame("trial" = rep(c(1,1,2,2,3,3), each = 4), 
                             "run" = rep(c(1,2,1,2,1,2), each = 4), 
                             "hole" = rep(c(1,2,3,4),  times = 6),
                             "chemical" = factor(c("2E1H", "Blank", "bromooctane", "methyl_benzoate",
                                            "Ammonia", "Blank", "bromooctane", "methyl_benzoate",
                                            "methyl_benzoate", "bromooctane", "Blank", "2E1H",
                                            "methyl_benzoate", "bromooctane", "Blank", "Ammonia",
                                            "2E1H", "Blank", "methyl_benzoate", "bromooctane", 
                                            "bromooctane", "methyl_benzoate", "Blank", "Ammonia"), 
                                            levels = c("2E1H", "Blank", "bromooctane", "methyl_benzoate",
                                                       "Ammonia"))
)
