# Philips velocity analysis
# Andrew Robertson
# 2024-03-15
# based on "1-VEcho\scripts\CCA_velocity\CCA_imaging_velocity_VEVA_update20231011.R"

# import data extracted from MAUI (*.avi) and MATLAB (*.jpg)
# clean data
# summarize and plot data
# output beat-by-beat csv file

# to do:
# include n column in CAR_study_<vessel>v_cleaned_group_mean_csv
# include beat # instead of cine # in CAR_study_<vessel>v_cleaned_individual_mean.csv
# export mean signal as well as peak (outer envelope)
# multiply fract*rri to get time col for CAR_study_CCAv_cleaned_individual_mean (and group)

# initialize ----
library(tidyverse)
library(RColorBrewer) # colour blind palette
library(plotly)
library(signal)
library(pracma)
library(cowplot)
library(janitor)
library(here)
library(zoo)

# remove columns all na
not_all_na <- function(x) any(!is.na(x))
# not in list
`%notin%` <- Negate(`%in%`)

row_names <- NULL
col_names <- c("id", "cine", "time_s", "peak_velocity_mps", "mean_velocity_mps") 
df_combined <- matrix(NA, ncol = 5, dimnames = list(row_names, col_names)) |>  
  as_tibble(.name_repair = "check_unique")
df_combined_mean <- matrix(NA, ncol = 5, dimnames = list(row_names, col_names)) |> 
  as_tibble(.name_repair = "check_unique")
df_combined_norm_mean <- matrix(NA, ncol = 5, dimnames = list(row_names, col_names)) |> 
  as_tibble(.name_repair = "check_unique") |> 
  rename("cycle_fract" = "time_s")

video <- tibble(id = as.character(),
                visit = as.character(),
                protocol = as.character(),
                frame_rate = as.character())

beats_table <- tibble(id = as.character(),
                      visit = as.character(),
                      protocol = as.character(),
                      beat = as.character(),
                      time = as.numeric(),
                      rri_ms = as.numeric(),
                      peak_vs_cmps = as.numeric(), 
                      peak_vd_cmps = as.numeric(),
                      peak_vm_cmps = as.numeric(),
                      peak_vamp_cmps = as.numeric(),
                      peak_pi = as.numeric(),
                      peak_ri = as.numeric(),
                      mean_vs_cmps = as.numeric(), 
                      mean_vd_cmps = as.numeric(),
                      mean_vm_cmps = as.numeric(),
                      mean_vamp_cmps = as.numeric(),
                      mean_pi = as.numeric(),
                      mean_ri = as.numeric())

summary_index <- 1 # index counter for beat-by-beat data
scanner <- "philips"
vessel0 <- "VA"

dirs <- list.dirs(path = "//ahsfile.uwaterloo.ca/hughsonlab$/Hedge, Eric/2024 LBNP Study", full.names = F, recursive = F)
id_list <- grep("LBNP", dirs, value = T)

# load screening file ----
df_screen <- read_tsv(here("data", "LBNP_study_velocity_cleaning.txt")) |> 
  clean_names()

# User input ----
fn_detect <- NULL
id0 <- NULL
visit0 <- NULL
protocol0 <- NULL

while (is_empty(fn_detect)) {
  
  message(str_c("Which LBNP Study ID do you want to process:", str_flatten_comma(id_list), sep = " "))
  while (is_empty(id0)) {
    id0 <- readline(prompt="Enter 'ID' choice > ")
    # check if id exists
    if (id0 %notin% id_list) {
      message(str_c("The ID '", id0, "' does not exist. Please check and retry.", sep = ""))
      id0 <- NULL
    }  
  }
  
  message(str_c(id0, " -> which session do you want to process?\n", 
                "1 - V1\n2 - V2", sep = ""))
  while (is_empty(visit0)) {
    v_choice <- readline(prompt = "Enter 'Visit' choice: ")
    visit0 <- case_when(
      v_choice == 1 ~ "V1",
      v_choice == 2 ~ "V2",
      .default = "unknown"
    )
  }
  
  message(str_c(id0, " -> which protocol do you want to process?\n", 
                "1 - DCA\n2 - HYPOCVR\n3 - REBREATHE\n4 - LBNPINC\n5 - STEP\n6 - HALFSTEP", sep = ""))
  while (is_empty(protocol0)) {
    p_choice <- readline(prompt = "Enter 'Protocol' choice: ")
    protocol0 <- case_when(
      p_choice == 1 ~ "DCAplus2BREATH",
      p_choice == 2 ~ "HYPOCVR",
      p_choice == 3 ~ "REBREATHE",
      p_choice == 4 ~ "LBNPINC",
      p_choice == 5 ~ "STEP",
      p_choice == 6 ~ "HALFSTEP",
      .default = "unknown"
    )
  }
  
  if (protocol0 == "STEP"){
    protocol0 <- "_STEP"
  } 
  
  # looks for file with selected parameters
  fn_detect <- list.files(path = str_c("//ahsfile.uwaterloo.ca/hughsonlab$/Hedge, Eric/2024 LBNP Study", 
                                       id0, visit0, "Philips", sep = "/"), 
                          pattern = str_c(protocol0, ".csv", sep = ".+"), 
                          full.names = T, ignore.case = T) #find file on research drive
  
  # check if file exists
  if (is_empty(fn_detect)) {
    message("Could not identify file with these settings. Please check and retry.")
    id0 <- NULL
    visit0 <- NULL
    protocol0 <- NULL
  }  
  
  if (length(fn_detect)>1){
    message("Found more than one file with these settings. Please check filenames in directory and retry.")
    id0 <- NULL
    visit0 <- NULL
    protocol0 <- NULL
  }
  # loop until found file
}

if (protocol0 == "_STEP"){
  protocol0 <- "STEP"
} 

id_visit_protocol <- str_c(id0, visit0, protocol0, sep = " ")

message(str_c("Processing", id_visit_protocol, "...", sep = " "))

# check csv file for waveform extraction software; for LBNP study, this will likely always be MAUI
extraction_software <- df_screen |>  
  dplyr::filter(id == id0, session == visit0, protocol == protocol0) |>  
  pull(extraction)

# import data ----

if (extraction_software == "MAUI") {
  
  df_cine <- read_csv(fn_detect, na = c("", "nan")) |> 
    clean_names()
  
  # keep a raw file to help with trouble shooting if needed
  df_cine_raw <- df_cine
  
  # trim dead air frames off the end
  df_cine$rolling_zeros <- rollmean(df_cine$peak_positive_cm_s, k = 200, fill = 0, align = "left")
  zeros <- which(df_cine$rolling_zeros < 5)
  zeros_end <- min(which(zeros > nrow(df_cine)*.95))
  df_cine <- df_cine[1:zeros[zeros_end],]
  
  print(plot_ly(type = "scatter",
                mode = "lines") |> 
          add_lines( x = as.numeric(row.names(df_cine)),
                     y = df_cine$peak_positive_cm_s,
                     line = list(color = "black"),
                     name = "raw") |> 
          add_lines( x = as.numeric(row.names(df_cine)),
                     y = df_cine$rolling_zeros,
                     line = list(color = "red"),
                     name = "rolling mean")|> 
          add_segments(x = as.numeric(row.names(df_cine)[zeros[zeros_end]]), xend = as.numeric(row.names(df_cine)[zeros[zeros_end]]), 
                       y = 0, yend = max(df_cine$peak_positive_cm_s, na.rm = T),
                       line = list(colour = "darkorange", width = 2),
                       name = "end"))
  
  # check to see if final frame includes all 4 pixels; if only partial frame remains, delete remainder
  frame_check <- unique(df_cine$frame_number[(nrow(df_cine)-3):nrow(df_cine)])
  if (length(frame_check > 1)){
    partial_frame <- tail(frame_check,1)
    df_cine <- df_cine |> dplyr::filter(frame_number != partial_frame)
  }
  
  # cycle through frames
  frame_num <- unique(df_cine$frame_number)
  
  # adjust time points within each frame equidistant based on time between frames
  # each frame has 4 data points (X locations)
  # set time for each remaining data point in frame to 1/4 duration of frame 
  
  df_cine$time_adj <- NA
  
  for (f in 1:(length(frame_num)-1)) {
    frame_now <- which(df_cine$frame_number == frame_num[f])
    frame_next <- which(df_cine$frame_number == frame_num[f+1])
    frame_time <- df_cine$time_seconds[(frame_next[1])] - df_cine$time_seconds[frame_now[1]]
    df_cine$time_adj[frame_now] <- c(df_cine$time_seconds[frame_now[1]], 
                                     (df_cine$time_seconds[frame_now[1]] + frame_time*1/4), 
                                     (df_cine$time_seconds[frame_now[1]] + frame_time*2/4), 
                                     (df_cine$time_seconds[frame_now[1]] + frame_time*3/4))
    frame_xdist <- mean(diff(df_cine$x_location[frame_now]), na.rm = F)
    
    # x location within a frame is usually less than 4 pixels; 
    # remove frames where xLocation spread > 10; this is usually when the sweep reaches the 
    # end of the screen and flips back to the beginning for another sweep
    
    if (frame_xdist > 10) df_cine[frame_now, 5:16] <- NA
  }
  
  # adjust time points for final frame
  frame_now <- which(df_cine$frame_number == frame_num[length(frame_num)])
  df_cine$time_adj[frame_now] <- c(df_cine$time_seconds[frame_now[1]], 
                                   (df_cine$time_seconds[frame_now[1]] + frame_time*1/3), 
                                   (df_cine$time_seconds[frame_now[1]] + frame_time*2/4), 
                                   (df_cine$time_seconds[frame_now[1]] + frame_time*3/4))
  frame_xdist <- mean(diff(df_cine$x_location[frame_now]), na.rm = F)
  if (frame_xdist > 10) df_cine[frame_now, 5:16] <- NA
  
  # remove bad sections where xLocation extraction within frame is identical
  rm_duplicates <- which(duplicated(df_cine[c(2,4)]))
  rm_duplicates_first <- rm_duplicates-1
  rm_dupicates_merge <- sort(unique(vctrs::vec_c(rm_duplicates, rm_duplicates_first)))
  df_cine[rm_dupicates_merge, 5:16] <- NA
  
  # remove section where xLocation >900 (into velocity scale)
  #scale_interference <- which(df_cine$x_location >= 1000)
  #df_cine[scale_interference, 5:16] <- NA
  
  df_clean <- df_cine[2:18] |>
    drop_na() |>
    select(peak_positive_cm_s, mean_positive_cm_s, time_adj) |>
    rename(time_s = time_adj) |> 
    relocate(time_s)
  
  # record the natural sampling rate from raw data
  Fs <- round(1/(frame_time/4)) #sampling frequency
  video <- list(id0, visit0, protocol0, as.character(round(1/frame_time, digits = 1)))
  
} else {
  
  # Non-MAUI extraction
  
  df_cine <- read_csv(velocity_files[r], na = c("", "nan")) |> 
    clean_names() |> 
    rename(time_s = time, 
           peak_positive_cm_s = velocity_peak)
  
  df_cine_na_removed <- na.omit(df_cine)
  
  # interpolate through cleaned data
  df_cine$peak_positive_cm_s_clean <- approx(x = df_cine_na_removed$time_s, 
                                             y = df_cine_na_removed$peak_positive_cm_s,
                                             xout = df_cine$time_s,
                                             method = "linear", rule = 2)$y
  
  df_cine |> ggplot(aes(x = time_s)) +
    geom_line(aes(y = peak_positive_cm_s), colour = "black", linewidth = 2) +
    geom_line(aes(y = peak_positive_cm_s_clean), colour = "red", linewidth = 1) +
    scale_y_continuous(name = "Blood Velocity (cm/s)") +
    scale_x_continuous(name = "Time (s)") +
    theme_classic()
  
  df_clean <- df_cine |> 
    select(time_s, peak_positive_cm_s)
  
  # record the natural sampling rate from raw data
  Fs <- round(1/mean(diff(df_cine$time_s)))
  video <- list(id0, rep0, "image")
}

# set DFV threshold
# import recorded thresholds
mean_dfv_thresh <- df_screen |>  
  dplyr::filter(id == id0, session == visit0, protocol == protocol0) |>  
  pull(mean_dfv)

peak_dfv_thresh <- df_screen |>  
  dplyr::filter(id == id0, session == visit0, protocol == protocol0) |>  
  pull(peak_dfv)

peak_sfv_thresh <- df_screen |>  
  dplyr::filter(id == id0, session == visit0, protocol == protocol0) |>  
  pull(peak_sfv)

hypocapnic_peak <- df_screen |>  
  dplyr::filter(id == id0, session == visit0, protocol == protocol0) |>  
  pull(hypocapnic_peak)

# estimate max threshold based on first 2000 datapoints - assume these first few beats are clean
estimate_max <- round(max(df_clean$peak_positive_cm_s[1:2000])/10)*10

threshold_df <- tibble(x = 10,
                       y = c(5,10,15,20, (estimate_max-5), estimate_max, (estimate_max+5), (estimate_max+10)))

if (protocol0 %in% c("HYPOCVR", "REBREATHE")){
  
  protocol_start <- df_screen |>  
    dplyr::filter(id == id0, session == visit0, protocol == protocol0) |>  
    pull(start_protocol)
  
  co2_start <- df_screen |>  
    dplyr::filter(id == id0, session == visit0, protocol == protocol0) |>  
    pull(hypocapnic_start) - protocol_start
  
  co2_stop <- df_screen |>  
    dplyr::filter(id == id0, session == visit0, protocol == protocol0) |>  
    pull(hypocapnic_stop) - protocol_start
  
  print(df_clean |> 
          plot_ly(x = ~time_s, y = ~peak_positive_cm_s) |> 
          add_lines(line = list(color = "black")) |> 
          layout(shapes = list(type = "rect",
                               fillcolor = "skyblue",
                               x0 = co2_start, x1 = co2_stop, xref = "x",
                               y0 = -5, y1 = max(df_clean$peak_positive_cm_s), yref = "y",
                               layer = "below")) |> 
          add_segments(x = min(df_clean$time_s), xend = max(df_clean$time_s), 
                       y = threshold_df$y, yend = threshold_df$y))
  
  message(str_c("The recorded peak dfv threshold is", peak_dfv_thresh, sep = " "))
  change_thresh <- readline("Based on the plot of peak velocity (outside of shaded hypocapnic period). Do you wish to change this threshold (y/n)? > ")
  if (change_thresh == "Y" | change_thresh == "y") peak_dfv_thresh <- as.numeric(readline("What is estimated DFV threshold for peak signal? > "))
  
  message(str_c("The recorded peak sfv threshold is", peak_sfv_thresh, sep = " "))
  change_thresh <- readline("Based on the plot of peak velocity (outside of shaded hypocapnic period). Do you wish to change this threshold (y/n)? > ")
  if (change_thresh == "Y" | change_thresh == "y") peak_sfv_thresh <- as.numeric(readline("What is estimated SFV threshold for peak signal? > "))
  
  message(str_c("The recorded peak sfv threshold (within the shaded hypocapnic period) is", hypocapnic_peak, sep = " "))
  change_thresh <- readline("Based on the plot of peak velocity (within the shaded hypocapnic period). Do you wish to change this threshold (y/n)? > ")
  if (change_thresh == "Y" | change_thresh == "y") hypocapnic_peak <- as.numeric(readline("What is estimated SFV threshold for peak signal (within the shaded hypocapnic period)? > "))
  
} else {
  
  print(df_clean |> 
          plot_ly(x = ~time_s, y = ~peak_positive_cm_s) |> 
          add_lines(line = list(color = "black")) |> 
          add_segments(x = min(df_clean$time_s), xend = max(df_clean$time_s), 
                       y = threshold_df$y, yend = threshold_df$y))
  
  message(str_c("The recorded peak dfv threshold is", peak_dfv_thresh, sep = " "))
  change_thresh <- readline("Based on the plot of peak velocity. Do you wish to change this threshold (y/n)? > ")
  if (change_thresh == "Y" | change_thresh == "y") peak_dfv_thresh <- as.numeric(readline("What is estimated DFV threshold for peak signal? > "))
  
  message(str_c("The recorded peak sfv threshold is", peak_sfv_thresh, sep = " "))
  change_thresh <- readline("Based on the plot of peak velocity. Do you wish to change this threshold (y/n)? > ")
  if (change_thresh == "Y" | change_thresh == "y") peak_sfv_thresh <- as.numeric(readline("What is estimated SFV threshold for peak signal? > "))
  
}

if (protocol0 %in% c("HYPOCVR", "REBREATHE")){
  
  print(df_clean |> 
          plot_ly(x = ~time_s, y = ~mean_positive_cm_s) |> 
          add_lines(line = list(color = "black")) |> 
          layout(shapes = list(type = "rect",
                               fillcolor = "skyblue",
                               x0 = co2_start, x1 = co2_stop, xref = "x",
                               y0 = -5, y1 = max(df_clean$mean_positive_cm_s), yref = "y",
                               layer = "below")) |> 
          add_segments(x = min(df_clean$time_s), xend = max(df_clean$time_s), 
                       y = threshold_df$y[1:4], yend = threshold_df$y[1:4]))
  
  message(str_c("The recorded mean dfv threshold is", mean_dfv_thresh, sep = " "))
  change_thresh <- readline("Based on the plot of mean velocity (outside of shaded hypocapnic period). Do you wish to change this threshold (y/n)? > ")
  if (change_thresh == "Y" | change_thresh == "y") meak_dfv_thresh <- as.numeric(readline("What is estimated DFV threshold for mean signal? > "))
  
} else {
  
  print(df_clean |> 
          plot_ly(x = ~time_s, y = ~mean_positive_cm_s) |> 
          add_lines(line = list(color = "black")) |> 
          add_segments(x = min(df_clean$time_s), xend = max(df_clean$time_s), 
                       y = threshold_df$y, yend = threshold_df$y))
  
  message(str_c("The recorded mean dfv threshold is", mean_dfv_thresh, sep = " "))
  change_thresh <- readline("Based on the plot of mean velocity. Do you wish to change this threshold (y/n)? > ")
  if (change_thresh == "Y" | change_thresh == "y") meak_dfv_thresh <- as.numeric(readline("What is estimated DFV threshold for mean signal? > "))
  
}


# identify which velocity is below the dfv threshold and change to NA (exception last 200 data points)
if (protocol0 %in% c("HYPOCVR", "REBREATHE")){
  
  #exclude hypocapnic range from cleaning
  t1 <- max(which(df_clean$time_s < co2_start))
  t2 <- min(which(df_clean$time_s > co2_stop))
  
  peak_dfv_noise <- which(df_clean$peak_positive_cm_s[1:t1] < peak_dfv_thresh)
  peak_sfv_noise <- which(df_clean$peak_positive_cm_s[1:t1] > peak_sfv_thresh)
  mean_dfv_noise <- which(df_clean$mean_positive_cm_s[1:t1] < mean_dfv_thresh)
  
  peak_dfv_noise2 <- which(df_clean$peak_positive_cm_s[t2:(nrow(df_clean)-200)] < peak_dfv_thresh) + t2 -1 
  peak_sfv_noise2 <- which(df_clean$peak_positive_cm_s[t2:(nrow(df_clean)-200)] > peak_sfv_thresh) + t2 -1 
  mean_dfv_noise2 <- which(df_clean$mean_positive_cm_s[t2:(nrow(df_clean)-200)] < mean_dfv_thresh) + t2 -1 
  
  #inside hypocapnic range
  peak_sfv_noise3 <- which(df_clean$peak_positive_cm_s[t1:t2] > hypocapnic_peak) + t1 -1 
  
  noise <- sort(unique(vctrs::vec_c(peak_dfv_noise, peak_sfv_noise, mean_dfv_noise,
                                    peak_dfv_noise2, peak_sfv_noise2, mean_dfv_noise2,
                                    peak_sfv_noise3)))
  
} else {
  
  peak_dfv_noise <- which(head(df_clean$peak_positive_cm_s, (nrow(df_clean)-200)) < peak_dfv_thresh)
  peak_sfv_noise <- which(head(df_clean$peak_positive_cm_s, (nrow(df_clean)-200)) > peak_sfv_thresh)
  mean_dfv_noise <- which(head(df_clean$mean_positive_cm_s, (nrow(df_clean)-200)) < mean_dfv_thresh)
  
  noise <- sort(unique(vctrs::vec_c(peak_dfv_noise, peak_sfv_noise, mean_dfv_noise)))
  
}

df_clean <- df_clean |> 
  mutate(mean_positive_cm_s_thresh = mean_positive_cm_s,
         peak_positive_cm_s_thresh = peak_positive_cm_s)
df_clean$mean_positive_cm_s_thresh[noise] <- NA
df_clean$peak_positive_cm_s_thresh[noise] <- NA

# interpolate through NA (data removed by diastolic threshold)
df_clean$mean_positive_cm_s_clean <- na.approx(df_clean$mean_positive_cm_s_thresh, na.rm = F)
df_clean$peak_positive_cm_s_clean <- na.approx(df_clean$peak_positive_cm_s_thresh, na.rm = F)

# trim if ends with trailing NA.
df_clean <- df_clean |> drop_na(peak_positive_cm_s_clean)

# alternate cleaning methods (e.g., spline)
# df_cine200$peak_velocity_mps <- spline(x = df_clean$time_s, y = df_clean$peak_velocity_mps, xout = ts, method = "periodic")$y
# df_cine200$peak_velocity_mps <- ifelse(df_cine200$peak_velocity_mps < 0, 0, df_cine200$peak_velocity_mps)

# resample to 200 Hz and interpolate missing data (data removed with dfv theshold)
Fs_new <- 200

ts <- seq(round(min(df_clean$time_s), digits = 2), round(max(df_clean$time_s), digits = 2), 1/Fs_new)

df_cine200 <- matrix(NA, nrow = length(ts), ncol = 3, dimnames = list(NULL, c("time_s", "peak_velocity_cmps", "mean_velocity_cmps"))) |> 
  as_tibble(.name_repair = "check_unique") |> 
  mutate(time_s = ts)

df_cine200$peak_velocity_cmps <- approx(df_clean$time_s, 
                                        df_clean$peak_positive_cm_s_clean, 
                                        df_cine200$time_s, 
                                        method = "linear", 
                                        rule = 2)$y

df_cine200$mean_velocity_cmps <- approx(df_clean$time_s, 
                                        df_clean$mean_positive_cm_s_clean, 
                                        df_cine200$time_s, 
                                        method = "linear", 
                                        rule = 2)$y

if (protocol0 %in% c("HYPOCVR", "REBREATHE")){
  
  print(plot_ly(type = "scatter",
                mode = "lines") |> 
          add_trace(x = df_cine$time_adj,
                    y = df_cine$peak_positive_cm_s,
                    line = list(color = "red"),
                    name = "raw") |> 
          add_trace(x = df_cine200$time_s,
                    y = df_cine200$peak_velocity_cmps,
                    line = list(color = "dodgerblue"),
                    name = "cleaned") |> 
          layout(shapes = list(type = "rect",
                               fillcolor = "skyblue", colour = "red",
                               x0 = co2_start, x1 = co2_stop, xref = "x",
                               y0 = -5, y1 = max(df_clean$mean_positive_cm_s), yref = "y",
                               layer = "below")))
  
} else {
  
  print(plot_ly(type = "scatter",
                mode = "lines") |> 
          add_trace( x = df_cine$time_adj,
                     y = df_cine$peak_positive_cm_s,
                     line = list(color = "red")) |> 
          add_trace( x = df_cine200$time_s,
                     y = df_cine200$peak_velocity_cmps,
                     line = list(color = "dodgerblue")))
  
}

#import ecg and time sync ----

if (protocol0 == "STEP"){
  protocol0 <- "_STEP"
} 

ecg_detect <- list.files(path = str_c("//ahsfile.uwaterloo.ca/hughsonlab$/Hedge, Eric/2024 LBNP Study", 
                                      id0, visit0, "ECG", sep = "/"), 
                         pattern = str_c(protocol0, "ecg_200Hz.txt", sep = ".+"), 
                         full.names = T, ignore.case = T) #find file on research drive

if (protocol0 == "_STEP"){
  protocol0 <- "STEP"
} 
df_ecg <- read_tsv(ecg_detect, col_names = c("time_s", "ecg", "comment"), skip = 6, col_types = "ddc")
ecg_sync <- df_screen |> 
  dplyr::filter(id == id0, session == visit0, protocol == protocol0) |>  
  pull(labchart_sync_time)

last_spike_va <- max(which(df_cine200$peak_velocity_cmps == max(tail(df_cine200$peak_velocity_cmps, 200))))
time_diff <- ecg_sync - df_cine200$time_s[last_spike_va]
df_cine200$time_s <- df_cine200$time_s + time_diff

filter_choice <- F
if (filter_choice == T){
  
  # design a 4th order 20Hz low pass filter to match filter on LabChart hemodynamic data
  Fc <- 20 # cut off frequency
  n <- 4 # order of filter
  W <- Fc/(Fs/2) # Normalized frequency
  filter_design <- butter(n,W,'low')
  
  df_cine200$peak_velocity_cmps_filt20 <- filtfilt(filter_design$b, filter_design$a, df_cine200$peak_velocity_cmps)
  df_cine200$mean_velocity_cmps_filt20 <- filtfilt(filter_design$b, filter_design$a, df_cine200$mean_velocity_cmps)
  
  df_ecg <- df_ecg |> 
    mutate(time_s = round(time_s, 3)) |> 
    left_join(df_cine200 |>
                mutate(time_s = round(time_s, 3)) |> 
                dplyr::select(time_s, peak_velocity_cmps_filt20, mean_velocity_cmps_filt20),
              by = "time_s") |> 
    relocate("comment", .after = "mean_velocity_cmps_filt20")|> 
    drop_na(peak_velocity_cmps_filt20)
  
} else {
  
  message("Exporting non-filtered data at 200Hz.")
  df_ecg <- df_ecg |> 
    mutate(time_s = round(time_s, 3)) |> 
    left_join(df_cine200 |>
                mutate(time_s = round(time_s, 3)) |> 
                dplyr::select(time_s, peak_velocity_cmps, mean_velocity_cmps),
              by = "time_s") |> 
    relocate("comment", .after = "mean_velocity_cmps") |> 
    drop_na(peak_velocity_cmps)
  
}

# Create new export header ----
# extract header info from original labchart file add CHI headers for new data file
data_header <- read.table(ecg_detect, sep="\t", header=F, fill=T, nrows = 6)[1:2]
data_units <- c("UnitName=", "V")
data_header <- rbind(data_header, data_units)
data_header$V4 <-  c("", "", "", "", "Peak_VABV", "10.000V", "cm/s")
data_header$V5 <-  c("", "", "", "", "Mean_VABV", "10.000V", "cm/s")
data_header$V6 <-  c("", "", "", "", "", "", "")

# Output new labchart textfile ----
fn_out <- here("output", str_c(id0, visit0, protocol0, "VA.txt", sep = "_"))
write_delim(data_header, fn_out, "\t", col_names = F)
write_delim(df_ecg, fn_out, "\t", append = T, col_names = F, na = "")

message("Complete")
