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

# initialize plots
plotlist_trace <- list() # plot full sweep

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
df_screen <- read_tsv(here("data", "LBNP_study_velocity_peakcounts.txt")) |> 
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

combined_beats_peak <- matrix(as.numeric(), nrow = 500, ncol = 1) |> 
  as_tibble(.name_repair = "minimal")

id_visit_protocol <- str_c(id0, visit0, protocol0, sep = " ")

message(str_c("Processing", id_visit_protocol, "...", sep = " "))

# check csv file to see if we want to skip cine
omit_cine <- df_screen |>  
  dplyr::filter(id == id0, session == visit0, protocol == protocol0) |>  
  pull(skip)

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
  # this noise interferes with foot detection
  df_cine$rolling_zeros <- rollmean(df_cine$peak_positive_cm_s, k = 100, fill = 0, align = "left")
  zeros <- which(df_cine$rolling_zeros < 5)
  df_cine <- df_cine[1:min(zeros),]
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
           peak_velocity_mps = velocity_peak)
  
  df_cine_na_removed <- na.omit(df_cine)
  
  # interpolate through cleaned data
  df_cine$peak_velocity_mps_clean <- approx(x = df_cine_na_removed$time_s, 
                                            y = df_cine_na_removed$peak_velocity_mps,
                                            xout = df_cine$time_s,
                                            method = "linear", rule = 2)$y
  
  df_cine |> ggplot(aes(x = time_s)) +
    geom_line(aes(y = peak_velocity_mps), colour = "black", linewidth = 2) +
    geom_line(aes(y = peak_velocity_mps_clean), colour = "red", linewidth = 1) +
    scale_y_continuous(name = "Blood Velocity (cm/s)") +
    scale_x_continuous(name = "Time (s)") +
    theme_classic()
  
  df_clean <- df_cine |> 
    select(time_s, peak_velocity_mps)
  
  # record the natural sampling rate from raw data
  Fs <- round(1/mean(diff(df_cine$time_s)))
  video <- list(id0, rep0, "image")
}


# resample to 200 Hz and interpolate missing data (data removed with dfv theshold)
Fs_new <- 200

ts <- seq(round(min(df_clean$time_s), digits = 2), round(max(df_clean$time_s), digits = 2), 1/Fs_new)

df_cine200 <- matrix(NA, nrow = length(ts), ncol = 3, dimnames = list(NULL, c("time_s", "peak_velocity_cmps", "mean_velocity_cmps"))) |> 
  as_tibble(.name_repair = "check_unique") |> 
  mutate(time_s = ts)

df_cine200$peak_velocity_cmps <- approx(df_clean$time_s, 
                                        df_clean$peak_positive_cm_s, 
                                        df_cine200$time_s, 
                                        method = "linear", 
                                        rule = 2)$y

df_cine200$mean_velocity_cmps <- approx(df_clean$time_s, 
                                        df_clean$mean_positive_cm_s, 
                                        df_cine200$time_s, 
                                        method = "linear", 
                                        rule = 2)$y

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

# estimate max threshold based on first 2000 datapoints - assume these first few beats are clean
estimate_max <- round(max(df_cine200$peak_velocity_cmps[1:2000])/10)*10

threshold_df <- tibble(x = 10,
                       y = c(5, 10, 15, 20, (estimate_max-5), estimate_max, (estimate_max+5), (estimate_max+10)))

df_cine200 |> ggplot() +
  geom_line(aes(x = time_s, y = peak_velocity_cmps), colour = "black", linewidth = 0.75) +
  geom_hline(data = threshold_df, aes(yintercept = y), colour = "red", linewidth = 1) +
  geom_text(data = threshold_df, aes(x = x, y = y-2.5, label = paste0('<b>',y,'</b>')), size = 4, colour = "cyan4", fontface = "bold") +
  theme_bw()
ggplotly()

message(str_c("The recorded peak dfv threshold is", peak_dfv_thresh, sep = " "))
change_thresh <- readline("Based on the plot of peak velocity. Do you wish to change this threshold (y/n)? > ")
if (change_thresh == "Y" | change_thresh == "y") peak_dfv_thresh <- as.numeric(readline("What is estimaged DFV threshold for peak signal? > "))

message(str_c("The recorded peak sfv threshold is", peak_sfv_thresh, sep = " "))
change_thresh <- readline("Based on the plot of peak velocity. Do you wish to change this threshold (y/n)? > ")
if (change_thresh == "Y" | change_thresh == "y") peak_sfv_thresh <- as.numeric(readline("What is estimaged SFV threshold for peak signal? > "))

df_cine200 |> ggplot() +
  geom_line(aes(x = time_s, y = mean_velocity_cmps), colour = "black", linewidth = 0.75) +
  geom_hline(data = threshold_df[1:4,], aes(yintercept = y), colour = "red", linewidth = 1) +
  geom_text(data = threshold_df[1:4,], aes(x = x, y = y-2.5, label = paste0('<b>',y,'</b>')), size = 4, colour = "cyan4", fontface = "bold") +
  theme_bw()
ggplotly()

message(str_c("The recorded mean dfv threshold is", mean_dfv_thresh, sep = " "))
change_thresh <- readline("Based on the plot of mean velocity. Do you wish to change this threshold (y/n)? > ")
if (change_thresh == "Y" | change_thresh == "y") peak_dfv_thresh <- as.numeric(readline("What is estimaged DFV threshold for mean signal? > "))

# identify which velocity is below the dfv threshold and change to NA (exception last 500 data points)
peak_dfv_noise <- which(head(df_cine200$peak_velocity_cmps, (nrow(df_cine200)-500)) < peak_dfv_thresh)
peak_sfv_noise <- which(head(df_cine200$peak_velocity_cmps, (nrow(df_cine200)-500)) > peak_sfv_thresh)
mean_dfv_noise <- which(head(df_cine200$mean_velocity_cmps, (nrow(df_cine200)-500)) < mean_dfv_thresh)

noise <- unique(vctrs::vec_c(peak_dfv_noise, peak_sfv_noise, mean_dfv_noise))

df_cine200 <- df_cine200 |> 
  mutate(mean_velocity_cmps_thresh = mean_velocity_cmps,
         peak_velocity_cmps_thresh = peak_velocity_cmps)
df_cine200$mean_velocity_cmps_thresh[noise] <- NA
df_cine200$peak_velocity_cmps_thresh[noise] <- NA

# interpolate through NA (data removed by diastolic threshold)
df_cine200$mean_velocity_cmps_clean <- na.approx(df_cine200$mean_velocity_cmps_thresh, na.rm = F)
df_cine200$peak_velocity_cmps_clean <- na.approx(df_cine200$peak_velocity_cmps_thresh, na.rm = F)

# trim if ends with trailing NA.
df_cine200 <- df_cine200 |> drop_na(peak_velocity_cmps_clean)

# alternate cleaning methods (e.g., spline)
# df_cine200$peak_velocity_mps <- spline(x = df_clean$time_s, y = df_clean$peak_velocity_mps, xout = ts, method = "periodic")$y
# df_cine200$peak_velocity_mps <- ifelse(df_cine200$peak_velocity_mps < 0, 0, df_cine200$peak_velocity_mps)

df_cine200 |> ggplot() +
  geom_line(aes(x = time_s, y = peak_velocity_cmps), colour = "dodgerblue", linewidth = 1.5) +
  geom_line(aes(x = time_s, y = peak_velocity_cmps_clean), colour = "red", linewidth = 0.5) +
  theme_classic(); print(ggplotly())

#import ecg and time sync ----
ecg_detect <- list.files(path = str_c("//ahsfile.uwaterloo.ca/hughsonlab$/Hedge, Eric/2024 LBNP Study", 
                                      id0, visit0, sep = "/"), 
                         pattern = str_c(protocol0, "ecg.txt", sep = ".+"), 
                         full.names = T, ignore.case = T) #find file on research drive
df_ecg <- read_tsv(ecg_detect, col_names = c("time_s", "ecg", "comment"), skip = 6, col_types = "ddc")
ecg_sync <- df_screen |> 
  dplyr::filter(id == id0, session == visit0, protocol == protocol0) |>  
  pull(labchart_sync_time)

last_spike_va <- which.max(df_cine200$peak_velocity_cmps == max(tail(df_cine200$peak_velocity_cmps, 1000)))
time_diff <- ecg_sync - df_cine200$time_s[last_spike_va]
df_cine200$time_s <- df_cine200$time_s + time_diff

# nfeet <- as.numeric(readline("How many waveform feet should be detected? > "))
nfeet <- df_screen |>  
  dplyr::filter(id == id0, session == visit0, protocol == protocol0) |>  
  pull(foot_count)

npeaks <- findpeaks(df_cine200$peak_velocity_cmps_clean, minpeakheight = peak_sfv_thresh/2)
nfeet <- round(nrow(npeaks)+.1*nrow(npeaks))

# design a 4th order 20Hz low pass filter to match filter on LabChart hemodynamic data
Fc <- 20 # cut off frequency
n <- 4 # order of filter
W <- Fc/(Fs/2) # Normalized frequency
filter_design <- butter(n,W,'low')

df_cine200$peak_velocity_cmps_filt20 <- filtfilt(filter_design$b, filter_design$a, df_cine200$peak_velocity_cmps_clean)
df_cine200$mean_velocity_cmps_filt20 <- filtfilt(filter_design$b, filter_design$a, df_cine200$mean_velocity_cmps_clean)

df_ecg <- df_ecg |> 
  mutate(time_s = round(time_s, 3)) |> 
  left_join(df_cine200 |>
              mutate(time_s = round(time_s, 3)) |> 
              dplyr::select(time_s, peak_velocity_cmps_filt20, mean_velocity_cmps_filt20),
            by = "time_s") |> 
  relocate("comment", .after = "mean_velocity_cmps_clean")

# Create new export header ----
# extract header info from original labchart file add CHI headers for new data file
data_header <- read.table(ecg_detect, sep="\t", header=F, fill=T, nrows = 6)[1:2]
data_header$V4 <-  c("", "", "", "", "Peak_VABV", "10.000V")
data_header$V5 <-  c("", "", "", "", "Mean_VABV", "10.000V")
data_header$V6 <-  c("", "", "", "", "", "")

# Output new data files ----
# new labchart textfile
fn_out <- here("output", str_c(id0, visit0, protocol0, "VA_cleaned.txt", sep = "_"))
write_delim(data_header, fn_out, "\t", col_names = F)
write_delim(df_ecg, fn_out, "\t", append = T, col_names = F, na = "")

## foot detection ----

# design a 4th order 8Hz low pass filter for foot detection
Fc <- 8 # cut off frequency
n <- 4 # order of filter
W <- Fc/(Fs/2) # Normalized frequency
filter_design <- butter(n,W,'low')

df_cine200$peak_velocity_cmps_filt8 <- filtfilt(filter_design$b, filter_design$a, df_cine200$peak_velocity_cmps_clean)

# find 2nd derivative peak for foot detection
diff_lag = 1

df_cine200 <- df_cine200 |> 
  mutate(peak_velocity_cmps_dif1 = vctrs::vec_c(rep(NA, 1*diff_lag), diff(df_cine200$peak_velocity_cmps_filt8, diff_lag, 1)),
         peak_velocity_cmps_dif2 = vctrs::vec_c(rep(NA, 2*diff_lag), diff(df_cine200$peak_velocity_cmps_filt8, diff_lag, 2)))

rri_est <- (max(df_cine200$time_s)-min(df_cine200$time_s))/nfeet # in sec

wf_feet <- findpeaks(df_cine200$peak_velocity_cmps_dif2, minpeakheight = max(df_cine200$peak_velocity_cmps_dif2, na.rm = T)/7) |> 
  as_tibble(.name_repair = "universal") |> 
  arrange(...2) |> 
  mutate(time_s = df_cine200$time_s[...2],
         peak_velocity_cmps = df_cine200$peak_velocity_cmps_clean[...2])

# remove 1st peak if detected during diastole of 1st beat
if (id_visit_protocol %in% c("LBNP01")) {
  wf_feet <- wf_feet[-1, ]
} 
if (id_visit_protocol %in% c("LBNP01_V1_NA")) {
  wf_feet <- wf_feet[-c(1:2), ]
} 

# check for double peak and remove 2nd peak
rr_dist <- rri_est*Fs_new #with buffer
peak_check = 2
#while (nrow(wf_feet) > nfeet){
while (wf_feet$time_s[peak_check] < tail(wf_feet$time_s, 1)){ 
  if ((wf_feet$...2[peak_check] - wf_feet$...2[(peak_check-1)]) < rr_dist){
    wf_feet <- wf_feet[-peak_check,]
  } else {
    peak_check <- peak_check + 1
  }
}

# confirm peak detection
df_cine200 |> ggplot() +
  geom_line(aes(x = time_s, y = peak_velocity_cmps_clean), colour = "black") +
  geom_path(aes(x = time_s, y = peak_velocity_cmps_dif2), colour = "darkcyan") +
  geom_vline(data = wf_feet, aes(xintercept = time_s), colour = "red", linewidth = .5) +
  labs(title = id_visit_protocol) +
  theme_bw(); ggplotly()

# find foot tangent (omit 1st and last beat)
wf_feet$vel_min <- NA
wf_feet$tangfoot <- NA
wf_feet$time_tang <- NA
wf_feet$vel_tang <- NA

for (i in 1:(nrow(wf_feet))){
  
  if (wf_feet$...2[i] <= 20) {
    wf_feet$vel_min[i] <- min(df_cine200$peak_velocity_cmps_clean[1:wf_feet$...2[i]], na.rm = T)
  } else {
    wf_feet$vel_min[i] <- min(df_cine200$peak_velocity_cmps_clean[(wf_feet$...2[i]-20):wf_feet$...2[i]], na.rm = T)
  }
  
  slope <- df_cine200$peak_velocity_cmps_dif1[wf_feet$...2[i]]
  
  if (slope <= 0) {
    slope <- df_cine200$peak_velocity_cmps_dif1[(wf_feet$...2[i]+2)]
  } 
  
  wf_feet$tangfoot[i] <- wf_feet$...2[i] - floor((wf_feet$peak_velocity_cmps[i] - wf_feet$vel_min[i])/slope) # x0 = x1 - (y1-y0)/m
  
}

wf_feet$time_tang <- df_cine200$time_s[wf_feet$tangfoot]
wf_feet$vel_tang <- df_cine200$peak_velocity_cmps_clean[wf_feet$tangfoot]

df_cine200 |> ggplot() +
  geom_path(aes(x = time_s, y = peak_velocity_cmps_filt20), colour = "black") +
  geom_point(data = wf_feet, aes(x = time_s, y = peak_velocity_cmps), colour = "blue") +
  geom_point(data = wf_feet, aes(x = time_tang, y = vel_tang), size = 2, shape = 21, fill = "red") +
  #scale_y_continuous(limits = c(0,1.5), expand = c(0,0)) +
  labs(title = id_visit_protocol) +
  theme_bw()
ggplotly()

# check <- readline("Check waveform foot detection.")

## separate beats ----

# check screen file to omit specific beats
skip <- df_screen |>  
  dplyr::filter(id == id0, session == visit0, protocol == protocol0) |>  
  pull(omit)

if (is_empty(skip)) {
  
  skip_beats <- as.numeric(NA)
  
} else {
  
  if (str_detect(skip, "-") & !is.na(skip)) { # checks if range marker exists
    
    if (str_detect(skip, ",")) { # checks if range and comma delineation exists
      
      skip_temp <- str_split(skip, ", ", simplify = T)
      range_finder <- which(str_detect(skip_temp, "-"))
      skip_expand <- as.numeric(str_split(skip_temp[range_finder], "-", simplify = T)) |> sort()
      skip_beats <- as.numeric(skip_temp[-range_finder])
      
      for (x in seq_along(range_finder)){
        
        skip_expansion <- seq(skip_expand[x*2-1],skip_expand[x*2])
        skip_beats <- vctrs::vec_c(skip_beats, skip_expansion) |> sort()
        
      }
      
    } else {
      
      skip_expand <- as.numeric(str_split(skip, "-", simplify = T))
      skip_beats <- seq(skip_expand[1], skip_expand[2])
      
    }
    
  } else {
    
    skip_beats <- as.numeric(str_split(skip, ", ", simplify = T))
    
  }
  
}

col_head <- c("time", seq(1, (nrow(wf_feet)-1)))
pw_raw_keep <- matrix(as.numeric(NA), nrow = Fs_new*2, ncol =  nrow(wf_feet), dimnames = list(NULL, col_head)) |> 
  as_tibble() |> 
  mutate(time = seq(1/Fs_new,2,1/Fs_new))

pw_raw_remove <- matrix(as.numeric(NA), nrow = Fs_new*2, ncol =  nrow(wf_feet), dimnames = list(NULL, col_head)) |> 
  as_tibble() |> 
  mutate(time = seq(1/Fs_new,2,1/Fs_new))

pw_raw_all <- matrix(as.numeric(NA), nrow = Fs_new*2, ncol =  nrow(wf_feet), dimnames = list(NULL, col_head)) |> 
  as_tibble() |> 
  mutate(time = seq(1/Fs_new,2,1/Fs_new))

df_cine200$bad <- 0
before_beat1 <- which(df_cine200$time_s < min(wf_feet$time_tang))
df_cine200$bad[before_beat1] <- 1
after_lastbeat <- which(df_cine200$time_s > max(wf_feet$time_tang))
df_cine200$bad[after_lastbeat] <- 1

for (b in seq(1:(nrow(wf_feet)-1))){
  
  beat <- df_cine200 |> dplyr::filter(time_s >= wf_feet$time_tang[b], time_s < wf_feet$time_tang[(b+1)])
  
  pw_raw_all[1:nrow(beat),(b+1)] <- beat$peak_velocity_cmps_filt20
  
  if (b %in% skip_beats) { # skips bad beat and moves to next
    
    continuous_timing <- which(df_cine200$time_s >= wf_feet$time_tang[b] & df_cine200$time_s < wf_feet$time_tang[(b+1)])
    df_cine200$bad[continuous_timing] <- 1
    
    pw_raw_remove[1:nrow(beat),(b+1)] <- beat$peak_velocity_cmps_filt20
    
    # add single beat variables to beat table
    beats_table[summary_index, 1:6] <- list(id0, 
                                            visit0,
                                            protocol0,
                                            as.character(b),
                                            beat$time_s[1],
                                            tail(beat$time_s,1) - beat$time_s[1])
    beats_table[summary_index, 7:18] <- NA
    
    summary_index <- summary_index+1
    
    print(b)
    next
  }
  
  pw_raw_keep[1:nrow(beat),(b+1)] <- beat$peak_velocity_cmps_filt20
  
  # add single beat variables to beat table
  beats_table[summary_index, 1:9] <- list(id0, 
                                          visit0,
                                          protocol0,
                                          as.character(b),
                                          beat$time_s[1],
                                          tail(beat$time_s,1) - beat$time_s[1],
                                          max(beat$peak_velocity_cmps_filt20, na.rm = T), # vs_cmps
                                          min(beat$peak_velocity_cmps_filt20, na.rm = T), # vd_cmps
                                          mean(beat$peak_velocity_cmps_filt20, na.rm = T)) # vm_cmps
  
  beats_table$peak_vamp_cmps[summary_index] <- beats_table$peak_vs_cmps[summary_index] - beats_table$peak_vd_cmps[summary_index] # Vamp
  beats_table$peak_pi[summary_index] <- beats_table$peak_vamp_cmps[summary_index] / beats_table$peak_vm_cmps[summary_index] # PI = Vamp/Vm
  beats_table$peak_ri[summary_index] <- beats_table$peak_vamp_cmps[summary_index] / beats_table$peak_vs_cmps[summary_index] # RI = Vamp/Vs
  
  beats_table[summary_index, 13:15] <- list(max(beat$mean_velocity_cmps_filt20, na.rm = T), # vs_cmps
                                            min(beat$mean_velocity_cmps_filt20, na.rm = T), # vd_cmps
                                            mean(beat$mean_velocity_cmps_filt20, na.rm = T)) # vm_cmps
  
  beats_table$mean_vamp_cmps[summary_index] <- beats_table$mean_vs_cmps[summary_index] - beats_table$mean_vd_cmps[summary_index] # Vamp
  beats_table$mean_pi[summary_index] <- beats_table$mean_vamp_cmps[summary_index] / beats_table$mean_vm_cmps[summary_index] # PI = Vamp/Vm
  beats_table$mean_ri[summary_index] <- beats_table$mean_vamp_cmps[summary_index] / beats_table$mean_vs_cmps[summary_index] # RI = Vamp/Vs
  
  summary_index <- summary_index+1
  
  print(b)
}

pw_raw_keep_long <- pivot_longer(pw_raw_keep, cols = 2:ncol(pw_raw_keep), values_to = "peak_velocity", names_to = "beat") |> 
  mutate(beat = as.numeric(beat),
         beat_pct = beat/max(beat)*100, 
         beat_c = as.factor(beat))

pw_raw_rem_long <- pivot_longer(pw_raw_remove, cols = 2:ncol(pw_raw_remove), values_to = "peak_velocity", names_to = "beat") |> 
  mutate(beat = as.numeric(beat),
         beat_pct = beat/max(beat)*100)

pw_raw_all_long <- pivot_longer(pw_raw_all, cols = 2:ncol(pw_raw_all), values_to = "peak_velocity", names_to = "beat") |> 
  mutate(beat = as.numeric(beat),
         beat_pct = beat/max(beat)*100)

pw_raw_keep_long_mean_sd <- pw_raw_keep_long |> 
  group_by(time) |> 
  summarize(across(.cols = "peak_velocity", 
                   .fns = list(mean = ~mean(., na.rm = T),
                               sd = ~sd(., na.rm = T)),
                   .names = "{.col}_{.fn}"), .groups = "drop")

x_limit <- round(mean(diff(wf_feet$time_tang,1)) + mean(diff(wf_feet$time_tang,1))*.25, 1)

pw_raw_keep_long |> 
  dplyr::filter(beat > 0, beat < 201) |> 
  ggplot() +
  geom_line(aes(x = time, y = peak_velocity, group = beat_c, colour = beat_c), #beat_pct), 
            alpha = 0.7, linewidth = 0.75, show.legend = T) +
  scale_colour_discrete(name = "Beat #") +
  scale_y_continuous(name = "Peak Velocity") +
  scale_x_continuous(name = "Time (s)", limits = c(0, x_limit), expand = c(0,0), breaks = seq(0, 3, .4), minor_breaks = NULL) +
  ggtitle(str_c(id0, ": individual beats", sep = "")) +
  theme_bw() +
  theme(text = element_text(size = 18, colour = "black"),
        axis.title = element_text(size = 20), 
        line = element_line(colour = "black"), 
        panel.spacing = unit(1, "cm"),
        legend.key.width = unit(2,"cm"))

ggplotly()

# split to individuals beats
# indy_beats <- matrix(as.numeric(), nrow = 500, ncol = nrow(wf_feet)-1, dimnames = list(NULL, seq(1:(nrow(wf_feet)-1)))) |> 
#   as_tibble(.name_repair = "check_unique")
# 
# for (i in 1:ncol(indy_beats)){
#   
#   if (any(!is.na(skip_beats)) & (i %in% skip_beats)) { # remove beat
#     #df_cine200$peak_velocity_mps[wf_feet$tangfoot[i]:(wf_feet$tangfoot[(i+1)]-1)] <- NA
#     df_cine200$peak_velocity_mps_filt20[wf_feet$tangfoot[i]:(wf_feet$tangfoot[(i+1)]-1)] <- NA
#     next
#   }
#   
#   if (!((wf_feet$time_tang[(i+1)]-wf_feet$time_tang[i]) > (1.5*rri_est))) { # ignore feet far away, likely two beats between
#     indy_beats[1:(wf_feet$tangfoot[(i+1)]-wf_feet$tangfoot[i]),i] <- df_cine200$peak_velocity_mps_filt20[wf_feet$tangfoot[i]:(wf_feet$tangfoot[(i+1)]-1)]
#   }
#   
# }

df_cine200$peak_velocity_cmps_filt20[1:(wf_feet$tangfoot[1] - 1)] <- NA # convert head to NA
df_cine200$peak_velocity_cmps_filt20[wf_feet$tangfoot[nrow(wf_feet)]:nrow(df_cine200)] <- NA # convert tail to NA

plotlist_trace[[pr]] <- ggplot() +
  geom_path(data = df_cine200, aes(x = time_s, y = peak_velocity_cmps), colour = "black") +
  geom_path(data = df_cine200, aes(x = time_s, y = peak_velocity_cmps_filt20), colour = "red") +
  geom_point(data = wf_feet, aes(x = time_tang, y = vel_tang), size = 2, shape = 21, fill = "dodgerblue") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(limits = c(0,100), expand = c(0,0)) +
  labs(title = id_visit_protocol) +
  theme_classic(); print(plotlist_trace[[pr]])

df_cine200 <- df_cine200[wf_feet$tangfoot[1]:(wf_feet$tangfoot[nrow(wf_feet)] - 1), ] # remove head and tail

# remove columns all na
indy_beats <- indy_beats |> select_if(not_all_na)

# add cine # to beat # as column name and add to data from other cines for same participant
indy_beats <- indy_beats |> rename_with(~paste(r, .x, sep = "_"))
combined_beats_peak <- bind_cols(combined_beats_peak, indy_beats)

pr <- pr+1

# if all data is NA, move to next participant
if (all(is.na(combined_beats_peak))) next

combined_beats_peak <- combined_beats_peak |> 
  select_if(not_all_na) |> # removed any columns that are all NA
  filter_all(any_vars(!is.na(.))) # remove any rows that are all NA

# set start time for each beat to zero
temp_seq <- seq(0, nrow(combined_beats_peak)/Fs_new, 1/Fs_new)
combined_beats_peak$time_s <- temp_seq[1:nrow(combined_beats_peak)]

## normalize beats to cardiac cycle ----

# create new tibble for time normalized waveforms, replace time with cycle fraction
combined_beats_peak_norm <- matrix(as.numeric(), nrow = 101, ncol = ncol(combined_beats_peak), 
                                   dimnames = list(NULL, colnames(combined_beats_peak))) |> 
  as_tibble(.name_repair = "check_unique") |> 
  select(-time_s) |> 
  mutate(cycle_fract = seq(0, 1, 0.01))

# process one beat-by-beat
for (i in 1:(ncol(combined_beats_peak)-1)){
  
  # select single beat
  single <- combined_beats_peak |> 
    select(all_of(i), time_s) |> 
    drop_na()
  
  # normalize time
  single$time_normalized <- (single$time_s - single$time_s[1])/(max(single$time_s)-single$time_s[1])
  
  # linear interpolation for cycle fraction from time
  combined_beats_peak_norm[ ,i] <- approx(x = single$time_normalized, 
                                          y = single[[1]], 
                                          xout = combined_beats_peak_norm$cycle_fract, 
                                          rule = 2, 
                                          na.rm = F)$y
  
  
  
}

## calculate mean waveform and plot ----
peak_velocity_mean_sd <- combined_beats_peak |>
  rowwise() |> 
  mutate(mean = mean(c_across(!starts_with("time")), na.rm = T),
         sd = sd(c_across(!starts_with("time")), na.rm = T))

keep <- tibble(id = id0,
               cine = rep0,
               time_s = peak_velocity_mean_sd$time_s,
               peak_velocity_mps = peak_velocity_mean_sd$mean)

df_combined_mean <- bind_rows(df_combined_mean, keep)

combined_beats_long <- pivot_longer(combined_beats_peak, 
                                    cols = 1:(ncol(combined_beats_peak)-1), 
                                    names_sep = "_", 
                                    names_to = c("cine", "beat"), 
                                    values_to = "peak_velocity") |> 
  arrange(cine, beat) |> 
  mutate(beat = as_factor(beat),
         cine = as_factor(cine))

indy_t <- combined_beats_long |> ggplot(aes(x = time_s, y = peak_velocity, group = interaction(beat, cine))) +
  geom_path(aes(colour = cine), linewidth = 1.5, show.legend = T)+
  scale_x_continuous(name = "Time (s)", limits = c(0, 1.6), breaks = seq(0, 1.6, by = 0.2), expand = c(0,0), minor_breaks = NULL) +
  scale_y_continuous(name = "Blood velocity (m/s)", limits = c(0, 1.4), breaks = seq(0, 1.4, by = .2), expand = c(0,0), minor_breaks = NULL) +
  labs(title = str_c(id0, vessel0, sep = " ")) +
  annotate("text", x = 1.3, y = 1.3, label = paste("# beats: ", (ncol(combined_beats_peak)-1)), size = 28/.pt) +
  theme_bw() +
  theme(line = element_line(colour = "black", lineend = "square", linewidth = 1),
        text = element_text(size = 28),
        axis.text = element_text(size = 24),
        legend.position = c(.8,.9), 
        legend.justification = c("right", "top"), 
        plot.margin = margin(1,1,1,1,"cm"))

mean_t <- peak_velocity_mean_sd |> ggplot(aes(x = time_s, y = mean)) +
  geom_ribbon(aes(ymin = mean-sd, ymax = mean+sd), alpha = 0.2) +
  geom_path(colour = "blue", linewidth = 2) +
  scale_x_continuous(name = "Time (s)", limits = c(0, 1.6), breaks = seq(0, 1.6, by = 0.2), expand = c(0,0), minor_breaks = NULL) +
  scale_y_continuous(name = "Blood velocity (m/s)", limits = c(0, 1.4), breaks = seq(0, 1.4, by = .2), expand = c(0,0), minor_breaks = NULL) +
  #labs(title = str_c(id0, vessel0, sep = " ")) +
  theme_bw() +
  theme(line = element_line(colour = "black", lineend = "square", linewidth = 1),
        text = element_text(size = 28),
        axis.text = element_text(size = 24),
        plot.margin = margin(1,1,1,1,"cm"))

plots_append <- plot_grid(indy_t, mean_t, align = "hv", axis = "tblr", nrow = 1)

fn <- here("output", "qc", "PW_Dopp", str_c("CAR_study_", vessel0, "v_indy_time_", id0, ".png", sep = ""))
ggsave(fn, plots_append, width = 700, height = 400, units = "mm")


fn <- here("output", str_c(id0, visit0, protocol0, "VA_cleaned_individual_bxb.csv", sep = "_"))
write_csv(beats_table, fn, na = "")

fn <-  here("output", str_c("LBNP_study_VA_frame_rates.csv", sep = ""))
write_csv(video, fn, append = T)

### plot full trace
velocity_plots_grid <- list()
p <- 1
for (pp in seq(1, length(plotlist_trace), by = 4)){
  if (pp+3 > length(plotlist_trace)){
    rem <- length(plotlist_trace) - pp
    velocity_plots_grid[[p]] <- do.call("plot_grid", c(plotlist_trace[pp:(pp+rem)], nrow = 4, align = "hv", axis = "tblr"))
  } else {
    velocity_plots_grid[[p]] <- do.call("plot_grid", c(plotlist_trace[pp:(pp+3)], nrow = 4, align = "hv", axis = "tblr"))
  }
  
  temp_plot <- add_sub(velocity_plots_grid[[p]], "Black: original data\nRed: retained waveforms (20Hz Butterworth low pass filter)\nBlue dot: waveform foot detection (tangent method)", x = 0.01, hjust = 0)
  
  fn <- here("output", "PW_Dopp", str_c("LBNP_study_", vessel0, "v_envelopes", p, ".png", sep = ""))
  ggsave(fn, temp_plot, width = 350, height = 200, units = "mm")
  p <- p+1
}

message("Complete")
