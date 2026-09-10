# global.R
library(shiny)
library(ross.wq.tools)
library(DT)

#Installing and loading all packages
invisible(
  lapply(c(
    "tidyverse", # Data manipulation
    "janitor", # Clean dirty data
    "lubridate", # Date-Time Manipulation
    "rvest", # HTML Retrieval and Manipulation
    "readxl", # Reading excel files
    "here", # Easy, cross platform file referencing
    "ggplot2", # Plotting libraries
    "ggpubr",
    "plotly",
    "devtools", # For downloading GitHub packages
    "remotes",
    "yaml",
    "arrow", # reading parquet files
    "groupdata2"
  ),
  function(x) {
    if (x %in% installed.packages()) {
      suppressMessages({
        library(x, character.only = TRUE)
      })
    } else {
      suppressMessages({
        install.packages(x)
        library(x, character.only = TRUE)
      })
    }
  })
)

# Define file paths for tracking system
# EDIT: Update the data_dir path to point to the data folder in the parent poudre_sonde_network dir.
data_dir <- file.path(year_cycle_path, "hydro_vu_pull","back_calibration")
tracking_file <- file.path(data_dir,"calibrated_sensor_field_data_tracking.rds")
finalized_file <- file.path(data_dir, "calibrated_sensor_field_data_finalized.rds")
original_file <- file.path(data_dir,paste0("calibrated_sensor_data_", year, ".rds"))

# Initialize tracking data (unverified calibrations)
if (file.exists(tracking_file)) {
  calibrated_data_tracking <- readr::read_rds(tracking_file)
  #If this year's data is empty, stop the app and display a message
  if(length(calibrated_data_tracking[[year]]) == 0){
    stop("Tracking data is empty. Either there is an error with the file or corrections are completed")
  }
} else {
  # First time running - copy original to tracking
  calibrated_data_tracking <- readr::read_rds(original_file)
  readr::write_rds(calibrated_data_tracking, tracking_file)
}

# Initialize finalized data (verified calibrations)
if (file.exists(finalized_file)) {
  final_calibrated_data <- readr::read_rds(finalized_file)
} else {
  # First time running - create empty nested list structure matching tracking data
  final_calibrated_data <- lapply(calibrated_data_tracking, function(year_data) {
    list()  # Empty list for each year
  })
  readr::write_rds(final_calibrated_data, finalized_file)
}

# Read in the sensor specific calibration data
sensor_calibration_data <- readr::read_rds(here::here("data", "collated", "sensor", "cal_reports",
                                                      "sensor_calibration_data.RDS"))

ross.wq.tools::load_calibration_data(
  cal_data_file_path = here::here("data", "collated", "sensor", "cal_reports",
                                  "munged_calibration_data.RDS")
)

# Read in field notes

field_notes_data <- ross.wq.tools::load_mWater(creds = yaml::read_yaml(here::here("creds", "mWaterCreds.yml")))%>%
  mutate(DT_round = floor_date(start_dt, unit = "15 minutes"))%>%
  filter(str_detect(visit_type, "Sensor Calibration, Cleaning or Check"))%>%
  pivot_longer(cols = contains(c("post", "pre" )),
               names_to = "col",
               values_to = "value")%>%
  select(site, DT_round, col, value)%>%
  mutate(parameter = str_split(col, "_", simplify = TRUE)[,1],
         parameter_clean = case_when(
           parameter == "rdo" ~ "DO",
           parameter == "cond" ~ "Specific Conductivity",
           parameter == "turb" ~ "Turbidity",
           parameter == "fdom" ~ "FDOM Fluorescence",
           parameter == "chla" ~ "Chl-a Fluorescence",
           parameter == "orp" ~ "ORP",
           parameter == "ph" ~ "pH",
           TRUE ~ parameter
         ),
         type = str_split(col, paste0(parameter, "_"), simplify = T)[,2])%>%
  ross.wq.tools::fix_site_names()
#add old field notes if year is 2023 or earlier
if(year <= 2023){
  old_field_notes <- readxl::read_excel(here("data","raw","field_notes","sensor_field_notes.xlsx")) %>%
    mutate(DT = (paste0(date, " ", start_time_mst))) %>%
    mutate(DT = ymd_hm(DT) + hours(7)) %>%
    arrange(DT) %>%
    mutate(DT_round = round_date(DT, "15 minutes")) %>%
    mutate(DT_round = with_tz(DT_round, tzone = "UTC"),
           last_site_visit = with_tz(DT_round, tzone = "UTC"),
           DT_join = as.character(DT_round),
           sonde_moved = NA,
           sonde_employed = case_when(is.na(sensor_deployed) & is.na(sensor_pulled) ~ NA,
                                      sensor_deployed == "x" ~ 0,
                                      sensor_pulled == "x" ~ 1),
           site = ifelse(site == "rist", "bellvue", site)) %>%
    ross.wq.tools::fix_site_names()%>%
    #add extra cols
    mutate(rdo_pre_clean = NA,
           rdo_post_clean = NA,
           rdo_post_cal = NA,
           cond_pre_clean = NA,
           cond_post_clean = NA,
           cond_post_cal = NA,
           ph_pre_clean = NA,
           ph_post_clean = NA,
           ph_post_cal = NA,
           orp_pre_clean = NA,
           orp_post_clean = NA,
           orp_post_cal = NA,
           turb_pre_clean = NA,
           turb_post_clean = NA,
           turb_post_cal = NA,
           fdom_pre_clean = NA,
           fdom_post_clean = NA,
           fdom_post_cal = NA,
           chla_pre_clean = NA,
           chla_post_clean = NA,
           chla_post_cal = NA)%>%
    select(-sensors_cleaned, -wiper_working, -cal_report_collected)

  all_field_notes <- bind_rows(field_notes_data, old_field_notes)
}
# Source helper functions
# source(here::here("R", "generate_final_df.R"))
# source(here::here("R", "cal_plot.R"))
# source(here::here("R", "final_calibration_plot.R"))
# source(here::here("R", "update_backend.R"))
