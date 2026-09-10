#Installing and loading all packages
invisible(
  lapply(c(
    "tidyverse", # Data manipulation
    "lubridate", # Date-Time Manipulation
    "readxl", # Reading excel files
    "here", # Easy, cross platform file referencing
    "ggplot2", # Plotting libraries
    "patchwork",
    "ggpubr",
    "bslib",
    "shinyWidgets",
    "shinyFiles",
    "fs",
    "plotly",
    "devtools", # For downloading GitHub packages
    "remotes",
    "arrow", # Reading parquet files
    "groupdata2",
    "shiny",
    "DT",
    "ross.wq.tools", # In house functions
    "writexl"
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

# EDIT: Year and cycle values set in the Process_2022.rmd
# This prevents the app from needing to be updated every year when new data is added
####---- SET YEAR -----###
# year = "2022"
# year_cycle = paste0( year, "_cycle")
#### ----------------- ###

in_progress_dir = in_progress_path

#Check to see if post verification folder exists
post_ver_dir = file.path(year_cycle_path, "post_verification")
if(!dir.exists(post_ver_dir)){
  dir.create(post_ver_dir)
}
excel_path <- file.path(post_ver_dir, "drift_corrections.xlsx")

# Template to find downstream and upstream data
site_order <- list(
  clp = c("joei", "cbri", "chd", "pfal", "pbr", "pman", "pbd", "bellvue",
          "salyer", "udall", "riverbend", "cottonwood", "elc", "archery", "riverbluffs"),
  springcreek = c("riverbend", "springcreek", "cottonwood"),
  boxcreek = c("elc", "boxcreek", "archery"),
  sfm = c("sfm"),
  mtn_campus = c("mtn_campus")
)

# Load source scripts containing the functions from your RMD
source("R/drift_functions.R")
field_notes <- grab_mWater_sensor_notes(load_mWater(creds = yaml::read_yaml(here::here("creds", "mWaterCreds.yml"))))

# Load data
verified_files <- list.files(file.path(in_progress_dir, "verified_directory"),
                             pattern = "Turbidity|FDOM", ignore.case = TRUE, full.names = TRUE)
ver_site_params <- str_remove(basename(verified_files), "_FINAL_.*\\.parquet$")
#Post Verified files
post_ver_files <- list.files(post_ver_dir,
                             pattern = "Turbidity|FDOM", ignore.case = TRUE, full.names = TRUE)
post_ver_site_params <- str_remove(basename(post_ver_files), ".parquet$")

prepped_data <- map(verified_files, arrow::read_parquet) %>%
  set_names(str_remove(basename(verified_files), "_FINAL_.*\\.parquet$")) %>%
  map(function(df) {
    df %>%
      add_field_notes(df = ., notes = field_notes) %>%
      mutate(
        drift = ifelse(stringr::str_detect(user_flag, "drift"), TRUE, FALSE),
        drift = ifelse(is.na(drift), FALSE, drift),
        mean_analysis = ifelse(verification_status %in% c("FLAGGED", "PASS"), mean, NA),
        DT_round = with_tz(DT_round, "MST")
      ) %>%
      fill(c("turb_pre_clean", "turb_post_clean", "fdom_pre_clean", "fdom_post_clean"),
           .direction = "down")
  })

# Get site options
# Extract valid site-parameter combinations from your actual loaded keys
# If keys look like "bellvue-Turbidity", "sfm-FDOM", etc.
all_combos <- names(prepped_data)
#Filter out files that have already been post-verified
remaining_sites <- setdiff(all_combos, post_ver_site_params)

# Auto-pass no-drift combos ----
# If a site-param has no drift flags whatsoever there is nothing for the user
# to correct. Write it straight to post_verification now so it never appears
# in the app queue. post_ver_site_params and remaining_sites are refreshed
# afterwards so the server sees a fully accurate picture at startup.
no_drift_combos <- remaining_sites %>%
  keep(~ {
    df <- prepped_data[[.x]] %>% filter(!is.na(mean_analysis))
    nrow(df) > 0 && !any(df$drift, na.rm = TRUE)
  })

if (length(no_drift_combos) > 0) {
  message("Auto-passing ", length(no_drift_combos),
          " combo(s) with no drift flags to post_verification...")

  walk(no_drift_combos, function(combo) {
    prepped_data[[combo]] %>%
      mutate(mean_drift_trans = mean_analysis, correction_type = "raw")%>%
      arrow::write_parquet(file.path(post_ver_dir, paste0(combo, ".parquet")))

    message("  Auto-passed: ", combo)
  })

  # Refresh so everything downstream in this file — and the server's
  # completed_combos seed — reflects the combos just written.
  post_ver_files       <- list.files(post_ver_dir,
                                     pattern = "Turbidity|FDOM",
                                     ignore.case = TRUE, full.names = TRUE)
  post_ver_site_params <- str_remove(basename(post_ver_files), ".parquet$")
  remaining_sites      <- setdiff(all_combos, post_ver_site_params)
}

# Dynamically extract unique sites and parameters directly from your dataset keys
site_choices <- map_chr(str_split(all_combos, "-"), 1) %>% unique()

param_choices <- c("Turbidity", "FDOM Fluorescence")

# Sites that still have at least one unprocessed parameter — used to seed the UI
# on startup so already-finished combos are excluded from the beginning
remaining_site_names <- map_chr(str_split(remaining_sites, "-"), 1) %>% unique()
if(length(remaining_site_names) == 0 ){
  message("All site-parameter combinations have been drift corrected. Move to the publishing step or next year's data.")
  stop()
}
