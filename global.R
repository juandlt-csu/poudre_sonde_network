##### GLOBAL OPTIONS #####
# This file sets global options for all necessary tools within the PSN Repository #

### Variables ###
# Set the year & data cycle
year <- 2022
year_cycle <- paste0(year, "_cycle")
usgs_flow_dt_start <- paste0(year, "-01-01")
usgs_flow_dt_end <- paste0(year, "-12-31")

# Set year_cycle path (the directory that contains all the data for a specific year)
year_cycle_path <- here::here("data", "raw", "sensor", "manual_data_verification", year_cycle)
# In progress path
in_progress_path <- here::here(year_cycle_path, "in_progress")
# Raw data path
raw_data_path <- here::here(in_progress_path, "raw_data")
# All data path
all_data_path <- here::here(in_progress_path, "all_data_directory")
# Pre-verification path
pre_verification_path <- here::here(in_progress_path, "pre_verification_directory")
# Intermediary path
intermediary_path <- here::here(in_progress_path, "intermediary_directory")
# Verified data path
verified_path <- here::here(in_progress_path, "verified_directory")

# Meta path
meta_path <- here::here(in_progress_path, "meta")
# Log dir path
log_dir_path <- here::here("data","raw", "sensor", "log_download", year)
