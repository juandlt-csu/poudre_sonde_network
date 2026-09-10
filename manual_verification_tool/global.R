#global variables

library(shiny)
library(bslib)
library(DT)
library(tidyverse)
library(lubridate)
library(here)
library(ggpubr)
library(gridExtra)
library(plotly)
library(keys)
library(patchwork)
library(digest)
library(fs)
library(shinyFiles)
library(shinyWidgets)
library(glue)
library(anytime)
library(arrow)
library(yaml)
library(cdssr) # Devtools packge, code to install:
# install.packages("devtools")
#devtools::install_github("anguswg-ucsb/cdssr")

options(shiny.maxRequestSize = 10000 * 1024^2)
options(shiny.autoload.r = TRUE)

`%nin%` = Negate(`%in%`)

#Source global root config first
source(here("global.R"))

#Source Helper files
walk(list.files(here("manual_verification_tool", "R"), pattern = "\\.R$", full.names = TRUE), source)
##### Colors + parameters #####

site_color_combo <- tibble(site = c("joei", "cbri", "chd", "pfal", "sfm", "pbr", "pman", "pbd","bellvue","salyer", "udall", "riverbend_virridy", "riverbend",
                                    "cottonwood_virridy", "cottonwood","elc",  "archery_virridy", "archery", "boxcreek", "springcreek", "riverbluffs"),
                           color = c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC",
                                     "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44","#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788"))


final_status_colors <- c("PASS" = "green",
                         "OMIT" = "red",
                         "FLAGGED" = "orange")



available_parameters <- get_filenames()%>%mutate(
  parameter = map_chr(filename, ~ split_filename(.x)$parameter))%>%
  pull(parameter)%>%
  unique()


available_sites <- get_filenames()%>%mutate(
  site = map_chr(filename, ~ split_filename(.x)$site))%>%
  pull(site)%>%
  unique()

#TODO: Automate for public version or ask for user input
available_flags <- read_csv(here("data", "raw", "sensor", "manual_data_verification", year_cycle, "in_progress", "meta", "available_flags.csv"), show_col_types = F)%>%
  pull(flags)%>%
  unique()

#### USGS STREAMFLOW API QUERY
cdwr_creds <- read_yaml(here("creds","CDWRCreds.yml"))
stations <- c("CLAFTCCO", "CLAFORCO", "CLABOXCO", "CLASRKCO", "JWCCHACO")

global_usgs_flow_data <- tryCatch({

  purrr::map_df(stations, function(station) {

    # Fetch telemetry time series for the current station
    data <- tryCatch({
      get_telemetry_ts(
        abbrev = station,
        start_date = usgs_flow_dt_start,
        end_date = usgs_flow_dt_end,
        timescale = "raw",
        api_key = cdwr_creds$api_key
      )
    }, error = function(e) {
      warning("Failed to fetch data for station: ", station, " - ", e$message)
      return(NULL)
    })

    # Add abbreviation column to keep track of locations
    if (!is.null(data) && nrow(data) > 0) {
      data$station_abbrev <- station
      return(data)
    }

    return(NULL)
  })

}, error = function(e) {
  warning("Failed to fetch global USGS flow data: ", e$message)
  return(NULL) # Explicitly returns NULL to global_usgs_flow_data
})

###### End Helper Functions ######
