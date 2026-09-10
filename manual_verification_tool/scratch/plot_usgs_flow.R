
#pull in data from 9/1/19 - 12/31-19 for 06752260
#usgs lincoln: CLAFORCO
#usgs boxelder: CLABOXCO
#devtools::install_github("anguswg-ucsb/cdssr")
library(cdssr)
library(tidyverse)
library(plotly)

start_dt = "2022-10-29"
end_dt = "2022-11-12"

canyon <- get_telemetry_ts(abbrev = "CLAFTCCO", start_date = start_dt, end_date = end_dt, timescale = "raw") # for PBD, Bellvue
lincoln <- get_telemetry_ts(abbrev = "CLAFORCO", start_date = start_dt, end_date = end_dt, timescale = "raw") # for Salyer, udall, riverbend
boxelder <- get_telemetry_ts(abbrev = "CLABOXCO", start_date = start_dt, end_date = end_dt, timescale = "raw") # for Cottonwood, elc, archery

ggplotly(
  ggplot() +
    geom_line(data = canyon, aes(x = datetime, y = meas_value)) +
    labs(x = "Date", y = "Streamflow (cfs)")
)

ggplotly(
  ggplot() +
    geom_line(data = lincoln, aes(x = datetime, y = meas_value)) +
    labs(x = "Date", y = "Streamflow (cfs)")
)

ggplotly(
  ggplot() +
    geom_line(data = boxelder, aes(x = datetime, y = meas_value)) +
    labs(x = "Date", y = "Streamflow (cfs)")
)
