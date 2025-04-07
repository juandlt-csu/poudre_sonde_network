library(tidyverse)
library(here)
# Temperature and turbidity

# Salyer natural area sonde
# river bend ponds/timberline

# Pull in the data
data <- readRDS("data/virridy_verification/all_data_flagged_complete.RDS")

# filter the data
data_temp <- data %>%
  keep(str_detect(names(.), "Temperature")) %>%
  keep(!str_detect(names(.), "virridy")) %>%
  keep(str_detect(names(.), "timberline|legacy")) %>%
  names()

data_turb <- data %>%
  keep(str_detect(names(.), "Turbidity")) %>%
  keep(!str_detect(names(.), "virridy")) %>%
  keep(str_detect(names(.), "timberline|legacy"))

# average the data by week
data_temp_week <- data_temp %>%
  map(~.x %>%
      mutate(
        week = week(DT_round),
        y_w = paste(year, "-", week)
      ) %>%
      group_by(y_w) %>%
      mutate(weekly_mean = round(mean(mean, na.rm = T), 3)) %>%
      select(site, parameter, year_week = y_w, weekly_mean) %>%
      distinct() %>%
      mutate(site = ifelse(site == "Timberline", "River Bend Ponds", "Salyer Natural Area"))
    )

data_turb_week <- data_turb %>%
  map(~.x %>%
        mutate(
          week = week(DT_round),
          y_w = paste(year, "-", week)
        ) %>%
        group_by(y_w) %>%
        mutate(weekly_mean = round(mean(mean, na.rm = T), 3)) %>%
        select(site, parameter, year_week = y_w, weekly_mean) %>%
        distinct() %>%
        mutate(site = ifelse(site == "Timberline", "River Bend Ponds", "Salyer Natural Area"))
  )

# Save the data

data_temp_week_output <- data_temp_week %>%
  set_names(c("salyer_natural_area_temperature", "river_bend_ponds_temperature"))

iwalk(data_temp_week_output, ~write_csv(.x, here(paste0(.y, ".csv"))))

data_turb_week_output <- data_turb_week %>%
  set_names(c("salyer_natural_area_turbidity", "river_bend_ponds_turbidity"))

iwalk(data_turb_week_output, ~write_csv(.x, here(paste0(.y, ".csv"))))
