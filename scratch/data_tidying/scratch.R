---
  title: "Organizing raw data"
author: "ROSSyndicate"
date: "`r Sys.Date()`"
output: html_document
editor_options:
  markdown:
  wrap: 90
---

  ```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, error = FALSE, message = FALSE)
```

Load necessary packages:

  ```{r}
source("src/package_loader.R")
lapply(c("data.table", "tidyverse", "rvest", "readxl", "lubridate", "zoo", "padr","plotly", "feather", "RcppRoll", "yaml"), package_loader)
```

# Import and collate data

Load and munge field calibration files:

  ```{r}
cal_files <- list.files("data/calibration_reports", pattern=".html")

cal_tabler <- function(cal_files){

  cal <- read_html(file.path("data/calibration_reports/", cal_files)) %>%
    html_nodes("div") %>%
    html_text() %>%
    as_tibble()

  rdo <- cal %>% filter(grepl("RDO", value)) %>% pull() %>% str_replace_all(., " ", "") %>% tolower()

  ph_orp <- cal %>% filter(grepl("pH/ORP", value)) %>% pull() %>% str_replace_all(., " ", "") %>% tolower()

  conductivity <- cal %>% filter(grepl("Conductivity",value)) %>% pull() %>% str_replace_all(., " ", "") %>% tolower()

  turbidity <- cal %>% filter(grepl("Turbidity",value)) %>% pull() %>% str_replace_all(., " ", "") %>% tolower()

  # Always the fifth sensor when depth is available:
  try(depth <- cal %>% .[5,] %>% pull() %>% str_replace_all(., " ", "") %>% tolower())

  time_mst <- paste0(str_sub(cal_files, -13, -12),":", str_sub(cal_files, -11, -10))

  date <- paste0(str_sub(cal_files, -22, -19),"-", str_sub(cal_files, -18, -17),"-", str_sub(cal_files, -16, -15))

  cal_table <- tibble(site = sub("\\_.*", "", cal_files),

                      DT = ymd_hm(paste(date, time_mst, tz = "MST")),

                      # Depth
                      depth_cal_date = "None",
                      depth_offset = "None",
                      depth_ref_depth = "None",
                      depth_ref_offset = "None",
                      depth_pre_psi = "None",
                      depth_post_psi = "None",

                      # Dissolved Oxygen
                      rdo_slope = str_match(rdo, "slope\\s*(.*?)\\s*offset")[,2],
                      rdo_offset = str_match(rdo, "offset\\s*(.*?)\\s*mg/l")[,2],
                      rdo_100 = str_match(rdo, "premeasurement\\s*(.*?)\\s*%satpost")[,2],
                      rdo_conc = str_match(rdo, "concentration\\s*(.*?)\\s*mg/lpremeasurement")[,2],
                      rdo_temp = str_match(rdo, "temperature\\s*(.*?)\\s*°c")[,2],
                      rdo_pressure = str_match(rdo, "pressure\\s*(.*?)\\s*mbar")[,2],

                      # pH
                      ph_slope_pre = str_match(ph_orp, "offset1slope\\s*(.*?)\\s*mv/ph")[,2],
                      ph_offset_pre = str_match(ph_orp, "mv/phoffset\\s*(.*?)\\s*mvslopeandoffset2")[,2],
                      ph_slope_post = str_match(ph_orp, "offset2slope\\s*(.*?)\\s*mv/ph")[,2],
                      ph_offset_post = str_match(ph_orp, paste0(ph_slope_post,"mv/phoffset\\s*(.*?)\\s*mvorporp"))[,2],
                      # Sometimes, the post value can actually be in the high 6 pH... therefore the post measurement regex matching text is conditional
                      ph_7_nice = str_sub(str_match(ph_orp, "postmeasurementph7\\s*(.*?)\\s*mvcal")[,2], 10, nchar(str_match(ph_orp, "postmeasurementph7\\s*(.*?)\\s*mvcal")[,2])),
                      ph_7_other = str_sub(str_match(ph_orp, "postmeasurementph6\\s*(.*?)\\s*mvcal")[,2], 10, nchar(str_match(ph_orp, "postmeasurementph6\\s*(.*?)\\s*mvcal")[,2])),
                      ph_7 = ifelse(is.na(ph_7_nice), ph_7_other, ph_7_nice),

                      # ORP
                      #Newly encountered thing: sometimes the calibration report calls the ORP standard Zobell's, sometimes it's just called "ORP Standard":
                      orp_offset = ifelse(is.na(str_match(ph_orp, "zobell'soffset\\s*(.*?)\\s*mvtemperature")[,2]),
                                          str_match(ph_orp, "orpstandardoffset\\s*(.*?)\\s*mvtemperature")[,2],
                                          str_match(ph_orp, "zobell'soffset\\s*(.*?)\\s*mvtemperature")[,2]),

                      # Conductivity
                      tds_conversion_ppm = str_sub(str_match(conductivity, "tdsconversionfactor\\s*(.*?)\\s*cellconstant")[,2], 6, nchar(str_match(conductivity, "tdsconversionfactor\\s*(.*?)\\s*cellconstant")[,2])),
                      cond_cell_constant = str_match(conductivity, "cellconstant\\s*(.*?)\\s*referencetemperature")[,2],
                      cond_pre = str_match(conductivity,paste0(str_match(conductivity,
                                                                         "premeasurementactual\\s*(.*?)\\s*specificconductivity")[,2],"specificconductivity\\s*(.*?)\\s*µs/cmpost"))[,2],
                      cond_post = str_match(conductivity,paste0(str_match(conductivity,
                                                                          "postmeasurementactual\\s*(.*?)\\s*specificconductivity")[,2],"specificconductivity\\s*(.*?)\\s*µs/cm"))[,2],

                      # Turbidity
                      ntu_slope = "None",
                      ntu_offset = "None",
                      ntu_10 = "None",
                      ntu_100 = "None") %>%

    select(-c(ph_7_nice, ph_7_other))

  # Not all sondes have depth. So, we "try" to get the values.
  try(cal_table <- cal_table %>%
        mutate(
          # Depth
          depth_cal_date = str_match(depth, "lastcalibrated\\s*(.*?)\\s*calibrationdetails")[,2],
          depth_offset = str_match(depth, "zerooffset\\s*(.*?)\\s*psireferencedepth")[,2],
          depth_ref_depth = str_match(depth, "psireferencedepth\\s*(.*?)\\s*ftreferenceoffset")[,2],
          depth_ref_offset = str_match(depth, "ftreferenceoffset\\s*(.*?)\\s*psipremeasurement")[,2],
          depth_pre_psi = str_match(depth, "psipremeasurement\\s*(.*?)\\s*psipostmeasurement")[,2],
          depth_post_psi = str_match(depth, "psipostmeasurement\\s*(.*?)\\s*psi")[,2]))



  # Not all sondes have turbidity. So, we "try" to get the values.
  try(cal_table <- cal_table %>%
        mutate(
          # Turbidity
          ntu_slope = str_match(turbidity, "slope\\s*(.*?)\\s*offset")[,2],
          ntu_offset = str_match(turbidity, "offset\\s*(.*?)\\s*ntu")[,2],
          ntu_10 = str_match(turbidity, "calibrationpoint1premeasurement\\s*(.*?)\\s*ntupost")[,2],
          ntu_100 = str_match(turbidity, "calibrationpoint2premeasurement\\s*(.*?)\\s*ntupost")[,2]))

  cal_table <- cal_table %>%
    mutate(
      #Factory Defaults
      factory_defaults = paste0(ifelse(is.na(ntu_slope), "Turbidity ", ""),
                                ifelse(is.na(rdo_slope), "RDO ", ""),
                                ifelse(is.na(ph_slope_post), "pH ", ""),
                                ifelse(is.na(orp_offset), "ORP ", ""),
                                ifelse(is.na(cond_post), "Conductivity ", ""),
                                ifelse(is.na(depth_cal_date), "Depth ", "")))
  cal_table

}

cal_table <- map_dfr(cal_files, cal_tabler) %>%
  distinct(.keep_all = TRUE) %>%
  group_by(site) %>%
  mutate(DT = round_date(DT, "15 minutes")) %>%
  # filter for years that are 2022 and greater
  filter(year(DT) >= 2022)
# this creates a lot of errors, is it okay to ignore them all? Just make a note that these will create a lot of errors.

rm(cal_files, cal_tabler, package_loader)
```

Load field notes and define the start time as the 15 minute preceding the field time

```{r}
# Pulling in field notes and adding relevant datetime columns
field_notes <- read_excel("data/sensor_field_notes.xlsx") %>%
  mutate(start_DT = ymd_hm(paste(date, start_time_mst, tzone = "MST"))) %>%
  mutate(#start_DT = with_tz(start_DT, tzone = "MST"),
    DT_round = floor_date(start_DT, "15 minutes"),
    DT_join = as.character(DT_round),
    site = tolower(site),
    season = year(DT_round))

# Determining when the sonde was employed (SE) based on `sensor_pulled` (SP) and `sensor_deployed` (SD) columns.
## when (SP & SD) columns are NOT empty, SE = 0 (employed).
## when (SP) column is NOT empty, but (SD) column is empty, SE = 1 (NOT employed).
## when (SP) column is empty, but (SD) column is NOT empty, SE = 0 (employed).
## when (SP & SD) columns are empty, SE = NA (unknown).

## Downstream, fill() (fill in missing values with previous value)
## is used on the sonde_employed column after it has been joined
## to the pulled API data to determine if the sonde was employed for that data.
deployment_record <- field_notes %>%
  # filter for years that are 2022 and greater
  filter(year(DT_round) >= 2022) %>%
  arrange(site, DT_round) %>%
  group_by(site) %>%
  # `sonde_employed` determines if the sonde is deployed or not. 0 = sonde deployed, 1 = sonde is not deployed
  mutate(sonde_employed = case_when(!is.na(sensor_pulled) & !is.na(sensor_deployed) ~ 0,
                                    !is.na(sensor_pulled) & is.na(sensor_deployed) ~ 1,
                                    is.na(sensor_pulled) & !is.na(sensor_deployed) ~ 0,
                                    is.na(sensor_pulled) & is.na(sensor_deployed) ~ NA),
         last_site_visit = DT_round)
```

Munge field notes

```{r}
# field_notes_clean <- field_notes %>%
#   mutate(sonde_impact = case_when(grepl('sonde', visit_comments, ignore.case = T) ~ 'y',
#                                   grepl('calib', visit_comments, ignore.case = T) ~ 'y',
#                                   sensor_pulled == 'x' ~ 'y',
#                                   sensor_deployed == 'x' ~ 'y',
#                                   grepl('yes', sensors_cleaned, ignore.case = T) ~ 'y'
#   ),
#   end_DT = case_when(sensor_pulled == 'x' ~ ymd_hm(paste0(season, '-12-31 23:59'), tz = 'MST'),
#                      sensor_deployed == 'x' ~ start_DT,
#                      TRUE ~ start_DT + minutes(15)),
#   start_DT = case_when(sensor_deployed == 'x' ~ ymd_hm(paste0(season, '-01-01 00:00'), tz = 'MST'),
#                        TRUE ~ start_DT)) %>%
#   filter(sonde_impact == 'y', !is.na(start_DT))
```

Merge the data sets from all API pulls:

  ```{r}
all_data <- list.files(path = "data/api/", full.names = TRUE, pattern = "*.csv") %>%
  map_dfr(~data.table::fread(.) %>% select(-id)) %>%
  # remove overlapping API-pull data
  distinct() %>%
  # remove VuLink data
  filter(!grepl("vulink", name, ignore.case = TRUE)) %>%
  # Convert DT to MST:
  mutate(DT = as_datetime(timestamp, tz = "UTC")) %>%
  mutate(DT = with_tz(DT, tzone = "MST"),
         DT_round = round_date(DT, "15 minutes"),
         DT_join = as.character(DT_round),
         site = tolower(site)) %>%
  # filter for years that are 2022 and greater
  filter(year(DT_round) >= 2022) %>%
  # Lastly, we swapped Boxelder's sonde out for Rist's late in 2022:
  mutate(site = ifelse(site == "rist" & DT > "2022-09-20" & DT < "2023-01-01", "boxelder", site))
```

### Export collated raw file

```{r}
# This will be a parquet file in the future?
#write_feather(all_data, paste0('data/SOME_FOLDER_FOR_POSTERITY/collated_raw_sonde_v', Sys.Date(), '.feather'))
```

# Level 1 QA-QC

### Recoding data where sonde was out of water

Filter instances in which a sonde was pulled out of field mid-season (basically a backup filter in case a sonde was pulled out of water but the log didn't get stopped)

```{r}
#this function is not working... not going to debug right now, but this is meant to do the lifting of this code block, but in an automated way...
#
# recode_for_maintenance = function(start, end, site) {
#   # start_DT = ymd_hm(start, tz = 'MST')
#   # end_DT = ymd_hm(end, tz = 'MST')
#   all_data %>%
#     setDT(.) %>%
#     mutate(value = if_else(ymd_hms(DT) >= ymd_hms(start) &
#                              ymd_hms(DT) <= ymd_hms(end) &
#                              site == site,
#                            NA_real_,
#                            value),
#            flag = if_else(ymd_hms(DT) >= ymd_hms(start) &
#                             ymd_hms(DT) <= ymd_hms(end) &
#                             site == site,
#                           'recoded for sensor maintenance',
#                           NA_character_))
# }

# preserving this for now.
# all_data <- all_data %>%
#   setDT(.) %>%
#   # Rist
#   filter(!(ymd_hms(DT_round) >= ymd_hms("2021-09-11 09:00:00") & ymd_hms(DT_round) <= ymd_hms("2022-04-22 14:00:00") & site == "rist"),
#          !(ymd_hms(DT_round) <= ymd_hms("2022-05-30 09:00:00") & site == "rist"),
#          !(ymd_hms(DT_round) >= ymd_hms('2022-05-06 12:00:00') & ymd_hms(DT_round) <= ymd_hms('2022-05-09 14:30:00') & site == "rist"),
#          # Legacy
#          !(ymd_hms(DT_round) >= ymd_hms('2021-12-04 19:30:00') & ymd_hms(DT_round) < ymd_hms('2022-04-06 17:30:00') & site == "legacy"),
#          !(ymd_hms(DT_round) >= ymd_hms('2022-05-24 09:30:00') & ymd_hms(DT_round) < ymd_hms('2022-06-01 13:30:00') & site == "legacy"),
#          !(ymd_hms(DT_round) > ymd_hms('2022-07-08 14:00:00') & ymd_hms(DT_round) <= ymd_hms('2022-07-12 10:00:00') & site == "legacy"),
#          !(ymd_hms(DT_round) >= ymd_hms('2022-08-04 09:50:00') & ymd_hms(DT_round) <= ymd_hms('2022-08-25 16:15:00') & site == "legacy"),
#          !(ymd_hms(DT_round) > ymd_hms('2022-09-07 06:57:00') & ymd_hms(DT_round) <= ymd_hms('2022-09-18 07:00:00') & site == "legacy"),
#          # Timberline
#          !(ymd_hms(DT_round) > ymd_hms('2022-01-01 08:15:00') & ymd_hms(DT_round) <= ymd_hms('2022-04-06 08:15:00') & site == "timberline"),
#          # Archery
#          !(ymd_hms(DT_round) > ymd_hms('2022-10-04 15:00:00') & ymd_hms(DT_round) <= ymd_hms('2022-10-07 16:00:00') & site == "archery"))
```

### Temperature at Archery

I think we will want to develop pipelines specific for each site and parameter for
processing speed (targets?). Here I am starting a pipeline for Archery's temperature:

  #### Format data

  ```{r}
# Determine each site and parameter in all_data
sites <- unique(all_data$site)
params <- unique(all_data$parameter)
# Constructing a df to iterate over each site-parameter combination
combinations <- crossing(sites, params)

# Function to summarize site-parameter combinations
summarize_site_param <- function(site_arg, parameter_arg, api_data) {

  # removing print statement
  # print(paste0('Trying site ', site_arg, " with parameter ", parameter_arg, "."))

  # filtering the data and generating results
  summary <- tryCatch({
    api_data %>%
      filter(site == site_arg & parameter == parameter_arg) %>%
      select(-name) %>%
      distinct() %>%
      group_by(DT_round) %>% # site & parameter does not need to be here anymore
      # to do: preserve values used with nest()
      summarize(mean = as.numeric(mean(value, na.rm = T)),
                diff = abs(min(value, na.rm = T) - max(value, na.rm = T)),
                n_obs = n()) %>%
      ungroup() %>%
      arrange(DT_round) %>%
      # pad the dataset so that all 15-min timestamps are present
      pad(by = "DT_round", interval = "15 min") %>%
      mutate(DT_join = as.character(DT_round),
             site = site_arg,
             parameter = parameter_arg,
             flag = NA) %>%
      full_join(filter(dplyr::select(deployment_record, sonde_employed, last_site_visit, visit_comments, sensor_malfunction, DT_join, site)),
                by = c('DT_join', 'site')) %>%
      # Use fill() to determine when sonde was in the field, and when the last site visit was.
      fill(c(sonde_employed, last_site_visit))
  },

  error = function(err) {
    # error message
    cat("An error occurred with site ", site_arg, " parameter ", parameter_arg, ".\n")
    cat("Error message:", conditionMessage(err), "\n")
    flush.console() # Immediately print the error messages
    NULL  # Return NULL in case of an error
  })

  return(summary)
}

# Make a list of the summarized data
all_data_summary_list <- map2(combinations$sites,
                              combinations$params,
                              summarize_site_param,
                              api_data = all_data) %>%
  # set the names for the dfs in the list
  set_names(paste0(combinations$sites, "-", combinations$params)) %>%
  # remove NULL values from the list
  keep(~ !is.null(.))

# Bind rows for each df in list
all_data_summary_df <- bind_rows(all_data_summary_list)
```

```{r}
generate_summary_statistics <- function(site_param_df) {
  # This should get wrapped up in a function to map over all_data_summary_list
  summary_stats_df <- site_param_df %>%
    # ... so that we can get the proper leading/lagging values across our entire timeseries:
    mutate(
      # Add the next value and previous value for mean.
      front1 = lead(mean, n = 1),
      back1 = lag(mean, n = 1),
      # Add the median for a point centered in a rolling median of 7 points.
      rollmed = roll_median(mean, n = 7, align = 'center', na.rm = F, fill = NA_real_),
      # Add the mean for a point centered in a rolling mean of 7 points.
      rollavg = roll_mean(mean, n = 7, align = 'center', na.rm = F, fill = NA_real_),
      # Add the standard deviation for a point centered in a rolling mean of 7 points.
      rollsd = roll_sd(mean, n = 7, align = 'center', na.rm = F, fill = NA_real_),
      # Determine the slope of a point in relation to the point ahead and behind.
      slope_ahead = abs(front1 - mean)/15,
      slope_behind = abs(mean - back1)/15,
      # add some summary info for future us
      month = month(DT_round),
      year = year(DT_round),
      y_m = paste(year, '-', month)
    ) %>%
    group_by(y_m) %>%
    mutate(ym_sd = sd(mean, na.rm = T)) %>%
    ungroup() # does this do anything for downstream analysis?

  return(summary_stats_df)
}
all_data_summary_stats_df <- map(all_data_summary_list, generate_summary_statistics)
```

#### Set thresholds (these thresholds need to become a yaml table)

#### Add the flags
Add flagging functions for each df in all_data_summary_list

```{r}
# flag addition function
# This function will be called inside of add_x_flags() functions (ex. add_field_flags())
# This function will be used to simply add flags into the flag column
add_flag <- function(df, condition_arg, description_arg) {
  df <- df %>% mutate(flag = case_when(
    {{condition_arg}} ~ if_else(is.na(flag), description_arg, paste(flag, description_arg, sep = "; ")),
    TRUE ~ flag))
  return(df)
}

# Flag alter function
# This function can be used to alter flags that have been established
alter_flag <- function(df, condition_arg, old_description_arg, new_description_arg) {
  df <- df %>% mutate(flag = case_when(
    {{condition_arg}} & str_detect(flag, old_description_arg) ~ str_replace(flag, old_description_arg, new_description_arg),
    TRUE ~ flag))
  return(df)
}

# Example with field flags

# creating add_field_flags function
add_field_flags <- function(df) {
  df <- df %>%
    # To use add_flag in a pipeline just input a condition for a flag and your description for the flag.
    add_flag(sonde_employed == 1, "sonde not employed") %>%
    add_flag(!is.na(visit_comments), "site visit")
  return(df)
}

# applying add_field_flags() function to temperature_archery_summary_stats df
temperature_archery_field_flags <- add_field_flags(temperature_archery_summary_stats)
```

```{r}
parameter_ranges <- yaml::read_yaml("src/qaqc/parameter_thresholds.yml")

# sensor spec ranges flagging function
# Some parameters are derived from a combination of other parameters. If any of those are wrong then they should be flagged.
add_range_flags <- function(df) {

  # get the parameter from the parameter column in the df of interest
  parameter_name <- unique(na.omit(df$parameter))
  # Pull the sensor specification range from the yaml file
  sensor_min <- parameter_ranges[[parameter_name]]$sensor_specifications$min
  sensor_max <- parameter_ranges[[parameter_name]]$sensor_specifications$max
  # Pull the lab bound range from the yaml file
  lab_min <- parameter_ranges[[parameter_name]]$lab_bounds$min
  lab_max <- parameter_ranges[[parameter_name]]$lab_bounds$max

  df <- df %>%
    # adding sensor range flags
    add_flag(parameter == parameter_name & (mean < sensor_min | mean > sensor_max),
             paste0("out of (", parameter_name,") sensor specification range (TEST)")) %>%
    # adding lab bound flags
    add_flag(parameter == parameter_name & (mean < lab_min | mean > lab_max),
             paste0("out of (", parameter_name,") lab bounds (TEST)"))

  return(df)

}

# applying add_sensor_spec_range_flags() function to temperature_archery_field_flags df
temperature_archery_range_flags <- add_sensor_spec_range_flags(temperature_archery_field_flags)
```

For wonky data next to each other (defined by roll sd) - here we make some seasonal summaries to help with this:

  ```{r}
# flag data points that are outside the range of 3 sds from the mean of previous + future 3 data points:
# Do we want to play around with using the monthly sd instead of the rolling sd?

# keep this df name temperature_archery_bound_flags bc the bounding flagging
# function was the last one that was applied

temperature_archery_bound_flags <- temperature_archery_bound_flags %>%
  add_flag((mean <= rollavg - (3 * ym_sd) | mean >= rollavg + (3 * ym_sd)), "Outside SD range")
```

For data that repeats:

  ```{r}
# flag data that repeat across time:
temperature_archery_bound_flags <- temperature_archery_bound_flags %>%
  add_flag((mean == front1 | mean == back1), "Repeated value")
```

For instances where change in a 15m period is > 1degC (for temperature data) (this should also be a look up table)

```{r}
# This is here to explore data easily

test_slope_df_2 <- all_data_summary_stats_df[["legacy-pH"]] %>%
  add_field_flags() %>%
  add_range_flags() %>%
  add_flag((mean <= rollavg - (3 * ym_sd) | mean >= rollavg + (3 * ym_sd)), "Outside SD range") %>%
  add_flag((mean == front1 | mean == back1), "Repeated value") %>%
  add_flag((slope_ahead >= slope_thresh | slope_behind >= slope_thresh), "slope flag suspect") %>%
  alter_flag((slope_ahead >= slope_thresh & slope_behind >= slope_thresh), "slope flag suspect", "slope flag actual") %>%
  add_flag(is.na(mean), "missing data")
```


```{r}
# flag data where change in a 15m period is > 1degC
temperature_archery_bound_flags <- temperature_archery_bound_flags %>%
  add_flag((slope_ahead >= slope_thresh | slope_behind >= slope_thresh), "Slope exceedance")

# Testing
```

For missing data:
  ```{r}
temperature_archery_bound_flags <- temperature_archery_bound_flags %>%
  add_flag(is.na(mean), "missing data")
```

#### Visualize the flags

```{r}
# Current vis:
ggplot() +
  geom_point(data = filter(temperature_archery_bound_flags, is.na(flag)),
             aes(x=DT_round, y = mean)) +
  geom_point(data = filter(temperature_archery_bound_flags, !is.na(flag)),
             aes(x=DT_round, y = mean, color = flag)) +
  theme_bw() +
  theme(legend.position = 'bottom')
```

```{r}
# Visualizing slope exceedance
plot_list <- list()

# flag_dates <- test_slope_df_2 %>%
#   filter(str_detect(flag, "slope flag suspect"))

flag_dates <- test_slope_df_2 %>%
  filter(str_detect(flag, "site visit")) %>%
  # filter(str_detect(flag, "slope flag suspect")) %>%
  # filter(str_detect(flag, "Outside SD range")) %>%
  # filter(!is.na(flag)) %>%
  group_by(day(DT_join)) %>%
  slice(1)

for (i in 1:nrow(flag_dates)) {

  flag_year <- flag_dates$year[i]
  flag_month <- flag_dates$month[i]
  flag_day <- flag_dates$DT_round[i]


  day_data <- test_slope_df_2 %>%
    filter(year == flag_year,
           month == flag_month,
           day(DT_round) == day(flag_day))

  plot <- ggplot(data = day_data, aes(x=DT_round, y = mean, color = flag)) +
    geom_point() +
    # exceeding sd visualized
    geom_ribbon(aes(ymin = rollavg - ym_sd, ymax = rollavg + ym_sd), alpha = 0.1, color = NA) +
    geom_ribbon(aes(ymin = rollavg - (ym_sd*2), ymax = rollavg + (ym_sd*2)), alpha = 0.1, color = NA) +
    geom_ribbon(aes(ymin = rollavg - (ym_sd*3), ymax = rollavg + (ym_sd*3)), alpha = 0.1, color = NA) +
    geom_line(aes(x = DT_round, y = rollavg, color = "mean"), show.legend = TRUE) +
    # exceeding slope visualized
    geom_vline(data = (day_data %>% filter(is.na(mean))), aes(xintercept = DT_round, color = flag)) +
    # geom_smooth(method='loess', formula= y~x, color = "blue", level = .99, na.rm = TRUE) +
    theme_bw() +
    theme(legend.position = 'bottom') +
    ggtitle(paste("Site Visit on", as.character(flag_day)))

  plot_list[[as.character(flag_day)]] <- plot

}

# for loop here because map was acting strange
for (i in 1:length(plot_list)) {
  print(plot_list[[i]])
}
```

```{r}
# This can be used to develop site specific flags for specific parameters
# (later)

# Use grep to find elements with "archery" in their title
archery_objects <- all_data_summary_list[grep("archery", names(all_data_summary_list), ignore.case = TRUE)]

# Resulting list of objects with "archery" in the title
archery_objects_df <- bind_rows(archery_objects)
```
