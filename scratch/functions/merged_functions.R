add_depth_shift_flag <- function(df){


  cal_flag <- function(df, condition_arg, description_arg) {
    df <- df %>% mutate(depth_change = case_when(
      {{condition_arg}} ~ if_else(is.na(depth_change), paste(description_arg),
                                  ifelse(!grepl(description_arg, depth_change), paste(depth_change, description_arg, sep = ";\n"), depth_change)),
      TRUE ~ depth_change))
    return(df)
  }

  # Function to add a column if it doesn't already exist
  add_column_if_not_exists <- function(df, column_name, default_value) {
    if (!column_name %in% colnames(df)) {
      df <- df %>% mutate(!!sym(column_name) := default_value)
    }
    return(df)
  }

    shift_dates <- read_csv('data/qaqc/level_shifts.csv') %>%
      filter(type == "sonde moved") %>%
      add_column_if_not_exists(., "depth_change", NA) %>%
      mutate(DT_join = as.character(ymd_hms(DT_join))) %>%
      left_join(df, ., by = c("site", "DT_join")) %>%
      cal_flag(type == "sonde moved", "sonde moved") %>%
      select(-type)

    return(shift_dates)

}

add_do_drops <- function(df){

if("DO" %in% df$parameter){

  df <- df %>%
    add_flag(back1 - mean >= 0.5 & front1 - mean >= 0.5, "DO interference")

  return(df)

} else {

  return(df)

}

  }
add_drift_flag <- function(df){

  # Only test for biofilm growth on turbidity sensors
  if("Turbidity" %in% df$parameter){

    # subset data to turbidity and conductivity only
    sub <- df %>%
      dplyr::filter(parameter %in% c("Turbidity", "Specific Conductivity")) %>%
      dplyr::select(DT_round, DT_join, parameter, mean) %>%
      tidyr::pivot_wider(names_from = parameter, values_from = mean)
    names(sub) <- make.names(names(sub))

    # Check if a given window's worth of sensor data has a linear
    # relationship with time (indicates a "steady" drift)
    progressive_drift <- function(x) {
      # Only assess time series with less than 90% missing data in
      # the rolling window:
      if(length(x[!is.na(x)]) > (length(x) - (length(x)*0.1))){

        # Fit linear model
        model <- lm(x ~ c(1:length(x)), na.action = na.omit)

        # Extract R-squared value
        r_squared <- summary(model)$r.squared

        # Return R-squared value
        return(r_squared)

      } else {

        # If not enough data, set slope to 0
        no_slope <- 0

        return(no_slope)

      }
    }

    # Function to check if a selected time window's mean R-squared value is
    # at least 60% (i.e., fairly high linear relationship with time indicates a "steady" drift)
    check_too_steady <- function(x) {
      mean(x) >= 0.60
    }

    # Function that uses all the functions above to see if a given time window's R-squared with time is strong. If the one-day OR the three-day slope
    # for a selected parameter is greater than 60%, we determine it has "failed" (i.e., drift seems to exist)

    biofilm_tester <- function(data = sub, col){

      data %>%
        data.table::data.table() %>%
        #bind_rows(flagged_data_dfs[["prospect-Turbidity"]]) %>%
        # dplyr::filter(lubridate::as_date(DT_round) >= "2022-09-10" & lubridate::as_date(DT_round) <= "2022-10-15") %>%
        dplyr::mutate(r2_s_right = data.table::frollapply(!!sym(col), n = 96, FUN = progressive_drift, align = "right", fill = NA),
                      r2_s_center = data.table::frollapply(!!sym(col), n = 96, FUN = progressive_drift, align = "left", fill = NA),
                      r2_l_right = data.table::frollapply(!!sym(col), n = 288, FUN = progressive_drift, align = "right", fill = NA),
                      r2_l_center = data.table::frollapply(!!sym(col), n = 288, FUN = progressive_drift, align = "left", fill = NA),
                      tightest_r = pmax(r2_s_center, r2_s_right, r2_l_center, r2_l_right, na.rm = TRUE),
                      failed = data.table::frollapply(tightest_r, n = 96, FUN = check_too_steady, align = "right", fill = NA)) %>%
        dplyr::select(!!(col) := "failed",
                      DT_join)
    }

    biofouling <- c("Turbidity", "Specific.Conductivity") %>%
      purrr::map(~biofilm_tester(data = sub, col = .)) %>%
      dplyr::bind_cols() %>%
      dplyr::rename(DT_join = DT_join...2) %>%
      dplyr::mutate(parameter = "Turbidity") %>%
      dplyr::right_join(., df, by = c("DT_join", "parameter")) %>%
      # If a steady slope is existing through time for turbidity, but NOT for conductivity, it is likely
      # sensor bio-fouling
      add_flag(., Turbidity == 1 & Specific.Conductivity != 1 & !grepl("drift", flag), "drift") %>%
      dplyr::select(-c(Turbidity, DT_join...4, Specific.Conductivity))

    return(biofouling)

  } else {

    return(df)}

}
#
# a=ggplot(drift_tested)+
#   geom_point(aes(x = DT_round, y = r_squaredSC_short ), color = "blue") +
#   geom_point(aes(x = DT_round, y = r_squaredT_short ), color = "red") +
#   geom_point(aes(x = DT_round, y = r_squaredDO_short), color = "black")
#
# b=ggplot(drift_tested)+
#   geom_point(aes(x = DT_round, y = r_squaredSC_long ), color = "blue") +
#   geom_point(aes(x = DT_round, y = r_squaredT_long ), color = "red") +
#   geom_point(aes(x = DT_round, y = r_squaredDO_long), color = "black")
#
# c=ggplot(drift_tested %>% filter(as_date(DT_round) >= "2022-10-10")) +
#   #geom_point(aes(x = DT_round, y = Specific.Conductivity)) +
#   geom_point(aes(x = DT_round, y = Turbidity, color = as.character(failed))) #+
# #geom_point(aes(x = DT_round, y = DO))
#
# ggpubr::ggarrange(a, b, c, ncol = 1)
#

#' @title Add field related flags to a data frame based on field notes.
#'
#' @description
#' A function that checks 3 different and separate conditions related to the
#' field notes and adds the corresponding flags accordingly. The "sonde not
#' employed" flag is added if the sonde was not employed at the site.
#' "site visit" flag is added if the last site visit date is the same as the
#' date of the data. "sv window" flag is added if site visit flag is detected
#' within a 45 minute window. Note that this function will only work on data
#' that has already been joined to field notes.
#'
#' @param df A data frame with a `flag` column.
#'
#' @return A data frame with a `flag` column that has been updated with the
#' relevant field note flags.
#'
#' @examples
#' add_range_flags(df = all_data_flagged$`archery-Actual Conductivity`)
#' add_range_flags(df = all_data_flagged$`boxelder-Temperature`)

add_field_flag <- function(df) {

  df <- df %>%
    # flag when sonde was not employed in the river
    add_flag(sonde_employed == 1, "sonde not employed") %>% # removing the instances where we check for the flag description in the flag column when we are flagging. This is already done in `add_flag()`
    # flag when sonde was handled in a site visit
    add_flag(as.character(last_site_visit) == as.character(DT_round), "site visit")

  # Add flags for the next 60 minutes after a site visit
  for (i in 1:4) {
    df <- df %>%
      add_flag(lag(str_detect(flag, "site visit"), n = i), "sv window")
  }

  # ... and the first 15 minutes before:
  df <- df %>%
    add_flag(lead(str_detect(flag, "site visit"), n = 1), "sv window")

  return(df)

}
#' @title Underlying function for flagging data.
#'
#' @description
#' This function adds a flag to the `flag` column of a given data frame based on
#' specified conditions for each row. The name of the flag is provided
#' by the user.
#'
#' @param df A data frame with a `flag` column.
#'
#' @param condition_arg A logical statement that is evaluated in the context of
#' the data frame.
#'
#' @param description_arg Flag added to the `flag` column.
#'
#' @returns
#' An identical data frame with a `flag` column that has been updated with the
#' flag description provided.
#'
#' @examples
# add_flag(df = all_data_flagged$`archery-Actual Conductivity`, condition_arg = mean >= 100, description_arg = "exceeds 100")
# add_flag(df = all_data_flagged$`archery-Actual Conductivity`, condition_arg = mean <= 10, description_arg = "below 10")

add_flag <- function(df, condition_arg, description_arg) {
  df <- df %>% mutate(flag = case_when(
    {{condition_arg}} ~ if_else(is.na(flag), paste(description_arg),
                                ifelse(!grepl(description_arg, flag), paste(flag, description_arg, sep = ";\n"), flag)),
    TRUE ~ flag))
  return(df)
}
#' @title Add a flag if the water temperature is freezing.
#'
#' @description
#' A function designed to append the 'frozen' flag to a row if the value
#' in the `mean` column is less than or equal to 0.
#'
#' @param df A data frame with a `flag` column.
#'
#' @return A data frame with a `flag` column that has been updated with the
#' 'frozen' flag.
#'
#' @examples
#' add_frozen_flag(df = all_data_flagged$`archery-Actual Conductivity`)
#' add_frozen_flag(df = all_data_flagged$`boxelder-Temperature`)
#'
#' @seealso [flag_all_data()]

add_frozen_flag <- function(df){

  # create a df of temperature for each site
  temperature <- df %>%
    data.table::data.table() %>%
    dplyr::select(DT_round, DT_join, parameter, mean) %>%
    dplyr::filter(parameter == "Temperature") %>%
    dplyr::select(DT_join, Temperature = mean)

    # add "temperature" column to df:
    temperature_checked <- df %>%
      dplyr::left_join(., temperature, by = "DT_join") %>%
      # If water temperature is freezing, flag all parameters
      add_flag(., Temperature <= 0, "frozen") %>%
      # remove the temp column so df is identical in structure to OG df
      dplyr::select(-Temperature)

    return(temperature_checked)

}
#' @title Add malfunction flags to a data frame based on field notes and HydroVu.
#'
#' @description
#' This function adds the 'sensor malfunction' flag to a data frame by the dates,
#' sites, and parameters that members of the lab know contain erroneous data due
#' to sensor malfunctions. Note that this flag is used in two instances: first
#' when removing erroneous data from our statistics calculations and again during
#' the actual flagging step in the QAQC pipeline.
#'
#' @param df A data frame with a `flag` column.
#'
#' @param malfunction_records The malfunction records pulled from mWater
#'
#' @return A data frame with a `flag` column that has been updated with the
#' 'sensor malfunction' flag.
#'
#' @examples
#' add_range_flags(df = all_data_flagged$`archery-Actual Conductivity`)
#' add_range_flags(df = all_data_flagged$`boxelder-Temperature`)
#'
#' @seealso [grab_mWater_malfunction_notes()]

add_malfunction_flag <- function(df, malfunction_records){

  # Filter malfunction records for relevant information
  df_site <- unique(df$site)
  df_parameter <- unique(df$parameter)

  malfunction_records_filtered <- malfunction_records %>%
    filter(site == df_site) %>%
    # if parameter is NA, that means entire sonde was malfunctioning. Otherwise,
    # just the parameter listed:
    filter(is.na(parameter) | parameter == df_parameter) %>%
    mutate(end_DT = ifelse(is.na(end_DT), ymd_hms("9999-12-31 23:59:59", tz = "MST"), end_DT)) %>%
    mutate(end_DT = as.POSIXct(end_DT, tz = "MST")) %>%
    rowid_to_column()

  if (nrow(malfunction_records_filtered > 0)) {

    drift <- malfunction_records_filtered %>%
      filter(grepl("grime|gunk|drift|biofoul|biofilm", notes, ignore.case = TRUE))

    burial <- malfunction_records_filtered %>%
      filter(grepl("buried|burial|bury|sediment|roots", notes, ignore.case = TRUE))

    depth_funk <- malfunction_records_filtered %>%
      filter(grepl("improper level calibration", notes, ignore.case = TRUE))

    unsubmerged <- malfunction_records_filtered %>%
      filter(grepl("not submerged", notes, ignore.case = TRUE))

    general_malfunction <- malfunction_records_filtered %>%
      filter(!rowid %in% drift$rowid & !rowid %in% burial$rowid &
               !rowid %in% depth_funk$rowid & !rowid %in% unsubmerged$rowid)

    # make a list of date intervals where sensor was malfunctioning
    drift_interval_list <- map2(
      .x = drift$start_DT,
      .y = drift$end_DT,
      .f = ~interval(.x, .y, tz = "MST"))

    burial_interval_list <- map2(
      .x = burial$start_DT,
      .y = burial$end_DT,
      .f = ~interval(.x, .y, tz = "MST"))

    depth_interval_list <- map2(
      .x = depth_funk$start_DT,
      .y = depth_funk$end_DT,
      .f = ~interval(.x, .y, tz = "MST"))

    unsubmerged_interval_list <- map2(
      .x = unsubmerged$start_DT,
      .y = unsubmerged$end_DT,
      .f = ~interval(.x, .y, tz = "MST"))

    general_interval_list <- map2(
      .x = general_malfunction$start_DT,
      .y = general_malfunction$end_DT,
      .f = ~interval(.x, .y, tz = "MST"))

    # If DT_round is within any of the intervals where a sensor was malfunctioning, flag it
    try(df <- df %>%
          add_flag(DT_round %within% burial_interval_list, "sonde burial"))

    try(df <- df %>%
          add_flag(DT_round %within% drift_interval_list, "sensor biofouling"))

    try(df <- df %>%
          add_flag(DT_round %within% depth_interval_list, "depth calibration malfunction"))

    try(df <- df %>%
          add_flag(DT_round %within% unsubmerged_interval_list, "sonde unsubmerged"))

    try(df <- df %>%
          add_flag(DT_round %within% general_interval_list, "sensor malfunction"))

  }


  return(df)

}
#' @title Add NA flags to a data frame based on `mean` column.
#'
#' @description
#' A function designed to append the 'missing data' flag to a row if the `mean`
#' column in that row contains NA values.
#'
#' @param df A data frame with a `flag` column.
#'
#' @return A data frame with a `flag` column that has been updated with the
#' 'missing data' flag.
#'
#' @examples
#' add_range_flags(df = all_data_flagged$`archery-Actual Conductivity`)
#' add_range_flags(df = all_data_flagged$`boxelder-Temperature`)

add_na_flag <- function(df){
  df <- df %>%
    add_flag(is.na(mean), "missing data")
}
#' @title Add a flag if the value of a parameter exceeds its realistic range.
#'
#' @description
#' A function designed to append the 'outside of sensor realistic range' flag to
#' a row if the `mean` column in that row exceeds the realistic ranges that are
#' set in the `src/qaqc/sensor_real_thresholds.yml` file. These are instances
#' where a value exceeds the expected ranges for the Poudre. Note that this is
#' currently only aplied to pH, SC, and Temperature.
#'
#' @param df A data frame with a `flag` column.
#'
#' @return A data frame with a `flag` column that has been updated with the
#' 'outside of sensor realistic range' flag.
#'
#' add_range_flags(df = all_data_flagged$`archery-Actual Conductivity`)
#' add_range_flags(df = all_data_flagged$`boxelder-Temperature`)

add_realistic_flag <- function(df){

  sensor_realistic_ranges <- read_yaml("data/qaqc/sensor_real_thresholds.yml")

  # make this a non yaml solution and add it to the threshold table
  # get the parameter from the parameter column in the df of interest
  parameter_name <- unique(na.omit(df$parameter))
  # Pull the sensor specification range from the yaml file
  sensor_min <- eval(parse(text = sensor_realistic_ranges[[parameter_name]]$min))
  sensor_max <- eval(parse(text = sensor_realistic_ranges[[parameter_name]]$max))

  df <- df %>%
    # adding sensor range flags
    add_flag(parameter == parameter_name & (mean < sensor_min | mean > sensor_max),
             "outside of sensor realistic range")

    return(df)

}
#' @title Add a flag if the value in the `mean` column repeats sequentially.
#'
#' @description
#' A function designed to append the 'repeated value' flag to a row if the value
#' in the `mean` column is equal to the previous or next value in the `mean`
#' column.
#'
#' @param df A data frame with a `flag` column.
#'
#' @return A data frame with a `flag` column that has been updated with the
#' 'repeated value' flag.
#'
#' add_range_flags(df = all_data_flagged$`archery-Actual Conductivity`)
#' add_range_flags(df = all_data_flagged$`boxelder-Temperature`)

add_repeat_flag <- function(df){
  df <- df %>%
    add_flag((mean == front1 | mean == back1), "repeated value")
}
#' @title Add flags related to calculated parameter seasonal ranges.
#'
#' @description
#' A function that checks 2 different and separate conditions related to the
#' calculated parameter seasonal ranges. The 'outside of seasonal range' flag is
#' added if the `mean` value is outside the seasonal 1st - 99th percentile. The
#' 'slope violation' flag is added if the `slope_ahead` or `slope_behind` value
#' is greater than or equal to the `t_slope_behind_99` value, or if the
#' `slope_ahead` or `slope_behind` value is less than or equal to the
#' `t_slope_behind_01` value. The `t_slope_behind_99` and `t_slope_behind_01`
#' values are specific to each site-parameter. Note that this is the version of
#' the function used in `flag_all_data()`.
#'
#' @param df A data frame with a `flag` column.
#'
#' @return A data frame with a `flag` column that has been updated with the
#' relevant calculated seasonal range flags.
#'
#' @examples
#' add_range_flags(df = all_data_flagged$`archery-Actual Conductivity`)
#' add_range_flags(df = all_data_flagged$`boxelder-Temperature`)
#'
#' @seealso [flag_all_data()]

add_seasonal_flag <- function(df) {

  # get the site name from the site column in the df of interest
  site_name <- unique(na.omit(df$site))
  # get the parameter from the parameter column in the df of interest
  parameter_name <- unique(na.omit(df$parameter))

  lookup <- threshold_lookup %>%
    filter(site == site_name & parameter == parameter_name) %>%
    select(!c(site, parameter))

  df <- df %>%
    # Using seasonal cut-offs...
    left_join(lookup, by = "season") %>%
    # ... flag obs that are outside the seasonal 1st - 99th percentile range:
    add_flag((mean < t_mean01 | mean > t_mean99), "outside of seasonal range") %>%
    # flag obs whose slope is outside the 1st - 99th percentile range:
    add_flag(((slope_ahead >= t_slope_behind_99 | slope_behind >= t_slope_behind_99) |
              (slope_ahead <= t_slope_behind_01 | slope_behind <= t_slope_behind_01)), "slope violation")

  return(df)

}
# Adding flags related to sensor specification ranges; these are instances where
# a value exceeds the manufacturer specified limits.

#' @title Add a flag if the value of a parameter exceeds its sensor specification range.
#'
#' @description
#' A function designed to append the 'outside of sensor specification range' flag to
#' a row if the `mean` column in that row exceeds the sensor specification ranges that are
#' set in the `src/qaqc/sensor_spec_thresholds.yml` file. These are instances
#' where a value exceeds the expected ranges for the Poudre based on the sensor
#' manufacturer's specifications.
#'
#' @param df A data frame with a `flag` column.
#'
#' @return A data frame with a `flag` column that has been updated with the
#' 'outside of sensor specification range' flag.
#'
#' @examples
#' add_spec_flag(df = all_data_flagged$`archery-Actual Conductivity`)
#' add_spec_flag(df = all_data_flagged$`boxelder-Temperature`)
#'
#' @seealso [flag_all_data()]

add_spec_flag <- function(df){

  sensor_spec_ranges <- read_yaml("data/qaqc/sensor_spec_thresholds.yml")

  # make this a non yaml solution and add it to the threshold table
  # get the parameter from the parameter column in the df of interest
  parameter_name <- unique(na.omit(df$parameter))
  # Pull the sensor specification range from the yaml file
  sensor_min <- eval(parse(text = sensor_spec_ranges[[parameter_name]]$min))
  sensor_max <- eval(parse(text = sensor_spec_ranges[[parameter_name]]$max))

  df <- df %>%
    # adding sensor range flags
    add_flag(parameter == parameter_name & (mean < sensor_min | mean > sensor_max) & !grepl("outside of sensor specification range", flag),
             "outside of sensor specification range") %>%

    return(df)

}
# Add flags to the `flag` column of a dataframe based on large swaths of suspect data.
# "24hr anomaly flag" is added if more than 50% of the data points in a 24 hour window are flagged.
# "anomaly window" flag is added if the point is included in a 24hr anomaly.
# @param df A dataframe with a `auto_cleaned_flag` column.
# @return A dataframe with a `auto_cleaned_flag` column that has been updated with the relevant large anomaly flags.
# @examples
# add_large_anomaly_flags(df = all_data_flagged$`archery-Actual Conductivity`)
# add_large_anomaly_flags(df = all_data_flagged$`boxelder-Temperature`)

add_suspect_flag_full <- function(df) {
  # these are the flags that we don't want to perform this exercise across
  auto_cleaned_flag_string <- "sonde not employed|missing data|site visit|sv window|sonde burial|sensor biofouling|depth calibration malfunction|sonde unsubmerged|sensor malfunction"

  # Define a function to check if a given 3-hour window has >= 50% fails
  check_3_hour_window_fail <- function(x) {
    sum(x) / length(x) >= 0.5
  }

  df_test <- df %>%
    dplyr::mutate(auto_cleaned_flag_binary = ifelse((is.na(auto_cleaned_flag) | grepl(auto_cleaned_flag_string, auto_cleaned_flag) | auto_cleaned_flag_string == "suspect data"), 0, 1)) %>%
    #arrange(timestamp) %>%
    dplyr::mutate(over_50_percent_fail_window = zoo::rollapply(auto_cleaned_flag_binary, width = 12, FUN = check_3_hour_window_fail, fill = NA, align = "right")) %>%
    add_flag(over_50_percent_fail_window == TRUE & is.na(auto_cleaned_flag), "suspect data")

  return(df_test)

}


#' @title Add flags to the `flag` column of a dataframe based on large swaths of suspect data.
#'
#' @description
#' "24hr anomaly flag" is added if more than 50% of the data points in a 24 hour window are flagged.
#' "anomaly window" flag is added if the point is included in a 24hr anomaly.
#'
#' @param df A data frame with a `flag` column.
#'
#' @return A data frame with a `flag` column that has been updated with the relevant calculated seasonal range flags.
#'
#' @examples
#' add_range_flags(df = all_data_flagged$`archery-Actual Conductivity`)
#' add_range_flags(df = all_data_flagged$`boxelder-Temperature`)
#'
#' @seealso [flag_all_data()]

add_suspect_flag <- function(df) {

  flag_string <- "sonde not employed|missing data|site visit|sv window|suspect data" # include the flag added by this function

  # Define a function to check if a given 3-hour window has >= 50% fails
  check_3_hour_window_fail <- function(x) {
    sum(x) / length(x) >= 0.5
  }

  df_test <- df %>%
    dplyr::mutate(flag_binary = ifelse((is.na(flag) | grepl(flag_string, flag)), 0, 1)) %>%
    #arrange(timestamp) %>%
    dplyr::mutate(over_50_percent_fail_window = ifelse(is.na(over_50_percent_fail_window),
                                                zoo::rollapply(flag_binary, width = 12, FUN = check_3_hour_window_fail, fill = NA, align = "right"),
                                                over_50_percent_fail_window)) %>%
    add_flag(over_50_percent_fail_window == TRUE & is.na(flag), "suspect data")

  return(df_test)

}
add_threshold_lines <- function(plot, plot_data, site_arg, parameter_arg) {

  # pull in threshold data (i don't like that I do this everytime the function is run)
  real_thresholds <- read_csv("data/qaqc/realistic_thresholds.csv", show_col_types = FALSE) %>%
    filter(parameter == parameter_arg)
  sensor_thresholds <- yaml::read_yaml("data/qaqc/sensor_spec_thresholds.yml")[[parameter_arg]]%>% #filter for parameter_arg
    #turn into tibble with min/max
    bind_rows()

  seasonal_thresholds <- read_csv('data/qaqc/seasonal_thresholds_virridy.csv', show_col_types = FALSE) %>%
    #to do: Check to make sure seasonal thresholds csv is not necessary
    #bind_rows(read_csv('data/qaqc/seasonal_thresholds.csv', show_col_type = FALSE),
    distinct(site, parameter, season, .keep_all = TRUE) %>%
    #read_csv("data/qaqc/seasonal_thresholds_virridy.csv", show_col_types = FALSE) %>%
    filter(parameter == parameter_arg,
           site == site_arg)

  # Determine the unique seasons in plot_data
  unique_seasons <- unique(plot_data$season)

  # Filter seasonal_thresholds for the seasons present in plot_data
  seasonal_thresholds <- seasonal_thresholds %>%
    filter(season %in% unique_seasons)

  if (nrow(seasonal_thresholds) > 1) { # make sure this works

    season_1 <- case_when(
      all(c("winter_baseflow", "snowmelt") %in% unique_seasons) ~ "winter_baseflow",
      all(c("snowmelt", "monsoon") %in% unique_seasons) ~ "snowmelt",
      all(c("monsoon", "fall_baseflow") %in% unique_seasons) ~ "monsoon",
      all(c("fall_baseflow", "winter_baseflow") %in% unique_seasons) ~ "fall_baseflow",
      TRUE ~ NA_character_
    )
    season_2 <- case_when(
      all(c("winter_baseflow", "snowmelt") %in% unique_seasons) ~ "snowmelt",
      all(c("snowmelt", "monsoon") %in% unique_seasons) ~ "monsoon",
      all(c("monsoon", "fall_baseflow") %in% unique_seasons) ~ "fall_baseflow",
      all(c("fall_baseflow", "winter_baseflow") %in% unique_seasons) ~ "winter_baseflow",
      TRUE ~ NA_character_
    )

    seasonal_threshold_s1_quantiles <- unname(quantile(c(seasonal_thresholds %>% filter(season == season_1) %>% pull(t_mean01),
                                                         seasonal_thresholds %>% filter(season == season_1) %>% pull(t_mean99)),
                                                       c(0.1, 0.9)))

    seasonal_threshold_s2_quantiles <- unname(quantile(c(seasonal_thresholds %>% filter(season == season_2) %>% pull(t_mean01),
                                                         seasonal_thresholds %>% filter(season == season_2) %>% pull(t_mean99)),
                                                       c(0.1, 0.9)))

    slice_data <- plot_data %>%
      group_by(season) %>%
      slice(1) %>%
      ungroup()

    site_data <- filter(plot_data, site == site_arg)

    if (!all(is.na(site_data$mean))) {
      # Lower bound
      if ( # !is.infinite(min(site_data$mean, na.rm = TRUE)) &
        (min(site_data$mean, na.rm = TRUE) <= seasonal_threshold_s1_quantiles[1] | min(site_data$mean, na.rm = TRUE) <= seasonal_threshold_s2_quantiles[1])) {

        # Xs
        start_x_s1 <- min(slice_data$DT_round)

        transition_date <- slice_data %>%
          filter(season == season_2) %>%
          pull(DT_round)

        end_x_s2 <- ceiling_date(max(plot_data$DT_round), "day")

        #Ys
        y_lower_s1 <- seasonal_thresholds %>%
          filter(season == season_1) %>%
          pull(t_mean01)

        y_lower_s2 <- seasonal_thresholds %>%
          filter(season == season_2) %>%
          pull(t_mean01)

        plot <- plot +
          # season 1
          ## lower bound
          geom_segment(aes(x = start_x_s1, y = y_lower_s1,
                           xend = transition_date, yend = y_lower_s1,
                           color = "Seasonal Min", linetype = "Seasonal Min")) +
          # season 2
          ## lower bound
          geom_segment(aes(x = transition_date, y = y_lower_s2,
                           xend = end_x_s2, yend = y_lower_s2,
                           color = "Seasonal Min", linetype = "Seasonal Min"))

      }

      # Upper bound
      if (# !is.infinite(max(site_data$mean, na.rm = TRUE)) &
        (max(site_data$mean, na.rm = TRUE) >= seasonal_threshold_s1_quantiles[2] | max(site_data$mean, na.rm = TRUE) >= seasonal_threshold_s2_quantiles[2])) {

        # Xs
        start_x_s1 <- min(slice_data$DT_round)


        transition_date <- slice_data %>%
          filter(season == season_2) %>%
          pull(DT_round)

        end_x_s2 <- ceiling_date(max(plot_data$DT_round), "day")

        #Ys
        y_upper_s1 <- seasonal_thresholds %>%
          filter(season == season_1) %>%
          pull(t_mean99)

        y_upper_s2 <- seasonal_thresholds %>%
          filter(season == season_2) %>%
          pull(t_mean99)

        plot <- plot +
          # season 1
          ## lower bound
          geom_segment(aes(x = start_x_s1, y = y_upper_s1,
                           xend = transition_date, yend = y_upper_s1,
                           color = "Seasonal Max", linetype = "Seasonal Max")) +
          # season 2
          ## lower bound
          geom_segment(aes(x = transition_date, y = y_upper_s2,
                           xend = end_x_s2, yend = y_upper_s2,
                           color = "Seasonal Max", linetype = "Seasonal Max"))

      }

      # real thresholds
      real_thresholds_quantiles <- unname(quantile(c(real_thresholds$min, real_thresholds$max), c(0.1, 0.9)))

      if (#!is.infinite(min(site_data$mean, na.rm = TRUE))  &
        (min(site_data$mean, na.rm = TRUE) <= real_thresholds_quantiles[1])) {
        plot <- plot +
          geom_hline(aes(yintercept = real_thresholds$min,
                         color = "Real Min",
                         linetype = "Real"))
      }
      if (#!is.infinite(max(site_data$mean, na.rm = TRUE)) &
        (max(site_data$mean, na.rm = TRUE) >= real_thresholds_quantiles[2])) {
        plot <- plot +
          geom_hline(aes(yintercept = real_thresholds$max,
                         color = "Real Max",
                         linetype = "Real"))
      }

      # sensor thresholds
      sensor_thresholds_quantiles <- unname(quantile(c(sensor_thresholds$min, sensor_thresholds$max), c(0.1, 0.9)))

      # if ((min(site_data$mean, na.rm = TRUE) <= sensor_thresholds_quantiles[1])) {
      #   plot <- plot +
      #     geom_hline(aes(yintercept = sensor_thresholds$mix, # *** this needs to be min
      #                    color = "Sensor Min",
      #                    linetype = "Sensor"))
      # }

      if ((max(site_data$mean, na.rm = TRUE) >= sensor_thresholds_quantiles[2])) {
        plot <- plot +
          geom_hline(aes(yintercept = sensor_thresholds$max,
                         color = "Sensor Max",
                         linetype = "Sensor"))
      }
      return(plot)
    }

  } else if (nrow(seasonal_thresholds) == 1){

    site_data <- filter(plot_data, site == site_arg)

    # Filter seasonal_thresholds for the seasons present in site_data
    seasonal_thresholds <- seasonal_thresholds %>%
      filter(season %in% unique(site_data$season))

    if (!all(is.na(site_data$mean))) {
      # seasonal thresholds
      seasonal_thresholds_quantiles <- unname(quantile(c(seasonal_thresholds$t_mean01, seasonal_thresholds$t_mean99), c(0.1, 0.9)))

      if ((min(site_data$mean, na.rm = TRUE) <= seasonal_thresholds_quantiles[1]) == TRUE){
        plot <- plot +
          geom_hline(aes(yintercept = seasonal_thresholds$t_mean01,
                         color = "Seasonal Min",
                         linetype = "Seasonal"))
      }

      if ((max(site_data$mean, na.rm = TRUE) >= seasonal_thresholds_quantiles[2]) == TRUE) {
        plot <- plot +
          geom_hline(aes(yintercept = seasonal_thresholds$t_mean99,
                         color = "Seasonal Max",
                         linetype = "Seasonal"))
      }

      # real thresholds
      real_thresholds_quantiles <- unname(quantile(c(real_thresholds$min, real_thresholds$max), c(0.1, 0.9)))

      if (min(site_data$mean, na.rm = TRUE) <= real_thresholds_quantiles[1]) {
        plot <- plot +
          geom_hline(aes(yintercept = real_thresholds$min,
                         color = "Real Min",
                         linetype = "Real"))
      }
      if (max(site_data$mean, na.rm = TRUE) >= real_thresholds_quantiles[2]) {
        plot <- plot +
          geom_hline(aes(yintercept = real_thresholds$min,
                         color = "Real Max",
                         linetype = "Real"))
      }

      # sensor thresholds
      sensor_thresholds_quantiles <- unname(quantile(c(sensor_thresholds$min, sensor_thresholds$max), c(0.1, 0.9)))

      # if (min(site_data$mean, na.rm = TRUE) <= sensor_thresholds_quantiles[1]) {
      #   plot <- plot +
      #     geom_hline(aes(yintercept = sensor_thresholds$mix,
      #                color = "Sensor Mine",
      #                linetype = "Sensor"))
      # }

      if (max(site_data$mean, na.rm = TRUE) >= sensor_thresholds_quantiles[2]) {
        plot <- plot +
          geom_hline(aes(yintercept = sensor_thresholds$max,
                         color = "Sensor Max",
                         linetype = "Sensor"))
      }
      return(plot)
    }
    return(plot)
  }
  return(plot)
}
#' @title Add a flag if the sonde was not fully submerged by water.
#'
#' @description
#' A function designed to append the 'unsubmerged' flag to a row if the value
#' in the `relative_depth` column is less than or equal to 0.
#'
#' @param df A data frame with a `flag` column.
#'
#' @return A data frame with a `flag` column that has been updated with the
#' 'sonde unsubmerged' flag.
#'
#' @examples
#' add_unsubmerged_flag(df = all_data_flagged$`archery-Actual Conductivity`)
#' add_unsubmerged_flag(df = all_data_flagged$`boxelder-Temperature`)
#'
#' @seealso [flag_all_data()]

add_unsubmerged_flag <- function(df){

  # create a df of temperature for each site
  depth <- df %>%
    data.table::data.table() %>%
    dplyr::select(DT_round, DT_join, parameter, mean) %>%
    dplyr::filter(parameter == "Depth") %>%
    dplyr::select(DT_join, Depth = mean)

  # add "depth" column to df:
  depth_checked <- df %>%
    dplyr::left_join(., depth, by = "DT_join") %>%
    # If water temperature is freezing, flag all parameters
    add_flag(., Depth <= 0, "sonde unsubmerged") %>%
    # remove the temp column so df is identical in structure to OG df
    dplyr::select(-Depth)

  return(depth_checked)

}
# Add a verification column to the data frame and automatically fail data points that meet certain criteria.
# @param df A dataframe with a `flag` column.
# @return A dataframe with a `verification` column.
# @examples
# add_verification_column(df = all_data_flagged$`archery-Actual Conductivity`)
# add_verification_column(df = all_data_flagged$`boxelder-Temperature`)

add_verification_column <- function(df) {
  df <- df %>%
    mutate(verification = case_when(
      # are there other conditions where a data point should automatically fail?
      # str_detect(flag, "outside sd range") ~ "fail",
      str_detect(flag, "missing data") ~ "fail",
      # need to add more situations where data points fail
      is.na(flag) ~ "pass",
      TRUE ~ NA))
  return(df)
}
alter_verification_column <- function(verification_choice, df_data) {

  df_name <- paste(unique(df_data$site), unique(df_data$parameter), head(df_data$DT_round, 1))

  # PASS/FAIL statements ----
  # these cases are inverses of each other, both are included
  # in the same if statements to reduce confusion

  # PASS ALL/FAIL NONE
  if (verification_choice %in% c("PASS ALL", "FAIL NONE")) {

    altered_df <- df_data %>%
      mutate(mean_verified = if_else(is.na(flag), mean, NA),
             is_verified = TRUE,
             verification_status = "PASS")

    cat("All points for ", df_name, "have PASSED and are accurate information.\n")

    return(altered_df)

  }

  # PASS VALID/FAIL FLAGGED
  if (verification_choice %in% c("PASS VALID", "FAIL FLAGGED")) {

    altered_df <- df_data %>%
      mutate(mean_verified = mean,
             is_verified = TRUE,
             verification_status = if_else(is.na(flag), "PASS", "FAIL")) # pass valid data, fail flagged data

    cat("All valid points for ", df_name, "have PASSED\n")
    cat("All flagged points for ", df_name, "have FAILED\n")

    return(altered_df)

  }

  # PASS FLAGGED, FAIL VALID
  # if all the flagged data is correct, and none of the valid data is correct,
  # none of it should be shown
  if (verification_choice %in% c("PASS FLAGGED", "FAIL VALID")) {

    altered_df <- df_data %>%
      mutate(mean_verified = NA,
             is_verified = TRUE,
             verification_status = if_else(!is.na(flag), "PASS", "FAIL"))

    cat("All valid points for ", df_name, "have FAILED\n")
    cat("All flagged points for ", df_name, "have PASSED.\n")

    return(altered_df)

  }

  # PASS NONE, FAIL ALL
  # This means that flagged data would actually be accurate data, and valid data would
  # actually be inaccurate data. this is the inverse of pass all
  if (verification_choice %in% c("PASS NONE", "FAIL ALL")) {

    altered_df <- df_data %>%
      mutate(mean_verified = if_else(!is.na(flag), mean, NA),
             is_verified = TRUE,
             verification_status = "FAIL")

    cat("All points for ", df_name, "have FAILED and are inaccurate information.\n")

    return(altered_df)

  }

  # skip statements
  if (verification_choice == "SKIP") {

    altered_df <- df_data %>%
      mutate(mean_verified = if_else(is.na(flag), mean, NA),
             is_verified = TRUE, #false?
             verification_status = "SKIP")

    cat("The accuracy of all points for ", df_name, "cannot be determined as\n")
    cat("accurate or inaccurate information at this time.")

    return(altered_df)

  }

  # inspect statements
  if (verification_choice %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME")) {
    stop("You can't inspect here")
  }

  # quit statements ***
  if (verification_choice == "QUIT") {

    cat("Quitting out of verification mode...")
    cat("\n")

    QUIT <<- TRUE

    return(NULL)

  }

}

# to do (j): document this function
api_puller <- function(site, start_dt, end_dt = Sys.time(), api_token, dump_dir, require = NULL) {

  locs <- hv_locations_all(hv_token)

  # make a list of site names

  options(scipen = 999)
  for(i in 1:length(site)){

    site_loc <- locs %>%
      dplyr::mutate(name = tolower(name)) %>%
      dplyr::filter(grepl(site[i], name, ignore.case = TRUE))

    site_loc_list <- site_loc$id

    # Get data for each site location. Note this maps over the entire list of locations,
    # many of which are unlikely to be active during the time you specify. Don't freak out if
    # you see a bunch of '404 Not Found' errors, you're just seeing the list of locations
    # that are not active. The data frame 'alldata' should contain your data from all applicable
    # sites during the time frame indicated. Note that this will take some time (one month of
    # data for 5 sites takes ~10 mins. Be patient!

    # Add date range you are interested in; data are stored in HydroVu in UTC
    # Here, a way to find the most recent download of the data. Use this as the start date to
    # reduce overlapping data

    # tz weirdness
    # utc_start_date <- format(as.POSIXct(start_dt, tz = "UTC") + hours(7), format = "%Y-%m-%d %H:%M:%S")
    #
    # utc_end_date <-   format(as.POSIXct(end_dt, tz = "UTC") + hours(7), format = "%Y-%m-%d %H:%M:%S")

    # doing this fixes the mismatch in date times during the combined_data step - jd
    utc_start_date <- format(as.POSIXct(start_dt, tz = "UTC") + lubridate::hours(0), format = "%Y-%m-%d %H:%M:%S")

    utc_end_date <-   format(as.POSIXct(end_dt, tz = "UTC") + lubridate::hours(0), format = "%Y-%m-%d %H:%M:%S")

    timezone = "UTC"

    # Map over the location ids
    alldata <- site_loc_list %>% purrr::map(~hv_data_id(.,
                                                        start_time = utc_start_date,
                                                        end_time = utc_end_date,
                                                        token = api_token,
                                                        tz = timezone))

    # grab only locations with data (stored as a data frame) / drop 404 errors
    filtered <- purrr::keep(alldata, is.data.frame)

    if(length(filtered) == 0){

      print(paste0("No data at ", site[i], " during this time frame"))

    } else {

      # bind lists together (now that all are dataframes, we can just collate quickly)
      one_df <- dplyr::bind_rows(filtered) %>%
        dplyr::rename(id = Location,
                      parameter = Parameter,
                      units = Units) %>%
        dplyr::left_join(., site_loc, by = "id") %>%
        dplyr::mutate(site = tolower(site[i])) %>%
        dplyr::select(site, id, name, timestamp, parameter, value, units)

      # if site contains both a csu and a virridy sonde, split them up:

      if(site[i] %in% c("Timberline", "Prospect", "Archery")){

        try(virridy_df <- one_df %>%
              filter(grepl("virridy", name, ignore.case = TRUE)) %>%
              mutate(site = paste0(site[i], " virridy")))

        try(readr::write_csv(virridy_df,
                             paste0(dump_dir, "/", site[i], " virridy_", stringr::str_replace(stringr::str_replace(substr(end_dt, 1, 16), "[:\\s]", "_"), ":", ""), ".csv")))

        csu_df <- one_df %>%
          filter(!grepl("virridy", name, ignore.case = TRUE))

        readr::write_csv(csu_df,
                         paste0(dump_dir, "/", site[i], "_", stringr::str_replace(stringr::str_replace(substr(end_dt, 1, 16), "[:\\s]", "_"), ":", ""), ".csv"))

      } else {

        # otherwise, save the full data set

        readr::write_csv(one_df,
                         paste0(dump_dir, "/", site[i], "_", stringr::str_replace(stringr::str_replace(substr(end_dt, 1, 16), "[:\\s]", "_"), ":", ""), ".csv"))
      }
    }
  }

}
append_historical_api_data <- function(hist_dir, inc_dir) {

  # list of site names
  site_names <- list(
    "Archery",
    "Boxelder",
    "Legacy",
    "Lincoln",
    "Prospect",
    "River Bluffs",
    "Tamasag",
    "Timberline"
  )

  # get the full file names
  inc_full_file_names <- list.files(inc_dir, pattern = "*.csv", full.names = TRUE)
  hist_full_file_names <- list.files(hist_dir, pattern = "*.csv", full.names = TRUE)

  # find the files that match using the site names
  walk(site_names, function(site_name) {

    # find the index of the matching site names in the file names
    inc_site_name_full_path <- grepl(site_name, inc_full_file_names, ignore.case = TRUE)
    hist_site_name_full_path <- grepl(site_name, hist_full_file_names, ignore.case = TRUE)

    # get the file names
    inc_site_path <- inc_full_file_names[inc_site_name_full_path]
    hist_site_path <- hist_full_file_names[hist_site_name_full_path]

    # read in the files
    inc_file <- read_csv(inc_site_path)
    hist_file <- read_csv(hist_site_path)

    # combine the files
    new_file <- bind_rows(hist_file, inc_file) %>%
      distinct()

    # write the file (this will replace the old file)
    write_csv(new_file, paste0(hist_dir, site_name, "_historical.csv"))
  })

  # delete the files in the inc_dir
  walk(inc_full_file_names, unlink)

}



basic_plot <- function(){
  max_1 <- max(wq_tl[[parameters[1]]], na.rm = T)

 param_1 <- ggplot(wq_tl, aes(x = DT_round, y = .data[[parameters[1]]], ymin = 0, ymax = .data[[parameters[1]]]))+
   geom_ribbon(color = "blue", fill = "blue" )+
   theme_bw()+
   labs(x = "Date", y = parameters[1])

 param_2 <- ggplot(wq_tl, aes(x = DT_round, y = .data[[parameters[2]]]))+
   geom_path(color = "red")+
   theme_bw()+
   labs(x = "Date", y = parameters[2])

 ggarrange(param_1, param_2,ncol = 2, nrow = 1)

}
check_incoming_api_dir <- function(incoming_dir, archive_dir) {
  # Check if data/api/incoming_api_data exists
  if (dir.exists(incoming_dir)) {
    # Check if incoming directory is empty
    if (length(list.files(incoming_dir)) == 0) {
      print(paste0(incoming_dir, " exists and is empty."))
      print("Incoming API data directory is configured properly.")
    } else {
      print(paste0(incoming_dir, " exists but is not empty..."))
      print("Please ensure previous workflow ran properly...")
      stop("Pipeline halted due to non-empty incoming directory.")
    }
  } else {
    print(paste0(incoming_dir, " does not exist..."))
    print("Creating incoming directory...")
    dir.create(incoming_dir, recursive = TRUE)
    print("Incoming directory created.")
  }

  # Check if data/api/archive_api_data exists
  if (dir.exists(archive_dir)) {
    print(paste0(archive_dir, " exists."))
    print("Directory is configured properly.")
  } else {
    print(paste0(archive_dir, " does not exist..."))
    print("Creating archive directory...")
    dir.create(archive_dir, recursive = TRUE)
    print("Archive directory created.")
  }
}
clean_directories <- function() {
  pre_dir_path <- pre_verification_path
  int_dir_path <- intermediary_path
  ver_dir_path <- verified_path

  pre_dir_names <- list.files(pre_dir_path)
  int_dir_names <- list.files(int_dir_path)
  ver_dir_names <- list.files(ver_dir_path)

  # Pre-verification Directory
  # keep if not in Verified Directory, else if in Verified Directory delete from this Directory
  for (i in pre_dir_names) {
    if (length(ver_dir_names) != 0 & i %in% ver_dir_names) {
      file.remove(paste0(pre_dir_path, i))
      cat("removed file ", i, " from ",  pre_dir_path, "\n")
    } else {
      next
    }
  }

  # Intermediary Directory
  # keep if skips in the data OR any(!is_verified), else move to Verified Directory
  for (i in int_dir_names) {
    df <- readRDS(paste0(int_dir_path, i))
    if (any(df$verification_status == 'SKIP') | any(!df$is_verified)) {
      next
    } else if (all(df$verification_status != 'SKIP') & all(df$is_verified)) {
      source_path <- paste0(int_dir_path, i)
      destination_path <- paste0(ver_dir_path, i)
      file.copy(source_path, destination_path, overwrite = TRUE)
      file.remove(source_path)
      cat("moved file ", i, " from ", int_dir_path, "to", ver_dir_path, "\n")
    }
  }

  # Verified Directory
  # final destination for data, check that the data that is here should be here
  for (i in ver_dir_names) {
    df <- readRDS(paste0(ver_dir_path, i))
    if (any(df$verification_status == 'SKIP') | any(!df$is_verified)) {
      cat(i, "should not be in verified directory\n")
      stop("THERE IS DATA WITH SKIPS OR NON-VERIFIED DATA IN THE VERIFIED DIRECTORY")
    } else {
      cat("data in verified directory is correct\n")
      next
    }
  }
}
clear_incoming_data_dir <- function(incoming_dir, archive_dir, require = NULL) {
  # Check if the previous step in the targets pipeline ran
  # this is done via the require arg

  # List files in the incoming and archive directories
  incoming_files <- list.files(incoming_dir, full.names = FALSE)
  archive_files <- list.files(archive_dir, full.names = FALSE)

  # Find files that are not already in the archive directory
  files_to_copy <- setdiff(incoming_files, archive_files)

  # Copy only the files that are not already in the archive directory
  if (length(files_to_copy) > 0) {
    for (file in files_to_copy) {
      full_file_name <- file.path(incoming_dir, file)
      file.copy(full_file_name, archive_dir)
      print(paste0(file, " has been moved to archive API data folder."))
    }
    print("Files have been copied from the incoming directory to the archive directory.")
  } else {
    print("All files are already present in the archive directory. Nothing to copy.")
  }
Sys.sleep(5)

# refresh archive_files now that the incoming_files have been copied over
archive_files <- list.files(archive_dir, full.names = FALSE)
  # Check if all files from incoming directory are now in the archive directory
  if (all(incoming_files %in% archive_files)) {
    print("All files in the incoming directory have been successfully copied to the archive directory.")
  } else {
    print("Not all files from the incoming directory have been successfully copied to the archive directory.")
    # Should this halt the pipeline?
    # Right now it seems to always print this and I haven't figured out why so
    # I don't think that it should not halt the pipeline in the state that its in right now -JD.
  }

  # Delete the copied files from the incoming directory
  if (length(files_to_copy) > 0) {
    for (file in files_to_copy) {
      full_file_name <- file.path(incoming_dir, file)
      file.remove(full_file_name)
    }
    print("Copied files have been removed from the incoming directory.")
  }

  # Delete any files in the incoming directory
  for (file in list.files(incoming_dir, full.names = TRUE)) {
    file.remove(file)
  }
  print("All files removed from incoming directory.")
}
#' @title Combine historical and incoming data
#'
#' @description
#' A function that combines the subset of historically flagged data and creates historical columns.
#' NOTE: historical columns should be a part of the historically flagged data, but it is not at this point.
#'
#' @param incoming_data_list A list of dataframes with the incoming data.
#'
#' @param historical_data_list A list of dataframes with the historical data.
#'
#' @return A list of dataframes with the combined historical and incoming data.
#'
#' @examples
#' combine_hist_inc_data(incoming_data_list = incoming_data_list, historical_data_list = historical_data_list)

combine_hist_inc_data <- function(incoming_data_list, historical_data_list) {

  # Get the matching index names
  matching_indexes <- intersect(names(incoming_data_list), names(historical_data_list))

  # Get the last 24 hours of the historically flagged data based on earliest
  # DT listed in the api pull for that site/data combo
  last_24_hours <- map(historical_data_list,
                      function(data) {
                        data %>%
                          # most recent day of already-flagged data
                          filter(DT_round >= ymd_hms(max(DT_round) - hours(24), tz = "MST")) %>%
                          # mark it as "historic"
                          mutate(historical = TRUE,
                                 # ensure last_site_visit column has proper DT:
                                 last_site_visit = force_tz(last_site_visit, tzone = "MST"))
                      })

  # bind summarized_incoming_data and last_24_hours together
  combined_hist_inc_data <- map(matching_indexes, function(index) {
    last_24_hours[[index]] %>%
      select(-c(t_mean01, # why is this being removed? -jd
                t_mean99,
                t_slope_behind_01,
                t_slope_behind_99,
                t_sd_0199)) %>%
      bind_rows(., mutate(incoming_data_list[[index]], historical = FALSE)[-1,]) %>% # why are we removing the first row? -jd # there is a datetime issue on this bind
      # make sure new data has info about last site visit, etc. filled out
      fill(c(sonde_employed, last_site_visit, sensor_malfunction))

  }) %>%
    set_names(matching_indexes)

  return(combined_hist_inc_data)
}


#@param folder_name: input the raw photo folder that is named the site you wish to compile
compile_files <- function(folder_name){
  site <- str_extract(folder_name, "tamasag|legacy|timberline|prospect|boxelder|archery|riverbluffs")
  photo_files <- list.files(path = folder_name, full.names = TRUE, recursive = TRUE)
  #files in folder
  new_folder_files <- list.files(path = paste0('data/timelapse_photos/2023_compiled/',site), full.names = TRUE)

  photo_renamer <- function(file) {
    #grab dt from file
    dt <- read_exif(path = file,
                    tags = c("DateTimeOriginal")) %>%
      pull(DateTimeOriginal) %>%
      parse_date_time("YmdHMS", tz="MST")%>%
      format("%Y%m%d_%H%M")
    #create new file name from dt
    new_file_name <- paste0('data/timelapse_photos/2023_compiled/',site,"/", dt,'.JPG')

    #check to see if this file is already in the folder
    if(new_file_name %nin% new_folder_files){
      #if it is not, copy it over
      file.copy(file,to =  new_file_name)
    }

  }
  map(photo_files, photo_renamer)
  print(paste0("Finished ", folder_name))

}
## Photos

# Goals:
# Download all user created photos ( upstream, downstream, clarity, filter and other pictures)
# Label according to site, date, description in the format site_YYYYMMDD_descriptor.jpg
# Only download photos which have not yet been downloaded




download_pictures <- function(){
  #source to grab all notes cleaned
  source("src/load_mWater_notes.R")

  all_notes_cleaned <- load_mWater_notes()
  # Find all the downloaded pictures
  all_file_names <- tolower(list.files(path = "data/field_pics/", recursive = TRUE))
  #grab notes
  sampling_photos <- all_notes_cleaned%>%
    #grab needed columns
    select(site, start_dt,photos_downloaded,upstream_pic,downstream_pic,clarity,filter_pic,other_pic,other_pic_descriptor)%>%
    mutate(
      #correct names if it is in our upper sites (upper case acronyms)
      site = tolower(site),
      #Date format for pictures
      yyyymmdd = format(start_dt, "%Y%m%d"),
      #create filenames ONLY if there is a URL associated with the site visit
      upstream_filename = case_when(
        !is.na(upstream_pic) ~ paste0(site, "_", yyyymmdd, "_upstream.jpg"),
        TRUE ~ NA_character_
      ),
      downstream_filename = case_when(
        !is.na(downstream_pic) ~ paste0(site, "_", yyyymmdd, "_downstream.jpg"),
        TRUE ~ NA_character_
      ),
      clarity_filename = case_when(
        !is.na(clarity) ~ paste0(site, "_", yyyymmdd, "_clarity.jpg"),
        TRUE ~ NA_character_
      ),
      filter_filename = case_when(
        !is.na(filter_pic) ~ paste0(site, "_", yyyymmdd, "_filter.jpg"),
        TRUE ~ NA_character_
      ),
      # check to see if the photo is already downloaded to folder
      upstream_downloaded = case_when(
        is.na(upstream_filename) ~ NA,
        upstream_filename %in% all_file_names ~ TRUE,
        TRUE ~ FALSE
      ),
      clarity_downloaded = case_when(
        is.na(clarity_filename) ~ NA,
        clarity_filename %in% all_file_names ~ TRUE,
        TRUE ~ FALSE
      ),
      filter_downloaded = case_when(
        is.na(filter_filename) ~ NA,
        filter_filename %in% all_file_names ~ TRUE,
        TRUE ~ FALSE
      ),
      downstream_downloaded = case_when(
        is.na(downstream_filename) ~ NA,
        downstream_filename %in% all_file_names ~ TRUE,
        TRUE ~ FALSE
      )
    )
# basic path to field pics
path <- "data/field_pics/"


  # loop thru dataset and download the photo ONLY if it is not yet downloaded and not NA
  for (i in 1:nrow(sampling_photos)) {
    if (!is.na(sampling_photos$upstream_downloaded[i]) && !sampling_photos$upstream_downloaded[i]) {
      print(sampling_photos$upstream_filename[i])
      download.file(sampling_photos$upstream_pic[i], destfile = paste0(path,sampling_photos$upstream_filename[i]))
    }

    if (!is.na(sampling_photos$downstream_downloaded[i]) && !sampling_photos$downstream_downloaded[i]) {
      print(sampling_photos$downstream_filename[i])
      download.file(sampling_photos$downstream_pic[i], destfile = paste0(path, sampling_photos$downstream_filename[i]))
    }
    if (!is.na(sampling_photos$clarity_downloaded[i]) && !sampling_photos$clarity_downloaded[i]) {
      print(sampling_photos$clarity_filename[i])
      download.file(sampling_photos$clarity[i], destfile = paste0(path, sampling_photos$clarity_filename[i]))
    }
    if (!is.na(sampling_photos$filter_downloaded[i]) && !sampling_photos$filter_downloaded[i]) {
      print(sampling_photos$filter_filename[i])
      download.file(sampling_photos$filter_pic[i], destfile = paste0(path, sampling_photos$filter_filename[i]))
    }
  }

  #grab notes for sites with other pictures
  other_photos <- all_notes_cleaned%>%
    #grab needed columns
    select(site, start_dt,other_pic,other_pic_descriptor)%>%
    #get rid of instances with no other pics
    filter(!is.na(other_pic))%>%
    mutate(
      site = tolower(site),
      #Date format for pictures
      yyyymmdd = format(start_dt, "%Y%m%d"),
      # separate multiple URLs in other pic column
      other_pic_sep = str_split(other_pic, "; "),
      #seperate multiple descriptors in the descriptor column
      other_descriptor_sep = str_split(other_pic_descriptor, ","))%>%
    #for rows with multiple pictures, create a new row for each picture

    unnest(cols = c(other_pic_sep, other_descriptor_sep))%>%
    #remove excess columns and rename sep columns to match old columns
    select(site, start_dt,yyyymmdd, other_pic = other_pic_sep, other_pic_descriptor = other_descriptor_sep)%>%
    # make descriptor lower case and remove any spaces in the name
    mutate(other_pic_descriptor = tolower(str_replace_all(other_pic_descriptor, " ", "")),
           other_filename = case_when(!is.na(other_pic) ~ paste0(site, "_", yyyymmdd, "_", other_pic_descriptor, ".jpg")),
           # Check to see if photo has already been downloaded
           other_downloaded = case_when(
             is.na(other_filename) ~ NA,
             other_filename %in% all_file_names ~ TRUE,
             TRUE ~ FALSE
           ))

  # loop thru dataset and download the photo ONLY if it is not yet downloaded and not NA
  for (i in 1:nrow(other_photos)) {
    if (!is.na(other_photos$other_downloaded[i]) && !other_photos$other_downloaded[i]) {
      print(other_photos$other_filename[i])
      download.file(other_photos$other_pic[i], destfile = paste0(path,other_photos$other_filename[i]))
    }}

cat("\nAll Available Pictures Downloaded\n")
}

#download_pictures()
#RUN TO DOWNLOAD NEW PICTURES
# It takes about 2-5 minutes to download ~25-50 photos
# Sometimes the request to mWater time out, just re run the function below if that happens
## Determining uploads

#This function looks at the user inputs for calibration report collect and logs collected.
#Based on these inputs, it looks at all the uploaded logs or calibration reports
#and will print out what logs are missing and who to contact to get those files uploaded.




files_missing <- function(){

  `%nin%` = Negate(`%in%`)
  # #source clean mwater script for all notes cleaned
  #
  # source("src/mWater_collate/load_mWater_notes.R")

  #grab context metadata
  site_meta <- read_csv("data/metadata/water_sampling_sites.csv",show_col_types = FALSE)%>%
    select(site = site_code, Site_Name, site_label_rmrs)
  # sort for sites in upper network (ie. acronyms rather than street names)
  upper_sites <- read_csv("data/metadata/water_sampling_sites.csv",show_col_types = FALSE)%>%
    filter(watershed != "CLP  Mainstem-Fort Collins")%>%
    #this is to help match with user input
    mutate(site_code = tolower(site_code))

  field_season <- year(Sys.Date())

  #grab cal reports from folder
  cal_reports_simple <- str_extract(list.files(path = "data/calibration_reports/"), ".*_\\d{8}" )%>%tolower()
  logs_simple <- str_extract(list.files(path = paste0("data/sensor_data/", field_season), recursive = TRUE), "\\w+_\\d{8}_(vulink|troll)")%>%tolower()




  #grab sensor notes that have logs or cal reports that should be  downloaded
  sensor_files <- all_notes_cleaned%>%
    filter(year(DT_round) == field_season)%>%
    filter(grepl("Sensor",visit_type, ignore.case = TRUE))%>%
    filter(cal_report_collected|log_downloaded)%>%
    select(site, crew, start_DT,end_dt, cal_report_collected, cals_performed, log_downloaded, log1_type,log1_mmdd,  log2_type, log2_mmdd)%>%
    mutate(
      #make all site names lower
      site = tolower(site),
      # Create basis for calibration report name
      # this will be used to check for calibration report in data files and then b
      cal_report_name = case_when(cal_report_collected == TRUE ~ paste0(site, "_", format(start_DT, "%Y%m%d")),
                                  cal_report_collected == NA ~ NA),
      full_cal_name = case_when(cal_report_collected == TRUE ~ paste0(site, "_", format(end_dt, "%Y%m%d_%H%M_mst")),
                                cal_report_collected == NA ~ NA),
      log1_mmdd = case_when(nchar(as.character(log1_mmdd)) == 3 ~ paste0("0",log1_mmdd),
                            TRUE ~ as.character(log1_mmdd)),
      log1_type = case_when( grepl("aquatroll", log1_type,ignore.case = TRUE) ~  "troll",
                             grepl("at", log1_type,ignore.case = TRUE) ~  "troll",
                             grepl("aqua troll", log1_type,ignore.case = TRUE) ~  "troll",
                             grepl("vulink", log1_type,ignore.case = TRUE) ~  "vulink",
                             grepl("vu link", log1_type,ignore.case = TRUE) ~  "vulink",
                             TRUE ~ tolower(log1_type)),
      log2_type = case_when( grepl("aquatroll", log2_type,ignore.case = TRUE) ~  "troll",
                             grepl("at", log2_type,ignore.case = TRUE) ~  "troll",
                             grepl("aqua troll", log2_type,ignore.case = TRUE) ~  "troll",
                             grepl("vulink", log2_type,ignore.case = TRUE) ~  "vulink",
                             grepl("vu link", log2_type,ignore.case = TRUE) ~  "vulink",
                             TRUE ~ tolower(log2_type)),
      log2_mmdd = case_when(nchar(as.character(log2_mmdd)) == 3 ~ paste0("0",log2_mmdd),
                            TRUE ~ as.character(log2_mmdd)),
      log1_user_error = case_when( is.na(log_downloaded)~ FALSE,
                                   log_downloaded == FALSE ~ FALSE,
                                   is.na(log1_mmdd) | is.na(log1_type) ~ TRUE,
                                   TRUE ~ FALSE ),
      log2_user_error = case_when(is.na(log_downloaded)~ FALSE,
                                  log_downloaded == FALSE ~ FALSE,
                                  is.na(log2_mmdd) & is.na(log2_type) ~ FALSE,
                                  is.na(log2_mmdd) | is.na(log2_type) ~ TRUE,
                                  TRUE ~ FALSE),
      log1_name = case_when(!(is.na(log1_mmdd) | is.na(log1_type)) ~ paste0(site,"_",year(start_DT),log1_mmdd,
                                                                            "_", format(start_DT, "%Y%m%d"), "_", log1_type),
                            (is.na(log1_mmdd) | is.na(log1_type))  ~ NA),
      log2_name = case_when(!(is.na(log2_mmdd) | is.na(log2_type)) ~ paste0(site,"_",year(start_DT),log2_mmdd,
                                                                            "_", format(start_DT, "%Y%m%d"), "_", log2_type),
                            (is.na(log2_mmdd) | is.na(log2_type))  ~ NA),
      log_missing1 = case_when(
        is.na(log_downloaded)| log_downloaded == FALSE ~ FALSE,
        is.na(log1_name) ~ FALSE,
        log1_name %nin% logs_simple ~ TRUE,
        TRUE ~ FALSE
      ),
      log_missing2 = case_when(
        is.na(log_downloaded)| log_downloaded == FALSE ~ FALSE,
        is.na(log2_name) ~ FALSE,
        log2_name %nin% logs_simple ~ TRUE,
        TRUE ~ FALSE
      ),
      log_missing = case_when(
        log_missing1 | log_missing2 ~ TRUE,
        TRUE ~ FALSE
      ),
      cal_missing = case_when(
        is.na(cal_report_collected) ~ FALSE,
        cal_report_name %nin% cal_reports_simple ~ TRUE,
        TRUE ~ FALSE)
      )



  for (i in 1:nrow(sensor_files)) {

    #if log missing print out missing logs or cal reports
    if(sensor_files$log_missing[i]){
      cat("\nLog Missing: ", sensor_files$log1_name[i], " and/or ", sensor_files$log2_name[i], "\nContact: ", sensor_files$crew[i], "\n")

    }
    #if log missing print out missing logs or cal reports
    if(sensor_files$cal_missing[i]){
      cat("\nCal Missing: ", sensor_files$full_cal_name[i]," \nContact: ", sensor_files$crew[i], "\n")
    }


  }
  cat("\nFile Check Complete")

}

#' @title Add malfunction flags to a data frame based on field notes and HydroVu.
#'
#' @description
#' This function adds the 'sensor malfunction' flag to a data frame by the dates,
#' sites, and parameters that members of the lab know contain erroneous data due
#' to sensor malfunctions. Note that this flag is used in two instances: first
#' when removing erroneous data from our statistics calculations and again during
#' the actual flagging step in the QAQC pipeline.
#'
#' @param df A data frame with a `flag` column.
#'
#' @param malfunction_records The malfunction records pulled from mWater
#'
#' @return A data frame with a `flag` column that has been updated with the
#' 'sensor malfunction' flag.
#'
#' @examples
#' add_range_flags(df = all_data_flagged$`archery-Actual Conductivity`)
#' add_range_flags(df = all_data_flagged$`boxelder-Temperature`)
#'
#' @seealso [grab_mWater_malfunction_notes()]

fix_calibration <- function(df){

  # Filter records for relevant site-param information
  df_site <- unique(df$site)
  df_parameter <- unique(df$parameter)

  # Depth calibration requires its own fix. It's hacky and I don't like it, but
  # the way that the calibration report reports pressure calibrations makes it
  # impossible to actually back calibrate. Therefore, I'm just assuming that level
  # doesn't actually change after bad calibrations, and hard code it so that the
  # first "bad" depth is forced to equal the last "good" depth and I offset all
  # subsequent "bad" depth values by the difference between them.

  if(df_parameter == "Depth"){

    if(!"Depth" %in% df$parameter){

      nope <- df %>% mutate(relative_depth = NA)

      return(nope)

    }

    if("Depth" %in% df$parameter & "archery" %in% df$site){

      #df <- all_data_flagged[["archery-Depth"]]

      site_depth <- df %>%
        mutate(relative_depth = ifelse(DT_round >= as_datetime('2022-05-21 15:45:00', "MST") & DT_round <= as_datetime('2022-05-24 15:45:00', "MST"),
                                       mean +
                                         abs(
                                           dplyr::filter(df, as.character(DT_round) == ('2022-05-24 15:15:00'))$mean -
                                             dplyr::filter(df, as.character(DT_round) == ('2022-05-24 15:45:00'))$mean
                                         ),
                                       mean))
      #return(archery_depth)

    } else if("Depth" %in% df$parameter & "timberline" %in% df$site){

      # df <- all_data_flagged[["timberline-Depth"]]

      site_depth <- df %>%
        dplyr::mutate(relative_depth = ifelse(year == "2022" & DT_round <= lubridate::as_datetime('2022-04-07 17:00:00', "MST"),
                                              as.numeric(mean) +
                                                abs(
                                                  dplyr::filter(df, as.character(DT_round) == ('2022-04-07 16:15:00'))$mean -
                                                    dplyr::filter(df, as.character(DT_round) == ('2022-04-07 17:15:00'))$mean
                                                ),
                                              as.numeric(mean)))
      # return(timberline_depth)

    } else if ("Depth" %in% df$parameter & "legacy" %in% df$site) {

      # df <- all_data_flagged[["legacy-Depth"]]

      site_depth <- df %>%
        mutate(relative_depth = ifelse(DT_round >= as_datetime('2022-04-06 06:00:00', "MST") & DT_round <= as_datetime('2022-04-12 09:15:00', "MST"),
                                       as.numeric(mean) +
                                         abs(
                                           dplyr::filter(df, as.character(DT_round) == ('2022-04-12 09:15:00'))$mean -
                                             dplyr::filter(df, as.character(DT_round) == ('2022-04-12 09:30:00'))$mean
                                         ),
                                       mean))

      site_depth <-  site_depth %>%
        mutate(relative_depth = ifelse(DT_round >= as_datetime('2022-07-08 17:00:00', "MST") & DT_round <= as_datetime('2022-07-12 09:00:00', "MST"),
                                       relative_depth +
                                         abs(
                                           dplyr::filter(site_depth, as.character(DT_round) == ('2022-07-08 14:15:00'))$relative_depth -
                                             dplyr::filter(site_depth, as.character(DT_round) == ('2022-07-08 17:00:00'))$relative_depth
                                         ),
                                       relative_depth))

      site_depth <-  site_depth %>%
        mutate(relative_depth = ifelse(DT_round >= as_datetime('2022-07-22 11:30:00', "MST") & DT_round <= as_datetime('2022-07-25 14:15:00', "MST"),
                                       relative_depth +
                                         abs(
                                           dplyr::filter(site_depth, as.character(DT_round) == ('2022-07-22 09:45:00'))$relative_depth -
                                             dplyr::filter(site_depth, as.character(DT_round) == ('2022-07-22 11:30:00'))$relative_depth
                                         ),
                                       relative_depth))



      #return(legacy_depth)

    } else if ("Depth" %in% df$parameter & "tamasag" %in% df$site) {

      # df <- all_data_flagged[["tamasag-Depth"]]

      site_depth <- df %>%
        mutate(relative_depth = ifelse(DT_round <= "2022-04-24 07:15:00" & year == "2022",
                                       as.numeric(mean) +
                                         abs(
                                           dplyr::filter(df, as.character(DT_round) == ('2022-04-24 07:15:00'))$mean -
                                             dplyr::filter(df, as.character(DT_round) == ('2022-04-24 07:30:00'))$mean
                                         ),
                                       mean))

      # return(tamasag_depth)

    } else {

      site_depth <- df %>%
        mutate(relative_depth = mean,
               cal_fix = NA)
    }

    depth_flagged <- site_depth %>%
      dplyr::mutate(raw = mean,
                    mean = relative_depth,
                    cal_fix = ifelse(raw != mean, "calibration fix", NA)) %>%
      dplyr::select(-relative_depth)

    return(depth_flagged)

  }

  # PLACEHOLDER UNTIL OTHER CALIBRATIONS ARE GOOD:
 return(df %>%
    mutate(raw = mean,
           cal_fix = NA))

}

  # # For non-depth parameters, we can refer to the calibration report:
  # bad_cal <- readxl::read_excel("data/calibration_error_log.xlsx") %>%
  #   dplyr::mutate(start_DT = as.character(lubridate::as_datetime(start_DT)),
  #                 end_DT = as.character(lubridate::as_datetime(end_DT)),
  #                 #report_to_correct = as.character(lubridate::as_datetime(report_to_correct))
  #   )
  #
  # bad_cal_records_filtered <- bad_cal %>%
  #   filter(site == df_site) %>%
  #   filter(grepl(df_parameter, parameter, ignore.case = TRUE)) %>%
  #   mutate(end_DT = ifelse(is.na(end_DT), ymd_hms("9999-12-31 23:59:59", tz = "MST"), end_DT)) %>%
  #   mutate(end_DT = as.character(as.POSIXct(end_DT, tz = "MST"))) %>%
  #   rowid_to_column()
  #
  # # If there are no bad calibrations listed for that site-param, return original
  # # dataframe, filling our updated "mean" column with the old unmodified values:
  # if(nrow(bad_cal_records_filtered == 0)){
  #
  #   df <- df %>%
  #     mutate(raw = mean,
  #            mean = mean)
  #
  #   return(df)
  #
  # } else {
  #
  #   cal_tabler <- function(cal_files){
  #
  #     #cal_files <- list.files("data/calibration_reports")[3] # for testing
  #
  #     cal <- read_html(paste0(getwd(), "/data/calibration_reports/", cal_files)) %>%
  #       html_nodes("div") %>%
  #       html_text() %>%
  #       as_tibble()
  #
  #     rdo <- cal %>% filter(grepl("RDO", value)) %>% pull() %>% str_replace_all(., " ", "") %>% tolower()
  #
  #     ph_orp <- cal %>% filter(grepl("pH/ORP", value)) %>% pull() %>% str_replace_all(., " ", "") %>% tolower()
  #
  #     conductivity <- cal %>% filter(grepl("Conductivity",value)) %>% pull() %>% str_replace_all(., " ", "") %>% tolower()
  #
  #     if(length(cal %>% filter(grepl("Turbidity",value)) %>% pull() %>% str_replace_all(., " ", "") %>% tolower()) != 0){
  #
  #       turbidity <- cal %>% filter(grepl("Turbidity",value)) %>% pull() %>% str_replace_all(., " ", "") %>% tolower()
  #
  #     } else {
  #
  #       turbidity <- "No Turbidity Sensor"
  #
  #     }
  #
  #     # Always the last sensor when depth is available:
  #     depth <- ifelse(str_detect(cal %>% .[nrow(.),] %>% pull() %>% str_replace_all(., " ", "") %>% tolower(), "pressure"),#"psireferencedepth"),
  #                     cal %>% .[nrow(.),] %>% pull() %>% str_replace_all(., " ", "") %>% tolower(),
  #                     "No Depth Sensor")
  #
  #     time_mst <- paste0(str_sub(str_match(cal_files, "(\\d+)_mst")[, 2:1][1], 1, 2), ":",
  #                        str_sub(str_match(cal_files, "(\\d+)_mst")[, 2:1][1], 3, 4))
  #     #str_sub(cal_files, -13, -12),":", str_sub(cal_files, -11, -10))
  #
  #     date <- str_match(cal_files, "^[^_]+_([0-9]{8})_")[, 2]
  #
  #     #str_sub(cal_files, -22, -19),"-", str_sub(cal_files, -18, -17),"-", str_sub(cal_files, -16, -15))
  #
  #     cal_table <- tibble(site = sub("\\_.*", "", cal_files) %>% tolower(),
  #
  #                         DT = ymd_hm(paste(date, time_mst, tz = "MST")),
  #
  #                         # Dissolved Oxygen
  #                         rdo_cal_date = as.character(mdy(str_match(rdo, "lastcalibrated\\s*(.*?)\\s*calibrationdetails")[,2])),
  #                         rdo_slope = str_match(rdo, "slope\\s*(.*?)\\s*offset")[,2],
  #                         rdo_offset = str_match(rdo, "offset\\s*(.*?)\\s*mg/l")[,2],
  #                         rdo_100 = str_match(rdo, "premeasurement\\s*(.*?)\\s*%satpost")[,2],
  #                         rdo_conc = str_match(rdo, "concentration\\s*(.*?)\\s*mg/lpremeasurement")[,2],
  #                         rdo_temp = str_match(rdo, "temperature\\s*(.*?)\\s*°c")[,2],
  #                         rdo_pressure = str_match(rdo, "pressure\\s*(.*?)\\s*mbar")[,2],
  #
  #                         # pH
  #                         ph_cal_date = as.character(mdy(str_match(ph_orp, "lastcalibrated\\s*(.*?)\\s*calibrationdetails")[,2])),
  #                         ph_slope_pre = str_match(ph_orp, "offset1slope\\s*(.*?)\\s*mv/ph")[,2],
  #                         ph_offset_pre = str_match(ph_orp, "mv/phoffset\\s*(.*?)\\s*mvslopeandoffset2")[,2],
  #                         ph_slope_post = str_match(ph_orp, "offset2slope\\s*(.*?)\\s*mv/ph")[,2],
  #                         ph_offset_post = str_match(ph_orp, paste0(ph_slope_post,"mv/phoffset\\s*(.*?)\\s*mvorporp"))[,2],
  #                         # Sometimes, the post value can actually be in the high 6 pH... therefore the post measurement regex matching text is conditional
  #                         ph_7_nice = str_sub(str_match(ph_orp, "postmeasurementph7\\s*(.*?)\\s*mvcal")[,2], 10, nchar(str_match(ph_orp, "postmeasurementph7\\s*(.*?)\\s*mvcal")[,2])),
  #                         ph_7_high = str_sub(str_match(ph_orp, "postmeasurementph8\\s*(.*?)\\s*mvcal")[,2], 10, nchar(str_match(ph_orp, "postmeasurementph8\\s*(.*?)\\s*mvcal")[,2])),
  #                         ph_7_low = str_sub(str_match(ph_orp, "postmeasurementph6\\s*(.*?)\\s*mvcal")[,2], 10, nchar(str_match(ph_orp, "postmeasurementph6\\s*(.*?)\\s*mvcal")[,2])),
  #                         ph_7 = ifelse(!is.na(ph_7_nice), ph_7_nice,
  #                                       ifelse(!is.na(ph_7_high), ph_7_high, ph_7_low)),
  #
  #                         # ORP
  #                         #Newly encountered thing: sometimes the calibration report calls the ORP standard Zobell's, sometimes it's just called "ORP Standard":
  #                         orp_offset = ifelse(is.na(str_match(ph_orp, "zobell'soffset\\s*(.*?)\\s*mvtemperature")[,2]) & is.na(str_match(ph_orp, "quickcal\\s*(.*?)\\s*mvtemperature")[,2]),
  #                                             str_match(ph_orp, "orpstandardoffset\\s*(.*?)\\s*mvtemperature")[,2],
  #                                             ifelse(is.na(str_match(ph_orp, "zobell'soffset\\s*(.*?)\\s*mvtemperature")[,2]) & is.na(str_match(ph_orp, "orpstandardoffset\\s*(.*?)\\s*mvtemperature")[,2]),
  #                                                    str_match(ph_orp, "quickcal\\s*(.*?)\\s*mvtemperature")[,2],
  #                                                    str_match(ph_orp, "zobell'soffset\\s*(.*?)\\s*mvtemperature")[,2])),
  #
  #                         # Conductivity
  #                         cond_cal_date = as.character(mdy(str_match(conductivity, "lastcalibrated\\s*(.*?)\\s*calibrationdetails")[,2])),
  #                         tds_conversion_ppm = str_sub(str_match(conductivity, "tdsconversionfactor\\s*(.*?)\\s*cellconstant")[,2], 6, nchar(str_match(conductivity, "tdsconversionfactor\\s*(.*?)\\s*cellconstant")[,2])),
  #
  #                         # calibration report formatting has changed in 2024 for this variable. Therefore a post-2024 correction must occur
  #                         cond_cell_constant = ifelse(year(DT) < 2024, str_match(conductivity, "cellconstant\\s*(.*?)\\s*referencetemperature")[,2],
  #                                                     str_match(conductivity, "cellconstant\\s*(.*?)\\s*offset")[,2]),
  #
  #                         cond_offset = ifelse(year(DT) < 2024, NA,
  #                                              str_match(conductivity, "offset\\s*(.*?)\\s*µs/cm")[,2]),
  #
  #                         cond_pre = str_match(conductivity,paste0(str_match(conductivity,
  #                                                                            "premeasurementactual\\s*(.*?)\\s*specificconductivity")[,2],"specificconductivity\\s*(.*?)\\s*µs/cmpost"))[,2],
  #                         cond_post = str_match(conductivity,paste0(str_match(conductivity,
  #                                                                             "postmeasurementactual\\s*(.*?)\\s*specificconductivity")[,2],"specificconductivity\\s*(.*?)\\s*µs/cm"))[,2]) %>%
  #                         # if(turbidity == "No Turbidity Sensor"){
  #                         # # Turbidity
  #                         # turb_cal_date = "None",
  #                         # ntu_slope = "None",
  #                         # ntu_offset = "None",
  #                         # ntu_10 = "None",
  #                         # ntu_100 = "None") %>%
  #
  #     select(-c(ph_7_nice, ph_7_high, ph_7_low))
  #
  #     # Not all sondes have depth.
  #     if(!is.na(str_match(depth, "lastcalibrated"))){#\\s*(.*?)\\s*calibrationdetails")[,2])){
  #       cal_table <- cal_table %>%
  #         mutate(
  #           # Depth
  #           depth_cal_date = as.character(mdy(str_match(depth, "lastcalibrated\\s*(.*?)\\s*calibrationdetails")[,2])),
  #           depth_offset = str_match(depth, "zerooffset\\s*(.*?)\\s*psireferencedepth")[,2],
  #           depth_ref_depth = str_match(depth, "psireferencedepth\\s*(.*?)\\s*ftreferenceoffset")[,2],
  #           depth_ref_offset = str_match(depth, "ftreferenceoffset\\s*(.*?)\\s*psipremeasurement")[,2],
  #           depth_pre_psi = str_match(depth, "psipremeasurement\\s*(.*?)\\s*psipostmeasurement")[,2],
  #           depth_post_psi = str_match(depth, "psipostmeasurement\\s*(.*?)\\s*psi")[,2])
  #     }
  #
  #     if(depth == "No Depth Sensor"){
  #
  #       cal_table <- cal_table %>%
  #         mutate(# Depth
  #       depth_cal_date = "No Depth Sensor",
  #       depth_offset = "No Depth Sensor",
  #       depth_ref_depth = "No Depth Sensor",
  #       depth_ref_offset = "No Depth Sensor",
  #       depth_pre_psi = "No Depth Sensor",
  #       depth_post_psi = "No Depth Sensor")
  #     }
  #
  #
  #     if(!is.na(str_match(turbidity, "lastcalibrated"))){#calibrationpoint1premeasurement\\s*(.*?)\\s*ntupost")[,2])){
  #       # Not all sondes have turbidity.
  #       cal_table <- cal_table %>%
  #         mutate(
  #           # Turbidity
  #           turb_cal_date = as.character(mdy(str_match(turbidity, "lastcalibrated\\s*(.*?)\\s*calibrationdetails")[,2])),
  #           ntu_slope = str_match(turbidity, "slope\\s*(.*?)\\s*offset")[,2],
  #           ntu_offset = str_match(turbidity, "offset\\s*(.*?)\\s*ntu")[,2],
  #           ntu_10 = str_match(turbidity, "calibrationpoint1premeasurement\\s*(.*?)\\s*ntupost")[,2],
  #           ntu_100 = str_match(turbidity, "calibrationpoint2premeasurement\\s*(.*?)\\s*ntupost")[,2])
  #     }
  #
  #     if(turbidity == "No Turbidity Sensor"){
  #       cal_table <- cal_table %>%
  #         mutate(
  #           # Turbidity
  #           turb_cal_date = "No Turbidity Sensor",
  #           ntu_slope = "No Turbidity Sensor",
  #           ntu_offset = "No Turbidity Sensor",
  #           ntu_10 = "No Turbidity Sensor",
  #           ntu_100 = "No Turbidity Sensor")
  #
  #
  #
  #
  #     }
  #
  #
  #     cal_table <- cal_table %>%
  #       mutate(
  #         #Factory Defaults
  #         factory_defaults = paste0(ifelse(grepl("factorydefault", turbidity), "Turbidity ", ""),
  #                                   ifelse(grepl("factorydefault", rdo), "RDO ", ""),
  #                                   ifelse(is.na(ph_slope_post), "pH ", ""),
  #                                   ifelse(is.na(orp_offset), "ORP ", ""),
  #                                   ifelse(grepl("factorydefault", conductivity), "Conductivity ", ""),
  #                                   ifelse(grepl("factorydefaults", depth), "Depth ", ""))) %>%
  #       # convert all columns to character values to preserve info
  #       mutate(across(.cols = everything(), .fns = as.character)) %>%
  #       # remove "," from big numbers
  #       mutate(across(everything(), ~str_replace_all(., ",", "")))
  #
  #   }
  #
  #   bad_cal_interval_list <- map2(
  #     .x = bad_cal_records_filtered$start_DT,
  #     .y = bad_cal_records_filtered$end_DT,
  #     .f = ~interval(.x, .y, tz = "MST"))
  #
  #   if(df_parameter == "DO"){
  #
  #     cal_table <- list.files("data/calibration_reports", pattern=".html") %>%
  #       .[grepl(df_site, ., ignore.case = TRUE)] %>%
  #       map_dfr(., cal_tabler) %>%
  #       distinct(.keep_all = TRUE) %>%
  #       mutate(DT = as.character(round_date(ymd_hms(DT, tz = "MST"), "15 minutes"))) %>%
  #       # mutate(across(-matches("date|site|DT|factory"), as.numeric)) %>%
  #       dplyr::select(DT, site, rdo_slope, rdo_offset)
  #
  #     df_mod <- df %>%
  #       left_join(., cal_table, by = c("DT_join" = "DT", "site")) %>%
  #       fill(names(cal_table)[!grepl("\\b(site|DT)\\b", names(cal_table))], .direction = "down") %>%
  #       mutate(rdo_slope_pre = as.numeric(rdo_slope), rdo_offset_pre = as.numeric(rdo_offset)) %>%
  #       select(names(df), contains(c("pre"))) %>%
  #       left_join(., cal_table, by = c("DT_join" = "DT", "site")) %>%
  #       fill(names(cal_table)[!grepl("\\b(site|DT)\\b", names(cal_table))], .direction = "up") %>%
  #       mutate(rdo_slope_post = as.numeric(rdo_slope), rdo_offset_post = as.numeric(rdo_offset)) %>%
  #       select(names(df), contains(c("pre", "post"))) %>%
  #       #mutate(raw = (mean -rdo_offset_pre)/rdo_offset_pre)
  #       mutate(raw = case_when(DT_round %within% bad_cal_interval_list & is.na(rdo_slope_pre) ~ mean,
  #                              DT_round %within% bad_cal_interval_list & !is.na(rdo_slope_pre) ~ ((mean - rdo_offset_pre)/rdo_slope_pre),
  #                              .default = mean),
  #              cal_fix = case_when(DT_round %within% bad_cal_interval_list ~ (raw*rdo_slope_post) + rdo_offset_post,
  #                                  .default = mean)) %>%
  #       add_flag(mean != cal_fix, "calibration fix") %>%
  #       mutate(raw = mean,
  #              mean = cal_fixed)
  #
  #
  #   }
  #
  #   if(df_parameter == "pH"){
  #
  #     cal_table <- list.files("data/calibration_reports", pattern=".html") %>%
  #       .[grepl(df_site, ., ignore.case = TRUE)] %>%
  #       map_dfr(., cal_tabler) %>%
  #       distinct(.keep_all = TRUE) %>%
  #       mutate(DT = as.character(round_date(ymd_hms(DT, tz = "MST"), "15 minutes"))) %>%
  #       # mutate(across(-matches("date|site|DT|factory"), as.numeric)) %>%
  #       dplyr::select(DT, site, ph_slope_pre, ph_offset_pre, ph_slope_post, ph_offset_post, factory_defaults)
  #
  #     df_mod <- df %>%
  #       left_join(., cal_table, by = c("DT_join" = "DT", "site")) %>%
  #       fill(names(cal_table)[!grepl("\\b(site|DT)\\b", names(cal_table))], .direction = "down") %>%
  #       mutate(ph_slope_pre = as.numeric(ph_slope_pre), ph_offset_pre = as.numeric(ph_offset_pre),
  #              ph_slope_post = as.numeric(ph_slope_post), ph_offset_post = as.numeric(ph_offset_post)) %>%
  #       mutate(raw = case_when(DT_round %within% bad_cal_interval_list & is.na(ph_slope_pre) & grepl("pH", factory_defaults, ignore.case = TRUE) ~ mean,
  #                              DT_round %within% bad_cal_interval_list & !is.na(ph_slope_pre) & !grepl("pH", factory_defaults, ignore.case = TRUE) ~ ((mean - ph_offset_pre)/ph_slope_pre),
  #                              .default = mean),
  #              cal_fix = case_when(DT_round %within% bad_cal_interval_list ~ (raw*ph_slope_post) + ph_offset_post,
  #                                  .default = mean)) %>%
  #       add_flag(mean != cal_fix, "calibration fix")
  #
  #   }
  #
  #   if(df_parameter == "")
  #
  # }
# }
fix_depth_cal <- function(df){


  if(!"Depth" %in% df$parameter){

    nope <- df %>% mutate(relative_depth = NA)

   return(nope)

  }

  if("Depth" %in% df$parameter & "archery" %in% df$site){

     #df <- all_data_flagged[["archery-Depth"]]

    site_depth <- df %>%
      mutate(relative_depth = ifelse(DT_round >= as_datetime('2022-05-21 15:45:00', "MST") & DT_round <= as_datetime('2022-05-24 15:45:00', "MST"),
                                     mean +
                                       abs(
                                         dplyr::filter(df, as.character(DT_round) == ('2022-05-24 15:15:00'))$mean -
                                           dplyr::filter(df, as.character(DT_round) == ('2022-05-24 15:45:00'))$mean
                                       ),
                                     mean))
    #return(archery_depth)

  } else if("Depth" %in% df$parameter & "timberline" %in% df$site){

    # df <- all_data_flagged[["timberline-Depth"]]

    site_depth <- df %>%
      dplyr::mutate(relative_depth = ifelse(year == "2022" & DT_round <= lubridate::as_datetime('2022-04-07 17:00:00', "MST"),
                                            as.numeric(mean) +
                                              abs(
                                                dplyr::filter(df, as.character(DT_round) == ('2022-04-07 16:15:00'))$mean -
                                                  dplyr::filter(df, as.character(DT_round) == ('2022-04-07 17:15:00'))$mean
                                              ),
                                            as.numeric(mean)))
    # return(timberline_depth)

   } else if ("Depth" %in% df$parameter & "legacy" %in% df$site) {

    # df <- all_data_flagged[["legacy-Depth"]]

    site_depth <- df %>%
      mutate(relative_depth = ifelse(DT_round >= as_datetime('2022-04-06 06:00:00', "MST") & DT_round <= as_datetime('2022-04-12 09:15:00', "MST"),
                                     as.numeric(mean) +
                                       abs(
                                         dplyr::filter(df, as.character(DT_round) == ('2022-04-12 09:15:00'))$mean -
                                           dplyr::filter(df, as.character(DT_round) == ('2022-04-12 09:30:00'))$mean
                                       ),
                                     mean))

    site_depth <-  site_depth %>%
      mutate(relative_depth = ifelse(DT_round >= as_datetime('2022-07-08 17:00:00', "MST") & DT_round <= as_datetime('2022-07-12 09:00:00', "MST"),
                                     relative_depth +
                                       abs(
                                         dplyr::filter(site_depth, as.character(DT_round) == ('2022-07-08 14:15:00'))$relative_depth -
                                           dplyr::filter(site_depth, as.character(DT_round) == ('2022-07-08 17:00:00'))$relative_depth
                                       ),
                                     relative_depth))

    site_depth <-  site_depth %>%
      mutate(relative_depth = ifelse(DT_round >= as_datetime('2022-07-22 11:30:00', "MST") & DT_round <= as_datetime('2022-07-25 14:15:00', "MST"),
                                     relative_depth +
                                       abs(
                                         dplyr::filter(site_depth, as.character(DT_round) == ('2022-07-22 09:45:00'))$relative_depth -
                                           dplyr::filter(site_depth, as.character(DT_round) == ('2022-07-22 11:30:00'))$relative_depth
                                       ),
                                     relative_depth))

       #return(legacy_depth)

  } else if ("Depth" %in% df$parameter & "tamasag" %in% df$site) {

    # df <- all_data_flagged[["tamasag-Depth"]]

    site_depth <- df %>%
      mutate(relative_depth = ifelse(DT_round <= "2022-04-24 07:15:00" & year == "2022",
                                     as.numeric(mean) +
                                       abs(
                                         dplyr::filter(df, as.character(DT_round) == ('2022-04-24 07:15:00'))$mean -
                                           dplyr::filter(df, as.character(DT_round) == ('2022-04-24 07:30:00'))$mean
                                       ),
                                     mean))

    # return(tamasag_depth)

  } else if("Depth" %in% df$parameter) {

    site_depth <- df %>%
             mutate(relative_depth = mean)
}

  depth_flagged <- site_depth %>%
    add_flag(., relative_depth != mean, "post-calibration") %>%
    dplyr::mutate(mean = relative_depth)

  return(depth_flagged)

}
#' @title Apply all flags to a data frame
#'
#' @description
#' A function that applies all flags to a data frame. This function is used to
#' apply all flags to a data frame in one step.
#'
#' @param data A data frame with a `flag` column.
#'
#' @param require An upstream dependency (i.e., alternative to `{necessary object}`).
#'
#' @return A data frame with a `flag` column that has been updated with all flags.
#'
#' @examples
#' flag_all_data(data = all_data_flagged$`archery-Actual Conductivity`)
#' flag_all_data(data = all_data_flagged$`boxelder-Temperature`)
#'
#' @seealso [add_field_flag()]
#' @seealso [add_spec_flag()]
#' @seealso [add_seasonal_flag()]
#' @seealso [add_na_flag()]
#' @seealso [add_repeat_flag()]
#' @seealso [add_suspect_flag()]
#' @seealso [add_malfunction_flag()]

flag_all_data <- function(data, require = NULL) {
  flagged_data <- data %>%
    add_field_flag() %>%
    add_spec_flag() %>%
    add_seasonal_flag() %>%
    add_na_flag() %>%
    add_repeat_flag() %>%
    add_suspect_flag() %>%
    add_malfunction_flag(malfunction_records = mWater_malfunction_records)
    # mutate(mean_public = ifelse(is.na(flag), mean, NA))
  return(flagged_data)
}
generate_daily_plot <- function(plot_data_arg, df_list_arg, site_arg, parameter_arg) {

  # site_param_df <- df_list_arg[[paste0(site_arg, "-", parameter_arg)]]

  start_date <- min(plot_data_arg$DT_round)
  end_date <- max(plot_data_arg$DT_round)

  #default is lower network
  site_vector <- c("tamasag", # rist
                   "legacy",
                   "lincoln",
                   "timberline",
                   "prospect",
                   "boxelder", # elc
                   "archery",
                   "river bluffs")

  if(network == "virridy"){ # this will be the new default
    # establish order for all the non-tributary sites
    sites_order <-  c("joei","cbri","chd","pfal","sfm","pbd","tamasag",
                      "legacy","lincoln","timberline","prospect","boxelder",
                      "archery","riverbluffs")
    # establish the order for the tributary sites
    trib_sites_order <- c("boxcreek", "archery", NA, "springcreek", "prospect",
                          NA, "penn", "sfm", "lbea")
  }

  # determining the index for the site of interest.
  if (site_arg %in% sites_order) {

    plot_filter <- tibble(site = c("joei","cbri","chd","pfal","sfm","pbd",
                                   "tamasag","legacy", "lincoln","timberline",
                                   "timberline virridy","prospect",
                                   "prospect virridy","boxelder","archery",
                                   "archery virridy","riverbluffs"))

    site_index <- which(sites_order == site_arg)
    site_list <- as.vector(na.omit(sites_order[max(1, site_index - 1):min(length(sites_order), site_index + 1)]))

    plot_filter <- plot_filter %>%
      filter(grepl(paste(site_list, collapse = "|"), site, ignore.case = TRUE),
             site != site_arg) %>%
      pull(site)

  } else {

    plot_filter <- tibble(site = c("boxcreek", "archery", "archery virridy",
                                   "springcreek", "prospect", "prospect virridy",
                                   "penn", "sfm", "lbea"))

    site_index <- which(trib_sites_order == site_arg)
    site_list <- as.vector(na.omit(trib_sites_order[max(1, site_index - 1):min(length(trib_sites_order), site_index + 1)]))

    plot_filter <- plot_filter %>%
      filter(grepl(paste(site_list, collapse = "|"), site, ignore.case = TRUE),
             site != site_arg) %>%
      pull(site)

  }

  # Get the relevant sonde data
  # relevant_sondes <- map(plot_filter,
  #                        ~ {
  #                          sonde_name <- paste0(.x,"-",parameter_arg)
  #                          tryCatch({
  #                            sonde_df <- df_list_arg[[sonde_name]]  %>%
  #                              filter(DT_round %within% interval(start_date, end_date))},
  #                            error = function(err) {
  #                              cat("Sonde ", sonde_name," not found.\n")})
  #                        })

  # get the relevant sonde data source
  relevant_sonde_source <- map(plot_filter, ~ {
    sonde_name <- paste0(.x, "-", parameter_arg)
    # Determine which directory to pull data from
    tryCatch({
      retrieve_relevant_data_name(sonde_name, interval_arg = interval(start_date, end_date))
    }, error = function(err) {
      return("all_data")
    })
  })

  # Get the relevant data
  relevant_sondes <- map2(plot_filter,
                          relevant_sonde_source,
                          function(name, source) {
                            sonde_name <- paste0(name, "-", parameter_arg)
                            # try to pull in the data
                            tryCatch({
                              get(source)[[sonde_name]] %>%
                                filter(DT_round %within% interval(start_date, end_date))
                            }, error = function(err) {
                              return(NULL)
                            })
                          })

  # combine the lists
  sonde_info <- map2(relevant_sondes, relevant_sonde_source, list)

  # Remove any NULL results from the list
  sonde_info <- keep(sonde_info, ~!is.null(.x[[1]]))

  # append site_df to relevant sonde list, clean list, and bind dfs
  # to find plot info
  relevant_dfs <- map(sonde_info, ~.x[[1]])

  # append plot_data_arg to relevant sonde list, clean list, and bind dfs
  daily_plot_data <- append(relevant_dfs, list(plot_data_arg)) %>%
    keep(~ !is.null(.)) %>%
    keep(~ nrow(.)>0) %>%
    bind_rows()

  # use the daily flag data day as flag_day
  flag_day <- min(plot_data_arg$DT_round)

  plot <- ggplot(data = daily_plot_data) +
    geom_point(data = filter(daily_plot_data, (site == site_arg)),
               aes(x=DT_round, y=mean, color=flag)) +
    # geom_line(data = filter(daily_plot_data, (site != site_arg)),
    #           aes(x=DT_round, y=mean, color=site)) +
    map(sonde_info, function(sonde_data) {
      data <- sonde_data[[1]]
      data_source <- sonde_data[[2]]

      y_column <- if (data_source == "all_data") "mean" else "mean_verified"

      geom_line(data = data, aes(x = DT_round, y = .data[[y_column]], color = site))
    }) +
    ggtitle(paste0(str_to_title(site_arg), " ", str_to_title(parameter_arg), " (", format(flag_day, "%B %d, %Y"), ")")) +
    labs(x = "Time",
         y = "Mean")

  plot <- add_threshold_lines(plot = plot,
                              plot_data = plot_data_arg,
                              site_arg = site_arg,
                              parameter_arg = parameter_arg)

  plot <- plot +
    geom_vline(xintercept = seq(min(plot_data_arg$DT_round), max(plot_data_arg$DT_round), by = "hour"),
               color = "lightgrey", size = 0.5) +
    scale_x_datetime(breaks = seq(min(plot_data_arg$DT_round), max(plot_data_arg$DT_round), by = "hour"),
                     labels = format(seq(min(plot_data_arg$DT_round), max(plot_data_arg$DT_round), by = "hour"), "%H")) +
    theme_bw() +
    theme(legend.position = 'right',
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1)) + # make the background white
    guides(color = guide_legend(nrow = 10, byrow = TRUE))

  return(plot)

}
# Generate histograms of a site parameter data frame based on the mean column.
# @param df A dataframe with a `mean` column.
# @param df_index The index of the dataframe.
# @return A histogram of the mean column.
# @examples
# generate_general_histogram(df = all_data_flagged$`archery-Actual Conductivity`, df_index = "archery-Actual Conductivity")
# generate_general_histogram(df = all_data_flagged$`boxelder-Temperature`, df_index = "boxelder-Temperature")
generate_general_histogram <- function(df, df_index) {

  site_param <- toupper(sub("-", " ", df_index, fixed = TRUE))
  n <- as.character(sum(!is.na(df$mean)))
  title <- paste0("Histogram of ", site_param, " (n = ", n, ")")

  # there should be checks for data here
  minimum <- floor(min(df$mean, na.rm = TRUE))
  maximum <- ceiling(max(df$mean, na.rm = TRUE))

  histogram <- ggplot(data = df, aes(x = mean)) +
    geom_histogram(
      breaks = seq(minimum, maximum, by = 1)
    ) +
    labs(title = title)

  return(histogram)
}
generate_initial_weekly_plots <- function(all_df_list, pending_df_list, site_arg, parameter_arg, flag_arg = NULL) {

  site_param <- paste0(site_arg, "-", parameter_arg)

  site_flag_dates <- pending_df_list[[site_param]] %>%
    group_by(y_w) %>%
    filter(any(verification_status == "SKIP") | any(!is_verified)) %>%
    ungroup()

  if (!is.null(site_flag_dates)){
    # vector of sites in the order that they are in spatially ----
    # some sites have some funkiness going on (not all of the sites are present in the final plot)
    #default is lower network
    sites_order <- c("tamasag", # rist
                     "legacy",
                     "lincoln",
                     "timberline",
                     "prospect",
                     "boxelder", # elc
                     "archery",
                     "river bluffs")

    if(network == "virridy"){ # this will be the new default
      # establish order for all the non-tributary sites
      sites_order <-  c("joei","cbri","chd","pfal","sfm","pbd","tamasag",
                        "legacy","lincoln","timberline","prospect","boxelder",
                        "archery","riverbluffs")
      # establish order for all the tributary sites
      trib_sites_order <- c("boxcreek", "archery", NA, "springcreek", "prospect",
                            NA, "penn", "sfm", "lbea")
    }

    # determining the sites relevant to the site of interest.
    if (site_arg %in% sites_order) {

      plot_filter <- tibble(site = c("joei","cbri","chd","pfal","sfm","pbd",
                                     "tamasag","legacy", "lincoln","timberline",
                                     "timberline virridy","prospect",
                                     "prospect virridy","boxelder","archery",
                                     "archery virridy","riverbluffs"))

      site_index <- which(sites_order == site_arg)
      site_list <- as.vector(na.omit(sites_order[max(1, site_index - 1):min(length(sites_order), site_index + 1)]))

      plot_filter <- plot_filter %>%
        filter(grepl(paste(site_list, collapse = "|"), site, ignore.case = TRUE),
               site != site_arg) %>%
        pull(site)

    } else {

      plot_filter <- tibble(site = c("boxcreek", "archery", "archery virridy",
                                     "springcreek", "prospect", "prospect virridy",
                                     "penn", "sfm", "lbea"))

      site_index <- which(trib_sites_order == site_arg)
      site_list <- as.vector(na.omit(trib_sites_order[max(1, site_index - 1):min(length(trib_sites_order), site_index + 1)]))

      plot_filter <- plot_filter %>%
        filter(grepl(paste(site_list, collapse = "|"), site, ignore.case = TRUE),
               site != site_arg) %>%
        pull(site)

    }

    if (nrow(site_flag_dates >0)) {

      if (is.null(flag_arg)) {
        # This for loop generates an overlayed plot of weekly data for the site of
        # interest sandwiched by the site above and below it for each day that was
        # tagged by a flag of interest
        plot_list <- list()

        grouped_data <- site_flag_dates %>%
          group_by(y_w) %>% #group_by(week, year) %>% # group_by(week, month, year) %>%
          group_split()

        for(i in 1:length(grouped_data)) {

          group_data <- grouped_data[[i]]

          year_week <- unique(group_data$y_w)

          # filtering dfs of interest for the week of interest
          site_df <- site_flag_dates %>%
            filter(y_w == year_week)

          # Get the relevant sonde data
          relevant_sondes <- map(plot_filter, ~ {
            sonde_name <- paste0(.x, "-", parameter_arg)
            data_source <- NULL
            sonde_df <- NULL

            # Determine which directory to pull data from
            tryCatch({
              data_source <- retrieve_relevant_data_name(sonde_name, year_week)
              # cat("Data for",sonde_name,"will be pulled from",data_source,"\n")
            }, error = function(err) {
              # cat("Data for",sonde_name,"not found.\n")
              return(NULL)  # Return NULL if data source can't be determined
            })

            # Only try to pull in the data if data_source was successfully determined
            if (!is.null(data_source)) {
              tryCatch({
                sonde_df <- get(data_source)[[sonde_name]] %>%
                  filter(y_w == group_data$y_w)
              }, error = function(err) {
                cat("Sonde", sonde_name, "not found.\n")
                return(NULL)  # Return NULL if sonde data can't be retrieved
              })
            }

            # Only return a list if both data_source and sonde_df are available
            if (!is.null(data_source) & !is.null(sonde_df)) {
              return(list(sonde_df = sonde_df, data_source = data_source))
            } else {
              return(NULL)  # Return NULL if either data_source or sonde_df is NULL
            }
          })

          # Remove any NULL results from the list
          relevant_sondes <- compact(relevant_sondes)

          # append site_df to relevant sonde list, clean list, and bind dfs
          # to find plot info
          relevant_dfs <- map(relevant_sondes, ~.x[[1]])
          week_plot_data <- append(relevant_dfs, list(site_df)) %>% # how to relevant sondes here
            keep(~ !is.null(.)) %>%
            keep(~ nrow(.)>0) %>%
            bind_rows() %>%
            arrange(day)

          # Create a sequence of dates for the vertical lines
          start_date <- floor_date(min(week_plot_data$DT_round), "day")
          end_date <- ceiling_date(max(week_plot_data$DT_round), "day")
          vline_dates <- seq(start_date, end_date, by = "day")

          date_seq <- seq(min(vline_dates), max(vline_dates) - days(1), by = "1 day") + hours(12)

          # Use the first day of the group as flag_day
          flag_day <- min(group_data$DT_round)

          week_plot <- ggplot(data = week_plot_data) +
            geom_point(data = filter(week_plot_data, (site == site_arg)),
                       aes(x=DT_round, y=mean, color=flag)) +
            map(relevant_sondes, function(sonde_data) {
              data <- sonde_data[[1]]
              data_source <- sonde_data[[2]]

              y_column <- if (data_source == "all_data") "mean" else "mean_verified"

              geom_line(data = data, aes(x = DT_round, y = .data[[y_column]], color = site))
            })+
            geom_vline(xintercept = vline_dates, color = "black") +
            ggtitle(paste0(str_to_title(site_arg), " ", parameter_arg, " (", format(flag_day, "%B %d, %Y"), ")")) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            labs(x = "Day",
                 y = "Mean")

          week_plot <- add_threshold_lines(plot = week_plot,
                                           plot_data = week_plot_data,
                                           site_arg = site_arg,
                                           parameter_arg = parameter_arg)

          week_plot <- week_plot +
            theme_bw() +
            scale_x_datetime(date_breaks = "1 day",
                             date_labels = "%b %d",
                             minor_breaks = date_seq,
                             sec.axis = sec_axis(~., breaks = date_seq, labels = unique(week_plot_data$weekday))) +
            theme(legend.position = 'bottom',
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()) +
            guides(color = guide_legend(nrow = 4, byrow = TRUE))

          plot_list[[paste(site_param, as.character(flag_day))]] <- week_plot
          sorted_plot_names <- names(plot_list)[order(names(plot_list))]

          plot_list <- plot_list[sorted_plot_names]

        }
      }
      # ---- if flag arg is not null ----
      if (!is.null(flag_arg)) {
        print("under construction...")
        # # This for loop generates an overlayed plot of weekly data for the site of
        # # interest sandwiched by the site above and below it for each day that was
        # # tagged by a flag of interest
        # plot_list <- list()
        #
        # site_flag_weeks <- site_flag_dates %>%
        #   filter(str_detect(flag, flag_arg)) %>%
        #   group_by(week, year) %>%
        #   slice(1)
        #
        # grouped_data <- site_flag_dates %>%
        #   filter(y_w %in% site_flag_weeks$y_w) %>%
        #   group_by(week, year) %>% # group_by(week, month, year) %>%
        #   group_split()
        #
        #
        # for(i in 1:length(grouped_data)) {
        #
        #   group_data <- grouped_data[[i]]
        #
        #   # flag_title <- site_flag_dates$flag[i] # no flag title ***
        #
        #   # filtering dfs of interest for the weeks where a flag was detected
        #   site_df <- site_flag_dates %>%
        #     filter(y_w == group_data$y_w)
        #
        #   # TryCatch used here to avoid erroring out on the first and last values of
        #   # sites_order object (there is no prior/next record after the first/last record).
        #   # Return df as NULL in case of an error
        #   prev_site_df <- NULL
        #   next_site_df <- NULL
        #
        #   tryCatch({
        #     previous_site <- paste0(sites_order[site_index-1],"-",parameter_arg)
        #     prev_site_df <- all_df_list[[previous_site]] %>%
        #       filter(y_w == group_data$y_w)},
        #     error = function(err) {
        #       cat("No previous site.\n")})
        #
        #   tryCatch({
        #     next_site <- paste0(sites_order[site_index+1],"-",parameter_arg)
        #     next_site_df <- all_df_list[[next_site]] %>%
        #       filter(y_w == group_data$y_w)},
        #     error = function(err) {
        #       cat("No next site.\n")})
        #
        #   # Bind all three dfs
        #   week_plot_data <- list(site_df, prev_site_df, next_site_df) %>%
        #     # remove NULL values from the list
        #     keep(~ !is.null(.)) %>%
        #     bind_rows()
        #
        #   # Create a sequence of dates for the vertical lines
        #   start_date <- floor_date(min(week_plot_data$DT_round), "day")
        #   end_date <- ceiling_date(max(week_plot_data$DT_round), "day")
        #   vline_dates <- seq(start_date, end_date, by = "day")
        #
        #   date_seq <- seq(min(vline_dates), max(vline_dates) - days(1), by = "1 day") + hours(12) ## here ----
        #
        #   # Use the first day of the group as flag_day
        #   flag_day <- min(group_data$DT_round)
        #
        #   week_plot <- ggplot(data = week_plot_data) +
        #     geom_point(data = filter(week_plot_data, (site == site_arg)),
        #                aes(x=DT_round, y=mean, color=flag)) +
        #     geom_line(data = filter(week_plot_data, (site != site_arg)),
        #               aes(x=DT_round, y=mean, color=site)) +
        #     geom_vline(xintercept = vline_dates, color = "black") +
        #     ggtitle(paste0(str_to_title(site_arg), " ", parameter_arg, " (", format(flag_day, "%B %d, %Y"), ")")) +
        #     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        #     labs(x = "Day",
        #          y = "Mean")
        #
        #   week_plot <- add_threshold_lines(plot = week_plot,
        #                                    plot_data = week_plot_data,
        #                                    site_arg = site_arg,
        #                                    parameter_arg = parameter_arg)
        #
        #   week_plot <- week_plot +
        #     theme_bw() +
        #     scale_x_datetime(date_breaks = "1 day",
        #                      date_labels = "%b %d",
        #                      minor_breaks = date_seq) +
        #     theme(legend.position = 'bottom',
        #           panel.grid.major = element_blank(),
        #           panel.grid.minor = element_blank()) +
        #     annotate("text", x = date_seq, y = min(week_plot_data$mean, na.rm = TRUE) - 1, label = 1:length(date_seq), hjust = 0) +
        #     guides(color = guide_legend(nrow = 4, byrow = TRUE))
        #
        #   plot_list[[paste(site_param, as.character(flag_day))]] <- week_plot
        #   sorted_plot_names <- names(plot_list)[order(names(plot_list))]
        #
        #   plot_list <- plot_list[sorted_plot_names]
        # }
      }
      return(plot_list)
    } else {
      return(paste(flag_arg, "not detected.\n"))
    }
  } else {
    return(paste(site_arg, parameter_arg, "combination not available.\n"))
  }

}
# Generate seasonal histograms for a given site parameter data frame based on the mean column.
# @param df A dataframe with a `mean` column.
# @param df_index The index of the dataframe.
# @return A plot with a histogram of the mean column for each season.
# @examples
# generate_seasonal_histogram(df = all_data_flagged$`archery-Actual Conductivity`, df_index = "archery-Actual Conductivity")
# generate_seasonal_histogram(df = all_data_flagged$`boxelder-Temperature`, df_index = "boxelder-Temperature")

generate_seasonal_histogram <- function(df, df_index) {

  winter_baseflow <- c(12,1,2,3,4)
  snowmelt <- c(5,6,NA,NA,NA)
  monsoon <- c(7,8,9,NA,NA)
  fall_baseflow <- c(10,11,NA,NA,NA)

  # do water years

  seasons <- data.frame(winter_baseflow, snowmelt, monsoon, fall_baseflow)

  site_param <- toupper(sub("-", " ", df_index, fixed = TRUE))

  param <- unique(na.omit(df$parameter))

  hist_list <- list()
  for (i in colnames(seasons)){

    filtered_df <- df %>%
      filter(month %in% seasons[[i]],
             !str_detect(flag, "sensor specification range"))

    n <- as.character(sum(!is.na(filtered_df$mean)))

    title <- paste0(i," (n = ",n ,")")

    histogram <- ggplot() +
      geom_histogram() +
      labs(title = title)

    tryCatch({
      minimum <- floor(min(filtered_df$mean, na.rm = TRUE))
      maximum <- ceiling(max(filtered_df$mean, na.rm = TRUE))

      if (param %in% c("Specific Conductivity", "Actual Conductivity", "Turbidity")) {
        bins <- seq(minimum, maximum, by = 10)
      } else if (param %in% c("ORP", "pH", "Depth")) {
        bins <- seq(minimum, maximum, by = 0.05)
      } else {
        bins <- seq(minimum, maximum)
      }

      x_min <- filtered_df$m_mean05[1]
      x_max <- filtered_df$m_mean99[1]

      histogram <- ggplot(data = filtered_df, aes(x = mean)) +
        geom_histogram(breaks = bins) +
        geom_vline(xintercept = x_min, color = "red", linetype = "dashed") +
        geom_vline(xintercept = x_max, color = "red", linetype = "dashed") +
        facet_wrap(~ year, nrow = 1) +
        labs(title = title)
    },
    error = function(err) {
      cat("No finite values for", site_param, i, "\n")})

    hist_list[[i]] <- histogram

  }

  collated_hist <- ggarrange(plotlist = hist_list, nrow=2, ncol=2) %>%
    annotate_figure(top = site_param)

  return(collated_hist)

}
#' @title Generate site metaparameter data from the HydroVu API.
#'
#' @description
#' A function that generates a site metaparameter dataframe from the HydroVu API.
#' A metaparameter is a parameter that is used to generate flags for parameters
#' other than itself.
#'
#' @param api_data A dataframe with the munged API data.
#'
#' @param require An upstream dependency (i.e., alternative to `{necessary object}`).
#'
#' @return A dataframe with site metaparameter data that will be joined to the other
#' site-parameter data.
#'
#' @examples
#' generate_site_metaparam(site_arg = "archery", metaparameter_arg = "Battery Level", api_data = incoming_data_collated_csvs)
#'
#' @seealso [summarize_site_param()]

generate_site_metaparam <- function(api_data, require = NULL) {

    sites <- unique(api_data$site)
    metaparameters <- c("Temperature", "Battery Level", "Baro", "External Voltage")

    df_list <- list()
    for (i in sites) {
        metaparameter_data <- api_data %>%
            data.table::data.table() %>%
            dplyr::select(DT_join, site, parameter, value) %>%
            dplyr::filter(site == i & (parameter %in% metaparameters)) %>%
            dplyr::select(-site) %>%
            tidyr::pivot_wider(names_from = parameter, values_from = value)
        df_list[[i]] <- metaparameter_data
    }
    return(df_list)
}
#' @title Generate Summary Statistics
#'
#' @description
#' A function that generates summary statistics for a given site parameter data frame.
#' The generated statistics include:
#'   - The next value and previous value for the mean.
#'   - The rolling 7-point median of the mean.
#'   - The rolling 7-point mean of the mean.
#'   - The rolling 7-point standard deviation of the mean.
#'   - The slope of a point in relation to the point ahead and behind.
#'   - The rolling 7-point slope of the mean.
#'   - The month and year of each data point.
#'   - The year-month combination.
#'
#' @param site_param_df A data frame with a `mean` column retrieved from HydroVu API.
#'
#' @return A data frame with summary statistics for a given site parameter data frame.
#'
#' @examples
#' generate_summary_statistics(site_param_df = all_data_flagged$`archery-Actual Conductivity`)
#' generate_summary_statistics(site_param_df = all_data_flagged$`boxelder-Temperature`)

generate_summary_statistics_full <- function(site_param_df) {

  summary_stats_df <- site_param_df %>%
    # ... so that we can get the proper leading/lagging values across our entire timeseries:
    dplyr::mutate(
      # Add the next value and previous value for mean.
      front1 = dplyr::lead(mean, n = 1),
      back1 = dplyr::lag(mean, n = 1),
      # Add the rolling 7-point median (using itself + data of the past).
      rollmed = RcppRoll::roll_median(mean, n = 7, align = 'right', na.rm = F, fill = NA_real_), #ifelse(is.na(rollmed), roll_median(mean, n = 7, align = 'right', na.rm = F, fill = NA_real_), rollmed), # to go (j): check_na() function for when we append data
      # Add the rolling 7-point mean (using itself + data of the past).
      rollavg = RcppRoll::roll_mean(mean, n = 7, align = 'right', na.rm = F, fill = NA_real_), #ifelse(is.na(rollavg), roll_mean(mean, n = 7, align = 'right', na.rm = F, fill = NA_real_), rollavg),
      # Add the rolling 7-point standard deviation (using itself + data of the past).
      rollsd = RcppRoll::roll_sd(mean, n = 7, align = 'right', na.rm = F, fill = NA_real_), #ifelse(is.na(rollsd), roll_sd(mean, n = 7, align = 'right', na.rm = F, fill = NA_real_), rollsd),
      # Determine the slope of a point in relation to the point ahead and behind
      slope_ahead = (front1 - mean)/15,
      slope_behind = (mean - back1)/15,
      # Get the rolling 7-point slope (using itself + data of the past).
      rollslope = RcppRoll::roll_mean(slope_behind, n = 7, align = 'right', na.rm = F, fill = NA_real_), #ifelse(is.na(rollslope), roll_mean(slope_behind, n = 7, align = 'right', na.rm = F, fill = NA_real_), rollslope),
      # add some summary info for future steps
      month = lubridate::month(DT_round),
      year = lubridate::year(DT_round),
      y_m = paste0(year, '-', month),
      # Define our seasons:
      season = dplyr::case_when(month %in% c(12, 1, 2, 3, 4) ~ "winter_baseflow",
                                month %in% c(5, 6) ~ "snowmelt",
                                month %in% c(7, 8, 9) ~ "monsoon",
                                month %in% c(10, 11) ~ "fall_baseflow",
                                TRUE ~ NA)
    )

  return(summary_stats_df)

}
# Generate summary statistics for a given site parameter data frame.
#' @param site_param_df A data frame with a `mean` column retrieved from HydroVu API.
#' @return A data frame with summary statistics for a given site parameter data frame.
#' @examples
# generate_summary_statistics(site_param_df = all_data_flagged$`archery-Actual Conductivity`)
# generate_summary_statistics(site_param_df = all_data_flagged$`boxelder-Temperature`)

generate_summary_statistics <- function(site_param_df) {

  summary_stats_df <- site_param_df %>%
    # ... so that we can get the proper leading/lagging values across our entire timeseries:
    dplyr::mutate(
      # Add the next value and previous value for mean.
      # Only do this for newest data (i.e., our appended historical
      # data already has these filled out and we don't want to over-
      # write them)
      front1 = ifelse(is.na(front1), dplyr::lead(mean, n = 1), front1),
      back1 = ifelse(is.na(back1), dplyr::lag(mean, n = 1), back1),
      # Add the median for a point and 6 points behind it:
      rollmed = ifelse(is.na(rollmed), RcppRoll::roll_median(mean, n = 7, align = 'right', na.rm = F, fill = NA_real_), rollmed),
      # Add the mean for a point and 6 points behind it:
      rollavg = ifelse(is.na(rollavg), RcppRoll::roll_mean(mean, n = 7, align = 'right', na.rm = F, fill = NA_real_), rollavg),
      # Add the standard deviation for a point and 6 points behind it:
      rollsd = ifelse(is.na(rollsd), RcppRoll::roll_sd(mean, n = 7, align = 'right', na.rm = F, fill = NA_real_), rollsd),
      # Determine the slope of a point in relation to the point ahead and behind.
      slope_ahead = ifelse(is.na(slope_ahead), (front1 - mean)/15, slope_ahead),
      slope_behind = ifelse(is.na(slope_behind), (mean - back1)/15, slope_behind),
      # Get the rolling 7-point slope (using itself + data of the past).
      rollslope = ifelse(is.na(rollslope), RcppRoll::roll_mean(slope_behind, n = 7, align = 'right', na.rm = F, fill = NA_real_), rollslope),
      # add some summary info for future us
      month = lubridate::month(DT_round),
      year = lubridate::year(DT_round),
      y_m = paste(year, '-', month),
      season = dplyr::case_when(month %in% c(12, 1, 2, 3, 4) ~ "winter_baseflow",
                                month %in% c(5, 6) ~ "snowmelt",
                                month %in% c(7, 8, 9) ~ "monsoon",
                                month %in% c(10, 11) ~ "fall_baseflow",
                                TRUE ~ NA)
    )

  return(summary_stats_df)

}
generate_supplemental_weekly_plot <- function(daily_plot_data_arg, df_list_arg, site_arg, parameter_arg) {

  site_param_df <- df_list_arg[[paste0(site_arg, "-", parameter_arg)]]

  start_date <- min(daily_plot_data_arg$DT_round) - days(3) # TODO: replace all daily_plot_data_arg instances with something more efficient
  end_date <- max(daily_plot_data_arg$DT_round) + days(3)

  week_plot_data <- site_param_df %>%
    filter(DT_round %within% interval(start_date, end_date))

  #default is lower network
  site_vector <- c("tamasag", # rist
                   "legacy",
                   "lincoln",
                   "timberline",
                   "prospect",
                   "boxelder", # elc
                   "archery",
                   "river bluffs")

  if(network == "virridy"){ # this will be the new default
    # establish order for all the non-tributary sites
    sites_order <-  c("joei","cbri","chd","pfal","sfm","pbd","tamasag",
                      "legacy","lincoln","timberline","prospect","boxelder",
                      "archery","riverbluffs")
    # establish the order for the tributary sites
    trib_sites_order <- c("boxcreek", "archery", NA, "springcreek", "prospect",
                          NA, "penn", "sfm", "lbea")
  }

  # determining the index for the site of interest.
  if (site_arg %in% sites_order) {

    plot_filter <- tibble(site = c("joei","cbri","chd","pfal","sfm","pbd",
                                   "tamasag","legacy", "lincoln","timberline",
                                   "timberline virridy","prospect",
                                   "prospect virridy","boxelder","archery",
                                   "archery virridy","riverbluffs"))

    site_index <- which(sites_order == site_arg)
    site_list <- as.vector(na.omit(sites_order[max(1, site_index - 1):min(length(sites_order), site_index + 1)]))

    plot_filter <- plot_filter %>%
      filter(grepl(paste(site_list, collapse = "|"), site, ignore.case = TRUE),
             site != site_arg) %>%
      pull(site)

  } else {

    plot_filter <- tibble(site = c("boxcreek", "archery", "archery virridy",
                                   "springcreek", "prospect", "prospect virridy",
                                   "penn", "sfm", "lbea"))

    site_index <- which(trib_sites_order == site_arg)
    site_list <- as.vector(na.omit(trib_sites_order[max(1, site_index - 1):min(length(trib_sites_order), site_index + 1)]))

    plot_filter <- plot_filter %>%
      filter(grepl(paste(site_list, collapse = "|"), site, ignore.case = TRUE),
             site != site_arg) %>%
      pull(site)

  }

  # get the relevant sonde data source
  relevant_sonde_source <- map(plot_filter, ~ {
    sonde_name <- paste0(.x, "-", parameter_arg)
    # Determine which directory to pull data from
    tryCatch({
      retrieve_relevant_data_name(sonde_name, interval_arg = interval(start_date, end_date))
    }, error = function(err) {
      return("all_data")
    })
  })

  # Get the relevant data
  relevant_sondes <- map2(plot_filter,
                          relevant_sonde_source,
    function(name, source) {
    sonde_name <- paste0(name, "-", parameter_arg)
    # try to pull in the data
    tryCatch({
      get(source)[[sonde_name]] %>%
        filter(DT_round %within% interval(start_date, end_date))
    }, error = function(err) {
      return(NULL)
    })
  })

  # combine the lists
  sonde_info <- map2(relevant_sondes, relevant_sonde_source, list)

  # Remove any NULL results from the list
  sonde_info <- keep(sonde_info, ~!is.null(.x[[1]]))

  # append site_df to relevant sonde list, clean list, and bind dfs
  # to find plot info
  relevant_dfs <- map(sonde_info, ~.x[[1]])

  week_plot_data <- append(relevant_dfs, list(week_plot_data)) %>%
    keep(~ !is.null(.)) %>%
    keep(~ nrow(.)>0) %>%
    bind_rows() %>%
    arrange(day)

  # Create a sequence of dates for the vertical lines
  start_date <- floor_date(min(week_plot_data$DT_round), "day")
  end_date <- ceiling_date(max(week_plot_data$DT_round), "day")
  vline_dates <- seq(start_date, end_date, by = "day")

  date_seq <- seq(min(vline_dates), max(vline_dates) - days(1), by = "1 day") + hours(12)

  # use the daily flag data day as flag_day
  flag_day <- unique(daily_plot_data_arg$DT_round)

  week_plot <- ggplot(data = week_plot_data) +
    geom_point(data = filter(week_plot_data, (site == unique(daily_plot_data_arg$site))),
               aes(x=DT_round, y=mean, color=flag)) +
    map(sonde_info, function(sonde_data) {

      data <- sonde_data[[1]]
      data_source <- sonde_data[[2]]

      y_column <- if (data_source == "all_data") "mean" else "mean_verified"

      geom_line(data = data, aes(x = DT_round, y = .data[[y_column]], color = site))
    }) +
    geom_rect(data = daily_plot_data_arg, aes(xmin = min(DT_round), xmax = max(DT_round),
                                              ymin = -Inf, ymax = Inf),
              fill = "grey",
              alpha = 0.01,
              color = NA) +
    geom_vline(xintercept = vline_dates, color = "black") +
    ggtitle(paste0(str_to_title(unique(daily_plot_data_arg$site)), " ", unique(daily_plot_data_arg$parameter), " (", format(flag_day, "%B %d, %Y"), ")")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Day",
         y = "Mean")

  week_plot <- add_threshold_lines(plot = week_plot,
                                   plot_data = week_plot_data,
                                   site_arg = unique(daily_plot_data_arg$site),
                                   parameter_arg = unique(daily_plot_data_arg$parameter))

  week_plot <- week_plot +
    theme_bw() +
    scale_x_datetime(date_breaks = "1 day",
                     date_labels = "%b %d",
                     minor_breaks = date_seq) +
    theme(legend.position = 'right',
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) + # make the background white
    guides(color = guide_legend(nrow = 10, byrow = TRUE))

  return(week_plot)
}
# Get day decision.
# @param prompt_text A string for the prompt text.
# @return A string of "pass", "fail", or "inspect" depending on input from the user.
# @examples
# get_day_decision(prompt_text = "Would you like to (pass/fail/inspect) all data points for: ")

get_day_decision <- function(prompt_text) {
  while (TRUE) {
    user_input <- readline(prompt = paste(prompt_text))
    user_input <- tolower(user_input)

    if (user_input %in% c("pass", "p")) {
      return("pass")
    } else if (user_input %in% c("fail", "f")) {
      return("fail")
    } else if (user_input %in% c("inspect", "i")){
      return("inspect")
    } else {
      cat("Invalid input. Please enter 'yes', 'no', or 'inspect'.\n")
    }
  }
}
get_dt_inspection_decisions <- function(daily_plot_data) {
  # Create a sequence of date times from the daily plot data
  dt_list <- sort(unique(pull(daily_plot_data, DT_round)))

  # Get the date
  plot_date <- as_date(daily_plot_data$DT_round[1])

  # Prompt for the number of intervals
  prompt_text <- paste("How many intervals would you like to inspect?\n(Must be less than", length(dt_list), ")\n")

  while (TRUE) {
    user_input <- readline(prompt = prompt_text)

    # The user input needs to check for:
    # Only numbers between 1 and length(dt_list)
    num_intervals <- suppressWarnings(as.numeric(user_input))

    if (!is.na(num_intervals) && num_intervals > 0 && num_intervals < length(dt_list)) {
      break
    }

    cat("Invalid input. Please enter a valid number of intervals.\n")
  }

  selected_intervals <- list()

  for (i in 1:num_intervals) {
    prompt_text <- paste("Enter the time range for interval", i, ".\n(format 'HH:MM:SS-HH:MM:SS'):\n")

    while (TRUE) {
      user_input <- readline(prompt = prompt_text) # these need to be 15 minute intervals ***
      interval <- unlist(strsplit(user_input, "-"))

      interval_set <- paste(plot_date, interval)
      interval_set <- as_datetime(interval_set, tz = "MST")

      if (length(interval_set) == 2 && all(interval_set %in% dt_list)) {
        start_time <- interval_set[1]
        end_time <- interval_set[2]

        # Check if selected_intervals is empty
        if (start_time <= end_time & length(selected_intervals) == 0) {
          selected_intervals[[i]] <- interval(start_time, end_time)
          break
        }

        if (start_time <= end_time & !suppressWarnings(any(map(selected_intervals, ~int_overlaps(interval(start_time, end_time), .x))))) {
          selected_intervals[[i]] <- interval(start_time, end_time)
          break
        }
      }

      cat("Invalid input.\nPlease enter a valid time range that doesn't overlap with previous intervals.\n")
    }
  }

  return(selected_intervals)
}
#' #' Import Fort Collins Floodwarning Rain Gage Data
#' #'
#' #' This function pulls rain data for sites within the FC Floodwarning System
#' #'
#' #'
#' #' @param start Start date of when you want data.
#' #' @param end End date of when you want data; default is the current date.
#' #' @param save Whether to save (TRUE) the resulting table or not (FALSE)
#' #' @param path If `save = TRUE`, the file path to save the shapefile
#' #'
#' #' @return A table of time series flow data across the FC Floodwarning System
#' #'
#' get_fc_rain <- function(start = '2018-10-01', end = Sys.Date(), save = TRUE, path = 'data/context_data/'){
#'
#'   call <- "https://opendata.fcgov.com/api/views/g87z-rviz/rows.csv?accessType=DOWNLOAD"
#'
#'   #download dataset from FC
#'   temp1 <- tempfile()
#'   download.file(paste0(call), destfile = temp1, method = "curl")
#'   data <- read_csv(temp1) %>%
#'     dplyr::mutate(date = as_date(mdy_hms(Timestamp))) %>%
#'     dplyr::filter(date >= start & date <= end)
#'
#'   s_platte_ws s_platte_flowlines <- nhdplusTools::get_nhdplus(AOI = s_platte_ws,
#'                                                               realization = "flowline") <- nhdplusTools::get_nldi_basin(list(featureSource = "nwissite", featureID = "USGS-06754000"))
#'   s_platte_flowlines <- nhdplusTools::get_nhdplus(AOI = s_platte_ws,
#'                                                   realization = "flowline")
#'   s_platte_catchments <- nhdplusTools::get_nhdplus(AOI = s_platte_ws, realization = "catchment")
#'
#'   poudre <- s_platte_flowlines %>%
#'     filter(grepl("Cache la Poudre", gnis_name, ignore.case = TRUE))
#'
#'   poudre <- s_platte_catchments %>%
#'     filter(featureid %in% poudre$comid) %>%
#'     .[sites,] %>%
#'     rowid_to_column()
#'
#'   sites <- read_csv("data/metadata/sonde_location_metadata.csv") %>%
#'     separate(lat_long, into = c("lat", "long"), sep = ",") %>%
#'     sf::st_as_sf(coords = c("long", "lat"), crs = 4326) %>%
#'     filter(Site %in% c("Tamasag","Legacy","Lincoln","Timberline","Prospect","Boxelder","Archery","River Bluffs")) %>%
#'     mutate(rowid = sf::st_nearest_feature(., poudre)) %>%
#'     left_join(sf::st_drop_geometry(poudre), by = "rowid")
#'
#'   ws_maker <- function(site_names){
#'
#'     df <- sites %>%
#'       dplyr::filter(Site == site_names)
#'
#'     # get_UT() creates a vector of all flowlines upstream of the comid of interest ...
#'     UT_comids <- nhdplusTools::get_UT(network = s_platte_flowlines,
#'                                       comid = df$featureid)
#'
#'     catchments <- filter(s_platte_catchments, featureid %in% c(UT_comids)) %>%
#'       summarize() %>%
#'       mutate(site = df$Site)
#'
#'     return(catchments)
#'
#'   }
#'
#'   watersheds <- sites$Site %>%
#'     map_dfr(~ws_maker(.))
#'
#'   rain_sites <- fread("data/context_data/fc_rain.csv") %>%
#'     distinct(`Sensor Name`, Latitude, Longitude) %>%
#'     sf::st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>%
#'     .[watersheds,] %>%
#'     sf::st_join(watersheds) %>%
#'     mutate(site = tolower(site)) %>%
#'     write_csv("data/context_data/site_link.csv")
#'
#'   if(save == TRUE){
#'     write_csv(data, paste0(path, "/fc_rain.csv"))
#'   }
#'
#'   return(data)
#' }
#'
#'
#'
#'
# Get flag decision.
# @param prompt_text A string for the prompt text.
# @return A boolean of TRUE or FALSE depending on input from the user.
# @examples
# get_flag_decision(prompt_text = "Would you like to (pass/fail) this data point: ")

get_flag_decision <- function(prompt_text) {
  while (TRUE) {
    user_input <- readline(prompt = paste(prompt_text))
    user_input <- tolower(user_input)

    if (user_input %in% c("pass", "p")) {
      return(TRUE)
    } else if (user_input %in% c("fail", "f")) {
      return(FALSE)
    } else {
      cat("Invalid input. Please enter 'pass' or 'fail'.\n")
    }
  }
}
#' @title Get the start dates for the HydroVu API pull from the historically
#' flagged data
#'
#' @description
#' This function finds the max datetime in the temperature data
#' frames for each site's historically flagged data. We use temperature because
#' it is always being tracked in the sensors and thus is the most accurate
#' representation of the last API pull for a given sensor.
#'
#' @param incoming_historically_flagged_data_list This is the output from
#' `flagged_data_dfs` and is the historically flagged data.
#'
#' @returns
#' A data frame that has the sites and their corresponding starting
#' date-times' for the HydroVu API pull. This data frame will then get passed into
#' `incoming_data_csvs_upload`.

get_start_dates_df <- function(incoming_historically_flagged_data_list) {

  temperature_subset <- grep("Temperature", names(incoming_historically_flagged_data_list))

  start_dates_df <- map_dfr(incoming_historically_flagged_data_list[temperature_subset],
                            function(temperature_df){
                              temperature_df %>%
                                select(DT_round, site) %>%
                                filter(DT_round == max(DT_round))
                              })

  return(start_dates_df)

}
# @ param site: string, name of site to get filenames for, select from "tamasag" "legacy", "timberline" "prospect" "boxelder" "archery" "riverbluffs"
get_tl_photo_filenames <- function(site = "legacy",start_dt = "2023-05-20 12:00" ,end_dt = "2023-07-20 12:00"){
  #convert dates to DT objects
  start_dt <- ymd_hm(start_dt, tz = "MST")
  end_dt <- ymd_hm(end_dt, tz = "MST")
  if(site == "river bluffs"){
    site <- "riverbluffs"
  }

  if(site %nin% c("tamasag", "legacy", "timberline", "prospect", "boxelder", "archery", "riverbluffs")){
    print("Invalid site name, please select from: 'tamasag' 'legacy' 'timberline' 'prospect' 'boxelder' 'archery' 'riverbluffs'")

  nah <- tibble()
    return(nah)
  }

  #look through all the files in the site's folder
  all_files <- tibble(filename = list.files(paste0("data/timelapse_photos/2023_compiled/",site), full.names = TRUE, recursive = TRUE))%>%
    mutate(DT = str_extract(filename, "[0-9]{8}_[0-9]{4}"),
           DT = parse_date_time(DT, "Ymd_HM", tz = "MST"),
           #round DT to 30 min to match with sensor data
           DT_round = round_date(DT, "30 minutes"))%>%
    # if there are multiple photos with the same DT_round value, only keep one row with that DT_round
    distinct(DT_round, .keep_all = TRUE)%>%
    select(-DT)%>%
    filter(between(DT_round, start_dt, end_dt))

  if(nrow(all_files) == 0){
    print("No files found for this site and date range")
    nah <- tibble()
    return(nah)
  }else{
   return(all_files)
  }


}

#test <- get_tl_photo_filenames(site = "legacy", start_dt = "2023-7-20 12:00", end_dt = "2023-12-06 12:00")




#' @title Get verification decision
#'
#' @description
#' Get decisions from user input
#'
#' @param prompt_text The prompt that the user is replying to
#'
#' @examples
#' # get_day_decision(prompt_text = "Would you like to (pass/fail/inspect) all data points for: ")

get_verification_decision <- function(prompt_text) { # should add a layer level to this to prevent user from inspect sub-daily data

  while (TRUE) {
    user_input <- readline(prompt = paste(prompt_text))
    user_input <- tolower(user_input)

    # pass statements
    if (user_input %in% c("pass all", "pass", "pa", "p")) {
      return("PASS ALL")
    }
    if (user_input %in% c("pass valid", "pv")) {
      return("PASS VALID")
    }
    if (user_input %in% c("pass flagged", "pf")) {
      return("PASS FLAGGED")
    }
    if (user_input %in% c("pass none", "pn")) {
      return("PASS NONE")
    }

    # fail statements
    if (user_input %in% c("fail all", "fail", "fa", "f")) {
      return("FAIL ALL")
    }
    if (user_input %in% c("fail valid", "fv")) {
      return("FAIL VALID")
    }
    if (user_input %in% c("fail flagged", "ff")) {
      return("FAIL FLAGGED")
    }
    if (user_input %in% c("fail none", "fn")) {
      return("FAIL NONE")
    }

    # skip statements
    if (user_input %in% c("skip", "skip all", "s", "sa")) {
      return("SKIP")
    }

    # inspect statements ***
    if (user_input %in% c("inspect all", "inspect", "ia", "i")){
      return("INSPECT ALL")
    }

    if (user_input %in% c("inspect valid", "iv")){
      return("INSPECT VALID")
    }

    if (user_input %in% c("inspect flagged", "if")){
      return("INSPECT FLAGGED")
    }

    if(user_input %in% c("inspect some", "is")){ # if we are at daily level we can't go further
      return("INSPECT SOME")
    }

    # quit statement
    if (user_input %in% c("quit", "q")) {
      return("QUIT")
    }

    cat("Invalid input. Please enter one of the options from the following decision matrix:\n\n")
    cat("          | PASS | FAIL | SKIP | INSPECT  |  | QUIT |\n")
    cat("----------+------+------+------+----------|  |------|\n")
    cat(" ALL      | pa   | fa   | sa   | ia       |  | q    |\n")
    cat(" VALID    | pv   | fv   |      | iv       |  +------+\n")
    cat(" FLAGGED  | pf   | ff   |      | if       |\n")
    cat(" NONE     | pn   | fn   |      |          |\n")
    cat(" SOME     |      |      |      | is       |\n")
    cat("----------+------+------+------+----------|\n")
  }

}

get_watersheds <- function(sf = NULL, coordinates = c(-112.2124, 37.63281)){

  sf::sf_use_s2(FALSE)

  if(is.null(sf)){

    # Create a data frame with a column named 'geometry'
    df <- tibble::tibble(long = coordinates[1],
                         lat = coordinates[2])

    aoi <- sf::st_as_sf(df, coords = c("long", "lat"), crs = 4326)
  }

  if(is.null(coordinates)){

    aoi <- sf %>% sf::st_transform(4326)

  }

  aoi <- aoi %>%
    nhdplusTools::get_nhdplus(AOI = .)

  # Use the NHDPlus digital elevation model to find the nearest downslope
  # NHD flowline for any point in space (here, our point of interest)
  # Not working anymore -
  # trace <- get_raindrop_trace(aoi, direction = "up")
  #
  # "Snap" our site to the nearest NHD flowline feature
  # snap_point <- sf::st_sfc(sf::st_point(trace$intersection_point[[1]][1:2]),
  #                          crs=4326)
  #
  # Clip/split our catchment to only include the portion of the
  # catchment upstream of our site:
  # better_termination <- get_split_catchment(snap_point, upstream = F)[2,]

  # read in the complete NHD (in tabular form) to make for much more efficient nhd crawling.
  # This data in tabular form doesn't exist anywhere online that I know of... -_-
  nhd <- readr::read_csv('data/context_data/nhd_flow_network.csv')

  upstream <- nhdplusTools::get_UT(nhd, aoi$comid) %>% #upstream trace function in nhdplusTools
    tibble::as_tibble() %>%
    dplyr::rename(comid_list = value)  %>%
    dplyr::distinct(comid_list, .keep_all = TRUE)

  nhd_catch <-  upstream$comid_list %>%
    map(~nhdplusTools::get_nhdplus(comid = .,
                                   realization='catchment',
                                   t_srs = 4269)) %>%
    dplyr::bind_rows() %>%
    dplyr::distinct(featureid,.keep_all=TRUE) %>%
    dplyr::group_by(featureid) %>%
    dplyr::summarize()

  return(nhd_catch)

}


get_weekly_inspection_decision <- function(weekly_plot_data) {

  prompt_text <- "Which days do you want to verify (1/2/3/4/5/6/7)? \nCan be any combination of weekdays, with no other characters.\nex: 456\n"
  day_min <- min(weekly_plot_data$weekday)
  day_max <- max(weekly_plot_data$weekday)

  while (TRUE) {

    user_input_og <- readline(prompt = paste(prompt_text))

    user_input <- suppressWarnings(
      sort(na.omit(as.numeric(unique(unlist(strsplit(user_input_og, ""))))) [as.numeric(unique(unlist(strsplit(user_input_og, "")))) >= day_min & as.numeric(unique(unlist(strsplit(user_input_og, "")))) <= day_max])
      ) # 1 and 7 will be set to data min and max

    if (length(user_input) != 0) {
      return(user_input)
    }

    cat(user_input_og, " is not a valid input")
  }
}
#' @title Get working data decision
#'
#' @description
#' Get working directory and data decisions from user input
#'
#' @param prompt_text The prompt that the user is replying to
#'
#' @examples
#' # get_day_decision(prompt_text = "Would you like to (pass/fail/inspect) all data points for: ")

get_working_data_decision <- function() { # should add a layer level to this to prevent user from inspect sub-daily data

  prompt_text <- "Which directory are you working from? (pre/int): "

  while (TRUE) {
    user_input <- readline(prompt = paste(prompt_text))
    user_input <- tolower(user_input)

    if (user_input == "pre") {
      working_data <<- set_names(map(list.files(pre_verification_path, full.names = TRUE), readRDS), list.files(pre_verification_path))
      break()
    }

    if (user_input == "int") {
      working_data <<- set_names(map(list.files(intermediary_path, full.names = TRUE), readRDS), list.files(intermediary_path))
      break()
    }

    cat("Invalid input. Please enter one of the options from the following:\n\n")
    cat("pre = pre_verification_dir\n")
    cat("int = intermediary_dir\n")
  }

}
grab_confirmed_mWater_malfunctions <- function(field_notes){

    # Grab notes about sensor malfunction
    malfunction_records <- all_notes_cleaned %>%
      dplyr::filter(grepl("Sensor malfunction", visit_type, ignore.case = TRUE)) %>%
      dplyr::select(start_DT, site, crew, which_sensor_malfunction, malfunction_end_dt, notes = visit_comments)

    # write to csv
    write_csv(malfunction_records, "data/mWater_malfunciton_records.csv")

  }
grab_mWater_malfunction_notes <- function(mWater_api_data){

  # Grab notes about sensor malfunction
  malfunction_records <- mWater_api_data %>%
    filter(grepl("Sensor malfunction", visit_type, ignore.case = TRUE)) %>%
    select(start_DT, site, crew, which_sensor_malfunction, malfunction_end_dt, notes = visit_comments)

  #write to csv
  write_csv(malfunction_records, "data/mWater_malfunction_records.csv")

  parameters <- c("Battery Level",
                  "Baro",
                  "Chl-a Fluorescence",
                  "Depth",
                  "DO",
                  "External Voltage",
                  "ORP",
                  "pH",
                  "Specific Conductivity",
                  "Temperature",
                  "Turbidity")

  malfunction_records <- malfunction_records %>%
    # keep records relevant to {target} analysis
    select(start_DT, end_DT = malfunction_end_dt, site, parameter = which_sensor_malfunction, notes) %>%
    # match the text in the sensor column to the text in the target analysis
    separate_rows(parameter, sep = ", ") %>%
    mutate(
      parameter = case_when(
      parameter == "Chlorophyll a" ~ "Chl-a Fluorescence",
      parameter == "RDO" ~ "DO",
      parameter == "Conductivity" ~ "Specific Conductivity",
      .default = parameter
    ),
    site = case_when(
      site == "riverbluffs" ~ "river bluffs",
      .default = site
    )) %>%
    filter((is.na(parameter)) | (parameter %in% parameters))

  return(malfunction_records)

}
grab_mWater_sensor_malfunctions <- function(){

  # API Pull of mWater submitted notes

  # Grab API url from yml
  # Contact Sam Struthers if you need access
  creds <- yaml::read_yaml("src/mWater_collate/mWater_API.yml")
  api_url <- as.character(creds["url"])

  # Read in from API and tidy for downstream use

  # This is basic tidying of data set to:
  # correct datetime from UTC to Denver time (always MST)
  # correct columns where Other input is allowed (Site, visit type, photos downloaded, sensor malfunction)
  # Add rounded date time

  mal_notes <- readr::read_csv(url(api_url), show_col_types = FALSE) %>%
    dplyr::mutate(
      # start and end dt comes in as UTC -> to MST
      start_DT = lubridate::with_tz(lubridate::parse_date_time(start_dt, orders = c("%Y%m%d %H:%M:%S", "%m%d%y %H:%M", "%m%d%Y %H:%M", "%b%d%y %H:%M")), tz = "MST"),
      end_dt = lubridate::with_tz(lubridate::parse_date_time(end_dt, orders = c("%Y%m%d %H:%M:%S", "%m%d%y %H:%M", "%m%d%Y %H:%M", "%b%d%y %H:%M" )), tz = "MST"),
      malfunction_end_dt = with_tz(lubridate::parse_date_time(malfunction_end_dt, orders = c("%Y%m%d %H:%M:%S", "%m%d%y %H:%M", "%m%d%Y %H:%M", "%b%d%y %H:%M" )), tz = "MST"),
      date = as.Date(start_DT, tz = "MST"),
      start_time_mst = format(start_DT, "%H:%M"),
      sensor_pulled = as.character(sn_removed),
      sensor_deployed = as.character(sn_deployed),
      # If other is chosen, make site == other response
      site = ifelse(site == "Other (please specify)", tolower(stringr::str_replace_all(site_other, " ", "")), site),
      # When I changed the mWater survey, I accidentally introduced ??? in the place of Sensor Calibration option, fixing that here
      visit_type = dplyr::case_when(stringr::str_detect(visit_type, "\\?\\?\\?") ~ stringr::str_replace(string = visit_type,
                                                                                                        pattern =  "\\?\\?\\?",
                                                                                                        replacement = "Sensor Calibration or Check"),
                                    TRUE ~ visit_type),
      # Merging visit_type and visit type other
      visit_type = dplyr::case_when(stringr::str_detect(visit_type, "Other") ~ stringr::str_replace(string = visit_type,
                                                                                                    pattern =  "Other \\(please specify\\)",
                                                                                                    replacement = visit_type_other),
                                    TRUE ~ visit_type),
      # Merge sensor malfunction and sensor malfunction other
      which_sensor_malfunction = dplyr::case_when(stringr::str_detect(which_sensor_malfunction, "Other") ~ stringr::str_replace(string = which_sensor_malfunction,
                                                                                                                                pattern =  "Other \\(please specify\\)",
                                                                                                                                replacement = as.character(other_which_sensor_malfunction)),
                                                  TRUE ~ which_sensor_malfunction),
      # If other is chosen, make photos downloaded equal to response
      photos_downloaded = ifelse(photos_downloaded == "Other (please specify)", photos_downloaded_other, photos_downloaded),
      # Rounded start date time
      DT_round = lubridate::floor_date(start_DT, "15 minutes")) %>%
    # arrange by most recent visit
    dplyr::arrange(DT_round) %>%
    # Remove other columns
    dplyr::select(-c(photos_downloaded_other, visit_type_other, site_other, other_which_sensor_malfunction)) %>%
    dplyr::filter(grepl("sensor malfunction", visit_type, ignore.case = TRUE)) %>%
    dplyr::select(malfunction_start_dt = DT_round,  malfunction_end_dt, which_sensor_malfunction)

  return(mal_notes)

}
grab_mWater_sensor_notes <- function(mWater_api_data){

  # Sensor Notes

  # These are the notes that will be added to the QAQC workflow notes Most of the code in this chunk is to get the df to
  # match the one in the QAQC workflow It can be saved as a CSV or pulled directly into QAQC workflow
  # grab only notes where technician is interacting with sensor on site (excludes sensor malfunction notes)

  mWater_field_notes <- mWater_api_data %>%
    filter(grepl("Sensor",visit_type, ignore.case = TRUE) & !grepl("Sensor malfunction",visit_type, ignore.case = TRUE)) %>%
    # determining sonde employed status based on sensor_change
    mutate(sonde_employed = case_when(is.na(sensor_change)  ~ NA,
                                      sensor_change == "Swapped" ~ NA,
                                      sensor_change == "Pulled" ~ 1,
                                      sensor_change == "Deployed" ~ 0),
                                      #sensor_change %in% c("Swapped", "Deployed") ~ 1),

           #Sensor swapped notes
           sensor_swapped_notes = case_when(is.na(sensor_change)  ~ NA,
                                            sensor_change == "Pulled" &!is.na(sensor_pulled) ~ paste0("SN Removed: ", sensor_pulled),
                                            sensor_change == "Swapped" ~ paste0("SN Removed: ", sensor_pulled, " SN Deployed: ", sensor_deployed),
                                            sensor_change == "Deployed" ~ sensor_deployed),
           #Date/field season columns to match QAQC workflow
           DT_join = as.character(DT_round),
           field_season = year(DT_round),
           last_site_visit = DT_round,
           date = as.character(date)
    )%>%
    arrange(desc(DT_round))%>%
    #order columns in easily readable ordering
    select(site, crew, DT_round,sonde_employed,  sensors_cleaned, wiper_working, rdo_cap_condition, rdo_cap_replaced , ph_junction_replaced ,
           cals_performed, cal_report_collected , sensor_malfunction,sensor_pulled,sensor_deployed, sensor_swapped_notes,
           visit_type,start_time_mst,DT_join,  start_DT, end_dt,date,  visit_comments,photos_downloaded, field_season, last_site_visit)

  # back up to CSV
  # write_csv(sensor_notes, "data/mWater_sensor_field_notes.csv")

  # rm(all_notes_cleaned)
  return(mWater_field_notes)

}

hv_data_id <- function(loc_id, start_time = startdate, end_time = enddate, tz = timezone, token) {

  # convert the time to timestamp, convert to UTC for lookup in HydroVu
  start <- as.numeric(lubridate::with_tz(lubridate::ymd_hms(start_time, tz = tz), tzone = "UTC"))
  end <- as.numeric(lubridate::with_tz(lubridate::ymd_hms(end_time, tz = tz), tzone = "UTC"))

  # build the url
  url = "https://www.hydrovu.com/public-api/v1/locations/"
  url <- paste0(url, loc_id, "/data?endTime=", end, '&startTime=', start)

  req <- httr2::request(url)
  print(paste0('Trying site ', loc_id))
  try({
    resp <-  req %>% httr2::req_oauth_client_credentials(token) %>% httr2::req_perform()
    data <- list(resp %>% httr2::resp_body_json())
    h <- resp %>% httr2::resp_headers()

    while (!is.null(h[["X-ISI-Next-Page"]]))
    {
      resp <- req %>% httr2::req_headers("X-ISI-Start-Page" = h[["X-ISI-Next-Page"]]) %>%
        httr2::req_oauth_client_credentials(token) %>% httr2::req_perform()
      data <- c(data, list(resp %>% httr2::resp_body_json()))
      h <- resp %>% httr2::resp_headers()
    }

    # get the params and units
    params <- hv_names(token, return = "params")
    units <- hv_names(token, return = "units")

    # collapse the paginated date and clean up
    df <- purrr::map_dfr(data, flatten_page_params) %>%
      dplyr::mutate(timestamp = lubridate::with_tz(lubridate::as_datetime(timestamp, tz = "UTC"), tzone = tz),
                    Location = loc_id) %>%
      dplyr::inner_join(params, by = "parameterId") %>%
      dplyr::inner_join(units, by = "unitId") %>%
      dplyr::select(-parameterId, -unitId) %>%
      dplyr::arrange(Parameter, timestamp)

      return(df)
  })

}
#' Return the list of locations for the given client from HydroVu
#'
#' @param client a valid OAuth2 token such as returned from \code{hv_auth()}
#' @param url HydroVu url that lists the locations
#'
#' @return a dataframe listing all the locations visible to the client
#' @export
#'
#' @examples
#' \dontrun{
#' locs <- hv_locations(client)
#' }

hv_locations_all <- function(client,
                         url = "https://www.hydrovu.com/public-api/v1/locations/list") {

  req <- httr2::request(url)

  try({
  resp <-  req %>% httr2::req_oauth_client_credentials(client) %>% httr2::req_perform()
  locs <- list(resp %>% httr2::resp_body_json())
  h <- resp %>% httr2::resp_headers()

  while (!is.null(h[["X-ISI-Next-Page"]]))
  {
    resp2 <- req %>%
      httr2::req_headers("X-ISI-Start-Page" = h[["X-ISI-Next-Page"]]) %>%
      httr2::req_oauth_client_credentials(client) %>%
      httr2::req_perform()
    locs <- c(locs, list(resp2 %>% httr2::resp_body_json()))
    h <- resp2 %>% httr2::resp_headers()
  }
  # collapse the paginated date and clean up
  df <- flatten_df(locs) %>%
    select(-gps) %>%
    filter(!duplicated(.))
  return(df)
  })
}

#' @title intersensor_check
#'
#' @description
#' A function designed to reduce overflagging: if a slope violation occurs at the
#' same time as a slope violation in either depth or temperature, it is likely
#' not a sensor malfunction and instead a real product of the river.
#'
#' @param df An updated data frame with a `flag` column.
#'
#' @return A data frame with a `flag` column that has been updated with the
#' inter-sensor flag reduction step.
#'
#' @seealso [add_flag()]

intersensor_check <- function(df){

  # create a df of temperature for each site
  temperature <- df %>%
    data.table::data.table() %>%
    dplyr::select(DT_round, DT_join, parameter, mean, flag) %>%
    dplyr::filter(parameter == "Temperature") %>%
    dplyr::select(DT_join, Temperature = parameter, Temperature_flag = flag) %>%
    dplyr:: mutate(Temperature_front1 = dplyr::lead(Temperature_flag, n = 1),
                   Temperature_back1 = dplyr::lag(Temperature_flag, n = 1))

  # create a df of depth for each site
  depth <- df %>%
    data.table::data.table() %>%
    dplyr::select(DT_round, DT_join, parameter, mean, flag) %>%
    dplyr::filter(parameter == "Depth") %>%
    dplyr::select(DT_join, Depth = parameter, Depth_flag = flag) %>%
    dplyr:: mutate(Depth_front1 = dplyr::lead(Depth_flag, n = 1),
                   Depth_back1 = dplyr::lag(Depth_flag, n = 1))

  # add "temperature" and "depth" data columns to df:
  intersensors_checked <- df %>%
    dplyr::filter(!parameter %in% c("Depth", "Temperature")) %>%
    dplyr::left_join(., temperature, by = "DT_join") %>%
    dplyr::left_join(., depth, by = "DT_join") %>%
    # If either the depth or temperature have the same flag as a given parameter
    # identified at the same time (or one obs before/after), tag it
    dplyr::mutate(intersensored = dplyr::case_when(grepl("slope violation", flag) &
                                                     (grepl("slope violation", Depth_flag)   | grepl("slope violation", Temperature_flag)   |
                                                        grepl("slope violation", Depth_front1) | grepl("slope violation", Temperature_front1) |
                                                        grepl("slope violation", Depth_back1)  | grepl("slope violation", Temperature_back1)
                                                     ) ~ TRUE)) %>%
    dplyr::mutate(flag = ifelse(is.na(intersensored), flag, stringr::str_replace(flag, "slope violation", "")))

  final_checked_data <- df %>%
    dplyr::filter(parameter %in% c("Depth", "Temperature")) %>%
    # After using the temp and depth slope flags, remove that flagging entirely
    # from those parameters. We have yet to find an instance of the slope flag
    # capturing "fake" spikes in either of those data sets:
    dplyr::mutate(flag = stringr::str_replace(flag, "slope violation", "")) %>%
    dplyr::bind_rows(., intersensors_checked) %>%
    dplyr::select(-c(Depth, Depth_flag, Temperature, Temperature_flag))

  return(final_checked_data)

}
#' @title Load and tidy mWater field notes
#'
#' @description A function that uploads and cleans the field notes submitted to mWater.
#'
#' @param creds A .yml file with necessary credentials for accessing the field notes. Contact Sam Struthers if you need access.
#'
#' @return A dataframe with the field notes.

load_mWater_notes <- function(creds = yaml::read_yaml("creds/mWaterCreds.yml")){

  # API Pull of mWater submitted notes

  # Grab API url from yml
  # Contact Sam Struthers if you need access
  api_url <- as.character(creds["url"])

  # Read in from API and tidy for downstream use

  # This is basic tidying of data set to:
  # correct datetime from UTC to Denver time (MST)
  # correct columns where Other input is allowed (Site, visit type, photos downloaded, sensor malfunction)
  # Add rounded date time

  all_notes_cleaned <- readr::read_csv(url(api_url), show_col_types = FALSE) %>%
    dplyr::mutate(
      # start and end dt comes in as UTC -> to MST
      start_DT = lubridate::with_tz(lubridate::parse_date_time(start_dt, orders = c("%Y%m%d %H:%M:%S", "%m%d%y %H:%M", "%m%d%Y %H:%M", "%b%d%y %H:%M" )), tz = "MST"),
      end_dt = lubridate::with_tz(lubridate::parse_date_time(end_dt, orders = c("%Y%m%d %H:%M:%S", "%m%d%y %H:%M", "%m%d%Y %H:%M", "%b%d%y %H:%M" )), tz = "MST"),
      malfunction_end_dt = lubridate::with_tz(lubridate::parse_date_time(malfunction_end_dt, orders = c("%Y%m%d %H:%M:%S", "%m%d%y %H:%M", "%m%d%Y %H:%M", "%b%d%y %H:%M" )), tz = "MST"),
      date = as.Date(start_DT, tz = "MST"),
      start_time_mst = format(start_DT, "%H:%M"),
      sensor_pulled = as.character(sn_removed),
      sensor_deployed = as.character(sn_deployed),
      # If other is chosen, make site == other response
      site = ifelse(site == "Other (please specify)", tolower(stringr::str_replace_all(site_other, " ", "")), site),
      # When I changed the mWater survey, I accidentally introduced ??? in the place of Sensor Calibration option, fixing that here
      visit_type = dplyr::case_when(stringr::str_detect(visit_type, "\\?\\?\\?") ~ stringr::str_replace(string = visit_type,
                                                                               pattern =  "\\?\\?\\?",
                                                                               replacement = "Sensor Calibration or Check"),
                             TRUE ~ visit_type),
      # Merging visit_type and visit type other
      visit_type = dplyr::case_when(stringr::str_detect(visit_type, "Other") ~ stringr::str_replace(string = visit_type,
                                                                           pattern =  "Other \\(please specify\\)",
                                                                           replacement = visit_type_other),
                             TRUE ~ visit_type),
      # Merge sensor malfunction and sensor malfunction other
      which_sensor_malfunction = dplyr::case_when(stringr::str_detect(which_sensor_malfunction, "Other") ~ stringr::str_replace(string = which_sensor_malfunction,
                                                                                                       pattern =  "Other \\(please specify\\)",
                                                                                                       replacement = as.character(other_which_sensor_malfunction)),
                                           TRUE ~ which_sensor_malfunction),
      # If other is chosen, make photos downloaded equal to response
      photos_downloaded = ifelse(photos_downloaded == "Other (please specify)", photos_downloaded_other, photos_downloaded),
      # Rounded start date time
      DT_round = lubridate::floor_date(start_DT, "15 minutes")) %>%
    # arrange by most recent visit
    dplyr::arrange(DT_round)%>%
    # Remove other columns
    dplyr::select(-c(photos_downloaded_other,visit_type_other, site_other, other_which_sensor_malfunction ))

  return(all_notes_cleaned)

}
#' @title Load and tidy old field notes
#'
#' @description A function that uploads and cleans the field notes excel file. This function adds datetime
#' columns to the field notes dataframe and filters out field notes where the sensor
#' was not handled.
#'
#' @param filepath A file path to the raw field notes.
#'
#' @return A dataframe with the field notes.
#'
#' @examples
#' clean_old_field_notes(filepath = "data/sensor_field_notes.xlsx")

load_old_field_notes <- function(filepath){

  raw_field_notes <- readxl::read_excel(filepath)

  field_notes <- raw_field_notes %>%
    dplyr::mutate(start_DT = lubridate::ymd_hm(paste(date, start_time_mst), tz = "MST")) %>%
    dplyr::mutate(
      DT_round = lubridate::floor_date(start_DT, "15 minutes"),
      DT_join = as.character(DT_round),
      site = tolower(site),
      field_season = lubridate::year(DT_round),
      last_site_visit = DT_round) %>%
    dplyr::arrange(site, DT_round) %>%
    # rename instances of old names:
    dplyr::mutate(site = ifelse(site == "rist", "tamasag",
                         ifelse(site == "elc", "boxelder", site))) %>%
    # `sonde_employed` determines if the sonde is deployed or not. 0 = sonde deployed, 1 = sonde is not deployed
    mutate(sonde_employed = dplyr::case_when(!is.na(sensor_pulled) & !is.na(sensor_deployed) ~ 0,
                                             !is.na(sensor_pulled) & is.na(sensor_deployed) ~ 1,
                                             is.na(sensor_pulled) & !is.na(sensor_deployed) ~ 0,
                                             is.na(sensor_pulled) & is.na(sensor_deployed) ~ NA),
           end_dt  = as.POSIXct(NA, tz = "MST")) %>%
    # remove field dates where sensor was not handled:
    dplyr::filter(grepl("Sensor Cleaning or Check|Sensor Calibration", visit_type, ignore.case = TRUE))

  return(field_notes)

}
#' @title Generate Threshold Table for QAQC
#'
#' @description
#' A function designed to generate a threshold table for QAQC. This table
#' contains the thresholds for the mean, slope_behind, and standard deviation
#' of the mean for each site and season.
#'
#' @param df A data frame with a `flag` column.
#'
#' @return A data frame with the thresholds for the mean, slope_behind, and
#' standard deviation of the mean for each site and season.
#'
#' @examples
#' make_threshold_table(df = all_data_flagged$`archery-Actual Conductivity`)

make_threshold_table <- function(df){

  # sensor_malfunction_notes <- grab_mWater_malfunction_notes(mWater_api_data = load_mWater_notes())

  slope_down <- df %>%
    # REMOVE DATA WE KNOW TO BE ERRONEOUS:
    add_field_flag() %>%
    # incoroporate back-calibrated values for development of thresholds:
    fix_calibration() %>%
    add_malfunction_flag(df = ., malfunction_records = sensor_malfunction_notes) %>%
    add_spec_flag() %>%
    # Filter to remove times when the sonde was under freexzing conditions or not submerged (depth <= 0)
    add_frozen_flag() %>%
    # fix_depth_cal() %>%
    add_unsubmerged_flag() %>%
    # remove times when sonde was moved (to avoid crazy slopes)
    add_depth_shift_flag() %>%
    # Filter to only un-flagged data, and data that was collected when sonde
    # was actually deployed
    filter(is.na(flag) | sonde_employed != 1) %>%
    #... and when the sonde wasn't moved
    filter(is.na(depth_change) | depth_change != "sonde moved") %>%
    # Get threshold for negative slope data
    filter(slope_behind < 0) %>%
    group_by(season) %>%
    summarize(f_slope_behind_01 = quantile(slope_behind, 0.01, na.rm = TRUE))

  slope_up <- df %>%
    # REMOVE DATA WE KNOW TO BE ERRONEOUS:
    add_field_flag() %>%
    fix_calibration() %>%
    add_malfunction_flag(df = ., malfunction_records = sensor_malfunction_notes) %>%
    add_spec_flag() %>%
    # Filter to remove times when the sonde was under freexzing conditions or not submerged (depth <= 0)
    add_frozen_flag() %>%
    # fix_depth_cal() %>%
    add_unsubmerged_flag() %>%
    # remove times when sonde was moved (to avoid crazy slopes)
    add_depth_shift_flag() %>%
    # Filter to only un-flagged data, and data that was collected when sonde
    # was actually deployed
    filter(is.na(flag) | sonde_employed != 1) %>%
    #... and when the sonde wasn't moved
    filter(is.na(depth_change) | depth_change != "sonde moved") %>%
    # Get threshold for positive slope data
    filter(slope_behind > 0) %>%
    group_by(season) %>%
    summarize(f_slope_behind_99 = quantile(slope_behind, 0.99, na.rm = TRUE))

  good_data_stats <- df %>%
    # REMOVE DATA WE KNOW TO BE ERRONEOUS:
    add_field_flag() %>%
    fix_calibration() %>%
    add_malfunction_flag(df = ., malfunction_records = sensor_malfunction_notes) %>%
    add_spec_flag() %>%
    # Filter to remove times when the sonde was under freexzing conditions or not submerged (depth <= 0)
    add_frozen_flag() %>%
    # fix_depth_cal() %>%
    add_unsubmerged_flag() %>%
    # remove times when sonde was moved (to avoid crazy slopes)
    add_depth_shift_flag() %>%
    # Filter to only un-flagged data, and data that was collected when sonde
    # was actually deployed
    filter(is.na(flag) | sonde_employed != 1) %>%
    #... and when the sonde wasn't moved
    filter(is.na(depth_change) | depth_change != "sonde moved") %>%
    group_by(season) %>%
    # join our slope data thresholds:
    left_join(slope_up, by = "season") %>%
    left_join(slope_down, by = "season") %>%
    # develop other thresholds across all data
    mutate(f01 = quantile(mean, 0.01, na.rm = TRUE),
           f99 = quantile(mean, 0.99, na.rm = TRUE)) %>%
           # f_slope_behind_01 = slope_down, #quantile(slope_behind, 0.01, na.rm = TRUE),
           # f_slope_behind_99 = slope_up) %>% #quantile(slope_behind, 0.99, na.rm = TRUE)) %>%
    # THEN, GET STANDARD DEVIATION OF ONLYYYY VALUES WITHIN THE 1-99th PERCENTILE OF THAT GOOD DATA:
    filter(mean > f01 & mean < f99) %>%
    # SD is the ONLY statistic that uses this winnowed-down data set in its development.
    # All else use the full, "good" data set.
    summarize(site = paste0(unique(site)),
              parameter = paste0(unique(parameter)),
              t_mean01 = as.numeric(paste0(unique(f01))),
              t_mean99 = as.numeric(paste0(unique(f99))),
              t_slope_behind_01 = as.numeric(paste0(unique(f_slope_behind_01))),
              t_slope_behind_99 = as.numeric(paste0(unique(f_slope_behind_99))),
              # This stat is useless. Should remove eventually.
              t_sd_0199 = sd(mean, na.rm = T))

  return(good_data_stats)

}
#' @title Munge API data for QAQC workflow
#'
#' @description
#' A function designed to munge the raw API data for the QAQC workflow.
#'
#' @param api_path Path where the raw API data lives.
#' @param network Options include "csu", "virridy"
#' @param require An upstream dependency (i.e., alternative to `{necessary object}`).
#'
#' @return A dataframe with the munged API data.
#'
#' @examples
# munge_api_data(api_path = "data/api/incoming/")

munge_api_data <- function(api_path, network = "csu", require = NULL) {

  api_data <- list.files(path = api_path, full.names = TRUE, pattern = "*.csv") %>%
    purrr::map_dfr(~data.table::fread(.) %>%
                     dplyr::select(-id)) %>%
    # remove overlapping API-pull data
    dplyr::distinct()

  if(network %in% c("csu", "CSU")){

  api_data <- api_data %>%
    # remove VuLink data
    dplyr::filter(!grepl("vulink", name, ignore.case = TRUE)) %>% # these will throw us errors if there is none?
    # remove Virridy data
    dplyr::filter(!grepl("virridy", name, ignore.case = TRUE)) %>%
    dplyr::select(-name) %>%
    # Convert UTC (as it is sent from HydroVU API) to MST:
    dplyr::mutate(DT = lubridate::as_datetime(timestamp, tz = "UTC")) %>%
    dplyr::mutate(DT = lubridate::with_tz(DT, tzone = "MST"),
           DT_round = lubridate::round_date(DT, "15 minutes"),
           DT_join = as.character(DT_round),
           site = tolower(site)) %>%
    # These sites will be considered the same site for this workflow
    dplyr::mutate(site = ifelse(site == "rist", "tamasag",
                  ifelse(site == "elc", "boxelder", site))) %>%
    # Lastly, we swapped Boxelder's sonde out for Rist's late in 2022:
    dplyr::mutate(site = ifelse(site == "tamasag" & DT > lubridate::ymd("2022-09-20", tz = "MST") & DT < lubridate::ymd("2023-01-01", tz = "MST"), "boxelder", site)) %>%
    dplyr::distinct(.keep_all = TRUE)
  }

  if(network %in% c("virridy", "Virridy")){

    api_data <- api_data %>%
      # remove VuLink data
      dplyr::filter(!grepl("vulink", name, ignore.case = TRUE)) %>% # these will throw us errors if there is none?
      dplyr::select(-name) %>%
      # Convert UTC (as it is sent from HydroVU API) to MST:
      dplyr::mutate(DT = lubridate::as_datetime(timestamp, tz = "UTC")) %>%
      dplyr::mutate(DT = lubridate::with_tz(DT, tzone = "MST"),
                    DT_round = lubridate::round_date(DT, "15 minutes"),
                    DT_join = as.character(DT_round),
                    site = tolower(site)) %>%
      dplyr::distinct(.keep_all = TRUE)
    }

  return(api_data)

}
#' @title Network Check
#'
#' @description
#' This function performs a network check on a given data frame, flagging potential
#' issues in the data based on upstream and downstream sites.
#'
#' @param df A site-parameter data frame that has gone through the initial flagging
#' process.
#' @param network Whether the network check is happening across the Virridy sites or CSU sites.
#' @return A modified data frame flags that have been altered based on the network check.
#'
#' @examples
#' network_check(df = all_data_flagged$`archery-Actual Conductivity`)

network_check <- function(df, network = "csu") {

  df <- df

  # get the site name from the site column in the df of interest
  site_name <- unique(na.omit(df$site))
  # get the parameter from the parameter column in the df of interest
  parameter_name <- unique(na.omit(df$parameter))

  # vector of sites in the order that they are in spatially
  # some sites have some funkiness going on

  sites_order <- c("tamasag", # rist
                   "legacy",
                   "lincoln",
                   "timberline",
                   "prospect",
                   "boxelder", # elc
                   "archery",
                   "river bluffs")

  width_fun = ifelse(site_name == "tamasag", 17, # 2 hours before/after
              ifelse(site_name == "legacy", 17,
              ifelse(site_name == "lincoln", 17,
              ifelse(site_name == "timberline", 17,
              ifelse(site_name == "prospect", 17,
              ifelse(site_name == "boxelder", 17,
              ifelse(site_name == "archery", 17,
              ifelse(site_name == "river bluffs", 17, NA))))))))

  if(network %in% c("virridy", "Virridy")){

    sites_order <-  c("joei",
                      "cbri",
                      "chd",
                      "pfal",
                      "pbd",
                      "tamasag",
                      "legacy",
                      "lincoln",
                      "timberline",
                      #"springcreek",
                      "prospect",
                      "boxelder",
                      #boxcreek,"
                      "archery",
                      "river bluffs")

    width_fun = ifelse(site_name == "joei", 17, # 2 hours before/after
                ifelse(site_name == "cbri", 17,
                ifelse(site_name == "chd", 17,
                ifelse(site_name == "pfal", 17,
                ifelse(site_name == "pbd", 17,
                ifelse(site_name == "sfm", 17,
                ifelse(site_name == "lbea", 17,
                ifelse(site_name == "penn", 17,
                ifelse(site_name == "tamasag", 17,
                ifelse(site_name == "legacy", 17,
                ifelse(site_name == "lincoln", 17,
                ifelse(site_name == "timberline", 17,
                ifelse(site_name == "timberline virridy", 17,
                ifelse(site_name == "springcreek", 17,
                ifelse(site_name == "prospect", 17,
                ifelse(site_name == "prospect virridy", 17,
                ifelse(site_name == "boxelder", 17,
                ifelse(site_name == "boxcreek", 17,
                ifelse(site_name == "archery", 17,
                ifelse(site_name == "archery virridy", 17,
                ifelse(site_name == "river bluffs", 17, NA)))))))))))))))))))))

    if(site_name %in% c("penn", "sfm", "lbea")){

      sites_order <- c("penn",
        "sfm",
        "lbea")

    }

    if(site_name == "springcreek"){

      sites_order <- c("timberline virridy",
                       "springcreek",
                       "prospect virridy")

    }

    if(site_name == "boxcreek"){

      sites_order <- c("boxelder virridy",
                       "boxcreek",
                       "archery virridy")

    }

  }

  # determining the index for the site of interest.
  site_index <- which(sites_order == sites_order[grep(gsub(" virridy", "", site_name), sites_order, ignore.case = TRUE)])

  # Generating df name to pull from df_list list
  site_param <- paste0(site_name, "-", parameter_name)

  prev_site_df <- tibble(DT_round = NA)
  next_site_df <- tibble(DT_round = NA)

  tryCatch({
    previous_site <- paste0(sites_order[site_index-1],"-",parameter_name)
    prev_site_df <- intersensor_checks[[previous_site]] %>%
      select(DT_round, site_up = site, flag_up = flag) %>%
      data.table()},
    error = function(err) {
      cat("No upstream site.\n")})

  tryCatch({
    next_site <- paste0(sites_order[site_index+1],"-",parameter_name)
    next_site_df <- intersensor_checks[[next_site]] %>%
      select(DT_round, site_down = site, flag_down = flag) %>%
      data.table()},
    error = function(err) {
      cat("No downstream site.\n")})


  join <- df %>%
    left_join(., prev_site_df, by = "DT_round") %>%
    left_join(., next_site_df, by = "DT_round")

  if(!("flag_down" %in% colnames(join))) {join$flag_down <- NA}
  if(!("flag_up" %in% colnames(join))) {join$flag_up <- NA}
  if(!("site_down" %in% colnames(join))) {join$site_down <- NA}
  if(!("site_up" %in% colnames(join))) {join$site_up <- NA}


  # Define a function to check if a given 2-hour window has any instances of the same word
  check_2_hour_window_fail <- function(x) {
    sum(x) >= 1
  }

  df_test <- join %>%
    # No upstream/downstream flag = 0
    mutate(flag_binary = ifelse(
      (is.na(flag_up) | grepl("drift|DO interference|repeat|sonde not employed|frozen|unsubmerged|missing data|site visit|sv window|sensor malfunction|sonde burial|sensor biofouling|improper level cal|sonde moved", flag_up)) &
        (is.na(flag_down) | grepl("drift|DO interference|repeat|sonde not employed|frozen|unsubmerged|missing data|site visit|sv window|sensor malfunction|sonde burial|sensor biofouling|improper level cal|sonde moved", flag_down)), 0, 1)) %>%
    mutate(overlapping_flag = zoo::rollapply(flag_binary, width = width_fun, FUN = check_2_hour_window_fail, fill = NA, align = "center")) %>%
    # If there is a flag (flags associated with spikes in concentration or funkiness like that), and there is also a flag up/downstream at the same time (2 hour window) it is likely a real
    # WQ event and should therefore not be considered "poor" data:
    mutate(auto_cleaned_flag = ifelse(!is.na(flag) & !grepl("drift|DO interference|repeat|sonde not employed|frozen|unsubmerged|missing data|site visit|sv window|sensor malfunction|sonde burial|sensor biofouling|improper level cal|sonde moved", flag) & overlapping_flag == TRUE, NA, flag)) %>%
    select(-c(flag_up, flag_down, site_up, site_down, flag_binary, overlapping_flag))

  # df_test <- join %>%
  #   mutate(flag_binary = ifelse(#grepl("slope|suspect", flag) &
  #     (is.na(flag_up) | grepl("seasonal range|repeat|sonde not employed|missing data|site visit|sv window", flag_up)) &
  #       (is.na(flag_down) | grepl("seasonal range|repeat|sonde not employed|missing data|site visit|sv window", flag_down)), 0, 1)) %>%
  #   #arrange(timestamp) %>%
  #   mutate(overlapping_flag = zoo::rollapply(flag_binary, width = width_fun, FUN = check_2_hour_window_fail, fill = NA, align = "center")) %>%
  #   mutate(cleaner_flag = ifelse(!is.na(flag) & !grepl("seasonal range|repeat|sonde not employed|missing data|site visit|sv window", flag) & overlapping_flag == TRUE, NA, flag)) %>%
  #   select(-c(flag_up, flag_down, site_up, site_down, flag_binary, overlapping_flag))

  return(df_test)

}
# Function to check for package installation, then install (if necessary) and load libraries.
# Adapted from code developed by Caitlin Mothes, PhD.


# fill in with packages that need to be loaded:
# packages <- c('tidyverse',
#               'sf')

package_loader <- function(x) {
  if (x %in% installed.packages()) {
    library(x, character.only = TRUE)
  } else {
    install.packages(x)
    library(x, character.only = TRUE)
  }
}

# lapply(packages, package_loader)
# generate_flag_report <- function(df) {
#
#   # Extract the site and parameter from the df_name
#   site <- unique(na.omit(df$site))
#   parameter <- unique(na.omit(df$parameter))
#
#   list_of_flags <- c("sonde not employed", # add_field_flag()
#                      "site visit", # add_field_flag()
#                      "sv window", # add_field_flag()
#                      "sensor malfunction", # add_malfunction_flag()
#                      "outside of sensor specification range", # add_spec_flag()
#                      "outside of seasonal range", # add_seasonal_flag()
#                      "slope violation", # add_seasonal_flag()
#                      "outside sd range", # add_seasonal_flag()
#                      "repeated value", # add_repeat_flag()
#                      "missing data", # add_na_flag()
#                      "suspect data") # add_suspect_flag()
#
#   # check these
#   sans_na_flags <- "^(missing data|
#                       sonde not employed;\\nmissing data|
#                       missing data;\\nsuspect data|
#                       sonde not employed;\\nmissing data;\\nsuspect data|
#                       site visit;\\nmissing data;\\nsuspect data|
#                       sv window;\\nmissing data;\\nsuspect data)$"
#
#   # summarize total data points
#   total_observations <- df %>%
#     summarise(n_total = n_distinct(DT_round)) %>%
#     pull(n_total)
#
#   # summarize total data points sans missing data
#   total_observations_1 <- df %>%
#     # filter out when flag has only missing data or only sonde not employed and missing data
#     filter(!str_detect(flag, sans_na_flags)) %>%
#     summarise(n_total = n_distinct(DT_round)) %>%
#     pull(n_total)
#
#   # summarize total days
#   total_observations_dates <- df %>%
#     group_by(date = format(DT_round, "%m-%d-%Y")) %>%
#     summarize(n_total = nrow(date)) %>%
#     nrow()
#
#   # summarize total days sans missing data
#   total_observations_dates_1 <- df %>%
#     filter(!str_detect(flag, sans_na_flags)) %>%
#     group_by(date = format(DT_round, "%m-%d-%Y")) %>%
#     summarize(n_total = nrow(date)) %>%
#     nrow()
#
#   row_list <- list()
#   for (i in list_of_flags) {
#
#     # summarize flagged data points
#     flagged_observations <- df %>%
#       filter(str_detect(flag, i)) %>%
#       summarise(n_flag = n_distinct(DT_round)) %>%
#       pull(n_flag)
#     # summarize flagged data points
#     flagged_observations_1 <- df %>%
#       filter(!str_detect(flag, sans_na_flags)) %>%
#       filter(str_detect(flag, i)) %>%
#       summarise(n_flag = n_distinct(DT_round)) %>%
#       pull(n_flag)
#     # summarize percent data points that are flagged
#     percent_flagged <- flagged_observations/total_observations
#     # summarize percent data points that are flagged sans missing data
#     percent_flagged_1 <- flagged_observations_1/total_observations_1
#
#     # summarize flagged days
#     flagged_observations_dates <- df %>%
#       filter(str_detect(flag, i)) %>%
#       group_by(date = format(DT_round, "%m-%d-%Y")) %>%
#       summarize(n_total = nrow(date)) %>%
#       nrow()
#     # summarize flagged days
#     flagged_observations_dates_1 <- df %>%
#       filter(!str_detect(flag, sans_na_flags)) %>%
#       filter(str_detect(flag, i)) %>%
#       group_by(date = format(DT_round, "%m-%d-%Y")) %>%
#       summarize(n_total = nrow(date)) %>%
#       nrow()
#     # summarize percent days that are flagged
#     percent_flagged_dates <- flagged_observations_dates/total_observations_dates
#     # summarize percent days that are flagged
#     percent_flagged_dates_1 <- flagged_observations_dates_1/total_observations_dates_1
#
#     # creating a row with the information
#     calculated_values <- tibble(
#       # metadata
#       site = site,
#       parameter = parameter,
#       flag = i,
#       # data points
#       data_points_flagged_percentage = format(round(percent_flagged * 100, 2), nsmall = 2),
#       data_points_flagged = flagged_observations,
#       total_data_points = total_observations,
#       data_points_flagged_percentage_sans_na = format(round(percent_flagged_1 * 100, 2), nsmall = 2),
#       data_points_flagged_sans_na = flagged_observations_1,
#       total_data_points_sans_na = total_observations_1,
#       # dates
#       dates_flagged_percentage = format(round(percent_flagged_dates * 100, 2), nsmall=2),
#       dates_flagged = flagged_observations_dates,
#       total_dates = total_observations_dates,
#       dates_flagged_percentage_sans_na = format(round(percent_flagged_dates_1 * 100, 2), nsmall=2),
#       dates_flagged_sans_na = flagged_observations_dates_1,
#       total_dates_sans_na = total_observations_dates_1
#       )
#
#     row_list[[i]] <- calculated_values
#   }
#
#   #calculated_df <- bind_cols(row_list)
#   return(bind_rows(row_list))
#
# }

photo_plotter <- function(index = 200, output_folder){

  param_1_label <- paste0(parameters[1], " (", param_unit[1],")")
  param_2_label <- paste0(parameters[2], " (", param_unit[2],")")

# get the adjustment, breaks etc so axes look nice
  param_1_max <- max(wq_tl[[parameters[1]]], na.rm = T)
  param_1_max_int <- as.integer(max(wq_tl[[parameters[1]]], na.rm = T))+1

  adjustment <- max(wq_tl[[parameters[2]]], na.rm = T) / param_1_max

  brk <- ifelse(param_1_max_int < 6, 1,2)

# Bounds for the y-axis

  lower_bound <- 0
  upper_bound <-  max(wq_tl[[parameters[2]]], na.rm = T)


  #This is the index of the image for the background for the individual photo
  simul = wq_tl[index,]
  # this is all the data before the image (ie previous photos from the timelapse)
  upto = wq_tl[1:index,]
  #read the image for the background
  photo_bg <- readJPEG(simul$filename)
  #create an individual image
  back <- ggplot() +
    #plot the photo
    annotation_custom(rasterGrob(photo_bg,
                                 width = unit(1,"npc"),
                                 height = unit(1,"npc")),
                      -Inf, Inf, -Inf, Inf)

  #plot the data for this image (includes all the preivous data)
  inset <- ggplot() +
    geom_ribbon(data = upto, aes(x = DT_round,
#to do: need to figure out adjustment factors for each parameter to plot correctly
                                 y = .data[[parameters[1]]] *adjustment,
                                 ymin = 0,
                                 ymax = .data[[parameters[1]]] * adjustment),
                color = "white",
                fill = "white",
                #linetype = "dash",
                alpha = 0.75)+
    geom_path(data = upto, aes(x = DT_round, y = .data[[parameters[2]]]),
              color = "#F34646", size=2) +
    geom_point(data = simul, aes(x = DT_round, y = .data[[parameters[2]]]),
               color = "#F34646")+
    # ylim(min(wq_tl[[parameters[2]]], na.rm = T),
    #      max(wq_tl[[parameters[2]]], na.rm = T))+
    xlim(min(wq_tl$DT_round, na.rm = T),
         max(wq_tl$DT_round, na.rm = T))+
    scale_y_continuous(name = param_2_label, limits = c(0, upper_bound),
                       sec.axis = sec_axis(trans = ~./adjustment ,
                                           name = param_1_label,
                                           breaks = seq(0,param_1_max_int,brk)))+
    dark_theme_light(base_size = 10) +
    theme(axis.title.y.right = element_text(color="white"),axis.text.y.right = element_text(color = "white"), axis.title.y.left = element_text(color="#F34646"), axis.text.y.left = element_text(color = "#F34646")) +
    theme(
      panel.background = element_rect(fill = "transparent"),
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.background = element_rect(fill = "transparent"),
      legend.box.background = element_rect(fill = "transparent"),
      panel.border = element_blank(),
      axis.line = element_line(color = 'gray80'),
    ) +
    ylab(param_2_label) +
    xlab('')

  gp1 <- back +
    inset_element(inset,
                  left = 0.16,
                  bottom = 0.15,
                  right = 0.9,
                  top = 0.6)
  #print(gp1)
  ggsave(paste0(folder_path, index, ".png"),
         width = 1920,
         height = 1080,
         units = "px")
  #dev.copy(gp1,paste0('data/timelapse_photos/vid_image/', max(upto$rowid), ".png"))
  #dev.off()
}

#test functions
#photo_plotter(2, "data/timelapse_photos/sjs_test/")
#map(1:nrow(wq_tl), ~photo_plotter(.x, "data/timelapse_photos/sjs_test/"))
plotter <- function(site_name = "boxelder",
                    start_date = "2023-03-04",
                    end_date = "2023-07-07",
                    title = "wack data title",
                    parameter = "all") {

  # Plotting options: all, Specific Conductivity, Turbidity, Chl-a Fluorescence, ORP, Depth, pH, DO

  sites_order <- c("tamasag", # rist
                   "legacy",
                   "lincoln",
                   "timberline",
                   "prospect",
                   "boxelder", # elc
                   "archery",
                   "river bluffs")

  # determining the index for the site of interest.
  site_index <- which(sites_order == site_name)
  previous_site <- sites_order[site_index-1]
  next_site <- sites_order[site_index+1]

  # Pull in WQ data from Manners Bridge (in the canyon):
  contrail <- list.files("data/context_data/contrail", full.names = TRUE) %>%
    map_dfr(~fread(.)) %>%
    mutate(date = as_date(ymd_hms(Reading)),
           Reading = floor_date(ymd_hms(Reading), "hour")) %>%
    filter(date >= start_date & date <= end_date)

  # Pull in a list of precip monitoring sites in each sonde's watershed:
  rain_list <- fread("data/context_data/site_link.csv") %>%
    filter(site == site_name) %>% pull(`Sensor Name`)

  # Pull in raw precip data
  rain_data <- fread("data/context_data/fc_rain.csv") %>%
    filter(`Sensor Name` %in% rain_list) %>%
    mutate(date = as_date(mdy_hms(Timestamp)),
           DT = floor_date(mdy_hms(Timestamp), unit = "hour")) %>%
    filter(date >= start_date & date <= end_date) %>%
    group_by(DT) %>%
    summarize(total_ws = sum(`Incremental Rainfall (in)`, na.rm = TRUE))

  # Grab the largest amount of rain across the sonde's watershed at an
  # hourly timestep:
  rain_stat <- fread("data/context_data/fc_rain.csv") %>%
    filter(`Sensor Name` %in% rain_list) %>%
    mutate(date = as_date(mdy_hms(Timestamp)),
           DT = floor_date(mdy_hms(Timestamp), unit = "hour")) %>%
    group_by(DT) %>%
    summarize(total_ws = sum(`Incremental Rainfall (in)`, na.rm = TRUE)) %>%
    filter(total_ws == max(total_ws))

  # Pull in flagged PSN WQ data:
  df <- fread('data/flagged/all_data_flagged_plotter.csv') %>%
    filter(!grepl("sonde not employed", final_cleaned_flag)) %>%
    mutate(date = as_date(DT_round)) %>%
    # filter to dates of interest:
    dplyr::filter((date) >= as_date(start_date) & (date) <= as_date(end_date))

  # site of interest data:
  site_df <- df %>%
    filter(!is.na(mean)) %>%
    filter(site == site_name)

  # downstream of site of interest (if there is one):
  downstream_df <- site_df
  try(downstream_df <- df %>%
        filter(!grepl("site visit|sv window", final_cleaned_flag)) %>%
        filter(site == next_site))

  # upstream of site of interest (if there is one):
  upstream_df <- site_df
  try(upstream_df <- df %>%
        filter(!grepl("site visit|sv window", final_cleaned_flag)) %>%
        filter(site == previous_site))

# All the plotz
  depth <- ggplot() +
    geom_line(data = filter(contrail, Unit == "ft"), aes(Reading, Value * 0.3048), color = "lightblue") +
    geom_line(data = filter(downstream_df, parameter == "Depth"), aes(DT_round, mean), cex = 0.2, alpha = 1, color = "lightgrey") +
    geom_line(data = filter(upstream_df, parameter == "Depth"), aes(DT_round, mean), cex = 0.2, alpha = 1, color = "darkgrey") +
    geom_line(data = filter(site_df, parameter == "Depth"), aes(DT_round, mean), color = "black") +
    geom_point(data = filter(site_df, parameter == "Depth" & grepl("site visit", final_cleaned_flag)), aes(DT_round, mean), color = "#E70870", cex = 1.5) +
    geom_point(data = filter(site_df, parameter == "Depth" & grepl("sonde moved", final_cleaned_flag)), aes(DT_round, mean), color = "#FFCA3A", cex = 0.75) +
    theme_bw() +
    ylim(min = min(filter(site_df, parameter == "Depth")$mean), max = max(filter(site_df, parameter == "Depth")$mean)) +
    ylab("Depth m") + xlab("") + ggtitle(title)

  temp <- ggplot() +
    geom_line(data = filter(contrail, Unit == "C"), aes(Reading, Value), color = "lightblue") +
    geom_line(data = filter(downstream_df, parameter == "Temperature"), aes(DT_round, mean), cex = 0.2, alpha = 1, color = "lightgrey") +
    geom_line(data = filter(upstream_df, parameter == "Temperature"), aes(DT_round, mean), cex = 0.2, alpha = 1, color = "darkgrey") +
    geom_line(data = filter(site_df, parameter == "Temperature"), aes(DT_round, mean), color = "black") +
    geom_point(data = filter(site_df, parameter == "Temperature" & grepl("site visit|sv window", final_cleaned_flag)), aes(DT_round, mean), color = "#E70870") +
    theme_bw() +
    ylim(min = min(filter(site_df, parameter == "Temperature")$mean), max = max(filter(site_df, parameter == "Temperature")$mean)) +
    ylab("Temperature C") + xlab("") + ggtitle("")

  ph <- ggplot() +
    geom_line(data = filter(contrail, Unit == "pH"), aes(Reading, Value), color = "lightblue") +
    geom_line(data = filter(downstream_df, parameter == "pH"), aes(DT_round, mean), cex = 0.2, alpha = 1, color = "lightgrey") +
    geom_line(data = filter(upstream_df, parameter == "pH"), aes(DT_round, mean), cex = 0.2, alpha = 1, color = "darkgrey") +
    geom_line(data = filter(site_df, parameter == "pH"), aes(DT_round, mean), color = "black") +
    geom_point(data = filter(site_df, parameter == "pH" & grepl("site visit|sv window", final_cleaned_flag)), aes(DT_round, mean), color = "#E70870") +
    theme_bw() +
    ylim(min = min(filter(site_df, parameter == "pH")$mean), max = max(filter(site_df, parameter == "pH")$mean)) +
    xlab("") + ylab("pH") + ggtitle("")

  orp <- ggplot() +
    geom_line(data = filter(downstream_df, parameter == "ORP"), aes(DT_round, mean), cex = 0.2, alpha = 1, color = "lightgrey") +
    geom_line(data = filter(upstream_df, parameter == "ORP"), aes(DT_round, mean), cex = 0.2, alpha = 1, color = "darkgrey") +
    geom_line(data = filter(site_df, parameter == "ORP"), aes(DT_round, mean), color = "black") +
    geom_point(data = filter(site_df, parameter == "ORP" & grepl("site visit|sv window", final_cleaned_flag)), aes(DT_round, mean), color = "#E70870") +
    theme_bw() +
    ylim(min = min(filter(site_df, parameter == "ORP")$mean), max = max(filter(site_df, parameter == "ORP")$mean)) +
    ylab("ORP") + xlab("") + ggtitle("")

  spc <- ggplot() +
    geom_line(data = filter(contrail, Unit == "uS/cm"), aes(Reading, Value), color = "lightblue") +
    geom_line(data = filter(downstream_df, parameter == "Specific Conductivity"), aes(DT_round, mean), cex = 0.2, alpha = 1, cex = 0.2, alpha = 1, color = "lightgrey") +
    geom_line(data = filter(upstream_df, parameter == "Specific Conductivity"), aes(DT_round, mean), cex = 0.2, alpha = 1, cex = 0.2, alpha = 1, color = "darkgrey") +
    geom_line(data = filter(site_df, parameter == "Specific Conductivity"), aes(DT_round, mean), color = "black", cex = 0.8) +
    geom_point(data = filter(site_df, parameter == "Specific Conductivity" & grepl("site visit|sv window", final_cleaned_flag)), aes(DT_round, mean), color = "#E70870", cex = 1.5) +
    theme_bw() +
    ylim(min = min(filter(site_df, parameter == "Specific Conductivity")$mean) - 1, max = max(filter(site_df, parameter == "Specific Conductivity")$mean) + 1) +
    ylab("SpC uS/cm") + xlab("") + ggtitle("")

  do <- ggplot() +
    geom_line(data = filter(contrail, Unit == "mg/L"), aes(Reading, Value), color = "lightblue") +
    geom_line(data = filter(downstream_df, parameter == "DO"), aes(DT_round, mean), cex = 0.2, alpha = 1, color = "lightgrey") +
    geom_line(data = filter(upstream_df, parameter == "DO"), aes(DT_round, mean), cex = 0.2, alpha = 1, color = "darkgrey") +
    geom_line(data = filter(site_df, parameter == "DO"), aes(DT_round, mean), color = "black") +
    geom_point(data = filter(site_df, parameter == "DO" & grepl("site visit|sv window", final_cleaned_flag)), aes(DT_round, mean), color = "#E70870") +
    theme_bw() +
    ylim(min = min(filter(site_df, parameter == "DO")$mean) - 1, max = max(filter(site_df, parameter == "DO")$mean) + 1) +
    ylab("DO mg/L") + xlab("") + ggtitle("")

  # Handle plotting for sondes with turbidity vs. chl-a data:

  turb <- ggplot() +
    geom_line(data = filter(contrail, Unit == "ntu"), aes(Reading, Value), color = "lightblue") +
    geom_line(data = filter(downstream_df, parameter == "Turbidity"), aes(DT_round, mean), cex = 0.2, alpha = 1, color = "lightgrey") +
    geom_line(data = filter(upstream_df, parameter == "Turbidity"), aes(DT_round, mean), cex = 0.2, alpha = 1, color = "darkgrey") +
    geom_line(data = filter(site_df, parameter == "Turbidity"), aes(DT_round, mean), color = "black") +
    geom_point(data = filter(site_df, parameter == "Turbidity" & grepl("site visit|sv window", final_cleaned_flag)), aes(DT_round, mean), color = "#E70870") +
    theme_bw() +
    ylim(min = min(filter(site_df, parameter == "Turbidity")$mean) + 10, max = max(filter(site_df, parameter == "Turbidity")$mean) + 10) +
    ylab("Turbidity NTU") + xlab("") + ggtitle("")

  if(nrow(filter(site_df, parameter == "Chl-a Fluorescence")) > 0){
    chla <- ggplot() +
      geom_line(data = filter(downstream_df, parameter == "Chl-a Fluorescence"), aes(DT_round, mean), cex = 0.2, alpha = 1, color = "lightgrey") +
      geom_line(data = filter(upstream_df, parameter == "Chl-a Fluorescence"), aes(DT_round, mean), cex = 0.2, alpha = 1, color = "darkgrey") +
      geom_line(data = filter(site_df, parameter == "Chl-a Fluorescence"), aes(DT_round, mean), color = "black") +
      geom_point(data = filter(site_df, parameter == "Chl-a Fluorescence" & grepl("site visit|sv window", final_cleaned_flag)), aes(DT_round, mean), color = "#E70870") +
      theme_bw() +
      ylab("Chl-a RFU") + xlab("") + ggtitle("")
  }

  if(nrow(filter(site_df, parameter == "Chl-a Fluorescence")) == 0){ chla <- NULL }

  # rain proxy plot. this just plots the sum of at an hourly timestep of
  # rain at all sites in the watershed. not rooted in anything scientific,
  # but just shows in a crude way whether or not there was rain happening
  # at any given time (and at any given point) during the time of interest.
  rain <- ggplot(data = rain_data) +
    geom_col(aes(x = DT,
                 y = total_ws), color = "#002EA3") +
    ylab("Rain (PROXY)") + xlab("") +
    ylim(min = 0, max = rain_stat$total_ws) +
    theme_bw()

  # handling what gets plotted based on user input:

  if(parameter %in% c("All", "all", "all parameters", "all params")){
    if(is.null(chla)) {

      plot <- ggarrange(depth, temp, ph, orp, spc, do, turb,
                         rain, nrow = 4, ncol = 2)
    }

    if(!is.null(chla)){
      plot <- ggarrange(depth, temp, ph, orp, spc, do, chla,
                         rain, nrow = 4, ncol = 2)
    }
  }

  if(parameter %in% c("pH", "ph")){
    plot <- ggarrange(depth, temp, ph, rain,
                      ncol = 1, nrow = 4)
  }

  if(parameter %in% c("orp", "ORP")){
    plot <- ggarrange(depth, temp, orp, rain,
                       ncol = 1, nrow = 4)
  }

  if(parameter %in% c("do","DO")){
    plot <- ggarrange(depth, temp, do, rain,
                      ncol = 1, nrow = 4)
  }

  if(parameter %in% c("turbidity", "turb", "Turb", "Turbidity")){
    plot <- ggarrange(depth, temp, turb, rain,
                      ncol = 1, nrow = 4)
  }

  if(parameter %in% c("Specific Conductivity", "conductivity", "specific conductivity", "sc", "spc", "SC","SpC")){
    plot <- ggarrange(depth, temp, spc, rain,
                      ncol = 1, nrow = 4)
  }

  if(parameter %in% c("chla", "Chla", "chl-a", "chlorophyll a", "Chl-a Fluorescence")){
    plot <- ggarrange(depth, temp, chla, rain,
                      ncol = 1, nrow = 4)
  }

  return(plot)
}
retrieve_relevant_data_name <- function(df_name_arg, year_week_arg = NULL, interval_arg = NULL) {
  if(is.null(interval_arg)){
    if (df_name_arg %in% names(verified_data) & any(year_week_arg %in% verified_data[[df_name_arg]]$y_w)) {
      return("verified_data")
    }
    if (df_name_arg %in% names(intermediary_data) & any(year_week_arg %in% intermediary_data[[df_name_arg]]$y_w)) {
      return("intermediary_data")
    }
    if (df_name_arg %in% names(all_data) & any(year_week_arg %in% all_data[[df_name_arg]]$y_w)) {
      return("all_data")
    }
  }
  if(is.null(year_week_arg)){
    if (df_name_arg %in% names(verified_data) & any(verified_data[[df_name_arg]]$DT_round %within% interval_arg)) {
      return("verified_data")
    }
    if (df_name_arg %in% names(intermediary_data) & any(intermediary_data[[df_name_arg]]$DT_round %within% interval_arg)) {
      return("intermediary_data")
    }
    if (df_name_arg %in% names(all_data) & any(all_data[[df_name_arg]]$DT_round %within% interval_arg)) {
      return("all_data")
    }
  }
}
## Water Sampling Data:

# Goal:
#Save data in the correct format for RMRS spreadsheet
#Save all water sampling probe values in a spreadsheet
# To get the RMRS style data for a specfic date of sampling,
# Input the date of interest in sampling_spreadsheet_creator
#sampling_spreadsheet_creator(date_oi = "2023-11-17")

# To get all the water sampling data and save to CSV in sampling notes
# This also returns the df sampling_notes in case you want to review in R
#sampling_spreadsheet_creator(all_dates = TRUE)


sampling_spreadsheet_creator <- function(date_oi = "2023-10-16", all_dates = FALSE ){

  #source clean mwater script for all notes cleaned

  source("src/load_mWater_notes.R")

  #pull in site meta data
  site_meta <- read_csv("data/metadata/water_sampling_sites.csv",show_col_types = FALSE)%>%
    select(site = site_code, Site_Name, site_label_rmrs)
  # sort for sites in upper network (ie. acronyms rather than street names)
  upper_sites <- read_csv("data/metadata/water_sampling_sites.csv",show_col_types = FALSE)%>%
    filter(watershed != "CLP  Mainstem-Fort Collins")%>%
    #this is to help match with user input
    mutate(site_code = tolower(site_code))

  # create df of all water samples and save DT, handheld probe and chla volume data
  sampling_notes <- load_mWater_notes()%>%
    filter(grepl("Sampling",visit_type))%>%
    mutate(all_pics_taken = case_when(!is.na(downstream_pic)&!is.na(upstream_pic)&!is.na(clarity)&!is.na(filter_pic) ~ TRUE, TRUE ~ FALSE),
           #correct names if it is in our upper sites (acronyms)
           site = ifelse(site %in% upper_sites$site_code, toupper(site), site),
           DT_round = round_date(start_DT, "15 minutes"))%>%
    select(site,crew, DT_round, date, time = start_time_mst, sample_collected, chla_volume_ml, vol_filtered_blank_dup, do_mgl, cond_ms_cm, temp_c, visit_comments, all_pics_taken, q_cfs)

  # Distinguish BLANK and DUPLICATE values
  blanks_dups <- sampling_notes %>%
    #find all values that have blank or dup
    filter(grepl("DUPLICATE|BLANK", sample_collected)) %>%
    # change sample collected to match BLANK/DUP
    mutate(sample_collected = ifelse(grepl("DUPLICATE", sample_collected), "DUPLICATE", "BLANK"),
           # Volume filtered blank dup becomes chla volume
           chla_volume_ml = vol_filtered_blank_dup,
           #drop vol_filtered_blank/dup
           vol_filtered_blank_dup = NULL)

  # Add blank and duplicate values back to main
  sampling_notes <- sampling_notes%>%
    #get rid of blank/dup in sample collcected
    mutate(sample_collected = gsub("DUPLICATE|BLANK| |,", "", sample_collected),
           #drop vol_filtered_blank/dup
           vol_filtered_blank_dup = NULL)%>%
    #bring in blank and dup rows
    rbind(blanks_dups)%>%
    #arrange by datetime and site (Blanks and dups go second)
    arrange(DT_round, site)%>%
    # join with RMRS friendly metadata
    left_join(site_meta, by = "site")



 if (all_dates == TRUE) {

   sampling_notes_output <- sampling_notes%>%
     # select only the needed columns, saved in the correct order and fix column names
     select(site_code = site, Date = date, SampleType = sample_collected, time_mst = time,chla_volume_ml,  do_mgl, cond_ms_cm, temp_c, visit_comments)

   # write to csv
   write_csv(x = sampling_notes_output, file = paste0("data/sampling_notes/all_samples_as_of_",as.character(Sys.Date()),".csv" ))
 }else{
   #grab desired date
   date_oi_clean <- as.Date(date_oi, tz = "America/Denver")
   # filter sampling notes df by desired date
   samples_of_day <- filter(sampling_notes,date == date_oi_clean )%>%
     #match date to RMRS style
     mutate(Date = format(date, "%d-%b-%y"),
            #create SiteDescr column for RMRS sheet
            SiteDescr = paste0(site, "_",format(date, "%m%d%y")))%>%
     # select only the needed columns, saved in the correct order and fix column names
     select(Site = site ,SiteDescr, SiteLabel = site_label_rmrs , Date, SampleType = sample_collected, q_cfs,
            time, do_mgl, cond_ms_cm, temp_c, notes = visit_comments  )
   # write to csv
   write_csv(x = samples_of_day, file = paste0("data/sampling_notes/samples_of_",date_oi,".csv" ))
 }
}

save_intermediary_dir <- function(updated_df, df_name) {
  # if updated_site_param_df is from all_data, but there is already intermediary data for it, refuse to move on
  if (all(!updated_df$is_verified) & df_name %in% list.files(intermediary_path)) {
    while (TRUE) {
      cat("!!!YOU ARE ABOUT TO OVERWRITE INTERMEDIARY DATA WITH PRE-VERIFICATION DATA!!!")
      user_input <- readline(prompt = "Are you sure you want to continue with this action? (yes/no):")
      user_input <- tolower(user_input)

      if (user_input == "yes") {
        saveRDS(updated_df, paste0(intermediary_path, df_name))
        break
      }

      if (user_input == "no") {
        working_data <<- set_names(map(list.files(intermediary_path, full.names = TRUE), readRDS), list.files(intermediary_path))
        updated_site_param_df <<- working_data[[site_param_name]] # working dir df
        saveRDS(updated_site_param_df, paste0(intermediary_path, df_name))
        break
      }

      cat("Invalid input. Please enter one of the options from the following:\n\n")
      cat("yes\n")
      cat("no\n")
    }
  }

  # if updated_site_param_df is from verified data, refuse to move on
  if (df_name %in% list.files(verified_path)) { # test to see if this works
    while (TRUE) {
      cat("!!!YOU ARE ABOUT TO WORK ON DATA THAT HAS BEEN FINALIZED!!!")
      user_input <- readline(prompt = "Are you sure you want to continue with this action? (yes/no):")
      user_input <- tolower(user_input)

      if (user_input == "yes") {
        saveRDS(updated_df, paste0(intermediary_path, df_name))
        break
      }

      if (user_input == "no") {
        cat("Please update your `site` and `parameter` objects to data that has not been finalized.")
        break
      }

      cat("Invalid input. Please enter one of the options from the following:\n\n")
      cat("yes\n")
      cat("no\n")
    }
  }

  saveRDS(updated_df, paste0(intermediary_path, df_name))
}
# Generate plots with both weekly and daily flagged data.
# This function will generate a list of plots with both weekly and daily flagged data.
# @param site_arg A string of the site name.
# @param parameter_arg A string of the parameter name.
# @param flag_arg A string of the flag name.
# @return A list of plots with both weekly and daily flagged data.
# @examples
# stack_flag_plots(site_arg = "archery", parameter_arg = "Actual Conductivity", flag_arg = "outside of Actual Conductivity sensor specification range")
# stack_flag_plots(site_arg = "boxelder", parameter_arg = "Temperature", flag_arg = "outside of Temperature sensor specification range")

# stack_flag_plots <- function(site_arg, parameter_arg, flag_arg, df_list) {
#   # Call on the weekly and daily functions and fill their args with this
#   # functions args
#   weekly_plot_list <- generate_weekly_flag_plots(site_arg = site_arg, parameter_arg = parameter_arg, flag_arg = flag_arg, df_list = df_list)
#   daily_plot_list <- generate_daily_flag_plots(site_arg = site_arg, parameter_arg = parameter_arg, flag_arg = flag_arg, df_list = df_list)
#
#   # These two functions should always return the same amount of plots, so we can
#   # use map2() to stack them with combine_plots()
#   weekly_daily_plots <- map2(.x = weekly_plot_list, .y = daily_plot_list, ~ggarrange(.x, .y, nrow = 2, ncol = 1, heights = 5, widths = 12))
#   return(weekly_daily_plots)
# }
#' @title Summarize site parameter data from the API and field notes data frames.
#' @description
#' A short description...
#' @param site_arg A site name.
#' @param parameter_arg A parameter name.
#' @param api_data A dataframe with the munged API data.
#' @param notes The munged field notes
#' @return A dataframe with summary statistics for a given site parameter data frame.
#' @examples
# summarize_site_param(site_arg = "archery", parameter_arg = "Actual Conductivity", api_data = incoming_data_collated_csvs)
# summarize_site_param(site_arg = "boxelder", parameter_arg = "Temperature", api_data = incoming_data_collated_csvs)

summarize_site_param_full <- function(site_arg, parameter_arg, api_data, notes = field_notes) {

  # filter deployment records for the full join
  site_field_notes <- notes %>%
    filter(grepl(paste(unlist(str_split(site_arg, " ")), collapse = "|"), site, ignore.case = TRUE))

  # filtering the data and generating results
  summary <- tryCatch({
    api_data %>%
      # subset to single site-parameter combo:
      dplyr::filter(site == site_arg & parameter == parameter_arg) %>%
      # safety step of removing any erroneous dupes
      dplyr::distinct() %>%
      # across each 15 timestep, get the average value, spread, and count of obs
      dplyr::group_by(DT_round, site, parameter) %>%
      dplyr::summarize(mean = as.numeric(mean(value, na.rm = T)),
                       diff = abs(min(value, na.rm = T) - max(value, na.rm = T)),
                       n_obs = n()) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(DT_round) %>%
      # pad the dataset so that all 15-min time stamps are present
      padr::pad(by = "DT_round", interval = "15 min") %>%
      # add a DT_join column to join field notes to (make DT_round character string, so no
      # funky DT issues occur during the join):
      dplyr::mutate(DT_join = as.character(DT_round),
                    site = site_arg,
                    parameter = parameter_arg,
                    flag = NA) %>% # add "flag" column for future processing
      # join our tidied data frame with our field notes data:
      dplyr::left_join(dplyr::filter(dplyr::select(site_field_notes, sonde_employed, last_site_visit, DT_join, visit_comments, sensor_malfunction, cals_performed)),
                       by = c('DT_join')) %>%
      # make sure DT_join is still correct:
      dplyr::mutate(DT_round = lubridate::as_datetime(DT_join, tz = "MST")) %>%
      # Use fill() to determine when sonde was in the field, and when the last site visit was
      # Necessary step for FULL dataset only (this step occurs in combine_hist_inc_data.R for auto)
      tidyr::fill(c(sonde_employed, last_site_visit, sensor_malfunction)) %>%
      dplyr::mutate(sonde_employed = ifelse(is.na(sonde_employed), 0, sonde_employed))
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
#' @title Summarize site parameter data from the HydroVu API and field notes data frames.
#'
#' @description
#' A function that summarizes site parameter data from the API and field notes data frames.
#'
#' @param site_arg A site name.
#'
#' @param parameter_arg A parameter name.
#'
#' @param api_data A dataframe with the munged API data.
#'
#' @param notes The munged field notes
#'
#' @param require An upstream dependency (i.e., alternative to `{necessary object}`).
#'
#' @return A dataframe with summary statistics and field notes for a given site parameter data frame.
#'
#' @examples
# summarize_site_param(site_arg = "archery", parameter_arg = "Actual Conductivity", api_data = incoming_data_collated_csvs)
# summarize_site_param(site_arg = "boxelder", parameter_arg = "Temperature", api_data = incoming_data_collated_csvs)

summarize_site_param <- function(site_arg, parameter_arg, api_data, notes, require = NULL) {

  # filter deployment records for the full join
  site_field_notes <- notes %>%
    dplyr::filter(site == site_arg)

  # filtering the data and generating results
  summary <- tryCatch({
    api_data %>%
      # summarize the data
      dplyr::filter(site == site_arg & parameter == parameter_arg) %>%
      dplyr::distinct() %>%
      dplyr::group_by(DT_round) %>%
      dplyr::summarize(mean = as.numeric(mean(value, na.rm = T)),
                       diff = abs(min(value, na.rm = T) - max(value, na.rm = T)),
                       n_obs = n()) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(DT_round) %>%
      # pad the dataset so that all 15-min timestamps are present
      padr::pad(by = "DT_round", interval = "15 min") %>%
      # join the field notes
      dplyr::mutate(DT_join = as.character(DT_round),
                    site = site_arg,
                    parameter = parameter_arg,
                    flag = NA) %>%
      dplyr::left_join(dplyr::filter(dplyr::select(site_field_notes, sonde_employed,
                                                   last_site_visit, DT_join, site,
                                                   visit_comments, sensor_malfunction,
                                                   cals_performed)),
                       by = c('DT_join', 'site')) %>%
      # join the metaparameter dfs
      dplyr::left_join(site_metaparams_list[[site_arg]], by = "DT_join") %>%
      dplyr::mutate(DT_round = lubridate::as_datetime(DT_join, tz = "MST"))
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
#' @title Update Historical Flag List
#'
#' @description
#' A function that updates the historical flag list with new flagged data.
#'
#' @param new_flagged_data A list of data frames that have been flagged.
#'
#' @param historical_flagged_data A list of data frames that have been flagged.
#'
#' @return A list of data frames that have been updated with new flagged data.
#'
#' @examples
#' update_historical_flag_list(new_flagged_data = all_data_flagged$`archery-Actual Conductivity`, historical_flagged_data = all_data_flagged$`archery-Actual Conductivity`)
#' update_historical_flag_list(new_flagged_data = all_data_flagged$`boxelder-Temperature`, historical_flagged_data = all_data_flagged$`boxelder-Temperature`)

update_historical_flag_list <- function(new_flagged_data, historical_flagged_data){

  # Get the matching index names
  matching_indexes <- intersect(names(new_flagged_data), names(historical_flagged_data))

  # bind new_flagged_data and historical_flagged_data together
  updated_historical_flag_list <- map(matching_indexes, function(index) {

    old <- historical_flagged_data[[index]] %>%
      filter(DT_round < ymd_hms(max(DT_round) - hours(24), tz = "MST")) %>% # this is the antijoin step (but i thought we were doing 24 hours?) -jd changing this to 24 removes the duplicate problem
      mutate(last_site_visit = force_tz(last_site_visit, tzone = "MST"))

    bind_rows(old, new_flagged_data[[index]]) %>%
      arrange(DT_round) %>%
      select(-historical) # if we are removing this then I think we should remove that step from combined_data -jd
  }) %>%
    set_names(matching_indexes) %>%
    keep(~ !is.null(.))

  # need to make sure that when we are doing this we are not getting rid of dfs in the RDS list, still needs to be checked -jd
    # if we are matching up the historical list(88) with the incoming list(72) and only returning the matches then we will miss
    # be returning a list that is being written that is shorter than the list that we started with... Definitely needs to be fixed.

  return(updated_historical_flag_list)

}

# error 1: tz discrepancy at the join between incoming data and HFD
  # this messes up the future pulls because the start times df will be wrong

# time_diffs <- diff(test_data$DT_round)
#
# time_diffs # these should all be 15 and there is one that is not (where we joined incoming data to the HFD)
verify_flag_data <- function(df_list_arg,
                             site_arg,
                             parameter_arg, # do we want this to be able to be null eventually?
                             flag_arg = NULL,
                             weekly_plot_object) {

  # set layer
  layer <- "weekly"

  #quit

  if (!QUIT) {
    # verification decision lists ----
    simple_decision_list <- c("PASS ALL", "FAIL NONE",
                              "PASS VALID", "FAIL FLAGGED",
                              "PASS FLAGGED", "FAIL VALID",
                              "PASS NONE", "FAIL ALL",
                              "SKIP", "QUIT")

    inspect_decision_list <- c("INSPECT ALL", "INSPECT VALID",
                               "INSPECT FLAGGED", "INSPECT SOME")

    # initialize site_param_df ----
    site_param_df <- df_list_arg[[paste0(site_arg, "-", parameter_arg)]]

    # get data from initial weekly plot object ----
    df_weekly_data <- ggplot_build(weekly_plot_object)$plot$data %>%
      filter(site == site_arg) %>%
      select(intersect(names(.), names(site_param_df)))

    df_weekly_name <- ggplot_build(weekly_plot_object)$plot$labels$title

    # initialize list to contain updated dfs OUTPUT, VERY IMPORTANT ----
    altered_df_list <- list()

    # initialize week list ----
    weekday_list <- c(1,2,3,4,5,6,7)

    # Make sure that non-removable flags are still in place ----
    non_removable_flags <- paste0(c("sensor biofouling","sonde burial",
                                    "sensor malfunction","frozen","missing data",
                                    "sonde not employed","sonde unsubmerged",
                                    "depth calibration malfunction"),
                                  collapse = "|")

    # Which data are you working on: ----
    cat("Currently verifying: ", df_weekly_name, "data. \n")

    print(weekly_plot_object)

    # Get weekly decision(s) ----
    verification_decision <- get_verification_decision(
      paste("Would you like to (pass/fail/skip/inspect/quit) ", df_weekly_name, "?\n(layer~", layer, ")(px/fx/sx/ix/q): ") # this should be weekly data
    )

    if (verification_decision %in% simple_decision_list){
      altered_df_list <- alter_verification_column(verification_decision, df_weekly_data)
    }

    if (verification_decision %in% inspect_decision_list){
      # set layer
      layer <- "daily"

      # read in the data relevant to the current site-parameter that we are doing
      # this will be used to make the inspect decision plots. pull this data from the working directory
      # site_param_df <<- df_list[[paste0(site_arg, "-", parameter_arg)]] # if parameter arg is null then pull it all in

      # inspect all weekly ----
      if (verification_decision == "INSPECT ALL"){

        for (i in unique(df_weekly_data$weekday)) {

          daily_plot_data <- site_param_df %>%
            filter((y_w %in% df_weekly_data$y_w) & (weekday == i))

          daily_plot_object <- generate_daily_plot(plot_data_arg = daily_plot_data,
                                                   df_list_arg = df_list_arg,
                                                   site_arg = site_arg,
                                                   parameter_arg = parameter_arg)

          # weekly_plot_object
          weekly_plot_object <- generate_supplemental_weekly_plot(daily_plot_data_arg = daily_plot_data,
                                                                  df_list_arg = df_list_arg,
                                                                  site_arg = site_arg,
                                                                  parameter_arg = parameter_arg) # need to make sure that these functions can read in information from outside of them

          # print ggarrange daily and weekly plots
          print(ggarrange(daily_plot_object, weekly_plot_object, nrow = 2, ncol = 1))

          # Generate the prompt for the user to answer
          day_dt <- as.character(head(daily_plot_data$DT_round, 1)) # reformat so that this reads as a date
          day_choice_prompt <- paste("Would you like to (pass/fail/skip/inspect) ", day_dt, "?\n(layer~", layer, ")(px/fx/sx/ix): ")
          daily_verification_decision <- get_verification_decision(day_choice_prompt)

          while (layer == "daily" & daily_verification_decision == "QUIT") {
            cat("You cannot quit during this step. Please choose a different option.\n")
            daily_verification_decision <- get_verification_decision(day_choice_prompt)
          }

          if (daily_verification_decision %in% simple_decision_list){
            altered_df_list[[i]] <- alter_verification_column(daily_verification_decision, daily_plot_data) # this should get put in a list
          }

          if (daily_verification_decision %in% inspect_decision_list){
            altered_df_list[[i]] <- verify_subdaily(daily_verification_decision_arg = daily_verification_decision,
                                                    daily_plot_data_arg = daily_plot_data,
                                                    daily_plot_object_arg = daily_plot_object,
                                                    weekly_plot_object_arg = weekly_plot_object
            )
          }
        }
      }

      # inspect valid weekly ----
      if (verification_decision == "INSPECT VALID") {

        # get which days to inspect
        days_to_inspect <- unique(df_weekly_data %>%
                                    filter(is.na(flag)) %>%
                                    pull(weekday))

        new_weekday_list <- weekday_list[weekday_list >= min(df_weekly_data$weekday) & weekday_list <= max(df_weekly_data$weekday)]
        non_days_to_inspect <- new_weekday_list[!new_weekday_list %in% days_to_inspect]
        # for each non-inspect day decision:
        if (length(non_days_to_inspect) > 0) {

          # get non-inspect day decision:
          noninspect_prompt_text = paste("Would you like to (pass/fail/skip) non-inspect days?\n(layer~", layer, ")(px/fx/sx): ")
          non_inspect_verification_decision <- get_verification_decision(prompt_text = noninspect_prompt_text)

          # Handle the case when the user selects an "INSPECT" or "QUIT" option for non-inspect days
          while (layer == "daily" & non_inspect_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME", "QUIT")) {
            if(non_inspect_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME")) {
              cat("You cannot inspect non-inspect days. Please choose a different option.\n")
              non_inspect_verification_decision <- get_verification_decision(prompt_text = noninspect_prompt_text)
            }
            if(non_inspect_verification_decision == "QUIT") {
              cat("You cannot quit during this step. Please choose a different option.\n")
              non_inspect_verification_decision <- get_verification_decision(prompt_text = noninspect_prompt_text)
            }
          }

          for (i in non_days_to_inspect) {
            daily_plot_data <- site_param_df %>%
              filter((y_w %in% df_weekly_data$y_w) & (weekday == i))

            altered_df_list[[i]] <- alter_verification_column(non_inspect_verification_decision, daily_plot_data)
          }
        }

        for (i in days_to_inspect) {
          daily_plot_data <- site_param_df %>%
            filter((y_w %in% df_weekly_data$y_w) & (weekday == i))

          daily_plot_object <- generate_daily_plot(plot_data_arg = daily_plot_data,
                                                   df_list_arg = df_list_arg,
                                                   site_arg = site_arg,
                                                   parameter_arg = parameter_arg)

          # weekly_plot_object
          weekly_plot_object <- generate_supplemental_weekly_plot(daily_plot_data_arg = daily_plot_data,
                                                                  df_list_arg = df_list_arg,
                                                                  site_arg = site_arg,
                                                                  parameter_arg = parameter_arg) # need to make sure that these functions can read in information from outside of them

          # print ggarrange daily and weekly plots
          print(ggarrange(daily_plot_object, weekly_plot_object, nrow = 2, ncol = 1))

          # Generate the prompt for the user to answer
          day_dt <- as.character(head(daily_plot_data$DT_round, 1)) # reformat so that this reads as a date
          day_choice_prompt <- paste("Would you like to (pass/fail/skip/inspect) ", day_dt, "?\n(layer~", layer, ")(px/fx/sx/ix): ")
          daily_verification_decision <- get_verification_decision(day_choice_prompt)

          while (daily_verification_decision == "QUIT") {
            cat("You cannot quit during this step. Please choose a different option.\n")
            daily_verification_decision <- get_verification_decision(day_choice_prompt)
          }

          if (daily_verification_decision %in% simple_decision_list){
            altered_df_list[[i]] <- alter_verification_column(daily_verification_decision, daily_plot_data) # this should get put in a list
          }

          if (daily_verification_decision %in% inspect_decision_list){
            altered_df_list[[i]] <- verify_subdaily(daily_verification_decision_arg = daily_verification_decision,
                                                    daily_plot_data_arg = daily_plot_data,
                                                    daily_plot_object_arg = daily_plot_object,
                                                    weekly_plot_object_arg = weekly_plot_object
            )
          }
        }
      }

      # inspect flagged weekly ----
      if (verification_decision == "INSPECT FLAGGED") {

        # get which days to inspect
        days_to_inspect <- unique(df_weekly_data %>%
                                    filter(!is.na(flag)) %>%
                                    pull(weekday))

        new_weekday_list <- weekday_list[weekday_list >= min(df_weekly_data$weekday) & weekday_list <= max(df_weekly_data$weekday)]
        non_days_to_inspect <- new_weekday_list[!new_weekday_list %in% days_to_inspect]

        # for each non-inspect day decision:
        if (length(non_days_to_inspect) > 0) {

          # get non-inspect day decision:
          noninspect_prompt_text = paste("Would you like to (pass/fail/skip) non-inspect days?\n(layer~", layer, ")(px/fx/sx): ")
          non_inspect_verification_decision <- get_verification_decision(prompt_text = noninspect_prompt_text)

          # Handle the case when the user selects an "INSPECT" or "QUIT" option for non-inspect days
          while (layer == "daily" & non_inspect_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME", "QUIT")) {
            if(non_inspect_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME")) {
              cat("You cannot inspect non-inspect days. Please choose a different option.\n")
              non_inspect_verification_decision <- get_verification_decision(prompt_text = noninspect_prompt_text)
            }
            if(non_inspect_verification_decision == "QUIT") {
              cat("You cannot quit during this step. Please choose a different option.\n")
              non_inspect_verification_decision <- get_verification_decision(prompt_text = noninspect_prompt_text)
            }
          }

          for (i in non_days_to_inspect) {
            daily_plot_data <- site_param_df %>%
              filter((y_w %in% df_weekly_data$y_w) & (weekday == i))

            altered_df_list[[i]] <- alter_verification_column(non_inspect_verification_decision, daily_plot_data)
          }
        }

        for (i in days_to_inspect) {
          daily_plot_data <- site_param_df %>%
            filter((y_w %in% df_weekly_data$y_w) & (weekday == i))

          daily_plot_object <- generate_daily_plot(plot_data_arg = daily_plot_data,
                                                   df_list_arg = df_list_arg,
                                                   site_arg = site_arg,
                                                   parameter_arg = parameter_arg)


          # weekly_plot_object
          weekly_plot_object <- generate_supplemental_weekly_plot(daily_plot_data_arg = daily_plot_data,
                                                                  df_list_arg = df_list_arg,
                                                                  site_arg = site_arg,
                                                                  parameter_arg = parameter_arg) # need to make sure that these functions can read in information from outside of them

          # print ggarrange daily and weekly plots
          print(ggarrange(daily_plot_object, weekly_plot_object, nrow = 2, ncol = 1))

          # Generate the prompt for the user to answer
          day_dt <- as.character(head(daily_plot_data$DT_round, 1)) # reformat so that this reads as a date
          day_choice_prompt <- paste("Would you like to (pass/fail/skip/inspect) ", day_dt, "?\n(layer~", layer, ")(px/fx/sx/ix): ")
          daily_verification_decision <- get_verification_decision(day_choice_prompt)

          while (layer == "daily" & daily_verification_decision == "QUIT") {
            cat("You cannot quit during this step. Please choose a different option.\n")
            daily_verification_decision <- get_verification_decision(day_choice_prompt)
          }

          if (daily_verification_decision %in% simple_decision_list){
            altered_df_list[[i]] <- alter_verification_column(daily_verification_decision, daily_plot_data) # this should get put in a list
          }

          if (daily_verification_decision %in% inspect_decision_list){
            altered_df_list[[i]] <- verify_subdaily(daily_verification_decision_arg = daily_verification_decision,
                                                    daily_plot_data_arg = daily_plot_data,
                                                    daily_plot_object_arg = daily_plot_object,
                                                    weekly_plot_object_arg = weekly_plot_object
            )
          }
        }
      }

      # inspect some weekly ----
      if (verification_decision == "INSPECT SOME") {

        # get which days to inspect
        days_to_inspect <- get_weekly_inspection_decision(weekly_plot_data = df_weekly_data) # rename to daily ***

        new_weekday_list <- weekday_list[weekday_list >= min(df_weekly_data$weekday) & weekday_list <= max(df_weekly_data$weekday)]
        non_days_to_inspect <- new_weekday_list[!new_weekday_list %in% days_to_inspect]

        # for each non-inspect day decision:
        if (length(non_days_to_inspect)>0) {

          # get non-inspect day decision:
          noninspect_prompt_text = paste("Would you like to (pass/fail/skip) non-inspect days?\n(layer~", layer, ")(px/fx/sx): ")
          non_inspect_verification_decision <- get_verification_decision(prompt_text = noninspect_prompt_text)

          # Handle the case when the user selects an "INSPECT" or "QUIT" option for non-inspect days
          while (layer == "daily" & non_inspect_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME", "QUIT")) {
            if(non_inspect_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME")) {
              cat("You cannot inspect non-inspect days. Please choose a different option.\n")
              non_inspect_verification_decision <- get_verification_decision(prompt_text = noninspect_prompt_text)
            }
            if(non_inspect_verification_decision == "QUIT") {
              cat("You cannot quit during this step. Please choose a different option.\n")
              non_inspect_verification_decision <- get_verification_decision(prompt_text = noninspect_prompt_text)
            }
          }

          for (i in non_days_to_inspect) {
            daily_plot_data <- site_param_df %>%
              filter((y_w %in% df_weekly_data$y_w) & (weekday == i))

            altered_df_list[[i]] <- alter_verification_column(non_inspect_verification_decision, daily_plot_data)
          }
        }

        for (i in days_to_inspect) {

          daily_plot_data <- site_param_df %>%
            filter((y_w %in% df_weekly_data$y_w) & (weekday == i))

          daily_plot_object <- generate_daily_plot(plot_data_arg = daily_plot_data,
                                                   df_list_arg = df_list_arg,
                                                   site_arg = site_arg,
                                                   parameter_arg = parameter_arg)

          # weekly_plot_object

          weekly_plot_object <- generate_supplemental_weekly_plot(daily_plot_data_arg = daily_plot_data,
                                                                  df_list_arg = df_list_arg,
                                                                  site_arg = site_arg,
                                                                  parameter_arg = parameter_arg)

          # print ggarrange daily and weekly plots
          print(ggarrange(daily_plot_object, weekly_plot_object, nrow = 2, ncol = 1))

          # Generate the prompt for the user to answer
          day_dt <- as.character(head(daily_plot_data$DT_round, 1)) # reformat so that this reads as a date
          day_choice_prompt <- paste("Would you like to (pass/fail/skip/inspect) ", day_dt, "?\n(layer~", layer, ")(px/fx/sx/ix): ")
          daily_verification_decision <- get_verification_decision(day_choice_prompt)

          while (layer == "daily" & daily_verification_decision == "QUIT") {
            cat("You cannot quit during this step. Please choose a different option.\n")
            daily_verification_decision <- get_verification_decision(day_choice_prompt)
          }

          if (daily_verification_decision %in% simple_decision_list){
            # browser()
            altered_df_list[[i]] <- alter_verification_column(daily_verification_decision, daily_plot_data) # this should get put in a list
          }

          if (daily_verification_decision %in% inspect_decision_list){
            altered_df_list[[i]] <- verify_subdaily(daily_verification_decision_arg = daily_verification_decision,
                                                    daily_plot_data_arg = daily_plot_data,
                                                    daily_plot_object_arg = daily_plot_object,
                                                    weekly_plot_object_arg = weekly_plot_object
            )
          }
        }
      }

    }


    if (is_tibble(altered_df_list)) {
      altered_df_list <- altered_df_list %>%
        mutate(mean_verified = if_else(grepl(non_removable_flags, flag, ignore.case = TRUE), NA, mean_verified),
               is_verified = TRUE,
               verification_status = ifelse(grepl(non_removable_flags, flag, ignore.case = TRUE), "PASS", verification_status))

    } else if (is.list(altered_df_list)) {
      if (!length(altered_df_list) %in% c(51,49,46)){
        # browser()
        altered_df_list <- map(altered_df_list, function(df) {
          if(!is.null(df)) {
            # browser()
            df %>%
              mutate(mean_verified = if_else(grepl(non_removable_flags, flag, ignore.case = TRUE), NA, mean_verified),
                     is_verified = TRUE,
                     verification_status = ifelse(grepl(non_removable_flags, flag, ignore.case = TRUE), "PASS", verification_status))
          }
        })
      } else {
        # browser()
        altered_df_list <- altered_df_list %>%
          mutate(mean_verified = if_else(grepl(non_removable_flags, flag, ignore.case = TRUE), NA, mean_verified),
                 is_verified = TRUE,
                 verification_status = ifelse(grepl(non_removable_flags, flag, ignore.case = TRUE), "PASS", verification_status))
      }

    }

    return(altered_df_list) # this should not be an altered df, but rather a list of dfs that have been altered.
    # this means that the output of this method will be a list that contains a list of lists

    # If decision in quit (this can only happen before they start verifying a week of data)
    # save the data into the intermediary directory
  } else if (QUIT) {
    df_weekly_name <- ggplot_build(weekly_plot_object)$plot$labels$title
    cat("quitting out of", df_weekly_name, "\n")
    return(NULL)
  }

}

verify_subdaily <- function(daily_verification_decision_arg,
                            daily_plot_data_arg,
                            daily_plot_object_arg,
                            weekly_plot_object_arg){

  # set layer
  layer <- "subdaily"

  # inspect all
  if (daily_verification_decision_arg == "INSPECT ALL") {

    for (j in 1:nrow(daily_plot_data_arg)) {

      # create a data frame with a single row containing the desired DT_round value
      # and add it to the week plot and day plot
      flag_demarcation_data <- tibble(DT_round = daily_plot_data_arg$DT_round[j])

      # ignore flags that you cannot alter
      daily_plot_object_new <- daily_plot_object_arg +
        geom_vline(data = flag_demarcation_data,
                   aes(xintercept = DT_round,
                       color = "data point of interest",
                       # size = 2, # fix this later
                       alpha = 0.05))

      weekly_plot_object_new <- weekly_plot_object_arg +
        geom_vline(data = flag_demarcation_data,
                   aes(xintercept = DT_round,
                       color = "data point of interest",
                       # size = 2, # fix this later
                       alpha = 0.05))

      # print ggarrange daily and weekly plots
      print(ggarrange(daily_plot_object_new, weekly_plot_object_new, nrow = 2, ncol = 1))

      # generate the prompt for the user to answer
      # the user should not be able to inspect again
      day_dt <- daily_plot_data_arg$DT_join[j] # reformat so that this reads as a date
      dt_choice_prompt <- paste("Would you like to (pass/fail/skip) ", day_dt, "?\n(layer~", layer, ")(px/fx/sx)")
      dt_verification_decision <- get_verification_decision(dt_choice_prompt)

      # Handle the case when the user selects an "INSPECT" or "QUIT" option inappropriately
      while (layer == "subdaily" & dt_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME", "QUIT")) {
        if (dt_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME")){
          cat("You cannot inspect at this level. Please choose a different option.\n")
          dt_verification_decision <- get_verification_decision(dt_choice_prompt)
        }
        if (dt_verification_decision == "QUIT"){
          cat("You cannot quit during this step. Please choose a different option.\n")
          dt_verification_decision <- get_verification_decision(dt_choice_prompt)
        }
      }

      # update the specific row in df_daily_data based on the index
      j_index <- which(daily_plot_data_arg$DT_join == day_dt)

      altered_row <- alter_verification_column(dt_verification_decision, daily_plot_data_arg[j_index, ])

      daily_plot_data_arg <- daily_plot_data_arg %>%
        mutate(
          mean_verified = if_else(row_number() == j_index, altered_row$mean_verified, mean_verified),
          is_verified = if_else(row_number() == j_index, altered_row$is_verified, is_verified),
          verification_status = if_else(row_number() == j_index, altered_row$verification_status, verification_status)
        )
    }
    # inform the user that they are done verifying which ever df
    cat("Finished verifying: ", as.character(min(daily_plot_data_arg$DT_round)), " data.\n")

  }

  # inspect valid
  if (daily_verification_decision_arg == "INSPECT VALID") {

    # get the dt of the valid data to inspect to get indexes of that data
    valid_daily_plot_data_dt_list <- daily_plot_data_arg %>%
      filter(is.na(flag)) %>%
      pull(DT_join)

    # get the dt of the flagged data to get the non-inspect decision
    new_dt_list <- pull(daily_plot_data_arg, DT_join)
    non_dts_to_inspect <- new_dt_list[!new_dt_list %in% valid_daily_plot_data_dt_list]

    # decide what to do with the flagged data ----
    # for each non-inspect dt decision:
    if (length(non_dts_to_inspect) > 0) {

      # get non-inspect day decision:
      noninspect_prompt_text = paste("Would you like to (pass/fail/skip) non-inspect dts?\n(layer~", layer, ")(px/fx/sx): ")
      non_inspect_verification_decision <- get_verification_decision(noninspect_prompt_text)

      # Handle the case when the user selects an "INSPECT" or "QUIT" option inappropriately
      while (non_inspect_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME", "QUIT")) {
        if (non_inspect_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME")){
          cat("You cannot inspect non-inspect dts. Please choose a different option.\n")
          non_inspect_verification_decision <- get_verification_decision(noninspect_prompt_text)
        }
        if (non_inspect_verification_decision == "QUIT"){
          cat("You cannot quit during this step. Please choose a different option.\n")
          non_inspect_verification_decision <- get_verification_decision(noninspect_prompt_text)
        }
      }

      for (j in non_dts_to_inspect) {
        # j <- format(as.POSIXct(j, tz = "MST"), "%Y-%m-%d %H:%M:%S")

        j_index <- which((daily_plot_data_arg$DT_join) == j)

        altered_row <- alter_verification_column(non_inspect_verification_decision, daily_plot_data_arg[j_index, ])

        daily_plot_data_arg <- daily_plot_data_arg %>%
          mutate(
            mean_verified = if_else(row_number() == j_index, altered_row$mean_verified, mean_verified),
            is_verified = if_else(row_number() == j_index, altered_row$is_verified, is_verified),
            verification_status = if_else(row_number() == j_index, altered_row$verification_status, verification_status)
          )

        # some sort of print statement to show what happened to the non-target dts
      }
    }

    # update those flagged data in the df

    for (j in valid_daily_plot_data_dt_list) {

      # j <- format(as.POSIXct(j, tz = "MST"), "%Y-%m-%d %H:%M:%S")

      j_index <- which((daily_plot_data_arg$DT_join) == j)

      # create a data frame with a single row containing the desired DT_round value
      # and add it to the week plot and day plot
      flag_demarcation_data <- tibble(DT_round = daily_plot_data_arg$DT_round[j_index])

      daily_plot_object_new <- daily_plot_object_arg +
        geom_vline(data = flag_demarcation_data,
                   aes(xintercept = DT_round,
                       color = "data point of interest",
                       # size = 2, # fix this later
                       alpha = 0.05))

      weekly_plot_object_new <- weekly_plot_object_arg +
        geom_vline(data = flag_demarcation_data,
                   aes(xintercept = DT_round,
                       color = "data point of interest",
                       # size = 2, # fix this later
                       alpha = 0.05))

      # print ggarrange daily and weekly plots
      print(ggarrange(daily_plot_object_new, weekly_plot_object_new, nrow = 2, ncol = 1))

      # generate the prompt for the user to answer
      # the user should not be able to inspect again
      day_dt <- as.character(daily_plot_data_arg$DT_round[j_index]) # reformat so that this reads as a date
      dt_choice_prompt <- paste("Would you like to (pass/fail/skip) ", day_dt, "?\n(layer~", layer, ")(px/fx/sx)")
      dt_verification_decision <- get_verification_decision(dt_choice_prompt)

      # Handle the case when the user selects an "INSPECT" or "QUIT" option inappropriately
      while (layer == "subdaily" & dt_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME", "QUIT")) {
        if (dt_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME")){
          cat("You cannot inspect at this level. Please choose a different option.\n")
          dt_verification_decision <- get_verification_decision(dt_choice_prompt)
        }
        if (dt_verification_decision == "QUIT"){
          cat("You cannot quit during this step. Please choose a different option.\n")
          dt_verification_decision <- get_verification_decision(dt_choice_prompt)
        }
      }

      # update the specific row in df_daily_data based on the index
      # row_index <- which(as.character(daily_plot_data_arg$DT_round) == day_dt)
      altered_row <- alter_verification_column(dt_verification_decision, daily_plot_data_arg[j_index, ])

      daily_plot_data_arg <- daily_plot_data_arg %>%
        mutate(
          mean_verified = if_else(row_number() == j_index, altered_row$mean_verified, mean_verified),
          is_verified = if_else(row_number() == j_index, altered_row$is_verified, is_verified),
          verification_status = if_else(row_number() == j_index, altered_row$verification_status, verification_status)
        )
    }
    # inform the user that they are done verifying which ever df
    cat("Finished verifying: ", as.character(min(daily_plot_data_arg$DT_round)), " data.\n")
    # return the altered df
    # altered_df_list_arg[[i]] <- daily_plot_data_arg # make sure that this is the updated data frame

  }

  # inspect flagged
  if (daily_verification_decision_arg == "INSPECT FLAGGED") {

    # get the dt of the flagged data to inspect to get indexes of that data
    flagged_daily_plot_data_dt_list <- daily_plot_data_arg %>%
      filter(!is.na(flag)) %>%
      pull(DT_join)

    # get the dt of the flagged data to get the non-inspect decision
    new_dt_list <- pull(daily_plot_data_arg, DT_join)
    non_dts_to_inspect <- new_dt_list[!new_dt_list %in% flagged_daily_plot_data_dt_list]

    # decide what to do with the flagged data ----
    # for each non-inspect dt decision:
    if (length(non_dts_to_inspect) > 0) {

      # get non-inspect day decision:
      noninspect_prompt_text = paste("Would you like to (pass/fail/skip) non-inspect dts?\n(layer~", layer, ")(px/fx/sx): ")
      non_inspect_verification_decision <- get_verification_decision(noninspect_prompt_text)

      # Handle the case when the user selects an "INSPECT" or "QUIT" option inappropriately
      while (non_inspect_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME", "QUIT")) {
        if (non_inspect_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME")){
          cat("You cannot inspect non-inspect dts. Please choose a different option.\n")
          non_inspect_verification_decision <- get_verification_decision(noninspect_prompt_text)
        }
        if (non_inspect_verification_decision == "QUIT"){
          cat("You cannot quit during this step. Please choose a different option.\n")
          non_inspect_verification_decision <- get_verification_decision(noninspect_prompt_text)
        }
      }

      for (j in non_dts_to_inspect) {
        # j <- format(as.POSIXct(j, tz = "MST"), "%Y-%m-%d %H:%M:%S") # time zone

        j_index <- which((daily_plot_data_arg$DT_join) == j)

        altered_row <- alter_verification_column(non_inspect_verification_decision, daily_plot_data_arg[j_index, ])

        daily_plot_data_arg <- daily_plot_data_arg %>%
          mutate(
            mean_verified = if_else(row_number() == j_index, altered_row$mean_verified, mean_verified),
            is_verified = if_else(row_number() == j_index, altered_row$is_verified, is_verified),
            verification_status = if_else(row_number() == j_index, altered_row$verification_status, verification_status)
            )
        # some sort of print statement to show what happened to the non-target dts
      }
    }

    # update those flagged data in the df

    for (j in flagged_daily_plot_data_dt_list) {

      # j <- format(as.POSIXct(j, tz = "MST"), "%Y-%m-%d %H:%M:%S")

      j_index <- which((daily_plot_data_arg$DT_join) == j)

      # create a data frame with a single row containing the desired DT_round value
      # and add it to the week plot and day plot
      flag_demarcation_data <- tibble(DT_round = daily_plot_data_arg$DT_round[j_index])

      daily_plot_object_new <- daily_plot_object_arg +
        geom_vline(data = flag_demarcation_data,
                   aes(xintercept = DT_round,
                       color = "data point of interest",
                       # size = 2, # fix this later
                       alpha = 0.05))

      weekly_plot_object_new <- weekly_plot_object_arg +
        geom_vline(data = flag_demarcation_data,
                   aes(xintercept = DT_round,
                       color = "data point of interest",
                       # size = 2, # fix this later
                       alpha = 0.05))


      # print ggarrange daily and weekly plots
      print(ggarrange(daily_plot_object_new, weekly_plot_object_new, nrow = 2, ncol = 1))

      # generate the prompt for the user to answer
      # the user should not be able to inspect again
      day_dt <- as.character(daily_plot_data_arg$DT_round[j_index]) # reformat so that this reads as a date
      dt_choice_prompt <- paste("Would you like to (pass/fail/skip) ", day_dt, "?\n(layer~", layer, ")(px/fx/sx)")
      dt_verification_decision <- get_verification_decision(dt_choice_prompt)

      # Handle the case when the user selects an "INSPECT" or "QUIT" option inappropriately
      while (layer == "subdaily" & dt_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME", "QUIT")) {
        if (dt_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME")){
          cat("You cannot inspect at this level. Please choose a different option.\n")
          dt_verification_decision <- get_verification_decision(dt_choice_prompt)
        }
        if (dt_verification_decision == "QUIT"){
          cat("You cannot quit during this step. Please choose a different option.\n")
          dt_verification_decision <- get_verification_decision(dt_choice_prompt)
        }
      }

      # update the specific row in df_daily_data based on the index
      altered_row <- alter_verification_column(dt_verification_decision, daily_plot_data_arg[j_index, ])

      daily_plot_data_arg <- daily_plot_data_arg %>%
        mutate(
          mean_verified = if_else(row_number() == j_index, altered_row$mean_verified, mean_verified),
          is_verified = if_else(row_number() == j_index, altered_row$is_verified, is_verified),
          verification_status = if_else(row_number() == j_index, altered_row$verification_status, verification_status)
        )

    }
    # inform the user that they are done verifying which ever df
    cat("Finished verifying: ", as.character(min(daily_plot_data_arg$DT_round)), " data.\n")

  }

  # inspect some
  if (daily_verification_decision_arg == "INSPECT SOME") {

    # get which dts to inspect
    dt_intervals_to_inspect <- get_dt_inspection_decisions(daily_plot_data = daily_plot_data_arg)

    # decided what to do with the non-inspect data
    new_dt_list <- pull(daily_plot_data_arg, DT_round)
    non_dts_to_inspect <- new_dt_list[!new_dt_list %within% dt_intervals_to_inspect]

    # decide what to do with the flagged data ----
    # for each non-inspect dt decision:
    if (length(non_dts_to_inspect) > 0) {

      # get non-inspect day decision:
      noninspect_prompt_text = paste("Would you like to (pass/fail/skip) non-inspect dts?\n(layer~", layer, ")(px/fx/sx): ")
      non_inspect_verification_decision <- get_verification_decision(noninspect_prompt_text)

      # Handle the case when the user selects an "INSPECT" or "QUIT" option inappropriately
      while (non_inspect_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME", "QUIT")) {
        if (non_inspect_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME")){
          cat("You cannot inspect non-inspect dts. Please choose a different option.\n")
          non_inspect_verification_decision <- get_verification_decision(noninspect_prompt_text)
        }
        if (non_inspect_verification_decision == "QUIT"){
          cat("You cannot quit during this step. Please choose a different option.\n")
          non_inspect_verification_decision <- get_verification_decision(noninspect_prompt_text)
        }
      }

      for (j in non_dts_to_inspect) {
        j_time <- format(as.POSIXct(j, tz = "MST"), "%H:%M:%S")
        j <- if_else(
          j_time == "00:00:00",
          format(as.POSIXct(j, tz = "MST"), "%Y-%m-%d"),
          format(as.POSIXct(j, tz = "MST"), "%Y-%m-%d %H:%M:%S")
        )

        j_index <- which((daily_plot_data_arg$DT_join) == j)

        # daily_plot_data_arg[j_index, ] <- alter_verification_column(non_inspect_verification_decision, daily_plot_data_arg[j_index, ])

        altered_row <- alter_verification_column(non_inspect_verification_decision, daily_plot_data_arg[j_index, ])

        daily_plot_data_arg <- daily_plot_data_arg %>%
          mutate(
            mean_verified = if_else(row_number() == j_index, altered_row$mean_verified, mean_verified),
            is_verified = if_else(row_number() == j_index, altered_row$is_verified, is_verified),
            verification_status = if_else(row_number() == j_index, altered_row$verification_status, verification_status)
          )
        # some sort of print statement to show what happened to the non-target dts
      }
    }

    # update those flagged data in the df

    for (j in dt_intervals_to_inspect) { # make sure that j can not be 0

      interval_df <- daily_plot_data_arg %>%
        filter(DT_round %within% j)

      for (k in 1:nrow(interval_df)) {

        k_index <- which((daily_plot_data_arg$DT_round) == interval_df$DT_round[k])

        # create a data frame with a single row containing the desired DT_round value
        # and add it to the week plot and day plot
        flag_demarcation_data <- tibble(DT_round = daily_plot_data_arg$DT_round[k_index])

        daily_plot_object_new <- daily_plot_object_arg +
          geom_vline(data = flag_demarcation_data,
                     aes(xintercept = DT_round,
                         color = "data point of interest",
                         # size = 2, # fix this later
                         alpha = 0.05))

        weekly_plot_object_new <- weekly_plot_object_arg +
          geom_vline(data = flag_demarcation_data,
                     aes(xintercept = DT_round,
                         color = "data point of interest",
                         # size = 2, # fix this later
                         alpha = 0.05))


        # print ggarrange daily and weekly plots
        print(ggarrange(daily_plot_object_new, weekly_plot_object_new, nrow = 2, ncol = 1))

        # generate the prompt for the user to answer
        # the user should not be able to inspect again
        day_dt <- as.character(daily_plot_data_arg$DT_round[k_index]) # reformat so that this reads as a date
        dt_choice_prompt <- paste("Would you like to (pass/fail/skip) ", day_dt, "?\n(layer~", layer, ")(px/fx/sx)")
        dt_verification_decision <- get_verification_decision(dt_choice_prompt)

        # Handle the case when the user selects an "INSPECT" or "QUIT" option inappropriately
        while (layer == "subdaily" & dt_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME", "QUIT")) {
          if (dt_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME")){
            cat("You cannot inspect at this level. Please choose a different option.\n")
            dt_verification_decision <- get_verification_decision(dt_choice_prompt)
          }
          if (dt_verification_decision == "QUIT"){
            cat("You cannot quit during this step. Please choose a different option.\n")
            dt_verification_decision <- get_verification_decision(dt_choice_prompt)
          }
        }

        # update the specific row in df_daily_data based on the index
        # row_index <- which(as.character(daily_plot_data_arg$DT_round) == day_dt)
        # daily_plot_data_arg[k_index, ] <- alter_verification_column(dt_verification_decision, daily_plot_data_arg[k_index, ])

        altered_row <- alter_verification_column(dt_verification_decision, daily_plot_data_arg[k_index, ])

        daily_plot_data_arg <- daily_plot_data_arg %>%
          mutate(
            mean_verified = if_else(row_number() == k_index, altered_row$mean_verified, mean_verified),
            is_verified = if_else(row_number() == k_index, altered_row$is_verified, is_verified),
            verification_status = if_else(row_number() == k_index, altered_row$verification_status, verification_status)
          )
      }
    }
    # inform the user that they are done verifying which ever df
    cat("Finished verifying: ", as.character(min(daily_plot_data_arg$DT_round)), " data.\n")
  }

  return(daily_plot_data_arg)

}

add_depth_shift_flag <- function(df){


  cal_flag <- function(df, condition_arg, description_arg) {
    df <- df %>% mutate(depth_change = case_when(
      {{condition_arg}} ~ if_else(is.na(depth_change), paste(description_arg),
                                  ifelse(!grepl(description_arg, depth_change), paste(depth_change, description_arg, sep = ";\n"), depth_change)),
      TRUE ~ depth_change))
    return(df)
  }

  # Function to add a column if it doesn't already exist
  add_column_if_not_exists <- function(df, column_name, default_value) {
    if (!column_name %in% colnames(df)) {
      df <- df %>% mutate(!!sym(column_name) := default_value)
    }
    return(df)
  }

    shift_dates <- read_csv('data/qaqc/level_shifts.csv') %>%
      filter(type == "sonde moved") %>%
      add_column_if_not_exists(., "depth_change", NA) %>%
      mutate(DT_join = as.character(ymd_hms(DT_join))) %>%
      left_join(df, ., by = c("site", "DT_join")) %>%
      cal_flag(type == "sonde moved", "sonde moved") %>%
      select(-type)

    return(shift_dates)

}


add_do_drops <- function(df){

if("DO" %in% df$parameter){

  df <- df %>%
    add_flag(back1 - mean >= 0.5 & front1 - mean >= 0.5, "DO interference")

  return(df)

} else {

  return(df)

}

  }

add_drift_flag <- function(df){

  # Only test for biofilm growth on turbidity sensors
  if("Turbidity" %in% df$parameter){

    # subset data to turbidity and conductivity only
    sub <- df %>%
      dplyr::filter(parameter %in% c("Turbidity", "Specific Conductivity")) %>%
      dplyr::select(DT_round, DT_join, parameter, mean) %>%
      tidyr::pivot_wider(names_from = parameter, values_from = mean)
    names(sub) <- make.names(names(sub))

    # Check if a given window's worth of sensor data has a linear
    # relationship with time (indicates a "steady" drift)
    progressive_drift <- function(x) {
      # Only assess time series with less than 90% missing data in
      # the rolling window:
      if(length(x[!is.na(x)]) > (length(x) - (length(x)*0.1))){

        # Fit linear model
        model <- lm(x ~ c(1:length(x)), na.action = na.omit)

        # Extract R-squared value
        r_squared <- summary(model)$r.squared

        # Return R-squared value
        return(r_squared)

      } else {

        # If not enough data, set slope to 0
        no_slope <- 0

        return(no_slope)

      }
    }

    # Function to check if a selected time window's mean R-squared value is
    # at least 60% (i.e., fairly high linear relationship with time indicates a "steady" drift)
    check_too_steady <- function(x) {
      mean(x) >= 0.60
    }

    # Function that uses all the functions above to see if a given time window's R-squared with time is strong. If the one-day OR the three-day slope
    # for a selected parameter is greater than 60%, we determine it has "failed" (i.e., drift seems to exist)

    biofilm_tester <- function(data = sub, col){

      data %>%
        data.table::data.table() %>%
        #bind_rows(flagged_data_dfs[["prospect-Turbidity"]]) %>%
        # dplyr::filter(lubridate::as_date(DT_round) >= "2022-09-10" & lubridate::as_date(DT_round) <= "2022-10-15") %>%
        dplyr::mutate(r2_s_right = data.table::frollapply(!!sym(col), n = 96, FUN = progressive_drift, align = "right", fill = NA),
                      r2_s_center = data.table::frollapply(!!sym(col), n = 96, FUN = progressive_drift, align = "left", fill = NA),
                      r2_l_right = data.table::frollapply(!!sym(col), n = 288, FUN = progressive_drift, align = "right", fill = NA),
                      r2_l_center = data.table::frollapply(!!sym(col), n = 288, FUN = progressive_drift, align = "left", fill = NA),
                      tightest_r = pmax(r2_s_center, r2_s_right, r2_l_center, r2_l_right, na.rm = TRUE),
                      failed = data.table::frollapply(tightest_r, n = 96, FUN = check_too_steady, align = "right", fill = NA)) %>%
        dplyr::select(!!(col) := "failed",
                      DT_join)
    }

    biofouling <- c("Turbidity", "Specific.Conductivity") %>%
      purrr::map(~biofilm_tester(data = sub, col = .)) %>%
      dplyr::bind_cols() %>%
      dplyr::rename(DT_join = DT_join...2) %>%
      dplyr::mutate(parameter = "Turbidity") %>%
      dplyr::right_join(., df, by = c("DT_join", "parameter")) %>%
      # If a steady slope is existing through time for turbidity, but NOT for conductivity, it is likely
      # sensor bio-fouling
      add_flag(., Turbidity == 1 & Specific.Conductivity != 1 & !grepl("drift", flag), "drift") %>%
      dplyr::select(-c(Turbidity, DT_join...4, Specific.Conductivity))

    return(biofouling)

  } else {

    return(df)}

}
#
# a=ggplot(drift_tested)+
#   geom_point(aes(x = DT_round, y = r_squaredSC_short ), color = "blue") +
#   geom_point(aes(x = DT_round, y = r_squaredT_short ), color = "red") +
#   geom_point(aes(x = DT_round, y = r_squaredDO_short), color = "black")
#
# b=ggplot(drift_tested)+
#   geom_point(aes(x = DT_round, y = r_squaredSC_long ), color = "blue") +
#   geom_point(aes(x = DT_round, y = r_squaredT_long ), color = "red") +
#   geom_point(aes(x = DT_round, y = r_squaredDO_long), color = "black")
#
# c=ggplot(drift_tested %>% filter(as_date(DT_round) >= "2022-10-10")) +
#   #geom_point(aes(x = DT_round, y = Specific.Conductivity)) +
#   geom_point(aes(x = DT_round, y = Turbidity, color = as.character(failed))) #+
# #geom_point(aes(x = DT_round, y = DO))
#
# ggpubr::ggarrange(a, b, c, ncol = 1)
#


#' @title Add field related flags to a data frame based on field notes.
#'
#' @description
#' A function that checks 3 different and separate conditions related to the
#' field notes and adds the corresponding flags accordingly. The "sonde not
#' employed" flag is added if the sonde was not employed at the site.
#' "site visit" flag is added if the last site visit date is the same as the
#' date of the data. "sv window" flag is added if site visit flag is detected
#' within a 45 minute window. Note that this function will only work on data
#' that has already been joined to field notes.
#'
#' @param df A data frame with a `flag` column.
#'
#' @return A data frame with a `flag` column that has been updated with the
#' relevant field note flags.
#'
#' @examples
#' add_range_flags(df = all_data_flagged$`archery-Actual Conductivity`)
#' add_range_flags(df = all_data_flagged$`boxelder-Temperature`)

add_field_flag <- function(df) {

  df <- df %>%
    # flag when sonde was not employed in the river
    add_flag(sonde_employed == 1, "sonde not employed") %>% # removing the instances where we check for the flag description in the flag column when we are flagging. This is already done in `add_flag()`
    # flag when sonde was handled in a site visit
    add_flag(as.character(last_site_visit) == as.character(DT_round), "site visit")

  # Add flags for the next 60 minutes after a site visit
  for (i in 1:4) {
    df <- df %>%
      add_flag(lag(str_detect(flag, "site visit"), n = i), "sv window")
  }

  # ... and the first 15 minutes before:
  df <- df %>%
    add_flag(lead(str_detect(flag, "site visit"), n = 1), "sv window")

  return(df)

}

#' @title Underlying function for flagging data.
#'
#' @description
#' This function adds a flag to the `flag` column of a given data frame based on
#' specified conditions for each row. The name of the flag is provided
#' by the user.
#'
#' @param df A data frame with a `flag` column.
#'
#' @param condition_arg A logical statement that is evaluated in the context of
#' the data frame.
#'
#' @param description_arg Flag added to the `flag` column.
#'
#' @returns
#' An identical data frame with a `flag` column that has been updated with the
#' flag description provided.
#'
#' @examples
# add_flag(df = all_data_flagged$`archery-Actual Conductivity`, condition_arg = mean >= 100, description_arg = "exceeds 100")
# add_flag(df = all_data_flagged$`archery-Actual Conductivity`, condition_arg = mean <= 10, description_arg = "below 10")

add_flag <- function(df, condition_arg, description_arg) {
  df <- df %>% mutate(flag = case_when(
    {{condition_arg}} ~ if_else(is.na(flag), paste(description_arg),
                                ifelse(!grepl(description_arg, flag), paste(flag, description_arg, sep = ";\n"), flag)),
    TRUE ~ flag))
  return(df)
}

#' @title Add a flag if the water temperature is freezing.
#'
#' @description
#' A function designed to append the 'frozen' flag to a row if the value
#' in the `mean` column is less than or equal to 0.
#'
#' @param df A data frame with a `flag` column.
#'
#' @return A data frame with a `flag` column that has been updated with the
#' 'frozen' flag.
#'
#' @examples
#' add_frozen_flag(df = all_data_flagged$`archery-Actual Conductivity`)
#' add_frozen_flag(df = all_data_flagged$`boxelder-Temperature`)
#'
#' @seealso [flag_all_data()]

add_frozen_flag <- function(df){

  # create a df of temperature for each site
  temperature <- df %>%
    data.table::data.table() %>%
    dplyr::select(DT_round, DT_join, parameter, mean) %>%
    dplyr::filter(parameter == "Temperature") %>%
    dplyr::select(DT_join, Temperature = mean)

    # add "temperature" column to df:
    temperature_checked <- df %>%
      dplyr::left_join(., temperature, by = "DT_join") %>%
      # If water temperature is freezing, flag all parameters
      add_flag(., Temperature <= 0, "frozen") %>%
      # remove the temp column so df is identical in structure to OG df
      dplyr::select(-Temperature)

    return(temperature_checked)

}

#' @title Add malfunction flags to a data frame based on field notes and HydroVu.
#'
#' @description
#' This function adds the 'sensor malfunction' flag to a data frame by the dates,
#' sites, and parameters that members of the lab know contain erroneous data due
#' to sensor malfunctions. Note that this flag is used in two instances: first
#' when removing erroneous data from our statistics calculations and again during
#' the actual flagging step in the QAQC pipeline.
#'
#' @param df A data frame with a `flag` column.
#'
#' @param malfunction_records The malfunction records pulled from mWater
#'
#' @return A data frame with a `flag` column that has been updated with the
#' 'sensor malfunction' flag.
#'
#' @examples
#' add_range_flags(df = all_data_flagged$`archery-Actual Conductivity`)
#' add_range_flags(df = all_data_flagged$`boxelder-Temperature`)
#'
#' @seealso [grab_mWater_malfunction_notes()]

add_malfunction_flag <- function(df, malfunction_records){

  # Filter malfunction records for relevant information
  df_site <- unique(df$site)
  df_parameter <- unique(df$parameter)

  malfunction_records_filtered <- malfunction_records %>%
    filter(site == df_site) %>%
    # if parameter is NA, that means entire sonde was malfunctioning. Otherwise,
    # just the parameter listed:
    filter(is.na(parameter) | parameter == df_parameter) %>%
    mutate(end_DT = ifelse(is.na(end_DT), ymd_hms("9999-12-31 23:59:59", tz = "MST"), end_DT)) %>%
    mutate(end_DT = as.POSIXct(end_DT, tz = "MST")) %>%
    rowid_to_column()

  if (nrow(malfunction_records_filtered > 0)) {

    drift <- malfunction_records_filtered %>%
      filter(grepl("grime|gunk|drift|biofoul|biofilm", notes, ignore.case = TRUE))

    burial <- malfunction_records_filtered %>%
      filter(grepl("buried|burial|bury|sediment|roots", notes, ignore.case = TRUE))

    depth_funk <- malfunction_records_filtered %>%
      filter(grepl("improper level calibration", notes, ignore.case = TRUE))

    unsubmerged <- malfunction_records_filtered %>%
      filter(grepl("not submerged", notes, ignore.case = TRUE))

    general_malfunction <- malfunction_records_filtered %>%
      filter(!rowid %in% drift$rowid & !rowid %in% burial$rowid &
               !rowid %in% depth_funk$rowid & !rowid %in% unsubmerged$rowid)

    # make a list of date intervals where sensor was malfunctioning
    drift_interval_list <- map2(
      .x = drift$start_DT,
      .y = drift$end_DT,
      .f = ~interval(.x, .y, tz = "MST"))

    burial_interval_list <- map2(
      .x = burial$start_DT,
      .y = burial$end_DT,
      .f = ~interval(.x, .y, tz = "MST"))

    depth_interval_list <- map2(
      .x = depth_funk$start_DT,
      .y = depth_funk$end_DT,
      .f = ~interval(.x, .y, tz = "MST"))

    unsubmerged_interval_list <- map2(
      .x = unsubmerged$start_DT,
      .y = unsubmerged$end_DT,
      .f = ~interval(.x, .y, tz = "MST"))

    general_interval_list <- map2(
      .x = general_malfunction$start_DT,
      .y = general_malfunction$end_DT,
      .f = ~interval(.x, .y, tz = "MST"))

    # If DT_round is within any of the intervals where a sensor was malfunctioning, flag it
    try(df <- df %>%
          add_flag(DT_round %within% burial_interval_list, "sonde burial"))

    try(df <- df %>%
          add_flag(DT_round %within% drift_interval_list, "sensor biofouling"))

    try(df <- df %>%
          add_flag(DT_round %within% depth_interval_list, "depth calibration malfunction"))

    try(df <- df %>%
          add_flag(DT_round %within% unsubmerged_interval_list, "sonde unsubmerged"))

    try(df <- df %>%
          add_flag(DT_round %within% general_interval_list, "sensor malfunction"))

  }


  return(df)

}

#' @title Add NA flags to a data frame based on `mean` column.
#'
#' @description
#' A function designed to append the 'missing data' flag to a row if the `mean`
#' column in that row contains NA values.
#'
#' @param df A data frame with a `flag` column.
#'
#' @return A data frame with a `flag` column that has been updated with the
#' 'missing data' flag.
#'
#' @examples
#' add_range_flags(df = all_data_flagged$`archery-Actual Conductivity`)
#' add_range_flags(df = all_data_flagged$`boxelder-Temperature`)

add_na_flag <- function(df){
  df <- df %>%
    add_flag(is.na(mean), "missing data")
}

#' @title Add a flag if the value of a parameter exceeds its realistic range.
#'
#' @description
#' A function designed to append the 'outside of sensor realistic range' flag to
#' a row if the `mean` column in that row exceeds the realistic ranges that are
#' set in the `src/qaqc/sensor_real_thresholds.yml` file. These are instances
#' where a value exceeds the expected ranges for the Poudre. Note that this is
#' currently only aplied to pH, SC, and Temperature.
#'
#' @param df A data frame with a `flag` column.
#'
#' @return A data frame with a `flag` column that has been updated with the
#' 'outside of sensor realistic range' flag.
#'
#' add_range_flags(df = all_data_flagged$`archery-Actual Conductivity`)
#' add_range_flags(df = all_data_flagged$`boxelder-Temperature`)

add_realistic_flag <- function(df){

  sensor_realistic_ranges <- read_yaml("data/qaqc/sensor_real_thresholds.yml")

  # make this a non yaml solution and add it to the threshold table
  # get the parameter from the parameter column in the df of interest
  parameter_name <- unique(na.omit(df$parameter))
  # Pull the sensor specification range from the yaml file
  sensor_min <- eval(parse(text = sensor_realistic_ranges[[parameter_name]]$min))
  sensor_max <- eval(parse(text = sensor_realistic_ranges[[parameter_name]]$max))

  df <- df %>%
    # adding sensor range flags
    add_flag(parameter == parameter_name & (mean < sensor_min | mean > sensor_max),
             "outside of sensor realistic range")

    return(df)

}

#' @title Add a flag if the value in the `mean` column repeats sequentially.
#'
#' @description
#' A function designed to append the 'repeated value' flag to a row if the value
#' in the `mean` column is equal to the previous or next value in the `mean`
#' column.
#'
#' @param df A data frame with a `flag` column.
#'
#' @return A data frame with a `flag` column that has been updated with the
#' 'repeated value' flag.
#'
#' add_range_flags(df = all_data_flagged$`archery-Actual Conductivity`)
#' add_range_flags(df = all_data_flagged$`boxelder-Temperature`)

add_repeat_flag <- function(df){
  df <- df %>%
    add_flag((mean == front1 | mean == back1), "repeated value")
}

#' @title Add flags related to calculated parameter seasonal ranges.
#'
#' @description
#' A function that checks 2 different and separate conditions related to the
#' calculated parameter seasonal ranges. The 'outside of seasonal range' flag is
#' added if the `mean` value is outside the seasonal 1st - 99th percentile. The
#' 'slope violation' flag is added if the `slope_ahead` or `slope_behind` value
#' is greater than or equal to the `t_slope_behind_99` value, or if the
#' `slope_ahead` or `slope_behind` value is less than or equal to the
#' `t_slope_behind_01` value. The `t_slope_behind_99` and `t_slope_behind_01`
#' values are specific to each site-parameter. Note that this is the version of
#' the function used in `flag_all_data()`.
#'
#' @param df A data frame with a `flag` column.
#'
#' @return A data frame with a `flag` column that has been updated with the
#' relevant calculated seasonal range flags.
#'
#' @examples
#' add_range_flags(df = all_data_flagged$`archery-Actual Conductivity`)
#' add_range_flags(df = all_data_flagged$`boxelder-Temperature`)
#'
#' @seealso [flag_all_data()]

add_seasonal_flag <- function(df) {

  # get the site name from the site column in the df of interest
  site_name <- unique(na.omit(df$site))
  # get the parameter from the parameter column in the df of interest
  parameter_name <- unique(na.omit(df$parameter))

  lookup <- threshold_lookup %>%
    filter(site == site_name & parameter == parameter_name) %>%
    select(!c(site, parameter))

  df <- df %>%
    # Using seasonal cut-offs...
    left_join(lookup, by = "season") %>%
    # ... flag obs that are outside the seasonal 1st - 99th percentile range:
    add_flag((mean < t_mean01 | mean > t_mean99), "outside of seasonal range") %>%
    # flag obs whose slope is outside the 1st - 99th percentile range:
    add_flag(((slope_ahead >= t_slope_behind_99 | slope_behind >= t_slope_behind_99) |
              (slope_ahead <= t_slope_behind_01 | slope_behind <= t_slope_behind_01)), "slope violation")

  return(df)

}

# Adding flags related to sensor specification ranges; these are instances where
# a value exceeds the manufacturer specified limits.

#' @title Add a flag if the value of a parameter exceeds its sensor specification range.
#'
#' @description
#' A function designed to append the 'outside of sensor specification range' flag to
#' a row if the `mean` column in that row exceeds the sensor specification ranges that are
#' set in the `src/qaqc/sensor_spec_thresholds.yml` file. These are instances
#' where a value exceeds the expected ranges for the Poudre based on the sensor
#' manufacturer's specifications.
#'
#' @param df A data frame with a `flag` column.
#'
#' @return A data frame with a `flag` column that has been updated with the
#' 'outside of sensor specification range' flag.
#'
#' @examples
#' add_spec_flag(df = all_data_flagged$`archery-Actual Conductivity`)
#' add_spec_flag(df = all_data_flagged$`boxelder-Temperature`)
#'
#' @seealso [flag_all_data()]

add_spec_flag <- function(df){

  sensor_spec_ranges <- read_yaml("data/qaqc/sensor_spec_thresholds.yml")

  # make this a non yaml solution and add it to the threshold table
  # get the parameter from the parameter column in the df of interest
  parameter_name <- unique(na.omit(df$parameter))
  # Pull the sensor specification range from the yaml file
  sensor_min <- eval(parse(text = sensor_spec_ranges[[parameter_name]]$min))
  sensor_max <- eval(parse(text = sensor_spec_ranges[[parameter_name]]$max))

  df <- df %>%
    # adding sensor range flags
    add_flag(parameter == parameter_name & (mean < sensor_min | mean > sensor_max) & !grepl("outside of sensor specification range", flag),
             "outside of sensor specification range") %>%

    return(df)

}

# Add flags to the `flag` column of a dataframe based on large swaths of suspect data.
# "24hr anomaly flag" is added if more than 50% of the data points in a 24 hour window are flagged.
# "anomaly window" flag is added if the point is included in a 24hr anomaly.
# @param df A dataframe with a `auto_cleaned_flag` column.
# @return A dataframe with a `auto_cleaned_flag` column that has been updated with the relevant large anomaly flags.
# @examples
# add_large_anomaly_flags(df = all_data_flagged$`archery-Actual Conductivity`)
# add_large_anomaly_flags(df = all_data_flagged$`boxelder-Temperature`)

add_suspect_flag_full <- function(df) {
  # these are the flags that we don't want to perform this exercise across
  auto_cleaned_flag_string <- "sonde not employed|missing data|site visit|sv window|sonde burial|sensor biofouling|depth calibration malfunction|sonde unsubmerged|sensor malfunction"

  # Define a function to check if a given 3-hour window has >= 50% fails
  check_3_hour_window_fail <- function(x) {
    sum(x) / length(x) >= 0.5
  }

  df_test <- df %>%
    dplyr::mutate(auto_cleaned_flag_binary = ifelse((is.na(auto_cleaned_flag) | grepl(auto_cleaned_flag_string, auto_cleaned_flag) | auto_cleaned_flag_string == "suspect data"), 0, 1)) %>%
    #arrange(timestamp) %>%
    dplyr::mutate(over_50_percent_fail_window = zoo::rollapply(auto_cleaned_flag_binary, width = 12, FUN = check_3_hour_window_fail, fill = NA, align = "right")) %>%
    add_flag(over_50_percent_fail_window == TRUE & is.na(auto_cleaned_flag), "suspect data")

  return(df_test)

}



#' @title Add flags to the `flag` column of a dataframe based on large swaths of suspect data.
#'
#' @description
#' "24hr anomaly flag" is added if more than 50% of the data points in a 24 hour window are flagged.
#' "anomaly window" flag is added if the point is included in a 24hr anomaly.
#'
#' @param df A data frame with a `flag` column.
#'
#' @return A data frame with a `flag` column that has been updated with the relevant calculated seasonal range flags.
#'
#' @examples
#' add_range_flags(df = all_data_flagged$`archery-Actual Conductivity`)
#' add_range_flags(df = all_data_flagged$`boxelder-Temperature`)
#'
#' @seealso [flag_all_data()]

add_suspect_flag <- function(df) {

  flag_string <- "sonde not employed|missing data|site visit|sv window|suspect data" # include the flag added by this function

  # Define a function to check if a given 3-hour window has >= 50% fails
  check_3_hour_window_fail <- function(x) {
    sum(x) / length(x) >= 0.5
  }

  df_test <- df %>%
    dplyr::mutate(flag_binary = ifelse((is.na(flag) | grepl(flag_string, flag)), 0, 1)) %>%
    #arrange(timestamp) %>%
    dplyr::mutate(over_50_percent_fail_window = ifelse(is.na(over_50_percent_fail_window),
                                                zoo::rollapply(flag_binary, width = 12, FUN = check_3_hour_window_fail, fill = NA, align = "right"),
                                                over_50_percent_fail_window)) %>%
    add_flag(over_50_percent_fail_window == TRUE & is.na(flag), "suspect data")

  return(df_test)

}

add_threshold_lines <- function(plot, plot_data, site_arg, parameter_arg) {

  # pull in threshold data (i don't like that I do this everytime the function is run)
  real_thresholds <- read_csv("data/qaqc/realistic_thresholds.csv", show_col_types = FALSE) %>%
    filter(parameter == parameter_arg)
  sensor_thresholds <- yaml::read_yaml("data/qaqc/sensor_spec_thresholds.yml")[[parameter_arg]]%>% #filter for parameter_arg
    #turn into tibble with min/max
    bind_rows()

  seasonal_thresholds <- read_csv('data/qaqc/seasonal_thresholds_virridy.csv', show_col_types = FALSE) %>%
    #to do: Check to make sure seasonal thresholds csv is not necessary
    #bind_rows(read_csv('data/qaqc/seasonal_thresholds.csv', show_col_type = FALSE),
    distinct(site, parameter, season, .keep_all = TRUE) %>%
    #read_csv("data/qaqc/seasonal_thresholds_virridy.csv", show_col_types = FALSE) %>%
    filter(parameter == parameter_arg,
           site == site_arg)

  # Determine the unique seasons in plot_data
  unique_seasons <- unique(plot_data$season)

  # Filter seasonal_thresholds for the seasons present in plot_data
  seasonal_thresholds <- seasonal_thresholds %>%
    filter(season %in% unique_seasons)

  if (nrow(seasonal_thresholds) > 1) { # make sure this works

    season_1 <- case_when(
      all(c("winter_baseflow", "snowmelt") %in% unique_seasons) ~ "winter_baseflow",
      all(c("snowmelt", "monsoon") %in% unique_seasons) ~ "snowmelt",
      all(c("monsoon", "fall_baseflow") %in% unique_seasons) ~ "monsoon",
      all(c("fall_baseflow", "winter_baseflow") %in% unique_seasons) ~ "fall_baseflow",
      TRUE ~ NA_character_
    )
    season_2 <- case_when(
      all(c("winter_baseflow", "snowmelt") %in% unique_seasons) ~ "snowmelt",
      all(c("snowmelt", "monsoon") %in% unique_seasons) ~ "monsoon",
      all(c("monsoon", "fall_baseflow") %in% unique_seasons) ~ "fall_baseflow",
      all(c("fall_baseflow", "winter_baseflow") %in% unique_seasons) ~ "winter_baseflow",
      TRUE ~ NA_character_
    )

    seasonal_threshold_s1_quantiles <- unname(quantile(c(seasonal_thresholds %>% filter(season == season_1) %>% pull(t_mean01),
                                                         seasonal_thresholds %>% filter(season == season_1) %>% pull(t_mean99)),
                                                       c(0.1, 0.9)))

    seasonal_threshold_s2_quantiles <- unname(quantile(c(seasonal_thresholds %>% filter(season == season_2) %>% pull(t_mean01),
                                                         seasonal_thresholds %>% filter(season == season_2) %>% pull(t_mean99)),
                                                       c(0.1, 0.9)))

    slice_data <- plot_data %>%
      group_by(season) %>%
      slice(1) %>%
      ungroup()

    site_data <- filter(plot_data, site == site_arg)

    if (!all(is.na(site_data$mean))) {
      # Lower bound
      if ( # !is.infinite(min(site_data$mean, na.rm = TRUE)) &
        (min(site_data$mean, na.rm = TRUE) <= seasonal_threshold_s1_quantiles[1] | min(site_data$mean, na.rm = TRUE) <= seasonal_threshold_s2_quantiles[1])) {

        # Xs
        start_x_s1 <- min(slice_data$DT_round)

        transition_date <- slice_data %>%
          filter(season == season_2) %>%
          pull(DT_round)

        end_x_s2 <- ceiling_date(max(plot_data$DT_round), "day")

        #Ys
        y_lower_s1 <- seasonal_thresholds %>%
          filter(season == season_1) %>%
          pull(t_mean01)

        y_lower_s2 <- seasonal_thresholds %>%
          filter(season == season_2) %>%
          pull(t_mean01)

        plot <- plot +
          # season 1
          ## lower bound
          geom_segment(aes(x = start_x_s1, y = y_lower_s1,
                           xend = transition_date, yend = y_lower_s1,
                           color = "Seasonal Min", linetype = "Seasonal Min")) +
          # season 2
          ## lower bound
          geom_segment(aes(x = transition_date, y = y_lower_s2,
                           xend = end_x_s2, yend = y_lower_s2,
                           color = "Seasonal Min", linetype = "Seasonal Min"))

      }

      # Upper bound
      if (# !is.infinite(max(site_data$mean, na.rm = TRUE)) &
        (max(site_data$mean, na.rm = TRUE) >= seasonal_threshold_s1_quantiles[2] | max(site_data$mean, na.rm = TRUE) >= seasonal_threshold_s2_quantiles[2])) {

        # Xs
        start_x_s1 <- min(slice_data$DT_round)


        transition_date <- slice_data %>%
          filter(season == season_2) %>%
          pull(DT_round)

        end_x_s2 <- ceiling_date(max(plot_data$DT_round), "day")

        #Ys
        y_upper_s1 <- seasonal_thresholds %>%
          filter(season == season_1) %>%
          pull(t_mean99)

        y_upper_s2 <- seasonal_thresholds %>%
          filter(season == season_2) %>%
          pull(t_mean99)

        plot <- plot +
          # season 1
          ## lower bound
          geom_segment(aes(x = start_x_s1, y = y_upper_s1,
                           xend = transition_date, yend = y_upper_s1,
                           color = "Seasonal Max", linetype = "Seasonal Max")) +
          # season 2
          ## lower bound
          geom_segment(aes(x = transition_date, y = y_upper_s2,
                           xend = end_x_s2, yend = y_upper_s2,
                           color = "Seasonal Max", linetype = "Seasonal Max"))

      }

      # real thresholds
      real_thresholds_quantiles <- unname(quantile(c(real_thresholds$min, real_thresholds$max), c(0.1, 0.9)))

      if (#!is.infinite(min(site_data$mean, na.rm = TRUE))  &
        (min(site_data$mean, na.rm = TRUE) <= real_thresholds_quantiles[1])) {
        plot <- plot +
          geom_hline(aes(yintercept = real_thresholds$min,
                         color = "Real Min",
                         linetype = "Real"))
      }
      if (#!is.infinite(max(site_data$mean, na.rm = TRUE)) &
        (max(site_data$mean, na.rm = TRUE) >= real_thresholds_quantiles[2])) {
        plot <- plot +
          geom_hline(aes(yintercept = real_thresholds$max,
                         color = "Real Max",
                         linetype = "Real"))
      }

      # sensor thresholds
      sensor_thresholds_quantiles <- unname(quantile(c(sensor_thresholds$min, sensor_thresholds$max), c(0.1, 0.9)))

      # if ((min(site_data$mean, na.rm = TRUE) <= sensor_thresholds_quantiles[1])) {
      #   plot <- plot +
      #     geom_hline(aes(yintercept = sensor_thresholds$mix, # *** this needs to be min
      #                    color = "Sensor Min",
      #                    linetype = "Sensor"))
      # }

      if ((max(site_data$mean, na.rm = TRUE) >= sensor_thresholds_quantiles[2])) {
        plot <- plot +
          geom_hline(aes(yintercept = sensor_thresholds$max,
                         color = "Sensor Max",
                         linetype = "Sensor"))
      }
      return(plot)
    }

  } else if (nrow(seasonal_thresholds) == 1){

    site_data <- filter(plot_data, site == site_arg)

    # Filter seasonal_thresholds for the seasons present in site_data
    seasonal_thresholds <- seasonal_thresholds %>%
      filter(season %in% unique(site_data$season))

    if (!all(is.na(site_data$mean))) {
      # seasonal thresholds
      seasonal_thresholds_quantiles <- unname(quantile(c(seasonal_thresholds$t_mean01, seasonal_thresholds$t_mean99), c(0.1, 0.9)))

      if ((min(site_data$mean, na.rm = TRUE) <= seasonal_thresholds_quantiles[1]) == TRUE){
        plot <- plot +
          geom_hline(aes(yintercept = seasonal_thresholds$t_mean01,
                         color = "Seasonal Min",
                         linetype = "Seasonal"))
      }

      if ((max(site_data$mean, na.rm = TRUE) >= seasonal_thresholds_quantiles[2]) == TRUE) {
        plot <- plot +
          geom_hline(aes(yintercept = seasonal_thresholds$t_mean99,
                         color = "Seasonal Max",
                         linetype = "Seasonal"))
      }

      # real thresholds
      real_thresholds_quantiles <- unname(quantile(c(real_thresholds$min, real_thresholds$max), c(0.1, 0.9)))

      if (min(site_data$mean, na.rm = TRUE) <= real_thresholds_quantiles[1]) {
        plot <- plot +
          geom_hline(aes(yintercept = real_thresholds$min,
                         color = "Real Min",
                         linetype = "Real"))
      }
      if (max(site_data$mean, na.rm = TRUE) >= real_thresholds_quantiles[2]) {
        plot <- plot +
          geom_hline(aes(yintercept = real_thresholds$min,
                         color = "Real Max",
                         linetype = "Real"))
      }

      # sensor thresholds
      sensor_thresholds_quantiles <- unname(quantile(c(sensor_thresholds$min, sensor_thresholds$max), c(0.1, 0.9)))

      # if (min(site_data$mean, na.rm = TRUE) <= sensor_thresholds_quantiles[1]) {
      #   plot <- plot +
      #     geom_hline(aes(yintercept = sensor_thresholds$mix,
      #                color = "Sensor Mine",
      #                linetype = "Sensor"))
      # }

      if (max(site_data$mean, na.rm = TRUE) >= sensor_thresholds_quantiles[2]) {
        plot <- plot +
          geom_hline(aes(yintercept = sensor_thresholds$max,
                         color = "Sensor Max",
                         linetype = "Sensor"))
      }
      return(plot)
    }
    return(plot)
  }
  return(plot)
}

#' @title Add a flag if the sonde was not fully submerged by water.
#'
#' @description
#' A function designed to append the 'unsubmerged' flag to a row if the value
#' in the `relative_depth` column is less than or equal to 0.
#'
#' @param df A data frame with a `flag` column.
#'
#' @return A data frame with a `flag` column that has been updated with the
#' 'sonde unsubmerged' flag.
#'
#' @examples
#' add_unsubmerged_flag(df = all_data_flagged$`archery-Actual Conductivity`)
#' add_unsubmerged_flag(df = all_data_flagged$`boxelder-Temperature`)
#'
#' @seealso [flag_all_data()]

add_unsubmerged_flag <- function(df){

  # create a df of temperature for each site
  depth <- df %>%
    data.table::data.table() %>%
    dplyr::select(DT_round, DT_join, parameter, mean) %>%
    dplyr::filter(parameter == "Depth") %>%
    dplyr::select(DT_join, Depth = mean)

  # add "depth" column to df:
  depth_checked <- df %>%
    dplyr::left_join(., depth, by = "DT_join") %>%
    # If water temperature is freezing, flag all parameters
    add_flag(., Depth <= 0, "sonde unsubmerged") %>%
    # remove the temp column so df is identical in structure to OG df
    dplyr::select(-Depth)

  return(depth_checked)

}

# Add a verification column to the data frame and automatically fail data points that meet certain criteria.
# @param df A dataframe with a `flag` column.
# @return A dataframe with a `verification` column.
# @examples
# add_verification_column(df = all_data_flagged$`archery-Actual Conductivity`)
# add_verification_column(df = all_data_flagged$`boxelder-Temperature`)

add_verification_column <- function(df) {
  df <- df %>%
    mutate(verification = case_when(
      # are there other conditions where a data point should automatically fail?
      # str_detect(flag, "outside sd range") ~ "fail",
      str_detect(flag, "missing data") ~ "fail",
      # need to add more situations where data points fail
      is.na(flag) ~ "pass",
      TRUE ~ NA))
  return(df)
}

alter_verification_column <- function(verification_choice, df_data) {

  df_name <- paste(unique(df_data$site), unique(df_data$parameter), head(df_data$DT_round, 1))

  # PASS/FAIL statements ----
  # these cases are inverses of each other, both are included
  # in the same if statements to reduce confusion

  # PASS ALL/FAIL NONE
  if (verification_choice %in% c("PASS ALL", "FAIL NONE")) {

    altered_df <- df_data %>%
      mutate(mean_verified = if_else(is.na(flag), mean, NA),
             is_verified = TRUE,
             verification_status = "PASS")

    cat("All points for ", df_name, "have PASSED and are accurate information.\n")

    return(altered_df)

  }

  # PASS VALID/FAIL FLAGGED
  if (verification_choice %in% c("PASS VALID", "FAIL FLAGGED")) {

    altered_df <- df_data %>%
      mutate(mean_verified = mean,
             is_verified = TRUE,
             verification_status = if_else(is.na(flag), "PASS", "FAIL")) # pass valid data, fail flagged data

    cat("All valid points for ", df_name, "have PASSED\n")
    cat("All flagged points for ", df_name, "have FAILED\n")

    return(altered_df)

  }

  # PASS FLAGGED, FAIL VALID
  # if all the flagged data is correct, and none of the valid data is correct,
  # none of it should be shown
  if (verification_choice %in% c("PASS FLAGGED", "FAIL VALID")) {

    altered_df <- df_data %>%
      mutate(mean_verified = NA,
             is_verified = TRUE,
             verification_status = if_else(!is.na(flag), "PASS", "FAIL"))

    cat("All valid points for ", df_name, "have FAILED\n")
    cat("All flagged points for ", df_name, "have PASSED.\n")

    return(altered_df)

  }

  # PASS NONE, FAIL ALL
  # This means that flagged data would actually be accurate data, and valid data would
  # actually be inaccurate data. this is the inverse of pass all
  if (verification_choice %in% c("PASS NONE", "FAIL ALL")) {

    altered_df <- df_data %>%
      mutate(mean_verified = if_else(!is.na(flag), mean, NA),
             is_verified = TRUE,
             verification_status = "FAIL")

    cat("All points for ", df_name, "have FAILED and are inaccurate information.\n")

    return(altered_df)

  }

  # skip statements
  if (verification_choice == "SKIP") {

    altered_df <- df_data %>%
      mutate(mean_verified = if_else(is.na(flag), mean, NA),
             is_verified = TRUE, #false?
             verification_status = "SKIP")

    cat("The accuracy of all points for ", df_name, "cannot be determined as\n")
    cat("accurate or inaccurate information at this time.")

    return(altered_df)

  }

  # inspect statements
  if (verification_choice %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME")) {
    stop("You can't inspect here")
  }

  # quit statements ***
  if (verification_choice == "QUIT") {

    cat("Quitting out of verification mode...")
    cat("\n")

    QUIT <<- TRUE

    return(NULL)

  }

}


# to do (j): document this function
api_puller <- function(site, start_dt, end_dt = Sys.time(), api_token, dump_dir, require = NULL) {

  locs <- hv_locations_all(hv_token)

  # make a list of site names

  options(scipen = 999)
  for(i in 1:length(site)){

    site_loc <- locs %>%
      dplyr::mutate(name = tolower(name)) %>%
      dplyr::filter(grepl(site[i], name, ignore.case = TRUE))

    site_loc_list <- site_loc$id

    # Get data for each site location. Note this maps over the entire list of locations,
    # many of which are unlikely to be active during the time you specify. Don't freak out if
    # you see a bunch of '404 Not Found' errors, you're just seeing the list of locations
    # that are not active. The data frame 'alldata' should contain your data from all applicable
    # sites during the time frame indicated. Note that this will take some time (one month of
    # data for 5 sites takes ~10 mins. Be patient!

    # Add date range you are interested in; data are stored in HydroVu in UTC
    # Here, a way to find the most recent download of the data. Use this as the start date to
    # reduce overlapping data

    # tz weirdness
    # utc_start_date <- format(as.POSIXct(start_dt, tz = "UTC") + hours(7), format = "%Y-%m-%d %H:%M:%S")
    #
    # utc_end_date <-   format(as.POSIXct(end_dt, tz = "UTC") + hours(7), format = "%Y-%m-%d %H:%M:%S")

    # doing this fixes the mismatch in date times during the combined_data step - jd
    utc_start_date <- format(as.POSIXct(start_dt, tz = "UTC") + lubridate::hours(0), format = "%Y-%m-%d %H:%M:%S")

    utc_end_date <-   format(as.POSIXct(end_dt, tz = "UTC") + lubridate::hours(0), format = "%Y-%m-%d %H:%M:%S")

    timezone = "UTC"

    # Map over the location ids
    alldata <- site_loc_list %>% purrr::map(~hv_data_id(.,
                                                        start_time = utc_start_date,
                                                        end_time = utc_end_date,
                                                        token = api_token,
                                                        tz = timezone))

    # grab only locations with data (stored as a data frame) / drop 404 errors
    filtered <- purrr::keep(alldata, is.data.frame)

    if(length(filtered) == 0){

      print(paste0("No data at ", site[i], " during this time frame"))

    } else {

      # bind lists together (now that all are dataframes, we can just collate quickly)
      one_df <- dplyr::bind_rows(filtered) %>%
        dplyr::rename(id = Location,
                      parameter = Parameter,
                      units = Units) %>%
        dplyr::left_join(., site_loc, by = "id") %>%
        dplyr::mutate(site = tolower(site[i])) %>%
        dplyr::select(site, id, name, timestamp, parameter, value, units)

      # if site contains both a csu and a virridy sonde, split them up:

      if(site[i] %in% c("Timberline", "Prospect", "Archery")){

        try(virridy_df <- one_df %>%
              filter(grepl("virridy", name, ignore.case = TRUE)) %>%
              mutate(site = paste0(site[i], " virridy")))

        try(readr::write_csv(virridy_df,
                             paste0(dump_dir, "/", site[i], " virridy_", stringr::str_replace(stringr::str_replace(substr(end_dt, 1, 16), "[:\\s]", "_"), ":", ""), ".csv")))

        csu_df <- one_df %>%
          filter(!grepl("virridy", name, ignore.case = TRUE))

        readr::write_csv(csu_df,
                         paste0(dump_dir, "/", site[i], "_", stringr::str_replace(stringr::str_replace(substr(end_dt, 1, 16), "[:\\s]", "_"), ":", ""), ".csv"))

      } else {

        # otherwise, save the full data set

        readr::write_csv(one_df,
                         paste0(dump_dir, "/", site[i], "_", stringr::str_replace(stringr::str_replace(substr(end_dt, 1, 16), "[:\\s]", "_"), ":", ""), ".csv"))
      }
    }
  }

}

append_historical_api_data <- function(hist_dir, inc_dir) {

  # list of site names
  site_names <- list(
    "Archery",
    "Boxelder",
    "Legacy",
    "Lincoln",
    "Prospect",
    "River Bluffs",
    "Tamasag",
    "Timberline"
  )

  # get the full file names
  inc_full_file_names <- list.files(inc_dir, pattern = "*.csv", full.names = TRUE)
  hist_full_file_names <- list.files(hist_dir, pattern = "*.csv", full.names = TRUE)

  # find the files that match using the site names
  walk(site_names, function(site_name) {

    # find the index of the matching site names in the file names
    inc_site_name_full_path <- grepl(site_name, inc_full_file_names, ignore.case = TRUE)
    hist_site_name_full_path <- grepl(site_name, hist_full_file_names, ignore.case = TRUE)

    # get the file names
    inc_site_path <- inc_full_file_names[inc_site_name_full_path]
    hist_site_path <- hist_full_file_names[hist_site_name_full_path]

    # read in the files
    inc_file <- read_csv(inc_site_path)
    hist_file <- read_csv(hist_site_path)

    # combine the files
    new_file <- bind_rows(hist_file, inc_file) %>%
      distinct()

    # write the file (this will replace the old file)
    write_csv(new_file, paste0(hist_dir, site_name, "_historical.csv"))
  })

  # delete the files in the inc_dir
  walk(inc_full_file_names, unlink)

}




basic_plot <- function(){
  max_1 <- max(wq_tl[[parameters[1]]], na.rm = T)

 param_1 <- ggplot(wq_tl, aes(x = DT_round, y = .data[[parameters[1]]], ymin = 0, ymax = .data[[parameters[1]]]))+
   geom_ribbon(color = "blue", fill = "blue" )+
   theme_bw()+
   labs(x = "Date", y = parameters[1])

 param_2 <- ggplot(wq_tl, aes(x = DT_round, y = .data[[parameters[2]]]))+
   geom_path(color = "red")+
   theme_bw()+
   labs(x = "Date", y = parameters[2])

 ggarrange(param_1, param_2,ncol = 2, nrow = 1)

}

check_incoming_api_dir <- function(incoming_dir, archive_dir) {
  # Check if data/api/incoming_api_data exists
  if (dir.exists(incoming_dir)) {
    # Check if incoming directory is empty
    if (length(list.files(incoming_dir)) == 0) {
      print(paste0(incoming_dir, " exists and is empty."))
      print("Incoming API data directory is configured properly.")
    } else {
      print(paste0(incoming_dir, " exists but is not empty..."))
      print("Please ensure previous workflow ran properly...")
      stop("Pipeline halted due to non-empty incoming directory.")
    }
  } else {
    print(paste0(incoming_dir, " does not exist..."))
    print("Creating incoming directory...")
    dir.create(incoming_dir, recursive = TRUE)
    print("Incoming directory created.")
  }

  # Check if data/api/archive_api_data exists
  if (dir.exists(archive_dir)) {
    print(paste0(archive_dir, " exists."))
    print("Directory is configured properly.")
  } else {
    print(paste0(archive_dir, " does not exist..."))
    print("Creating archive directory...")
    dir.create(archive_dir, recursive = TRUE)
    print("Archive directory created.")
  }
}

clean_directories <- function() {
  pre_dir_path <- pre_verification_path
  int_dir_path <- intermediary_path
  ver_dir_path <- verified_path

  pre_dir_names <- list.files(pre_dir_path)
  int_dir_names <- list.files(int_dir_path)
  ver_dir_names <- list.files(ver_dir_path)

  # Pre-verification Directory
  # keep if not in Verified Directory, else if in Verified Directory delete from this Directory
  for (i in pre_dir_names) {
    if (length(ver_dir_names) != 0 & i %in% ver_dir_names) {
      file.remove(paste0(pre_dir_path, i))
      cat("removed file ", i, " from ",  pre_dir_path, "\n")
    } else {
      next
    }
  }

  # Intermediary Directory
  # keep if skips in the data OR any(!is_verified), else move to Verified Directory
  for (i in int_dir_names) {
    df <- readRDS(paste0(int_dir_path, i))
    if (any(df$verification_status == 'SKIP') | any(!df$is_verified)) {
      next
    } else if (all(df$verification_status != 'SKIP') & all(df$is_verified)) {
      source_path <- paste0(int_dir_path, i)
      destination_path <- paste0(ver_dir_path, i)
      file.copy(source_path, destination_path, overwrite = TRUE)
      file.remove(source_path)
      cat("moved file ", i, " from ", int_dir_path, "to", ver_dir_path, "\n")
    }
  }

  # Verified Directory
  # final destination for data, check that the data that is here should be here
  for (i in ver_dir_names) {
    df <- readRDS(paste0(ver_dir_path, i))
    if (any(df$verification_status == 'SKIP') | any(!df$is_verified)) {
      cat(i, "should not be in verified directory\n")
      stop("THERE IS DATA WITH SKIPS OR NON-VERIFIED DATA IN THE VERIFIED DIRECTORY")
    } else {
      cat("data in verified directory is correct\n")
      next
    }
  }
}

clear_incoming_data_dir <- function(incoming_dir, archive_dir, require = NULL) {
  # Check if the previous step in the targets pipeline ran
  # this is done via the require arg

  # List files in the incoming and archive directories
  incoming_files <- list.files(incoming_dir, full.names = FALSE)
  archive_files <- list.files(archive_dir, full.names = FALSE)

  # Find files that are not already in the archive directory
  files_to_copy <- setdiff(incoming_files, archive_files)

  # Copy only the files that are not already in the archive directory
  if (length(files_to_copy) > 0) {
    for (file in files_to_copy) {
      full_file_name <- file.path(incoming_dir, file)
      file.copy(full_file_name, archive_dir)
      print(paste0(file, " has been moved to archive API data folder."))
    }
    print("Files have been copied from the incoming directory to the archive directory.")
  } else {
    print("All files are already present in the archive directory. Nothing to copy.")
  }
Sys.sleep(5)

# refresh archive_files now that the incoming_files have been copied over
archive_files <- list.files(archive_dir, full.names = FALSE)
  # Check if all files from incoming directory are now in the archive directory
  if (all(incoming_files %in% archive_files)) {
    print("All files in the incoming directory have been successfully copied to the archive directory.")
  } else {
    print("Not all files from the incoming directory have been successfully copied to the archive directory.")
    # Should this halt the pipeline?
    # Right now it seems to always print this and I haven't figured out why so
    # I don't think that it should not halt the pipeline in the state that its in right now -JD.
  }

  # Delete the copied files from the incoming directory
  if (length(files_to_copy) > 0) {
    for (file in files_to_copy) {
      full_file_name <- file.path(incoming_dir, file)
      file.remove(full_file_name)
    }
    print("Copied files have been removed from the incoming directory.")
  }

  # Delete any files in the incoming directory
  for (file in list.files(incoming_dir, full.names = TRUE)) {
    file.remove(file)
  }
  print("All files removed from incoming directory.")
}

#' @title Combine historical and incoming data
#'
#' @description
#' A function that combines the subset of historically flagged data and creates historical columns.
#' NOTE: historical columns should be a part of the historically flagged data, but it is not at this point.
#'
#' @param incoming_data_list A list of dataframes with the incoming data.
#'
#' @param historical_data_list A list of dataframes with the historical data.
#'
#' @return A list of dataframes with the combined historical and incoming data.
#'
#' @examples
#' combine_hist_inc_data(incoming_data_list = incoming_data_list, historical_data_list = historical_data_list)

combine_hist_inc_data <- function(incoming_data_list, historical_data_list) {

  # Get the matching index names
  matching_indexes <- intersect(names(incoming_data_list), names(historical_data_list))

  # Get the last 24 hours of the historically flagged data based on earliest
  # DT listed in the api pull for that site/data combo
  last_24_hours <- map(historical_data_list,
                      function(data) {
                        data %>%
                          # most recent day of already-flagged data
                          filter(DT_round >= ymd_hms(max(DT_round) - hours(24), tz = "MST")) %>%
                          # mark it as "historic"
                          mutate(historical = TRUE,
                                 # ensure last_site_visit column has proper DT:
                                 last_site_visit = force_tz(last_site_visit, tzone = "MST"))
                      })

  # bind summarized_incoming_data and last_24_hours together
  combined_hist_inc_data <- map(matching_indexes, function(index) {
    last_24_hours[[index]] %>%
      select(-c(t_mean01, # why is this being removed? -jd
                t_mean99,
                t_slope_behind_01,
                t_slope_behind_99,
                t_sd_0199)) %>%
      bind_rows(., mutate(incoming_data_list[[index]], historical = FALSE)[-1,]) %>% # why are we removing the first row? -jd # there is a datetime issue on this bind
      # make sure new data has info about last site visit, etc. filled out
      fill(c(sonde_employed, last_site_visit, sensor_malfunction))

  }) %>%
    set_names(matching_indexes)

  return(combined_hist_inc_data)
}



#@param folder_name: input the raw photo folder that is named the site you wish to compile
compile_files <- function(folder_name){
  site <- str_extract(folder_name, "tamasag|legacy|timberline|prospect|boxelder|archery|riverbluffs")
  photo_files <- list.files(path = folder_name, full.names = TRUE, recursive = TRUE)
  #files in folder
  new_folder_files <- list.files(path = paste0('data/timelapse_photos/2023_compiled/',site), full.names = TRUE)

  photo_renamer <- function(file) {
    #grab dt from file
    dt <- read_exif(path = file,
                    tags = c("DateTimeOriginal")) %>%
      pull(DateTimeOriginal) %>%
      parse_date_time("YmdHMS", tz="MST")%>%
      format("%Y%m%d_%H%M")
    #create new file name from dt
    new_file_name <- paste0('data/timelapse_photos/2023_compiled/',site,"/", dt,'.JPG')

    #check to see if this file is already in the folder
    if(new_file_name %nin% new_folder_files){
      #if it is not, copy it over
      file.copy(file,to =  new_file_name)
    }

  }
  map(photo_files, photo_renamer)
  print(paste0("Finished ", folder_name))

}

## Photos

# Goals:
# Download all user created photos ( upstream, downstream, clarity, filter and other pictures)
# Label according to site, date, description in the format site_YYYYMMDD_descriptor.jpg
# Only download photos which have not yet been downloaded




download_pictures <- function(){
  #source to grab all notes cleaned
  source("src/load_mWater_notes.R")

  all_notes_cleaned <- load_mWater_notes()
  # Find all the downloaded pictures
  all_file_names <- tolower(list.files(path = "data/field_pics/", recursive = TRUE))
  #grab notes
  sampling_photos <- all_notes_cleaned%>%
    #grab needed columns
    select(site, start_dt,photos_downloaded,upstream_pic,downstream_pic,clarity,filter_pic,other_pic,other_pic_descriptor)%>%
    mutate(
      #correct names if it is in our upper sites (upper case acronyms)
      site = tolower(site),
      #Date format for pictures
      yyyymmdd = format(start_dt, "%Y%m%d"),
      #create filenames ONLY if there is a URL associated with the site visit
      upstream_filename = case_when(
        !is.na(upstream_pic) ~ paste0(site, "_", yyyymmdd, "_upstream.jpg"),
        TRUE ~ NA_character_
      ),
      downstream_filename = case_when(
        !is.na(downstream_pic) ~ paste0(site, "_", yyyymmdd, "_downstream.jpg"),
        TRUE ~ NA_character_
      ),
      clarity_filename = case_when(
        !is.na(clarity) ~ paste0(site, "_", yyyymmdd, "_clarity.jpg"),
        TRUE ~ NA_character_
      ),
      filter_filename = case_when(
        !is.na(filter_pic) ~ paste0(site, "_", yyyymmdd, "_filter.jpg"),
        TRUE ~ NA_character_
      ),
      # check to see if the photo is already downloaded to folder
      upstream_downloaded = case_when(
        is.na(upstream_filename) ~ NA,
        upstream_filename %in% all_file_names ~ TRUE,
        TRUE ~ FALSE
      ),
      clarity_downloaded = case_when(
        is.na(clarity_filename) ~ NA,
        clarity_filename %in% all_file_names ~ TRUE,
        TRUE ~ FALSE
      ),
      filter_downloaded = case_when(
        is.na(filter_filename) ~ NA,
        filter_filename %in% all_file_names ~ TRUE,
        TRUE ~ FALSE
      ),
      downstream_downloaded = case_when(
        is.na(downstream_filename) ~ NA,
        downstream_filename %in% all_file_names ~ TRUE,
        TRUE ~ FALSE
      )
    )
# basic path to field pics
path <- "data/field_pics/"


  # loop thru dataset and download the photo ONLY if it is not yet downloaded and not NA
  for (i in 1:nrow(sampling_photos)) {
    if (!is.na(sampling_photos$upstream_downloaded[i]) && !sampling_photos$upstream_downloaded[i]) {
      print(sampling_photos$upstream_filename[i])
      download.file(sampling_photos$upstream_pic[i], destfile = paste0(path,sampling_photos$upstream_filename[i]))
    }

    if (!is.na(sampling_photos$downstream_downloaded[i]) && !sampling_photos$downstream_downloaded[i]) {
      print(sampling_photos$downstream_filename[i])
      download.file(sampling_photos$downstream_pic[i], destfile = paste0(path, sampling_photos$downstream_filename[i]))
    }
    if (!is.na(sampling_photos$clarity_downloaded[i]) && !sampling_photos$clarity_downloaded[i]) {
      print(sampling_photos$clarity_filename[i])
      download.file(sampling_photos$clarity[i], destfile = paste0(path, sampling_photos$clarity_filename[i]))
    }
    if (!is.na(sampling_photos$filter_downloaded[i]) && !sampling_photos$filter_downloaded[i]) {
      print(sampling_photos$filter_filename[i])
      download.file(sampling_photos$filter_pic[i], destfile = paste0(path, sampling_photos$filter_filename[i]))
    }
  }

  #grab notes for sites with other pictures
  other_photos <- all_notes_cleaned%>%
    #grab needed columns
    select(site, start_dt,other_pic,other_pic_descriptor)%>%
    #get rid of instances with no other pics
    filter(!is.na(other_pic))%>%
    mutate(
      site = tolower(site),
      #Date format for pictures
      yyyymmdd = format(start_dt, "%Y%m%d"),
      # separate multiple URLs in other pic column
      other_pic_sep = str_split(other_pic, "; "),
      #seperate multiple descriptors in the descriptor column
      other_descriptor_sep = str_split(other_pic_descriptor, ","))%>%
    #for rows with multiple pictures, create a new row for each picture

    unnest(cols = c(other_pic_sep, other_descriptor_sep))%>%
    #remove excess columns and rename sep columns to match old columns
    select(site, start_dt,yyyymmdd, other_pic = other_pic_sep, other_pic_descriptor = other_descriptor_sep)%>%
    # make descriptor lower case and remove any spaces in the name
    mutate(other_pic_descriptor = tolower(str_replace_all(other_pic_descriptor, " ", "")),
           other_filename = case_when(!is.na(other_pic) ~ paste0(site, "_", yyyymmdd, "_", other_pic_descriptor, ".jpg")),
           # Check to see if photo has already been downloaded
           other_downloaded = case_when(
             is.na(other_filename) ~ NA,
             other_filename %in% all_file_names ~ TRUE,
             TRUE ~ FALSE
           ))

  # loop thru dataset and download the photo ONLY if it is not yet downloaded and not NA
  for (i in 1:nrow(other_photos)) {
    if (!is.na(other_photos$other_downloaded[i]) && !other_photos$other_downloaded[i]) {
      print(other_photos$other_filename[i])
      download.file(other_photos$other_pic[i], destfile = paste0(path,other_photos$other_filename[i]))
    }}

cat("\nAll Available Pictures Downloaded\n")
}

#download_pictures()
#RUN TO DOWNLOAD NEW PICTURES
# It takes about 2-5 minutes to download ~25-50 photos
# Sometimes the request to mWater time out, just re run the function below if that happens

## Determining uploads

#This function looks at the user inputs for calibration report collect and logs collected.
#Based on these inputs, it looks at all the uploaded logs or calibration reports
#and will print out what logs are missing and who to contact to get those files uploaded.




files_missing <- function(){

  `%nin%` = Negate(`%in%`)
  # #source clean mwater script for all notes cleaned
  #
  # source("src/mWater_collate/load_mWater_notes.R")

  #grab context metadata
  site_meta <- read_csv("data/metadata/water_sampling_sites.csv",show_col_types = FALSE)%>%
    select(site = site_code, Site_Name, site_label_rmrs)
  # sort for sites in upper network (ie. acronyms rather than street names)
  upper_sites <- read_csv("data/metadata/water_sampling_sites.csv",show_col_types = FALSE)%>%
    filter(watershed != "CLP  Mainstem-Fort Collins")%>%
    #this is to help match with user input
    mutate(site_code = tolower(site_code))

  field_season <- year(Sys.Date())

  #grab cal reports from folder
  cal_reports_simple <- str_extract(list.files(path = "data/calibration_reports/"), ".*_\\d{8}" )%>%tolower()
  logs_simple <- str_extract(list.files(path = paste0("data/sensor_data/", field_season), recursive = TRUE), "\\w+_\\d{8}_(vulink|troll)")%>%tolower()




  #grab sensor notes that have logs or cal reports that should be  downloaded
  sensor_files <- all_notes_cleaned%>%
    filter(year(DT_round) == field_season)%>%
    filter(grepl("Sensor",visit_type, ignore.case = TRUE))%>%
    filter(cal_report_collected|log_downloaded)%>%
    select(site, crew, start_DT,end_dt, cal_report_collected, cals_performed, log_downloaded, log1_type,log1_mmdd,  log2_type, log2_mmdd)%>%
    mutate(
      #make all site names lower
      site = tolower(site),
      # Create basis for calibration report name
      # this will be used to check for calibration report in data files and then b
      cal_report_name = case_when(cal_report_collected == TRUE ~ paste0(site, "_", format(start_DT, "%Y%m%d")),
                                  cal_report_collected == NA ~ NA),
      full_cal_name = case_when(cal_report_collected == TRUE ~ paste0(site, "_", format(end_dt, "%Y%m%d_%H%M_mst")),
                                cal_report_collected == NA ~ NA),
      log1_mmdd = case_when(nchar(as.character(log1_mmdd)) == 3 ~ paste0("0",log1_mmdd),
                            TRUE ~ as.character(log1_mmdd)),
      log1_type = case_when( grepl("aquatroll", log1_type,ignore.case = TRUE) ~  "troll",
                             grepl("at", log1_type,ignore.case = TRUE) ~  "troll",
                             grepl("aqua troll", log1_type,ignore.case = TRUE) ~  "troll",
                             grepl("vulink", log1_type,ignore.case = TRUE) ~  "vulink",
                             grepl("vu link", log1_type,ignore.case = TRUE) ~  "vulink",
                             TRUE ~ tolower(log1_type)),
      log2_type = case_when( grepl("aquatroll", log2_type,ignore.case = TRUE) ~  "troll",
                             grepl("at", log2_type,ignore.case = TRUE) ~  "troll",
                             grepl("aqua troll", log2_type,ignore.case = TRUE) ~  "troll",
                             grepl("vulink", log2_type,ignore.case = TRUE) ~  "vulink",
                             grepl("vu link", log2_type,ignore.case = TRUE) ~  "vulink",
                             TRUE ~ tolower(log2_type)),
      log2_mmdd = case_when(nchar(as.character(log2_mmdd)) == 3 ~ paste0("0",log2_mmdd),
                            TRUE ~ as.character(log2_mmdd)),
      log1_user_error = case_when( is.na(log_downloaded)~ FALSE,
                                   log_downloaded == FALSE ~ FALSE,
                                   is.na(log1_mmdd) | is.na(log1_type) ~ TRUE,
                                   TRUE ~ FALSE ),
      log2_user_error = case_when(is.na(log_downloaded)~ FALSE,
                                  log_downloaded == FALSE ~ FALSE,
                                  is.na(log2_mmdd) & is.na(log2_type) ~ FALSE,
                                  is.na(log2_mmdd) | is.na(log2_type) ~ TRUE,
                                  TRUE ~ FALSE),
      log1_name = case_when(!(is.na(log1_mmdd) | is.na(log1_type)) ~ paste0(site,"_",year(start_DT),log1_mmdd,
                                                                            "_", format(start_DT, "%Y%m%d"), "_", log1_type),
                            (is.na(log1_mmdd) | is.na(log1_type))  ~ NA),
      log2_name = case_when(!(is.na(log2_mmdd) | is.na(log2_type)) ~ paste0(site,"_",year(start_DT),log2_mmdd,
                                                                            "_", format(start_DT, "%Y%m%d"), "_", log2_type),
                            (is.na(log2_mmdd) | is.na(log2_type))  ~ NA),
      log_missing1 = case_when(
        is.na(log_downloaded)| log_downloaded == FALSE ~ FALSE,
        is.na(log1_name) ~ FALSE,
        log1_name %nin% logs_simple ~ TRUE,
        TRUE ~ FALSE
      ),
      log_missing2 = case_when(
        is.na(log_downloaded)| log_downloaded == FALSE ~ FALSE,
        is.na(log2_name) ~ FALSE,
        log2_name %nin% logs_simple ~ TRUE,
        TRUE ~ FALSE
      ),
      log_missing = case_when(
        log_missing1 | log_missing2 ~ TRUE,
        TRUE ~ FALSE
      ),
      cal_missing = case_when(
        is.na(cal_report_collected) ~ FALSE,
        cal_report_name %nin% cal_reports_simple ~ TRUE,
        TRUE ~ FALSE)
      )



  for (i in 1:nrow(sensor_files)) {

    #if log missing print out missing logs or cal reports
    if(sensor_files$log_missing[i]){
      cat("\nLog Missing: ", sensor_files$log1_name[i], " and/or ", sensor_files$log2_name[i], "\nContact: ", sensor_files$crew[i], "\n")

    }
    #if log missing print out missing logs or cal reports
    if(sensor_files$cal_missing[i]){
      cat("\nCal Missing: ", sensor_files$full_cal_name[i]," \nContact: ", sensor_files$crew[i], "\n")
    }


  }
  cat("\nFile Check Complete")

}


#' @title Add malfunction flags to a data frame based on field notes and HydroVu.
#'
#' @description
#' This function adds the 'sensor malfunction' flag to a data frame by the dates,
#' sites, and parameters that members of the lab know contain erroneous data due
#' to sensor malfunctions. Note that this flag is used in two instances: first
#' when removing erroneous data from our statistics calculations and again during
#' the actual flagging step in the QAQC pipeline.
#'
#' @param df A data frame with a `flag` column.
#'
#' @param malfunction_records The malfunction records pulled from mWater
#'
#' @return A data frame with a `flag` column that has been updated with the
#' 'sensor malfunction' flag.
#'
#' @examples
#' add_range_flags(df = all_data_flagged$`archery-Actual Conductivity`)
#' add_range_flags(df = all_data_flagged$`boxelder-Temperature`)
#'
#' @seealso [grab_mWater_malfunction_notes()]

fix_calibration <- function(df){

  # Filter records for relevant site-param information
  df_site <- unique(df$site)
  df_parameter <- unique(df$parameter)

  # Depth calibration requires its own fix. It's hacky and I don't like it, but
  # the way that the calibration report reports pressure calibrations makes it
  # impossible to actually back calibrate. Therefore, I'm just assuming that level
  # doesn't actually change after bad calibrations, and hard code it so that the
  # first "bad" depth is forced to equal the last "good" depth and I offset all
  # subsequent "bad" depth values by the difference between them.

  if(df_parameter == "Depth"){

    if(!"Depth" %in% df$parameter){

      nope <- df %>% mutate(relative_depth = NA)

      return(nope)

    }

    if("Depth" %in% df$parameter & "archery" %in% df$site){

      #df <- all_data_flagged[["archery-Depth"]]

      site_depth <- df %>%
        mutate(relative_depth = ifelse(DT_round >= as_datetime('2022-05-21 15:45:00', "MST") & DT_round <= as_datetime('2022-05-24 15:45:00', "MST"),
                                       mean +
                                         abs(
                                           dplyr::filter(df, as.character(DT_round) == ('2022-05-24 15:15:00'))$mean -
                                             dplyr::filter(df, as.character(DT_round) == ('2022-05-24 15:45:00'))$mean
                                         ),
                                       mean))
      #return(archery_depth)

    } else if("Depth" %in% df$parameter & "timberline" %in% df$site){

      # df <- all_data_flagged[["timberline-Depth"]]

      site_depth <- df %>%
        dplyr::mutate(relative_depth = ifelse(year == "2022" & DT_round <= lubridate::as_datetime('2022-04-07 17:00:00', "MST"),
                                              as.numeric(mean) +
                                                abs(
                                                  dplyr::filter(df, as.character(DT_round) == ('2022-04-07 16:15:00'))$mean -
                                                    dplyr::filter(df, as.character(DT_round) == ('2022-04-07 17:15:00'))$mean
                                                ),
                                              as.numeric(mean)))
      # return(timberline_depth)

    } else if ("Depth" %in% df$parameter & "legacy" %in% df$site) {

      # df <- all_data_flagged[["legacy-Depth"]]

      site_depth <- df %>%
        mutate(relative_depth = ifelse(DT_round >= as_datetime('2022-04-06 06:00:00', "MST") & DT_round <= as_datetime('2022-04-12 09:15:00', "MST"),
                                       as.numeric(mean) +
                                         abs(
                                           dplyr::filter(df, as.character(DT_round) == ('2022-04-12 09:15:00'))$mean -
                                             dplyr::filter(df, as.character(DT_round) == ('2022-04-12 09:30:00'))$mean
                                         ),
                                       mean))

      site_depth <-  site_depth %>%
        mutate(relative_depth = ifelse(DT_round >= as_datetime('2022-07-08 17:00:00', "MST") & DT_round <= as_datetime('2022-07-12 09:00:00', "MST"),
                                       relative_depth +
                                         abs(
                                           dplyr::filter(site_depth, as.character(DT_round) == ('2022-07-08 14:15:00'))$relative_depth -
                                             dplyr::filter(site_depth, as.character(DT_round) == ('2022-07-08 17:00:00'))$relative_depth
                                         ),
                                       relative_depth))

      site_depth <-  site_depth %>%
        mutate(relative_depth = ifelse(DT_round >= as_datetime('2022-07-22 11:30:00', "MST") & DT_round <= as_datetime('2022-07-25 14:15:00', "MST"),
                                       relative_depth +
                                         abs(
                                           dplyr::filter(site_depth, as.character(DT_round) == ('2022-07-22 09:45:00'))$relative_depth -
                                             dplyr::filter(site_depth, as.character(DT_round) == ('2022-07-22 11:30:00'))$relative_depth
                                         ),
                                       relative_depth))



      #return(legacy_depth)

    } else if ("Depth" %in% df$parameter & "tamasag" %in% df$site) {

      # df <- all_data_flagged[["tamasag-Depth"]]

      site_depth <- df %>%
        mutate(relative_depth = ifelse(DT_round <= "2022-04-24 07:15:00" & year == "2022",
                                       as.numeric(mean) +
                                         abs(
                                           dplyr::filter(df, as.character(DT_round) == ('2022-04-24 07:15:00'))$mean -
                                             dplyr::filter(df, as.character(DT_round) == ('2022-04-24 07:30:00'))$mean
                                         ),
                                       mean))

      # return(tamasag_depth)

    } else {

      site_depth <- df %>%
        mutate(relative_depth = mean,
               cal_fix = NA)
    }

    depth_flagged <- site_depth %>%
      dplyr::mutate(raw = mean,
                    mean = relative_depth,
                    cal_fix = ifelse(raw != mean, "calibration fix", NA)) %>%
      dplyr::select(-relative_depth)

    return(depth_flagged)

  }

  # PLACEHOLDER UNTIL OTHER CALIBRATIONS ARE GOOD:
 return(df %>%
    mutate(raw = mean,
           cal_fix = NA))

}

  # # For non-depth parameters, we can refer to the calibration report:
  # bad_cal <- readxl::read_excel("data/calibration_error_log.xlsx") %>%
  #   dplyr::mutate(start_DT = as.character(lubridate::as_datetime(start_DT)),
  #                 end_DT = as.character(lubridate::as_datetime(end_DT)),
  #                 #report_to_correct = as.character(lubridate::as_datetime(report_to_correct))
  #   )
  #
  # bad_cal_records_filtered <- bad_cal %>%
  #   filter(site == df_site) %>%
  #   filter(grepl(df_parameter, parameter, ignore.case = TRUE)) %>%
  #   mutate(end_DT = ifelse(is.na(end_DT), ymd_hms("9999-12-31 23:59:59", tz = "MST"), end_DT)) %>%
  #   mutate(end_DT = as.character(as.POSIXct(end_DT, tz = "MST"))) %>%
  #   rowid_to_column()
  #
  # # If there are no bad calibrations listed for that site-param, return original
  # # dataframe, filling our updated "mean" column with the old unmodified values:
  # if(nrow(bad_cal_records_filtered == 0)){
  #
  #   df <- df %>%
  #     mutate(raw = mean,
  #            mean = mean)
  #
  #   return(df)
  #
  # } else {
  #
  #   cal_tabler <- function(cal_files){
  #
  #     #cal_files <- list.files("data/calibration_reports")[3] # for testing
  #
  #     cal <- read_html(paste0(getwd(), "/data/calibration_reports/", cal_files)) %>%
  #       html_nodes("div") %>%
  #       html_text() %>%
  #       as_tibble()
  #
  #     rdo <- cal %>% filter(grepl("RDO", value)) %>% pull() %>% str_replace_all(., " ", "") %>% tolower()
  #
  #     ph_orp <- cal %>% filter(grepl("pH/ORP", value)) %>% pull() %>% str_replace_all(., " ", "") %>% tolower()
  #
  #     conductivity <- cal %>% filter(grepl("Conductivity",value)) %>% pull() %>% str_replace_all(., " ", "") %>% tolower()
  #
  #     if(length(cal %>% filter(grepl("Turbidity",value)) %>% pull() %>% str_replace_all(., " ", "") %>% tolower()) != 0){
  #
  #       turbidity <- cal %>% filter(grepl("Turbidity",value)) %>% pull() %>% str_replace_all(., " ", "") %>% tolower()
  #
  #     } else {
  #
  #       turbidity <- "No Turbidity Sensor"
  #
  #     }
  #
  #     # Always the last sensor when depth is available:
  #     depth <- ifelse(str_detect(cal %>% .[nrow(.),] %>% pull() %>% str_replace_all(., " ", "") %>% tolower(), "pressure"),#"psireferencedepth"),
  #                     cal %>% .[nrow(.),] %>% pull() %>% str_replace_all(., " ", "") %>% tolower(),
  #                     "No Depth Sensor")
  #
  #     time_mst <- paste0(str_sub(str_match(cal_files, "(\\d+)_mst")[, 2:1][1], 1, 2), ":",
  #                        str_sub(str_match(cal_files, "(\\d+)_mst")[, 2:1][1], 3, 4))
  #     #str_sub(cal_files, -13, -12),":", str_sub(cal_files, -11, -10))
  #
  #     date <- str_match(cal_files, "^[^_]+_([0-9]{8})_")[, 2]
  #
  #     #str_sub(cal_files, -22, -19),"-", str_sub(cal_files, -18, -17),"-", str_sub(cal_files, -16, -15))
  #
  #     cal_table <- tibble(site = sub("\\_.*", "", cal_files) %>% tolower(),
  #
  #                         DT = ymd_hm(paste(date, time_mst, tz = "MST")),
  #
  #                         # Dissolved Oxygen
  #                         rdo_cal_date = as.character(mdy(str_match(rdo, "lastcalibrated\\s*(.*?)\\s*calibrationdetails")[,2])),
  #                         rdo_slope = str_match(rdo, "slope\\s*(.*?)\\s*offset")[,2],
  #                         rdo_offset = str_match(rdo, "offset\\s*(.*?)\\s*mg/l")[,2],
  #                         rdo_100 = str_match(rdo, "premeasurement\\s*(.*?)\\s*%satpost")[,2],
  #                         rdo_conc = str_match(rdo, "concentration\\s*(.*?)\\s*mg/lpremeasurement")[,2],
  #                         rdo_temp = str_match(rdo, "temperature\\s*(.*?)\\s*°c")[,2],
  #                         rdo_pressure = str_match(rdo, "pressure\\s*(.*?)\\s*mbar")[,2],
  #
  #                         # pH
  #                         ph_cal_date = as.character(mdy(str_match(ph_orp, "lastcalibrated\\s*(.*?)\\s*calibrationdetails")[,2])),
  #                         ph_slope_pre = str_match(ph_orp, "offset1slope\\s*(.*?)\\s*mv/ph")[,2],
  #                         ph_offset_pre = str_match(ph_orp, "mv/phoffset\\s*(.*?)\\s*mvslopeandoffset2")[,2],
  #                         ph_slope_post = str_match(ph_orp, "offset2slope\\s*(.*?)\\s*mv/ph")[,2],
  #                         ph_offset_post = str_match(ph_orp, paste0(ph_slope_post,"mv/phoffset\\s*(.*?)\\s*mvorporp"))[,2],
  #                         # Sometimes, the post value can actually be in the high 6 pH... therefore the post measurement regex matching text is conditional
  #                         ph_7_nice = str_sub(str_match(ph_orp, "postmeasurementph7\\s*(.*?)\\s*mvcal")[,2], 10, nchar(str_match(ph_orp, "postmeasurementph7\\s*(.*?)\\s*mvcal")[,2])),
  #                         ph_7_high = str_sub(str_match(ph_orp, "postmeasurementph8\\s*(.*?)\\s*mvcal")[,2], 10, nchar(str_match(ph_orp, "postmeasurementph8\\s*(.*?)\\s*mvcal")[,2])),
  #                         ph_7_low = str_sub(str_match(ph_orp, "postmeasurementph6\\s*(.*?)\\s*mvcal")[,2], 10, nchar(str_match(ph_orp, "postmeasurementph6\\s*(.*?)\\s*mvcal")[,2])),
  #                         ph_7 = ifelse(!is.na(ph_7_nice), ph_7_nice,
  #                                       ifelse(!is.na(ph_7_high), ph_7_high, ph_7_low)),
  #
  #                         # ORP
  #                         #Newly encountered thing: sometimes the calibration report calls the ORP standard Zobell's, sometimes it's just called "ORP Standard":
  #                         orp_offset = ifelse(is.na(str_match(ph_orp, "zobell'soffset\\s*(.*?)\\s*mvtemperature")[,2]) & is.na(str_match(ph_orp, "quickcal\\s*(.*?)\\s*mvtemperature")[,2]),
  #                                             str_match(ph_orp, "orpstandardoffset\\s*(.*?)\\s*mvtemperature")[,2],
  #                                             ifelse(is.na(str_match(ph_orp, "zobell'soffset\\s*(.*?)\\s*mvtemperature")[,2]) & is.na(str_match(ph_orp, "orpstandardoffset\\s*(.*?)\\s*mvtemperature")[,2]),
  #                                                    str_match(ph_orp, "quickcal\\s*(.*?)\\s*mvtemperature")[,2],
  #                                                    str_match(ph_orp, "zobell'soffset\\s*(.*?)\\s*mvtemperature")[,2])),
  #
  #                         # Conductivity
  #                         cond_cal_date = as.character(mdy(str_match(conductivity, "lastcalibrated\\s*(.*?)\\s*calibrationdetails")[,2])),
  #                         tds_conversion_ppm = str_sub(str_match(conductivity, "tdsconversionfactor\\s*(.*?)\\s*cellconstant")[,2], 6, nchar(str_match(conductivity, "tdsconversionfactor\\s*(.*?)\\s*cellconstant")[,2])),
  #
  #                         # calibration report formatting has changed in 2024 for this variable. Therefore a post-2024 correction must occur
  #                         cond_cell_constant = ifelse(year(DT) < 2024, str_match(conductivity, "cellconstant\\s*(.*?)\\s*referencetemperature")[,2],
  #                                                     str_match(conductivity, "cellconstant\\s*(.*?)\\s*offset")[,2]),
  #
  #                         cond_offset = ifelse(year(DT) < 2024, NA,
  #                                              str_match(conductivity, "offset\\s*(.*?)\\s*µs/cm")[,2]),
  #
  #                         cond_pre = str_match(conductivity,paste0(str_match(conductivity,
  #                                                                            "premeasurementactual\\s*(.*?)\\s*specificconductivity")[,2],"specificconductivity\\s*(.*?)\\s*µs/cmpost"))[,2],
  #                         cond_post = str_match(conductivity,paste0(str_match(conductivity,
  #                                                                             "postmeasurementactual\\s*(.*?)\\s*specificconductivity")[,2],"specificconductivity\\s*(.*?)\\s*µs/cm"))[,2]) %>%
  #                         # if(turbidity == "No Turbidity Sensor"){
  #                         # # Turbidity
  #                         # turb_cal_date = "None",
  #                         # ntu_slope = "None",
  #                         # ntu_offset = "None",
  #                         # ntu_10 = "None",
  #                         # ntu_100 = "None") %>%
  #
  #     select(-c(ph_7_nice, ph_7_high, ph_7_low))
  #
  #     # Not all sondes have depth.
  #     if(!is.na(str_match(depth, "lastcalibrated"))){#\\s*(.*?)\\s*calibrationdetails")[,2])){
  #       cal_table <- cal_table %>%
  #         mutate(
  #           # Depth
  #           depth_cal_date = as.character(mdy(str_match(depth, "lastcalibrated\\s*(.*?)\\s*calibrationdetails")[,2])),
  #           depth_offset = str_match(depth, "zerooffset\\s*(.*?)\\s*psireferencedepth")[,2],
  #           depth_ref_depth = str_match(depth, "psireferencedepth\\s*(.*?)\\s*ftreferenceoffset")[,2],
  #           depth_ref_offset = str_match(depth, "ftreferenceoffset\\s*(.*?)\\s*psipremeasurement")[,2],
  #           depth_pre_psi = str_match(depth, "psipremeasurement\\s*(.*?)\\s*psipostmeasurement")[,2],
  #           depth_post_psi = str_match(depth, "psipostmeasurement\\s*(.*?)\\s*psi")[,2])
  #     }
  #
  #     if(depth == "No Depth Sensor"){
  #
  #       cal_table <- cal_table %>%
  #         mutate(# Depth
  #       depth_cal_date = "No Depth Sensor",
  #       depth_offset = "No Depth Sensor",
  #       depth_ref_depth = "No Depth Sensor",
  #       depth_ref_offset = "No Depth Sensor",
  #       depth_pre_psi = "No Depth Sensor",
  #       depth_post_psi = "No Depth Sensor")
  #     }
  #
  #
  #     if(!is.na(str_match(turbidity, "lastcalibrated"))){#calibrationpoint1premeasurement\\s*(.*?)\\s*ntupost")[,2])){
  #       # Not all sondes have turbidity.
  #       cal_table <- cal_table %>%
  #         mutate(
  #           # Turbidity
  #           turb_cal_date = as.character(mdy(str_match(turbidity, "lastcalibrated\\s*(.*?)\\s*calibrationdetails")[,2])),
  #           ntu_slope = str_match(turbidity, "slope\\s*(.*?)\\s*offset")[,2],
  #           ntu_offset = str_match(turbidity, "offset\\s*(.*?)\\s*ntu")[,2],
  #           ntu_10 = str_match(turbidity, "calibrationpoint1premeasurement\\s*(.*?)\\s*ntupost")[,2],
  #           ntu_100 = str_match(turbidity, "calibrationpoint2premeasurement\\s*(.*?)\\s*ntupost")[,2])
  #     }
  #
  #     if(turbidity == "No Turbidity Sensor"){
  #       cal_table <- cal_table %>%
  #         mutate(
  #           # Turbidity
  #           turb_cal_date = "No Turbidity Sensor",
  #           ntu_slope = "No Turbidity Sensor",
  #           ntu_offset = "No Turbidity Sensor",
  #           ntu_10 = "No Turbidity Sensor",
  #           ntu_100 = "No Turbidity Sensor")
  #
  #
  #
  #
  #     }
  #
  #
  #     cal_table <- cal_table %>%
  #       mutate(
  #         #Factory Defaults
  #         factory_defaults = paste0(ifelse(grepl("factorydefault", turbidity), "Turbidity ", ""),
  #                                   ifelse(grepl("factorydefault", rdo), "RDO ", ""),
  #                                   ifelse(is.na(ph_slope_post), "pH ", ""),
  #                                   ifelse(is.na(orp_offset), "ORP ", ""),
  #                                   ifelse(grepl("factorydefault", conductivity), "Conductivity ", ""),
  #                                   ifelse(grepl("factorydefaults", depth), "Depth ", ""))) %>%
  #       # convert all columns to character values to preserve info
  #       mutate(across(.cols = everything(), .fns = as.character)) %>%
  #       # remove "," from big numbers
  #       mutate(across(everything(), ~str_replace_all(., ",", "")))
  #
  #   }
  #
  #   bad_cal_interval_list <- map2(
  #     .x = bad_cal_records_filtered$start_DT,
  #     .y = bad_cal_records_filtered$end_DT,
  #     .f = ~interval(.x, .y, tz = "MST"))
  #
  #   if(df_parameter == "DO"){
  #
  #     cal_table <- list.files("data/calibration_reports", pattern=".html") %>%
  #       .[grepl(df_site, ., ignore.case = TRUE)] %>%
  #       map_dfr(., cal_tabler) %>%
  #       distinct(.keep_all = TRUE) %>%
  #       mutate(DT = as.character(round_date(ymd_hms(DT, tz = "MST"), "15 minutes"))) %>%
  #       # mutate(across(-matches("date|site|DT|factory"), as.numeric)) %>%
  #       dplyr::select(DT, site, rdo_slope, rdo_offset)
  #
  #     df_mod <- df %>%
  #       left_join(., cal_table, by = c("DT_join" = "DT", "site")) %>%
  #       fill(names(cal_table)[!grepl("\\b(site|DT)\\b", names(cal_table))], .direction = "down") %>%
  #       mutate(rdo_slope_pre = as.numeric(rdo_slope), rdo_offset_pre = as.numeric(rdo_offset)) %>%
  #       select(names(df), contains(c("pre"))) %>%
  #       left_join(., cal_table, by = c("DT_join" = "DT", "site")) %>%
  #       fill(names(cal_table)[!grepl("\\b(site|DT)\\b", names(cal_table))], .direction = "up") %>%
  #       mutate(rdo_slope_post = as.numeric(rdo_slope), rdo_offset_post = as.numeric(rdo_offset)) %>%
  #       select(names(df), contains(c("pre", "post"))) %>%
  #       #mutate(raw = (mean -rdo_offset_pre)/rdo_offset_pre)
  #       mutate(raw = case_when(DT_round %within% bad_cal_interval_list & is.na(rdo_slope_pre) ~ mean,
  #                              DT_round %within% bad_cal_interval_list & !is.na(rdo_slope_pre) ~ ((mean - rdo_offset_pre)/rdo_slope_pre),
  #                              .default = mean),
  #              cal_fix = case_when(DT_round %within% bad_cal_interval_list ~ (raw*rdo_slope_post) + rdo_offset_post,
  #                                  .default = mean)) %>%
  #       add_flag(mean != cal_fix, "calibration fix") %>%
  #       mutate(raw = mean,
  #              mean = cal_fixed)
  #
  #
  #   }
  #
  #   if(df_parameter == "pH"){
  #
  #     cal_table <- list.files("data/calibration_reports", pattern=".html") %>%
  #       .[grepl(df_site, ., ignore.case = TRUE)] %>%
  #       map_dfr(., cal_tabler) %>%
  #       distinct(.keep_all = TRUE) %>%
  #       mutate(DT = as.character(round_date(ymd_hms(DT, tz = "MST"), "15 minutes"))) %>%
  #       # mutate(across(-matches("date|site|DT|factory"), as.numeric)) %>%
  #       dplyr::select(DT, site, ph_slope_pre, ph_offset_pre, ph_slope_post, ph_offset_post, factory_defaults)
  #
  #     df_mod <- df %>%
  #       left_join(., cal_table, by = c("DT_join" = "DT", "site")) %>%
  #       fill(names(cal_table)[!grepl("\\b(site|DT)\\b", names(cal_table))], .direction = "down") %>%
  #       mutate(ph_slope_pre = as.numeric(ph_slope_pre), ph_offset_pre = as.numeric(ph_offset_pre),
  #              ph_slope_post = as.numeric(ph_slope_post), ph_offset_post = as.numeric(ph_offset_post)) %>%
  #       mutate(raw = case_when(DT_round %within% bad_cal_interval_list & is.na(ph_slope_pre) & grepl("pH", factory_defaults, ignore.case = TRUE) ~ mean,
  #                              DT_round %within% bad_cal_interval_list & !is.na(ph_slope_pre) & !grepl("pH", factory_defaults, ignore.case = TRUE) ~ ((mean - ph_offset_pre)/ph_slope_pre),
  #                              .default = mean),
  #              cal_fix = case_when(DT_round %within% bad_cal_interval_list ~ (raw*ph_slope_post) + ph_offset_post,
  #                                  .default = mean)) %>%
  #       add_flag(mean != cal_fix, "calibration fix")
  #
  #   }
  #
  #   if(df_parameter == "")
  #
  # }
# }

fix_depth_cal <- function(df){


  if(!"Depth" %in% df$parameter){

    nope <- df %>% mutate(relative_depth = NA)

   return(nope)

  }

  if("Depth" %in% df$parameter & "archery" %in% df$site){

     #df <- all_data_flagged[["archery-Depth"]]

    site_depth <- df %>%
      mutate(relative_depth = ifelse(DT_round >= as_datetime('2022-05-21 15:45:00', "MST") & DT_round <= as_datetime('2022-05-24 15:45:00', "MST"),
                                     mean +
                                       abs(
                                         dplyr::filter(df, as.character(DT_round) == ('2022-05-24 15:15:00'))$mean -
                                           dplyr::filter(df, as.character(DT_round) == ('2022-05-24 15:45:00'))$mean
                                       ),
                                     mean))
    #return(archery_depth)

  } else if("Depth" %in% df$parameter & "timberline" %in% df$site){

    # df <- all_data_flagged[["timberline-Depth"]]

    site_depth <- df %>%
      dplyr::mutate(relative_depth = ifelse(year == "2022" & DT_round <= lubridate::as_datetime('2022-04-07 17:00:00', "MST"),
                                            as.numeric(mean) +
                                              abs(
                                                dplyr::filter(df, as.character(DT_round) == ('2022-04-07 16:15:00'))$mean -
                                                  dplyr::filter(df, as.character(DT_round) == ('2022-04-07 17:15:00'))$mean
                                              ),
                                            as.numeric(mean)))
    # return(timberline_depth)

   } else if ("Depth" %in% df$parameter & "legacy" %in% df$site) {

    # df <- all_data_flagged[["legacy-Depth"]]

    site_depth <- df %>%
      mutate(relative_depth = ifelse(DT_round >= as_datetime('2022-04-06 06:00:00', "MST") & DT_round <= as_datetime('2022-04-12 09:15:00', "MST"),
                                     as.numeric(mean) +
                                       abs(
                                         dplyr::filter(df, as.character(DT_round) == ('2022-04-12 09:15:00'))$mean -
                                           dplyr::filter(df, as.character(DT_round) == ('2022-04-12 09:30:00'))$mean
                                       ),
                                     mean))

    site_depth <-  site_depth %>%
      mutate(relative_depth = ifelse(DT_round >= as_datetime('2022-07-08 17:00:00', "MST") & DT_round <= as_datetime('2022-07-12 09:00:00', "MST"),
                                     relative_depth +
                                       abs(
                                         dplyr::filter(site_depth, as.character(DT_round) == ('2022-07-08 14:15:00'))$relative_depth -
                                           dplyr::filter(site_depth, as.character(DT_round) == ('2022-07-08 17:00:00'))$relative_depth
                                       ),
                                     relative_depth))

    site_depth <-  site_depth %>%
      mutate(relative_depth = ifelse(DT_round >= as_datetime('2022-07-22 11:30:00', "MST") & DT_round <= as_datetime('2022-07-25 14:15:00', "MST"),
                                     relative_depth +
                                       abs(
                                         dplyr::filter(site_depth, as.character(DT_round) == ('2022-07-22 09:45:00'))$relative_depth -
                                           dplyr::filter(site_depth, as.character(DT_round) == ('2022-07-22 11:30:00'))$relative_depth
                                       ),
                                     relative_depth))

       #return(legacy_depth)

  } else if ("Depth" %in% df$parameter & "tamasag" %in% df$site) {

    # df <- all_data_flagged[["tamasag-Depth"]]

    site_depth <- df %>%
      mutate(relative_depth = ifelse(DT_round <= "2022-04-24 07:15:00" & year == "2022",
                                     as.numeric(mean) +
                                       abs(
                                         dplyr::filter(df, as.character(DT_round) == ('2022-04-24 07:15:00'))$mean -
                                           dplyr::filter(df, as.character(DT_round) == ('2022-04-24 07:30:00'))$mean
                                       ),
                                     mean))

    # return(tamasag_depth)

  } else if("Depth" %in% df$parameter) {

    site_depth <- df %>%
             mutate(relative_depth = mean)
}

  depth_flagged <- site_depth %>%
    add_flag(., relative_depth != mean, "post-calibration") %>%
    dplyr::mutate(mean = relative_depth)

  return(depth_flagged)

}

#' @title Apply all flags to a data frame
#'
#' @description
#' A function that applies all flags to a data frame. This function is used to
#' apply all flags to a data frame in one step.
#'
#' @param data A data frame with a `flag` column.
#'
#' @param require An upstream dependency (i.e., alternative to `{necessary object}`).
#'
#' @return A data frame with a `flag` column that has been updated with all flags.
#'
#' @examples
#' flag_all_data(data = all_data_flagged$`archery-Actual Conductivity`)
#' flag_all_data(data = all_data_flagged$`boxelder-Temperature`)
#'
#' @seealso [add_field_flag()]
#' @seealso [add_spec_flag()]
#' @seealso [add_seasonal_flag()]
#' @seealso [add_na_flag()]
#' @seealso [add_repeat_flag()]
#' @seealso [add_suspect_flag()]
#' @seealso [add_malfunction_flag()]

flag_all_data <- function(data, require = NULL) {
  flagged_data <- data %>%
    add_field_flag() %>%
    add_spec_flag() %>%
    add_seasonal_flag() %>%
    add_na_flag() %>%
    add_repeat_flag() %>%
    add_suspect_flag() %>%
    add_malfunction_flag(malfunction_records = mWater_malfunction_records)
    # mutate(mean_public = ifelse(is.na(flag), mean, NA))
  return(flagged_data)
}

generate_daily_plot <- function(plot_data_arg, df_list_arg, site_arg, parameter_arg) {

  # site_param_df <- df_list_arg[[paste0(site_arg, "-", parameter_arg)]]

  start_date <- min(plot_data_arg$DT_round)
  end_date <- max(plot_data_arg$DT_round)

  #default is lower network
  site_vector <- c("tamasag", # rist
                   "legacy",
                   "lincoln",
                   "timberline",
                   "prospect",
                   "boxelder", # elc
                   "archery",
                   "river bluffs")

  if(network == "virridy"){ # this will be the new default
    # establish order for all the non-tributary sites
    sites_order <-  c("joei","cbri","chd","pfal","sfm","pbd","tamasag",
                      "legacy","lincoln","timberline","prospect","boxelder",
                      "archery","riverbluffs")
    # establish the order for the tributary sites
    trib_sites_order <- c("boxcreek", "archery", NA, "springcreek", "prospect",
                          NA, "penn", "sfm", "lbea")
  }

  # determining the index for the site of interest.
  if (site_arg %in% sites_order) {

    plot_filter <- tibble(site = c("joei","cbri","chd","pfal","sfm","pbd",
                                   "tamasag","legacy", "lincoln","timberline",
                                   "timberline virridy","prospect",
                                   "prospect virridy","boxelder","archery",
                                   "archery virridy","riverbluffs"))

    site_index <- which(sites_order == site_arg)
    site_list <- as.vector(na.omit(sites_order[max(1, site_index - 1):min(length(sites_order), site_index + 1)]))

    plot_filter <- plot_filter %>%
      filter(grepl(paste(site_list, collapse = "|"), site, ignore.case = TRUE),
             site != site_arg) %>%
      pull(site)

  } else {

    plot_filter <- tibble(site = c("boxcreek", "archery", "archery virridy",
                                   "springcreek", "prospect", "prospect virridy",
                                   "penn", "sfm", "lbea"))

    site_index <- which(trib_sites_order == site_arg)
    site_list <- as.vector(na.omit(trib_sites_order[max(1, site_index - 1):min(length(trib_sites_order), site_index + 1)]))

    plot_filter <- plot_filter %>%
      filter(grepl(paste(site_list, collapse = "|"), site, ignore.case = TRUE),
             site != site_arg) %>%
      pull(site)

  }

  # Get the relevant sonde data
  # relevant_sondes <- map(plot_filter,
  #                        ~ {
  #                          sonde_name <- paste0(.x,"-",parameter_arg)
  #                          tryCatch({
  #                            sonde_df <- df_list_arg[[sonde_name]]  %>%
  #                              filter(DT_round %within% interval(start_date, end_date))},
  #                            error = function(err) {
  #                              cat("Sonde ", sonde_name," not found.\n")})
  #                        })

  # get the relevant sonde data source
  relevant_sonde_source <- map(plot_filter, ~ {
    sonde_name <- paste0(.x, "-", parameter_arg)
    # Determine which directory to pull data from
    tryCatch({
      retrieve_relevant_data_name(sonde_name, interval_arg = interval(start_date, end_date))
    }, error = function(err) {
      return("all_data")
    })
  })

  # Get the relevant data
  relevant_sondes <- map2(plot_filter,
                          relevant_sonde_source,
                          function(name, source) {
                            sonde_name <- paste0(name, "-", parameter_arg)
                            # try to pull in the data
                            tryCatch({
                              get(source)[[sonde_name]] %>%
                                filter(DT_round %within% interval(start_date, end_date))
                            }, error = function(err) {
                              return(NULL)
                            })
                          })

  # combine the lists
  sonde_info <- map2(relevant_sondes, relevant_sonde_source, list)

  # Remove any NULL results from the list
  sonde_info <- keep(sonde_info, ~!is.null(.x[[1]]))

  # append site_df to relevant sonde list, clean list, and bind dfs
  # to find plot info
  relevant_dfs <- map(sonde_info, ~.x[[1]])

  # append plot_data_arg to relevant sonde list, clean list, and bind dfs
  daily_plot_data <- append(relevant_dfs, list(plot_data_arg)) %>%
    keep(~ !is.null(.)) %>%
    keep(~ nrow(.)>0) %>%
    bind_rows()

  # use the daily flag data day as flag_day
  flag_day <- min(plot_data_arg$DT_round)

  plot <- ggplot(data = daily_plot_data) +
    geom_point(data = filter(daily_plot_data, (site == site_arg)),
               aes(x=DT_round, y=mean, color=flag)) +
    # geom_line(data = filter(daily_plot_data, (site != site_arg)),
    #           aes(x=DT_round, y=mean, color=site)) +
    map(sonde_info, function(sonde_data) {
      data <- sonde_data[[1]]
      data_source <- sonde_data[[2]]

      y_column <- if (data_source == "all_data") "mean" else "mean_verified"

      geom_line(data = data, aes(x = DT_round, y = .data[[y_column]], color = site))
    }) +
    ggtitle(paste0(str_to_title(site_arg), " ", str_to_title(parameter_arg), " (", format(flag_day, "%B %d, %Y"), ")")) +
    labs(x = "Time",
         y = "Mean")

  plot <- add_threshold_lines(plot = plot,
                              plot_data = plot_data_arg,
                              site_arg = site_arg,
                              parameter_arg = parameter_arg)

  plot <- plot +
    geom_vline(xintercept = seq(min(plot_data_arg$DT_round), max(plot_data_arg$DT_round), by = "hour"),
               color = "lightgrey", size = 0.5) +
    scale_x_datetime(breaks = seq(min(plot_data_arg$DT_round), max(plot_data_arg$DT_round), by = "hour"),
                     labels = format(seq(min(plot_data_arg$DT_round), max(plot_data_arg$DT_round), by = "hour"), "%H")) +
    theme_bw() +
    theme(legend.position = 'right',
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1)) + # make the background white
    guides(color = guide_legend(nrow = 10, byrow = TRUE))

  return(plot)

}

# Generate histograms of a site parameter data frame based on the mean column.
# @param df A dataframe with a `mean` column.
# @param df_index The index of the dataframe.
# @return A histogram of the mean column.
# @examples
# generate_general_histogram(df = all_data_flagged$`archery-Actual Conductivity`, df_index = "archery-Actual Conductivity")
# generate_general_histogram(df = all_data_flagged$`boxelder-Temperature`, df_index = "boxelder-Temperature")
generate_general_histogram <- function(df, df_index) {

  site_param <- toupper(sub("-", " ", df_index, fixed = TRUE))
  n <- as.character(sum(!is.na(df$mean)))
  title <- paste0("Histogram of ", site_param, " (n = ", n, ")")

  # there should be checks for data here
  minimum <- floor(min(df$mean, na.rm = TRUE))
  maximum <- ceiling(max(df$mean, na.rm = TRUE))

  histogram <- ggplot(data = df, aes(x = mean)) +
    geom_histogram(
      breaks = seq(minimum, maximum, by = 1)
    ) +
    labs(title = title)

  return(histogram)
}

generate_initial_weekly_plots <- function(all_df_list, pending_df_list, site_arg, parameter_arg, flag_arg = NULL) {

  site_param <- paste0(site_arg, "-", parameter_arg)

  site_flag_dates <- pending_df_list[[site_param]] %>%
    group_by(y_w) %>%
    filter(any(verification_status == "SKIP") | any(!is_verified)) %>%
    ungroup()

  if (!is.null(site_flag_dates)){
    # vector of sites in the order that they are in spatially ----
    # some sites have some funkiness going on (not all of the sites are present in the final plot)
    #default is lower network
    sites_order <- c("tamasag", # rist
                     "legacy",
                     "lincoln",
                     "timberline",
                     "prospect",
                     "boxelder", # elc
                     "archery",
                     "river bluffs")

    if(network == "virridy"){ # this will be the new default
      # establish order for all the non-tributary sites
      sites_order <-  c("joei","cbri","chd","pfal","sfm","pbd","tamasag",
                        "legacy","lincoln","timberline","prospect","boxelder",
                        "archery","riverbluffs")
      # establish order for all the tributary sites
      trib_sites_order <- c("boxcreek", "archery", NA, "springcreek", "prospect",
                            NA, "penn", "sfm", "lbea")
    }

    # determining the sites relevant to the site of interest.
    if (site_arg %in% sites_order) {

      plot_filter <- tibble(site = c("joei","cbri","chd","pfal","sfm","pbd",
                                     "tamasag","legacy", "lincoln","timberline",
                                     "timberline virridy","prospect",
                                     "prospect virridy","boxelder","archery",
                                     "archery virridy","riverbluffs"))

      site_index <- which(sites_order == site_arg)
      site_list <- as.vector(na.omit(sites_order[max(1, site_index - 1):min(length(sites_order), site_index + 1)]))

      plot_filter <- plot_filter %>%
        filter(grepl(paste(site_list, collapse = "|"), site, ignore.case = TRUE),
               site != site_arg) %>%
        pull(site)

    } else {

      plot_filter <- tibble(site = c("boxcreek", "archery", "archery virridy",
                                     "springcreek", "prospect", "prospect virridy",
                                     "penn", "sfm", "lbea"))

      site_index <- which(trib_sites_order == site_arg)
      site_list <- as.vector(na.omit(trib_sites_order[max(1, site_index - 1):min(length(trib_sites_order), site_index + 1)]))

      plot_filter <- plot_filter %>%
        filter(grepl(paste(site_list, collapse = "|"), site, ignore.case = TRUE),
               site != site_arg) %>%
        pull(site)

    }

    if (nrow(site_flag_dates >0)) {

      if (is.null(flag_arg)) {
        # This for loop generates an overlayed plot of weekly data for the site of
        # interest sandwiched by the site above and below it for each day that was
        # tagged by a flag of interest
        plot_list <- list()

        grouped_data <- site_flag_dates %>%
          group_by(y_w) %>% #group_by(week, year) %>% # group_by(week, month, year) %>%
          group_split()

        for(i in 1:length(grouped_data)) {

          group_data <- grouped_data[[i]]

          year_week <- unique(group_data$y_w)

          # filtering dfs of interest for the week of interest
          site_df <- site_flag_dates %>%
            filter(y_w == year_week)

          # Get the relevant sonde data
          relevant_sondes <- map(plot_filter, ~ {
            sonde_name <- paste0(.x, "-", parameter_arg)
            data_source <- NULL
            sonde_df <- NULL

            # Determine which directory to pull data from
            tryCatch({
              data_source <- retrieve_relevant_data_name(sonde_name, year_week)
              # cat("Data for",sonde_name,"will be pulled from",data_source,"\n")
            }, error = function(err) {
              # cat("Data for",sonde_name,"not found.\n")
              return(NULL)  # Return NULL if data source can't be determined
            })

            # Only try to pull in the data if data_source was successfully determined
            if (!is.null(data_source)) {
              tryCatch({
                sonde_df <- get(data_source)[[sonde_name]] %>%
                  filter(y_w == group_data$y_w)
              }, error = function(err) {
                cat("Sonde", sonde_name, "not found.\n")
                return(NULL)  # Return NULL if sonde data can't be retrieved
              })
            }

            # Only return a list if both data_source and sonde_df are available
            if (!is.null(data_source) & !is.null(sonde_df)) {
              return(list(sonde_df = sonde_df, data_source = data_source))
            } else {
              return(NULL)  # Return NULL if either data_source or sonde_df is NULL
            }
          })

          # Remove any NULL results from the list
          relevant_sondes <- compact(relevant_sondes)

          # append site_df to relevant sonde list, clean list, and bind dfs
          # to find plot info
          relevant_dfs <- map(relevant_sondes, ~.x[[1]])
          week_plot_data <- append(relevant_dfs, list(site_df)) %>% # how to relevant sondes here
            keep(~ !is.null(.)) %>%
            keep(~ nrow(.)>0) %>%
            bind_rows() %>%
            arrange(day)

          # Create a sequence of dates for the vertical lines
          start_date <- floor_date(min(week_plot_data$DT_round), "day")
          end_date <- ceiling_date(max(week_plot_data$DT_round), "day")
          vline_dates <- seq(start_date, end_date, by = "day")

          date_seq <- seq(min(vline_dates), max(vline_dates) - days(1), by = "1 day") + hours(12)

          # Use the first day of the group as flag_day
          flag_day <- min(group_data$DT_round)

          week_plot <- ggplot(data = week_plot_data) +
            geom_point(data = filter(week_plot_data, (site == site_arg)),
                       aes(x=DT_round, y=mean, color=flag)) +
            map(relevant_sondes, function(sonde_data) {
              data <- sonde_data[[1]]
              data_source <- sonde_data[[2]]

              y_column <- if (data_source == "all_data") "mean" else "mean_verified"

              geom_line(data = data, aes(x = DT_round, y = .data[[y_column]], color = site))
            })+
            geom_vline(xintercept = vline_dates, color = "black") +
            ggtitle(paste0(str_to_title(site_arg), " ", parameter_arg, " (", format(flag_day, "%B %d, %Y"), ")")) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            labs(x = "Day",
                 y = "Mean")

          week_plot <- add_threshold_lines(plot = week_plot,
                                           plot_data = week_plot_data,
                                           site_arg = site_arg,
                                           parameter_arg = parameter_arg)

          week_plot <- week_plot +
            theme_bw() +
            scale_x_datetime(date_breaks = "1 day",
                             date_labels = "%b %d",
                             minor_breaks = date_seq,
                             sec.axis = sec_axis(~., breaks = date_seq, labels = unique(week_plot_data$weekday))) +
            theme(legend.position = 'bottom',
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()) +
            guides(color = guide_legend(nrow = 4, byrow = TRUE))

          plot_list[[paste(site_param, as.character(flag_day))]] <- week_plot
          sorted_plot_names <- names(plot_list)[order(names(plot_list))]

          plot_list <- plot_list[sorted_plot_names]

        }
      }
      # ---- if flag arg is not null ----
      if (!is.null(flag_arg)) {
        print("under construction...")
        # # This for loop generates an overlayed plot of weekly data for the site of
        # # interest sandwiched by the site above and below it for each day that was
        # # tagged by a flag of interest
        # plot_list <- list()
        #
        # site_flag_weeks <- site_flag_dates %>%
        #   filter(str_detect(flag, flag_arg)) %>%
        #   group_by(week, year) %>%
        #   slice(1)
        #
        # grouped_data <- site_flag_dates %>%
        #   filter(y_w %in% site_flag_weeks$y_w) %>%
        #   group_by(week, year) %>% # group_by(week, month, year) %>%
        #   group_split()
        #
        #
        # for(i in 1:length(grouped_data)) {
        #
        #   group_data <- grouped_data[[i]]
        #
        #   # flag_title <- site_flag_dates$flag[i] # no flag title ***
        #
        #   # filtering dfs of interest for the weeks where a flag was detected
        #   site_df <- site_flag_dates %>%
        #     filter(y_w == group_data$y_w)
        #
        #   # TryCatch used here to avoid erroring out on the first and last values of
        #   # sites_order object (there is no prior/next record after the first/last record).
        #   # Return df as NULL in case of an error
        #   prev_site_df <- NULL
        #   next_site_df <- NULL
        #
        #   tryCatch({
        #     previous_site <- paste0(sites_order[site_index-1],"-",parameter_arg)
        #     prev_site_df <- all_df_list[[previous_site]] %>%
        #       filter(y_w == group_data$y_w)},
        #     error = function(err) {
        #       cat("No previous site.\n")})
        #
        #   tryCatch({
        #     next_site <- paste0(sites_order[site_index+1],"-",parameter_arg)
        #     next_site_df <- all_df_list[[next_site]] %>%
        #       filter(y_w == group_data$y_w)},
        #     error = function(err) {
        #       cat("No next site.\n")})
        #
        #   # Bind all three dfs
        #   week_plot_data <- list(site_df, prev_site_df, next_site_df) %>%
        #     # remove NULL values from the list
        #     keep(~ !is.null(.)) %>%
        #     bind_rows()
        #
        #   # Create a sequence of dates for the vertical lines
        #   start_date <- floor_date(min(week_plot_data$DT_round), "day")
        #   end_date <- ceiling_date(max(week_plot_data$DT_round), "day")
        #   vline_dates <- seq(start_date, end_date, by = "day")
        #
        #   date_seq <- seq(min(vline_dates), max(vline_dates) - days(1), by = "1 day") + hours(12) ## here ----
        #
        #   # Use the first day of the group as flag_day
        #   flag_day <- min(group_data$DT_round)
        #
        #   week_plot <- ggplot(data = week_plot_data) +
        #     geom_point(data = filter(week_plot_data, (site == site_arg)),
        #                aes(x=DT_round, y=mean, color=flag)) +
        #     geom_line(data = filter(week_plot_data, (site != site_arg)),
        #               aes(x=DT_round, y=mean, color=site)) +
        #     geom_vline(xintercept = vline_dates, color = "black") +
        #     ggtitle(paste0(str_to_title(site_arg), " ", parameter_arg, " (", format(flag_day, "%B %d, %Y"), ")")) +
        #     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        #     labs(x = "Day",
        #          y = "Mean")
        #
        #   week_plot <- add_threshold_lines(plot = week_plot,
        #                                    plot_data = week_plot_data,
        #                                    site_arg = site_arg,
        #                                    parameter_arg = parameter_arg)
        #
        #   week_plot <- week_plot +
        #     theme_bw() +
        #     scale_x_datetime(date_breaks = "1 day",
        #                      date_labels = "%b %d",
        #                      minor_breaks = date_seq) +
        #     theme(legend.position = 'bottom',
        #           panel.grid.major = element_blank(),
        #           panel.grid.minor = element_blank()) +
        #     annotate("text", x = date_seq, y = min(week_plot_data$mean, na.rm = TRUE) - 1, label = 1:length(date_seq), hjust = 0) +
        #     guides(color = guide_legend(nrow = 4, byrow = TRUE))
        #
        #   plot_list[[paste(site_param, as.character(flag_day))]] <- week_plot
        #   sorted_plot_names <- names(plot_list)[order(names(plot_list))]
        #
        #   plot_list <- plot_list[sorted_plot_names]
        # }
      }
      return(plot_list)
    } else {
      return(paste(flag_arg, "not detected.\n"))
    }
  } else {
    return(paste(site_arg, parameter_arg, "combination not available.\n"))
  }

}

# Generate seasonal histograms for a given site parameter data frame based on the mean column.
# @param df A dataframe with a `mean` column.
# @param df_index The index of the dataframe.
# @return A plot with a histogram of the mean column for each season.
# @examples
# generate_seasonal_histogram(df = all_data_flagged$`archery-Actual Conductivity`, df_index = "archery-Actual Conductivity")
# generate_seasonal_histogram(df = all_data_flagged$`boxelder-Temperature`, df_index = "boxelder-Temperature")

generate_seasonal_histogram <- function(df, df_index) {

  winter_baseflow <- c(12,1,2,3,4)
  snowmelt <- c(5,6,NA,NA,NA)
  monsoon <- c(7,8,9,NA,NA)
  fall_baseflow <- c(10,11,NA,NA,NA)

  # do water years

  seasons <- data.frame(winter_baseflow, snowmelt, monsoon, fall_baseflow)

  site_param <- toupper(sub("-", " ", df_index, fixed = TRUE))

  param <- unique(na.omit(df$parameter))

  hist_list <- list()
  for (i in colnames(seasons)){

    filtered_df <- df %>%
      filter(month %in% seasons[[i]],
             !str_detect(flag, "sensor specification range"))

    n <- as.character(sum(!is.na(filtered_df$mean)))

    title <- paste0(i," (n = ",n ,")")

    histogram <- ggplot() +
      geom_histogram() +
      labs(title = title)

    tryCatch({
      minimum <- floor(min(filtered_df$mean, na.rm = TRUE))
      maximum <- ceiling(max(filtered_df$mean, na.rm = TRUE))

      if (param %in% c("Specific Conductivity", "Actual Conductivity", "Turbidity")) {
        bins <- seq(minimum, maximum, by = 10)
      } else if (param %in% c("ORP", "pH", "Depth")) {
        bins <- seq(minimum, maximum, by = 0.05)
      } else {
        bins <- seq(minimum, maximum)
      }

      x_min <- filtered_df$m_mean05[1]
      x_max <- filtered_df$m_mean99[1]

      histogram <- ggplot(data = filtered_df, aes(x = mean)) +
        geom_histogram(breaks = bins) +
        geom_vline(xintercept = x_min, color = "red", linetype = "dashed") +
        geom_vline(xintercept = x_max, color = "red", linetype = "dashed") +
        facet_wrap(~ year, nrow = 1) +
        labs(title = title)
    },
    error = function(err) {
      cat("No finite values for", site_param, i, "\n")})

    hist_list[[i]] <- histogram

  }

  collated_hist <- ggarrange(plotlist = hist_list, nrow=2, ncol=2) %>%
    annotate_figure(top = site_param)

  return(collated_hist)

}

#' @title Generate site metaparameter data from the HydroVu API.
#'
#' @description
#' A function that generates a site metaparameter dataframe from the HydroVu API.
#' A metaparameter is a parameter that is used to generate flags for parameters
#' other than itself.
#'
#' @param api_data A dataframe with the munged API data.
#'
#' @param require An upstream dependency (i.e., alternative to `{necessary object}`).
#'
#' @return A dataframe with site metaparameter data that will be joined to the other
#' site-parameter data.
#'
#' @examples
#' generate_site_metaparam(site_arg = "archery", metaparameter_arg = "Battery Level", api_data = incoming_data_collated_csvs)
#'
#' @seealso [summarize_site_param()]

generate_site_metaparam <- function(api_data, require = NULL) {

    sites <- unique(api_data$site)
    metaparameters <- c("Temperature", "Battery Level", "Baro", "External Voltage")

    df_list <- list()
    for (i in sites) {
        metaparameter_data <- api_data %>%
            data.table::data.table() %>%
            dplyr::select(DT_join, site, parameter, value) %>%
            dplyr::filter(site == i & (parameter %in% metaparameters)) %>%
            dplyr::select(-site) %>%
            tidyr::pivot_wider(names_from = parameter, values_from = value)
        df_list[[i]] <- metaparameter_data
    }
    return(df_list)
}

#' @title Generate Summary Statistics
#'
#' @description
#' A function that generates summary statistics for a given site parameter data frame.
#' The generated statistics include:
#'   - The next value and previous value for the mean.
#'   - The rolling 7-point median of the mean.
#'   - The rolling 7-point mean of the mean.
#'   - The rolling 7-point standard deviation of the mean.
#'   - The slope of a point in relation to the point ahead and behind.
#'   - The rolling 7-point slope of the mean.
#'   - The month and year of each data point.
#'   - The year-month combination.
#'
#' @param site_param_df A data frame with a `mean` column retrieved from HydroVu API.
#'
#' @return A data frame with summary statistics for a given site parameter data frame.
#'
#' @examples
#' generate_summary_statistics(site_param_df = all_data_flagged$`archery-Actual Conductivity`)
#' generate_summary_statistics(site_param_df = all_data_flagged$`boxelder-Temperature`)

generate_summary_statistics_full <- function(site_param_df) {

  summary_stats_df <- site_param_df %>%
    # ... so that we can get the proper leading/lagging values across our entire timeseries:
    dplyr::mutate(
      # Add the next value and previous value for mean.
      front1 = dplyr::lead(mean, n = 1),
      back1 = dplyr::lag(mean, n = 1),
      # Add the rolling 7-point median (using itself + data of the past).
      rollmed = RcppRoll::roll_median(mean, n = 7, align = 'right', na.rm = F, fill = NA_real_), #ifelse(is.na(rollmed), roll_median(mean, n = 7, align = 'right', na.rm = F, fill = NA_real_), rollmed), # to go (j): check_na() function for when we append data
      # Add the rolling 7-point mean (using itself + data of the past).
      rollavg = RcppRoll::roll_mean(mean, n = 7, align = 'right', na.rm = F, fill = NA_real_), #ifelse(is.na(rollavg), roll_mean(mean, n = 7, align = 'right', na.rm = F, fill = NA_real_), rollavg),
      # Add the rolling 7-point standard deviation (using itself + data of the past).
      rollsd = RcppRoll::roll_sd(mean, n = 7, align = 'right', na.rm = F, fill = NA_real_), #ifelse(is.na(rollsd), roll_sd(mean, n = 7, align = 'right', na.rm = F, fill = NA_real_), rollsd),
      # Determine the slope of a point in relation to the point ahead and behind
      slope_ahead = (front1 - mean)/15,
      slope_behind = (mean - back1)/15,
      # Get the rolling 7-point slope (using itself + data of the past).
      rollslope = RcppRoll::roll_mean(slope_behind, n = 7, align = 'right', na.rm = F, fill = NA_real_), #ifelse(is.na(rollslope), roll_mean(slope_behind, n = 7, align = 'right', na.rm = F, fill = NA_real_), rollslope),
      # add some summary info for future steps
      month = lubridate::month(DT_round),
      year = lubridate::year(DT_round),
      y_m = paste0(year, '-', month),
      # Define our seasons:
      season = dplyr::case_when(month %in% c(12, 1, 2, 3, 4) ~ "winter_baseflow",
                                month %in% c(5, 6) ~ "snowmelt",
                                month %in% c(7, 8, 9) ~ "monsoon",
                                month %in% c(10, 11) ~ "fall_baseflow",
                                TRUE ~ NA)
    )

  return(summary_stats_df)

}

# Generate summary statistics for a given site parameter data frame.
#' @param site_param_df A data frame with a `mean` column retrieved from HydroVu API.
#' @return A data frame with summary statistics for a given site parameter data frame.
#' @examples
# generate_summary_statistics(site_param_df = all_data_flagged$`archery-Actual Conductivity`)
# generate_summary_statistics(site_param_df = all_data_flagged$`boxelder-Temperature`)

generate_summary_statistics <- function(site_param_df) {

  summary_stats_df <- site_param_df %>%
    # ... so that we can get the proper leading/lagging values across our entire timeseries:
    dplyr::mutate(
      # Add the next value and previous value for mean.
      # Only do this for newest data (i.e., our appended historical
      # data already has these filled out and we don't want to over-
      # write them)
      front1 = ifelse(is.na(front1), dplyr::lead(mean, n = 1), front1),
      back1 = ifelse(is.na(back1), dplyr::lag(mean, n = 1), back1),
      # Add the median for a point and 6 points behind it:
      rollmed = ifelse(is.na(rollmed), RcppRoll::roll_median(mean, n = 7, align = 'right', na.rm = F, fill = NA_real_), rollmed),
      # Add the mean for a point and 6 points behind it:
      rollavg = ifelse(is.na(rollavg), RcppRoll::roll_mean(mean, n = 7, align = 'right', na.rm = F, fill = NA_real_), rollavg),
      # Add the standard deviation for a point and 6 points behind it:
      rollsd = ifelse(is.na(rollsd), RcppRoll::roll_sd(mean, n = 7, align = 'right', na.rm = F, fill = NA_real_), rollsd),
      # Determine the slope of a point in relation to the point ahead and behind.
      slope_ahead = ifelse(is.na(slope_ahead), (front1 - mean)/15, slope_ahead),
      slope_behind = ifelse(is.na(slope_behind), (mean - back1)/15, slope_behind),
      # Get the rolling 7-point slope (using itself + data of the past).
      rollslope = ifelse(is.na(rollslope), RcppRoll::roll_mean(slope_behind, n = 7, align = 'right', na.rm = F, fill = NA_real_), rollslope),
      # add some summary info for future us
      month = lubridate::month(DT_round),
      year = lubridate::year(DT_round),
      y_m = paste(year, '-', month),
      season = dplyr::case_when(month %in% c(12, 1, 2, 3, 4) ~ "winter_baseflow",
                                month %in% c(5, 6) ~ "snowmelt",
                                month %in% c(7, 8, 9) ~ "monsoon",
                                month %in% c(10, 11) ~ "fall_baseflow",
                                TRUE ~ NA)
    )

  return(summary_stats_df)

}

generate_supplemental_weekly_plot <- function(daily_plot_data_arg, df_list_arg, site_arg, parameter_arg) {

  site_param_df <- df_list_arg[[paste0(site_arg, "-", parameter_arg)]]

  start_date <- min(daily_plot_data_arg$DT_round) - days(3) # TODO: replace all daily_plot_data_arg instances with something more efficient
  end_date <- max(daily_plot_data_arg$DT_round) + days(3)

  week_plot_data <- site_param_df %>%
    filter(DT_round %within% interval(start_date, end_date))

  #default is lower network
  site_vector <- c("tamasag", # rist
                   "legacy",
                   "lincoln",
                   "timberline",
                   "prospect",
                   "boxelder", # elc
                   "archery",
                   "river bluffs")

  if(network == "virridy"){ # this will be the new default
    # establish order for all the non-tributary sites
    sites_order <-  c("joei","cbri","chd","pfal","sfm","pbd","tamasag",
                      "legacy","lincoln","timberline","prospect","boxelder",
                      "archery","riverbluffs")
    # establish the order for the tributary sites
    trib_sites_order <- c("boxcreek", "archery", NA, "springcreek", "prospect",
                          NA, "penn", "sfm", "lbea")
  }

  # determining the index for the site of interest.
  if (site_arg %in% sites_order) {

    plot_filter <- tibble(site = c("joei","cbri","chd","pfal","sfm","pbd",
                                   "tamasag","legacy", "lincoln","timberline",
                                   "timberline virridy","prospect",
                                   "prospect virridy","boxelder","archery",
                                   "archery virridy","riverbluffs"))

    site_index <- which(sites_order == site_arg)
    site_list <- as.vector(na.omit(sites_order[max(1, site_index - 1):min(length(sites_order), site_index + 1)]))

    plot_filter <- plot_filter %>%
      filter(grepl(paste(site_list, collapse = "|"), site, ignore.case = TRUE),
             site != site_arg) %>%
      pull(site)

  } else {

    plot_filter <- tibble(site = c("boxcreek", "archery", "archery virridy",
                                   "springcreek", "prospect", "prospect virridy",
                                   "penn", "sfm", "lbea"))

    site_index <- which(trib_sites_order == site_arg)
    site_list <- as.vector(na.omit(trib_sites_order[max(1, site_index - 1):min(length(trib_sites_order), site_index + 1)]))

    plot_filter <- plot_filter %>%
      filter(grepl(paste(site_list, collapse = "|"), site, ignore.case = TRUE),
             site != site_arg) %>%
      pull(site)

  }

  # get the relevant sonde data source
  relevant_sonde_source <- map(plot_filter, ~ {
    sonde_name <- paste0(.x, "-", parameter_arg)
    # Determine which directory to pull data from
    tryCatch({
      retrieve_relevant_data_name(sonde_name, interval_arg = interval(start_date, end_date))
    }, error = function(err) {
      return("all_data")
    })
  })

  # Get the relevant data
  relevant_sondes <- map2(plot_filter,
                          relevant_sonde_source,
    function(name, source) {
    sonde_name <- paste0(name, "-", parameter_arg)
    # try to pull in the data
    tryCatch({
      get(source)[[sonde_name]] %>%
        filter(DT_round %within% interval(start_date, end_date))
    }, error = function(err) {
      return(NULL)
    })
  })

  # combine the lists
  sonde_info <- map2(relevant_sondes, relevant_sonde_source, list)

  # Remove any NULL results from the list
  sonde_info <- keep(sonde_info, ~!is.null(.x[[1]]))

  # append site_df to relevant sonde list, clean list, and bind dfs
  # to find plot info
  relevant_dfs <- map(sonde_info, ~.x[[1]])

  week_plot_data <- append(relevant_dfs, list(week_plot_data)) %>%
    keep(~ !is.null(.)) %>%
    keep(~ nrow(.)>0) %>%
    bind_rows() %>%
    arrange(day)

  # Create a sequence of dates for the vertical lines
  start_date <- floor_date(min(week_plot_data$DT_round), "day")
  end_date <- ceiling_date(max(week_plot_data$DT_round), "day")
  vline_dates <- seq(start_date, end_date, by = "day")

  date_seq <- seq(min(vline_dates), max(vline_dates) - days(1), by = "1 day") + hours(12)

  # use the daily flag data day as flag_day
  flag_day <- unique(daily_plot_data_arg$DT_round)

  week_plot <- ggplot(data = week_plot_data) +
    geom_point(data = filter(week_plot_data, (site == unique(daily_plot_data_arg$site))),
               aes(x=DT_round, y=mean, color=flag)) +
    map(sonde_info, function(sonde_data) {

      data <- sonde_data[[1]]
      data_source <- sonde_data[[2]]

      y_column <- if (data_source == "all_data") "mean" else "mean_verified"

      geom_line(data = data, aes(x = DT_round, y = .data[[y_column]], color = site))
    }) +
    geom_rect(data = daily_plot_data_arg, aes(xmin = min(DT_round), xmax = max(DT_round),
                                              ymin = -Inf, ymax = Inf),
              fill = "grey",
              alpha = 0.01,
              color = NA) +
    geom_vline(xintercept = vline_dates, color = "black") +
    ggtitle(paste0(str_to_title(unique(daily_plot_data_arg$site)), " ", unique(daily_plot_data_arg$parameter), " (", format(flag_day, "%B %d, %Y"), ")")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Day",
         y = "Mean")

  week_plot <- add_threshold_lines(plot = week_plot,
                                   plot_data = week_plot_data,
                                   site_arg = unique(daily_plot_data_arg$site),
                                   parameter_arg = unique(daily_plot_data_arg$parameter))

  week_plot <- week_plot +
    theme_bw() +
    scale_x_datetime(date_breaks = "1 day",
                     date_labels = "%b %d",
                     minor_breaks = date_seq) +
    theme(legend.position = 'right',
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) + # make the background white
    guides(color = guide_legend(nrow = 10, byrow = TRUE))

  return(week_plot)
}

# Get day decision.
# @param prompt_text A string for the prompt text.
# @return A string of "pass", "fail", or "inspect" depending on input from the user.
# @examples
# get_day_decision(prompt_text = "Would you like to (pass/fail/inspect) all data points for: ")

get_day_decision <- function(prompt_text) {
  while (TRUE) {
    user_input <- readline(prompt = paste(prompt_text))
    user_input <- tolower(user_input)

    if (user_input %in% c("pass", "p")) {
      return("pass")
    } else if (user_input %in% c("fail", "f")) {
      return("fail")
    } else if (user_input %in% c("inspect", "i")){
      return("inspect")
    } else {
      cat("Invalid input. Please enter 'yes', 'no', or 'inspect'.\n")
    }
  }
}

get_dt_inspection_decisions <- function(daily_plot_data) {
  # Create a sequence of date times from the daily plot data
  dt_list <- sort(unique(pull(daily_plot_data, DT_round)))

  # Get the date
  plot_date <- as_date(daily_plot_data$DT_round[1])

  # Prompt for the number of intervals
  prompt_text <- paste("How many intervals would you like to inspect?\n(Must be less than", length(dt_list), ")\n")

  while (TRUE) {
    user_input <- readline(prompt = prompt_text)

    # The user input needs to check for:
    # Only numbers between 1 and length(dt_list)
    num_intervals <- suppressWarnings(as.numeric(user_input))

    if (!is.na(num_intervals) && num_intervals > 0 && num_intervals < length(dt_list)) {
      break
    }

    cat("Invalid input. Please enter a valid number of intervals.\n")
  }

  selected_intervals <- list()

  for (i in 1:num_intervals) {
    prompt_text <- paste("Enter the time range for interval", i, ".\n(format 'HH:MM:SS-HH:MM:SS'):\n")

    while (TRUE) {
      user_input <- readline(prompt = prompt_text) # these need to be 15 minute intervals ***
      interval <- unlist(strsplit(user_input, "-"))

      interval_set <- paste(plot_date, interval)
      interval_set <- as_datetime(interval_set, tz = "MST")

      if (length(interval_set) == 2 && all(interval_set %in% dt_list)) {
        start_time <- interval_set[1]
        end_time <- interval_set[2]

        # Check if selected_intervals is empty
        if (start_time <= end_time & length(selected_intervals) == 0) {
          selected_intervals[[i]] <- interval(start_time, end_time)
          break
        }

        if (start_time <= end_time & !suppressWarnings(any(map(selected_intervals, ~int_overlaps(interval(start_time, end_time), .x))))) {
          selected_intervals[[i]] <- interval(start_time, end_time)
          break
        }
      }

      cat("Invalid input.\nPlease enter a valid time range that doesn't overlap with previous intervals.\n")
    }
  }

  return(selected_intervals)
}

#' #' Import Fort Collins Floodwarning Rain Gage Data
#' #'
#' #' This function pulls rain data for sites within the FC Floodwarning System
#' #'
#' #'
#' #' @param start Start date of when you want data.
#' #' @param end End date of when you want data; default is the current date.
#' #' @param save Whether to save (TRUE) the resulting table or not (FALSE)
#' #' @param path If `save = TRUE`, the file path to save the shapefile
#' #'
#' #' @return A table of time series flow data across the FC Floodwarning System
#' #'
#' get_fc_rain <- function(start = '2018-10-01', end = Sys.Date(), save = TRUE, path = 'data/context_data/'){
#'
#'   call <- "https://opendata.fcgov.com/api/views/g87z-rviz/rows.csv?accessType=DOWNLOAD"
#'
#'   #download dataset from FC
#'   temp1 <- tempfile()
#'   download.file(paste0(call), destfile = temp1, method = "curl")
#'   data <- read_csv(temp1) %>%
#'     dplyr::mutate(date = as_date(mdy_hms(Timestamp))) %>%
#'     dplyr::filter(date >= start & date <= end)
#'
#'   s_platte_ws s_platte_flowlines <- nhdplusTools::get_nhdplus(AOI = s_platte_ws,
#'                                                               realization = "flowline") <- nhdplusTools::get_nldi_basin(list(featureSource = "nwissite", featureID = "USGS-06754000"))
#'   s_platte_flowlines <- nhdplusTools::get_nhdplus(AOI = s_platte_ws,
#'                                                   realization = "flowline")
#'   s_platte_catchments <- nhdplusTools::get_nhdplus(AOI = s_platte_ws, realization = "catchment")
#'
#'   poudre <- s_platte_flowlines %>%
#'     filter(grepl("Cache la Poudre", gnis_name, ignore.case = TRUE))
#'
#'   poudre <- s_platte_catchments %>%
#'     filter(featureid %in% poudre$comid) %>%
#'     .[sites,] %>%
#'     rowid_to_column()
#'
#'   sites <- read_csv("data/metadata/sonde_location_metadata.csv") %>%
#'     separate(lat_long, into = c("lat", "long"), sep = ",") %>%
#'     sf::st_as_sf(coords = c("long", "lat"), crs = 4326) %>%
#'     filter(Site %in% c("Tamasag","Legacy","Lincoln","Timberline","Prospect","Boxelder","Archery","River Bluffs")) %>%
#'     mutate(rowid = sf::st_nearest_feature(., poudre)) %>%
#'     left_join(sf::st_drop_geometry(poudre), by = "rowid")
#'
#'   ws_maker <- function(site_names){
#'
#'     df <- sites %>%
#'       dplyr::filter(Site == site_names)
#'
#'     # get_UT() creates a vector of all flowlines upstream of the comid of interest ...
#'     UT_comids <- nhdplusTools::get_UT(network = s_platte_flowlines,
#'                                       comid = df$featureid)
#'
#'     catchments <- filter(s_platte_catchments, featureid %in% c(UT_comids)) %>%
#'       summarize() %>%
#'       mutate(site = df$Site)
#'
#'     return(catchments)
#'
#'   }
#'
#'   watersheds <- sites$Site %>%
#'     map_dfr(~ws_maker(.))
#'
#'   rain_sites <- fread("data/context_data/fc_rain.csv") %>%
#'     distinct(`Sensor Name`, Latitude, Longitude) %>%
#'     sf::st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>%
#'     .[watersheds,] %>%
#'     sf::st_join(watersheds) %>%
#'     mutate(site = tolower(site)) %>%
#'     write_csv("data/context_data/site_link.csv")
#'
#'   if(save == TRUE){
#'     write_csv(data, paste0(path, "/fc_rain.csv"))
#'   }
#'
#'   return(data)
#' }
#'
#'
#'
#'

# Get flag decision.
# @param prompt_text A string for the prompt text.
# @return A boolean of TRUE or FALSE depending on input from the user.
# @examples
# get_flag_decision(prompt_text = "Would you like to (pass/fail) this data point: ")

get_flag_decision <- function(prompt_text) {
  while (TRUE) {
    user_input <- readline(prompt = paste(prompt_text))
    user_input <- tolower(user_input)

    if (user_input %in% c("pass", "p")) {
      return(TRUE)
    } else if (user_input %in% c("fail", "f")) {
      return(FALSE)
    } else {
      cat("Invalid input. Please enter 'pass' or 'fail'.\n")
    }
  }
}

#' @title Get the start dates for the HydroVu API pull from the historically
#' flagged data
#'
#' @description
#' This function finds the max datetime in the temperature data
#' frames for each site's historically flagged data. We use temperature because
#' it is always being tracked in the sensors and thus is the most accurate
#' representation of the last API pull for a given sensor.
#'
#' @param incoming_historically_flagged_data_list This is the output from
#' `flagged_data_dfs` and is the historically flagged data.
#'
#' @returns
#' A data frame that has the sites and their corresponding starting
#' date-times' for the HydroVu API pull. This data frame will then get passed into
#' `incoming_data_csvs_upload`.

get_start_dates_df <- function(incoming_historically_flagged_data_list) {

  temperature_subset <- grep("Temperature", names(incoming_historically_flagged_data_list))

  start_dates_df <- map_dfr(incoming_historically_flagged_data_list[temperature_subset],
                            function(temperature_df){
                              temperature_df %>%
                                select(DT_round, site) %>%
                                filter(DT_round == max(DT_round))
                              })

  return(start_dates_df)

}

# @ param site: string, name of site to get filenames for, select from "tamasag" "legacy", "timberline" "prospect" "boxelder" "archery" "riverbluffs"
get_tl_photo_filenames <- function(site = "legacy",start_dt = "2023-05-20 12:00" ,end_dt = "2023-07-20 12:00"){
  #convert dates to DT objects
  start_dt <- ymd_hm(start_dt, tz = "MST")
  end_dt <- ymd_hm(end_dt, tz = "MST")
  if(site == "river bluffs"){
    site <- "riverbluffs"
  }

  if(site %nin% c("tamasag", "legacy", "timberline", "prospect", "boxelder", "archery", "riverbluffs")){
    print("Invalid site name, please select from: 'tamasag' 'legacy' 'timberline' 'prospect' 'boxelder' 'archery' 'riverbluffs'")

  nah <- tibble()
    return(nah)
  }

  #look through all the files in the site's folder
  all_files <- tibble(filename = list.files(paste0("data/timelapse_photos/2023_compiled/",site), full.names = TRUE, recursive = TRUE))%>%
    mutate(DT = str_extract(filename, "[0-9]{8}_[0-9]{4}"),
           DT = parse_date_time(DT, "Ymd_HM", tz = "MST"),
           #round DT to 30 min to match with sensor data
           DT_round = round_date(DT, "30 minutes"))%>%
    # if there are multiple photos with the same DT_round value, only keep one row with that DT_round
    distinct(DT_round, .keep_all = TRUE)%>%
    select(-DT)%>%
    filter(between(DT_round, start_dt, end_dt))

  if(nrow(all_files) == 0){
    print("No files found for this site and date range")
    nah <- tibble()
    return(nah)
  }else{
   return(all_files)
  }


}

#test <- get_tl_photo_filenames(site = "legacy", start_dt = "2023-7-20 12:00", end_dt = "2023-12-06 12:00")





#' @title Get verification decision
#'
#' @description
#' Get decisions from user input
#'
#' @param prompt_text The prompt that the user is replying to
#'
#' @examples
#' # get_day_decision(prompt_text = "Would you like to (pass/fail/inspect) all data points for: ")

get_verification_decision <- function(prompt_text) { # should add a layer level to this to prevent user from inspect sub-daily data

  while (TRUE) {
    user_input <- readline(prompt = paste(prompt_text))
    user_input <- tolower(user_input)

    # pass statements
    if (user_input %in% c("pass all", "pass", "pa", "p")) {
      return("PASS ALL")
    }
    if (user_input %in% c("pass valid", "pv")) {
      return("PASS VALID")
    }
    if (user_input %in% c("pass flagged", "pf")) {
      return("PASS FLAGGED")
    }
    if (user_input %in% c("pass none", "pn")) {
      return("PASS NONE")
    }

    # fail statements
    if (user_input %in% c("fail all", "fail", "fa", "f")) {
      return("FAIL ALL")
    }
    if (user_input %in% c("fail valid", "fv")) {
      return("FAIL VALID")
    }
    if (user_input %in% c("fail flagged", "ff")) {
      return("FAIL FLAGGED")
    }
    if (user_input %in% c("fail none", "fn")) {
      return("FAIL NONE")
    }

    # skip statements
    if (user_input %in% c("skip", "skip all", "s", "sa")) {
      return("SKIP")
    }

    # inspect statements ***
    if (user_input %in% c("inspect all", "inspect", "ia", "i")){
      return("INSPECT ALL")
    }

    if (user_input %in% c("inspect valid", "iv")){
      return("INSPECT VALID")
    }

    if (user_input %in% c("inspect flagged", "if")){
      return("INSPECT FLAGGED")
    }

    if(user_input %in% c("inspect some", "is")){ # if we are at daily level we can't go further
      return("INSPECT SOME")
    }

    # quit statement
    if (user_input %in% c("quit", "q")) {
      return("QUIT")
    }

    cat("Invalid input. Please enter one of the options from the following decision matrix:\n\n")
    cat("          | PASS | FAIL | SKIP | INSPECT  |  | QUIT |\n")
    cat("----------+------+------+------+----------|  |------|\n")
    cat(" ALL      | pa   | fa   | sa   | ia       |  | q    |\n")
    cat(" VALID    | pv   | fv   |      | iv       |  +------+\n")
    cat(" FLAGGED  | pf   | ff   |      | if       |\n")
    cat(" NONE     | pn   | fn   |      |          |\n")
    cat(" SOME     |      |      |      | is       |\n")
    cat("----------+------+------+------+----------|\n")
  }

}


get_watersheds <- function(sf = NULL, coordinates = c(-112.2124, 37.63281)){

  sf::sf_use_s2(FALSE)

  if(is.null(sf)){

    # Create a data frame with a column named 'geometry'
    df <- tibble::tibble(long = coordinates[1],
                         lat = coordinates[2])

    aoi <- sf::st_as_sf(df, coords = c("long", "lat"), crs = 4326)
  }

  if(is.null(coordinates)){

    aoi <- sf %>% sf::st_transform(4326)

  }

  aoi <- aoi %>%
    nhdplusTools::get_nhdplus(AOI = .)

  # Use the NHDPlus digital elevation model to find the nearest downslope
  # NHD flowline for any point in space (here, our point of interest)
  # Not working anymore -
  # trace <- get_raindrop_trace(aoi, direction = "up")
  #
  # "Snap" our site to the nearest NHD flowline feature
  # snap_point <- sf::st_sfc(sf::st_point(trace$intersection_point[[1]][1:2]),
  #                          crs=4326)
  #
  # Clip/split our catchment to only include the portion of the
  # catchment upstream of our site:
  # better_termination <- get_split_catchment(snap_point, upstream = F)[2,]

  # read in the complete NHD (in tabular form) to make for much more efficient nhd crawling.
  # This data in tabular form doesn't exist anywhere online that I know of... -_-
  nhd <- readr::read_csv('data/context_data/nhd_flow_network.csv')

  upstream <- nhdplusTools::get_UT(nhd, aoi$comid) %>% #upstream trace function in nhdplusTools
    tibble::as_tibble() %>%
    dplyr::rename(comid_list = value)  %>%
    dplyr::distinct(comid_list, .keep_all = TRUE)

  nhd_catch <-  upstream$comid_list %>%
    map(~nhdplusTools::get_nhdplus(comid = .,
                                   realization='catchment',
                                   t_srs = 4269)) %>%
    dplyr::bind_rows() %>%
    dplyr::distinct(featureid,.keep_all=TRUE) %>%
    dplyr::group_by(featureid) %>%
    dplyr::summarize()

  return(nhd_catch)

}



get_weekly_inspection_decision <- function(weekly_plot_data) {

  prompt_text <- "Which days do you want to verify (1/2/3/4/5/6/7)? \nCan be any combination of weekdays, with no other characters.\nex: 456\n"
  day_min <- min(weekly_plot_data$weekday)
  day_max <- max(weekly_plot_data$weekday)

  while (TRUE) {

    user_input_og <- readline(prompt = paste(prompt_text))

    user_input <- suppressWarnings(
      sort(na.omit(as.numeric(unique(unlist(strsplit(user_input_og, ""))))) [as.numeric(unique(unlist(strsplit(user_input_og, "")))) >= day_min & as.numeric(unique(unlist(strsplit(user_input_og, "")))) <= day_max])
      ) # 1 and 7 will be set to data min and max

    if (length(user_input) != 0) {
      return(user_input)
    }

    cat(user_input_og, " is not a valid input")
  }
}

#' @title Get working data decision
#'
#' @description
#' Get working directory and data decisions from user input
#'
#' @param prompt_text The prompt that the user is replying to
#'
#' @examples
#' # get_day_decision(prompt_text = "Would you like to (pass/fail/inspect) all data points for: ")

get_working_data_decision <- function() { # should add a layer level to this to prevent user from inspect sub-daily data

  prompt_text <- "Which directory are you working from? (pre/int): "

  while (TRUE) {
    user_input <- readline(prompt = paste(prompt_text))
    user_input <- tolower(user_input)

    if (user_input == "pre") {
      working_data <<- set_names(map(list.files(pre_verification_path, full.names = TRUE), readRDS), list.files(pre_verification_path))
      break()
    }

    if (user_input == "int") {
      working_data <<- set_names(map(list.files(intermediary_path, full.names = TRUE), readRDS), list.files(intermediary_path))
      break()
    }

    cat("Invalid input. Please enter one of the options from the following:\n\n")
    cat("pre = pre_verification_dir\n")
    cat("int = intermediary_dir\n")
  }

}

grab_confirmed_mWater_malfunctions <- function(field_notes){

    # Grab notes about sensor malfunction
    malfunction_records <- all_notes_cleaned %>%
      dplyr::filter(grepl("Sensor malfunction", visit_type, ignore.case = TRUE)) %>%
      dplyr::select(start_DT, site, crew, which_sensor_malfunction, malfunction_end_dt, notes = visit_comments)

    # write to csv
    write_csv(malfunction_records, "data/mWater_malfunciton_records.csv")

  }

grab_mWater_malfunction_notes <- function(mWater_api_data){

  # Grab notes about sensor malfunction
  malfunction_records <- mWater_api_data %>%
    filter(grepl("Sensor malfunction", visit_type, ignore.case = TRUE)) %>%
    select(start_DT, site, crew, which_sensor_malfunction, malfunction_end_dt, notes = visit_comments)

  #write to csv
  write_csv(malfunction_records, "data/mWater_malfunction_records.csv")

  parameters <- c("Battery Level",
                  "Baro",
                  "Chl-a Fluorescence",
                  "Depth",
                  "DO",
                  "External Voltage",
                  "ORP",
                  "pH",
                  "Specific Conductivity",
                  "Temperature",
                  "Turbidity")

  malfunction_records <- malfunction_records %>%
    # keep records relevant to {target} analysis
    select(start_DT, end_DT = malfunction_end_dt, site, parameter = which_sensor_malfunction, notes) %>%
    # match the text in the sensor column to the text in the target analysis
    separate_rows(parameter, sep = ", ") %>%
    mutate(
      parameter = case_when(
      parameter == "Chlorophyll a" ~ "Chl-a Fluorescence",
      parameter == "RDO" ~ "DO",
      parameter == "Conductivity" ~ "Specific Conductivity",
      .default = parameter
    ),
    site = case_when(
      site == "riverbluffs" ~ "river bluffs",
      .default = site
    )) %>%
    filter((is.na(parameter)) | (parameter %in% parameters))

  return(malfunction_records)

}

grab_mWater_sensor_malfunctions <- function(){

  # API Pull of mWater submitted notes

  # Grab API url from yml
  # Contact Sam Struthers if you need access
  creds <- yaml::read_yaml("src/mWater_collate/mWater_API.yml")
  api_url <- as.character(creds["url"])

  # Read in from API and tidy for downstream use

  # This is basic tidying of data set to:
  # correct datetime from UTC to Denver time (always MST)
  # correct columns where Other input is allowed (Site, visit type, photos downloaded, sensor malfunction)
  # Add rounded date time

  mal_notes <- readr::read_csv(url(api_url), show_col_types = FALSE) %>%
    dplyr::mutate(
      # start and end dt comes in as UTC -> to MST
      start_DT = lubridate::with_tz(lubridate::parse_date_time(start_dt, orders = c("%Y%m%d %H:%M:%S", "%m%d%y %H:%M", "%m%d%Y %H:%M", "%b%d%y %H:%M")), tz = "MST"),
      end_dt = lubridate::with_tz(lubridate::parse_date_time(end_dt, orders = c("%Y%m%d %H:%M:%S", "%m%d%y %H:%M", "%m%d%Y %H:%M", "%b%d%y %H:%M" )), tz = "MST"),
      malfunction_end_dt = with_tz(lubridate::parse_date_time(malfunction_end_dt, orders = c("%Y%m%d %H:%M:%S", "%m%d%y %H:%M", "%m%d%Y %H:%M", "%b%d%y %H:%M" )), tz = "MST"),
      date = as.Date(start_DT, tz = "MST"),
      start_time_mst = format(start_DT, "%H:%M"),
      sensor_pulled = as.character(sn_removed),
      sensor_deployed = as.character(sn_deployed),
      # If other is chosen, make site == other response
      site = ifelse(site == "Other (please specify)", tolower(stringr::str_replace_all(site_other, " ", "")), site),
      # When I changed the mWater survey, I accidentally introduced ??? in the place of Sensor Calibration option, fixing that here
      visit_type = dplyr::case_when(stringr::str_detect(visit_type, "\\?\\?\\?") ~ stringr::str_replace(string = visit_type,
                                                                                                        pattern =  "\\?\\?\\?",
                                                                                                        replacement = "Sensor Calibration or Check"),
                                    TRUE ~ visit_type),
      # Merging visit_type and visit type other
      visit_type = dplyr::case_when(stringr::str_detect(visit_type, "Other") ~ stringr::str_replace(string = visit_type,
                                                                                                    pattern =  "Other \\(please specify\\)",
                                                                                                    replacement = visit_type_other),
                                    TRUE ~ visit_type),
      # Merge sensor malfunction and sensor malfunction other
      which_sensor_malfunction = dplyr::case_when(stringr::str_detect(which_sensor_malfunction, "Other") ~ stringr::str_replace(string = which_sensor_malfunction,
                                                                                                                                pattern =  "Other \\(please specify\\)",
                                                                                                                                replacement = as.character(other_which_sensor_malfunction)),
                                                  TRUE ~ which_sensor_malfunction),
      # If other is chosen, make photos downloaded equal to response
      photos_downloaded = ifelse(photos_downloaded == "Other (please specify)", photos_downloaded_other, photos_downloaded),
      # Rounded start date time
      DT_round = lubridate::floor_date(start_DT, "15 minutes")) %>%
    # arrange by most recent visit
    dplyr::arrange(DT_round) %>%
    # Remove other columns
    dplyr::select(-c(photos_downloaded_other, visit_type_other, site_other, other_which_sensor_malfunction)) %>%
    dplyr::filter(grepl("sensor malfunction", visit_type, ignore.case = TRUE)) %>%
    dplyr::select(malfunction_start_dt = DT_round,  malfunction_end_dt, which_sensor_malfunction)

  return(mal_notes)

}

grab_mWater_sensor_notes <- function(mWater_api_data){

  # Sensor Notes

  # These are the notes that will be added to the QAQC workflow notes Most of the code in this chunk is to get the df to
  # match the one in the QAQC workflow It can be saved as a CSV or pulled directly into QAQC workflow
  # grab only notes where technician is interacting with sensor on site (excludes sensor malfunction notes)

  mWater_field_notes <- mWater_api_data %>%
    filter(grepl("Sensor",visit_type, ignore.case = TRUE) & !grepl("Sensor malfunction",visit_type, ignore.case = TRUE)) %>%
    # determining sonde employed status based on sensor_change
    mutate(sonde_employed = case_when(is.na(sensor_change)  ~ NA,
                                      sensor_change == "Swapped" ~ NA,
                                      sensor_change == "Pulled" ~ 1,
                                      sensor_change == "Deployed" ~ 0),
                                      #sensor_change %in% c("Swapped", "Deployed") ~ 1),

           #Sensor swapped notes
           sensor_swapped_notes = case_when(is.na(sensor_change)  ~ NA,
                                            sensor_change == "Pulled" &!is.na(sensor_pulled) ~ paste0("SN Removed: ", sensor_pulled),
                                            sensor_change == "Swapped" ~ paste0("SN Removed: ", sensor_pulled, " SN Deployed: ", sensor_deployed),
                                            sensor_change == "Deployed" ~ sensor_deployed),
           #Date/field season columns to match QAQC workflow
           DT_join = as.character(DT_round),
           field_season = year(DT_round),
           last_site_visit = DT_round,
           date = as.character(date)
    )%>%
    arrange(desc(DT_round))%>%
    #order columns in easily readable ordering
    select(site, crew, DT_round,sonde_employed,  sensors_cleaned, wiper_working, rdo_cap_condition, rdo_cap_replaced , ph_junction_replaced ,
           cals_performed, cal_report_collected , sensor_malfunction,sensor_pulled,sensor_deployed, sensor_swapped_notes,
           visit_type,start_time_mst,DT_join,  start_DT, end_dt,date,  visit_comments,photos_downloaded, field_season, last_site_visit)

  # back up to CSV
  # write_csv(sensor_notes, "data/mWater_sensor_field_notes.csv")

  # rm(all_notes_cleaned)
  return(mWater_field_notes)

}


hv_data_id <- function(loc_id, start_time = startdate, end_time = enddate, tz = timezone, token) {

  # convert the time to timestamp, convert to UTC for lookup in HydroVu
  start <- as.numeric(lubridate::with_tz(lubridate::ymd_hms(start_time, tz = tz), tzone = "UTC"))
  end <- as.numeric(lubridate::with_tz(lubridate::ymd_hms(end_time, tz = tz), tzone = "UTC"))

  # build the url
  url = "https://www.hydrovu.com/public-api/v1/locations/"
  url <- paste0(url, loc_id, "/data?endTime=", end, '&startTime=', start)

  req <- httr2::request(url)
  print(paste0('Trying site ', loc_id))
  try({
    resp <-  req %>% httr2::req_oauth_client_credentials(token) %>% httr2::req_perform()
    data <- list(resp %>% httr2::resp_body_json())
    h <- resp %>% httr2::resp_headers()

    while (!is.null(h[["X-ISI-Next-Page"]]))
    {
      resp <- req %>% httr2::req_headers("X-ISI-Start-Page" = h[["X-ISI-Next-Page"]]) %>%
        httr2::req_oauth_client_credentials(token) %>% httr2::req_perform()
      data <- c(data, list(resp %>% httr2::resp_body_json()))
      h <- resp %>% httr2::resp_headers()
    }

    # get the params and units
    params <- hv_names(token, return = "params")
    units <- hv_names(token, return = "units")

    # collapse the paginated date and clean up
    df <- purrr::map_dfr(data, flatten_page_params) %>%
      dplyr::mutate(timestamp = lubridate::with_tz(lubridate::as_datetime(timestamp, tz = "UTC"), tzone = tz),
                    Location = loc_id) %>%
      dplyr::inner_join(params, by = "parameterId") %>%
      dplyr::inner_join(units, by = "unitId") %>%
      dplyr::select(-parameterId, -unitId) %>%
      dplyr::arrange(Parameter, timestamp)

      return(df)
  })

}

#' Return the list of locations for the given client from HydroVu
#'
#' @param client a valid OAuth2 token such as returned from \code{hv_auth()}
#' @param url HydroVu url that lists the locations
#'
#' @return a dataframe listing all the locations visible to the client
#' @export
#'
#' @examples
#' \dontrun{
#' locs <- hv_locations(client)
#' }

hv_locations_all <- function(client,
                         url = "https://www.hydrovu.com/public-api/v1/locations/list") {

  req <- httr2::request(url)

  try({
  resp <-  req %>% httr2::req_oauth_client_credentials(client) %>% httr2::req_perform()
  locs <- list(resp %>% httr2::resp_body_json())
  h <- resp %>% httr2::resp_headers()

  while (!is.null(h[["X-ISI-Next-Page"]]))
  {
    resp2 <- req %>%
      httr2::req_headers("X-ISI-Start-Page" = h[["X-ISI-Next-Page"]]) %>%
      httr2::req_oauth_client_credentials(client) %>%
      httr2::req_perform()
    locs <- c(locs, list(resp2 %>% httr2::resp_body_json()))
    h <- resp2 %>% httr2::resp_headers()
  }
  # collapse the paginated date and clean up
  df <- flatten_df(locs) %>%
    select(-gps) %>%
    filter(!duplicated(.))
  return(df)
  })
}


#' @title intersensor_check
#'
#' @description
#' A function designed to reduce overflagging: if a slope violation occurs at the
#' same time as a slope violation in either depth or temperature, it is likely
#' not a sensor malfunction and instead a real product of the river.
#'
#' @param df An updated data frame with a `flag` column.
#'
#' @return A data frame with a `flag` column that has been updated with the
#' inter-sensor flag reduction step.
#'
#' @seealso [add_flag()]

intersensor_check <- function(df){

  # create a df of temperature for each site
  temperature <- df %>%
    data.table::data.table() %>%
    dplyr::select(DT_round, DT_join, parameter, mean, flag) %>%
    dplyr::filter(parameter == "Temperature") %>%
    dplyr::select(DT_join, Temperature = parameter, Temperature_flag = flag) %>%
    dplyr:: mutate(Temperature_front1 = dplyr::lead(Temperature_flag, n = 1),
                   Temperature_back1 = dplyr::lag(Temperature_flag, n = 1))

  # create a df of depth for each site
  depth <- df %>%
    data.table::data.table() %>%
    dplyr::select(DT_round, DT_join, parameter, mean, flag) %>%
    dplyr::filter(parameter == "Depth") %>%
    dplyr::select(DT_join, Depth = parameter, Depth_flag = flag) %>%
    dplyr:: mutate(Depth_front1 = dplyr::lead(Depth_flag, n = 1),
                   Depth_back1 = dplyr::lag(Depth_flag, n = 1))

  # add "temperature" and "depth" data columns to df:
  intersensors_checked <- df %>%
    dplyr::filter(!parameter %in% c("Depth", "Temperature")) %>%
    dplyr::left_join(., temperature, by = "DT_join") %>%
    dplyr::left_join(., depth, by = "DT_join") %>%
    # If either the depth or temperature have the same flag as a given parameter
    # identified at the same time (or one obs before/after), tag it
    dplyr::mutate(intersensored = dplyr::case_when(grepl("slope violation", flag) &
                                                     (grepl("slope violation", Depth_flag)   | grepl("slope violation", Temperature_flag)   |
                                                        grepl("slope violation", Depth_front1) | grepl("slope violation", Temperature_front1) |
                                                        grepl("slope violation", Depth_back1)  | grepl("slope violation", Temperature_back1)
                                                     ) ~ TRUE)) %>%
    dplyr::mutate(flag = ifelse(is.na(intersensored), flag, stringr::str_replace(flag, "slope violation", "")))

  final_checked_data <- df %>%
    dplyr::filter(parameter %in% c("Depth", "Temperature")) %>%
    # After using the temp and depth slope flags, remove that flagging entirely
    # from those parameters. We have yet to find an instance of the slope flag
    # capturing "fake" spikes in either of those data sets:
    dplyr::mutate(flag = stringr::str_replace(flag, "slope violation", "")) %>%
    dplyr::bind_rows(., intersensors_checked) %>%
    dplyr::select(-c(Depth, Depth_flag, Temperature, Temperature_flag))

  return(final_checked_data)

}

#' @title Load and tidy mWater field notes
#'
#' @description A function that uploads and cleans the field notes submitted to mWater.
#'
#' @param creds A .yml file with necessary credentials for accessing the field notes. Contact Sam Struthers if you need access.
#'
#' @return A dataframe with the field notes.

load_mWater_notes <- function(creds = yaml::read_yaml("creds/mWaterCreds.yml")){

  # API Pull of mWater submitted notes

  # Grab API url from yml
  # Contact Sam Struthers if you need access
  api_url <- as.character(creds["url"])

  # Read in from API and tidy for downstream use

  # This is basic tidying of data set to:
  # correct datetime from UTC to Denver time (MST)
  # correct columns where Other input is allowed (Site, visit type, photos downloaded, sensor malfunction)
  # Add rounded date time

  all_notes_cleaned <- readr::read_csv(url(api_url), show_col_types = FALSE) %>%
    dplyr::mutate(
      # start and end dt comes in as UTC -> to MST
      start_DT = lubridate::with_tz(lubridate::parse_date_time(start_dt, orders = c("%Y%m%d %H:%M:%S", "%m%d%y %H:%M", "%m%d%Y %H:%M", "%b%d%y %H:%M" )), tz = "MST"),
      end_dt = lubridate::with_tz(lubridate::parse_date_time(end_dt, orders = c("%Y%m%d %H:%M:%S", "%m%d%y %H:%M", "%m%d%Y %H:%M", "%b%d%y %H:%M" )), tz = "MST"),
      malfunction_end_dt = lubridate::with_tz(lubridate::parse_date_time(malfunction_end_dt, orders = c("%Y%m%d %H:%M:%S", "%m%d%y %H:%M", "%m%d%Y %H:%M", "%b%d%y %H:%M" )), tz = "MST"),
      date = as.Date(start_DT, tz = "MST"),
      start_time_mst = format(start_DT, "%H:%M"),
      sensor_pulled = as.character(sn_removed),
      sensor_deployed = as.character(sn_deployed),
      # If other is chosen, make site == other response
      site = ifelse(site == "Other (please specify)", tolower(stringr::str_replace_all(site_other, " ", "")), site),
      # When I changed the mWater survey, I accidentally introduced ??? in the place of Sensor Calibration option, fixing that here
      visit_type = dplyr::case_when(stringr::str_detect(visit_type, "\\?\\?\\?") ~ stringr::str_replace(string = visit_type,
                                                                               pattern =  "\\?\\?\\?",
                                                                               replacement = "Sensor Calibration or Check"),
                             TRUE ~ visit_type),
      # Merging visit_type and visit type other
      visit_type = dplyr::case_when(stringr::str_detect(visit_type, "Other") ~ stringr::str_replace(string = visit_type,
                                                                           pattern =  "Other \\(please specify\\)",
                                                                           replacement = visit_type_other),
                             TRUE ~ visit_type),
      # Merge sensor malfunction and sensor malfunction other
      which_sensor_malfunction = dplyr::case_when(stringr::str_detect(which_sensor_malfunction, "Other") ~ stringr::str_replace(string = which_sensor_malfunction,
                                                                                                       pattern =  "Other \\(please specify\\)",
                                                                                                       replacement = as.character(other_which_sensor_malfunction)),
                                           TRUE ~ which_sensor_malfunction),
      # If other is chosen, make photos downloaded equal to response
      photos_downloaded = ifelse(photos_downloaded == "Other (please specify)", photos_downloaded_other, photos_downloaded),
      # Rounded start date time
      DT_round = lubridate::floor_date(start_DT, "15 minutes")) %>%
    # arrange by most recent visit
    dplyr::arrange(DT_round)%>%
    # Remove other columns
    dplyr::select(-c(photos_downloaded_other,visit_type_other, site_other, other_which_sensor_malfunction ))

  return(all_notes_cleaned)

}

#' @title Load and tidy old field notes
#'
#' @description A function that uploads and cleans the field notes excel file. This function adds datetime
#' columns to the field notes dataframe and filters out field notes where the sensor
#' was not handled.
#'
#' @param filepath A file path to the raw field notes.
#'
#' @return A dataframe with the field notes.
#'
#' @examples
#' clean_old_field_notes(filepath = "data/sensor_field_notes.xlsx")

load_old_field_notes <- function(filepath){

  raw_field_notes <- readxl::read_excel(filepath)

  field_notes <- raw_field_notes %>%
    dplyr::mutate(start_DT = lubridate::ymd_hm(paste(date, start_time_mst), tz = "MST")) %>%
    dplyr::mutate(
      DT_round = lubridate::floor_date(start_DT, "15 minutes"),
      DT_join = as.character(DT_round),
      site = tolower(site),
      field_season = lubridate::year(DT_round),
      last_site_visit = DT_round) %>%
    dplyr::arrange(site, DT_round) %>%
    # rename instances of old names:
    dplyr::mutate(site = ifelse(site == "rist", "tamasag",
                         ifelse(site == "elc", "boxelder", site))) %>%
    # `sonde_employed` determines if the sonde is deployed or not. 0 = sonde deployed, 1 = sonde is not deployed
    mutate(sonde_employed = dplyr::case_when(!is.na(sensor_pulled) & !is.na(sensor_deployed) ~ 0,
                                             !is.na(sensor_pulled) & is.na(sensor_deployed) ~ 1,
                                             is.na(sensor_pulled) & !is.na(sensor_deployed) ~ 0,
                                             is.na(sensor_pulled) & is.na(sensor_deployed) ~ NA),
           end_dt  = as.POSIXct(NA, tz = "MST")) %>%
    # remove field dates where sensor was not handled:
    dplyr::filter(grepl("Sensor Cleaning or Check|Sensor Calibration", visit_type, ignore.case = TRUE))

  return(field_notes)

}

#' @title Generate Threshold Table for QAQC
#'
#' @description
#' A function designed to generate a threshold table for QAQC. This table
#' contains the thresholds for the mean, slope_behind, and standard deviation
#' of the mean for each site and season.
#'
#' @param df A data frame with a `flag` column.
#'
#' @return A data frame with the thresholds for the mean, slope_behind, and
#' standard deviation of the mean for each site and season.
#'
#' @examples
#' make_threshold_table(df = all_data_flagged$`archery-Actual Conductivity`)

make_threshold_table <- function(df){

  # sensor_malfunction_notes <- grab_mWater_malfunction_notes(mWater_api_data = load_mWater_notes())

  slope_down <- df %>%
    # REMOVE DATA WE KNOW TO BE ERRONEOUS:
    add_field_flag() %>%
    # incoroporate back-calibrated values for development of thresholds:
    fix_calibration() %>%
    add_malfunction_flag(df = ., malfunction_records = sensor_malfunction_notes) %>%
    add_spec_flag() %>%
    # Filter to remove times when the sonde was under freexzing conditions or not submerged (depth <= 0)
    add_frozen_flag() %>%
    # fix_depth_cal() %>%
    add_unsubmerged_flag() %>%
    # remove times when sonde was moved (to avoid crazy slopes)
    add_depth_shift_flag() %>%
    # Filter to only un-flagged data, and data that was collected when sonde
    # was actually deployed
    filter(is.na(flag) | sonde_employed != 1) %>%
    #... and when the sonde wasn't moved
    filter(is.na(depth_change) | depth_change != "sonde moved") %>%
    # Get threshold for negative slope data
    filter(slope_behind < 0) %>%
    group_by(season) %>%
    summarize(f_slope_behind_01 = quantile(slope_behind, 0.01, na.rm = TRUE))

  slope_up <- df %>%
    # REMOVE DATA WE KNOW TO BE ERRONEOUS:
    add_field_flag() %>%
    fix_calibration() %>%
    add_malfunction_flag(df = ., malfunction_records = sensor_malfunction_notes) %>%
    add_spec_flag() %>%
    # Filter to remove times when the sonde was under freexzing conditions or not submerged (depth <= 0)
    add_frozen_flag() %>%
    # fix_depth_cal() %>%
    add_unsubmerged_flag() %>%
    # remove times when sonde was moved (to avoid crazy slopes)
    add_depth_shift_flag() %>%
    # Filter to only un-flagged data, and data that was collected when sonde
    # was actually deployed
    filter(is.na(flag) | sonde_employed != 1) %>%
    #... and when the sonde wasn't moved
    filter(is.na(depth_change) | depth_change != "sonde moved") %>%
    # Get threshold for positive slope data
    filter(slope_behind > 0) %>%
    group_by(season) %>%
    summarize(f_slope_behind_99 = quantile(slope_behind, 0.99, na.rm = TRUE))

  good_data_stats <- df %>%
    # REMOVE DATA WE KNOW TO BE ERRONEOUS:
    add_field_flag() %>%
    fix_calibration() %>%
    add_malfunction_flag(df = ., malfunction_records = sensor_malfunction_notes) %>%
    add_spec_flag() %>%
    # Filter to remove times when the sonde was under freexzing conditions or not submerged (depth <= 0)
    add_frozen_flag() %>%
    # fix_depth_cal() %>%
    add_unsubmerged_flag() %>%
    # remove times when sonde was moved (to avoid crazy slopes)
    add_depth_shift_flag() %>%
    # Filter to only un-flagged data, and data that was collected when sonde
    # was actually deployed
    filter(is.na(flag) | sonde_employed != 1) %>%
    #... and when the sonde wasn't moved
    filter(is.na(depth_change) | depth_change != "sonde moved") %>%
    group_by(season) %>%
    # join our slope data thresholds:
    left_join(slope_up, by = "season") %>%
    left_join(slope_down, by = "season") %>%
    # develop other thresholds across all data
    mutate(f01 = quantile(mean, 0.01, na.rm = TRUE),
           f99 = quantile(mean, 0.99, na.rm = TRUE)) %>%
           # f_slope_behind_01 = slope_down, #quantile(slope_behind, 0.01, na.rm = TRUE),
           # f_slope_behind_99 = slope_up) %>% #quantile(slope_behind, 0.99, na.rm = TRUE)) %>%
    # THEN, GET STANDARD DEVIATION OF ONLYYYY VALUES WITHIN THE 1-99th PERCENTILE OF THAT GOOD DATA:
    filter(mean > f01 & mean < f99) %>%
    # SD is the ONLY statistic that uses this winnowed-down data set in its development.
    # All else use the full, "good" data set.
    summarize(site = paste0(unique(site)),
              parameter = paste0(unique(parameter)),
              t_mean01 = as.numeric(paste0(unique(f01))),
              t_mean99 = as.numeric(paste0(unique(f99))),
              t_slope_behind_01 = as.numeric(paste0(unique(f_slope_behind_01))),
              t_slope_behind_99 = as.numeric(paste0(unique(f_slope_behind_99))),
              # This stat is useless. Should remove eventually.
              t_sd_0199 = sd(mean, na.rm = T))

  return(good_data_stats)

}

#' @title Munge API data for QAQC workflow
#'
#' @description
#' A function designed to munge the raw API data for the QAQC workflow.
#'
#' @param api_path Path where the raw API data lives.
#' @param network Options include "csu", "virridy"
#' @param require An upstream dependency (i.e., alternative to `{necessary object}`).
#'
#' @return A dataframe with the munged API data.
#'
#' @examples
# munge_api_data(api_path = "data/api/incoming/")

munge_api_data <- function(api_path, network = "csu", require = NULL) {

  api_data <- list.files(path = api_path, full.names = TRUE, pattern = "*.csv") %>%
    purrr::map_dfr(~data.table::fread(.) %>%
                     dplyr::select(-id)) %>%
    # remove overlapping API-pull data
    dplyr::distinct()

  if(network %in% c("csu", "CSU")){

  api_data <- api_data %>%
    # remove VuLink data
    dplyr::filter(!grepl("vulink", name, ignore.case = TRUE)) %>% # these will throw us errors if there is none?
    # remove Virridy data
    dplyr::filter(!grepl("virridy", name, ignore.case = TRUE)) %>%
    dplyr::select(-name) %>%
    # Convert UTC (as it is sent from HydroVU API) to MST:
    dplyr::mutate(DT = lubridate::as_datetime(timestamp, tz = "UTC")) %>%
    dplyr::mutate(DT = lubridate::with_tz(DT, tzone = "MST"),
           DT_round = lubridate::round_date(DT, "15 minutes"),
           DT_join = as.character(DT_round),
           site = tolower(site)) %>%
    # These sites will be considered the same site for this workflow
    dplyr::mutate(site = ifelse(site == "rist", "tamasag",
                  ifelse(site == "elc", "boxelder", site))) %>%
    # Lastly, we swapped Boxelder's sonde out for Rist's late in 2022:
    dplyr::mutate(site = ifelse(site == "tamasag" & DT > lubridate::ymd("2022-09-20", tz = "MST") & DT < lubridate::ymd("2023-01-01", tz = "MST"), "boxelder", site)) %>%
    dplyr::distinct(.keep_all = TRUE)
  }

  if(network %in% c("virridy", "Virridy")){

    api_data <- api_data %>%
      # remove VuLink data
      dplyr::filter(!grepl("vulink", name, ignore.case = TRUE)) %>% # these will throw us errors if there is none?
      dplyr::select(-name) %>%
      # Convert UTC (as it is sent from HydroVU API) to MST:
      dplyr::mutate(DT = lubridate::as_datetime(timestamp, tz = "UTC")) %>%
      dplyr::mutate(DT = lubridate::with_tz(DT, tzone = "MST"),
                    DT_round = lubridate::round_date(DT, "15 minutes"),
                    DT_join = as.character(DT_round),
                    site = tolower(site)) %>%
      dplyr::distinct(.keep_all = TRUE)
    }

  return(api_data)

}

#' @title Network Check
#'
#' @description
#' This function performs a network check on a given data frame, flagging potential
#' issues in the data based on upstream and downstream sites.
#'
#' @param df A site-parameter data frame that has gone through the initial flagging
#' process.
#' @param network Whether the network check is happening across the Virridy sites or CSU sites.
#' @return A modified data frame flags that have been altered based on the network check.
#'
#' @examples
#' network_check(df = all_data_flagged$`archery-Actual Conductivity`)

network_check <- function(df, network = "csu") {

  df <- df

  # get the site name from the site column in the df of interest
  site_name <- unique(na.omit(df$site))
  # get the parameter from the parameter column in the df of interest
  parameter_name <- unique(na.omit(df$parameter))

  # vector of sites in the order that they are in spatially
  # some sites have some funkiness going on

  sites_order <- c("tamasag", # rist
                   "legacy",
                   "lincoln",
                   "timberline",
                   "prospect",
                   "boxelder", # elc
                   "archery",
                   "river bluffs")

  width_fun = ifelse(site_name == "tamasag", 17, # 2 hours before/after
              ifelse(site_name == "legacy", 17,
              ifelse(site_name == "lincoln", 17,
              ifelse(site_name == "timberline", 17,
              ifelse(site_name == "prospect", 17,
              ifelse(site_name == "boxelder", 17,
              ifelse(site_name == "archery", 17,
              ifelse(site_name == "river bluffs", 17, NA))))))))

  if(network %in% c("virridy", "Virridy")){

    sites_order <-  c("joei",
                      "cbri",
                      "chd",
                      "pfal",
                      "pbd",
                      "tamasag",
                      "legacy",
                      "lincoln",
                      "timberline",
                      #"springcreek",
                      "prospect",
                      "boxelder",
                      #boxcreek,"
                      "archery",
                      "river bluffs")

    width_fun = ifelse(site_name == "joei", 17, # 2 hours before/after
                ifelse(site_name == "cbri", 17,
                ifelse(site_name == "chd", 17,
                ifelse(site_name == "pfal", 17,
                ifelse(site_name == "pbd", 17,
                ifelse(site_name == "sfm", 17,
                ifelse(site_name == "lbea", 17,
                ifelse(site_name == "penn", 17,
                ifelse(site_name == "tamasag", 17,
                ifelse(site_name == "legacy", 17,
                ifelse(site_name == "lincoln", 17,
                ifelse(site_name == "timberline", 17,
                ifelse(site_name == "timberline virridy", 17,
                ifelse(site_name == "springcreek", 17,
                ifelse(site_name == "prospect", 17,
                ifelse(site_name == "prospect virridy", 17,
                ifelse(site_name == "boxelder", 17,
                ifelse(site_name == "boxcreek", 17,
                ifelse(site_name == "archery", 17,
                ifelse(site_name == "archery virridy", 17,
                ifelse(site_name == "river bluffs", 17, NA)))))))))))))))))))))

    if(site_name %in% c("penn", "sfm", "lbea")){

      sites_order <- c("penn",
        "sfm",
        "lbea")

    }

    if(site_name == "springcreek"){

      sites_order <- c("timberline virridy",
                       "springcreek",
                       "prospect virridy")

    }

    if(site_name == "boxcreek"){

      sites_order <- c("boxelder virridy",
                       "boxcreek",
                       "archery virridy")

    }

  }

  # determining the index for the site of interest.
  site_index <- which(sites_order == sites_order[grep(gsub(" virridy", "", site_name), sites_order, ignore.case = TRUE)])

  # Generating df name to pull from df_list list
  site_param <- paste0(site_name, "-", parameter_name)

  prev_site_df <- tibble(DT_round = NA)
  next_site_df <- tibble(DT_round = NA)

  tryCatch({
    previous_site <- paste0(sites_order[site_index-1],"-",parameter_name)
    prev_site_df <- intersensor_checks[[previous_site]] %>%
      select(DT_round, site_up = site, flag_up = flag) %>%
      data.table()},
    error = function(err) {
      cat("No upstream site.\n")})

  tryCatch({
    next_site <- paste0(sites_order[site_index+1],"-",parameter_name)
    next_site_df <- intersensor_checks[[next_site]] %>%
      select(DT_round, site_down = site, flag_down = flag) %>%
      data.table()},
    error = function(err) {
      cat("No downstream site.\n")})


  join <- df %>%
    left_join(., prev_site_df, by = "DT_round") %>%
    left_join(., next_site_df, by = "DT_round")

  if(!("flag_down" %in% colnames(join))) {join$flag_down <- NA}
  if(!("flag_up" %in% colnames(join))) {join$flag_up <- NA}
  if(!("site_down" %in% colnames(join))) {join$site_down <- NA}
  if(!("site_up" %in% colnames(join))) {join$site_up <- NA}


  # Define a function to check if a given 2-hour window has any instances of the same word
  check_2_hour_window_fail <- function(x) {
    sum(x) >= 1
  }

  df_test <- join %>%
    # No upstream/downstream flag = 0
    mutate(flag_binary = ifelse(
      (is.na(flag_up) | grepl("drift|DO interference|repeat|sonde not employed|frozen|unsubmerged|missing data|site visit|sv window|sensor malfunction|sonde burial|sensor biofouling|improper level cal|sonde moved", flag_up)) &
        (is.na(flag_down) | grepl("drift|DO interference|repeat|sonde not employed|frozen|unsubmerged|missing data|site visit|sv window|sensor malfunction|sonde burial|sensor biofouling|improper level cal|sonde moved", flag_down)), 0, 1)) %>%
    mutate(overlapping_flag = zoo::rollapply(flag_binary, width = width_fun, FUN = check_2_hour_window_fail, fill = NA, align = "center")) %>%
    # If there is a flag (flags associated with spikes in concentration or funkiness like that), and there is also a flag up/downstream at the same time (2 hour window) it is likely a real
    # WQ event and should therefore not be considered "poor" data:
    mutate(auto_cleaned_flag = ifelse(!is.na(flag) & !grepl("drift|DO interference|repeat|sonde not employed|frozen|unsubmerged|missing data|site visit|sv window|sensor malfunction|sonde burial|sensor biofouling|improper level cal|sonde moved", flag) & overlapping_flag == TRUE, NA, flag)) %>%
    select(-c(flag_up, flag_down, site_up, site_down, flag_binary, overlapping_flag))

  # df_test <- join %>%
  #   mutate(flag_binary = ifelse(#grepl("slope|suspect", flag) &
  #     (is.na(flag_up) | grepl("seasonal range|repeat|sonde not employed|missing data|site visit|sv window", flag_up)) &
  #       (is.na(flag_down) | grepl("seasonal range|repeat|sonde not employed|missing data|site visit|sv window", flag_down)), 0, 1)) %>%
  #   #arrange(timestamp) %>%
  #   mutate(overlapping_flag = zoo::rollapply(flag_binary, width = width_fun, FUN = check_2_hour_window_fail, fill = NA, align = "center")) %>%
  #   mutate(cleaner_flag = ifelse(!is.na(flag) & !grepl("seasonal range|repeat|sonde not employed|missing data|site visit|sv window", flag) & overlapping_flag == TRUE, NA, flag)) %>%
  #   select(-c(flag_up, flag_down, site_up, site_down, flag_binary, overlapping_flag))

  return(df_test)

}

# Function to check for package installation, then install (if necessary) and load libraries.
# Adapted from code developed by Caitlin Mothes, PhD.


# fill in with packages that need to be loaded:
# packages <- c('tidyverse',
#               'sf')

package_loader <- function(x) {
  if (x %in% installed.packages()) {
    library(x, character.only = TRUE)
  } else {
    install.packages(x)
    library(x, character.only = TRUE)
  }
}

# lapply(packages, package_loader)

# generate_flag_report <- function(df) {
#
#   # Extract the site and parameter from the df_name
#   site <- unique(na.omit(df$site))
#   parameter <- unique(na.omit(df$parameter))
#
#   list_of_flags <- c("sonde not employed", # add_field_flag()
#                      "site visit", # add_field_flag()
#                      "sv window", # add_field_flag()
#                      "sensor malfunction", # add_malfunction_flag()
#                      "outside of sensor specification range", # add_spec_flag()
#                      "outside of seasonal range", # add_seasonal_flag()
#                      "slope violation", # add_seasonal_flag()
#                      "outside sd range", # add_seasonal_flag()
#                      "repeated value", # add_repeat_flag()
#                      "missing data", # add_na_flag()
#                      "suspect data") # add_suspect_flag()
#
#   # check these
#   sans_na_flags <- "^(missing data|
#                       sonde not employed;\\nmissing data|
#                       missing data;\\nsuspect data|
#                       sonde not employed;\\nmissing data;\\nsuspect data|
#                       site visit;\\nmissing data;\\nsuspect data|
#                       sv window;\\nmissing data;\\nsuspect data)$"
#
#   # summarize total data points
#   total_observations <- df %>%
#     summarise(n_total = n_distinct(DT_round)) %>%
#     pull(n_total)
#
#   # summarize total data points sans missing data
#   total_observations_1 <- df %>%
#     # filter out when flag has only missing data or only sonde not employed and missing data
#     filter(!str_detect(flag, sans_na_flags)) %>%
#     summarise(n_total = n_distinct(DT_round)) %>%
#     pull(n_total)
#
#   # summarize total days
#   total_observations_dates <- df %>%
#     group_by(date = format(DT_round, "%m-%d-%Y")) %>%
#     summarize(n_total = nrow(date)) %>%
#     nrow()
#
#   # summarize total days sans missing data
#   total_observations_dates_1 <- df %>%
#     filter(!str_detect(flag, sans_na_flags)) %>%
#     group_by(date = format(DT_round, "%m-%d-%Y")) %>%
#     summarize(n_total = nrow(date)) %>%
#     nrow()
#
#   row_list <- list()
#   for (i in list_of_flags) {
#
#     # summarize flagged data points
#     flagged_observations <- df %>%
#       filter(str_detect(flag, i)) %>%
#       summarise(n_flag = n_distinct(DT_round)) %>%
#       pull(n_flag)
#     # summarize flagged data points
#     flagged_observations_1 <- df %>%
#       filter(!str_detect(flag, sans_na_flags)) %>%
#       filter(str_detect(flag, i)) %>%
#       summarise(n_flag = n_distinct(DT_round)) %>%
#       pull(n_flag)
#     # summarize percent data points that are flagged
#     percent_flagged <- flagged_observations/total_observations
#     # summarize percent data points that are flagged sans missing data
#     percent_flagged_1 <- flagged_observations_1/total_observations_1
#
#     # summarize flagged days
#     flagged_observations_dates <- df %>%
#       filter(str_detect(flag, i)) %>%
#       group_by(date = format(DT_round, "%m-%d-%Y")) %>%
#       summarize(n_total = nrow(date)) %>%
#       nrow()
#     # summarize flagged days
#     flagged_observations_dates_1 <- df %>%
#       filter(!str_detect(flag, sans_na_flags)) %>%
#       filter(str_detect(flag, i)) %>%
#       group_by(date = format(DT_round, "%m-%d-%Y")) %>%
#       summarize(n_total = nrow(date)) %>%
#       nrow()
#     # summarize percent days that are flagged
#     percent_flagged_dates <- flagged_observations_dates/total_observations_dates
#     # summarize percent days that are flagged
#     percent_flagged_dates_1 <- flagged_observations_dates_1/total_observations_dates_1
#
#     # creating a row with the information
#     calculated_values <- tibble(
#       # metadata
#       site = site,
#       parameter = parameter,
#       flag = i,
#       # data points
#       data_points_flagged_percentage = format(round(percent_flagged * 100, 2), nsmall = 2),
#       data_points_flagged = flagged_observations,
#       total_data_points = total_observations,
#       data_points_flagged_percentage_sans_na = format(round(percent_flagged_1 * 100, 2), nsmall = 2),
#       data_points_flagged_sans_na = flagged_observations_1,
#       total_data_points_sans_na = total_observations_1,
#       # dates
#       dates_flagged_percentage = format(round(percent_flagged_dates * 100, 2), nsmall=2),
#       dates_flagged = flagged_observations_dates,
#       total_dates = total_observations_dates,
#       dates_flagged_percentage_sans_na = format(round(percent_flagged_dates_1 * 100, 2), nsmall=2),
#       dates_flagged_sans_na = flagged_observations_dates_1,
#       total_dates_sans_na = total_observations_dates_1
#       )
#
#     row_list[[i]] <- calculated_values
#   }
#
#   #calculated_df <- bind_cols(row_list)
#   return(bind_rows(row_list))
#
# }


photo_plotter <- function(index = 200, output_folder){

  param_1_label <- paste0(parameters[1], " (", param_unit[1],")")
  param_2_label <- paste0(parameters[2], " (", param_unit[2],")")

# get the adjustment, breaks etc so axes look nice
  param_1_max <- max(wq_tl[[parameters[1]]], na.rm = T)
  param_1_max_int <- as.integer(max(wq_tl[[parameters[1]]], na.rm = T))+1

  adjustment <- max(wq_tl[[parameters[2]]], na.rm = T) / param_1_max

  brk <- ifelse(param_1_max_int < 6, 1,2)

# Bounds for the y-axis

  lower_bound <- 0
  upper_bound <-  max(wq_tl[[parameters[2]]], na.rm = T)


  #This is the index of the image for the background for the individual photo
  simul = wq_tl[index,]
  # this is all the data before the image (ie previous photos from the timelapse)
  upto = wq_tl[1:index,]
  #read the image for the background
  photo_bg <- readJPEG(simul$filename)
  #create an individual image
  back <- ggplot() +
    #plot the photo
    annotation_custom(rasterGrob(photo_bg,
                                 width = unit(1,"npc"),
                                 height = unit(1,"npc")),
                      -Inf, Inf, -Inf, Inf)

  #plot the data for this image (includes all the preivous data)
  inset <- ggplot() +
    geom_ribbon(data = upto, aes(x = DT_round,
#to do: need to figure out adjustment factors for each parameter to plot correctly
                                 y = .data[[parameters[1]]] *adjustment,
                                 ymin = 0,
                                 ymax = .data[[parameters[1]]] * adjustment),
                color = "white",
                fill = "white",
                #linetype = "dash",
                alpha = 0.75)+
    geom_path(data = upto, aes(x = DT_round, y = .data[[parameters[2]]]),
              color = "#F34646", size=2) +
    geom_point(data = simul, aes(x = DT_round, y = .data[[parameters[2]]]),
               color = "#F34646")+
    # ylim(min(wq_tl[[parameters[2]]], na.rm = T),
    #      max(wq_tl[[parameters[2]]], na.rm = T))+
    xlim(min(wq_tl$DT_round, na.rm = T),
         max(wq_tl$DT_round, na.rm = T))+
    scale_y_continuous(name = param_2_label, limits = c(0, upper_bound),
                       sec.axis = sec_axis(trans = ~./adjustment ,
                                           name = param_1_label,
                                           breaks = seq(0,param_1_max_int,brk)))+
    dark_theme_light(base_size = 10) +
    theme(axis.title.y.right = element_text(color="white"),axis.text.y.right = element_text(color = "white"), axis.title.y.left = element_text(color="#F34646"), axis.text.y.left = element_text(color = "#F34646")) +
    theme(
      panel.background = element_rect(fill = "transparent"),
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.background = element_rect(fill = "transparent"),
      legend.box.background = element_rect(fill = "transparent"),
      panel.border = element_blank(),
      axis.line = element_line(color = 'gray80'),
    ) +
    ylab(param_2_label) +
    xlab('')

  gp1 <- back +
    inset_element(inset,
                  left = 0.16,
                  bottom = 0.15,
                  right = 0.9,
                  top = 0.6)
  #print(gp1)
  ggsave(paste0(folder_path, index, ".png"),
         width = 1920,
         height = 1080,
         units = "px")
  #dev.copy(gp1,paste0('data/timelapse_photos/vid_image/', max(upto$rowid), ".png"))
  #dev.off()
}

#test functions
#photo_plotter(2, "data/timelapse_photos/sjs_test/")
#map(1:nrow(wq_tl), ~photo_plotter(.x, "data/timelapse_photos/sjs_test/"))

plotter <- function(site_name = "boxelder",
                    start_date = "2023-03-04",
                    end_date = "2023-07-07",
                    title = "wack data title",
                    parameter = "all") {

  # Plotting options: all, Specific Conductivity, Turbidity, Chl-a Fluorescence, ORP, Depth, pH, DO

  sites_order <- c("tamasag", # rist
                   "legacy",
                   "lincoln",
                   "timberline",
                   "prospect",
                   "boxelder", # elc
                   "archery",
                   "river bluffs")

  # determining the index for the site of interest.
  site_index <- which(sites_order == site_name)
  previous_site <- sites_order[site_index-1]
  next_site <- sites_order[site_index+1]

  # Pull in WQ data from Manners Bridge (in the canyon):
  contrail <- list.files("data/context_data/contrail", full.names = TRUE) %>%
    map_dfr(~fread(.)) %>%
    mutate(date = as_date(ymd_hms(Reading)),
           Reading = floor_date(ymd_hms(Reading), "hour")) %>%
    filter(date >= start_date & date <= end_date)

  # Pull in a list of precip monitoring sites in each sonde's watershed:
  rain_list <- fread("data/context_data/site_link.csv") %>%
    filter(site == site_name) %>% pull(`Sensor Name`)

  # Pull in raw precip data
  rain_data <- fread("data/context_data/fc_rain.csv") %>%
    filter(`Sensor Name` %in% rain_list) %>%
    mutate(date = as_date(mdy_hms(Timestamp)),
           DT = floor_date(mdy_hms(Timestamp), unit = "hour")) %>%
    filter(date >= start_date & date <= end_date) %>%
    group_by(DT) %>%
    summarize(total_ws = sum(`Incremental Rainfall (in)`, na.rm = TRUE))

  # Grab the largest amount of rain across the sonde's watershed at an
  # hourly timestep:
  rain_stat <- fread("data/context_data/fc_rain.csv") %>%
    filter(`Sensor Name` %in% rain_list) %>%
    mutate(date = as_date(mdy_hms(Timestamp)),
           DT = floor_date(mdy_hms(Timestamp), unit = "hour")) %>%
    group_by(DT) %>%
    summarize(total_ws = sum(`Incremental Rainfall (in)`, na.rm = TRUE)) %>%
    filter(total_ws == max(total_ws))

  # Pull in flagged PSN WQ data:
  df <- fread('data/flagged/all_data_flagged_plotter.csv') %>%
    filter(!grepl("sonde not employed", final_cleaned_flag)) %>%
    mutate(date = as_date(DT_round)) %>%
    # filter to dates of interest:
    dplyr::filter((date) >= as_date(start_date) & (date) <= as_date(end_date))

  # site of interest data:
  site_df <- df %>%
    filter(!is.na(mean)) %>%
    filter(site == site_name)

  # downstream of site of interest (if there is one):
  downstream_df <- site_df
  try(downstream_df <- df %>%
        filter(!grepl("site visit|sv window", final_cleaned_flag)) %>%
        filter(site == next_site))

  # upstream of site of interest (if there is one):
  upstream_df <- site_df
  try(upstream_df <- df %>%
        filter(!grepl("site visit|sv window", final_cleaned_flag)) %>%
        filter(site == previous_site))

# All the plotz
  depth <- ggplot() +
    geom_line(data = filter(contrail, Unit == "ft"), aes(Reading, Value * 0.3048), color = "lightblue") +
    geom_line(data = filter(downstream_df, parameter == "Depth"), aes(DT_round, mean), cex = 0.2, alpha = 1, color = "lightgrey") +
    geom_line(data = filter(upstream_df, parameter == "Depth"), aes(DT_round, mean), cex = 0.2, alpha = 1, color = "darkgrey") +
    geom_line(data = filter(site_df, parameter == "Depth"), aes(DT_round, mean), color = "black") +
    geom_point(data = filter(site_df, parameter == "Depth" & grepl("site visit", final_cleaned_flag)), aes(DT_round, mean), color = "#E70870", cex = 1.5) +
    geom_point(data = filter(site_df, parameter == "Depth" & grepl("sonde moved", final_cleaned_flag)), aes(DT_round, mean), color = "#FFCA3A", cex = 0.75) +
    theme_bw() +
    ylim(min = min(filter(site_df, parameter == "Depth")$mean), max = max(filter(site_df, parameter == "Depth")$mean)) +
    ylab("Depth m") + xlab("") + ggtitle(title)

  temp <- ggplot() +
    geom_line(data = filter(contrail, Unit == "C"), aes(Reading, Value), color = "lightblue") +
    geom_line(data = filter(downstream_df, parameter == "Temperature"), aes(DT_round, mean), cex = 0.2, alpha = 1, color = "lightgrey") +
    geom_line(data = filter(upstream_df, parameter == "Temperature"), aes(DT_round, mean), cex = 0.2, alpha = 1, color = "darkgrey") +
    geom_line(data = filter(site_df, parameter == "Temperature"), aes(DT_round, mean), color = "black") +
    geom_point(data = filter(site_df, parameter == "Temperature" & grepl("site visit|sv window", final_cleaned_flag)), aes(DT_round, mean), color = "#E70870") +
    theme_bw() +
    ylim(min = min(filter(site_df, parameter == "Temperature")$mean), max = max(filter(site_df, parameter == "Temperature")$mean)) +
    ylab("Temperature C") + xlab("") + ggtitle("")

  ph <- ggplot() +
    geom_line(data = filter(contrail, Unit == "pH"), aes(Reading, Value), color = "lightblue") +
    geom_line(data = filter(downstream_df, parameter == "pH"), aes(DT_round, mean), cex = 0.2, alpha = 1, color = "lightgrey") +
    geom_line(data = filter(upstream_df, parameter == "pH"), aes(DT_round, mean), cex = 0.2, alpha = 1, color = "darkgrey") +
    geom_line(data = filter(site_df, parameter == "pH"), aes(DT_round, mean), color = "black") +
    geom_point(data = filter(site_df, parameter == "pH" & grepl("site visit|sv window", final_cleaned_flag)), aes(DT_round, mean), color = "#E70870") +
    theme_bw() +
    ylim(min = min(filter(site_df, parameter == "pH")$mean), max = max(filter(site_df, parameter == "pH")$mean)) +
    xlab("") + ylab("pH") + ggtitle("")

  orp <- ggplot() +
    geom_line(data = filter(downstream_df, parameter == "ORP"), aes(DT_round, mean), cex = 0.2, alpha = 1, color = "lightgrey") +
    geom_line(data = filter(upstream_df, parameter == "ORP"), aes(DT_round, mean), cex = 0.2, alpha = 1, color = "darkgrey") +
    geom_line(data = filter(site_df, parameter == "ORP"), aes(DT_round, mean), color = "black") +
    geom_point(data = filter(site_df, parameter == "ORP" & grepl("site visit|sv window", final_cleaned_flag)), aes(DT_round, mean), color = "#E70870") +
    theme_bw() +
    ylim(min = min(filter(site_df, parameter == "ORP")$mean), max = max(filter(site_df, parameter == "ORP")$mean)) +
    ylab("ORP") + xlab("") + ggtitle("")

  spc <- ggplot() +
    geom_line(data = filter(contrail, Unit == "uS/cm"), aes(Reading, Value), color = "lightblue") +
    geom_line(data = filter(downstream_df, parameter == "Specific Conductivity"), aes(DT_round, mean), cex = 0.2, alpha = 1, cex = 0.2, alpha = 1, color = "lightgrey") +
    geom_line(data = filter(upstream_df, parameter == "Specific Conductivity"), aes(DT_round, mean), cex = 0.2, alpha = 1, cex = 0.2, alpha = 1, color = "darkgrey") +
    geom_line(data = filter(site_df, parameter == "Specific Conductivity"), aes(DT_round, mean), color = "black", cex = 0.8) +
    geom_point(data = filter(site_df, parameter == "Specific Conductivity" & grepl("site visit|sv window", final_cleaned_flag)), aes(DT_round, mean), color = "#E70870", cex = 1.5) +
    theme_bw() +
    ylim(min = min(filter(site_df, parameter == "Specific Conductivity")$mean) - 1, max = max(filter(site_df, parameter == "Specific Conductivity")$mean) + 1) +
    ylab("SpC uS/cm") + xlab("") + ggtitle("")

  do <- ggplot() +
    geom_line(data = filter(contrail, Unit == "mg/L"), aes(Reading, Value), color = "lightblue") +
    geom_line(data = filter(downstream_df, parameter == "DO"), aes(DT_round, mean), cex = 0.2, alpha = 1, color = "lightgrey") +
    geom_line(data = filter(upstream_df, parameter == "DO"), aes(DT_round, mean), cex = 0.2, alpha = 1, color = "darkgrey") +
    geom_line(data = filter(site_df, parameter == "DO"), aes(DT_round, mean), color = "black") +
    geom_point(data = filter(site_df, parameter == "DO" & grepl("site visit|sv window", final_cleaned_flag)), aes(DT_round, mean), color = "#E70870") +
    theme_bw() +
    ylim(min = min(filter(site_df, parameter == "DO")$mean) - 1, max = max(filter(site_df, parameter == "DO")$mean) + 1) +
    ylab("DO mg/L") + xlab("") + ggtitle("")

  # Handle plotting for sondes with turbidity vs. chl-a data:

  turb <- ggplot() +
    geom_line(data = filter(contrail, Unit == "ntu"), aes(Reading, Value), color = "lightblue") +
    geom_line(data = filter(downstream_df, parameter == "Turbidity"), aes(DT_round, mean), cex = 0.2, alpha = 1, color = "lightgrey") +
    geom_line(data = filter(upstream_df, parameter == "Turbidity"), aes(DT_round, mean), cex = 0.2, alpha = 1, color = "darkgrey") +
    geom_line(data = filter(site_df, parameter == "Turbidity"), aes(DT_round, mean), color = "black") +
    geom_point(data = filter(site_df, parameter == "Turbidity" & grepl("site visit|sv window", final_cleaned_flag)), aes(DT_round, mean), color = "#E70870") +
    theme_bw() +
    ylim(min = min(filter(site_df, parameter == "Turbidity")$mean) + 10, max = max(filter(site_df, parameter == "Turbidity")$mean) + 10) +
    ylab("Turbidity NTU") + xlab("") + ggtitle("")

  if(nrow(filter(site_df, parameter == "Chl-a Fluorescence")) > 0){
    chla <- ggplot() +
      geom_line(data = filter(downstream_df, parameter == "Chl-a Fluorescence"), aes(DT_round, mean), cex = 0.2, alpha = 1, color = "lightgrey") +
      geom_line(data = filter(upstream_df, parameter == "Chl-a Fluorescence"), aes(DT_round, mean), cex = 0.2, alpha = 1, color = "darkgrey") +
      geom_line(data = filter(site_df, parameter == "Chl-a Fluorescence"), aes(DT_round, mean), color = "black") +
      geom_point(data = filter(site_df, parameter == "Chl-a Fluorescence" & grepl("site visit|sv window", final_cleaned_flag)), aes(DT_round, mean), color = "#E70870") +
      theme_bw() +
      ylab("Chl-a RFU") + xlab("") + ggtitle("")
  }

  if(nrow(filter(site_df, parameter == "Chl-a Fluorescence")) == 0){ chla <- NULL }

  # rain proxy plot. this just plots the sum of at an hourly timestep of
  # rain at all sites in the watershed. not rooted in anything scientific,
  # but just shows in a crude way whether or not there was rain happening
  # at any given time (and at any given point) during the time of interest.
  rain <- ggplot(data = rain_data) +
    geom_col(aes(x = DT,
                 y = total_ws), color = "#002EA3") +
    ylab("Rain (PROXY)") + xlab("") +
    ylim(min = 0, max = rain_stat$total_ws) +
    theme_bw()

  # handling what gets plotted based on user input:

  if(parameter %in% c("All", "all", "all parameters", "all params")){
    if(is.null(chla)) {

      plot <- ggarrange(depth, temp, ph, orp, spc, do, turb,
                         rain, nrow = 4, ncol = 2)
    }

    if(!is.null(chla)){
      plot <- ggarrange(depth, temp, ph, orp, spc, do, chla,
                         rain, nrow = 4, ncol = 2)
    }
  }

  if(parameter %in% c("pH", "ph")){
    plot <- ggarrange(depth, temp, ph, rain,
                      ncol = 1, nrow = 4)
  }

  if(parameter %in% c("orp", "ORP")){
    plot <- ggarrange(depth, temp, orp, rain,
                       ncol = 1, nrow = 4)
  }

  if(parameter %in% c("do","DO")){
    plot <- ggarrange(depth, temp, do, rain,
                      ncol = 1, nrow = 4)
  }

  if(parameter %in% c("turbidity", "turb", "Turb", "Turbidity")){
    plot <- ggarrange(depth, temp, turb, rain,
                      ncol = 1, nrow = 4)
  }

  if(parameter %in% c("Specific Conductivity", "conductivity", "specific conductivity", "sc", "spc", "SC","SpC")){
    plot <- ggarrange(depth, temp, spc, rain,
                      ncol = 1, nrow = 4)
  }

  if(parameter %in% c("chla", "Chla", "chl-a", "chlorophyll a", "Chl-a Fluorescence")){
    plot <- ggarrange(depth, temp, chla, rain,
                      ncol = 1, nrow = 4)
  }

  return(plot)
}

retrieve_relevant_data_name <- function(df_name_arg, year_week_arg = NULL, interval_arg = NULL) {
  if(is.null(interval_arg)){
    if (df_name_arg %in% names(verified_data) & any(year_week_arg %in% verified_data[[df_name_arg]]$y_w)) {
      return("verified_data")
    }
    if (df_name_arg %in% names(intermediary_data) & any(year_week_arg %in% intermediary_data[[df_name_arg]]$y_w)) {
      return("intermediary_data")
    }
    if (df_name_arg %in% names(all_data) & any(year_week_arg %in% all_data[[df_name_arg]]$y_w)) {
      return("all_data")
    }
  }
  if(is.null(year_week_arg)){
    if (df_name_arg %in% names(verified_data) & any(verified_data[[df_name_arg]]$DT_round %within% interval_arg)) {
      return("verified_data")
    }
    if (df_name_arg %in% names(intermediary_data) & any(intermediary_data[[df_name_arg]]$DT_round %within% interval_arg)) {
      return("intermediary_data")
    }
    if (df_name_arg %in% names(all_data) & any(all_data[[df_name_arg]]$DT_round %within% interval_arg)) {
      return("all_data")
    }
  }
}

## Water Sampling Data:

# Goal:
#Save data in the correct format for RMRS spreadsheet
#Save all water sampling probe values in a spreadsheet
# To get the RMRS style data for a specfic date of sampling,
# Input the date of interest in sampling_spreadsheet_creator
#sampling_spreadsheet_creator(date_oi = "2023-11-17")

# To get all the water sampling data and save to CSV in sampling notes
# This also returns the df sampling_notes in case you want to review in R
#sampling_spreadsheet_creator(all_dates = TRUE)


sampling_spreadsheet_creator <- function(date_oi = "2023-10-16", all_dates = FALSE ){

  #source clean mwater script for all notes cleaned

  source("src/load_mWater_notes.R")

  #pull in site meta data
  site_meta <- read_csv("data/metadata/water_sampling_sites.csv",show_col_types = FALSE)%>%
    select(site = site_code, Site_Name, site_label_rmrs)
  # sort for sites in upper network (ie. acronyms rather than street names)
  upper_sites <- read_csv("data/metadata/water_sampling_sites.csv",show_col_types = FALSE)%>%
    filter(watershed != "CLP  Mainstem-Fort Collins")%>%
    #this is to help match with user input
    mutate(site_code = tolower(site_code))

  # create df of all water samples and save DT, handheld probe and chla volume data
  sampling_notes <- load_mWater_notes()%>%
    filter(grepl("Sampling",visit_type))%>%
    mutate(all_pics_taken = case_when(!is.na(downstream_pic)&!is.na(upstream_pic)&!is.na(clarity)&!is.na(filter_pic) ~ TRUE, TRUE ~ FALSE),
           #correct names if it is in our upper sites (acronyms)
           site = ifelse(site %in% upper_sites$site_code, toupper(site), site),
           DT_round = round_date(start_DT, "15 minutes"))%>%
    select(site,crew, DT_round, date, time = start_time_mst, sample_collected, chla_volume_ml, vol_filtered_blank_dup, do_mgl, cond_ms_cm, temp_c, visit_comments, all_pics_taken, q_cfs)

  # Distinguish BLANK and DUPLICATE values
  blanks_dups <- sampling_notes %>%
    #find all values that have blank or dup
    filter(grepl("DUPLICATE|BLANK", sample_collected)) %>%
    # change sample collected to match BLANK/DUP
    mutate(sample_collected = ifelse(grepl("DUPLICATE", sample_collected), "DUPLICATE", "BLANK"),
           # Volume filtered blank dup becomes chla volume
           chla_volume_ml = vol_filtered_blank_dup,
           #drop vol_filtered_blank/dup
           vol_filtered_blank_dup = NULL)

  # Add blank and duplicate values back to main
  sampling_notes <- sampling_notes%>%
    #get rid of blank/dup in sample collcected
    mutate(sample_collected = gsub("DUPLICATE|BLANK| |,", "", sample_collected),
           #drop vol_filtered_blank/dup
           vol_filtered_blank_dup = NULL)%>%
    #bring in blank and dup rows
    rbind(blanks_dups)%>%
    #arrange by datetime and site (Blanks and dups go second)
    arrange(DT_round, site)%>%
    # join with RMRS friendly metadata
    left_join(site_meta, by = "site")



 if (all_dates == TRUE) {

   sampling_notes_output <- sampling_notes%>%
     # select only the needed columns, saved in the correct order and fix column names
     select(site_code = site, Date = date, SampleType = sample_collected, time_mst = time,chla_volume_ml,  do_mgl, cond_ms_cm, temp_c, visit_comments)

   # write to csv
   write_csv(x = sampling_notes_output, file = paste0("data/sampling_notes/all_samples_as_of_",as.character(Sys.Date()),".csv" ))
 }else{
   #grab desired date
   date_oi_clean <- as.Date(date_oi, tz = "America/Denver")
   # filter sampling notes df by desired date
   samples_of_day <- filter(sampling_notes,date == date_oi_clean )%>%
     #match date to RMRS style
     mutate(Date = format(date, "%d-%b-%y"),
            #create SiteDescr column for RMRS sheet
            SiteDescr = paste0(site, "_",format(date, "%m%d%y")))%>%
     # select only the needed columns, saved in the correct order and fix column names
     select(Site = site ,SiteDescr, SiteLabel = site_label_rmrs , Date, SampleType = sample_collected, q_cfs,
            time, do_mgl, cond_ms_cm, temp_c, notes = visit_comments  )
   # write to csv
   write_csv(x = samples_of_day, file = paste0("data/sampling_notes/samples_of_",date_oi,".csv" ))
 }
}


save_intermediary_dir <- function(updated_df, df_name) {
  # if updated_site_param_df is from all_data, but there is already intermediary data for it, refuse to move on
  if (all(!updated_df$is_verified) & df_name %in% list.files(intermediary_path)) {
    while (TRUE) {
      cat("!!!YOU ARE ABOUT TO OVERWRITE INTERMEDIARY DATA WITH PRE-VERIFICATION DATA!!!")
      user_input <- readline(prompt = "Are you sure you want to continue with this action? (yes/no):")
      user_input <- tolower(user_input)

      if (user_input == "yes") {
        saveRDS(updated_df, paste0(intermediary_path, df_name))
        break
      }

      if (user_input == "no") {
        working_data <<- set_names(map(list.files(intermediary_path, full.names = TRUE), readRDS), list.files(intermediary_path))
        updated_site_param_df <<- working_data[[site_param_name]] # working dir df
        saveRDS(updated_site_param_df, paste0(intermediary_path, df_name))
        break
      }

      cat("Invalid input. Please enter one of the options from the following:\n\n")
      cat("yes\n")
      cat("no\n")
    }
  }

  # if updated_site_param_df is from verified data, refuse to move on
  if (df_name %in% list.files(verified_path)) { # test to see if this works
    while (TRUE) {
      cat("!!!YOU ARE ABOUT TO WORK ON DATA THAT HAS BEEN FINALIZED!!!")
      user_input <- readline(prompt = "Are you sure you want to continue with this action? (yes/no):")
      user_input <- tolower(user_input)

      if (user_input == "yes") {
        saveRDS(updated_df, paste0(intermediary_path, df_name))
        break
      }

      if (user_input == "no") {
        cat("Please update your `site` and `parameter` objects to data that has not been finalized.")
        break
      }

      cat("Invalid input. Please enter one of the options from the following:\n\n")
      cat("yes\n")
      cat("no\n")
    }
  }

  saveRDS(updated_df, paste0(intermediary_path, df_name))
}

# Generate plots with both weekly and daily flagged data.
# This function will generate a list of plots with both weekly and daily flagged data.
# @param site_arg A string of the site name.
# @param parameter_arg A string of the parameter name.
# @param flag_arg A string of the flag name.
# @return A list of plots with both weekly and daily flagged data.
# @examples
# stack_flag_plots(site_arg = "archery", parameter_arg = "Actual Conductivity", flag_arg = "outside of Actual Conductivity sensor specification range")
# stack_flag_plots(site_arg = "boxelder", parameter_arg = "Temperature", flag_arg = "outside of Temperature sensor specification range")

# stack_flag_plots <- function(site_arg, parameter_arg, flag_arg, df_list) {
#   # Call on the weekly and daily functions and fill their args with this
#   # functions args
#   weekly_plot_list <- generate_weekly_flag_plots(site_arg = site_arg, parameter_arg = parameter_arg, flag_arg = flag_arg, df_list = df_list)
#   daily_plot_list <- generate_daily_flag_plots(site_arg = site_arg, parameter_arg = parameter_arg, flag_arg = flag_arg, df_list = df_list)
#
#   # These two functions should always return the same amount of plots, so we can
#   # use map2() to stack them with combine_plots()
#   weekly_daily_plots <- map2(.x = weekly_plot_list, .y = daily_plot_list, ~ggarrange(.x, .y, nrow = 2, ncol = 1, heights = 5, widths = 12))
#   return(weekly_daily_plots)
# }

#' @title Summarize site parameter data from the API and field notes data frames.
#' @description
#' A short description...
#' @param site_arg A site name.
#' @param parameter_arg A parameter name.
#' @param api_data A dataframe with the munged API data.
#' @param notes The munged field notes
#' @return A dataframe with summary statistics for a given site parameter data frame.
#' @examples
# summarize_site_param(site_arg = "archery", parameter_arg = "Actual Conductivity", api_data = incoming_data_collated_csvs)
# summarize_site_param(site_arg = "boxelder", parameter_arg = "Temperature", api_data = incoming_data_collated_csvs)

summarize_site_param_full <- function(site_arg, parameter_arg, api_data, notes = field_notes) {

  # filter deployment records for the full join
  site_field_notes <- notes %>%
    filter(grepl(paste(unlist(str_split(site_arg, " ")), collapse = "|"), site, ignore.case = TRUE))

  # filtering the data and generating results
  summary <- tryCatch({
    api_data %>%
      # subset to single site-parameter combo:
      dplyr::filter(site == site_arg & parameter == parameter_arg) %>%
      # safety step of removing any erroneous dupes
      dplyr::distinct() %>%
      # across each 15 timestep, get the average value, spread, and count of obs
      dplyr::group_by(DT_round, site, parameter) %>%
      dplyr::summarize(mean = as.numeric(mean(value, na.rm = T)),
                       diff = abs(min(value, na.rm = T) - max(value, na.rm = T)),
                       n_obs = n()) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(DT_round) %>%
      # pad the dataset so that all 15-min time stamps are present
      padr::pad(by = "DT_round", interval = "15 min") %>%
      # add a DT_join column to join field notes to (make DT_round character string, so no
      # funky DT issues occur during the join):
      dplyr::mutate(DT_join = as.character(DT_round),
                    site = site_arg,
                    parameter = parameter_arg,
                    flag = NA) %>% # add "flag" column for future processing
      # join our tidied data frame with our field notes data:
      dplyr::left_join(dplyr::filter(dplyr::select(site_field_notes, sonde_employed, last_site_visit, DT_join, visit_comments, sensor_malfunction, cals_performed)),
                       by = c('DT_join')) %>%
      # make sure DT_join is still correct:
      dplyr::mutate(DT_round = lubridate::as_datetime(DT_join, tz = "MST")) %>%
      # Use fill() to determine when sonde was in the field, and when the last site visit was
      # Necessary step for FULL dataset only (this step occurs in combine_hist_inc_data.R for auto)
      tidyr::fill(c(sonde_employed, last_site_visit, sensor_malfunction)) %>%
      dplyr::mutate(sonde_employed = ifelse(is.na(sonde_employed), 0, sonde_employed))
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

#' @title Summarize site parameter data from the HydroVu API and field notes data frames.
#'
#' @description
#' A function that summarizes site parameter data from the API and field notes data frames.
#'
#' @param site_arg A site name.
#'
#' @param parameter_arg A parameter name.
#'
#' @param api_data A dataframe with the munged API data.
#'
#' @param notes The munged field notes
#'
#' @param require An upstream dependency (i.e., alternative to `{necessary object}`).
#'
#' @return A dataframe with summary statistics and field notes for a given site parameter data frame.
#'
#' @examples
# summarize_site_param(site_arg = "archery", parameter_arg = "Actual Conductivity", api_data = incoming_data_collated_csvs)
# summarize_site_param(site_arg = "boxelder", parameter_arg = "Temperature", api_data = incoming_data_collated_csvs)

summarize_site_param <- function(site_arg, parameter_arg, api_data, notes, require = NULL) {

  # filter deployment records for the full join
  site_field_notes <- notes %>%
    dplyr::filter(site == site_arg)

  # filtering the data and generating results
  summary <- tryCatch({
    api_data %>%
      # summarize the data
      dplyr::filter(site == site_arg & parameter == parameter_arg) %>%
      dplyr::distinct() %>%
      dplyr::group_by(DT_round) %>%
      dplyr::summarize(mean = as.numeric(mean(value, na.rm = T)),
                       diff = abs(min(value, na.rm = T) - max(value, na.rm = T)),
                       n_obs = n()) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(DT_round) %>%
      # pad the dataset so that all 15-min timestamps are present
      padr::pad(by = "DT_round", interval = "15 min") %>%
      # join the field notes
      dplyr::mutate(DT_join = as.character(DT_round),
                    site = site_arg,
                    parameter = parameter_arg,
                    flag = NA) %>%
      dplyr::left_join(dplyr::filter(dplyr::select(site_field_notes, sonde_employed,
                                                   last_site_visit, DT_join, site,
                                                   visit_comments, sensor_malfunction,
                                                   cals_performed)),
                       by = c('DT_join', 'site')) %>%
      # join the metaparameter dfs
      dplyr::left_join(site_metaparams_list[[site_arg]], by = "DT_join") %>%
      dplyr::mutate(DT_round = lubridate::as_datetime(DT_join, tz = "MST"))
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

#' @title Update Historical Flag List
#'
#' @description
#' A function that updates the historical flag list with new flagged data.
#'
#' @param new_flagged_data A list of data frames that have been flagged.
#'
#' @param historical_flagged_data A list of data frames that have been flagged.
#'
#' @return A list of data frames that have been updated with new flagged data.
#'
#' @examples
#' update_historical_flag_list(new_flagged_data = all_data_flagged$`archery-Actual Conductivity`, historical_flagged_data = all_data_flagged$`archery-Actual Conductivity`)
#' update_historical_flag_list(new_flagged_data = all_data_flagged$`boxelder-Temperature`, historical_flagged_data = all_data_flagged$`boxelder-Temperature`)

update_historical_flag_list <- function(new_flagged_data, historical_flagged_data){

  # Get the matching index names
  matching_indexes <- intersect(names(new_flagged_data), names(historical_flagged_data))

  # bind new_flagged_data and historical_flagged_data together
  updated_historical_flag_list <- map(matching_indexes, function(index) {

    old <- historical_flagged_data[[index]] %>%
      filter(DT_round < ymd_hms(max(DT_round) - hours(24), tz = "MST")) %>% # this is the antijoin step (but i thought we were doing 24 hours?) -jd changing this to 24 removes the duplicate problem
      mutate(last_site_visit = force_tz(last_site_visit, tzone = "MST"))

    bind_rows(old, new_flagged_data[[index]]) %>%
      arrange(DT_round) %>%
      select(-historical) # if we are removing this then I think we should remove that step from combined_data -jd
  }) %>%
    set_names(matching_indexes) %>%
    keep(~ !is.null(.))

  # need to make sure that when we are doing this we are not getting rid of dfs in the RDS list, still needs to be checked -jd
    # if we are matching up the historical list(88) with the incoming list(72) and only returning the matches then we will miss
    # be returning a list that is being written that is shorter than the list that we started with... Definitely needs to be fixed.

  return(updated_historical_flag_list)

}

# error 1: tz discrepancy at the join between incoming data and HFD
  # this messes up the future pulls because the start times df will be wrong

# time_diffs <- diff(test_data$DT_round)
#
# time_diffs # these should all be 15 and there is one that is not (where we joined incoming data to the HFD)

verify_flag_data <- function(df_list_arg,
                             site_arg,
                             parameter_arg, # do we want this to be able to be null eventually?
                             flag_arg = NULL,
                             weekly_plot_object) {

  # set layer
  layer <- "weekly"

  #quit

  if (!QUIT) {
    # verification decision lists ----
    simple_decision_list <- c("PASS ALL", "FAIL NONE",
                              "PASS VALID", "FAIL FLAGGED",
                              "PASS FLAGGED", "FAIL VALID",
                              "PASS NONE", "FAIL ALL",
                              "SKIP", "QUIT")

    inspect_decision_list <- c("INSPECT ALL", "INSPECT VALID",
                               "INSPECT FLAGGED", "INSPECT SOME")

    # initialize site_param_df ----
    site_param_df <- df_list_arg[[paste0(site_arg, "-", parameter_arg)]]

    # get data from initial weekly plot object ----
    df_weekly_data <- ggplot_build(weekly_plot_object)$plot$data %>%
      filter(site == site_arg) %>%
      select(intersect(names(.), names(site_param_df)))

    df_weekly_name <- ggplot_build(weekly_plot_object)$plot$labels$title

    # initialize list to contain updated dfs OUTPUT, VERY IMPORTANT ----
    altered_df_list <- list()

    # initialize week list ----
    weekday_list <- c(1,2,3,4,5,6,7)

    # Make sure that non-removable flags are still in place ----
    non_removable_flags <- paste0(c("sensor biofouling","sonde burial",
                                    "sensor malfunction","frozen","missing data",
                                    "sonde not employed","sonde unsubmerged",
                                    "depth calibration malfunction"),
                                  collapse = "|")

    # Which data are you working on: ----
    cat("Currently verifying: ", df_weekly_name, "data. \n")

    print(weekly_plot_object)

    # Get weekly decision(s) ----
    verification_decision <- get_verification_decision(
      paste("Would you like to (pass/fail/skip/inspect/quit) ", df_weekly_name, "?\n(layer~", layer, ")(px/fx/sx/ix/q): ") # this should be weekly data
    )

    if (verification_decision %in% simple_decision_list){
      altered_df_list <- alter_verification_column(verification_decision, df_weekly_data)
    }

    if (verification_decision %in% inspect_decision_list){
      # set layer
      layer <- "daily"

      # read in the data relevant to the current site-parameter that we are doing
      # this will be used to make the inspect decision plots. pull this data from the working directory
      # site_param_df <<- df_list[[paste0(site_arg, "-", parameter_arg)]] # if parameter arg is null then pull it all in

      # inspect all weekly ----
      if (verification_decision == "INSPECT ALL"){

        for (i in unique(df_weekly_data$weekday)) {

          daily_plot_data <- site_param_df %>%
            filter((y_w %in% df_weekly_data$y_w) & (weekday == i))

          daily_plot_object <- generate_daily_plot(plot_data_arg = daily_plot_data,
                                                   df_list_arg = df_list_arg,
                                                   site_arg = site_arg,
                                                   parameter_arg = parameter_arg)

          # weekly_plot_object
          weekly_plot_object <- generate_supplemental_weekly_plot(daily_plot_data_arg = daily_plot_data,
                                                                  df_list_arg = df_list_arg,
                                                                  site_arg = site_arg,
                                                                  parameter_arg = parameter_arg) # need to make sure that these functions can read in information from outside of them

          # print ggarrange daily and weekly plots
          print(ggarrange(daily_plot_object, weekly_plot_object, nrow = 2, ncol = 1))

          # Generate the prompt for the user to answer
          day_dt <- as.character(head(daily_plot_data$DT_round, 1)) # reformat so that this reads as a date
          day_choice_prompt <- paste("Would you like to (pass/fail/skip/inspect) ", day_dt, "?\n(layer~", layer, ")(px/fx/sx/ix): ")
          daily_verification_decision <- get_verification_decision(day_choice_prompt)

          while (layer == "daily" & daily_verification_decision == "QUIT") {
            cat("You cannot quit during this step. Please choose a different option.\n")
            daily_verification_decision <- get_verification_decision(day_choice_prompt)
          }

          if (daily_verification_decision %in% simple_decision_list){
            altered_df_list[[i]] <- alter_verification_column(daily_verification_decision, daily_plot_data) # this should get put in a list
          }

          if (daily_verification_decision %in% inspect_decision_list){
            altered_df_list[[i]] <- verify_subdaily(daily_verification_decision_arg = daily_verification_decision,
                                                    daily_plot_data_arg = daily_plot_data,
                                                    daily_plot_object_arg = daily_plot_object,
                                                    weekly_plot_object_arg = weekly_plot_object
            )
          }
        }
      }

      # inspect valid weekly ----
      if (verification_decision == "INSPECT VALID") {

        # get which days to inspect
        days_to_inspect <- unique(df_weekly_data %>%
                                    filter(is.na(flag)) %>%
                                    pull(weekday))

        new_weekday_list <- weekday_list[weekday_list >= min(df_weekly_data$weekday) & weekday_list <= max(df_weekly_data$weekday)]
        non_days_to_inspect <- new_weekday_list[!new_weekday_list %in% days_to_inspect]
        # for each non-inspect day decision:
        if (length(non_days_to_inspect) > 0) {

          # get non-inspect day decision:
          noninspect_prompt_text = paste("Would you like to (pass/fail/skip) non-inspect days?\n(layer~", layer, ")(px/fx/sx): ")
          non_inspect_verification_decision <- get_verification_decision(prompt_text = noninspect_prompt_text)

          # Handle the case when the user selects an "INSPECT" or "QUIT" option for non-inspect days
          while (layer == "daily" & non_inspect_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME", "QUIT")) {
            if(non_inspect_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME")) {
              cat("You cannot inspect non-inspect days. Please choose a different option.\n")
              non_inspect_verification_decision <- get_verification_decision(prompt_text = noninspect_prompt_text)
            }
            if(non_inspect_verification_decision == "QUIT") {
              cat("You cannot quit during this step. Please choose a different option.\n")
              non_inspect_verification_decision <- get_verification_decision(prompt_text = noninspect_prompt_text)
            }
          }

          for (i in non_days_to_inspect) {
            daily_plot_data <- site_param_df %>%
              filter((y_w %in% df_weekly_data$y_w) & (weekday == i))

            altered_df_list[[i]] <- alter_verification_column(non_inspect_verification_decision, daily_plot_data)
          }
        }

        for (i in days_to_inspect) {
          daily_plot_data <- site_param_df %>%
            filter((y_w %in% df_weekly_data$y_w) & (weekday == i))

          daily_plot_object <- generate_daily_plot(plot_data_arg = daily_plot_data,
                                                   df_list_arg = df_list_arg,
                                                   site_arg = site_arg,
                                                   parameter_arg = parameter_arg)

          # weekly_plot_object
          weekly_plot_object <- generate_supplemental_weekly_plot(daily_plot_data_arg = daily_plot_data,
                                                                  df_list_arg = df_list_arg,
                                                                  site_arg = site_arg,
                                                                  parameter_arg = parameter_arg) # need to make sure that these functions can read in information from outside of them

          # print ggarrange daily and weekly plots
          print(ggarrange(daily_plot_object, weekly_plot_object, nrow = 2, ncol = 1))

          # Generate the prompt for the user to answer
          day_dt <- as.character(head(daily_plot_data$DT_round, 1)) # reformat so that this reads as a date
          day_choice_prompt <- paste("Would you like to (pass/fail/skip/inspect) ", day_dt, "?\n(layer~", layer, ")(px/fx/sx/ix): ")
          daily_verification_decision <- get_verification_decision(day_choice_prompt)

          while (daily_verification_decision == "QUIT") {
            cat("You cannot quit during this step. Please choose a different option.\n")
            daily_verification_decision <- get_verification_decision(day_choice_prompt)
          }

          if (daily_verification_decision %in% simple_decision_list){
            altered_df_list[[i]] <- alter_verification_column(daily_verification_decision, daily_plot_data) # this should get put in a list
          }

          if (daily_verification_decision %in% inspect_decision_list){
            altered_df_list[[i]] <- verify_subdaily(daily_verification_decision_arg = daily_verification_decision,
                                                    daily_plot_data_arg = daily_plot_data,
                                                    daily_plot_object_arg = daily_plot_object,
                                                    weekly_plot_object_arg = weekly_plot_object
            )
          }
        }
      }

      # inspect flagged weekly ----
      if (verification_decision == "INSPECT FLAGGED") {

        # get which days to inspect
        days_to_inspect <- unique(df_weekly_data %>%
                                    filter(!is.na(flag)) %>%
                                    pull(weekday))

        new_weekday_list <- weekday_list[weekday_list >= min(df_weekly_data$weekday) & weekday_list <= max(df_weekly_data$weekday)]
        non_days_to_inspect <- new_weekday_list[!new_weekday_list %in% days_to_inspect]

        # for each non-inspect day decision:
        if (length(non_days_to_inspect) > 0) {

          # get non-inspect day decision:
          noninspect_prompt_text = paste("Would you like to (pass/fail/skip) non-inspect days?\n(layer~", layer, ")(px/fx/sx): ")
          non_inspect_verification_decision <- get_verification_decision(prompt_text = noninspect_prompt_text)

          # Handle the case when the user selects an "INSPECT" or "QUIT" option for non-inspect days
          while (layer == "daily" & non_inspect_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME", "QUIT")) {
            if(non_inspect_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME")) {
              cat("You cannot inspect non-inspect days. Please choose a different option.\n")
              non_inspect_verification_decision <- get_verification_decision(prompt_text = noninspect_prompt_text)
            }
            if(non_inspect_verification_decision == "QUIT") {
              cat("You cannot quit during this step. Please choose a different option.\n")
              non_inspect_verification_decision <- get_verification_decision(prompt_text = noninspect_prompt_text)
            }
          }

          for (i in non_days_to_inspect) {
            daily_plot_data <- site_param_df %>%
              filter((y_w %in% df_weekly_data$y_w) & (weekday == i))

            altered_df_list[[i]] <- alter_verification_column(non_inspect_verification_decision, daily_plot_data)
          }
        }

        for (i in days_to_inspect) {
          daily_plot_data <- site_param_df %>%
            filter((y_w %in% df_weekly_data$y_w) & (weekday == i))

          daily_plot_object <- generate_daily_plot(plot_data_arg = daily_plot_data,
                                                   df_list_arg = df_list_arg,
                                                   site_arg = site_arg,
                                                   parameter_arg = parameter_arg)


          # weekly_plot_object
          weekly_plot_object <- generate_supplemental_weekly_plot(daily_plot_data_arg = daily_plot_data,
                                                                  df_list_arg = df_list_arg,
                                                                  site_arg = site_arg,
                                                                  parameter_arg = parameter_arg) # need to make sure that these functions can read in information from outside of them

          # print ggarrange daily and weekly plots
          print(ggarrange(daily_plot_object, weekly_plot_object, nrow = 2, ncol = 1))

          # Generate the prompt for the user to answer
          day_dt <- as.character(head(daily_plot_data$DT_round, 1)) # reformat so that this reads as a date
          day_choice_prompt <- paste("Would you like to (pass/fail/skip/inspect) ", day_dt, "?\n(layer~", layer, ")(px/fx/sx/ix): ")
          daily_verification_decision <- get_verification_decision(day_choice_prompt)

          while (layer == "daily" & daily_verification_decision == "QUIT") {
            cat("You cannot quit during this step. Please choose a different option.\n")
            daily_verification_decision <- get_verification_decision(day_choice_prompt)
          }

          if (daily_verification_decision %in% simple_decision_list){
            altered_df_list[[i]] <- alter_verification_column(daily_verification_decision, daily_plot_data) # this should get put in a list
          }

          if (daily_verification_decision %in% inspect_decision_list){
            altered_df_list[[i]] <- verify_subdaily(daily_verification_decision_arg = daily_verification_decision,
                                                    daily_plot_data_arg = daily_plot_data,
                                                    daily_plot_object_arg = daily_plot_object,
                                                    weekly_plot_object_arg = weekly_plot_object
            )
          }
        }
      }

      # inspect some weekly ----
      if (verification_decision == "INSPECT SOME") {

        # get which days to inspect
        days_to_inspect <- get_weekly_inspection_decision(weekly_plot_data = df_weekly_data) # rename to daily ***

        new_weekday_list <- weekday_list[weekday_list >= min(df_weekly_data$weekday) & weekday_list <= max(df_weekly_data$weekday)]
        non_days_to_inspect <- new_weekday_list[!new_weekday_list %in% days_to_inspect]

        # for each non-inspect day decision:
        if (length(non_days_to_inspect)>0) {

          # get non-inspect day decision:
          noninspect_prompt_text = paste("Would you like to (pass/fail/skip) non-inspect days?\n(layer~", layer, ")(px/fx/sx): ")
          non_inspect_verification_decision <- get_verification_decision(prompt_text = noninspect_prompt_text)

          # Handle the case when the user selects an "INSPECT" or "QUIT" option for non-inspect days
          while (layer == "daily" & non_inspect_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME", "QUIT")) {
            if(non_inspect_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME")) {
              cat("You cannot inspect non-inspect days. Please choose a different option.\n")
              non_inspect_verification_decision <- get_verification_decision(prompt_text = noninspect_prompt_text)
            }
            if(non_inspect_verification_decision == "QUIT") {
              cat("You cannot quit during this step. Please choose a different option.\n")
              non_inspect_verification_decision <- get_verification_decision(prompt_text = noninspect_prompt_text)
            }
          }

          for (i in non_days_to_inspect) {
            daily_plot_data <- site_param_df %>%
              filter((y_w %in% df_weekly_data$y_w) & (weekday == i))

            altered_df_list[[i]] <- alter_verification_column(non_inspect_verification_decision, daily_plot_data)
          }
        }

        for (i in days_to_inspect) {

          daily_plot_data <- site_param_df %>%
            filter((y_w %in% df_weekly_data$y_w) & (weekday == i))

          daily_plot_object <- generate_daily_plot(plot_data_arg = daily_plot_data,
                                                   df_list_arg = df_list_arg,
                                                   site_arg = site_arg,
                                                   parameter_arg = parameter_arg)

          # weekly_plot_object

          weekly_plot_object <- generate_supplemental_weekly_plot(daily_plot_data_arg = daily_plot_data,
                                                                  df_list_arg = df_list_arg,
                                                                  site_arg = site_arg,
                                                                  parameter_arg = parameter_arg)

          # print ggarrange daily and weekly plots
          print(ggarrange(daily_plot_object, weekly_plot_object, nrow = 2, ncol = 1))

          # Generate the prompt for the user to answer
          day_dt <- as.character(head(daily_plot_data$DT_round, 1)) # reformat so that this reads as a date
          day_choice_prompt <- paste("Would you like to (pass/fail/skip/inspect) ", day_dt, "?\n(layer~", layer, ")(px/fx/sx/ix): ")
          daily_verification_decision <- get_verification_decision(day_choice_prompt)

          while (layer == "daily" & daily_verification_decision == "QUIT") {
            cat("You cannot quit during this step. Please choose a different option.\n")
            daily_verification_decision <- get_verification_decision(day_choice_prompt)
          }

          if (daily_verification_decision %in% simple_decision_list){
            # browser()
            altered_df_list[[i]] <- alter_verification_column(daily_verification_decision, daily_plot_data) # this should get put in a list
          }

          if (daily_verification_decision %in% inspect_decision_list){
            altered_df_list[[i]] <- verify_subdaily(daily_verification_decision_arg = daily_verification_decision,
                                                    daily_plot_data_arg = daily_plot_data,
                                                    daily_plot_object_arg = daily_plot_object,
                                                    weekly_plot_object_arg = weekly_plot_object
            )
          }
        }
      }

    }


    if (is_tibble(altered_df_list)) {
      altered_df_list <- altered_df_list %>%
        mutate(mean_verified = if_else(grepl(non_removable_flags, flag, ignore.case = TRUE), NA, mean_verified),
               is_verified = TRUE,
               verification_status = ifelse(grepl(non_removable_flags, flag, ignore.case = TRUE), "PASS", verification_status))

    } else if (is.list(altered_df_list)) {
      if (!length(altered_df_list) %in% c(51,49,46)){
        # browser()
        altered_df_list <- map(altered_df_list, function(df) {
          if(!is.null(df)) {
            # browser()
            df %>%
              mutate(mean_verified = if_else(grepl(non_removable_flags, flag, ignore.case = TRUE), NA, mean_verified),
                     is_verified = TRUE,
                     verification_status = ifelse(grepl(non_removable_flags, flag, ignore.case = TRUE), "PASS", verification_status))
          }
        })
      } else {
        # browser()
        altered_df_list <- altered_df_list %>%
          mutate(mean_verified = if_else(grepl(non_removable_flags, flag, ignore.case = TRUE), NA, mean_verified),
                 is_verified = TRUE,
                 verification_status = ifelse(grepl(non_removable_flags, flag, ignore.case = TRUE), "PASS", verification_status))
      }

    }

    return(altered_df_list) # this should not be an altered df, but rather a list of dfs that have been altered.
    # this means that the output of this method will be a list that contains a list of lists

    # If decision in quit (this can only happen before they start verifying a week of data)
    # save the data into the intermediary directory
  } else if (QUIT) {
    df_weekly_name <- ggplot_build(weekly_plot_object)$plot$labels$title
    cat("quitting out of", df_weekly_name, "\n")
    return(NULL)
  }

}


verify_subdaily <- function(daily_verification_decision_arg,
                            daily_plot_data_arg,
                            daily_plot_object_arg,
                            weekly_plot_object_arg){

  # set layer
  layer <- "subdaily"

  # inspect all
  if (daily_verification_decision_arg == "INSPECT ALL") {

    for (j in 1:nrow(daily_plot_data_arg)) {

      # create a data frame with a single row containing the desired DT_round value
      # and add it to the week plot and day plot
      flag_demarcation_data <- tibble(DT_round = daily_plot_data_arg$DT_round[j])

      # ignore flags that you cannot alter
      daily_plot_object_new <- daily_plot_object_arg +
        geom_vline(data = flag_demarcation_data,
                   aes(xintercept = DT_round,
                       color = "data point of interest",
                       # size = 2, # fix this later
                       alpha = 0.05))

      weekly_plot_object_new <- weekly_plot_object_arg +
        geom_vline(data = flag_demarcation_data,
                   aes(xintercept = DT_round,
                       color = "data point of interest",
                       # size = 2, # fix this later
                       alpha = 0.05))

      # print ggarrange daily and weekly plots
      print(ggarrange(daily_plot_object_new, weekly_plot_object_new, nrow = 2, ncol = 1))

      # generate the prompt for the user to answer
      # the user should not be able to inspect again
      day_dt <- daily_plot_data_arg$DT_join[j] # reformat so that this reads as a date
      dt_choice_prompt <- paste("Would you like to (pass/fail/skip) ", day_dt, "?\n(layer~", layer, ")(px/fx/sx)")
      dt_verification_decision <- get_verification_decision(dt_choice_prompt)

      # Handle the case when the user selects an "INSPECT" or "QUIT" option inappropriately
      while (layer == "subdaily" & dt_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME", "QUIT")) {
        if (dt_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME")){
          cat("You cannot inspect at this level. Please choose a different option.\n")
          dt_verification_decision <- get_verification_decision(dt_choice_prompt)
        }
        if (dt_verification_decision == "QUIT"){
          cat("You cannot quit during this step. Please choose a different option.\n")
          dt_verification_decision <- get_verification_decision(dt_choice_prompt)
        }
      }

      # update the specific row in df_daily_data based on the index
      j_index <- which(daily_plot_data_arg$DT_join == day_dt)

      altered_row <- alter_verification_column(dt_verification_decision, daily_plot_data_arg[j_index, ])

      daily_plot_data_arg <- daily_plot_data_arg %>%
        mutate(
          mean_verified = if_else(row_number() == j_index, altered_row$mean_verified, mean_verified),
          is_verified = if_else(row_number() == j_index, altered_row$is_verified, is_verified),
          verification_status = if_else(row_number() == j_index, altered_row$verification_status, verification_status)
        )
    }
    # inform the user that they are done verifying which ever df
    cat("Finished verifying: ", as.character(min(daily_plot_data_arg$DT_round)), " data.\n")

  }

  # inspect valid
  if (daily_verification_decision_arg == "INSPECT VALID") {

    # get the dt of the valid data to inspect to get indexes of that data
    valid_daily_plot_data_dt_list <- daily_plot_data_arg %>%
      filter(is.na(flag)) %>%
      pull(DT_join)

    # get the dt of the flagged data to get the non-inspect decision
    new_dt_list <- pull(daily_plot_data_arg, DT_join)
    non_dts_to_inspect <- new_dt_list[!new_dt_list %in% valid_daily_plot_data_dt_list]

    # decide what to do with the flagged data ----
    # for each non-inspect dt decision:
    if (length(non_dts_to_inspect) > 0) {

      # get non-inspect day decision:
      noninspect_prompt_text = paste("Would you like to (pass/fail/skip) non-inspect dts?\n(layer~", layer, ")(px/fx/sx): ")
      non_inspect_verification_decision <- get_verification_decision(noninspect_prompt_text)

      # Handle the case when the user selects an "INSPECT" or "QUIT" option inappropriately
      while (non_inspect_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME", "QUIT")) {
        if (non_inspect_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME")){
          cat("You cannot inspect non-inspect dts. Please choose a different option.\n")
          non_inspect_verification_decision <- get_verification_decision(noninspect_prompt_text)
        }
        if (non_inspect_verification_decision == "QUIT"){
          cat("You cannot quit during this step. Please choose a different option.\n")
          non_inspect_verification_decision <- get_verification_decision(noninspect_prompt_text)
        }
      }

      for (j in non_dts_to_inspect) {
        # j <- format(as.POSIXct(j, tz = "MST"), "%Y-%m-%d %H:%M:%S")

        j_index <- which((daily_plot_data_arg$DT_join) == j)

        altered_row <- alter_verification_column(non_inspect_verification_decision, daily_plot_data_arg[j_index, ])

        daily_plot_data_arg <- daily_plot_data_arg %>%
          mutate(
            mean_verified = if_else(row_number() == j_index, altered_row$mean_verified, mean_verified),
            is_verified = if_else(row_number() == j_index, altered_row$is_verified, is_verified),
            verification_status = if_else(row_number() == j_index, altered_row$verification_status, verification_status)
          )

        # some sort of print statement to show what happened to the non-target dts
      }
    }

    # update those flagged data in the df

    for (j in valid_daily_plot_data_dt_list) {

      # j <- format(as.POSIXct(j, tz = "MST"), "%Y-%m-%d %H:%M:%S")

      j_index <- which((daily_plot_data_arg$DT_join) == j)

      # create a data frame with a single row containing the desired DT_round value
      # and add it to the week plot and day plot
      flag_demarcation_data <- tibble(DT_round = daily_plot_data_arg$DT_round[j_index])

      daily_plot_object_new <- daily_plot_object_arg +
        geom_vline(data = flag_demarcation_data,
                   aes(xintercept = DT_round,
                       color = "data point of interest",
                       # size = 2, # fix this later
                       alpha = 0.05))

      weekly_plot_object_new <- weekly_plot_object_arg +
        geom_vline(data = flag_demarcation_data,
                   aes(xintercept = DT_round,
                       color = "data point of interest",
                       # size = 2, # fix this later
                       alpha = 0.05))

      # print ggarrange daily and weekly plots
      print(ggarrange(daily_plot_object_new, weekly_plot_object_new, nrow = 2, ncol = 1))

      # generate the prompt for the user to answer
      # the user should not be able to inspect again
      day_dt <- as.character(daily_plot_data_arg$DT_round[j_index]) # reformat so that this reads as a date
      dt_choice_prompt <- paste("Would you like to (pass/fail/skip) ", day_dt, "?\n(layer~", layer, ")(px/fx/sx)")
      dt_verification_decision <- get_verification_decision(dt_choice_prompt)

      # Handle the case when the user selects an "INSPECT" or "QUIT" option inappropriately
      while (layer == "subdaily" & dt_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME", "QUIT")) {
        if (dt_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME")){
          cat("You cannot inspect at this level. Please choose a different option.\n")
          dt_verification_decision <- get_verification_decision(dt_choice_prompt)
        }
        if (dt_verification_decision == "QUIT"){
          cat("You cannot quit during this step. Please choose a different option.\n")
          dt_verification_decision <- get_verification_decision(dt_choice_prompt)
        }
      }

      # update the specific row in df_daily_data based on the index
      # row_index <- which(as.character(daily_plot_data_arg$DT_round) == day_dt)
      altered_row <- alter_verification_column(dt_verification_decision, daily_plot_data_arg[j_index, ])

      daily_plot_data_arg <- daily_plot_data_arg %>%
        mutate(
          mean_verified = if_else(row_number() == j_index, altered_row$mean_verified, mean_verified),
          is_verified = if_else(row_number() == j_index, altered_row$is_verified, is_verified),
          verification_status = if_else(row_number() == j_index, altered_row$verification_status, verification_status)
        )
    }
    # inform the user that they are done verifying which ever df
    cat("Finished verifying: ", as.character(min(daily_plot_data_arg$DT_round)), " data.\n")
    # return the altered df
    # altered_df_list_arg[[i]] <- daily_plot_data_arg # make sure that this is the updated data frame

  }

  # inspect flagged
  if (daily_verification_decision_arg == "INSPECT FLAGGED") {

    # get the dt of the flagged data to inspect to get indexes of that data
    flagged_daily_plot_data_dt_list <- daily_plot_data_arg %>%
      filter(!is.na(flag)) %>%
      pull(DT_join)

    # get the dt of the flagged data to get the non-inspect decision
    new_dt_list <- pull(daily_plot_data_arg, DT_join)
    non_dts_to_inspect <- new_dt_list[!new_dt_list %in% flagged_daily_plot_data_dt_list]

    # decide what to do with the flagged data ----
    # for each non-inspect dt decision:
    if (length(non_dts_to_inspect) > 0) {

      # get non-inspect day decision:
      noninspect_prompt_text = paste("Would you like to (pass/fail/skip) non-inspect dts?\n(layer~", layer, ")(px/fx/sx): ")
      non_inspect_verification_decision <- get_verification_decision(noninspect_prompt_text)

      # Handle the case when the user selects an "INSPECT" or "QUIT" option inappropriately
      while (non_inspect_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME", "QUIT")) {
        if (non_inspect_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME")){
          cat("You cannot inspect non-inspect dts. Please choose a different option.\n")
          non_inspect_verification_decision <- get_verification_decision(noninspect_prompt_text)
        }
        if (non_inspect_verification_decision == "QUIT"){
          cat("You cannot quit during this step. Please choose a different option.\n")
          non_inspect_verification_decision <- get_verification_decision(noninspect_prompt_text)
        }
      }

      for (j in non_dts_to_inspect) {
        # j <- format(as.POSIXct(j, tz = "MST"), "%Y-%m-%d %H:%M:%S") # time zone

        j_index <- which((daily_plot_data_arg$DT_join) == j)

        altered_row <- alter_verification_column(non_inspect_verification_decision, daily_plot_data_arg[j_index, ])

        daily_plot_data_arg <- daily_plot_data_arg %>%
          mutate(
            mean_verified = if_else(row_number() == j_index, altered_row$mean_verified, mean_verified),
            is_verified = if_else(row_number() == j_index, altered_row$is_verified, is_verified),
            verification_status = if_else(row_number() == j_index, altered_row$verification_status, verification_status)
            )
        # some sort of print statement to show what happened to the non-target dts
      }
    }

    # update those flagged data in the df

    for (j in flagged_daily_plot_data_dt_list) {

      # j <- format(as.POSIXct(j, tz = "MST"), "%Y-%m-%d %H:%M:%S")

      j_index <- which((daily_plot_data_arg$DT_join) == j)

      # create a data frame with a single row containing the desired DT_round value
      # and add it to the week plot and day plot
      flag_demarcation_data <- tibble(DT_round = daily_plot_data_arg$DT_round[j_index])

      daily_plot_object_new <- daily_plot_object_arg +
        geom_vline(data = flag_demarcation_data,
                   aes(xintercept = DT_round,
                       color = "data point of interest",
                       # size = 2, # fix this later
                       alpha = 0.05))

      weekly_plot_object_new <- weekly_plot_object_arg +
        geom_vline(data = flag_demarcation_data,
                   aes(xintercept = DT_round,
                       color = "data point of interest",
                       # size = 2, # fix this later
                       alpha = 0.05))


      # print ggarrange daily and weekly plots
      print(ggarrange(daily_plot_object_new, weekly_plot_object_new, nrow = 2, ncol = 1))

      # generate the prompt for the user to answer
      # the user should not be able to inspect again
      day_dt <- as.character(daily_plot_data_arg$DT_round[j_index]) # reformat so that this reads as a date
      dt_choice_prompt <- paste("Would you like to (pass/fail/skip) ", day_dt, "?\n(layer~", layer, ")(px/fx/sx)")
      dt_verification_decision <- get_verification_decision(dt_choice_prompt)

      # Handle the case when the user selects an "INSPECT" or "QUIT" option inappropriately
      while (layer == "subdaily" & dt_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME", "QUIT")) {
        if (dt_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME")){
          cat("You cannot inspect at this level. Please choose a different option.\n")
          dt_verification_decision <- get_verification_decision(dt_choice_prompt)
        }
        if (dt_verification_decision == "QUIT"){
          cat("You cannot quit during this step. Please choose a different option.\n")
          dt_verification_decision <- get_verification_decision(dt_choice_prompt)
        }
      }

      # update the specific row in df_daily_data based on the index
      altered_row <- alter_verification_column(dt_verification_decision, daily_plot_data_arg[j_index, ])

      daily_plot_data_arg <- daily_plot_data_arg %>%
        mutate(
          mean_verified = if_else(row_number() == j_index, altered_row$mean_verified, mean_verified),
          is_verified = if_else(row_number() == j_index, altered_row$is_verified, is_verified),
          verification_status = if_else(row_number() == j_index, altered_row$verification_status, verification_status)
        )

    }
    # inform the user that they are done verifying which ever df
    cat("Finished verifying: ", as.character(min(daily_plot_data_arg$DT_round)), " data.\n")

  }

  # inspect some
  if (daily_verification_decision_arg == "INSPECT SOME") {

    # get which dts to inspect
    dt_intervals_to_inspect <- get_dt_inspection_decisions(daily_plot_data = daily_plot_data_arg)

    # decided what to do with the non-inspect data
    new_dt_list <- pull(daily_plot_data_arg, DT_round)
    non_dts_to_inspect <- new_dt_list[!new_dt_list %within% dt_intervals_to_inspect]

    # decide what to do with the flagged data ----
    # for each non-inspect dt decision:
    if (length(non_dts_to_inspect) > 0) {

      # get non-inspect day decision:
      noninspect_prompt_text = paste("Would you like to (pass/fail/skip) non-inspect dts?\n(layer~", layer, ")(px/fx/sx): ")
      non_inspect_verification_decision <- get_verification_decision(noninspect_prompt_text)

      # Handle the case when the user selects an "INSPECT" or "QUIT" option inappropriately
      while (non_inspect_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME", "QUIT")) {
        if (non_inspect_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME")){
          cat("You cannot inspect non-inspect dts. Please choose a different option.\n")
          non_inspect_verification_decision <- get_verification_decision(noninspect_prompt_text)
        }
        if (non_inspect_verification_decision == "QUIT"){
          cat("You cannot quit during this step. Please choose a different option.\n")
          non_inspect_verification_decision <- get_verification_decision(noninspect_prompt_text)
        }
      }

      for (j in non_dts_to_inspect) {
        j_time <- format(as.POSIXct(j, tz = "MST"), "%H:%M:%S")
        j <- if_else(
          j_time == "00:00:00",
          format(as.POSIXct(j, tz = "MST"), "%Y-%m-%d"),
          format(as.POSIXct(j, tz = "MST"), "%Y-%m-%d %H:%M:%S")
        )

        j_index <- which((daily_plot_data_arg$DT_join) == j)

        # daily_plot_data_arg[j_index, ] <- alter_verification_column(non_inspect_verification_decision, daily_plot_data_arg[j_index, ])

        altered_row <- alter_verification_column(non_inspect_verification_decision, daily_plot_data_arg[j_index, ])

        daily_plot_data_arg <- daily_plot_data_arg %>%
          mutate(
            mean_verified = if_else(row_number() == j_index, altered_row$mean_verified, mean_verified),
            is_verified = if_else(row_number() == j_index, altered_row$is_verified, is_verified),
            verification_status = if_else(row_number() == j_index, altered_row$verification_status, verification_status)
          )
        # some sort of print statement to show what happened to the non-target dts
      }
    }

    # update those flagged data in the df

    for (j in dt_intervals_to_inspect) { # make sure that j can not be 0

      interval_df <- daily_plot_data_arg %>%
        filter(DT_round %within% j)

      for (k in 1:nrow(interval_df)) {

        k_index <- which((daily_plot_data_arg$DT_round) == interval_df$DT_round[k])

        # create a data frame with a single row containing the desired DT_round value
        # and add it to the week plot and day plot
        flag_demarcation_data <- tibble(DT_round = daily_plot_data_arg$DT_round[k_index])

        daily_plot_object_new <- daily_plot_object_arg +
          geom_vline(data = flag_demarcation_data,
                     aes(xintercept = DT_round,
                         color = "data point of interest",
                         # size = 2, # fix this later
                         alpha = 0.05))

        weekly_plot_object_new <- weekly_plot_object_arg +
          geom_vline(data = flag_demarcation_data,
                     aes(xintercept = DT_round,
                         color = "data point of interest",
                         # size = 2, # fix this later
                         alpha = 0.05))


        # print ggarrange daily and weekly plots
        print(ggarrange(daily_plot_object_new, weekly_plot_object_new, nrow = 2, ncol = 1))

        # generate the prompt for the user to answer
        # the user should not be able to inspect again
        day_dt <- as.character(daily_plot_data_arg$DT_round[k_index]) # reformat so that this reads as a date
        dt_choice_prompt <- paste("Would you like to (pass/fail/skip) ", day_dt, "?\n(layer~", layer, ")(px/fx/sx)")
        dt_verification_decision <- get_verification_decision(dt_choice_prompt)

        # Handle the case when the user selects an "INSPECT" or "QUIT" option inappropriately
        while (layer == "subdaily" & dt_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME", "QUIT")) {
          if (dt_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME")){
            cat("You cannot inspect at this level. Please choose a different option.\n")
            dt_verification_decision <- get_verification_decision(dt_choice_prompt)
          }
          if (dt_verification_decision == "QUIT"){
            cat("You cannot quit during this step. Please choose a different option.\n")
            dt_verification_decision <- get_verification_decision(dt_choice_prompt)
          }
        }

        # update the specific row in df_daily_data based on the index
        # row_index <- which(as.character(daily_plot_data_arg$DT_round) == day_dt)
        # daily_plot_data_arg[k_index, ] <- alter_verification_column(dt_verification_decision, daily_plot_data_arg[k_index, ])

        altered_row <- alter_verification_column(dt_verification_decision, daily_plot_data_arg[k_index, ])

        daily_plot_data_arg <- daily_plot_data_arg %>%
          mutate(
            mean_verified = if_else(row_number() == k_index, altered_row$mean_verified, mean_verified),
            is_verified = if_else(row_number() == k_index, altered_row$is_verified, is_verified),
            verification_status = if_else(row_number() == k_index, altered_row$verification_status, verification_status)
          )
      }
    }
    # inform the user that they are done verifying which ever df
    cat("Finished verifying: ", as.character(min(daily_plot_data_arg$DT_round)), " data.\n")
  }

  return(daily_plot_data_arg)

}

add_depth_shift_flag <- function(df){


  cal_flag <- function(df, condition_arg, description_arg) {
    df <- df %>% mutate(depth_change = case_when(
      {{condition_arg}} ~ if_else(is.na(depth_change), paste(description_arg),
                                  ifelse(!grepl(description_arg, depth_change), paste(depth_change, description_arg, sep = ";\n"), depth_change)),
      TRUE ~ depth_change))
    return(df)
  }

  # Function to add a column if it doesn't already exist
  add_column_if_not_exists <- function(df, column_name, default_value) {
    if (!column_name %in% colnames(df)) {
      df <- df %>% mutate(!!sym(column_name) := default_value)
    }
    return(df)
  }

    shift_dates <- read_csv('data/qaqc/level_shifts.csv') %>%
      filter(type == "sonde moved") %>%
      add_column_if_not_exists(., "depth_change", NA) %>%
      mutate(DT_join = as.character(ymd_hms(DT_join))) %>%
      left_join(df, ., by = c("site", "DT_join")) %>%
      cal_flag(type == "sonde moved", "sonde moved") %>%
      select(-type)

    return(shift_dates)

}


add_do_drops <- function(df){

if("DO" %in% df$parameter){

  df <- df %>%
    add_flag(back1 - mean >= 0.5 & front1 - mean >= 0.5, "DO interference")

  return(df)

} else {

  return(df)

}

  }

add_drift_flag <- function(df){

  # Only test for biofilm growth on turbidity sensors
  if("Turbidity" %in% df$parameter){

    # subset data to turbidity and conductivity only
    sub <- df %>%
      dplyr::filter(parameter %in% c("Turbidity", "Specific Conductivity")) %>%
      dplyr::select(DT_round, DT_join, parameter, mean) %>%
      tidyr::pivot_wider(names_from = parameter, values_from = mean)
    names(sub) <- make.names(names(sub))

    # Check if a given window's worth of sensor data has a linear
    # relationship with time (indicates a "steady" drift)
    progressive_drift <- function(x) {
      # Only assess time series with less than 90% missing data in
      # the rolling window:
      if(length(x[!is.na(x)]) > (length(x) - (length(x)*0.1))){

        # Fit linear model
        model <- lm(x ~ c(1:length(x)), na.action = na.omit)

        # Extract R-squared value
        r_squared <- summary(model)$r.squared

        # Return R-squared value
        return(r_squared)

      } else {

        # If not enough data, set slope to 0
        no_slope <- 0

        return(no_slope)

      }
    }

    # Function to check if a selected time window's mean R-squared value is
    # at least 60% (i.e., fairly high linear relationship with time indicates a "steady" drift)
    check_too_steady <- function(x) {
      mean(x) >= 0.60
    }

    # Function that uses all the functions above to see if a given time window's R-squared with time is strong. If the one-day OR the three-day slope
    # for a selected parameter is greater than 60%, we determine it has "failed" (i.e., drift seems to exist)

    biofilm_tester <- function(data = sub, col){

      data %>%
        data.table::data.table() %>%
        #bind_rows(flagged_data_dfs[["prospect-Turbidity"]]) %>%
        # dplyr::filter(lubridate::as_date(DT_round) >= "2022-09-10" & lubridate::as_date(DT_round) <= "2022-10-15") %>%
        dplyr::mutate(r2_s_right = data.table::frollapply(!!sym(col), n = 96, FUN = progressive_drift, align = "right", fill = NA),
                      r2_s_center = data.table::frollapply(!!sym(col), n = 96, FUN = progressive_drift, align = "left", fill = NA),
                      r2_l_right = data.table::frollapply(!!sym(col), n = 288, FUN = progressive_drift, align = "right", fill = NA),
                      r2_l_center = data.table::frollapply(!!sym(col), n = 288, FUN = progressive_drift, align = "left", fill = NA),
                      tightest_r = pmax(r2_s_center, r2_s_right, r2_l_center, r2_l_right, na.rm = TRUE),
                      failed = data.table::frollapply(tightest_r, n = 96, FUN = check_too_steady, align = "right", fill = NA)) %>%
        dplyr::select(!!(col) := "failed",
                      DT_join)
    }

    biofouling <- c("Turbidity", "Specific.Conductivity") %>%
      purrr::map(~biofilm_tester(data = sub, col = .)) %>%
      dplyr::bind_cols() %>%
      dplyr::rename(DT_join = DT_join...2) %>%
      dplyr::mutate(parameter = "Turbidity") %>%
      dplyr::right_join(., df, by = c("DT_join", "parameter")) %>%
      # If a steady slope is existing through time for turbidity, but NOT for conductivity, it is likely
      # sensor bio-fouling
      add_flag(., Turbidity == 1 & Specific.Conductivity != 1 & !grepl("drift", flag), "drift") %>%
      dplyr::select(-c(Turbidity, DT_join...4, Specific.Conductivity))

    return(biofouling)

  } else {

    return(df)}

}
#
# a=ggplot(drift_tested)+
#   geom_point(aes(x = DT_round, y = r_squaredSC_short ), color = "blue") +
#   geom_point(aes(x = DT_round, y = r_squaredT_short ), color = "red") +
#   geom_point(aes(x = DT_round, y = r_squaredDO_short), color = "black")
#
# b=ggplot(drift_tested)+
#   geom_point(aes(x = DT_round, y = r_squaredSC_long ), color = "blue") +
#   geom_point(aes(x = DT_round, y = r_squaredT_long ), color = "red") +
#   geom_point(aes(x = DT_round, y = r_squaredDO_long), color = "black")
#
# c=ggplot(drift_tested %>% filter(as_date(DT_round) >= "2022-10-10")) +
#   #geom_point(aes(x = DT_round, y = Specific.Conductivity)) +
#   geom_point(aes(x = DT_round, y = Turbidity, color = as.character(failed))) #+
# #geom_point(aes(x = DT_round, y = DO))
#
# ggpubr::ggarrange(a, b, c, ncol = 1)
#


#' @title Add field related flags to a data frame based on field notes.
#'
#' @description
#' A function that checks 3 different and separate conditions related to the
#' field notes and adds the corresponding flags accordingly. The "sonde not
#' employed" flag is added if the sonde was not employed at the site.
#' "site visit" flag is added if the last site visit date is the same as the
#' date of the data. "sv window" flag is added if site visit flag is detected
#' within a 45 minute window. Note that this function will only work on data
#' that has already been joined to field notes.
#'
#' @param df A data frame with a `flag` column.
#'
#' @return A data frame with a `flag` column that has been updated with the
#' relevant field note flags.
#'
#' @examples
#' add_range_flags(df = all_data_flagged$`archery-Actual Conductivity`)
#' add_range_flags(df = all_data_flagged$`boxelder-Temperature`)

add_field_flag <- function(df) {

  df <- df %>%
    # flag when sonde was not employed in the river
    add_flag(sonde_employed == 1, "sonde not employed") %>% # removing the instances where we check for the flag description in the flag column when we are flagging. This is already done in `add_flag()`
    # flag when sonde was handled in a site visit
    add_flag(as.character(last_site_visit) == as.character(DT_round), "site visit")

  # Add flags for the next 60 minutes after a site visit
  for (i in 1:4) {
    df <- df %>%
      add_flag(lag(str_detect(flag, "site visit"), n = i), "sv window")
  }

  # ... and the first 15 minutes before:
  df <- df %>%
    add_flag(lead(str_detect(flag, "site visit"), n = 1), "sv window")

  return(df)

}

#' @title Underlying function for flagging data.
#'
#' @description
#' This function adds a flag to the `flag` column of a given data frame based on
#' specified conditions for each row. The name of the flag is provided
#' by the user.
#'
#' @param df A data frame with a `flag` column.
#'
#' @param condition_arg A logical statement that is evaluated in the context of
#' the data frame.
#'
#' @param description_arg Flag added to the `flag` column.
#'
#' @returns
#' An identical data frame with a `flag` column that has been updated with the
#' flag description provided.
#'
#' @examples
# add_flag(df = all_data_flagged$`archery-Actual Conductivity`, condition_arg = mean >= 100, description_arg = "exceeds 100")
# add_flag(df = all_data_flagged$`archery-Actual Conductivity`, condition_arg = mean <= 10, description_arg = "below 10")

add_flag <- function(df, condition_arg, description_arg) {
  df <- df %>% mutate(flag = case_when(
    {{condition_arg}} ~ if_else(is.na(flag), paste(description_arg),
                                ifelse(!grepl(description_arg, flag), paste(flag, description_arg, sep = ";\n"), flag)),
    TRUE ~ flag))
  return(df)
}

#' @title Add a flag if the water temperature is freezing.
#'
#' @description
#' A function designed to append the 'frozen' flag to a row if the value
#' in the `mean` column is less than or equal to 0.
#'
#' @param df A data frame with a `flag` column.
#'
#' @return A data frame with a `flag` column that has been updated with the
#' 'frozen' flag.
#'
#' @examples
#' add_frozen_flag(df = all_data_flagged$`archery-Actual Conductivity`)
#' add_frozen_flag(df = all_data_flagged$`boxelder-Temperature`)
#'
#' @seealso [flag_all_data()]

add_frozen_flag <- function(df){

  # create a df of temperature for each site
  temperature <- df %>%
    data.table::data.table() %>%
    dplyr::select(DT_round, DT_join, parameter, mean) %>%
    dplyr::filter(parameter == "Temperature") %>%
    dplyr::select(DT_join, Temperature = mean)

    # add "temperature" column to df:
    temperature_checked <- df %>%
      dplyr::left_join(., temperature, by = "DT_join") %>%
      # If water temperature is freezing, flag all parameters
      add_flag(., Temperature <= 0, "frozen") %>%
      # remove the temp column so df is identical in structure to OG df
      dplyr::select(-Temperature)

    return(temperature_checked)

}

#' @title Add malfunction flags to a data frame based on field notes and HydroVu.
#'
#' @description
#' This function adds the 'sensor malfunction' flag to a data frame by the dates,
#' sites, and parameters that members of the lab know contain erroneous data due
#' to sensor malfunctions. Note that this flag is used in two instances: first
#' when removing erroneous data from our statistics calculations and again during
#' the actual flagging step in the QAQC pipeline.
#'
#' @param df A data frame with a `flag` column.
#'
#' @param malfunction_records The malfunction records pulled from mWater
#'
#' @return A data frame with a `flag` column that has been updated with the
#' 'sensor malfunction' flag.
#'
#' @examples
#' add_range_flags(df = all_data_flagged$`archery-Actual Conductivity`)
#' add_range_flags(df = all_data_flagged$`boxelder-Temperature`)
#'
#' @seealso [grab_mWater_malfunction_notes()]

add_malfunction_flag <- function(df, malfunction_records){

  # Filter malfunction records for relevant information
  df_site <- unique(df$site)
  df_parameter <- unique(df$parameter)

  malfunction_records_filtered <- malfunction_records %>%
    filter(site == df_site) %>%
    # if parameter is NA, that means entire sonde was malfunctioning. Otherwise,
    # just the parameter listed:
    filter(is.na(parameter) | parameter == df_parameter) %>%
    mutate(end_DT = ifelse(is.na(end_DT), ymd_hms("9999-12-31 23:59:59", tz = "MST"), end_DT)) %>%
    mutate(end_DT = as.POSIXct(end_DT, tz = "MST")) %>%
    rowid_to_column()

  if (nrow(malfunction_records_filtered > 0)) {

    drift <- malfunction_records_filtered %>%
      filter(grepl("grime|gunk|drift|biofoul|biofilm", notes, ignore.case = TRUE))

    burial <- malfunction_records_filtered %>%
      filter(grepl("buried|burial|bury|sediment|roots", notes, ignore.case = TRUE))

    depth_funk <- malfunction_records_filtered %>%
      filter(grepl("improper level calibration", notes, ignore.case = TRUE))

    unsubmerged <- malfunction_records_filtered %>%
      filter(grepl("not submerged", notes, ignore.case = TRUE))

    general_malfunction <- malfunction_records_filtered %>%
      filter(!rowid %in% drift$rowid & !rowid %in% burial$rowid &
               !rowid %in% depth_funk$rowid & !rowid %in% unsubmerged$rowid)

    # make a list of date intervals where sensor was malfunctioning
    drift_interval_list <- map2(
      .x = drift$start_DT,
      .y = drift$end_DT,
      .f = ~interval(.x, .y, tz = "MST"))

    burial_interval_list <- map2(
      .x = burial$start_DT,
      .y = burial$end_DT,
      .f = ~interval(.x, .y, tz = "MST"))

    depth_interval_list <- map2(
      .x = depth_funk$start_DT,
      .y = depth_funk$end_DT,
      .f = ~interval(.x, .y, tz = "MST"))

    unsubmerged_interval_list <- map2(
      .x = unsubmerged$start_DT,
      .y = unsubmerged$end_DT,
      .f = ~interval(.x, .y, tz = "MST"))

    general_interval_list <- map2(
      .x = general_malfunction$start_DT,
      .y = general_malfunction$end_DT,
      .f = ~interval(.x, .y, tz = "MST"))

    # If DT_round is within any of the intervals where a sensor was malfunctioning, flag it
    try(df <- df %>%
          add_flag(DT_round %within% burial_interval_list, "sonde burial"))

    try(df <- df %>%
          add_flag(DT_round %within% drift_interval_list, "sensor biofouling"))

    try(df <- df %>%
          add_flag(DT_round %within% depth_interval_list, "depth calibration malfunction"))

    try(df <- df %>%
          add_flag(DT_round %within% unsubmerged_interval_list, "sonde unsubmerged"))

    try(df <- df %>%
          add_flag(DT_round %within% general_interval_list, "sensor malfunction"))

  }


  return(df)

}

#' @title Add NA flags to a data frame based on `mean` column.
#'
#' @description
#' A function designed to append the 'missing data' flag to a row if the `mean`
#' column in that row contains NA values.
#'
#' @param df A data frame with a `flag` column.
#'
#' @return A data frame with a `flag` column that has been updated with the
#' 'missing data' flag.
#'
#' @examples
#' add_range_flags(df = all_data_flagged$`archery-Actual Conductivity`)
#' add_range_flags(df = all_data_flagged$`boxelder-Temperature`)

add_na_flag <- function(df){
  df <- df %>%
    add_flag(is.na(mean), "missing data")
}

#' @title Add a flag if the value of a parameter exceeds its realistic range.
#'
#' @description
#' A function designed to append the 'outside of sensor realistic range' flag to
#' a row if the `mean` column in that row exceeds the realistic ranges that are
#' set in the `src/qaqc/sensor_real_thresholds.yml` file. These are instances
#' where a value exceeds the expected ranges for the Poudre. Note that this is
#' currently only aplied to pH, SC, and Temperature.
#'
#' @param df A data frame with a `flag` column.
#'
#' @return A data frame with a `flag` column that has been updated with the
#' 'outside of sensor realistic range' flag.
#'
#' add_range_flags(df = all_data_flagged$`archery-Actual Conductivity`)
#' add_range_flags(df = all_data_flagged$`boxelder-Temperature`)

add_realistic_flag <- function(df){

  sensor_realistic_ranges <- read_yaml("data/qaqc/sensor_real_thresholds.yml")

  # make this a non yaml solution and add it to the threshold table
  # get the parameter from the parameter column in the df of interest
  parameter_name <- unique(na.omit(df$parameter))
  # Pull the sensor specification range from the yaml file
  sensor_min <- eval(parse(text = sensor_realistic_ranges[[parameter_name]]$min))
  sensor_max <- eval(parse(text = sensor_realistic_ranges[[parameter_name]]$max))

  df <- df %>%
    # adding sensor range flags
    add_flag(parameter == parameter_name & (mean < sensor_min | mean > sensor_max),
             "outside of sensor realistic range")

    return(df)

}

#' @title Add a flag if the value in the `mean` column repeats sequentially.
#'
#' @description
#' A function designed to append the 'repeated value' flag to a row if the value
#' in the `mean` column is equal to the previous or next value in the `mean`
#' column.
#'
#' @param df A data frame with a `flag` column.
#'
#' @return A data frame with a `flag` column that has been updated with the
#' 'repeated value' flag.
#'
#' add_range_flags(df = all_data_flagged$`archery-Actual Conductivity`)
#' add_range_flags(df = all_data_flagged$`boxelder-Temperature`)

add_repeat_flag <- function(df){
  df <- df %>%
    add_flag((mean == front1 | mean == back1), "repeated value")
}

#' @title Add flags related to calculated parameter seasonal ranges.
#'
#' @description
#' A function that checks 2 different and separate conditions related to the
#' calculated parameter seasonal ranges. The 'outside of seasonal range' flag is
#' added if the `mean` value is outside the seasonal 1st - 99th percentile. The
#' 'slope violation' flag is added if the `slope_ahead` or `slope_behind` value
#' is greater than or equal to the `t_slope_behind_99` value, or if the
#' `slope_ahead` or `slope_behind` value is less than or equal to the
#' `t_slope_behind_01` value. The `t_slope_behind_99` and `t_slope_behind_01`
#' values are specific to each site-parameter. Note that this is the version of
#' the function used in `flag_all_data()`.
#'
#' @param df A data frame with a `flag` column.
#'
#' @return A data frame with a `flag` column that has been updated with the
#' relevant calculated seasonal range flags.
#'
#' @examples
#' add_range_flags(df = all_data_flagged$`archery-Actual Conductivity`)
#' add_range_flags(df = all_data_flagged$`boxelder-Temperature`)
#'
#' @seealso [flag_all_data()]

add_seasonal_flag <- function(df) {

  # get the site name from the site column in the df of interest
  site_name <- unique(na.omit(df$site))
  # get the parameter from the parameter column in the df of interest
  parameter_name <- unique(na.omit(df$parameter))

  lookup <- threshold_lookup %>%
    filter(site == site_name & parameter == parameter_name) %>%
    select(!c(site, parameter))

  df <- df %>%
    # Using seasonal cut-offs...
    left_join(lookup, by = "season") %>%
    # ... flag obs that are outside the seasonal 1st - 99th percentile range:
    add_flag((mean < t_mean01 | mean > t_mean99), "outside of seasonal range") %>%
    # flag obs whose slope is outside the 1st - 99th percentile range:
    add_flag(((slope_ahead >= t_slope_behind_99 | slope_behind >= t_slope_behind_99) |
              (slope_ahead <= t_slope_behind_01 | slope_behind <= t_slope_behind_01)), "slope violation")

  return(df)

}

# Adding flags related to sensor specification ranges; these are instances where
# a value exceeds the manufacturer specified limits.

#' @title Add a flag if the value of a parameter exceeds its sensor specification range.
#'
#' @description
#' A function designed to append the 'outside of sensor specification range' flag to
#' a row if the `mean` column in that row exceeds the sensor specification ranges that are
#' set in the `src/qaqc/sensor_spec_thresholds.yml` file. These are instances
#' where a value exceeds the expected ranges for the Poudre based on the sensor
#' manufacturer's specifications.
#'
#' @param df A data frame with a `flag` column.
#'
#' @return A data frame with a `flag` column that has been updated with the
#' 'outside of sensor specification range' flag.
#'
#' @examples
#' add_spec_flag(df = all_data_flagged$`archery-Actual Conductivity`)
#' add_spec_flag(df = all_data_flagged$`boxelder-Temperature`)
#'
#' @seealso [flag_all_data()]

add_spec_flag <- function(df){

  sensor_spec_ranges <- read_yaml("data/qaqc/sensor_spec_thresholds.yml")

  # make this a non yaml solution and add it to the threshold table
  # get the parameter from the parameter column in the df of interest
  parameter_name <- unique(na.omit(df$parameter))
  # Pull the sensor specification range from the yaml file
  sensor_min <- eval(parse(text = sensor_spec_ranges[[parameter_name]]$min))
  sensor_max <- eval(parse(text = sensor_spec_ranges[[parameter_name]]$max))

  df <- df %>%
    # adding sensor range flags
    add_flag(parameter == parameter_name & (mean < sensor_min | mean > sensor_max) & !grepl("outside of sensor specification range", flag),
             "outside of sensor specification range") %>%

    return(df)

}

# Add flags to the `flag` column of a dataframe based on large swaths of suspect data.
# "24hr anomaly flag" is added if more than 50% of the data points in a 24 hour window are flagged.
# "anomaly window" flag is added if the point is included in a 24hr anomaly.
# @param df A dataframe with a `auto_cleaned_flag` column.
# @return A dataframe with a `auto_cleaned_flag` column that has been updated with the relevant large anomaly flags.
# @examples
# add_large_anomaly_flags(df = all_data_flagged$`archery-Actual Conductivity`)
# add_large_anomaly_flags(df = all_data_flagged$`boxelder-Temperature`)

add_suspect_flag_full <- function(df) {
  # these are the flags that we don't want to perform this exercise across
  auto_cleaned_flag_string <- "sonde not employed|missing data|site visit|sv window|sonde burial|sensor biofouling|depth calibration malfunction|sonde unsubmerged|sensor malfunction"

  # Define a function to check if a given 3-hour window has >= 50% fails
  check_3_hour_window_fail <- function(x) {
    sum(x) / length(x) >= 0.5
  }

  df_test <- df %>%
    dplyr::mutate(auto_cleaned_flag_binary = ifelse((is.na(auto_cleaned_flag) | grepl(auto_cleaned_flag_string, auto_cleaned_flag) | auto_cleaned_flag_string == "suspect data"), 0, 1)) %>%
    #arrange(timestamp) %>%
    dplyr::mutate(over_50_percent_fail_window = zoo::rollapply(auto_cleaned_flag_binary, width = 12, FUN = check_3_hour_window_fail, fill = NA, align = "right")) %>%
    add_flag(over_50_percent_fail_window == TRUE & is.na(auto_cleaned_flag), "suspect data")

  return(df_test)

}



#' @title Add flags to the `flag` column of a dataframe based on large swaths of suspect data.
#'
#' @description
#' "24hr anomaly flag" is added if more than 50% of the data points in a 24 hour window are flagged.
#' "anomaly window" flag is added if the point is included in a 24hr anomaly.
#'
#' @param df A data frame with a `flag` column.
#'
#' @return A data frame with a `flag` column that has been updated with the relevant calculated seasonal range flags.
#'
#' @examples
#' add_range_flags(df = all_data_flagged$`archery-Actual Conductivity`)
#' add_range_flags(df = all_data_flagged$`boxelder-Temperature`)
#'
#' @seealso [flag_all_data()]

add_suspect_flag <- function(df) {

  flag_string <- "sonde not employed|missing data|site visit|sv window|suspect data" # include the flag added by this function

  # Define a function to check if a given 3-hour window has >= 50% fails
  check_3_hour_window_fail <- function(x) {
    sum(x) / length(x) >= 0.5
  }

  df_test <- df %>%
    dplyr::mutate(flag_binary = ifelse((is.na(flag) | grepl(flag_string, flag)), 0, 1)) %>%
    #arrange(timestamp) %>%
    dplyr::mutate(over_50_percent_fail_window = ifelse(is.na(over_50_percent_fail_window),
                                                zoo::rollapply(flag_binary, width = 12, FUN = check_3_hour_window_fail, fill = NA, align = "right"),
                                                over_50_percent_fail_window)) %>%
    add_flag(over_50_percent_fail_window == TRUE & is.na(flag), "suspect data")

  return(df_test)

}

add_threshold_lines <- function(plot, plot_data, site_arg, parameter_arg) {

  # pull in threshold data (i don't like that I do this everytime the function is run)
  real_thresholds <- read_csv("data/qaqc/realistic_thresholds.csv", show_col_types = FALSE) %>%
    filter(parameter == parameter_arg)
  sensor_thresholds <- yaml::read_yaml("data/qaqc/sensor_spec_thresholds.yml")[[parameter_arg]]%>% #filter for parameter_arg
    #turn into tibble with min/max
    bind_rows()

  seasonal_thresholds <- read_csv('data/qaqc/seasonal_thresholds_virridy.csv', show_col_types = FALSE) %>%
    #to do: Check to make sure seasonal thresholds csv is not necessary
    #bind_rows(read_csv('data/qaqc/seasonal_thresholds.csv', show_col_type = FALSE),
    distinct(site, parameter, season, .keep_all = TRUE) %>%
    #read_csv("data/qaqc/seasonal_thresholds_virridy.csv", show_col_types = FALSE) %>%
    filter(parameter == parameter_arg,
           site == site_arg)

  # Determine the unique seasons in plot_data
  unique_seasons <- unique(plot_data$season)

  # Filter seasonal_thresholds for the seasons present in plot_data
  seasonal_thresholds <- seasonal_thresholds %>%
    filter(season %in% unique_seasons)

  if (nrow(seasonal_thresholds) > 1) { # make sure this works

    season_1 <- case_when(
      all(c("winter_baseflow", "snowmelt") %in% unique_seasons) ~ "winter_baseflow",
      all(c("snowmelt", "monsoon") %in% unique_seasons) ~ "snowmelt",
      all(c("monsoon", "fall_baseflow") %in% unique_seasons) ~ "monsoon",
      all(c("fall_baseflow", "winter_baseflow") %in% unique_seasons) ~ "fall_baseflow",
      TRUE ~ NA_character_
    )
    season_2 <- case_when(
      all(c("winter_baseflow", "snowmelt") %in% unique_seasons) ~ "snowmelt",
      all(c("snowmelt", "monsoon") %in% unique_seasons) ~ "monsoon",
      all(c("monsoon", "fall_baseflow") %in% unique_seasons) ~ "fall_baseflow",
      all(c("fall_baseflow", "winter_baseflow") %in% unique_seasons) ~ "winter_baseflow",
      TRUE ~ NA_character_
    )

    seasonal_threshold_s1_quantiles <- unname(quantile(c(seasonal_thresholds %>% filter(season == season_1) %>% pull(t_mean01),
                                                         seasonal_thresholds %>% filter(season == season_1) %>% pull(t_mean99)),
                                                       c(0.1, 0.9)))

    seasonal_threshold_s2_quantiles <- unname(quantile(c(seasonal_thresholds %>% filter(season == season_2) %>% pull(t_mean01),
                                                         seasonal_thresholds %>% filter(season == season_2) %>% pull(t_mean99)),
                                                       c(0.1, 0.9)))

    slice_data <- plot_data %>%
      group_by(season) %>%
      slice(1) %>%
      ungroup()

    site_data <- filter(plot_data, site == site_arg)

    if (!all(is.na(site_data$mean))) {
      # Lower bound
      if ( # !is.infinite(min(site_data$mean, na.rm = TRUE)) &
        (min(site_data$mean, na.rm = TRUE) <= seasonal_threshold_s1_quantiles[1] | min(site_data$mean, na.rm = TRUE) <= seasonal_threshold_s2_quantiles[1])) {

        # Xs
        start_x_s1 <- min(slice_data$DT_round)

        transition_date <- slice_data %>%
          filter(season == season_2) %>%
          pull(DT_round)

        end_x_s2 <- ceiling_date(max(plot_data$DT_round), "day")

        #Ys
        y_lower_s1 <- seasonal_thresholds %>%
          filter(season == season_1) %>%
          pull(t_mean01)

        y_lower_s2 <- seasonal_thresholds %>%
          filter(season == season_2) %>%
          pull(t_mean01)

        plot <- plot +
          # season 1
          ## lower bound
          geom_segment(aes(x = start_x_s1, y = y_lower_s1,
                           xend = transition_date, yend = y_lower_s1,
                           color = "Seasonal Min", linetype = "Seasonal Min")) +
          # season 2
          ## lower bound
          geom_segment(aes(x = transition_date, y = y_lower_s2,
                           xend = end_x_s2, yend = y_lower_s2,
                           color = "Seasonal Min", linetype = "Seasonal Min"))

      }

      # Upper bound
      if (# !is.infinite(max(site_data$mean, na.rm = TRUE)) &
        (max(site_data$mean, na.rm = TRUE) >= seasonal_threshold_s1_quantiles[2] | max(site_data$mean, na.rm = TRUE) >= seasonal_threshold_s2_quantiles[2])) {

        # Xs
        start_x_s1 <- min(slice_data$DT_round)


        transition_date <- slice_data %>%
          filter(season == season_2) %>%
          pull(DT_round)

        end_x_s2 <- ceiling_date(max(plot_data$DT_round), "day")

        #Ys
        y_upper_s1 <- seasonal_thresholds %>%
          filter(season == season_1) %>%
          pull(t_mean99)

        y_upper_s2 <- seasonal_thresholds %>%
          filter(season == season_2) %>%
          pull(t_mean99)

        plot <- plot +
          # season 1
          ## lower bound
          geom_segment(aes(x = start_x_s1, y = y_upper_s1,
                           xend = transition_date, yend = y_upper_s1,
                           color = "Seasonal Max", linetype = "Seasonal Max")) +
          # season 2
          ## lower bound
          geom_segment(aes(x = transition_date, y = y_upper_s2,
                           xend = end_x_s2, yend = y_upper_s2,
                           color = "Seasonal Max", linetype = "Seasonal Max"))

      }

      # real thresholds
      real_thresholds_quantiles <- unname(quantile(c(real_thresholds$min, real_thresholds$max), c(0.1, 0.9)))

      if (#!is.infinite(min(site_data$mean, na.rm = TRUE))  &
        (min(site_data$mean, na.rm = TRUE) <= real_thresholds_quantiles[1])) {
        plot <- plot +
          geom_hline(aes(yintercept = real_thresholds$min,
                         color = "Real Min",
                         linetype = "Real"))
      }
      if (#!is.infinite(max(site_data$mean, na.rm = TRUE)) &
        (max(site_data$mean, na.rm = TRUE) >= real_thresholds_quantiles[2])) {
        plot <- plot +
          geom_hline(aes(yintercept = real_thresholds$max,
                         color = "Real Max",
                         linetype = "Real"))
      }

      # sensor thresholds
      sensor_thresholds_quantiles <- unname(quantile(c(sensor_thresholds$min, sensor_thresholds$max), c(0.1, 0.9)))

      # if ((min(site_data$mean, na.rm = TRUE) <= sensor_thresholds_quantiles[1])) {
      #   plot <- plot +
      #     geom_hline(aes(yintercept = sensor_thresholds$mix, # *** this needs to be min
      #                    color = "Sensor Min",
      #                    linetype = "Sensor"))
      # }

      if ((max(site_data$mean, na.rm = TRUE) >= sensor_thresholds_quantiles[2])) {
        plot <- plot +
          geom_hline(aes(yintercept = sensor_thresholds$max,
                         color = "Sensor Max",
                         linetype = "Sensor"))
      }
      return(plot)
    }

  } else if (nrow(seasonal_thresholds) == 1){

    site_data <- filter(plot_data, site == site_arg)

    # Filter seasonal_thresholds for the seasons present in site_data
    seasonal_thresholds <- seasonal_thresholds %>%
      filter(season %in% unique(site_data$season))

    if (!all(is.na(site_data$mean))) {
      # seasonal thresholds
      seasonal_thresholds_quantiles <- unname(quantile(c(seasonal_thresholds$t_mean01, seasonal_thresholds$t_mean99), c(0.1, 0.9)))

      if ((min(site_data$mean, na.rm = TRUE) <= seasonal_thresholds_quantiles[1]) == TRUE){
        plot <- plot +
          geom_hline(aes(yintercept = seasonal_thresholds$t_mean01,
                         color = "Seasonal Min",
                         linetype = "Seasonal"))
      }

      if ((max(site_data$mean, na.rm = TRUE) >= seasonal_thresholds_quantiles[2]) == TRUE) {
        plot <- plot +
          geom_hline(aes(yintercept = seasonal_thresholds$t_mean99,
                         color = "Seasonal Max",
                         linetype = "Seasonal"))
      }

      # real thresholds
      real_thresholds_quantiles <- unname(quantile(c(real_thresholds$min, real_thresholds$max), c(0.1, 0.9)))

      if (min(site_data$mean, na.rm = TRUE) <= real_thresholds_quantiles[1]) {
        plot <- plot +
          geom_hline(aes(yintercept = real_thresholds$min,
                         color = "Real Min",
                         linetype = "Real"))
      }
      if (max(site_data$mean, na.rm = TRUE) >= real_thresholds_quantiles[2]) {
        plot <- plot +
          geom_hline(aes(yintercept = real_thresholds$min,
                         color = "Real Max",
                         linetype = "Real"))
      }

      # sensor thresholds
      sensor_thresholds_quantiles <- unname(quantile(c(sensor_thresholds$min, sensor_thresholds$max), c(0.1, 0.9)))

      # if (min(site_data$mean, na.rm = TRUE) <= sensor_thresholds_quantiles[1]) {
      #   plot <- plot +
      #     geom_hline(aes(yintercept = sensor_thresholds$mix,
      #                color = "Sensor Mine",
      #                linetype = "Sensor"))
      # }

      if (max(site_data$mean, na.rm = TRUE) >= sensor_thresholds_quantiles[2]) {
        plot <- plot +
          geom_hline(aes(yintercept = sensor_thresholds$max,
                         color = "Sensor Max",
                         linetype = "Sensor"))
      }
      return(plot)
    }
    return(plot)
  }
  return(plot)
}

#' @title Add a flag if the sonde was not fully submerged by water.
#'
#' @description
#' A function designed to append the 'unsubmerged' flag to a row if the value
#' in the `relative_depth` column is less than or equal to 0.
#'
#' @param df A data frame with a `flag` column.
#'
#' @return A data frame with a `flag` column that has been updated with the
#' 'sonde unsubmerged' flag.
#'
#' @examples
#' add_unsubmerged_flag(df = all_data_flagged$`archery-Actual Conductivity`)
#' add_unsubmerged_flag(df = all_data_flagged$`boxelder-Temperature`)
#'
#' @seealso [flag_all_data()]

add_unsubmerged_flag <- function(df){

  # create a df of temperature for each site
  depth <- df %>%
    data.table::data.table() %>%
    dplyr::select(DT_round, DT_join, parameter, mean) %>%
    dplyr::filter(parameter == "Depth") %>%
    dplyr::select(DT_join, Depth = mean)

  # add "depth" column to df:
  depth_checked <- df %>%
    dplyr::left_join(., depth, by = "DT_join") %>%
    # If water temperature is freezing, flag all parameters
    add_flag(., Depth <= 0, "sonde unsubmerged") %>%
    # remove the temp column so df is identical in structure to OG df
    dplyr::select(-Depth)

  return(depth_checked)

}

# Add a verification column to the data frame and automatically fail data points that meet certain criteria.
# @param df A dataframe with a `flag` column.
# @return A dataframe with a `verification` column.
# @examples
# add_verification_column(df = all_data_flagged$`archery-Actual Conductivity`)
# add_verification_column(df = all_data_flagged$`boxelder-Temperature`)

add_verification_column <- function(df) {
  df <- df %>%
    mutate(verification = case_when(
      # are there other conditions where a data point should automatically fail?
      # str_detect(flag, "outside sd range") ~ "fail",
      str_detect(flag, "missing data") ~ "fail",
      # need to add more situations where data points fail
      is.na(flag) ~ "pass",
      TRUE ~ NA))
  return(df)
}

alter_verification_column <- function(verification_choice, df_data) {

  df_name <- paste(unique(df_data$site), unique(df_data$parameter), head(df_data$DT_round, 1))

  # PASS/FAIL statements ----
  # these cases are inverses of each other, both are included
  # in the same if statements to reduce confusion

  # PASS ALL/FAIL NONE
  if (verification_choice %in% c("PASS ALL", "FAIL NONE")) {

    altered_df <- df_data %>%
      mutate(mean_verified = if_else(is.na(flag), mean, NA),
             is_verified = TRUE,
             verification_status = "PASS")

    cat("All points for ", df_name, "have PASSED and are accurate information.\n")

    return(altered_df)

  }

  # PASS VALID/FAIL FLAGGED
  if (verification_choice %in% c("PASS VALID", "FAIL FLAGGED")) {

    altered_df <- df_data %>%
      mutate(mean_verified = mean,
             is_verified = TRUE,
             verification_status = if_else(is.na(flag), "PASS", "FAIL")) # pass valid data, fail flagged data

    cat("All valid points for ", df_name, "have PASSED\n")
    cat("All flagged points for ", df_name, "have FAILED\n")

    return(altered_df)

  }

  # PASS FLAGGED, FAIL VALID
  # if all the flagged data is correct, and none of the valid data is correct,
  # none of it should be shown
  if (verification_choice %in% c("PASS FLAGGED", "FAIL VALID")) {

    altered_df <- df_data %>%
      mutate(mean_verified = NA,
             is_verified = TRUE,
             verification_status = if_else(!is.na(flag), "PASS", "FAIL"))

    cat("All valid points for ", df_name, "have FAILED\n")
    cat("All flagged points for ", df_name, "have PASSED.\n")

    return(altered_df)

  }

  # PASS NONE, FAIL ALL
  # This means that flagged data would actually be accurate data, and valid data would
  # actually be inaccurate data. this is the inverse of pass all
  if (verification_choice %in% c("PASS NONE", "FAIL ALL")) {

    altered_df <- df_data %>%
      mutate(mean_verified = if_else(!is.na(flag), mean, NA),
             is_verified = TRUE,
             verification_status = "FAIL")

    cat("All points for ", df_name, "have FAILED and are inaccurate information.\n")

    return(altered_df)

  }

  # skip statements
  if (verification_choice == "SKIP") {

    altered_df <- df_data %>%
      mutate(mean_verified = if_else(is.na(flag), mean, NA),
             is_verified = TRUE, #false?
             verification_status = "SKIP")

    cat("The accuracy of all points for ", df_name, "cannot be determined as\n")
    cat("accurate or inaccurate information at this time.")

    return(altered_df)

  }

  # inspect statements
  if (verification_choice %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME")) {
    stop("You can't inspect here")
  }

  # quit statements ***
  if (verification_choice == "QUIT") {

    cat("Quitting out of verification mode...")
    cat("\n")

    QUIT <<- TRUE

    return(NULL)

  }

}


# to do (j): document this function
api_puller <- function(site, start_dt, end_dt = Sys.time(), api_token, dump_dir, require = NULL) {

  locs <- hv_locations_all(hv_token)

  # make a list of site names

  options(scipen = 999)
  for(i in 1:length(site)){

    site_loc <- locs %>%
      dplyr::mutate(name = tolower(name)) %>%
      dplyr::filter(grepl(site[i], name, ignore.case = TRUE))

    site_loc_list <- site_loc$id

    # Get data for each site location. Note this maps over the entire list of locations,
    # many of which are unlikely to be active during the time you specify. Don't freak out if
    # you see a bunch of '404 Not Found' errors, you're just seeing the list of locations
    # that are not active. The data frame 'alldata' should contain your data from all applicable
    # sites during the time frame indicated. Note that this will take some time (one month of
    # data for 5 sites takes ~10 mins. Be patient!

    # Add date range you are interested in; data are stored in HydroVu in UTC
    # Here, a way to find the most recent download of the data. Use this as the start date to
    # reduce overlapping data

    # tz weirdness
    # utc_start_date <- format(as.POSIXct(start_dt, tz = "UTC") + hours(7), format = "%Y-%m-%d %H:%M:%S")
    #
    # utc_end_date <-   format(as.POSIXct(end_dt, tz = "UTC") + hours(7), format = "%Y-%m-%d %H:%M:%S")

    # doing this fixes the mismatch in date times during the combined_data step - jd
    utc_start_date <- format(as.POSIXct(start_dt, tz = "UTC") + lubridate::hours(0), format = "%Y-%m-%d %H:%M:%S")

    utc_end_date <-   format(as.POSIXct(end_dt, tz = "UTC") + lubridate::hours(0), format = "%Y-%m-%d %H:%M:%S")

    timezone = "UTC"

    # Map over the location ids
    alldata <- site_loc_list %>% purrr::map(~hv_data_id(.,
                                                        start_time = utc_start_date,
                                                        end_time = utc_end_date,
                                                        token = api_token,
                                                        tz = timezone))

    # grab only locations with data (stored as a data frame) / drop 404 errors
    filtered <- purrr::keep(alldata, is.data.frame)

    if(length(filtered) == 0){

      print(paste0("No data at ", site[i], " during this time frame"))

    } else {

      # bind lists together (now that all are dataframes, we can just collate quickly)
      one_df <- dplyr::bind_rows(filtered) %>%
        dplyr::rename(id = Location,
                      parameter = Parameter,
                      units = Units) %>%
        dplyr::left_join(., site_loc, by = "id") %>%
        dplyr::mutate(site = tolower(site[i])) %>%
        dplyr::select(site, id, name, timestamp, parameter, value, units)

      # if site contains both a csu and a virridy sonde, split them up:

      if(site[i] %in% c("Timberline", "Prospect", "Archery")){

        try(virridy_df <- one_df %>%
              filter(grepl("virridy", name, ignore.case = TRUE)) %>%
              mutate(site = paste0(site[i], " virridy")))

        try(readr::write_csv(virridy_df,
                             paste0(dump_dir, "/", site[i], " virridy_", stringr::str_replace(stringr::str_replace(substr(end_dt, 1, 16), "[:\\s]", "_"), ":", ""), ".csv")))

        csu_df <- one_df %>%
          filter(!grepl("virridy", name, ignore.case = TRUE))

        readr::write_csv(csu_df,
                         paste0(dump_dir, "/", site[i], "_", stringr::str_replace(stringr::str_replace(substr(end_dt, 1, 16), "[:\\s]", "_"), ":", ""), ".csv"))

      } else {

        # otherwise, save the full data set

        readr::write_csv(one_df,
                         paste0(dump_dir, "/", site[i], "_", stringr::str_replace(stringr::str_replace(substr(end_dt, 1, 16), "[:\\s]", "_"), ":", ""), ".csv"))
      }
    }
  }

}

append_historical_api_data <- function(hist_dir, inc_dir) {

  # list of site names
  site_names <- list(
    "Archery",
    "Boxelder",
    "Legacy",
    "Lincoln",
    "Prospect",
    "River Bluffs",
    "Tamasag",
    "Timberline"
  )

  # get the full file names
  inc_full_file_names <- list.files(inc_dir, pattern = "*.csv", full.names = TRUE)
  hist_full_file_names <- list.files(hist_dir, pattern = "*.csv", full.names = TRUE)

  # find the files that match using the site names
  walk(site_names, function(site_name) {

    # find the index of the matching site names in the file names
    inc_site_name_full_path <- grepl(site_name, inc_full_file_names, ignore.case = TRUE)
    hist_site_name_full_path <- grepl(site_name, hist_full_file_names, ignore.case = TRUE)

    # get the file names
    inc_site_path <- inc_full_file_names[inc_site_name_full_path]
    hist_site_path <- hist_full_file_names[hist_site_name_full_path]

    # read in the files
    inc_file <- read_csv(inc_site_path)
    hist_file <- read_csv(hist_site_path)

    # combine the files
    new_file <- bind_rows(hist_file, inc_file) %>%
      distinct()

    # write the file (this will replace the old file)
    write_csv(new_file, paste0(hist_dir, site_name, "_historical.csv"))
  })

  # delete the files in the inc_dir
  walk(inc_full_file_names, unlink)

}




basic_plot <- function(){
  max_1 <- max(wq_tl[[parameters[1]]], na.rm = T)

 param_1 <- ggplot(wq_tl, aes(x = DT_round, y = .data[[parameters[1]]], ymin = 0, ymax = .data[[parameters[1]]]))+
   geom_ribbon(color = "blue", fill = "blue" )+
   theme_bw()+
   labs(x = "Date", y = parameters[1])

 param_2 <- ggplot(wq_tl, aes(x = DT_round, y = .data[[parameters[2]]]))+
   geom_path(color = "red")+
   theme_bw()+
   labs(x = "Date", y = parameters[2])

 ggarrange(param_1, param_2,ncol = 2, nrow = 1)

}

check_incoming_api_dir <- function(incoming_dir, archive_dir) {
  # Check if data/api/incoming_api_data exists
  if (dir.exists(incoming_dir)) {
    # Check if incoming directory is empty
    if (length(list.files(incoming_dir)) == 0) {
      print(paste0(incoming_dir, " exists and is empty."))
      print("Incoming API data directory is configured properly.")
    } else {
      print(paste0(incoming_dir, " exists but is not empty..."))
      print("Please ensure previous workflow ran properly...")
      stop("Pipeline halted due to non-empty incoming directory.")
    }
  } else {
    print(paste0(incoming_dir, " does not exist..."))
    print("Creating incoming directory...")
    dir.create(incoming_dir, recursive = TRUE)
    print("Incoming directory created.")
  }

  # Check if data/api/archive_api_data exists
  if (dir.exists(archive_dir)) {
    print(paste0(archive_dir, " exists."))
    print("Directory is configured properly.")
  } else {
    print(paste0(archive_dir, " does not exist..."))
    print("Creating archive directory...")
    dir.create(archive_dir, recursive = TRUE)
    print("Archive directory created.")
  }
}

clean_directories <- function() {
  pre_dir_path <- pre_verification_path
  int_dir_path <- intermediary_path
  ver_dir_path <- verified_path

  pre_dir_names <- list.files(pre_dir_path)
  int_dir_names <- list.files(int_dir_path)
  ver_dir_names <- list.files(ver_dir_path)

  # Pre-verification Directory
  # keep if not in Verified Directory, else if in Verified Directory delete from this Directory
  for (i in pre_dir_names) {
    if (length(ver_dir_names) != 0 & i %in% ver_dir_names) {
      file.remove(paste0(pre_dir_path, i))
      cat("removed file ", i, " from ",  pre_dir_path, "\n")
    } else {
      next
    }
  }

  # Intermediary Directory
  # keep if skips in the data OR any(!is_verified), else move to Verified Directory
  for (i in int_dir_names) {
    df <- readRDS(paste0(int_dir_path, i))
    if (any(df$verification_status == 'SKIP') | any(!df$is_verified)) {
      next
    } else if (all(df$verification_status != 'SKIP') & all(df$is_verified)) {
      source_path <- paste0(int_dir_path, i)
      destination_path <- paste0(ver_dir_path, i)
      file.copy(source_path, destination_path, overwrite = TRUE)
      file.remove(source_path)
      cat("moved file ", i, " from ", int_dir_path, "to", ver_dir_path, "\n")
    }
  }

  # Verified Directory
  # final destination for data, check that the data that is here should be here
  for (i in ver_dir_names) {
    df <- readRDS(paste0(ver_dir_path, i))
    if (any(df$verification_status == 'SKIP') | any(!df$is_verified)) {
      cat(i, "should not be in verified directory\n")
      stop("THERE IS DATA WITH SKIPS OR NON-VERIFIED DATA IN THE VERIFIED DIRECTORY")
    } else {
      cat("data in verified directory is correct\n")
      next
    }
  }
}

clear_incoming_data_dir <- function(incoming_dir, archive_dir, require = NULL) {
  # Check if the previous step in the targets pipeline ran
  # this is done via the require arg

  # List files in the incoming and archive directories
  incoming_files <- list.files(incoming_dir, full.names = FALSE)
  archive_files <- list.files(archive_dir, full.names = FALSE)

  # Find files that are not already in the archive directory
  files_to_copy <- setdiff(incoming_files, archive_files)

  # Copy only the files that are not already in the archive directory
  if (length(files_to_copy) > 0) {
    for (file in files_to_copy) {
      full_file_name <- file.path(incoming_dir, file)
      file.copy(full_file_name, archive_dir)
      print(paste0(file, " has been moved to archive API data folder."))
    }
    print("Files have been copied from the incoming directory to the archive directory.")
  } else {
    print("All files are already present in the archive directory. Nothing to copy.")
  }
Sys.sleep(5)

# refresh archive_files now that the incoming_files have been copied over
archive_files <- list.files(archive_dir, full.names = FALSE)
  # Check if all files from incoming directory are now in the archive directory
  if (all(incoming_files %in% archive_files)) {
    print("All files in the incoming directory have been successfully copied to the archive directory.")
  } else {
    print("Not all files from the incoming directory have been successfully copied to the archive directory.")
    # Should this halt the pipeline?
    # Right now it seems to always print this and I haven't figured out why so
    # I don't think that it should not halt the pipeline in the state that its in right now -JD.
  }

  # Delete the copied files from the incoming directory
  if (length(files_to_copy) > 0) {
    for (file in files_to_copy) {
      full_file_name <- file.path(incoming_dir, file)
      file.remove(full_file_name)
    }
    print("Copied files have been removed from the incoming directory.")
  }

  # Delete any files in the incoming directory
  for (file in list.files(incoming_dir, full.names = TRUE)) {
    file.remove(file)
  }
  print("All files removed from incoming directory.")
}

#' @title Combine historical and incoming data
#'
#' @description
#' A function that combines the subset of historically flagged data and creates historical columns.
#' NOTE: historical columns should be a part of the historically flagged data, but it is not at this point.
#'
#' @param incoming_data_list A list of dataframes with the incoming data.
#'
#' @param historical_data_list A list of dataframes with the historical data.
#'
#' @return A list of dataframes with the combined historical and incoming data.
#'
#' @examples
#' combine_hist_inc_data(incoming_data_list = incoming_data_list, historical_data_list = historical_data_list)

combine_hist_inc_data <- function(incoming_data_list, historical_data_list) {

  # Get the matching index names
  matching_indexes <- intersect(names(incoming_data_list), names(historical_data_list))

  # Get the last 24 hours of the historically flagged data based on earliest
  # DT listed in the api pull for that site/data combo
  last_24_hours <- map(historical_data_list,
                      function(data) {
                        data %>%
                          # most recent day of already-flagged data
                          filter(DT_round >= ymd_hms(max(DT_round) - hours(24), tz = "MST")) %>%
                          # mark it as "historic"
                          mutate(historical = TRUE,
                                 # ensure last_site_visit column has proper DT:
                                 last_site_visit = force_tz(last_site_visit, tzone = "MST"))
                      })

  # bind summarized_incoming_data and last_24_hours together
  combined_hist_inc_data <- map(matching_indexes, function(index) {
    last_24_hours[[index]] %>%
      select(-c(t_mean01, # why is this being removed? -jd
                t_mean99,
                t_slope_behind_01,
                t_slope_behind_99,
                t_sd_0199)) %>%
      bind_rows(., mutate(incoming_data_list[[index]], historical = FALSE)[-1,]) %>% # why are we removing the first row? -jd # there is a datetime issue on this bind
      # make sure new data has info about last site visit, etc. filled out
      fill(c(sonde_employed, last_site_visit, sensor_malfunction))

  }) %>%
    set_names(matching_indexes)

  return(combined_hist_inc_data)
}



#@param folder_name: input the raw photo folder that is named the site you wish to compile
compile_files <- function(folder_name){
  site <- str_extract(folder_name, "tamasag|legacy|timberline|prospect|boxelder|archery|riverbluffs")
  photo_files <- list.files(path = folder_name, full.names = TRUE, recursive = TRUE)
  #files in folder
  new_folder_files <- list.files(path = paste0('data/timelapse_photos/2023_compiled/',site), full.names = TRUE)

  photo_renamer <- function(file) {
    #grab dt from file
    dt <- read_exif(path = file,
                    tags = c("DateTimeOriginal")) %>%
      pull(DateTimeOriginal) %>%
      parse_date_time("YmdHMS", tz="MST")%>%
      format("%Y%m%d_%H%M")
    #create new file name from dt
    new_file_name <- paste0('data/timelapse_photos/2023_compiled/',site,"/", dt,'.JPG')

    #check to see if this file is already in the folder
    if(new_file_name %nin% new_folder_files){
      #if it is not, copy it over
      file.copy(file,to =  new_file_name)
    }

  }
  map(photo_files, photo_renamer)
  print(paste0("Finished ", folder_name))

}

## Photos

# Goals:
# Download all user created photos ( upstream, downstream, clarity, filter and other pictures)
# Label according to site, date, description in the format site_YYYYMMDD_descriptor.jpg
# Only download photos which have not yet been downloaded




download_pictures <- function(){
  #source to grab all notes cleaned
  source("src/load_mWater_notes.R")

  all_notes_cleaned <- load_mWater_notes()
  # Find all the downloaded pictures
  all_file_names <- tolower(list.files(path = "data/field_pics/", recursive = TRUE))
  #grab notes
  sampling_photos <- all_notes_cleaned%>%
    #grab needed columns
    select(site, start_dt,photos_downloaded,upstream_pic,downstream_pic,clarity,filter_pic,other_pic,other_pic_descriptor)%>%
    mutate(
      #correct names if it is in our upper sites (upper case acronyms)
      site = tolower(site),
      #Date format for pictures
      yyyymmdd = format(start_dt, "%Y%m%d"),
      #create filenames ONLY if there is a URL associated with the site visit
      upstream_filename = case_when(
        !is.na(upstream_pic) ~ paste0(site, "_", yyyymmdd, "_upstream.jpg"),
        TRUE ~ NA_character_
      ),
      downstream_filename = case_when(
        !is.na(downstream_pic) ~ paste0(site, "_", yyyymmdd, "_downstream.jpg"),
        TRUE ~ NA_character_
      ),
      clarity_filename = case_when(
        !is.na(clarity) ~ paste0(site, "_", yyyymmdd, "_clarity.jpg"),
        TRUE ~ NA_character_
      ),
      filter_filename = case_when(
        !is.na(filter_pic) ~ paste0(site, "_", yyyymmdd, "_filter.jpg"),
        TRUE ~ NA_character_
      ),
      # check to see if the photo is already downloaded to folder
      upstream_downloaded = case_when(
        is.na(upstream_filename) ~ NA,
        upstream_filename %in% all_file_names ~ TRUE,
        TRUE ~ FALSE
      ),
      clarity_downloaded = case_when(
        is.na(clarity_filename) ~ NA,
        clarity_filename %in% all_file_names ~ TRUE,
        TRUE ~ FALSE
      ),
      filter_downloaded = case_when(
        is.na(filter_filename) ~ NA,
        filter_filename %in% all_file_names ~ TRUE,
        TRUE ~ FALSE
      ),
      downstream_downloaded = case_when(
        is.na(downstream_filename) ~ NA,
        downstream_filename %in% all_file_names ~ TRUE,
        TRUE ~ FALSE
      )
    )
# basic path to field pics
path <- "data/field_pics/"


  # loop thru dataset and download the photo ONLY if it is not yet downloaded and not NA
  for (i in 1:nrow(sampling_photos)) {
    if (!is.na(sampling_photos$upstream_downloaded[i]) && !sampling_photos$upstream_downloaded[i]) {
      print(sampling_photos$upstream_filename[i])
      download.file(sampling_photos$upstream_pic[i], destfile = paste0(path,sampling_photos$upstream_filename[i]))
    }

    if (!is.na(sampling_photos$downstream_downloaded[i]) && !sampling_photos$downstream_downloaded[i]) {
      print(sampling_photos$downstream_filename[i])
      download.file(sampling_photos$downstream_pic[i], destfile = paste0(path, sampling_photos$downstream_filename[i]))
    }
    if (!is.na(sampling_photos$clarity_downloaded[i]) && !sampling_photos$clarity_downloaded[i]) {
      print(sampling_photos$clarity_filename[i])
      download.file(sampling_photos$clarity[i], destfile = paste0(path, sampling_photos$clarity_filename[i]))
    }
    if (!is.na(sampling_photos$filter_downloaded[i]) && !sampling_photos$filter_downloaded[i]) {
      print(sampling_photos$filter_filename[i])
      download.file(sampling_photos$filter_pic[i], destfile = paste0(path, sampling_photos$filter_filename[i]))
    }
  }

  #grab notes for sites with other pictures
  other_photos <- all_notes_cleaned%>%
    #grab needed columns
    select(site, start_dt,other_pic,other_pic_descriptor)%>%
    #get rid of instances with no other pics
    filter(!is.na(other_pic))%>%
    mutate(
      site = tolower(site),
      #Date format for pictures
      yyyymmdd = format(start_dt, "%Y%m%d"),
      # separate multiple URLs in other pic column
      other_pic_sep = str_split(other_pic, "; "),
      #seperate multiple descriptors in the descriptor column
      other_descriptor_sep = str_split(other_pic_descriptor, ","))%>%
    #for rows with multiple pictures, create a new row for each picture

    unnest(cols = c(other_pic_sep, other_descriptor_sep))%>%
    #remove excess columns and rename sep columns to match old columns
    select(site, start_dt,yyyymmdd, other_pic = other_pic_sep, other_pic_descriptor = other_descriptor_sep)%>%
    # make descriptor lower case and remove any spaces in the name
    mutate(other_pic_descriptor = tolower(str_replace_all(other_pic_descriptor, " ", "")),
           other_filename = case_when(!is.na(other_pic) ~ paste0(site, "_", yyyymmdd, "_", other_pic_descriptor, ".jpg")),
           # Check to see if photo has already been downloaded
           other_downloaded = case_when(
             is.na(other_filename) ~ NA,
             other_filename %in% all_file_names ~ TRUE,
             TRUE ~ FALSE
           ))

  # loop thru dataset and download the photo ONLY if it is not yet downloaded and not NA
  for (i in 1:nrow(other_photos)) {
    if (!is.na(other_photos$other_downloaded[i]) && !other_photos$other_downloaded[i]) {
      print(other_photos$other_filename[i])
      download.file(other_photos$other_pic[i], destfile = paste0(path,other_photos$other_filename[i]))
    }}

cat("\nAll Available Pictures Downloaded\n")
}

#download_pictures()
#RUN TO DOWNLOAD NEW PICTURES
# It takes about 2-5 minutes to download ~25-50 photos
# Sometimes the request to mWater time out, just re run the function below if that happens

## Determining uploads

#This function looks at the user inputs for calibration report collect and logs collected.
#Based on these inputs, it looks at all the uploaded logs or calibration reports
#and will print out what logs are missing and who to contact to get those files uploaded.




files_missing <- function(){

  `%nin%` = Negate(`%in%`)
  # #source clean mwater script for all notes cleaned
  #
  # source("src/mWater_collate/load_mWater_notes.R")

  #grab context metadata
  site_meta <- read_csv("data/metadata/water_sampling_sites.csv",show_col_types = FALSE)%>%
    select(site = site_code, Site_Name, site_label_rmrs)
  # sort for sites in upper network (ie. acronyms rather than street names)
  upper_sites <- read_csv("data/metadata/water_sampling_sites.csv",show_col_types = FALSE)%>%
    filter(watershed != "CLP  Mainstem-Fort Collins")%>%
    #this is to help match with user input
    mutate(site_code = tolower(site_code))

  field_season <- year(Sys.Date())

  #grab cal reports from folder
  cal_reports_simple <- str_extract(list.files(path = "data/calibration_reports/"), ".*_\\d{8}" )%>%tolower()
  logs_simple <- str_extract(list.files(path = paste0("data/sensor_data/", field_season), recursive = TRUE), "\\w+_\\d{8}_(vulink|troll)")%>%tolower()




  #grab sensor notes that have logs or cal reports that should be  downloaded
  sensor_files <- all_notes_cleaned%>%
    filter(year(DT_round) == field_season)%>%
    filter(grepl("Sensor",visit_type, ignore.case = TRUE))%>%
    filter(cal_report_collected|log_downloaded)%>%
    select(site, crew, start_DT,end_dt, cal_report_collected, cals_performed, log_downloaded, log1_type,log1_mmdd,  log2_type, log2_mmdd)%>%
    mutate(
      #make all site names lower
      site = tolower(site),
      # Create basis for calibration report name
      # this will be used to check for calibration report in data files and then b
      cal_report_name = case_when(cal_report_collected == TRUE ~ paste0(site, "_", format(start_DT, "%Y%m%d")),
                                  cal_report_collected == NA ~ NA),
      full_cal_name = case_when(cal_report_collected == TRUE ~ paste0(site, "_", format(end_dt, "%Y%m%d_%H%M_mst")),
                                cal_report_collected == NA ~ NA),
      log1_mmdd = case_when(nchar(as.character(log1_mmdd)) == 3 ~ paste0("0",log1_mmdd),
                            TRUE ~ as.character(log1_mmdd)),
      log1_type = case_when( grepl("aquatroll", log1_type,ignore.case = TRUE) ~  "troll",
                             grepl("at", log1_type,ignore.case = TRUE) ~  "troll",
                             grepl("aqua troll", log1_type,ignore.case = TRUE) ~  "troll",
                             grepl("vulink", log1_type,ignore.case = TRUE) ~  "vulink",
                             grepl("vu link", log1_type,ignore.case = TRUE) ~  "vulink",
                             TRUE ~ tolower(log1_type)),
      log2_type = case_when( grepl("aquatroll", log2_type,ignore.case = TRUE) ~  "troll",
                             grepl("at", log2_type,ignore.case = TRUE) ~  "troll",
                             grepl("aqua troll", log2_type,ignore.case = TRUE) ~  "troll",
                             grepl("vulink", log2_type,ignore.case = TRUE) ~  "vulink",
                             grepl("vu link", log2_type,ignore.case = TRUE) ~  "vulink",
                             TRUE ~ tolower(log2_type)),
      log2_mmdd = case_when(nchar(as.character(log2_mmdd)) == 3 ~ paste0("0",log2_mmdd),
                            TRUE ~ as.character(log2_mmdd)),
      log1_user_error = case_when( is.na(log_downloaded)~ FALSE,
                                   log_downloaded == FALSE ~ FALSE,
                                   is.na(log1_mmdd) | is.na(log1_type) ~ TRUE,
                                   TRUE ~ FALSE ),
      log2_user_error = case_when(is.na(log_downloaded)~ FALSE,
                                  log_downloaded == FALSE ~ FALSE,
                                  is.na(log2_mmdd) & is.na(log2_type) ~ FALSE,
                                  is.na(log2_mmdd) | is.na(log2_type) ~ TRUE,
                                  TRUE ~ FALSE),
      log1_name = case_when(!(is.na(log1_mmdd) | is.na(log1_type)) ~ paste0(site,"_",year(start_DT),log1_mmdd,
                                                                            "_", format(start_DT, "%Y%m%d"), "_", log1_type),
                            (is.na(log1_mmdd) | is.na(log1_type))  ~ NA),
      log2_name = case_when(!(is.na(log2_mmdd) | is.na(log2_type)) ~ paste0(site,"_",year(start_DT),log2_mmdd,
                                                                            "_", format(start_DT, "%Y%m%d"), "_", log2_type),
                            (is.na(log2_mmdd) | is.na(log2_type))  ~ NA),
      log_missing1 = case_when(
        is.na(log_downloaded)| log_downloaded == FALSE ~ FALSE,
        is.na(log1_name) ~ FALSE,
        log1_name %nin% logs_simple ~ TRUE,
        TRUE ~ FALSE
      ),
      log_missing2 = case_when(
        is.na(log_downloaded)| log_downloaded == FALSE ~ FALSE,
        is.na(log2_name) ~ FALSE,
        log2_name %nin% logs_simple ~ TRUE,
        TRUE ~ FALSE
      ),
      log_missing = case_when(
        log_missing1 | log_missing2 ~ TRUE,
        TRUE ~ FALSE
      ),
      cal_missing = case_when(
        is.na(cal_report_collected) ~ FALSE,
        cal_report_name %nin% cal_reports_simple ~ TRUE,
        TRUE ~ FALSE)
      )



  for (i in 1:nrow(sensor_files)) {

    #if log missing print out missing logs or cal reports
    if(sensor_files$log_missing[i]){
      cat("\nLog Missing: ", sensor_files$log1_name[i], " and/or ", sensor_files$log2_name[i], "\nContact: ", sensor_files$crew[i], "\n")

    }
    #if log missing print out missing logs or cal reports
    if(sensor_files$cal_missing[i]){
      cat("\nCal Missing: ", sensor_files$full_cal_name[i]," \nContact: ", sensor_files$crew[i], "\n")
    }


  }
  cat("\nFile Check Complete")

}


#' @title Add malfunction flags to a data frame based on field notes and HydroVu.
#'
#' @description
#' This function adds the 'sensor malfunction' flag to a data frame by the dates,
#' sites, and parameters that members of the lab know contain erroneous data due
#' to sensor malfunctions. Note that this flag is used in two instances: first
#' when removing erroneous data from our statistics calculations and again during
#' the actual flagging step in the QAQC pipeline.
#'
#' @param df A data frame with a `flag` column.
#'
#' @param malfunction_records The malfunction records pulled from mWater
#'
#' @return A data frame with a `flag` column that has been updated with the
#' 'sensor malfunction' flag.
#'
#' @examples
#' add_range_flags(df = all_data_flagged$`archery-Actual Conductivity`)
#' add_range_flags(df = all_data_flagged$`boxelder-Temperature`)
#'
#' @seealso [grab_mWater_malfunction_notes()]

fix_calibration <- function(df){

  # Filter records for relevant site-param information
  df_site <- unique(df$site)
  df_parameter <- unique(df$parameter)

  # Depth calibration requires its own fix. It's hacky and I don't like it, but
  # the way that the calibration report reports pressure calibrations makes it
  # impossible to actually back calibrate. Therefore, I'm just assuming that level
  # doesn't actually change after bad calibrations, and hard code it so that the
  # first "bad" depth is forced to equal the last "good" depth and I offset all
  # subsequent "bad" depth values by the difference between them.

  if(df_parameter == "Depth"){

    if(!"Depth" %in% df$parameter){

      nope <- df %>% mutate(relative_depth = NA)

      return(nope)

    }

    if("Depth" %in% df$parameter & "archery" %in% df$site){

      #df <- all_data_flagged[["archery-Depth"]]

      site_depth <- df %>%
        mutate(relative_depth = ifelse(DT_round >= as_datetime('2022-05-21 15:45:00', "MST") & DT_round <= as_datetime('2022-05-24 15:45:00', "MST"),
                                       mean +
                                         abs(
                                           dplyr::filter(df, as.character(DT_round) == ('2022-05-24 15:15:00'))$mean -
                                             dplyr::filter(df, as.character(DT_round) == ('2022-05-24 15:45:00'))$mean
                                         ),
                                       mean))
      #return(archery_depth)

    } else if("Depth" %in% df$parameter & "timberline" %in% df$site){

      # df <- all_data_flagged[["timberline-Depth"]]

      site_depth <- df %>%
        dplyr::mutate(relative_depth = ifelse(year == "2022" & DT_round <= lubridate::as_datetime('2022-04-07 17:00:00', "MST"),
                                              as.numeric(mean) +
                                                abs(
                                                  dplyr::filter(df, as.character(DT_round) == ('2022-04-07 16:15:00'))$mean -
                                                    dplyr::filter(df, as.character(DT_round) == ('2022-04-07 17:15:00'))$mean
                                                ),
                                              as.numeric(mean)))
      # return(timberline_depth)

    } else if ("Depth" %in% df$parameter & "legacy" %in% df$site) {

      # df <- all_data_flagged[["legacy-Depth"]]

      site_depth <- df %>%
        mutate(relative_depth = ifelse(DT_round >= as_datetime('2022-04-06 06:00:00', "MST") & DT_round <= as_datetime('2022-04-12 09:15:00', "MST"),
                                       as.numeric(mean) +
                                         abs(
                                           dplyr::filter(df, as.character(DT_round) == ('2022-04-12 09:15:00'))$mean -
                                             dplyr::filter(df, as.character(DT_round) == ('2022-04-12 09:30:00'))$mean
                                         ),
                                       mean))

      site_depth <-  site_depth %>%
        mutate(relative_depth = ifelse(DT_round >= as_datetime('2022-07-08 17:00:00', "MST") & DT_round <= as_datetime('2022-07-12 09:00:00', "MST"),
                                       relative_depth +
                                         abs(
                                           dplyr::filter(site_depth, as.character(DT_round) == ('2022-07-08 14:15:00'))$relative_depth -
                                             dplyr::filter(site_depth, as.character(DT_round) == ('2022-07-08 17:00:00'))$relative_depth
                                         ),
                                       relative_depth))

      site_depth <-  site_depth %>%
        mutate(relative_depth = ifelse(DT_round >= as_datetime('2022-07-22 11:30:00', "MST") & DT_round <= as_datetime('2022-07-25 14:15:00', "MST"),
                                       relative_depth +
                                         abs(
                                           dplyr::filter(site_depth, as.character(DT_round) == ('2022-07-22 09:45:00'))$relative_depth -
                                             dplyr::filter(site_depth, as.character(DT_round) == ('2022-07-22 11:30:00'))$relative_depth
                                         ),
                                       relative_depth))



      #return(legacy_depth)

    } else if ("Depth" %in% df$parameter & "tamasag" %in% df$site) {

      # df <- all_data_flagged[["tamasag-Depth"]]

      site_depth <- df %>%
        mutate(relative_depth = ifelse(DT_round <= "2022-04-24 07:15:00" & year == "2022",
                                       as.numeric(mean) +
                                         abs(
                                           dplyr::filter(df, as.character(DT_round) == ('2022-04-24 07:15:00'))$mean -
                                             dplyr::filter(df, as.character(DT_round) == ('2022-04-24 07:30:00'))$mean
                                         ),
                                       mean))

      # return(tamasag_depth)

    } else {

      site_depth <- df %>%
        mutate(relative_depth = mean,
               cal_fix = NA)
    }

    depth_flagged <- site_depth %>%
      dplyr::mutate(raw = mean,
                    mean = relative_depth,
                    cal_fix = ifelse(raw != mean, "calibration fix", NA)) %>%
      dplyr::select(-relative_depth)

    return(depth_flagged)

  }

  # PLACEHOLDER UNTIL OTHER CALIBRATIONS ARE GOOD:
 return(df %>%
    mutate(raw = mean,
           cal_fix = NA))

}

  # # For non-depth parameters, we can refer to the calibration report:
  # bad_cal <- readxl::read_excel("data/calibration_error_log.xlsx") %>%
  #   dplyr::mutate(start_DT = as.character(lubridate::as_datetime(start_DT)),
  #                 end_DT = as.character(lubridate::as_datetime(end_DT)),
  #                 #report_to_correct = as.character(lubridate::as_datetime(report_to_correct))
  #   )
  #
  # bad_cal_records_filtered <- bad_cal %>%
  #   filter(site == df_site) %>%
  #   filter(grepl(df_parameter, parameter, ignore.case = TRUE)) %>%
  #   mutate(end_DT = ifelse(is.na(end_DT), ymd_hms("9999-12-31 23:59:59", tz = "MST"), end_DT)) %>%
  #   mutate(end_DT = as.character(as.POSIXct(end_DT, tz = "MST"))) %>%
  #   rowid_to_column()
  #
  # # If there are no bad calibrations listed for that site-param, return original
  # # dataframe, filling our updated "mean" column with the old unmodified values:
  # if(nrow(bad_cal_records_filtered == 0)){
  #
  #   df <- df %>%
  #     mutate(raw = mean,
  #            mean = mean)
  #
  #   return(df)
  #
  # } else {
  #
  #   cal_tabler <- function(cal_files){
  #
  #     #cal_files <- list.files("data/calibration_reports")[3] # for testing
  #
  #     cal <- read_html(paste0(getwd(), "/data/calibration_reports/", cal_files)) %>%
  #       html_nodes("div") %>%
  #       html_text() %>%
  #       as_tibble()
  #
  #     rdo <- cal %>% filter(grepl("RDO", value)) %>% pull() %>% str_replace_all(., " ", "") %>% tolower()
  #
  #     ph_orp <- cal %>% filter(grepl("pH/ORP", value)) %>% pull() %>% str_replace_all(., " ", "") %>% tolower()
  #
  #     conductivity <- cal %>% filter(grepl("Conductivity",value)) %>% pull() %>% str_replace_all(., " ", "") %>% tolower()
  #
  #     if(length(cal %>% filter(grepl("Turbidity",value)) %>% pull() %>% str_replace_all(., " ", "") %>% tolower()) != 0){
  #
  #       turbidity <- cal %>% filter(grepl("Turbidity",value)) %>% pull() %>% str_replace_all(., " ", "") %>% tolower()
  #
  #     } else {
  #
  #       turbidity <- "No Turbidity Sensor"
  #
  #     }
  #
  #     # Always the last sensor when depth is available:
  #     depth <- ifelse(str_detect(cal %>% .[nrow(.),] %>% pull() %>% str_replace_all(., " ", "") %>% tolower(), "pressure"),#"psireferencedepth"),
  #                     cal %>% .[nrow(.),] %>% pull() %>% str_replace_all(., " ", "") %>% tolower(),
  #                     "No Depth Sensor")
  #
  #     time_mst <- paste0(str_sub(str_match(cal_files, "(\\d+)_mst")[, 2:1][1], 1, 2), ":",
  #                        str_sub(str_match(cal_files, "(\\d+)_mst")[, 2:1][1], 3, 4))
  #     #str_sub(cal_files, -13, -12),":", str_sub(cal_files, -11, -10))
  #
  #     date <- str_match(cal_files, "^[^_]+_([0-9]{8})_")[, 2]
  #
  #     #str_sub(cal_files, -22, -19),"-", str_sub(cal_files, -18, -17),"-", str_sub(cal_files, -16, -15))
  #
  #     cal_table <- tibble(site = sub("\\_.*", "", cal_files) %>% tolower(),
  #
  #                         DT = ymd_hm(paste(date, time_mst, tz = "MST")),
  #
  #                         # Dissolved Oxygen
  #                         rdo_cal_date = as.character(mdy(str_match(rdo, "lastcalibrated\\s*(.*?)\\s*calibrationdetails")[,2])),
  #                         rdo_slope = str_match(rdo, "slope\\s*(.*?)\\s*offset")[,2],
  #                         rdo_offset = str_match(rdo, "offset\\s*(.*?)\\s*mg/l")[,2],
  #                         rdo_100 = str_match(rdo, "premeasurement\\s*(.*?)\\s*%satpost")[,2],
  #                         rdo_conc = str_match(rdo, "concentration\\s*(.*?)\\s*mg/lpremeasurement")[,2],
  #                         rdo_temp = str_match(rdo, "temperature\\s*(.*?)\\s*°c")[,2],
  #                         rdo_pressure = str_match(rdo, "pressure\\s*(.*?)\\s*mbar")[,2],
  #
  #                         # pH
  #                         ph_cal_date = as.character(mdy(str_match(ph_orp, "lastcalibrated\\s*(.*?)\\s*calibrationdetails")[,2])),
  #                         ph_slope_pre = str_match(ph_orp, "offset1slope\\s*(.*?)\\s*mv/ph")[,2],
  #                         ph_offset_pre = str_match(ph_orp, "mv/phoffset\\s*(.*?)\\s*mvslopeandoffset2")[,2],
  #                         ph_slope_post = str_match(ph_orp, "offset2slope\\s*(.*?)\\s*mv/ph")[,2],
  #                         ph_offset_post = str_match(ph_orp, paste0(ph_slope_post,"mv/phoffset\\s*(.*?)\\s*mvorporp"))[,2],
  #                         # Sometimes, the post value can actually be in the high 6 pH... therefore the post measurement regex matching text is conditional
  #                         ph_7_nice = str_sub(str_match(ph_orp, "postmeasurementph7\\s*(.*?)\\s*mvcal")[,2], 10, nchar(str_match(ph_orp, "postmeasurementph7\\s*(.*?)\\s*mvcal")[,2])),
  #                         ph_7_high = str_sub(str_match(ph_orp, "postmeasurementph8\\s*(.*?)\\s*mvcal")[,2], 10, nchar(str_match(ph_orp, "postmeasurementph8\\s*(.*?)\\s*mvcal")[,2])),
  #                         ph_7_low = str_sub(str_match(ph_orp, "postmeasurementph6\\s*(.*?)\\s*mvcal")[,2], 10, nchar(str_match(ph_orp, "postmeasurementph6\\s*(.*?)\\s*mvcal")[,2])),
  #                         ph_7 = ifelse(!is.na(ph_7_nice), ph_7_nice,
  #                                       ifelse(!is.na(ph_7_high), ph_7_high, ph_7_low)),
  #
  #                         # ORP
  #                         #Newly encountered thing: sometimes the calibration report calls the ORP standard Zobell's, sometimes it's just called "ORP Standard":
  #                         orp_offset = ifelse(is.na(str_match(ph_orp, "zobell'soffset\\s*(.*?)\\s*mvtemperature")[,2]) & is.na(str_match(ph_orp, "quickcal\\s*(.*?)\\s*mvtemperature")[,2]),
  #                                             str_match(ph_orp, "orpstandardoffset\\s*(.*?)\\s*mvtemperature")[,2],
  #                                             ifelse(is.na(str_match(ph_orp, "zobell'soffset\\s*(.*?)\\s*mvtemperature")[,2]) & is.na(str_match(ph_orp, "orpstandardoffset\\s*(.*?)\\s*mvtemperature")[,2]),
  #                                                    str_match(ph_orp, "quickcal\\s*(.*?)\\s*mvtemperature")[,2],
  #                                                    str_match(ph_orp, "zobell'soffset\\s*(.*?)\\s*mvtemperature")[,2])),
  #
  #                         # Conductivity
  #                         cond_cal_date = as.character(mdy(str_match(conductivity, "lastcalibrated\\s*(.*?)\\s*calibrationdetails")[,2])),
  #                         tds_conversion_ppm = str_sub(str_match(conductivity, "tdsconversionfactor\\s*(.*?)\\s*cellconstant")[,2], 6, nchar(str_match(conductivity, "tdsconversionfactor\\s*(.*?)\\s*cellconstant")[,2])),
  #
  #                         # calibration report formatting has changed in 2024 for this variable. Therefore a post-2024 correction must occur
  #                         cond_cell_constant = ifelse(year(DT) < 2024, str_match(conductivity, "cellconstant\\s*(.*?)\\s*referencetemperature")[,2],
  #                                                     str_match(conductivity, "cellconstant\\s*(.*?)\\s*offset")[,2]),
  #
  #                         cond_offset = ifelse(year(DT) < 2024, NA,
  #                                              str_match(conductivity, "offset\\s*(.*?)\\s*µs/cm")[,2]),
  #
  #                         cond_pre = str_match(conductivity,paste0(str_match(conductivity,
  #                                                                            "premeasurementactual\\s*(.*?)\\s*specificconductivity")[,2],"specificconductivity\\s*(.*?)\\s*µs/cmpost"))[,2],
  #                         cond_post = str_match(conductivity,paste0(str_match(conductivity,
  #                                                                             "postmeasurementactual\\s*(.*?)\\s*specificconductivity")[,2],"specificconductivity\\s*(.*?)\\s*µs/cm"))[,2]) %>%
  #                         # if(turbidity == "No Turbidity Sensor"){
  #                         # # Turbidity
  #                         # turb_cal_date = "None",
  #                         # ntu_slope = "None",
  #                         # ntu_offset = "None",
  #                         # ntu_10 = "None",
  #                         # ntu_100 = "None") %>%
  #
  #     select(-c(ph_7_nice, ph_7_high, ph_7_low))
  #
  #     # Not all sondes have depth.
  #     if(!is.na(str_match(depth, "lastcalibrated"))){#\\s*(.*?)\\s*calibrationdetails")[,2])){
  #       cal_table <- cal_table %>%
  #         mutate(
  #           # Depth
  #           depth_cal_date = as.character(mdy(str_match(depth, "lastcalibrated\\s*(.*?)\\s*calibrationdetails")[,2])),
  #           depth_offset = str_match(depth, "zerooffset\\s*(.*?)\\s*psireferencedepth")[,2],
  #           depth_ref_depth = str_match(depth, "psireferencedepth\\s*(.*?)\\s*ftreferenceoffset")[,2],
  #           depth_ref_offset = str_match(depth, "ftreferenceoffset\\s*(.*?)\\s*psipremeasurement")[,2],
  #           depth_pre_psi = str_match(depth, "psipremeasurement\\s*(.*?)\\s*psipostmeasurement")[,2],
  #           depth_post_psi = str_match(depth, "psipostmeasurement\\s*(.*?)\\s*psi")[,2])
  #     }
  #
  #     if(depth == "No Depth Sensor"){
  #
  #       cal_table <- cal_table %>%
  #         mutate(# Depth
  #       depth_cal_date = "No Depth Sensor",
  #       depth_offset = "No Depth Sensor",
  #       depth_ref_depth = "No Depth Sensor",
  #       depth_ref_offset = "No Depth Sensor",
  #       depth_pre_psi = "No Depth Sensor",
  #       depth_post_psi = "No Depth Sensor")
  #     }
  #
  #
  #     if(!is.na(str_match(turbidity, "lastcalibrated"))){#calibrationpoint1premeasurement\\s*(.*?)\\s*ntupost")[,2])){
  #       # Not all sondes have turbidity.
  #       cal_table <- cal_table %>%
  #         mutate(
  #           # Turbidity
  #           turb_cal_date = as.character(mdy(str_match(turbidity, "lastcalibrated\\s*(.*?)\\s*calibrationdetails")[,2])),
  #           ntu_slope = str_match(turbidity, "slope\\s*(.*?)\\s*offset")[,2],
  #           ntu_offset = str_match(turbidity, "offset\\s*(.*?)\\s*ntu")[,2],
  #           ntu_10 = str_match(turbidity, "calibrationpoint1premeasurement\\s*(.*?)\\s*ntupost")[,2],
  #           ntu_100 = str_match(turbidity, "calibrationpoint2premeasurement\\s*(.*?)\\s*ntupost")[,2])
  #     }
  #
  #     if(turbidity == "No Turbidity Sensor"){
  #       cal_table <- cal_table %>%
  #         mutate(
  #           # Turbidity
  #           turb_cal_date = "No Turbidity Sensor",
  #           ntu_slope = "No Turbidity Sensor",
  #           ntu_offset = "No Turbidity Sensor",
  #           ntu_10 = "No Turbidity Sensor",
  #           ntu_100 = "No Turbidity Sensor")
  #
  #
  #
  #
  #     }
  #
  #
  #     cal_table <- cal_table %>%
  #       mutate(
  #         #Factory Defaults
  #         factory_defaults = paste0(ifelse(grepl("factorydefault", turbidity), "Turbidity ", ""),
  #                                   ifelse(grepl("factorydefault", rdo), "RDO ", ""),
  #                                   ifelse(is.na(ph_slope_post), "pH ", ""),
  #                                   ifelse(is.na(orp_offset), "ORP ", ""),
  #                                   ifelse(grepl("factorydefault", conductivity), "Conductivity ", ""),
  #                                   ifelse(grepl("factorydefaults", depth), "Depth ", ""))) %>%
  #       # convert all columns to character values to preserve info
  #       mutate(across(.cols = everything(), .fns = as.character)) %>%
  #       # remove "," from big numbers
  #       mutate(across(everything(), ~str_replace_all(., ",", "")))
  #
  #   }
  #
  #   bad_cal_interval_list <- map2(
  #     .x = bad_cal_records_filtered$start_DT,
  #     .y = bad_cal_records_filtered$end_DT,
  #     .f = ~interval(.x, .y, tz = "MST"))
  #
  #   if(df_parameter == "DO"){
  #
  #     cal_table <- list.files("data/calibration_reports", pattern=".html") %>%
  #       .[grepl(df_site, ., ignore.case = TRUE)] %>%
  #       map_dfr(., cal_tabler) %>%
  #       distinct(.keep_all = TRUE) %>%
  #       mutate(DT = as.character(round_date(ymd_hms(DT, tz = "MST"), "15 minutes"))) %>%
  #       # mutate(across(-matches("date|site|DT|factory"), as.numeric)) %>%
  #       dplyr::select(DT, site, rdo_slope, rdo_offset)
  #
  #     df_mod <- df %>%
  #       left_join(., cal_table, by = c("DT_join" = "DT", "site")) %>%
  #       fill(names(cal_table)[!grepl("\\b(site|DT)\\b", names(cal_table))], .direction = "down") %>%
  #       mutate(rdo_slope_pre = as.numeric(rdo_slope), rdo_offset_pre = as.numeric(rdo_offset)) %>%
  #       select(names(df), contains(c("pre"))) %>%
  #       left_join(., cal_table, by = c("DT_join" = "DT", "site")) %>%
  #       fill(names(cal_table)[!grepl("\\b(site|DT)\\b", names(cal_table))], .direction = "up") %>%
  #       mutate(rdo_slope_post = as.numeric(rdo_slope), rdo_offset_post = as.numeric(rdo_offset)) %>%
  #       select(names(df), contains(c("pre", "post"))) %>%
  #       #mutate(raw = (mean -rdo_offset_pre)/rdo_offset_pre)
  #       mutate(raw = case_when(DT_round %within% bad_cal_interval_list & is.na(rdo_slope_pre) ~ mean,
  #                              DT_round %within% bad_cal_interval_list & !is.na(rdo_slope_pre) ~ ((mean - rdo_offset_pre)/rdo_slope_pre),
  #                              .default = mean),
  #              cal_fix = case_when(DT_round %within% bad_cal_interval_list ~ (raw*rdo_slope_post) + rdo_offset_post,
  #                                  .default = mean)) %>%
  #       add_flag(mean != cal_fix, "calibration fix") %>%
  #       mutate(raw = mean,
  #              mean = cal_fixed)
  #
  #
  #   }
  #
  #   if(df_parameter == "pH"){
  #
  #     cal_table <- list.files("data/calibration_reports", pattern=".html") %>%
  #       .[grepl(df_site, ., ignore.case = TRUE)] %>%
  #       map_dfr(., cal_tabler) %>%
  #       distinct(.keep_all = TRUE) %>%
  #       mutate(DT = as.character(round_date(ymd_hms(DT, tz = "MST"), "15 minutes"))) %>%
  #       # mutate(across(-matches("date|site|DT|factory"), as.numeric)) %>%
  #       dplyr::select(DT, site, ph_slope_pre, ph_offset_pre, ph_slope_post, ph_offset_post, factory_defaults)
  #
  #     df_mod <- df %>%
  #       left_join(., cal_table, by = c("DT_join" = "DT", "site")) %>%
  #       fill(names(cal_table)[!grepl("\\b(site|DT)\\b", names(cal_table))], .direction = "down") %>%
  #       mutate(ph_slope_pre = as.numeric(ph_slope_pre), ph_offset_pre = as.numeric(ph_offset_pre),
  #              ph_slope_post = as.numeric(ph_slope_post), ph_offset_post = as.numeric(ph_offset_post)) %>%
  #       mutate(raw = case_when(DT_round %within% bad_cal_interval_list & is.na(ph_slope_pre) & grepl("pH", factory_defaults, ignore.case = TRUE) ~ mean,
  #                              DT_round %within% bad_cal_interval_list & !is.na(ph_slope_pre) & !grepl("pH", factory_defaults, ignore.case = TRUE) ~ ((mean - ph_offset_pre)/ph_slope_pre),
  #                              .default = mean),
  #              cal_fix = case_when(DT_round %within% bad_cal_interval_list ~ (raw*ph_slope_post) + ph_offset_post,
  #                                  .default = mean)) %>%
  #       add_flag(mean != cal_fix, "calibration fix")
  #
  #   }
  #
  #   if(df_parameter == "")
  #
  # }
# }

fix_depth_cal <- function(df){


  if(!"Depth" %in% df$parameter){

    nope <- df %>% mutate(relative_depth = NA)

   return(nope)

  }

  if("Depth" %in% df$parameter & "archery" %in% df$site){

     #df <- all_data_flagged[["archery-Depth"]]

    site_depth <- df %>%
      mutate(relative_depth = ifelse(DT_round >= as_datetime('2022-05-21 15:45:00', "MST") & DT_round <= as_datetime('2022-05-24 15:45:00', "MST"),
                                     mean +
                                       abs(
                                         dplyr::filter(df, as.character(DT_round) == ('2022-05-24 15:15:00'))$mean -
                                           dplyr::filter(df, as.character(DT_round) == ('2022-05-24 15:45:00'))$mean
                                       ),
                                     mean))
    #return(archery_depth)

  } else if("Depth" %in% df$parameter & "timberline" %in% df$site){

    # df <- all_data_flagged[["timberline-Depth"]]

    site_depth <- df %>%
      dplyr::mutate(relative_depth = ifelse(year == "2022" & DT_round <= lubridate::as_datetime('2022-04-07 17:00:00', "MST"),
                                            as.numeric(mean) +
                                              abs(
                                                dplyr::filter(df, as.character(DT_round) == ('2022-04-07 16:15:00'))$mean -
                                                  dplyr::filter(df, as.character(DT_round) == ('2022-04-07 17:15:00'))$mean
                                              ),
                                            as.numeric(mean)))
    # return(timberline_depth)

   } else if ("Depth" %in% df$parameter & "legacy" %in% df$site) {

    # df <- all_data_flagged[["legacy-Depth"]]

    site_depth <- df %>%
      mutate(relative_depth = ifelse(DT_round >= as_datetime('2022-04-06 06:00:00', "MST") & DT_round <= as_datetime('2022-04-12 09:15:00', "MST"),
                                     as.numeric(mean) +
                                       abs(
                                         dplyr::filter(df, as.character(DT_round) == ('2022-04-12 09:15:00'))$mean -
                                           dplyr::filter(df, as.character(DT_round) == ('2022-04-12 09:30:00'))$mean
                                       ),
                                     mean))

    site_depth <-  site_depth %>%
      mutate(relative_depth = ifelse(DT_round >= as_datetime('2022-07-08 17:00:00', "MST") & DT_round <= as_datetime('2022-07-12 09:00:00', "MST"),
                                     relative_depth +
                                       abs(
                                         dplyr::filter(site_depth, as.character(DT_round) == ('2022-07-08 14:15:00'))$relative_depth -
                                           dplyr::filter(site_depth, as.character(DT_round) == ('2022-07-08 17:00:00'))$relative_depth
                                       ),
                                     relative_depth))

    site_depth <-  site_depth %>%
      mutate(relative_depth = ifelse(DT_round >= as_datetime('2022-07-22 11:30:00', "MST") & DT_round <= as_datetime('2022-07-25 14:15:00', "MST"),
                                     relative_depth +
                                       abs(
                                         dplyr::filter(site_depth, as.character(DT_round) == ('2022-07-22 09:45:00'))$relative_depth -
                                           dplyr::filter(site_depth, as.character(DT_round) == ('2022-07-22 11:30:00'))$relative_depth
                                       ),
                                     relative_depth))

       #return(legacy_depth)

  } else if ("Depth" %in% df$parameter & "tamasag" %in% df$site) {

    # df <- all_data_flagged[["tamasag-Depth"]]

    site_depth <- df %>%
      mutate(relative_depth = ifelse(DT_round <= "2022-04-24 07:15:00" & year == "2022",
                                     as.numeric(mean) +
                                       abs(
                                         dplyr::filter(df, as.character(DT_round) == ('2022-04-24 07:15:00'))$mean -
                                           dplyr::filter(df, as.character(DT_round) == ('2022-04-24 07:30:00'))$mean
                                       ),
                                     mean))

    # return(tamasag_depth)

  } else if("Depth" %in% df$parameter) {

    site_depth <- df %>%
             mutate(relative_depth = mean)
}

  depth_flagged <- site_depth %>%
    add_flag(., relative_depth != mean, "post-calibration") %>%
    dplyr::mutate(mean = relative_depth)

  return(depth_flagged)

}

#' @title Apply all flags to a data frame
#'
#' @description
#' A function that applies all flags to a data frame. This function is used to
#' apply all flags to a data frame in one step.
#'
#' @param data A data frame with a `flag` column.
#'
#' @param require An upstream dependency (i.e., alternative to `{necessary object}`).
#'
#' @return A data frame with a `flag` column that has been updated with all flags.
#'
#' @examples
#' flag_all_data(data = all_data_flagged$`archery-Actual Conductivity`)
#' flag_all_data(data = all_data_flagged$`boxelder-Temperature`)
#'
#' @seealso [add_field_flag()]
#' @seealso [add_spec_flag()]
#' @seealso [add_seasonal_flag()]
#' @seealso [add_na_flag()]
#' @seealso [add_repeat_flag()]
#' @seealso [add_suspect_flag()]
#' @seealso [add_malfunction_flag()]

flag_all_data <- function(data, require = NULL) {
  flagged_data <- data %>%
    add_field_flag() %>%
    add_spec_flag() %>%
    add_seasonal_flag() %>%
    add_na_flag() %>%
    add_repeat_flag() %>%
    add_suspect_flag() %>%
    add_malfunction_flag(malfunction_records = mWater_malfunction_records)
    # mutate(mean_public = ifelse(is.na(flag), mean, NA))
  return(flagged_data)
}

generate_daily_plot <- function(plot_data_arg, df_list_arg, site_arg, parameter_arg) {

  # site_param_df <- df_list_arg[[paste0(site_arg, "-", parameter_arg)]]

  start_date <- min(plot_data_arg$DT_round)
  end_date <- max(plot_data_arg$DT_round)

  #default is lower network
  site_vector <- c("tamasag", # rist
                   "legacy",
                   "lincoln",
                   "timberline",
                   "prospect",
                   "boxelder", # elc
                   "archery",
                   "river bluffs")

  if(network == "virridy"){ # this will be the new default
    # establish order for all the non-tributary sites
    sites_order <-  c("joei","cbri","chd","pfal","sfm","pbd","tamasag",
                      "legacy","lincoln","timberline","prospect","boxelder",
                      "archery","riverbluffs")
    # establish the order for the tributary sites
    trib_sites_order <- c("boxcreek", "archery", NA, "springcreek", "prospect",
                          NA, "penn", "sfm", "lbea")
  }

  # determining the index for the site of interest.
  if (site_arg %in% sites_order) {

    plot_filter <- tibble(site = c("joei","cbri","chd","pfal","sfm","pbd",
                                   "tamasag","legacy", "lincoln","timberline",
                                   "timberline virridy","prospect",
                                   "prospect virridy","boxelder","archery",
                                   "archery virridy","riverbluffs"))

    site_index <- which(sites_order == site_arg)
    site_list <- as.vector(na.omit(sites_order[max(1, site_index - 1):min(length(sites_order), site_index + 1)]))

    plot_filter <- plot_filter %>%
      filter(grepl(paste(site_list, collapse = "|"), site, ignore.case = TRUE),
             site != site_arg) %>%
      pull(site)

  } else {

    plot_filter <- tibble(site = c("boxcreek", "archery", "archery virridy",
                                   "springcreek", "prospect", "prospect virridy",
                                   "penn", "sfm", "lbea"))

    site_index <- which(trib_sites_order == site_arg)
    site_list <- as.vector(na.omit(trib_sites_order[max(1, site_index - 1):min(length(trib_sites_order), site_index + 1)]))

    plot_filter <- plot_filter %>%
      filter(grepl(paste(site_list, collapse = "|"), site, ignore.case = TRUE),
             site != site_arg) %>%
      pull(site)

  }

  # Get the relevant sonde data
  # relevant_sondes <- map(plot_filter,
  #                        ~ {
  #                          sonde_name <- paste0(.x,"-",parameter_arg)
  #                          tryCatch({
  #                            sonde_df <- df_list_arg[[sonde_name]]  %>%
  #                              filter(DT_round %within% interval(start_date, end_date))},
  #                            error = function(err) {
  #                              cat("Sonde ", sonde_name," not found.\n")})
  #                        })

  # get the relevant sonde data source
  relevant_sonde_source <- map(plot_filter, ~ {
    sonde_name <- paste0(.x, "-", parameter_arg)
    # Determine which directory to pull data from
    tryCatch({
      retrieve_relevant_data_name(sonde_name, interval_arg = interval(start_date, end_date))
    }, error = function(err) {
      return("all_data")
    })
  })

  # Get the relevant data
  relevant_sondes <- map2(plot_filter,
                          relevant_sonde_source,
                          function(name, source) {
                            sonde_name <- paste0(name, "-", parameter_arg)
                            # try to pull in the data
                            tryCatch({
                              get(source)[[sonde_name]] %>%
                                filter(DT_round %within% interval(start_date, end_date))
                            }, error = function(err) {
                              return(NULL)
                            })
                          })

  # combine the lists
  sonde_info <- map2(relevant_sondes, relevant_sonde_source, list)

  # Remove any NULL results from the list
  sonde_info <- keep(sonde_info, ~!is.null(.x[[1]]))

  # append site_df to relevant sonde list, clean list, and bind dfs
  # to find plot info
  relevant_dfs <- map(sonde_info, ~.x[[1]])

  # append plot_data_arg to relevant sonde list, clean list, and bind dfs
  daily_plot_data <- append(relevant_dfs, list(plot_data_arg)) %>%
    keep(~ !is.null(.)) %>%
    keep(~ nrow(.)>0) %>%
    bind_rows()

  # use the daily flag data day as flag_day
  flag_day <- min(plot_data_arg$DT_round)

  plot <- ggplot(data = daily_plot_data) +
    geom_point(data = filter(daily_plot_data, (site == site_arg)),
               aes(x=DT_round, y=mean, color=flag)) +
    # geom_line(data = filter(daily_plot_data, (site != site_arg)),
    #           aes(x=DT_round, y=mean, color=site)) +
    map(sonde_info, function(sonde_data) {
      data <- sonde_data[[1]]
      data_source <- sonde_data[[2]]

      y_column <- if (data_source == "all_data") "mean" else "mean_verified"

      geom_line(data = data, aes(x = DT_round, y = .data[[y_column]], color = site))
    }) +
    ggtitle(paste0(str_to_title(site_arg), " ", str_to_title(parameter_arg), " (", format(flag_day, "%B %d, %Y"), ")")) +
    labs(x = "Time",
         y = "Mean")

  plot <- add_threshold_lines(plot = plot,
                              plot_data = plot_data_arg,
                              site_arg = site_arg,
                              parameter_arg = parameter_arg)

  plot <- plot +
    geom_vline(xintercept = seq(min(plot_data_arg$DT_round), max(plot_data_arg$DT_round), by = "hour"),
               color = "lightgrey", size = 0.5) +
    scale_x_datetime(breaks = seq(min(plot_data_arg$DT_round), max(plot_data_arg$DT_round), by = "hour"),
                     labels = format(seq(min(plot_data_arg$DT_round), max(plot_data_arg$DT_round), by = "hour"), "%H")) +
    theme_bw() +
    theme(legend.position = 'right',
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1)) + # make the background white
    guides(color = guide_legend(nrow = 10, byrow = TRUE))

  return(plot)

}

# Generate histograms of a site parameter data frame based on the mean column.
# @param df A dataframe with a `mean` column.
# @param df_index The index of the dataframe.
# @return A histogram of the mean column.
# @examples
# generate_general_histogram(df = all_data_flagged$`archery-Actual Conductivity`, df_index = "archery-Actual Conductivity")
# generate_general_histogram(df = all_data_flagged$`boxelder-Temperature`, df_index = "boxelder-Temperature")
generate_general_histogram <- function(df, df_index) {

  site_param <- toupper(sub("-", " ", df_index, fixed = TRUE))
  n <- as.character(sum(!is.na(df$mean)))
  title <- paste0("Histogram of ", site_param, " (n = ", n, ")")

  # there should be checks for data here
  minimum <- floor(min(df$mean, na.rm = TRUE))
  maximum <- ceiling(max(df$mean, na.rm = TRUE))

  histogram <- ggplot(data = df, aes(x = mean)) +
    geom_histogram(
      breaks = seq(minimum, maximum, by = 1)
    ) +
    labs(title = title)

  return(histogram)
}

generate_initial_weekly_plots <- function(all_df_list, pending_df_list, site_arg, parameter_arg, flag_arg = NULL) {

  site_param <- paste0(site_arg, "-", parameter_arg)

  site_flag_dates <- pending_df_list[[site_param]] %>%
    group_by(y_w) %>%
    filter(any(verification_status == "SKIP") | any(!is_verified)) %>%
    ungroup()

  if (!is.null(site_flag_dates)){
    # vector of sites in the order that they are in spatially ----
    # some sites have some funkiness going on (not all of the sites are present in the final plot)
    #default is lower network
    sites_order <- c("tamasag", # rist
                     "legacy",
                     "lincoln",
                     "timberline",
                     "prospect",
                     "boxelder", # elc
                     "archery",
                     "river bluffs")

    if(network == "virridy"){ # this will be the new default
      # establish order for all the non-tributary sites
      sites_order <-  c("joei","cbri","chd","pfal","sfm","pbd","tamasag",
                        "legacy","lincoln","timberline","prospect","boxelder",
                        "archery","riverbluffs")
      # establish order for all the tributary sites
      trib_sites_order <- c("boxcreek", "archery", NA, "springcreek", "prospect",
                            NA, "penn", "sfm", "lbea")
    }

    # determining the sites relevant to the site of interest.
    if (site_arg %in% sites_order) {

      plot_filter <- tibble(site = c("joei","cbri","chd","pfal","sfm","pbd",
                                     "tamasag","legacy", "lincoln","timberline",
                                     "timberline virridy","prospect",
                                     "prospect virridy","boxelder","archery",
                                     "archery virridy","riverbluffs"))

      site_index <- which(sites_order == site_arg)
      site_list <- as.vector(na.omit(sites_order[max(1, site_index - 1):min(length(sites_order), site_index + 1)]))

      plot_filter <- plot_filter %>%
        filter(grepl(paste(site_list, collapse = "|"), site, ignore.case = TRUE),
               site != site_arg) %>%
        pull(site)

    } else {

      plot_filter <- tibble(site = c("boxcreek", "archery", "archery virridy",
                                     "springcreek", "prospect", "prospect virridy",
                                     "penn", "sfm", "lbea"))

      site_index <- which(trib_sites_order == site_arg)
      site_list <- as.vector(na.omit(trib_sites_order[max(1, site_index - 1):min(length(trib_sites_order), site_index + 1)]))

      plot_filter <- plot_filter %>%
        filter(grepl(paste(site_list, collapse = "|"), site, ignore.case = TRUE),
               site != site_arg) %>%
        pull(site)

    }

    if (nrow(site_flag_dates >0)) {

      if (is.null(flag_arg)) {
        # This for loop generates an overlayed plot of weekly data for the site of
        # interest sandwiched by the site above and below it for each day that was
        # tagged by a flag of interest
        plot_list <- list()

        grouped_data <- site_flag_dates %>%
          group_by(y_w) %>% #group_by(week, year) %>% # group_by(week, month, year) %>%
          group_split()

        for(i in 1:length(grouped_data)) {

          group_data <- grouped_data[[i]]

          year_week <- unique(group_data$y_w)

          # filtering dfs of interest for the week of interest
          site_df <- site_flag_dates %>%
            filter(y_w == year_week)

          # Get the relevant sonde data
          relevant_sondes <- map(plot_filter, ~ {
            sonde_name <- paste0(.x, "-", parameter_arg)
            data_source <- NULL
            sonde_df <- NULL

            # Determine which directory to pull data from
            tryCatch({
              data_source <- retrieve_relevant_data_name(sonde_name, year_week)
              # cat("Data for",sonde_name,"will be pulled from",data_source,"\n")
            }, error = function(err) {
              # cat("Data for",sonde_name,"not found.\n")
              return(NULL)  # Return NULL if data source can't be determined
            })

            # Only try to pull in the data if data_source was successfully determined
            if (!is.null(data_source)) {
              tryCatch({
                sonde_df <- get(data_source)[[sonde_name]] %>%
                  filter(y_w == group_data$y_w)
              }, error = function(err) {
                cat("Sonde", sonde_name, "not found.\n")
                return(NULL)  # Return NULL if sonde data can't be retrieved
              })
            }

            # Only return a list if both data_source and sonde_df are available
            if (!is.null(data_source) & !is.null(sonde_df)) {
              return(list(sonde_df = sonde_df, data_source = data_source))
            } else {
              return(NULL)  # Return NULL if either data_source or sonde_df is NULL
            }
          })

          # Remove any NULL results from the list
          relevant_sondes <- compact(relevant_sondes)

          # append site_df to relevant sonde list, clean list, and bind dfs
          # to find plot info
          relevant_dfs <- map(relevant_sondes, ~.x[[1]])
          week_plot_data <- append(relevant_dfs, list(site_df)) %>% # how to relevant sondes here
            keep(~ !is.null(.)) %>%
            keep(~ nrow(.)>0) %>%
            bind_rows() %>%
            arrange(day)

          # Create a sequence of dates for the vertical lines
          start_date <- floor_date(min(week_plot_data$DT_round), "day")
          end_date <- ceiling_date(max(week_plot_data$DT_round), "day")
          vline_dates <- seq(start_date, end_date, by = "day")

          date_seq <- seq(min(vline_dates), max(vline_dates) - days(1), by = "1 day") + hours(12)

          # Use the first day of the group as flag_day
          flag_day <- min(group_data$DT_round)

          week_plot <- ggplot(data = week_plot_data) +
            geom_point(data = filter(week_plot_data, (site == site_arg)),
                       aes(x=DT_round, y=mean, color=flag)) +
            map(relevant_sondes, function(sonde_data) {
              data <- sonde_data[[1]]
              data_source <- sonde_data[[2]]

              y_column <- if (data_source == "all_data") "mean" else "mean_verified"

              geom_line(data = data, aes(x = DT_round, y = .data[[y_column]], color = site))
            })+
            geom_vline(xintercept = vline_dates, color = "black") +
            ggtitle(paste0(str_to_title(site_arg), " ", parameter_arg, " (", format(flag_day, "%B %d, %Y"), ")")) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            labs(x = "Day",
                 y = "Mean")

          week_plot <- add_threshold_lines(plot = week_plot,
                                           plot_data = week_plot_data,
                                           site_arg = site_arg,
                                           parameter_arg = parameter_arg)

          week_plot <- week_plot +
            theme_bw() +
            scale_x_datetime(date_breaks = "1 day",
                             date_labels = "%b %d",
                             minor_breaks = date_seq,
                             sec.axis = sec_axis(~., breaks = date_seq, labels = unique(week_plot_data$weekday))) +
            theme(legend.position = 'bottom',
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()) +
            guides(color = guide_legend(nrow = 4, byrow = TRUE))

          plot_list[[paste(site_param, as.character(flag_day))]] <- week_plot
          sorted_plot_names <- names(plot_list)[order(names(plot_list))]

          plot_list <- plot_list[sorted_plot_names]

        }
      }
      # ---- if flag arg is not null ----
      if (!is.null(flag_arg)) {
        print("under construction...")
        # # This for loop generates an overlayed plot of weekly data for the site of
        # # interest sandwiched by the site above and below it for each day that was
        # # tagged by a flag of interest
        # plot_list <- list()
        #
        # site_flag_weeks <- site_flag_dates %>%
        #   filter(str_detect(flag, flag_arg)) %>%
        #   group_by(week, year) %>%
        #   slice(1)
        #
        # grouped_data <- site_flag_dates %>%
        #   filter(y_w %in% site_flag_weeks$y_w) %>%
        #   group_by(week, year) %>% # group_by(week, month, year) %>%
        #   group_split()
        #
        #
        # for(i in 1:length(grouped_data)) {
        #
        #   group_data <- grouped_data[[i]]
        #
        #   # flag_title <- site_flag_dates$flag[i] # no flag title ***
        #
        #   # filtering dfs of interest for the weeks where a flag was detected
        #   site_df <- site_flag_dates %>%
        #     filter(y_w == group_data$y_w)
        #
        #   # TryCatch used here to avoid erroring out on the first and last values of
        #   # sites_order object (there is no prior/next record after the first/last record).
        #   # Return df as NULL in case of an error
        #   prev_site_df <- NULL
        #   next_site_df <- NULL
        #
        #   tryCatch({
        #     previous_site <- paste0(sites_order[site_index-1],"-",parameter_arg)
        #     prev_site_df <- all_df_list[[previous_site]] %>%
        #       filter(y_w == group_data$y_w)},
        #     error = function(err) {
        #       cat("No previous site.\n")})
        #
        #   tryCatch({
        #     next_site <- paste0(sites_order[site_index+1],"-",parameter_arg)
        #     next_site_df <- all_df_list[[next_site]] %>%
        #       filter(y_w == group_data$y_w)},
        #     error = function(err) {
        #       cat("No next site.\n")})
        #
        #   # Bind all three dfs
        #   week_plot_data <- list(site_df, prev_site_df, next_site_df) %>%
        #     # remove NULL values from the list
        #     keep(~ !is.null(.)) %>%
        #     bind_rows()
        #
        #   # Create a sequence of dates for the vertical lines
        #   start_date <- floor_date(min(week_plot_data$DT_round), "day")
        #   end_date <- ceiling_date(max(week_plot_data$DT_round), "day")
        #   vline_dates <- seq(start_date, end_date, by = "day")
        #
        #   date_seq <- seq(min(vline_dates), max(vline_dates) - days(1), by = "1 day") + hours(12) ## here ----
        #
        #   # Use the first day of the group as flag_day
        #   flag_day <- min(group_data$DT_round)
        #
        #   week_plot <- ggplot(data = week_plot_data) +
        #     geom_point(data = filter(week_plot_data, (site == site_arg)),
        #                aes(x=DT_round, y=mean, color=flag)) +
        #     geom_line(data = filter(week_plot_data, (site != site_arg)),
        #               aes(x=DT_round, y=mean, color=site)) +
        #     geom_vline(xintercept = vline_dates, color = "black") +
        #     ggtitle(paste0(str_to_title(site_arg), " ", parameter_arg, " (", format(flag_day, "%B %d, %Y"), ")")) +
        #     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        #     labs(x = "Day",
        #          y = "Mean")
        #
        #   week_plot <- add_threshold_lines(plot = week_plot,
        #                                    plot_data = week_plot_data,
        #                                    site_arg = site_arg,
        #                                    parameter_arg = parameter_arg)
        #
        #   week_plot <- week_plot +
        #     theme_bw() +
        #     scale_x_datetime(date_breaks = "1 day",
        #                      date_labels = "%b %d",
        #                      minor_breaks = date_seq) +
        #     theme(legend.position = 'bottom',
        #           panel.grid.major = element_blank(),
        #           panel.grid.minor = element_blank()) +
        #     annotate("text", x = date_seq, y = min(week_plot_data$mean, na.rm = TRUE) - 1, label = 1:length(date_seq), hjust = 0) +
        #     guides(color = guide_legend(nrow = 4, byrow = TRUE))
        #
        #   plot_list[[paste(site_param, as.character(flag_day))]] <- week_plot
        #   sorted_plot_names <- names(plot_list)[order(names(plot_list))]
        #
        #   plot_list <- plot_list[sorted_plot_names]
        # }
      }
      return(plot_list)
    } else {
      return(paste(flag_arg, "not detected.\n"))
    }
  } else {
    return(paste(site_arg, parameter_arg, "combination not available.\n"))
  }

}

# Generate seasonal histograms for a given site parameter data frame based on the mean column.
# @param df A dataframe with a `mean` column.
# @param df_index The index of the dataframe.
# @return A plot with a histogram of the mean column for each season.
# @examples
# generate_seasonal_histogram(df = all_data_flagged$`archery-Actual Conductivity`, df_index = "archery-Actual Conductivity")
# generate_seasonal_histogram(df = all_data_flagged$`boxelder-Temperature`, df_index = "boxelder-Temperature")

generate_seasonal_histogram <- function(df, df_index) {

  winter_baseflow <- c(12,1,2,3,4)
  snowmelt <- c(5,6,NA,NA,NA)
  monsoon <- c(7,8,9,NA,NA)
  fall_baseflow <- c(10,11,NA,NA,NA)

  # do water years

  seasons <- data.frame(winter_baseflow, snowmelt, monsoon, fall_baseflow)

  site_param <- toupper(sub("-", " ", df_index, fixed = TRUE))

  param <- unique(na.omit(df$parameter))

  hist_list <- list()
  for (i in colnames(seasons)){

    filtered_df <- df %>%
      filter(month %in% seasons[[i]],
             !str_detect(flag, "sensor specification range"))

    n <- as.character(sum(!is.na(filtered_df$mean)))

    title <- paste0(i," (n = ",n ,")")

    histogram <- ggplot() +
      geom_histogram() +
      labs(title = title)

    tryCatch({
      minimum <- floor(min(filtered_df$mean, na.rm = TRUE))
      maximum <- ceiling(max(filtered_df$mean, na.rm = TRUE))

      if (param %in% c("Specific Conductivity", "Actual Conductivity", "Turbidity")) {
        bins <- seq(minimum, maximum, by = 10)
      } else if (param %in% c("ORP", "pH", "Depth")) {
        bins <- seq(minimum, maximum, by = 0.05)
      } else {
        bins <- seq(minimum, maximum)
      }

      x_min <- filtered_df$m_mean05[1]
      x_max <- filtered_df$m_mean99[1]

      histogram <- ggplot(data = filtered_df, aes(x = mean)) +
        geom_histogram(breaks = bins) +
        geom_vline(xintercept = x_min, color = "red", linetype = "dashed") +
        geom_vline(xintercept = x_max, color = "red", linetype = "dashed") +
        facet_wrap(~ year, nrow = 1) +
        labs(title = title)
    },
    error = function(err) {
      cat("No finite values for", site_param, i, "\n")})

    hist_list[[i]] <- histogram

  }

  collated_hist <- ggarrange(plotlist = hist_list, nrow=2, ncol=2) %>%
    annotate_figure(top = site_param)

  return(collated_hist)

}

#' @title Generate site metaparameter data from the HydroVu API.
#'
#' @description
#' A function that generates a site metaparameter dataframe from the HydroVu API.
#' A metaparameter is a parameter that is used to generate flags for parameters
#' other than itself.
#'
#' @param api_data A dataframe with the munged API data.
#'
#' @param require An upstream dependency (i.e., alternative to `{necessary object}`).
#'
#' @return A dataframe with site metaparameter data that will be joined to the other
#' site-parameter data.
#'
#' @examples
#' generate_site_metaparam(site_arg = "archery", metaparameter_arg = "Battery Level", api_data = incoming_data_collated_csvs)
#'
#' @seealso [summarize_site_param()]

generate_site_metaparam <- function(api_data, require = NULL) {

    sites <- unique(api_data$site)
    metaparameters <- c("Temperature", "Battery Level", "Baro", "External Voltage")

    df_list <- list()
    for (i in sites) {
        metaparameter_data <- api_data %>%
            data.table::data.table() %>%
            dplyr::select(DT_join, site, parameter, value) %>%
            dplyr::filter(site == i & (parameter %in% metaparameters)) %>%
            dplyr::select(-site) %>%
            tidyr::pivot_wider(names_from = parameter, values_from = value)
        df_list[[i]] <- metaparameter_data
    }
    return(df_list)
}

#' @title Generate Summary Statistics
#'
#' @description
#' A function that generates summary statistics for a given site parameter data frame.
#' The generated statistics include:
#'   - The next value and previous value for the mean.
#'   - The rolling 7-point median of the mean.
#'   - The rolling 7-point mean of the mean.
#'   - The rolling 7-point standard deviation of the mean.
#'   - The slope of a point in relation to the point ahead and behind.
#'   - The rolling 7-point slope of the mean.
#'   - The month and year of each data point.
#'   - The year-month combination.
#'
#' @param site_param_df A data frame with a `mean` column retrieved from HydroVu API.
#'
#' @return A data frame with summary statistics for a given site parameter data frame.
#'
#' @examples
#' generate_summary_statistics(site_param_df = all_data_flagged$`archery-Actual Conductivity`)
#' generate_summary_statistics(site_param_df = all_data_flagged$`boxelder-Temperature`)

generate_summary_statistics_full <- function(site_param_df) {

  summary_stats_df <- site_param_df %>%
    # ... so that we can get the proper leading/lagging values across our entire timeseries:
    dplyr::mutate(
      # Add the next value and previous value for mean.
      front1 = dplyr::lead(mean, n = 1),
      back1 = dplyr::lag(mean, n = 1),
      # Add the rolling 7-point median (using itself + data of the past).
      rollmed = RcppRoll::roll_median(mean, n = 7, align = 'right', na.rm = F, fill = NA_real_), #ifelse(is.na(rollmed), roll_median(mean, n = 7, align = 'right', na.rm = F, fill = NA_real_), rollmed), # to go (j): check_na() function for when we append data
      # Add the rolling 7-point mean (using itself + data of the past).
      rollavg = RcppRoll::roll_mean(mean, n = 7, align = 'right', na.rm = F, fill = NA_real_), #ifelse(is.na(rollavg), roll_mean(mean, n = 7, align = 'right', na.rm = F, fill = NA_real_), rollavg),
      # Add the rolling 7-point standard deviation (using itself + data of the past).
      rollsd = RcppRoll::roll_sd(mean, n = 7, align = 'right', na.rm = F, fill = NA_real_), #ifelse(is.na(rollsd), roll_sd(mean, n = 7, align = 'right', na.rm = F, fill = NA_real_), rollsd),
      # Determine the slope of a point in relation to the point ahead and behind
      slope_ahead = (front1 - mean)/15,
      slope_behind = (mean - back1)/15,
      # Get the rolling 7-point slope (using itself + data of the past).
      rollslope = RcppRoll::roll_mean(slope_behind, n = 7, align = 'right', na.rm = F, fill = NA_real_), #ifelse(is.na(rollslope), roll_mean(slope_behind, n = 7, align = 'right', na.rm = F, fill = NA_real_), rollslope),
      # add some summary info for future steps
      month = lubridate::month(DT_round),
      year = lubridate::year(DT_round),
      y_m = paste0(year, '-', month),
      # Define our seasons:
      season = dplyr::case_when(month %in% c(12, 1, 2, 3, 4) ~ "winter_baseflow",
                                month %in% c(5, 6) ~ "snowmelt",
                                month %in% c(7, 8, 9) ~ "monsoon",
                                month %in% c(10, 11) ~ "fall_baseflow",
                                TRUE ~ NA)
    )

  return(summary_stats_df)

}

# Generate summary statistics for a given site parameter data frame.
#' @param site_param_df A data frame with a `mean` column retrieved from HydroVu API.
#' @return A data frame with summary statistics for a given site parameter data frame.
#' @examples
# generate_summary_statistics(site_param_df = all_data_flagged$`archery-Actual Conductivity`)
# generate_summary_statistics(site_param_df = all_data_flagged$`boxelder-Temperature`)

generate_summary_statistics <- function(site_param_df) {

  summary_stats_df <- site_param_df %>%
    # ... so that we can get the proper leading/lagging values across our entire timeseries:
    dplyr::mutate(
      # Add the next value and previous value for mean.
      # Only do this for newest data (i.e., our appended historical
      # data already has these filled out and we don't want to over-
      # write them)
      front1 = ifelse(is.na(front1), dplyr::lead(mean, n = 1), front1),
      back1 = ifelse(is.na(back1), dplyr::lag(mean, n = 1), back1),
      # Add the median for a point and 6 points behind it:
      rollmed = ifelse(is.na(rollmed), RcppRoll::roll_median(mean, n = 7, align = 'right', na.rm = F, fill = NA_real_), rollmed),
      # Add the mean for a point and 6 points behind it:
      rollavg = ifelse(is.na(rollavg), RcppRoll::roll_mean(mean, n = 7, align = 'right', na.rm = F, fill = NA_real_), rollavg),
      # Add the standard deviation for a point and 6 points behind it:
      rollsd = ifelse(is.na(rollsd), RcppRoll::roll_sd(mean, n = 7, align = 'right', na.rm = F, fill = NA_real_), rollsd),
      # Determine the slope of a point in relation to the point ahead and behind.
      slope_ahead = ifelse(is.na(slope_ahead), (front1 - mean)/15, slope_ahead),
      slope_behind = ifelse(is.na(slope_behind), (mean - back1)/15, slope_behind),
      # Get the rolling 7-point slope (using itself + data of the past).
      rollslope = ifelse(is.na(rollslope), RcppRoll::roll_mean(slope_behind, n = 7, align = 'right', na.rm = F, fill = NA_real_), rollslope),
      # add some summary info for future us
      month = lubridate::month(DT_round),
      year = lubridate::year(DT_round),
      y_m = paste(year, '-', month),
      season = dplyr::case_when(month %in% c(12, 1, 2, 3, 4) ~ "winter_baseflow",
                                month %in% c(5, 6) ~ "snowmelt",
                                month %in% c(7, 8, 9) ~ "monsoon",
                                month %in% c(10, 11) ~ "fall_baseflow",
                                TRUE ~ NA)
    )

  return(summary_stats_df)

}

generate_supplemental_weekly_plot <- function(daily_plot_data_arg, df_list_arg, site_arg, parameter_arg) {

  site_param_df <- df_list_arg[[paste0(site_arg, "-", parameter_arg)]]

  start_date <- min(daily_plot_data_arg$DT_round) - days(3) # TODO: replace all daily_plot_data_arg instances with something more efficient
  end_date <- max(daily_plot_data_arg$DT_round) + days(3)

  week_plot_data <- site_param_df %>%
    filter(DT_round %within% interval(start_date, end_date))

  #default is lower network
  site_vector <- c("tamasag", # rist
                   "legacy",
                   "lincoln",
                   "timberline",
                   "prospect",
                   "boxelder", # elc
                   "archery",
                   "river bluffs")

  if(network == "virridy"){ # this will be the new default
    # establish order for all the non-tributary sites
    sites_order <-  c("joei","cbri","chd","pfal","sfm","pbd","tamasag",
                      "legacy","lincoln","timberline","prospect","boxelder",
                      "archery","riverbluffs")
    # establish the order for the tributary sites
    trib_sites_order <- c("boxcreek", "archery", NA, "springcreek", "prospect",
                          NA, "penn", "sfm", "lbea")
  }

  # determining the index for the site of interest.
  if (site_arg %in% sites_order) {

    plot_filter <- tibble(site = c("joei","cbri","chd","pfal","sfm","pbd",
                                   "tamasag","legacy", "lincoln","timberline",
                                   "timberline virridy","prospect",
                                   "prospect virridy","boxelder","archery",
                                   "archery virridy","riverbluffs"))

    site_index <- which(sites_order == site_arg)
    site_list <- as.vector(na.omit(sites_order[max(1, site_index - 1):min(length(sites_order), site_index + 1)]))

    plot_filter <- plot_filter %>%
      filter(grepl(paste(site_list, collapse = "|"), site, ignore.case = TRUE),
             site != site_arg) %>%
      pull(site)

  } else {

    plot_filter <- tibble(site = c("boxcreek", "archery", "archery virridy",
                                   "springcreek", "prospect", "prospect virridy",
                                   "penn", "sfm", "lbea"))

    site_index <- which(trib_sites_order == site_arg)
    site_list <- as.vector(na.omit(trib_sites_order[max(1, site_index - 1):min(length(trib_sites_order), site_index + 1)]))

    plot_filter <- plot_filter %>%
      filter(grepl(paste(site_list, collapse = "|"), site, ignore.case = TRUE),
             site != site_arg) %>%
      pull(site)

  }

  # get the relevant sonde data source
  relevant_sonde_source <- map(plot_filter, ~ {
    sonde_name <- paste0(.x, "-", parameter_arg)
    # Determine which directory to pull data from
    tryCatch({
      retrieve_relevant_data_name(sonde_name, interval_arg = interval(start_date, end_date))
    }, error = function(err) {
      return("all_data")
    })
  })

  # Get the relevant data
  relevant_sondes <- map2(plot_filter,
                          relevant_sonde_source,
    function(name, source) {
    sonde_name <- paste0(name, "-", parameter_arg)
    # try to pull in the data
    tryCatch({
      get(source)[[sonde_name]] %>%
        filter(DT_round %within% interval(start_date, end_date))
    }, error = function(err) {
      return(NULL)
    })
  })

  # combine the lists
  sonde_info <- map2(relevant_sondes, relevant_sonde_source, list)

  # Remove any NULL results from the list
  sonde_info <- keep(sonde_info, ~!is.null(.x[[1]]))

  # append site_df to relevant sonde list, clean list, and bind dfs
  # to find plot info
  relevant_dfs <- map(sonde_info, ~.x[[1]])

  week_plot_data <- append(relevant_dfs, list(week_plot_data)) %>%
    keep(~ !is.null(.)) %>%
    keep(~ nrow(.)>0) %>%
    bind_rows() %>%
    arrange(day)

  # Create a sequence of dates for the vertical lines
  start_date <- floor_date(min(week_plot_data$DT_round), "day")
  end_date <- ceiling_date(max(week_plot_data$DT_round), "day")
  vline_dates <- seq(start_date, end_date, by = "day")

  date_seq <- seq(min(vline_dates), max(vline_dates) - days(1), by = "1 day") + hours(12)

  # use the daily flag data day as flag_day
  flag_day <- unique(daily_plot_data_arg$DT_round)

  week_plot <- ggplot(data = week_plot_data) +
    geom_point(data = filter(week_plot_data, (site == unique(daily_plot_data_arg$site))),
               aes(x=DT_round, y=mean, color=flag)) +
    map(sonde_info, function(sonde_data) {

      data <- sonde_data[[1]]
      data_source <- sonde_data[[2]]

      y_column <- if (data_source == "all_data") "mean" else "mean_verified"

      geom_line(data = data, aes(x = DT_round, y = .data[[y_column]], color = site))
    }) +
    geom_rect(data = daily_plot_data_arg, aes(xmin = min(DT_round), xmax = max(DT_round),
                                              ymin = -Inf, ymax = Inf),
              fill = "grey",
              alpha = 0.01,
              color = NA) +
    geom_vline(xintercept = vline_dates, color = "black") +
    ggtitle(paste0(str_to_title(unique(daily_plot_data_arg$site)), " ", unique(daily_plot_data_arg$parameter), " (", format(flag_day, "%B %d, %Y"), ")")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Day",
         y = "Mean")

  week_plot <- add_threshold_lines(plot = week_plot,
                                   plot_data = week_plot_data,
                                   site_arg = unique(daily_plot_data_arg$site),
                                   parameter_arg = unique(daily_plot_data_arg$parameter))

  week_plot <- week_plot +
    theme_bw() +
    scale_x_datetime(date_breaks = "1 day",
                     date_labels = "%b %d",
                     minor_breaks = date_seq) +
    theme(legend.position = 'right',
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) + # make the background white
    guides(color = guide_legend(nrow = 10, byrow = TRUE))

  return(week_plot)
}

# Get day decision.
# @param prompt_text A string for the prompt text.
# @return A string of "pass", "fail", or "inspect" depending on input from the user.
# @examples
# get_day_decision(prompt_text = "Would you like to (pass/fail/inspect) all data points for: ")

get_day_decision <- function(prompt_text) {
  while (TRUE) {
    user_input <- readline(prompt = paste(prompt_text))
    user_input <- tolower(user_input)

    if (user_input %in% c("pass", "p")) {
      return("pass")
    } else if (user_input %in% c("fail", "f")) {
      return("fail")
    } else if (user_input %in% c("inspect", "i")){
      return("inspect")
    } else {
      cat("Invalid input. Please enter 'yes', 'no', or 'inspect'.\n")
    }
  }
}

get_dt_inspection_decisions <- function(daily_plot_data) {
  # Create a sequence of date times from the daily plot data
  dt_list <- sort(unique(pull(daily_plot_data, DT_round)))

  # Get the date
  plot_date <- as_date(daily_plot_data$DT_round[1])

  # Prompt for the number of intervals
  prompt_text <- paste("How many intervals would you like to inspect?\n(Must be less than", length(dt_list), ")\n")

  while (TRUE) {
    user_input <- readline(prompt = prompt_text)

    # The user input needs to check for:
    # Only numbers between 1 and length(dt_list)
    num_intervals <- suppressWarnings(as.numeric(user_input))

    if (!is.na(num_intervals) && num_intervals > 0 && num_intervals < length(dt_list)) {
      break
    }

    cat("Invalid input. Please enter a valid number of intervals.\n")
  }

  selected_intervals <- list()

  for (i in 1:num_intervals) {
    prompt_text <- paste("Enter the time range for interval", i, ".\n(format 'HH:MM:SS-HH:MM:SS'):\n")

    while (TRUE) {
      user_input <- readline(prompt = prompt_text) # these need to be 15 minute intervals ***
      interval <- unlist(strsplit(user_input, "-"))

      interval_set <- paste(plot_date, interval)
      interval_set <- as_datetime(interval_set, tz = "MST")

      if (length(interval_set) == 2 && all(interval_set %in% dt_list)) {
        start_time <- interval_set[1]
        end_time <- interval_set[2]

        # Check if selected_intervals is empty
        if (start_time <= end_time & length(selected_intervals) == 0) {
          selected_intervals[[i]] <- interval(start_time, end_time)
          break
        }

        if (start_time <= end_time & !suppressWarnings(any(map(selected_intervals, ~int_overlaps(interval(start_time, end_time), .x))))) {
          selected_intervals[[i]] <- interval(start_time, end_time)
          break
        }
      }

      cat("Invalid input.\nPlease enter a valid time range that doesn't overlap with previous intervals.\n")
    }
  }

  return(selected_intervals)
}

#' #' Import Fort Collins Floodwarning Rain Gage Data
#' #'
#' #' This function pulls rain data for sites within the FC Floodwarning System
#' #'
#' #'
#' #' @param start Start date of when you want data.
#' #' @param end End date of when you want data; default is the current date.
#' #' @param save Whether to save (TRUE) the resulting table or not (FALSE)
#' #' @param path If `save = TRUE`, the file path to save the shapefile
#' #'
#' #' @return A table of time series flow data across the FC Floodwarning System
#' #'
#' get_fc_rain <- function(start = '2018-10-01', end = Sys.Date(), save = TRUE, path = 'data/context_data/'){
#'
#'   call <- "https://opendata.fcgov.com/api/views/g87z-rviz/rows.csv?accessType=DOWNLOAD"
#'
#'   #download dataset from FC
#'   temp1 <- tempfile()
#'   download.file(paste0(call), destfile = temp1, method = "curl")
#'   data <- read_csv(temp1) %>%
#'     dplyr::mutate(date = as_date(mdy_hms(Timestamp))) %>%
#'     dplyr::filter(date >= start & date <= end)
#'
#'   s_platte_ws s_platte_flowlines <- nhdplusTools::get_nhdplus(AOI = s_platte_ws,
#'                                                               realization = "flowline") <- nhdplusTools::get_nldi_basin(list(featureSource = "nwissite", featureID = "USGS-06754000"))
#'   s_platte_flowlines <- nhdplusTools::get_nhdplus(AOI = s_platte_ws,
#'                                                   realization = "flowline")
#'   s_platte_catchments <- nhdplusTools::get_nhdplus(AOI = s_platte_ws, realization = "catchment")
#'
#'   poudre <- s_platte_flowlines %>%
#'     filter(grepl("Cache la Poudre", gnis_name, ignore.case = TRUE))
#'
#'   poudre <- s_platte_catchments %>%
#'     filter(featureid %in% poudre$comid) %>%
#'     .[sites,] %>%
#'     rowid_to_column()
#'
#'   sites <- read_csv("data/metadata/sonde_location_metadata.csv") %>%
#'     separate(lat_long, into = c("lat", "long"), sep = ",") %>%
#'     sf::st_as_sf(coords = c("long", "lat"), crs = 4326) %>%
#'     filter(Site %in% c("Tamasag","Legacy","Lincoln","Timberline","Prospect","Boxelder","Archery","River Bluffs")) %>%
#'     mutate(rowid = sf::st_nearest_feature(., poudre)) %>%
#'     left_join(sf::st_drop_geometry(poudre), by = "rowid")
#'
#'   ws_maker <- function(site_names){
#'
#'     df <- sites %>%
#'       dplyr::filter(Site == site_names)
#'
#'     # get_UT() creates a vector of all flowlines upstream of the comid of interest ...
#'     UT_comids <- nhdplusTools::get_UT(network = s_platte_flowlines,
#'                                       comid = df$featureid)
#'
#'     catchments <- filter(s_platte_catchments, featureid %in% c(UT_comids)) %>%
#'       summarize() %>%
#'       mutate(site = df$Site)
#'
#'     return(catchments)
#'
#'   }
#'
#'   watersheds <- sites$Site %>%
#'     map_dfr(~ws_maker(.))
#'
#'   rain_sites <- fread("data/context_data/fc_rain.csv") %>%
#'     distinct(`Sensor Name`, Latitude, Longitude) %>%
#'     sf::st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>%
#'     .[watersheds,] %>%
#'     sf::st_join(watersheds) %>%
#'     mutate(site = tolower(site)) %>%
#'     write_csv("data/context_data/site_link.csv")
#'
#'   if(save == TRUE){
#'     write_csv(data, paste0(path, "/fc_rain.csv"))
#'   }
#'
#'   return(data)
#' }
#'
#'
#'
#'

# Get flag decision.
# @param prompt_text A string for the prompt text.
# @return A boolean of TRUE or FALSE depending on input from the user.
# @examples
# get_flag_decision(prompt_text = "Would you like to (pass/fail) this data point: ")

get_flag_decision <- function(prompt_text) {
  while (TRUE) {
    user_input <- readline(prompt = paste(prompt_text))
    user_input <- tolower(user_input)

    if (user_input %in% c("pass", "p")) {
      return(TRUE)
    } else if (user_input %in% c("fail", "f")) {
      return(FALSE)
    } else {
      cat("Invalid input. Please enter 'pass' or 'fail'.\n")
    }
  }
}

#' @title Get the start dates for the HydroVu API pull from the historically
#' flagged data
#'
#' @description
#' This function finds the max datetime in the temperature data
#' frames for each site's historically flagged data. We use temperature because
#' it is always being tracked in the sensors and thus is the most accurate
#' representation of the last API pull for a given sensor.
#'
#' @param incoming_historically_flagged_data_list This is the output from
#' `flagged_data_dfs` and is the historically flagged data.
#'
#' @returns
#' A data frame that has the sites and their corresponding starting
#' date-times' for the HydroVu API pull. This data frame will then get passed into
#' `incoming_data_csvs_upload`.

get_start_dates_df <- function(incoming_historically_flagged_data_list) {

  temperature_subset <- grep("Temperature", names(incoming_historically_flagged_data_list))

  start_dates_df <- map_dfr(incoming_historically_flagged_data_list[temperature_subset],
                            function(temperature_df){
                              temperature_df %>%
                                select(DT_round, site) %>%
                                filter(DT_round == max(DT_round))
                              })

  return(start_dates_df)

}

# @ param site: string, name of site to get filenames for, select from "tamasag" "legacy", "timberline" "prospect" "boxelder" "archery" "riverbluffs"
get_tl_photo_filenames <- function(site = "legacy",start_dt = "2023-05-20 12:00" ,end_dt = "2023-07-20 12:00"){
  #convert dates to DT objects
  start_dt <- ymd_hm(start_dt, tz = "MST")
  end_dt <- ymd_hm(end_dt, tz = "MST")
  if(site == "river bluffs"){
    site <- "riverbluffs"
  }

  if(site %nin% c("tamasag", "legacy", "timberline", "prospect", "boxelder", "archery", "riverbluffs")){
    print("Invalid site name, please select from: 'tamasag' 'legacy' 'timberline' 'prospect' 'boxelder' 'archery' 'riverbluffs'")

  nah <- tibble()
    return(nah)
  }

  #look through all the files in the site's folder
  all_files <- tibble(filename = list.files(paste0("data/timelapse_photos/2023_compiled/",site), full.names = TRUE, recursive = TRUE))%>%
    mutate(DT = str_extract(filename, "[0-9]{8}_[0-9]{4}"),
           DT = parse_date_time(DT, "Ymd_HM", tz = "MST"),
           #round DT to 30 min to match with sensor data
           DT_round = round_date(DT, "30 minutes"))%>%
    # if there are multiple photos with the same DT_round value, only keep one row with that DT_round
    distinct(DT_round, .keep_all = TRUE)%>%
    select(-DT)%>%
    filter(between(DT_round, start_dt, end_dt))

  if(nrow(all_files) == 0){
    print("No files found for this site and date range")
    nah <- tibble()
    return(nah)
  }else{
   return(all_files)
  }


}

#test <- get_tl_photo_filenames(site = "legacy", start_dt = "2023-7-20 12:00", end_dt = "2023-12-06 12:00")





#' @title Get verification decision
#'
#' @description
#' Get decisions from user input
#'
#' @param prompt_text The prompt that the user is replying to
#'
#' @examples
#' # get_day_decision(prompt_text = "Would you like to (pass/fail/inspect) all data points for: ")

get_verification_decision <- function(prompt_text) { # should add a layer level to this to prevent user from inspect sub-daily data

  while (TRUE) {
    user_input <- readline(prompt = paste(prompt_text))
    user_input <- tolower(user_input)

    # pass statements
    if (user_input %in% c("pass all", "pass", "pa", "p")) {
      return("PASS ALL")
    }
    if (user_input %in% c("pass valid", "pv")) {
      return("PASS VALID")
    }
    if (user_input %in% c("pass flagged", "pf")) {
      return("PASS FLAGGED")
    }
    if (user_input %in% c("pass none", "pn")) {
      return("PASS NONE")
    }

    # fail statements
    if (user_input %in% c("fail all", "fail", "fa", "f")) {
      return("FAIL ALL")
    }
    if (user_input %in% c("fail valid", "fv")) {
      return("FAIL VALID")
    }
    if (user_input %in% c("fail flagged", "ff")) {
      return("FAIL FLAGGED")
    }
    if (user_input %in% c("fail none", "fn")) {
      return("FAIL NONE")
    }

    # skip statements
    if (user_input %in% c("skip", "skip all", "s", "sa")) {
      return("SKIP")
    }

    # inspect statements ***
    if (user_input %in% c("inspect all", "inspect", "ia", "i")){
      return("INSPECT ALL")
    }

    if (user_input %in% c("inspect valid", "iv")){
      return("INSPECT VALID")
    }

    if (user_input %in% c("inspect flagged", "if")){
      return("INSPECT FLAGGED")
    }

    if(user_input %in% c("inspect some", "is")){ # if we are at daily level we can't go further
      return("INSPECT SOME")
    }

    # quit statement
    if (user_input %in% c("quit", "q")) {
      return("QUIT")
    }

    cat("Invalid input. Please enter one of the options from the following decision matrix:\n\n")
    cat("          | PASS | FAIL | SKIP | INSPECT  |  | QUIT |\n")
    cat("----------+------+------+------+----------|  |------|\n")
    cat(" ALL      | pa   | fa   | sa   | ia       |  | q    |\n")
    cat(" VALID    | pv   | fv   |      | iv       |  +------+\n")
    cat(" FLAGGED  | pf   | ff   |      | if       |\n")
    cat(" NONE     | pn   | fn   |      |          |\n")
    cat(" SOME     |      |      |      | is       |\n")
    cat("----------+------+------+------+----------|\n")
  }

}


get_watersheds <- function(sf = NULL, coordinates = c(-112.2124, 37.63281)){

  sf::sf_use_s2(FALSE)

  if(is.null(sf)){

    # Create a data frame with a column named 'geometry'
    df <- tibble::tibble(long = coordinates[1],
                         lat = coordinates[2])

    aoi <- sf::st_as_sf(df, coords = c("long", "lat"), crs = 4326)
  }

  if(is.null(coordinates)){

    aoi <- sf %>% sf::st_transform(4326)

  }

  aoi <- aoi %>%
    nhdplusTools::get_nhdplus(AOI = .)

  # Use the NHDPlus digital elevation model to find the nearest downslope
  # NHD flowline for any point in space (here, our point of interest)
  # Not working anymore -
  # trace <- get_raindrop_trace(aoi, direction = "up")
  #
  # "Snap" our site to the nearest NHD flowline feature
  # snap_point <- sf::st_sfc(sf::st_point(trace$intersection_point[[1]][1:2]),
  #                          crs=4326)
  #
  # Clip/split our catchment to only include the portion of the
  # catchment upstream of our site:
  # better_termination <- get_split_catchment(snap_point, upstream = F)[2,]

  # read in the complete NHD (in tabular form) to make for much more efficient nhd crawling.
  # This data in tabular form doesn't exist anywhere online that I know of... -_-
  nhd <- readr::read_csv('data/context_data/nhd_flow_network.csv')

  upstream <- nhdplusTools::get_UT(nhd, aoi$comid) %>% #upstream trace function in nhdplusTools
    tibble::as_tibble() %>%
    dplyr::rename(comid_list = value)  %>%
    dplyr::distinct(comid_list, .keep_all = TRUE)

  nhd_catch <-  upstream$comid_list %>%
    map(~nhdplusTools::get_nhdplus(comid = .,
                                   realization='catchment',
                                   t_srs = 4269)) %>%
    dplyr::bind_rows() %>%
    dplyr::distinct(featureid,.keep_all=TRUE) %>%
    dplyr::group_by(featureid) %>%
    dplyr::summarize()

  return(nhd_catch)

}



get_weekly_inspection_decision <- function(weekly_plot_data) {

  prompt_text <- "Which days do you want to verify (1/2/3/4/5/6/7)? \nCan be any combination of weekdays, with no other characters.\nex: 456\n"
  day_min <- min(weekly_plot_data$weekday)
  day_max <- max(weekly_plot_data$weekday)

  while (TRUE) {

    user_input_og <- readline(prompt = paste(prompt_text))

    user_input <- suppressWarnings(
      sort(na.omit(as.numeric(unique(unlist(strsplit(user_input_og, ""))))) [as.numeric(unique(unlist(strsplit(user_input_og, "")))) >= day_min & as.numeric(unique(unlist(strsplit(user_input_og, "")))) <= day_max])
      ) # 1 and 7 will be set to data min and max

    if (length(user_input) != 0) {
      return(user_input)
    }

    cat(user_input_og, " is not a valid input")
  }
}

#' @title Get working data decision
#'
#' @description
#' Get working directory and data decisions from user input
#'
#' @param prompt_text The prompt that the user is replying to
#'
#' @examples
#' # get_day_decision(prompt_text = "Would you like to (pass/fail/inspect) all data points for: ")

get_working_data_decision <- function() { # should add a layer level to this to prevent user from inspect sub-daily data

  prompt_text <- "Which directory are you working from? (pre/int): "

  while (TRUE) {
    user_input <- readline(prompt = paste(prompt_text))
    user_input <- tolower(user_input)

    if (user_input == "pre") {
      working_data <<- set_names(map(list.files(pre_verification_path, full.names = TRUE), readRDS), list.files(pre_verification_path))
      break()
    }

    if (user_input == "int") {
      working_data <<- set_names(map(list.files(intermediary_path, full.names = TRUE), readRDS), list.files(intermediary_path))
      break()
    }

    cat("Invalid input. Please enter one of the options from the following:\n\n")
    cat("pre = pre_verification_dir\n")
    cat("int = intermediary_dir\n")
  }

}

grab_confirmed_mWater_malfunctions <- function(field_notes){

    # Grab notes about sensor malfunction
    malfunction_records <- all_notes_cleaned %>%
      dplyr::filter(grepl("Sensor malfunction", visit_type, ignore.case = TRUE)) %>%
      dplyr::select(start_DT, site, crew, which_sensor_malfunction, malfunction_end_dt, notes = visit_comments)

    # write to csv
    write_csv(malfunction_records, "data/mWater_malfunciton_records.csv")

  }

grab_mWater_malfunction_notes <- function(mWater_api_data){

  # Grab notes about sensor malfunction
  malfunction_records <- mWater_api_data %>%
    filter(grepl("Sensor malfunction", visit_type, ignore.case = TRUE)) %>%
    select(start_DT, site, crew, which_sensor_malfunction, malfunction_end_dt, notes = visit_comments)

  #write to csv
  write_csv(malfunction_records, "data/mWater_malfunction_records.csv")

  parameters <- c("Battery Level",
                  "Baro",
                  "Chl-a Fluorescence",
                  "Depth",
                  "DO",
                  "External Voltage",
                  "ORP",
                  "pH",
                  "Specific Conductivity",
                  "Temperature",
                  "Turbidity")

  malfunction_records <- malfunction_records %>%
    # keep records relevant to {target} analysis
    select(start_DT, end_DT = malfunction_end_dt, site, parameter = which_sensor_malfunction, notes) %>%
    # match the text in the sensor column to the text in the target analysis
    separate_rows(parameter, sep = ", ") %>%
    mutate(
      parameter = case_when(
      parameter == "Chlorophyll a" ~ "Chl-a Fluorescence",
      parameter == "RDO" ~ "DO",
      parameter == "Conductivity" ~ "Specific Conductivity",
      .default = parameter
    ),
    site = case_when(
      site == "riverbluffs" ~ "river bluffs",
      .default = site
    )) %>%
    filter((is.na(parameter)) | (parameter %in% parameters))

  return(malfunction_records)

}

grab_mWater_sensor_malfunctions <- function(){

  # API Pull of mWater submitted notes

  # Grab API url from yml
  # Contact Sam Struthers if you need access
  creds <- yaml::read_yaml("src/mWater_collate/mWater_API.yml")
  api_url <- as.character(creds["url"])

  # Read in from API and tidy for downstream use

  # This is basic tidying of data set to:
  # correct datetime from UTC to Denver time (always MST)
  # correct columns where Other input is allowed (Site, visit type, photos downloaded, sensor malfunction)
  # Add rounded date time

  mal_notes <- readr::read_csv(url(api_url), show_col_types = FALSE) %>%
    dplyr::mutate(
      # start and end dt comes in as UTC -> to MST
      start_DT = lubridate::with_tz(lubridate::parse_date_time(start_dt, orders = c("%Y%m%d %H:%M:%S", "%m%d%y %H:%M", "%m%d%Y %H:%M", "%b%d%y %H:%M")), tz = "MST"),
      end_dt = lubridate::with_tz(lubridate::parse_date_time(end_dt, orders = c("%Y%m%d %H:%M:%S", "%m%d%y %H:%M", "%m%d%Y %H:%M", "%b%d%y %H:%M" )), tz = "MST"),
      malfunction_end_dt = with_tz(lubridate::parse_date_time(malfunction_end_dt, orders = c("%Y%m%d %H:%M:%S", "%m%d%y %H:%M", "%m%d%Y %H:%M", "%b%d%y %H:%M" )), tz = "MST"),
      date = as.Date(start_DT, tz = "MST"),
      start_time_mst = format(start_DT, "%H:%M"),
      sensor_pulled = as.character(sn_removed),
      sensor_deployed = as.character(sn_deployed),
      # If other is chosen, make site == other response
      site = ifelse(site == "Other (please specify)", tolower(stringr::str_replace_all(site_other, " ", "")), site),
      # When I changed the mWater survey, I accidentally introduced ??? in the place of Sensor Calibration option, fixing that here
      visit_type = dplyr::case_when(stringr::str_detect(visit_type, "\\?\\?\\?") ~ stringr::str_replace(string = visit_type,
                                                                                                        pattern =  "\\?\\?\\?",
                                                                                                        replacement = "Sensor Calibration or Check"),
                                    TRUE ~ visit_type),
      # Merging visit_type and visit type other
      visit_type = dplyr::case_when(stringr::str_detect(visit_type, "Other") ~ stringr::str_replace(string = visit_type,
                                                                                                    pattern =  "Other \\(please specify\\)",
                                                                                                    replacement = visit_type_other),
                                    TRUE ~ visit_type),
      # Merge sensor malfunction and sensor malfunction other
      which_sensor_malfunction = dplyr::case_when(stringr::str_detect(which_sensor_malfunction, "Other") ~ stringr::str_replace(string = which_sensor_malfunction,
                                                                                                                                pattern =  "Other \\(please specify\\)",
                                                                                                                                replacement = as.character(other_which_sensor_malfunction)),
                                                  TRUE ~ which_sensor_malfunction),
      # If other is chosen, make photos downloaded equal to response
      photos_downloaded = ifelse(photos_downloaded == "Other (please specify)", photos_downloaded_other, photos_downloaded),
      # Rounded start date time
      DT_round = lubridate::floor_date(start_DT, "15 minutes")) %>%
    # arrange by most recent visit
    dplyr::arrange(DT_round) %>%
    # Remove other columns
    dplyr::select(-c(photos_downloaded_other, visit_type_other, site_other, other_which_sensor_malfunction)) %>%
    dplyr::filter(grepl("sensor malfunction", visit_type, ignore.case = TRUE)) %>%
    dplyr::select(malfunction_start_dt = DT_round,  malfunction_end_dt, which_sensor_malfunction)

  return(mal_notes)

}

grab_mWater_sensor_notes <- function(mWater_api_data){

  # Sensor Notes

  # These are the notes that will be added to the QAQC workflow notes Most of the code in this chunk is to get the df to
  # match the one in the QAQC workflow It can be saved as a CSV or pulled directly into QAQC workflow
  # grab only notes where technician is interacting with sensor on site (excludes sensor malfunction notes)

  mWater_field_notes <- mWater_api_data %>%
    filter(grepl("Sensor",visit_type, ignore.case = TRUE) & !grepl("Sensor malfunction",visit_type, ignore.case = TRUE)) %>%
    # determining sonde employed status based on sensor_change
    mutate(sonde_employed = case_when(is.na(sensor_change)  ~ NA,
                                      sensor_change == "Swapped" ~ NA,
                                      sensor_change == "Pulled" ~ 1,
                                      sensor_change == "Deployed" ~ 0),
                                      #sensor_change %in% c("Swapped", "Deployed") ~ 1),

           #Sensor swapped notes
           sensor_swapped_notes = case_when(is.na(sensor_change)  ~ NA,
                                            sensor_change == "Pulled" &!is.na(sensor_pulled) ~ paste0("SN Removed: ", sensor_pulled),
                                            sensor_change == "Swapped" ~ paste0("SN Removed: ", sensor_pulled, " SN Deployed: ", sensor_deployed),
                                            sensor_change == "Deployed" ~ sensor_deployed),
           #Date/field season columns to match QAQC workflow
           DT_join = as.character(DT_round),
           field_season = year(DT_round),
           last_site_visit = DT_round,
           date = as.character(date)
    )%>%
    arrange(desc(DT_round))%>%
    #order columns in easily readable ordering
    select(site, crew, DT_round,sonde_employed,  sensors_cleaned, wiper_working, rdo_cap_condition, rdo_cap_replaced , ph_junction_replaced ,
           cals_performed, cal_report_collected , sensor_malfunction,sensor_pulled,sensor_deployed, sensor_swapped_notes,
           visit_type,start_time_mst,DT_join,  start_DT, end_dt,date,  visit_comments,photos_downloaded, field_season, last_site_visit)

  # back up to CSV
  # write_csv(sensor_notes, "data/mWater_sensor_field_notes.csv")

  # rm(all_notes_cleaned)
  return(mWater_field_notes)

}


hv_data_id <- function(loc_id, start_time = startdate, end_time = enddate, tz = timezone, token) {

  # convert the time to timestamp, convert to UTC for lookup in HydroVu
  start <- as.numeric(lubridate::with_tz(lubridate::ymd_hms(start_time, tz = tz), tzone = "UTC"))
  end <- as.numeric(lubridate::with_tz(lubridate::ymd_hms(end_time, tz = tz), tzone = "UTC"))

  # build the url
  url = "https://www.hydrovu.com/public-api/v1/locations/"
  url <- paste0(url, loc_id, "/data?endTime=", end, '&startTime=', start)

  req <- httr2::request(url)
  print(paste0('Trying site ', loc_id))
  try({
    resp <-  req %>% httr2::req_oauth_client_credentials(token) %>% httr2::req_perform()
    data <- list(resp %>% httr2::resp_body_json())
    h <- resp %>% httr2::resp_headers()

    while (!is.null(h[["X-ISI-Next-Page"]]))
    {
      resp <- req %>% httr2::req_headers("X-ISI-Start-Page" = h[["X-ISI-Next-Page"]]) %>%
        httr2::req_oauth_client_credentials(token) %>% httr2::req_perform()
      data <- c(data, list(resp %>% httr2::resp_body_json()))
      h <- resp %>% httr2::resp_headers()
    }

    # get the params and units
    params <- hv_names(token, return = "params")
    units <- hv_names(token, return = "units")

    # collapse the paginated date and clean up
    df <- purrr::map_dfr(data, flatten_page_params) %>%
      dplyr::mutate(timestamp = lubridate::with_tz(lubridate::as_datetime(timestamp, tz = "UTC"), tzone = tz),
                    Location = loc_id) %>%
      dplyr::inner_join(params, by = "parameterId") %>%
      dplyr::inner_join(units, by = "unitId") %>%
      dplyr::select(-parameterId, -unitId) %>%
      dplyr::arrange(Parameter, timestamp)

      return(df)
  })

}

#' Return the list of locations for the given client from HydroVu
#'
#' @param client a valid OAuth2 token such as returned from \code{hv_auth()}
#' @param url HydroVu url that lists the locations
#'
#' @return a dataframe listing all the locations visible to the client
#' @export
#'
#' @examples
#' \dontrun{
#' locs <- hv_locations(client)
#' }

hv_locations_all <- function(client,
                         url = "https://www.hydrovu.com/public-api/v1/locations/list") {

  req <- httr2::request(url)

  try({
  resp <-  req %>% httr2::req_oauth_client_credentials(client) %>% httr2::req_perform()
  locs <- list(resp %>% httr2::resp_body_json())
  h <- resp %>% httr2::resp_headers()

  while (!is.null(h[["X-ISI-Next-Page"]]))
  {
    resp2 <- req %>%
      httr2::req_headers("X-ISI-Start-Page" = h[["X-ISI-Next-Page"]]) %>%
      httr2::req_oauth_client_credentials(client) %>%
      httr2::req_perform()
    locs <- c(locs, list(resp2 %>% httr2::resp_body_json()))
    h <- resp2 %>% httr2::resp_headers()
  }
  # collapse the paginated date and clean up
  df <- flatten_df(locs) %>%
    select(-gps) %>%
    filter(!duplicated(.))
  return(df)
  })
}


#' @title intersensor_check
#'
#' @description
#' A function designed to reduce overflagging: if a slope violation occurs at the
#' same time as a slope violation in either depth or temperature, it is likely
#' not a sensor malfunction and instead a real product of the river.
#'
#' @param df An updated data frame with a `flag` column.
#'
#' @return A data frame with a `flag` column that has been updated with the
#' inter-sensor flag reduction step.
#'
#' @seealso [add_flag()]

intersensor_check <- function(df){

  # create a df of temperature for each site
  temperature <- df %>%
    data.table::data.table() %>%
    dplyr::select(DT_round, DT_join, parameter, mean, flag) %>%
    dplyr::filter(parameter == "Temperature") %>%
    dplyr::select(DT_join, Temperature = parameter, Temperature_flag = flag) %>%
    dplyr:: mutate(Temperature_front1 = dplyr::lead(Temperature_flag, n = 1),
                   Temperature_back1 = dplyr::lag(Temperature_flag, n = 1))

  # create a df of depth for each site
  depth <- df %>%
    data.table::data.table() %>%
    dplyr::select(DT_round, DT_join, parameter, mean, flag) %>%
    dplyr::filter(parameter == "Depth") %>%
    dplyr::select(DT_join, Depth = parameter, Depth_flag = flag) %>%
    dplyr:: mutate(Depth_front1 = dplyr::lead(Depth_flag, n = 1),
                   Depth_back1 = dplyr::lag(Depth_flag, n = 1))

  # add "temperature" and "depth" data columns to df:
  intersensors_checked <- df %>%
    dplyr::filter(!parameter %in% c("Depth", "Temperature")) %>%
    dplyr::left_join(., temperature, by = "DT_join") %>%
    dplyr::left_join(., depth, by = "DT_join") %>%
    # If either the depth or temperature have the same flag as a given parameter
    # identified at the same time (or one obs before/after), tag it
    dplyr::mutate(intersensored = dplyr::case_when(grepl("slope violation", flag) &
                                                     (grepl("slope violation", Depth_flag)   | grepl("slope violation", Temperature_flag)   |
                                                        grepl("slope violation", Depth_front1) | grepl("slope violation", Temperature_front1) |
                                                        grepl("slope violation", Depth_back1)  | grepl("slope violation", Temperature_back1)
                                                     ) ~ TRUE)) %>%
    dplyr::mutate(flag = ifelse(is.na(intersensored), flag, stringr::str_replace(flag, "slope violation", "")))

  final_checked_data <- df %>%
    dplyr::filter(parameter %in% c("Depth", "Temperature")) %>%
    # After using the temp and depth slope flags, remove that flagging entirely
    # from those parameters. We have yet to find an instance of the slope flag
    # capturing "fake" spikes in either of those data sets:
    dplyr::mutate(flag = stringr::str_replace(flag, "slope violation", "")) %>%
    dplyr::bind_rows(., intersensors_checked) %>%
    dplyr::select(-c(Depth, Depth_flag, Temperature, Temperature_flag))

  return(final_checked_data)

}

#' @title Load and tidy mWater field notes
#'
#' @description A function that uploads and cleans the field notes submitted to mWater.
#'
#' @param creds A .yml file with necessary credentials for accessing the field notes. Contact Sam Struthers if you need access.
#'
#' @return A dataframe with the field notes.

load_mWater_notes <- function(creds = yaml::read_yaml("creds/mWaterCreds.yml")){

  # API Pull of mWater submitted notes

  # Grab API url from yml
  # Contact Sam Struthers if you need access
  api_url <- as.character(creds["url"])

  # Read in from API and tidy for downstream use

  # This is basic tidying of data set to:
  # correct datetime from UTC to Denver time (MST)
  # correct columns where Other input is allowed (Site, visit type, photos downloaded, sensor malfunction)
  # Add rounded date time

  all_notes_cleaned <- readr::read_csv(url(api_url), show_col_types = FALSE) %>%
    dplyr::mutate(
      # start and end dt comes in as UTC -> to MST
      start_DT = lubridate::with_tz(lubridate::parse_date_time(start_dt, orders = c("%Y%m%d %H:%M:%S", "%m%d%y %H:%M", "%m%d%Y %H:%M", "%b%d%y %H:%M" )), tz = "MST"),
      end_dt = lubridate::with_tz(lubridate::parse_date_time(end_dt, orders = c("%Y%m%d %H:%M:%S", "%m%d%y %H:%M", "%m%d%Y %H:%M", "%b%d%y %H:%M" )), tz = "MST"),
      malfunction_end_dt = lubridate::with_tz(lubridate::parse_date_time(malfunction_end_dt, orders = c("%Y%m%d %H:%M:%S", "%m%d%y %H:%M", "%m%d%Y %H:%M", "%b%d%y %H:%M" )), tz = "MST"),
      date = as.Date(start_DT, tz = "MST"),
      start_time_mst = format(start_DT, "%H:%M"),
      sensor_pulled = as.character(sn_removed),
      sensor_deployed = as.character(sn_deployed),
      # If other is chosen, make site == other response
      site = ifelse(site == "Other (please specify)", tolower(stringr::str_replace_all(site_other, " ", "")), site),
      # When I changed the mWater survey, I accidentally introduced ??? in the place of Sensor Calibration option, fixing that here
      visit_type = dplyr::case_when(stringr::str_detect(visit_type, "\\?\\?\\?") ~ stringr::str_replace(string = visit_type,
                                                                               pattern =  "\\?\\?\\?",
                                                                               replacement = "Sensor Calibration or Check"),
                             TRUE ~ visit_type),
      # Merging visit_type and visit type other
      visit_type = dplyr::case_when(stringr::str_detect(visit_type, "Other") ~ stringr::str_replace(string = visit_type,
                                                                           pattern =  "Other \\(please specify\\)",
                                                                           replacement = visit_type_other),
                             TRUE ~ visit_type),
      # Merge sensor malfunction and sensor malfunction other
      which_sensor_malfunction = dplyr::case_when(stringr::str_detect(which_sensor_malfunction, "Other") ~ stringr::str_replace(string = which_sensor_malfunction,
                                                                                                       pattern =  "Other \\(please specify\\)",
                                                                                                       replacement = as.character(other_which_sensor_malfunction)),
                                           TRUE ~ which_sensor_malfunction),
      # If other is chosen, make photos downloaded equal to response
      photos_downloaded = ifelse(photos_downloaded == "Other (please specify)", photos_downloaded_other, photos_downloaded),
      # Rounded start date time
      DT_round = lubridate::floor_date(start_DT, "15 minutes")) %>%
    # arrange by most recent visit
    dplyr::arrange(DT_round)%>%
    # Remove other columns
    dplyr::select(-c(photos_downloaded_other,visit_type_other, site_other, other_which_sensor_malfunction ))

  return(all_notes_cleaned)

}

#' @title Load and tidy old field notes
#'
#' @description A function that uploads and cleans the field notes excel file. This function adds datetime
#' columns to the field notes dataframe and filters out field notes where the sensor
#' was not handled.
#'
#' @param filepath A file path to the raw field notes.
#'
#' @return A dataframe with the field notes.
#'
#' @examples
#' clean_old_field_notes(filepath = "data/sensor_field_notes.xlsx")

load_old_field_notes <- function(filepath){

  raw_field_notes <- readxl::read_excel(filepath)

  field_notes <- raw_field_notes %>%
    dplyr::mutate(start_DT = lubridate::ymd_hm(paste(date, start_time_mst), tz = "MST")) %>%
    dplyr::mutate(
      DT_round = lubridate::floor_date(start_DT, "15 minutes"),
      DT_join = as.character(DT_round),
      site = tolower(site),
      field_season = lubridate::year(DT_round),
      last_site_visit = DT_round) %>%
    dplyr::arrange(site, DT_round) %>%
    # rename instances of old names:
    dplyr::mutate(site = ifelse(site == "rist", "tamasag",
                         ifelse(site == "elc", "boxelder", site))) %>%
    # `sonde_employed` determines if the sonde is deployed or not. 0 = sonde deployed, 1 = sonde is not deployed
    mutate(sonde_employed = dplyr::case_when(!is.na(sensor_pulled) & !is.na(sensor_deployed) ~ 0,
                                             !is.na(sensor_pulled) & is.na(sensor_deployed) ~ 1,
                                             is.na(sensor_pulled) & !is.na(sensor_deployed) ~ 0,
                                             is.na(sensor_pulled) & is.na(sensor_deployed) ~ NA),
           end_dt  = as.POSIXct(NA, tz = "MST")) %>%
    # remove field dates where sensor was not handled:
    dplyr::filter(grepl("Sensor Cleaning or Check|Sensor Calibration", visit_type, ignore.case = TRUE))

  return(field_notes)

}

#' @title Generate Threshold Table for QAQC
#'
#' @description
#' A function designed to generate a threshold table for QAQC. This table
#' contains the thresholds for the mean, slope_behind, and standard deviation
#' of the mean for each site and season.
#'
#' @param df A data frame with a `flag` column.
#'
#' @return A data frame with the thresholds for the mean, slope_behind, and
#' standard deviation of the mean for each site and season.
#'
#' @examples
#' make_threshold_table(df = all_data_flagged$`archery-Actual Conductivity`)

make_threshold_table <- function(df){

  # sensor_malfunction_notes <- grab_mWater_malfunction_notes(mWater_api_data = load_mWater_notes())

  slope_down <- df %>%
    # REMOVE DATA WE KNOW TO BE ERRONEOUS:
    add_field_flag() %>%
    # incoroporate back-calibrated values for development of thresholds:
    fix_calibration() %>%
    add_malfunction_flag(df = ., malfunction_records = sensor_malfunction_notes) %>%
    add_spec_flag() %>%
    # Filter to remove times when the sonde was under freexzing conditions or not submerged (depth <= 0)
    add_frozen_flag() %>%
    # fix_depth_cal() %>%
    add_unsubmerged_flag() %>%
    # remove times when sonde was moved (to avoid crazy slopes)
    add_depth_shift_flag() %>%
    # Filter to only un-flagged data, and data that was collected when sonde
    # was actually deployed
    filter(is.na(flag) | sonde_employed != 1) %>%
    #... and when the sonde wasn't moved
    filter(is.na(depth_change) | depth_change != "sonde moved") %>%
    # Get threshold for negative slope data
    filter(slope_behind < 0) %>%
    group_by(season) %>%
    summarize(f_slope_behind_01 = quantile(slope_behind, 0.01, na.rm = TRUE))

  slope_up <- df %>%
    # REMOVE DATA WE KNOW TO BE ERRONEOUS:
    add_field_flag() %>%
    fix_calibration() %>%
    add_malfunction_flag(df = ., malfunction_records = sensor_malfunction_notes) %>%
    add_spec_flag() %>%
    # Filter to remove times when the sonde was under freexzing conditions or not submerged (depth <= 0)
    add_frozen_flag() %>%
    # fix_depth_cal() %>%
    add_unsubmerged_flag() %>%
    # remove times when sonde was moved (to avoid crazy slopes)
    add_depth_shift_flag() %>%
    # Filter to only un-flagged data, and data that was collected when sonde
    # was actually deployed
    filter(is.na(flag) | sonde_employed != 1) %>%
    #... and when the sonde wasn't moved
    filter(is.na(depth_change) | depth_change != "sonde moved") %>%
    # Get threshold for positive slope data
    filter(slope_behind > 0) %>%
    group_by(season) %>%
    summarize(f_slope_behind_99 = quantile(slope_behind, 0.99, na.rm = TRUE))

  good_data_stats <- df %>%
    # REMOVE DATA WE KNOW TO BE ERRONEOUS:
    add_field_flag() %>%
    fix_calibration() %>%
    add_malfunction_flag(df = ., malfunction_records = sensor_malfunction_notes) %>%
    add_spec_flag() %>%
    # Filter to remove times when the sonde was under freexzing conditions or not submerged (depth <= 0)
    add_frozen_flag() %>%
    # fix_depth_cal() %>%
    add_unsubmerged_flag() %>%
    # remove times when sonde was moved (to avoid crazy slopes)
    add_depth_shift_flag() %>%
    # Filter to only un-flagged data, and data that was collected when sonde
    # was actually deployed
    filter(is.na(flag) | sonde_employed != 1) %>%
    #... and when the sonde wasn't moved
    filter(is.na(depth_change) | depth_change != "sonde moved") %>%
    group_by(season) %>%
    # join our slope data thresholds:
    left_join(slope_up, by = "season") %>%
    left_join(slope_down, by = "season") %>%
    # develop other thresholds across all data
    mutate(f01 = quantile(mean, 0.01, na.rm = TRUE),
           f99 = quantile(mean, 0.99, na.rm = TRUE)) %>%
           # f_slope_behind_01 = slope_down, #quantile(slope_behind, 0.01, na.rm = TRUE),
           # f_slope_behind_99 = slope_up) %>% #quantile(slope_behind, 0.99, na.rm = TRUE)) %>%
    # THEN, GET STANDARD DEVIATION OF ONLYYYY VALUES WITHIN THE 1-99th PERCENTILE OF THAT GOOD DATA:
    filter(mean > f01 & mean < f99) %>%
    # SD is the ONLY statistic that uses this winnowed-down data set in its development.
    # All else use the full, "good" data set.
    summarize(site = paste0(unique(site)),
              parameter = paste0(unique(parameter)),
              t_mean01 = as.numeric(paste0(unique(f01))),
              t_mean99 = as.numeric(paste0(unique(f99))),
              t_slope_behind_01 = as.numeric(paste0(unique(f_slope_behind_01))),
              t_slope_behind_99 = as.numeric(paste0(unique(f_slope_behind_99))),
              # This stat is useless. Should remove eventually.
              t_sd_0199 = sd(mean, na.rm = T))

  return(good_data_stats)

}

#' @title Munge API data for QAQC workflow
#'
#' @description
#' A function designed to munge the raw API data for the QAQC workflow.
#'
#' @param api_path Path where the raw API data lives.
#' @param network Options include "csu", "virridy"
#' @param require An upstream dependency (i.e., alternative to `{necessary object}`).
#'
#' @return A dataframe with the munged API data.
#'
#' @examples
# munge_api_data(api_path = "data/api/incoming/")

munge_api_data <- function(api_path, network = "csu", require = NULL) {

  api_data <- list.files(path = api_path, full.names = TRUE, pattern = "*.csv") %>%
    purrr::map_dfr(~data.table::fread(.) %>%
                     dplyr::select(-id)) %>%
    # remove overlapping API-pull data
    dplyr::distinct()

  if(network %in% c("csu", "CSU")){

  api_data <- api_data %>%
    # remove VuLink data
    dplyr::filter(!grepl("vulink", name, ignore.case = TRUE)) %>% # these will throw us errors if there is none?
    # remove Virridy data
    dplyr::filter(!grepl("virridy", name, ignore.case = TRUE)) %>%
    dplyr::select(-name) %>%
    # Convert UTC (as it is sent from HydroVU API) to MST:
    dplyr::mutate(DT = lubridate::as_datetime(timestamp, tz = "UTC")) %>%
    dplyr::mutate(DT = lubridate::with_tz(DT, tzone = "MST"),
           DT_round = lubridate::round_date(DT, "15 minutes"),
           DT_join = as.character(DT_round),
           site = tolower(site)) %>%
    # These sites will be considered the same site for this workflow
    dplyr::mutate(site = ifelse(site == "rist", "tamasag",
                  ifelse(site == "elc", "boxelder", site))) %>%
    # Lastly, we swapped Boxelder's sonde out for Rist's late in 2022:
    dplyr::mutate(site = ifelse(site == "tamasag" & DT > lubridate::ymd("2022-09-20", tz = "MST") & DT < lubridate::ymd("2023-01-01", tz = "MST"), "boxelder", site)) %>%
    dplyr::distinct(.keep_all = TRUE)
  }

  if(network %in% c("virridy", "Virridy")){

    api_data <- api_data %>%
      # remove VuLink data
      dplyr::filter(!grepl("vulink", name, ignore.case = TRUE)) %>% # these will throw us errors if there is none?
      dplyr::select(-name) %>%
      # Convert UTC (as it is sent from HydroVU API) to MST:
      dplyr::mutate(DT = lubridate::as_datetime(timestamp, tz = "UTC")) %>%
      dplyr::mutate(DT = lubridate::with_tz(DT, tzone = "MST"),
                    DT_round = lubridate::round_date(DT, "15 minutes"),
                    DT_join = as.character(DT_round),
                    site = tolower(site)) %>%
      dplyr::distinct(.keep_all = TRUE)
    }

  return(api_data)

}

#' @title Network Check
#'
#' @description
#' This function performs a network check on a given data frame, flagging potential
#' issues in the data based on upstream and downstream sites.
#'
#' @param df A site-parameter data frame that has gone through the initial flagging
#' process.
#' @param network Whether the network check is happening across the Virridy sites or CSU sites.
#' @return A modified data frame flags that have been altered based on the network check.
#'
#' @examples
#' network_check(df = all_data_flagged$`archery-Actual Conductivity`)

network_check <- function(df, network = "csu") {

  df <- df

  # get the site name from the site column in the df of interest
  site_name <- unique(na.omit(df$site))
  # get the parameter from the parameter column in the df of interest
  parameter_name <- unique(na.omit(df$parameter))

  # vector of sites in the order that they are in spatially
  # some sites have some funkiness going on

  sites_order <- c("tamasag", # rist
                   "legacy",
                   "lincoln",
                   "timberline",
                   "prospect",
                   "boxelder", # elc
                   "archery",
                   "river bluffs")

  width_fun = ifelse(site_name == "tamasag", 17, # 2 hours before/after
              ifelse(site_name == "legacy", 17,
              ifelse(site_name == "lincoln", 17,
              ifelse(site_name == "timberline", 17,
              ifelse(site_name == "prospect", 17,
              ifelse(site_name == "boxelder", 17,
              ifelse(site_name == "archery", 17,
              ifelse(site_name == "river bluffs", 17, NA))))))))

  if(network %in% c("virridy", "Virridy")){

    sites_order <-  c("joei",
                      "cbri",
                      "chd",
                      "pfal",
                      "pbd",
                      "tamasag",
                      "legacy",
                      "lincoln",
                      "timberline",
                      #"springcreek",
                      "prospect",
                      "boxelder",
                      #boxcreek,"
                      "archery",
                      "river bluffs")

    width_fun = ifelse(site_name == "joei", 17, # 2 hours before/after
                ifelse(site_name == "cbri", 17,
                ifelse(site_name == "chd", 17,
                ifelse(site_name == "pfal", 17,
                ifelse(site_name == "pbd", 17,
                ifelse(site_name == "sfm", 17,
                ifelse(site_name == "lbea", 17,
                ifelse(site_name == "penn", 17,
                ifelse(site_name == "tamasag", 17,
                ifelse(site_name == "legacy", 17,
                ifelse(site_name == "lincoln", 17,
                ifelse(site_name == "timberline", 17,
                ifelse(site_name == "timberline virridy", 17,
                ifelse(site_name == "springcreek", 17,
                ifelse(site_name == "prospect", 17,
                ifelse(site_name == "prospect virridy", 17,
                ifelse(site_name == "boxelder", 17,
                ifelse(site_name == "boxcreek", 17,
                ifelse(site_name == "archery", 17,
                ifelse(site_name == "archery virridy", 17,
                ifelse(site_name == "river bluffs", 17, NA)))))))))))))))))))))

    if(site_name %in% c("penn", "sfm", "lbea")){

      sites_order <- c("penn",
        "sfm",
        "lbea")

    }

    if(site_name == "springcreek"){

      sites_order <- c("timberline virridy",
                       "springcreek",
                       "prospect virridy")

    }

    if(site_name == "boxcreek"){

      sites_order <- c("boxelder virridy",
                       "boxcreek",
                       "archery virridy")

    }

  }

  # determining the index for the site of interest.
  site_index <- which(sites_order == sites_order[grep(gsub(" virridy", "", site_name), sites_order, ignore.case = TRUE)])

  # Generating df name to pull from df_list list
  site_param <- paste0(site_name, "-", parameter_name)

  prev_site_df <- tibble(DT_round = NA)
  next_site_df <- tibble(DT_round = NA)

  tryCatch({
    previous_site <- paste0(sites_order[site_index-1],"-",parameter_name)
    prev_site_df <- intersensor_checks[[previous_site]] %>%
      select(DT_round, site_up = site, flag_up = flag) %>%
      data.table()},
    error = function(err) {
      cat("No upstream site.\n")})

  tryCatch({
    next_site <- paste0(sites_order[site_index+1],"-",parameter_name)
    next_site_df <- intersensor_checks[[next_site]] %>%
      select(DT_round, site_down = site, flag_down = flag) %>%
      data.table()},
    error = function(err) {
      cat("No downstream site.\n")})


  join <- df %>%
    left_join(., prev_site_df, by = "DT_round") %>%
    left_join(., next_site_df, by = "DT_round")

  if(!("flag_down" %in% colnames(join))) {join$flag_down <- NA}
  if(!("flag_up" %in% colnames(join))) {join$flag_up <- NA}
  if(!("site_down" %in% colnames(join))) {join$site_down <- NA}
  if(!("site_up" %in% colnames(join))) {join$site_up <- NA}


  # Define a function to check if a given 2-hour window has any instances of the same word
  check_2_hour_window_fail <- function(x) {
    sum(x) >= 1
  }

  df_test <- join %>%
    # No upstream/downstream flag = 0
    mutate(flag_binary = ifelse(
      (is.na(flag_up) | grepl("drift|DO interference|repeat|sonde not employed|frozen|unsubmerged|missing data|site visit|sv window|sensor malfunction|sonde burial|sensor biofouling|improper level cal|sonde moved", flag_up)) &
        (is.na(flag_down) | grepl("drift|DO interference|repeat|sonde not employed|frozen|unsubmerged|missing data|site visit|sv window|sensor malfunction|sonde burial|sensor biofouling|improper level cal|sonde moved", flag_down)), 0, 1)) %>%
    mutate(overlapping_flag = zoo::rollapply(flag_binary, width = width_fun, FUN = check_2_hour_window_fail, fill = NA, align = "center")) %>%
    # If there is a flag (flags associated with spikes in concentration or funkiness like that), and there is also a flag up/downstream at the same time (2 hour window) it is likely a real
    # WQ event and should therefore not be considered "poor" data:
    mutate(auto_cleaned_flag = ifelse(!is.na(flag) & !grepl("drift|DO interference|repeat|sonde not employed|frozen|unsubmerged|missing data|site visit|sv window|sensor malfunction|sonde burial|sensor biofouling|improper level cal|sonde moved", flag) & overlapping_flag == TRUE, NA, flag)) %>%
    select(-c(flag_up, flag_down, site_up, site_down, flag_binary, overlapping_flag))

  # df_test <- join %>%
  #   mutate(flag_binary = ifelse(#grepl("slope|suspect", flag) &
  #     (is.na(flag_up) | grepl("seasonal range|repeat|sonde not employed|missing data|site visit|sv window", flag_up)) &
  #       (is.na(flag_down) | grepl("seasonal range|repeat|sonde not employed|missing data|site visit|sv window", flag_down)), 0, 1)) %>%
  #   #arrange(timestamp) %>%
  #   mutate(overlapping_flag = zoo::rollapply(flag_binary, width = width_fun, FUN = check_2_hour_window_fail, fill = NA, align = "center")) %>%
  #   mutate(cleaner_flag = ifelse(!is.na(flag) & !grepl("seasonal range|repeat|sonde not employed|missing data|site visit|sv window", flag) & overlapping_flag == TRUE, NA, flag)) %>%
  #   select(-c(flag_up, flag_down, site_up, site_down, flag_binary, overlapping_flag))

  return(df_test)

}

# Function to check for package installation, then install (if necessary) and load libraries.
# Adapted from code developed by Caitlin Mothes, PhD.


# fill in with packages that need to be loaded:
# packages <- c('tidyverse',
#               'sf')

package_loader <- function(x) {
  if (x %in% installed.packages()) {
    library(x, character.only = TRUE)
  } else {
    install.packages(x)
    library(x, character.only = TRUE)
  }
}

# lapply(packages, package_loader)

# generate_flag_report <- function(df) {
#
#   # Extract the site and parameter from the df_name
#   site <- unique(na.omit(df$site))
#   parameter <- unique(na.omit(df$parameter))
#
#   list_of_flags <- c("sonde not employed", # add_field_flag()
#                      "site visit", # add_field_flag()
#                      "sv window", # add_field_flag()
#                      "sensor malfunction", # add_malfunction_flag()
#                      "outside of sensor specification range", # add_spec_flag()
#                      "outside of seasonal range", # add_seasonal_flag()
#                      "slope violation", # add_seasonal_flag()
#                      "outside sd range", # add_seasonal_flag()
#                      "repeated value", # add_repeat_flag()
#                      "missing data", # add_na_flag()
#                      "suspect data") # add_suspect_flag()
#
#   # check these
#   sans_na_flags <- "^(missing data|
#                       sonde not employed;\\nmissing data|
#                       missing data;\\nsuspect data|
#                       sonde not employed;\\nmissing data;\\nsuspect data|
#                       site visit;\\nmissing data;\\nsuspect data|
#                       sv window;\\nmissing data;\\nsuspect data)$"
#
#   # summarize total data points
#   total_observations <- df %>%
#     summarise(n_total = n_distinct(DT_round)) %>%
#     pull(n_total)
#
#   # summarize total data points sans missing data
#   total_observations_1 <- df %>%
#     # filter out when flag has only missing data or only sonde not employed and missing data
#     filter(!str_detect(flag, sans_na_flags)) %>%
#     summarise(n_total = n_distinct(DT_round)) %>%
#     pull(n_total)
#
#   # summarize total days
#   total_observations_dates <- df %>%
#     group_by(date = format(DT_round, "%m-%d-%Y")) %>%
#     summarize(n_total = nrow(date)) %>%
#     nrow()
#
#   # summarize total days sans missing data
#   total_observations_dates_1 <- df %>%
#     filter(!str_detect(flag, sans_na_flags)) %>%
#     group_by(date = format(DT_round, "%m-%d-%Y")) %>%
#     summarize(n_total = nrow(date)) %>%
#     nrow()
#
#   row_list <- list()
#   for (i in list_of_flags) {
#
#     # summarize flagged data points
#     flagged_observations <- df %>%
#       filter(str_detect(flag, i)) %>%
#       summarise(n_flag = n_distinct(DT_round)) %>%
#       pull(n_flag)
#     # summarize flagged data points
#     flagged_observations_1 <- df %>%
#       filter(!str_detect(flag, sans_na_flags)) %>%
#       filter(str_detect(flag, i)) %>%
#       summarise(n_flag = n_distinct(DT_round)) %>%
#       pull(n_flag)
#     # summarize percent data points that are flagged
#     percent_flagged <- flagged_observations/total_observations
#     # summarize percent data points that are flagged sans missing data
#     percent_flagged_1 <- flagged_observations_1/total_observations_1
#
#     # summarize flagged days
#     flagged_observations_dates <- df %>%
#       filter(str_detect(flag, i)) %>%
#       group_by(date = format(DT_round, "%m-%d-%Y")) %>%
#       summarize(n_total = nrow(date)) %>%
#       nrow()
#     # summarize flagged days
#     flagged_observations_dates_1 <- df %>%
#       filter(!str_detect(flag, sans_na_flags)) %>%
#       filter(str_detect(flag, i)) %>%
#       group_by(date = format(DT_round, "%m-%d-%Y")) %>%
#       summarize(n_total = nrow(date)) %>%
#       nrow()
#     # summarize percent days that are flagged
#     percent_flagged_dates <- flagged_observations_dates/total_observations_dates
#     # summarize percent days that are flagged
#     percent_flagged_dates_1 <- flagged_observations_dates_1/total_observations_dates_1
#
#     # creating a row with the information
#     calculated_values <- tibble(
#       # metadata
#       site = site,
#       parameter = parameter,
#       flag = i,
#       # data points
#       data_points_flagged_percentage = format(round(percent_flagged * 100, 2), nsmall = 2),
#       data_points_flagged = flagged_observations,
#       total_data_points = total_observations,
#       data_points_flagged_percentage_sans_na = format(round(percent_flagged_1 * 100, 2), nsmall = 2),
#       data_points_flagged_sans_na = flagged_observations_1,
#       total_data_points_sans_na = total_observations_1,
#       # dates
#       dates_flagged_percentage = format(round(percent_flagged_dates * 100, 2), nsmall=2),
#       dates_flagged = flagged_observations_dates,
#       total_dates = total_observations_dates,
#       dates_flagged_percentage_sans_na = format(round(percent_flagged_dates_1 * 100, 2), nsmall=2),
#       dates_flagged_sans_na = flagged_observations_dates_1,
#       total_dates_sans_na = total_observations_dates_1
#       )
#
#     row_list[[i]] <- calculated_values
#   }
#
#   #calculated_df <- bind_cols(row_list)
#   return(bind_rows(row_list))
#
# }


photo_plotter <- function(index = 200, output_folder){

  param_1_label <- paste0(parameters[1], " (", param_unit[1],")")
  param_2_label <- paste0(parameters[2], " (", param_unit[2],")")

# get the adjustment, breaks etc so axes look nice
  param_1_max <- max(wq_tl[[parameters[1]]], na.rm = T)
  param_1_max_int <- as.integer(max(wq_tl[[parameters[1]]], na.rm = T))+1

  adjustment <- max(wq_tl[[parameters[2]]], na.rm = T) / param_1_max

  brk <- ifelse(param_1_max_int < 6, 1,2)

# Bounds for the y-axis

  lower_bound <- 0
  upper_bound <-  max(wq_tl[[parameters[2]]], na.rm = T)


  #This is the index of the image for the background for the individual photo
  simul = wq_tl[index,]
  # this is all the data before the image (ie previous photos from the timelapse)
  upto = wq_tl[1:index,]
  #read the image for the background
  photo_bg <- readJPEG(simul$filename)
  #create an individual image
  back <- ggplot() +
    #plot the photo
    annotation_custom(rasterGrob(photo_bg,
                                 width = unit(1,"npc"),
                                 height = unit(1,"npc")),
                      -Inf, Inf, -Inf, Inf)

  #plot the data for this image (includes all the preivous data)
  inset <- ggplot() +
    geom_ribbon(data = upto, aes(x = DT_round,
#to do: need to figure out adjustment factors for each parameter to plot correctly
                                 y = .data[[parameters[1]]] *adjustment,
                                 ymin = 0,
                                 ymax = .data[[parameters[1]]] * adjustment),
                color = "white",
                fill = "white",
                #linetype = "dash",
                alpha = 0.75)+
    geom_path(data = upto, aes(x = DT_round, y = .data[[parameters[2]]]),
              color = "#F34646", size=2) +
    geom_point(data = simul, aes(x = DT_round, y = .data[[parameters[2]]]),
               color = "#F34646")+
    # ylim(min(wq_tl[[parameters[2]]], na.rm = T),
    #      max(wq_tl[[parameters[2]]], na.rm = T))+
    xlim(min(wq_tl$DT_round, na.rm = T),
         max(wq_tl$DT_round, na.rm = T))+
    scale_y_continuous(name = param_2_label, limits = c(0, upper_bound),
                       sec.axis = sec_axis(trans = ~./adjustment ,
                                           name = param_1_label,
                                           breaks = seq(0,param_1_max_int,brk)))+
    dark_theme_light(base_size = 10) +
    theme(axis.title.y.right = element_text(color="white"),axis.text.y.right = element_text(color = "white"), axis.title.y.left = element_text(color="#F34646"), axis.text.y.left = element_text(color = "#F34646")) +
    theme(
      panel.background = element_rect(fill = "transparent"),
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.background = element_rect(fill = "transparent"),
      legend.box.background = element_rect(fill = "transparent"),
      panel.border = element_blank(),
      axis.line = element_line(color = 'gray80'),
    ) +
    ylab(param_2_label) +
    xlab('')

  gp1 <- back +
    inset_element(inset,
                  left = 0.16,
                  bottom = 0.15,
                  right = 0.9,
                  top = 0.6)
  #print(gp1)
  ggsave(paste0(folder_path, index, ".png"),
         width = 1920,
         height = 1080,
         units = "px")
  #dev.copy(gp1,paste0('data/timelapse_photos/vid_image/', max(upto$rowid), ".png"))
  #dev.off()
}

#test functions
#photo_plotter(2, "data/timelapse_photos/sjs_test/")
#map(1:nrow(wq_tl), ~photo_plotter(.x, "data/timelapse_photos/sjs_test/"))

plotter <- function(site_name = "boxelder",
                    start_date = "2023-03-04",
                    end_date = "2023-07-07",
                    title = "wack data title",
                    parameter = "all") {

  # Plotting options: all, Specific Conductivity, Turbidity, Chl-a Fluorescence, ORP, Depth, pH, DO

  sites_order <- c("tamasag", # rist
                   "legacy",
                   "lincoln",
                   "timberline",
                   "prospect",
                   "boxelder", # elc
                   "archery",
                   "river bluffs")

  # determining the index for the site of interest.
  site_index <- which(sites_order == site_name)
  previous_site <- sites_order[site_index-1]
  next_site <- sites_order[site_index+1]

  # Pull in WQ data from Manners Bridge (in the canyon):
  contrail <- list.files("data/context_data/contrail", full.names = TRUE) %>%
    map_dfr(~fread(.)) %>%
    mutate(date = as_date(ymd_hms(Reading)),
           Reading = floor_date(ymd_hms(Reading), "hour")) %>%
    filter(date >= start_date & date <= end_date)

  # Pull in a list of precip monitoring sites in each sonde's watershed:
  rain_list <- fread("data/context_data/site_link.csv") %>%
    filter(site == site_name) %>% pull(`Sensor Name`)

  # Pull in raw precip data
  rain_data <- fread("data/context_data/fc_rain.csv") %>%
    filter(`Sensor Name` %in% rain_list) %>%
    mutate(date = as_date(mdy_hms(Timestamp)),
           DT = floor_date(mdy_hms(Timestamp), unit = "hour")) %>%
    filter(date >= start_date & date <= end_date) %>%
    group_by(DT) %>%
    summarize(total_ws = sum(`Incremental Rainfall (in)`, na.rm = TRUE))

  # Grab the largest amount of rain across the sonde's watershed at an
  # hourly timestep:
  rain_stat <- fread("data/context_data/fc_rain.csv") %>%
    filter(`Sensor Name` %in% rain_list) %>%
    mutate(date = as_date(mdy_hms(Timestamp)),
           DT = floor_date(mdy_hms(Timestamp), unit = "hour")) %>%
    group_by(DT) %>%
    summarize(total_ws = sum(`Incremental Rainfall (in)`, na.rm = TRUE)) %>%
    filter(total_ws == max(total_ws))

  # Pull in flagged PSN WQ data:
  df <- fread('data/flagged/all_data_flagged_plotter.csv') %>%
    filter(!grepl("sonde not employed", final_cleaned_flag)) %>%
    mutate(date = as_date(DT_round)) %>%
    # filter to dates of interest:
    dplyr::filter((date) >= as_date(start_date) & (date) <= as_date(end_date))

  # site of interest data:
  site_df <- df %>%
    filter(!is.na(mean)) %>%
    filter(site == site_name)

  # downstream of site of interest (if there is one):
  downstream_df <- site_df
  try(downstream_df <- df %>%
        filter(!grepl("site visit|sv window", final_cleaned_flag)) %>%
        filter(site == next_site))

  # upstream of site of interest (if there is one):
  upstream_df <- site_df
  try(upstream_df <- df %>%
        filter(!grepl("site visit|sv window", final_cleaned_flag)) %>%
        filter(site == previous_site))

# All the plotz
  depth <- ggplot() +
    geom_line(data = filter(contrail, Unit == "ft"), aes(Reading, Value * 0.3048), color = "lightblue") +
    geom_line(data = filter(downstream_df, parameter == "Depth"), aes(DT_round, mean), cex = 0.2, alpha = 1, color = "lightgrey") +
    geom_line(data = filter(upstream_df, parameter == "Depth"), aes(DT_round, mean), cex = 0.2, alpha = 1, color = "darkgrey") +
    geom_line(data = filter(site_df, parameter == "Depth"), aes(DT_round, mean), color = "black") +
    geom_point(data = filter(site_df, parameter == "Depth" & grepl("site visit", final_cleaned_flag)), aes(DT_round, mean), color = "#E70870", cex = 1.5) +
    geom_point(data = filter(site_df, parameter == "Depth" & grepl("sonde moved", final_cleaned_flag)), aes(DT_round, mean), color = "#FFCA3A", cex = 0.75) +
    theme_bw() +
    ylim(min = min(filter(site_df, parameter == "Depth")$mean), max = max(filter(site_df, parameter == "Depth")$mean)) +
    ylab("Depth m") + xlab("") + ggtitle(title)

  temp <- ggplot() +
    geom_line(data = filter(contrail, Unit == "C"), aes(Reading, Value), color = "lightblue") +
    geom_line(data = filter(downstream_df, parameter == "Temperature"), aes(DT_round, mean), cex = 0.2, alpha = 1, color = "lightgrey") +
    geom_line(data = filter(upstream_df, parameter == "Temperature"), aes(DT_round, mean), cex = 0.2, alpha = 1, color = "darkgrey") +
    geom_line(data = filter(site_df, parameter == "Temperature"), aes(DT_round, mean), color = "black") +
    geom_point(data = filter(site_df, parameter == "Temperature" & grepl("site visit|sv window", final_cleaned_flag)), aes(DT_round, mean), color = "#E70870") +
    theme_bw() +
    ylim(min = min(filter(site_df, parameter == "Temperature")$mean), max = max(filter(site_df, parameter == "Temperature")$mean)) +
    ylab("Temperature C") + xlab("") + ggtitle("")

  ph <- ggplot() +
    geom_line(data = filter(contrail, Unit == "pH"), aes(Reading, Value), color = "lightblue") +
    geom_line(data = filter(downstream_df, parameter == "pH"), aes(DT_round, mean), cex = 0.2, alpha = 1, color = "lightgrey") +
    geom_line(data = filter(upstream_df, parameter == "pH"), aes(DT_round, mean), cex = 0.2, alpha = 1, color = "darkgrey") +
    geom_line(data = filter(site_df, parameter == "pH"), aes(DT_round, mean), color = "black") +
    geom_point(data = filter(site_df, parameter == "pH" & grepl("site visit|sv window", final_cleaned_flag)), aes(DT_round, mean), color = "#E70870") +
    theme_bw() +
    ylim(min = min(filter(site_df, parameter == "pH")$mean), max = max(filter(site_df, parameter == "pH")$mean)) +
    xlab("") + ylab("pH") + ggtitle("")

  orp <- ggplot() +
    geom_line(data = filter(downstream_df, parameter == "ORP"), aes(DT_round, mean), cex = 0.2, alpha = 1, color = "lightgrey") +
    geom_line(data = filter(upstream_df, parameter == "ORP"), aes(DT_round, mean), cex = 0.2, alpha = 1, color = "darkgrey") +
    geom_line(data = filter(site_df, parameter == "ORP"), aes(DT_round, mean), color = "black") +
    geom_point(data = filter(site_df, parameter == "ORP" & grepl("site visit|sv window", final_cleaned_flag)), aes(DT_round, mean), color = "#E70870") +
    theme_bw() +
    ylim(min = min(filter(site_df, parameter == "ORP")$mean), max = max(filter(site_df, parameter == "ORP")$mean)) +
    ylab("ORP") + xlab("") + ggtitle("")

  spc <- ggplot() +
    geom_line(data = filter(contrail, Unit == "uS/cm"), aes(Reading, Value), color = "lightblue") +
    geom_line(data = filter(downstream_df, parameter == "Specific Conductivity"), aes(DT_round, mean), cex = 0.2, alpha = 1, cex = 0.2, alpha = 1, color = "lightgrey") +
    geom_line(data = filter(upstream_df, parameter == "Specific Conductivity"), aes(DT_round, mean), cex = 0.2, alpha = 1, cex = 0.2, alpha = 1, color = "darkgrey") +
    geom_line(data = filter(site_df, parameter == "Specific Conductivity"), aes(DT_round, mean), color = "black", cex = 0.8) +
    geom_point(data = filter(site_df, parameter == "Specific Conductivity" & grepl("site visit|sv window", final_cleaned_flag)), aes(DT_round, mean), color = "#E70870", cex = 1.5) +
    theme_bw() +
    ylim(min = min(filter(site_df, parameter == "Specific Conductivity")$mean) - 1, max = max(filter(site_df, parameter == "Specific Conductivity")$mean) + 1) +
    ylab("SpC uS/cm") + xlab("") + ggtitle("")

  do <- ggplot() +
    geom_line(data = filter(contrail, Unit == "mg/L"), aes(Reading, Value), color = "lightblue") +
    geom_line(data = filter(downstream_df, parameter == "DO"), aes(DT_round, mean), cex = 0.2, alpha = 1, color = "lightgrey") +
    geom_line(data = filter(upstream_df, parameter == "DO"), aes(DT_round, mean), cex = 0.2, alpha = 1, color = "darkgrey") +
    geom_line(data = filter(site_df, parameter == "DO"), aes(DT_round, mean), color = "black") +
    geom_point(data = filter(site_df, parameter == "DO" & grepl("site visit|sv window", final_cleaned_flag)), aes(DT_round, mean), color = "#E70870") +
    theme_bw() +
    ylim(min = min(filter(site_df, parameter == "DO")$mean) - 1, max = max(filter(site_df, parameter == "DO")$mean) + 1) +
    ylab("DO mg/L") + xlab("") + ggtitle("")

  # Handle plotting for sondes with turbidity vs. chl-a data:

  turb <- ggplot() +
    geom_line(data = filter(contrail, Unit == "ntu"), aes(Reading, Value), color = "lightblue") +
    geom_line(data = filter(downstream_df, parameter == "Turbidity"), aes(DT_round, mean), cex = 0.2, alpha = 1, color = "lightgrey") +
    geom_line(data = filter(upstream_df, parameter == "Turbidity"), aes(DT_round, mean), cex = 0.2, alpha = 1, color = "darkgrey") +
    geom_line(data = filter(site_df, parameter == "Turbidity"), aes(DT_round, mean), color = "black") +
    geom_point(data = filter(site_df, parameter == "Turbidity" & grepl("site visit|sv window", final_cleaned_flag)), aes(DT_round, mean), color = "#E70870") +
    theme_bw() +
    ylim(min = min(filter(site_df, parameter == "Turbidity")$mean) + 10, max = max(filter(site_df, parameter == "Turbidity")$mean) + 10) +
    ylab("Turbidity NTU") + xlab("") + ggtitle("")

  if(nrow(filter(site_df, parameter == "Chl-a Fluorescence")) > 0){
    chla <- ggplot() +
      geom_line(data = filter(downstream_df, parameter == "Chl-a Fluorescence"), aes(DT_round, mean), cex = 0.2, alpha = 1, color = "lightgrey") +
      geom_line(data = filter(upstream_df, parameter == "Chl-a Fluorescence"), aes(DT_round, mean), cex = 0.2, alpha = 1, color = "darkgrey") +
      geom_line(data = filter(site_df, parameter == "Chl-a Fluorescence"), aes(DT_round, mean), color = "black") +
      geom_point(data = filter(site_df, parameter == "Chl-a Fluorescence" & grepl("site visit|sv window", final_cleaned_flag)), aes(DT_round, mean), color = "#E70870") +
      theme_bw() +
      ylab("Chl-a RFU") + xlab("") + ggtitle("")
  }

  if(nrow(filter(site_df, parameter == "Chl-a Fluorescence")) == 0){ chla <- NULL }

  # rain proxy plot. this just plots the sum of at an hourly timestep of
  # rain at all sites in the watershed. not rooted in anything scientific,
  # but just shows in a crude way whether or not there was rain happening
  # at any given time (and at any given point) during the time of interest.
  rain <- ggplot(data = rain_data) +
    geom_col(aes(x = DT,
                 y = total_ws), color = "#002EA3") +
    ylab("Rain (PROXY)") + xlab("") +
    ylim(min = 0, max = rain_stat$total_ws) +
    theme_bw()

  # handling what gets plotted based on user input:

  if(parameter %in% c("All", "all", "all parameters", "all params")){
    if(is.null(chla)) {

      plot <- ggarrange(depth, temp, ph, orp, spc, do, turb,
                         rain, nrow = 4, ncol = 2)
    }

    if(!is.null(chla)){
      plot <- ggarrange(depth, temp, ph, orp, spc, do, chla,
                         rain, nrow = 4, ncol = 2)
    }
  }

  if(parameter %in% c("pH", "ph")){
    plot <- ggarrange(depth, temp, ph, rain,
                      ncol = 1, nrow = 4)
  }

  if(parameter %in% c("orp", "ORP")){
    plot <- ggarrange(depth, temp, orp, rain,
                       ncol = 1, nrow = 4)
  }

  if(parameter %in% c("do","DO")){
    plot <- ggarrange(depth, temp, do, rain,
                      ncol = 1, nrow = 4)
  }

  if(parameter %in% c("turbidity", "turb", "Turb", "Turbidity")){
    plot <- ggarrange(depth, temp, turb, rain,
                      ncol = 1, nrow = 4)
  }

  if(parameter %in% c("Specific Conductivity", "conductivity", "specific conductivity", "sc", "spc", "SC","SpC")){
    plot <- ggarrange(depth, temp, spc, rain,
                      ncol = 1, nrow = 4)
  }

  if(parameter %in% c("chla", "Chla", "chl-a", "chlorophyll a", "Chl-a Fluorescence")){
    plot <- ggarrange(depth, temp, chla, rain,
                      ncol = 1, nrow = 4)
  }

  return(plot)
}

retrieve_relevant_data_name <- function(df_name_arg, year_week_arg = NULL, interval_arg = NULL) {
  if(is.null(interval_arg)){
    if (df_name_arg %in% names(verified_data) & any(year_week_arg %in% verified_data[[df_name_arg]]$y_w)) {
      return("verified_data")
    }
    if (df_name_arg %in% names(intermediary_data) & any(year_week_arg %in% intermediary_data[[df_name_arg]]$y_w)) {
      return("intermediary_data")
    }
    if (df_name_arg %in% names(all_data) & any(year_week_arg %in% all_data[[df_name_arg]]$y_w)) {
      return("all_data")
    }
  }
  if(is.null(year_week_arg)){
    if (df_name_arg %in% names(verified_data) & any(verified_data[[df_name_arg]]$DT_round %within% interval_arg)) {
      return("verified_data")
    }
    if (df_name_arg %in% names(intermediary_data) & any(intermediary_data[[df_name_arg]]$DT_round %within% interval_arg)) {
      return("intermediary_data")
    }
    if (df_name_arg %in% names(all_data) & any(all_data[[df_name_arg]]$DT_round %within% interval_arg)) {
      return("all_data")
    }
  }
}

## Water Sampling Data:

# Goal:
#Save data in the correct format for RMRS spreadsheet
#Save all water sampling probe values in a spreadsheet
# To get the RMRS style data for a specfic date of sampling,
# Input the date of interest in sampling_spreadsheet_creator
#sampling_spreadsheet_creator(date_oi = "2023-11-17")

# To get all the water sampling data and save to CSV in sampling notes
# This also returns the df sampling_notes in case you want to review in R
#sampling_spreadsheet_creator(all_dates = TRUE)


sampling_spreadsheet_creator <- function(date_oi = "2023-10-16", all_dates = FALSE ){

  #source clean mwater script for all notes cleaned

  source("src/load_mWater_notes.R")

  #pull in site meta data
  site_meta <- read_csv("data/metadata/water_sampling_sites.csv",show_col_types = FALSE)%>%
    select(site = site_code, Site_Name, site_label_rmrs)
  # sort for sites in upper network (ie. acronyms rather than street names)
  upper_sites <- read_csv("data/metadata/water_sampling_sites.csv",show_col_types = FALSE)%>%
    filter(watershed != "CLP  Mainstem-Fort Collins")%>%
    #this is to help match with user input
    mutate(site_code = tolower(site_code))

  # create df of all water samples and save DT, handheld probe and chla volume data
  sampling_notes <- load_mWater_notes()%>%
    filter(grepl("Sampling",visit_type))%>%
    mutate(all_pics_taken = case_when(!is.na(downstream_pic)&!is.na(upstream_pic)&!is.na(clarity)&!is.na(filter_pic) ~ TRUE, TRUE ~ FALSE),
           #correct names if it is in our upper sites (acronyms)
           site = ifelse(site %in% upper_sites$site_code, toupper(site), site),
           DT_round = round_date(start_DT, "15 minutes"))%>%
    select(site,crew, DT_round, date, time = start_time_mst, sample_collected, chla_volume_ml, vol_filtered_blank_dup, do_mgl, cond_ms_cm, temp_c, visit_comments, all_pics_taken, q_cfs)

  # Distinguish BLANK and DUPLICATE values
  blanks_dups <- sampling_notes %>%
    #find all values that have blank or dup
    filter(grepl("DUPLICATE|BLANK", sample_collected)) %>%
    # change sample collected to match BLANK/DUP
    mutate(sample_collected = ifelse(grepl("DUPLICATE", sample_collected), "DUPLICATE", "BLANK"),
           # Volume filtered blank dup becomes chla volume
           chla_volume_ml = vol_filtered_blank_dup,
           #drop vol_filtered_blank/dup
           vol_filtered_blank_dup = NULL)

  # Add blank and duplicate values back to main
  sampling_notes <- sampling_notes%>%
    #get rid of blank/dup in sample collcected
    mutate(sample_collected = gsub("DUPLICATE|BLANK| |,", "", sample_collected),
           #drop vol_filtered_blank/dup
           vol_filtered_blank_dup = NULL)%>%
    #bring in blank and dup rows
    rbind(blanks_dups)%>%
    #arrange by datetime and site (Blanks and dups go second)
    arrange(DT_round, site)%>%
    # join with RMRS friendly metadata
    left_join(site_meta, by = "site")



 if (all_dates == TRUE) {

   sampling_notes_output <- sampling_notes%>%
     # select only the needed columns, saved in the correct order and fix column names
     select(site_code = site, Date = date, SampleType = sample_collected, time_mst = time,chla_volume_ml,  do_mgl, cond_ms_cm, temp_c, visit_comments)

   # write to csv
   write_csv(x = sampling_notes_output, file = paste0("data/sampling_notes/all_samples_as_of_",as.character(Sys.Date()),".csv" ))
 }else{
   #grab desired date
   date_oi_clean <- as.Date(date_oi, tz = "America/Denver")
   # filter sampling notes df by desired date
   samples_of_day <- filter(sampling_notes,date == date_oi_clean )%>%
     #match date to RMRS style
     mutate(Date = format(date, "%d-%b-%y"),
            #create SiteDescr column for RMRS sheet
            SiteDescr = paste0(site, "_",format(date, "%m%d%y")))%>%
     # select only the needed columns, saved in the correct order and fix column names
     select(Site = site ,SiteDescr, SiteLabel = site_label_rmrs , Date, SampleType = sample_collected, q_cfs,
            time, do_mgl, cond_ms_cm, temp_c, notes = visit_comments  )
   # write to csv
   write_csv(x = samples_of_day, file = paste0("data/sampling_notes/samples_of_",date_oi,".csv" ))
 }
}


save_intermediary_dir <- function(updated_df, df_name) {
  # if updated_site_param_df is from all_data, but there is already intermediary data for it, refuse to move on
  if (all(!updated_df$is_verified) & df_name %in% list.files(intermediary_path)) {
    while (TRUE) {
      cat("!!!YOU ARE ABOUT TO OVERWRITE INTERMEDIARY DATA WITH PRE-VERIFICATION DATA!!!")
      user_input <- readline(prompt = "Are you sure you want to continue with this action? (yes/no):")
      user_input <- tolower(user_input)

      if (user_input == "yes") {
        saveRDS(updated_df, paste0(intermediary_path, df_name))
        break
      }

      if (user_input == "no") {
        working_data <<- set_names(map(list.files(intermediary_path, full.names = TRUE), readRDS), list.files(intermediary_path))
        updated_site_param_df <<- working_data[[site_param_name]] # working dir df
        saveRDS(updated_site_param_df, paste0(intermediary_path, df_name))
        break
      }

      cat("Invalid input. Please enter one of the options from the following:\n\n")
      cat("yes\n")
      cat("no\n")
    }
  }

  # if updated_site_param_df is from verified data, refuse to move on
  if (df_name %in% list.files(verified_path)) { # test to see if this works
    while (TRUE) {
      cat("!!!YOU ARE ABOUT TO WORK ON DATA THAT HAS BEEN FINALIZED!!!")
      user_input <- readline(prompt = "Are you sure you want to continue with this action? (yes/no):")
      user_input <- tolower(user_input)

      if (user_input == "yes") {
        saveRDS(updated_df, paste0(intermediary_path, df_name))
        break
      }

      if (user_input == "no") {
        cat("Please update your `site` and `parameter` objects to data that has not been finalized.")
        break
      }

      cat("Invalid input. Please enter one of the options from the following:\n\n")
      cat("yes\n")
      cat("no\n")
    }
  }

  saveRDS(updated_df, paste0(intermediary_path, df_name))
}

# Generate plots with both weekly and daily flagged data.
# This function will generate a list of plots with both weekly and daily flagged data.
# @param site_arg A string of the site name.
# @param parameter_arg A string of the parameter name.
# @param flag_arg A string of the flag name.
# @return A list of plots with both weekly and daily flagged data.
# @examples
# stack_flag_plots(site_arg = "archery", parameter_arg = "Actual Conductivity", flag_arg = "outside of Actual Conductivity sensor specification range")
# stack_flag_plots(site_arg = "boxelder", parameter_arg = "Temperature", flag_arg = "outside of Temperature sensor specification range")

# stack_flag_plots <- function(site_arg, parameter_arg, flag_arg, df_list) {
#   # Call on the weekly and daily functions and fill their args with this
#   # functions args
#   weekly_plot_list <- generate_weekly_flag_plots(site_arg = site_arg, parameter_arg = parameter_arg, flag_arg = flag_arg, df_list = df_list)
#   daily_plot_list <- generate_daily_flag_plots(site_arg = site_arg, parameter_arg = parameter_arg, flag_arg = flag_arg, df_list = df_list)
#
#   # These two functions should always return the same amount of plots, so we can
#   # use map2() to stack them with combine_plots()
#   weekly_daily_plots <- map2(.x = weekly_plot_list, .y = daily_plot_list, ~ggarrange(.x, .y, nrow = 2, ncol = 1, heights = 5, widths = 12))
#   return(weekly_daily_plots)
# }

#' @title Summarize site parameter data from the API and field notes data frames.
#' @description
#' A short description...
#' @param site_arg A site name.
#' @param parameter_arg A parameter name.
#' @param api_data A dataframe with the munged API data.
#' @param notes The munged field notes
#' @return A dataframe with summary statistics for a given site parameter data frame.
#' @examples
# summarize_site_param(site_arg = "archery", parameter_arg = "Actual Conductivity", api_data = incoming_data_collated_csvs)
# summarize_site_param(site_arg = "boxelder", parameter_arg = "Temperature", api_data = incoming_data_collated_csvs)

summarize_site_param_full <- function(site_arg, parameter_arg, api_data, notes = field_notes) {

  # filter deployment records for the full join
  site_field_notes <- notes %>%
    filter(grepl(paste(unlist(str_split(site_arg, " ")), collapse = "|"), site, ignore.case = TRUE))

  # filtering the data and generating results
  summary <- tryCatch({
    api_data %>%
      # subset to single site-parameter combo:
      dplyr::filter(site == site_arg & parameter == parameter_arg) %>%
      # safety step of removing any erroneous dupes
      dplyr::distinct() %>%
      # across each 15 timestep, get the average value, spread, and count of obs
      dplyr::group_by(DT_round, site, parameter) %>%
      dplyr::summarize(mean = as.numeric(mean(value, na.rm = T)),
                       diff = abs(min(value, na.rm = T) - max(value, na.rm = T)),
                       n_obs = n()) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(DT_round) %>%
      # pad the dataset so that all 15-min time stamps are present
      padr::pad(by = "DT_round", interval = "15 min") %>%
      # add a DT_join column to join field notes to (make DT_round character string, so no
      # funky DT issues occur during the join):
      dplyr::mutate(DT_join = as.character(DT_round),
                    site = site_arg,
                    parameter = parameter_arg,
                    flag = NA) %>% # add "flag" column for future processing
      # join our tidied data frame with our field notes data:
      dplyr::left_join(dplyr::filter(dplyr::select(site_field_notes, sonde_employed, last_site_visit, DT_join, visit_comments, sensor_malfunction, cals_performed)),
                       by = c('DT_join')) %>%
      # make sure DT_join is still correct:
      dplyr::mutate(DT_round = lubridate::as_datetime(DT_join, tz = "MST")) %>%
      # Use fill() to determine when sonde was in the field, and when the last site visit was
      # Necessary step for FULL dataset only (this step occurs in combine_hist_inc_data.R for auto)
      tidyr::fill(c(sonde_employed, last_site_visit, sensor_malfunction)) %>%
      dplyr::mutate(sonde_employed = ifelse(is.na(sonde_employed), 0, sonde_employed))
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

#' @title Summarize site parameter data from the HydroVu API and field notes data frames.
#'
#' @description
#' A function that summarizes site parameter data from the API and field notes data frames.
#'
#' @param site_arg A site name.
#'
#' @param parameter_arg A parameter name.
#'
#' @param api_data A dataframe with the munged API data.
#'
#' @param notes The munged field notes
#'
#' @param require An upstream dependency (i.e., alternative to `{necessary object}`).
#'
#' @return A dataframe with summary statistics and field notes for a given site parameter data frame.
#'
#' @examples
# summarize_site_param(site_arg = "archery", parameter_arg = "Actual Conductivity", api_data = incoming_data_collated_csvs)
# summarize_site_param(site_arg = "boxelder", parameter_arg = "Temperature", api_data = incoming_data_collated_csvs)

summarize_site_param <- function(site_arg, parameter_arg, api_data, notes, require = NULL) {

  # filter deployment records for the full join
  site_field_notes <- notes %>%
    dplyr::filter(site == site_arg)

  # filtering the data and generating results
  summary <- tryCatch({
    api_data %>%
      # summarize the data
      dplyr::filter(site == site_arg & parameter == parameter_arg) %>%
      dplyr::distinct() %>%
      dplyr::group_by(DT_round) %>%
      dplyr::summarize(mean = as.numeric(mean(value, na.rm = T)),
                       diff = abs(min(value, na.rm = T) - max(value, na.rm = T)),
                       n_obs = n()) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(DT_round) %>%
      # pad the dataset so that all 15-min timestamps are present
      padr::pad(by = "DT_round", interval = "15 min") %>%
      # join the field notes
      dplyr::mutate(DT_join = as.character(DT_round),
                    site = site_arg,
                    parameter = parameter_arg,
                    flag = NA) %>%
      dplyr::left_join(dplyr::filter(dplyr::select(site_field_notes, sonde_employed,
                                                   last_site_visit, DT_join, site,
                                                   visit_comments, sensor_malfunction,
                                                   cals_performed)),
                       by = c('DT_join', 'site')) %>%
      # join the metaparameter dfs
      dplyr::left_join(site_metaparams_list[[site_arg]], by = "DT_join") %>%
      dplyr::mutate(DT_round = lubridate::as_datetime(DT_join, tz = "MST"))
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

#' @title Update Historical Flag List
#'
#' @description
#' A function that updates the historical flag list with new flagged data.
#'
#' @param new_flagged_data A list of data frames that have been flagged.
#'
#' @param historical_flagged_data A list of data frames that have been flagged.
#'
#' @return A list of data frames that have been updated with new flagged data.
#'
#' @examples
#' update_historical_flag_list(new_flagged_data = all_data_flagged$`archery-Actual Conductivity`, historical_flagged_data = all_data_flagged$`archery-Actual Conductivity`)
#' update_historical_flag_list(new_flagged_data = all_data_flagged$`boxelder-Temperature`, historical_flagged_data = all_data_flagged$`boxelder-Temperature`)

update_historical_flag_list <- function(new_flagged_data, historical_flagged_data){

  # Get the matching index names
  matching_indexes <- intersect(names(new_flagged_data), names(historical_flagged_data))

  # bind new_flagged_data and historical_flagged_data together
  updated_historical_flag_list <- map(matching_indexes, function(index) {

    old <- historical_flagged_data[[index]] %>%
      filter(DT_round < ymd_hms(max(DT_round) - hours(24), tz = "MST")) %>% # this is the antijoin step (but i thought we were doing 24 hours?) -jd changing this to 24 removes the duplicate problem
      mutate(last_site_visit = force_tz(last_site_visit, tzone = "MST"))

    bind_rows(old, new_flagged_data[[index]]) %>%
      arrange(DT_round) %>%
      select(-historical) # if we are removing this then I think we should remove that step from combined_data -jd
  }) %>%
    set_names(matching_indexes) %>%
    keep(~ !is.null(.))

  # need to make sure that when we are doing this we are not getting rid of dfs in the RDS list, still needs to be checked -jd
    # if we are matching up the historical list(88) with the incoming list(72) and only returning the matches then we will miss
    # be returning a list that is being written that is shorter than the list that we started with... Definitely needs to be fixed.

  return(updated_historical_flag_list)

}

# error 1: tz discrepancy at the join between incoming data and HFD
  # this messes up the future pulls because the start times df will be wrong

# time_diffs <- diff(test_data$DT_round)
#
# time_diffs # these should all be 15 and there is one that is not (where we joined incoming data to the HFD)

verify_flag_data <- function(df_list_arg,
                             site_arg,
                             parameter_arg, # do we want this to be able to be null eventually?
                             flag_arg = NULL,
                             weekly_plot_object) {

  # set layer
  layer <- "weekly"

  #quit

  if (!QUIT) {
    # verification decision lists ----
    simple_decision_list <- c("PASS ALL", "FAIL NONE",
                              "PASS VALID", "FAIL FLAGGED",
                              "PASS FLAGGED", "FAIL VALID",
                              "PASS NONE", "FAIL ALL",
                              "SKIP", "QUIT")

    inspect_decision_list <- c("INSPECT ALL", "INSPECT VALID",
                               "INSPECT FLAGGED", "INSPECT SOME")

    # initialize site_param_df ----
    site_param_df <- df_list_arg[[paste0(site_arg, "-", parameter_arg)]]

    # get data from initial weekly plot object ----
    df_weekly_data <- ggplot_build(weekly_plot_object)$plot$data %>%
      filter(site == site_arg) %>%
      select(intersect(names(.), names(site_param_df)))

    df_weekly_name <- ggplot_build(weekly_plot_object)$plot$labels$title

    # initialize list to contain updated dfs OUTPUT, VERY IMPORTANT ----
    altered_df_list <- list()

    # initialize week list ----
    weekday_list <- c(1,2,3,4,5,6,7)

    # Make sure that non-removable flags are still in place ----
    non_removable_flags <- paste0(c("sensor biofouling","sonde burial",
                                    "sensor malfunction","frozen","missing data",
                                    "sonde not employed","sonde unsubmerged",
                                    "depth calibration malfunction"),
                                  collapse = "|")

    # Which data are you working on: ----
    cat("Currently verifying: ", df_weekly_name, "data. \n")

    print(weekly_plot_object)

    # Get weekly decision(s) ----
    verification_decision <- get_verification_decision(
      paste("Would you like to (pass/fail/skip/inspect/quit) ", df_weekly_name, "?\n(layer~", layer, ")(px/fx/sx/ix/q): ") # this should be weekly data
    )

    if (verification_decision %in% simple_decision_list){
      altered_df_list <- alter_verification_column(verification_decision, df_weekly_data)
    }

    if (verification_decision %in% inspect_decision_list){
      # set layer
      layer <- "daily"

      # read in the data relevant to the current site-parameter that we are doing
      # this will be used to make the inspect decision plots. pull this data from the working directory
      # site_param_df <<- df_list[[paste0(site_arg, "-", parameter_arg)]] # if parameter arg is null then pull it all in

      # inspect all weekly ----
      if (verification_decision == "INSPECT ALL"){

        for (i in unique(df_weekly_data$weekday)) {

          daily_plot_data <- site_param_df %>%
            filter((y_w %in% df_weekly_data$y_w) & (weekday == i))

          daily_plot_object <- generate_daily_plot(plot_data_arg = daily_plot_data,
                                                   df_list_arg = df_list_arg,
                                                   site_arg = site_arg,
                                                   parameter_arg = parameter_arg)

          # weekly_plot_object
          weekly_plot_object <- generate_supplemental_weekly_plot(daily_plot_data_arg = daily_plot_data,
                                                                  df_list_arg = df_list_arg,
                                                                  site_arg = site_arg,
                                                                  parameter_arg = parameter_arg) # need to make sure that these functions can read in information from outside of them

          # print ggarrange daily and weekly plots
          print(ggarrange(daily_plot_object, weekly_plot_object, nrow = 2, ncol = 1))

          # Generate the prompt for the user to answer
          day_dt <- as.character(head(daily_plot_data$DT_round, 1)) # reformat so that this reads as a date
          day_choice_prompt <- paste("Would you like to (pass/fail/skip/inspect) ", day_dt, "?\n(layer~", layer, ")(px/fx/sx/ix): ")
          daily_verification_decision <- get_verification_decision(day_choice_prompt)

          while (layer == "daily" & daily_verification_decision == "QUIT") {
            cat("You cannot quit during this step. Please choose a different option.\n")
            daily_verification_decision <- get_verification_decision(day_choice_prompt)
          }

          if (daily_verification_decision %in% simple_decision_list){
            altered_df_list[[i]] <- alter_verification_column(daily_verification_decision, daily_plot_data) # this should get put in a list
          }

          if (daily_verification_decision %in% inspect_decision_list){
            altered_df_list[[i]] <- verify_subdaily(daily_verification_decision_arg = daily_verification_decision,
                                                    daily_plot_data_arg = daily_plot_data,
                                                    daily_plot_object_arg = daily_plot_object,
                                                    weekly_plot_object_arg = weekly_plot_object
            )
          }
        }
      }

      # inspect valid weekly ----
      if (verification_decision == "INSPECT VALID") {

        # get which days to inspect
        days_to_inspect <- unique(df_weekly_data %>%
                                    filter(is.na(flag)) %>%
                                    pull(weekday))

        new_weekday_list <- weekday_list[weekday_list >= min(df_weekly_data$weekday) & weekday_list <= max(df_weekly_data$weekday)]
        non_days_to_inspect <- new_weekday_list[!new_weekday_list %in% days_to_inspect]
        # for each non-inspect day decision:
        if (length(non_days_to_inspect) > 0) {

          # get non-inspect day decision:
          noninspect_prompt_text = paste("Would you like to (pass/fail/skip) non-inspect days?\n(layer~", layer, ")(px/fx/sx): ")
          non_inspect_verification_decision <- get_verification_decision(prompt_text = noninspect_prompt_text)

          # Handle the case when the user selects an "INSPECT" or "QUIT" option for non-inspect days
          while (layer == "daily" & non_inspect_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME", "QUIT")) {
            if(non_inspect_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME")) {
              cat("You cannot inspect non-inspect days. Please choose a different option.\n")
              non_inspect_verification_decision <- get_verification_decision(prompt_text = noninspect_prompt_text)
            }
            if(non_inspect_verification_decision == "QUIT") {
              cat("You cannot quit during this step. Please choose a different option.\n")
              non_inspect_verification_decision <- get_verification_decision(prompt_text = noninspect_prompt_text)
            }
          }

          for (i in non_days_to_inspect) {
            daily_plot_data <- site_param_df %>%
              filter((y_w %in% df_weekly_data$y_w) & (weekday == i))

            altered_df_list[[i]] <- alter_verification_column(non_inspect_verification_decision, daily_plot_data)
          }
        }

        for (i in days_to_inspect) {
          daily_plot_data <- site_param_df %>%
            filter((y_w %in% df_weekly_data$y_w) & (weekday == i))

          daily_plot_object <- generate_daily_plot(plot_data_arg = daily_plot_data,
                                                   df_list_arg = df_list_arg,
                                                   site_arg = site_arg,
                                                   parameter_arg = parameter_arg)

          # weekly_plot_object
          weekly_plot_object <- generate_supplemental_weekly_plot(daily_plot_data_arg = daily_plot_data,
                                                                  df_list_arg = df_list_arg,
                                                                  site_arg = site_arg,
                                                                  parameter_arg = parameter_arg) # need to make sure that these functions can read in information from outside of them

          # print ggarrange daily and weekly plots
          print(ggarrange(daily_plot_object, weekly_plot_object, nrow = 2, ncol = 1))

          # Generate the prompt for the user to answer
          day_dt <- as.character(head(daily_plot_data$DT_round, 1)) # reformat so that this reads as a date
          day_choice_prompt <- paste("Would you like to (pass/fail/skip/inspect) ", day_dt, "?\n(layer~", layer, ")(px/fx/sx/ix): ")
          daily_verification_decision <- get_verification_decision(day_choice_prompt)

          while (daily_verification_decision == "QUIT") {
            cat("You cannot quit during this step. Please choose a different option.\n")
            daily_verification_decision <- get_verification_decision(day_choice_prompt)
          }

          if (daily_verification_decision %in% simple_decision_list){
            altered_df_list[[i]] <- alter_verification_column(daily_verification_decision, daily_plot_data) # this should get put in a list
          }

          if (daily_verification_decision %in% inspect_decision_list){
            altered_df_list[[i]] <- verify_subdaily(daily_verification_decision_arg = daily_verification_decision,
                                                    daily_plot_data_arg = daily_plot_data,
                                                    daily_plot_object_arg = daily_plot_object,
                                                    weekly_plot_object_arg = weekly_plot_object
            )
          }
        }
      }

      # inspect flagged weekly ----
      if (verification_decision == "INSPECT FLAGGED") {

        # get which days to inspect
        days_to_inspect <- unique(df_weekly_data %>%
                                    filter(!is.na(flag)) %>%
                                    pull(weekday))

        new_weekday_list <- weekday_list[weekday_list >= min(df_weekly_data$weekday) & weekday_list <= max(df_weekly_data$weekday)]
        non_days_to_inspect <- new_weekday_list[!new_weekday_list %in% days_to_inspect]

        # for each non-inspect day decision:
        if (length(non_days_to_inspect) > 0) {

          # get non-inspect day decision:
          noninspect_prompt_text = paste("Would you like to (pass/fail/skip) non-inspect days?\n(layer~", layer, ")(px/fx/sx): ")
          non_inspect_verification_decision <- get_verification_decision(prompt_text = noninspect_prompt_text)

          # Handle the case when the user selects an "INSPECT" or "QUIT" option for non-inspect days
          while (layer == "daily" & non_inspect_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME", "QUIT")) {
            if(non_inspect_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME")) {
              cat("You cannot inspect non-inspect days. Please choose a different option.\n")
              non_inspect_verification_decision <- get_verification_decision(prompt_text = noninspect_prompt_text)
            }
            if(non_inspect_verification_decision == "QUIT") {
              cat("You cannot quit during this step. Please choose a different option.\n")
              non_inspect_verification_decision <- get_verification_decision(prompt_text = noninspect_prompt_text)
            }
          }

          for (i in non_days_to_inspect) {
            daily_plot_data <- site_param_df %>%
              filter((y_w %in% df_weekly_data$y_w) & (weekday == i))

            altered_df_list[[i]] <- alter_verification_column(non_inspect_verification_decision, daily_plot_data)
          }
        }

        for (i in days_to_inspect) {
          daily_plot_data <- site_param_df %>%
            filter((y_w %in% df_weekly_data$y_w) & (weekday == i))

          daily_plot_object <- generate_daily_plot(plot_data_arg = daily_plot_data,
                                                   df_list_arg = df_list_arg,
                                                   site_arg = site_arg,
                                                   parameter_arg = parameter_arg)


          # weekly_plot_object
          weekly_plot_object <- generate_supplemental_weekly_plot(daily_plot_data_arg = daily_plot_data,
                                                                  df_list_arg = df_list_arg,
                                                                  site_arg = site_arg,
                                                                  parameter_arg = parameter_arg) # need to make sure that these functions can read in information from outside of them

          # print ggarrange daily and weekly plots
          print(ggarrange(daily_plot_object, weekly_plot_object, nrow = 2, ncol = 1))

          # Generate the prompt for the user to answer
          day_dt <- as.character(head(daily_plot_data$DT_round, 1)) # reformat so that this reads as a date
          day_choice_prompt <- paste("Would you like to (pass/fail/skip/inspect) ", day_dt, "?\n(layer~", layer, ")(px/fx/sx/ix): ")
          daily_verification_decision <- get_verification_decision(day_choice_prompt)

          while (layer == "daily" & daily_verification_decision == "QUIT") {
            cat("You cannot quit during this step. Please choose a different option.\n")
            daily_verification_decision <- get_verification_decision(day_choice_prompt)
          }

          if (daily_verification_decision %in% simple_decision_list){
            altered_df_list[[i]] <- alter_verification_column(daily_verification_decision, daily_plot_data) # this should get put in a list
          }

          if (daily_verification_decision %in% inspect_decision_list){
            altered_df_list[[i]] <- verify_subdaily(daily_verification_decision_arg = daily_verification_decision,
                                                    daily_plot_data_arg = daily_plot_data,
                                                    daily_plot_object_arg = daily_plot_object,
                                                    weekly_plot_object_arg = weekly_plot_object
            )
          }
        }
      }

      # inspect some weekly ----
      if (verification_decision == "INSPECT SOME") {

        # get which days to inspect
        days_to_inspect <- get_weekly_inspection_decision(weekly_plot_data = df_weekly_data) # rename to daily ***

        new_weekday_list <- weekday_list[weekday_list >= min(df_weekly_data$weekday) & weekday_list <= max(df_weekly_data$weekday)]
        non_days_to_inspect <- new_weekday_list[!new_weekday_list %in% days_to_inspect]

        # for each non-inspect day decision:
        if (length(non_days_to_inspect)>0) {

          # get non-inspect day decision:
          noninspect_prompt_text = paste("Would you like to (pass/fail/skip) non-inspect days?\n(layer~", layer, ")(px/fx/sx): ")
          non_inspect_verification_decision <- get_verification_decision(prompt_text = noninspect_prompt_text)

          # Handle the case when the user selects an "INSPECT" or "QUIT" option for non-inspect days
          while (layer == "daily" & non_inspect_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME", "QUIT")) {
            if(non_inspect_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME")) {
              cat("You cannot inspect non-inspect days. Please choose a different option.\n")
              non_inspect_verification_decision <- get_verification_decision(prompt_text = noninspect_prompt_text)
            }
            if(non_inspect_verification_decision == "QUIT") {
              cat("You cannot quit during this step. Please choose a different option.\n")
              non_inspect_verification_decision <- get_verification_decision(prompt_text = noninspect_prompt_text)
            }
          }

          for (i in non_days_to_inspect) {
            daily_plot_data <- site_param_df %>%
              filter((y_w %in% df_weekly_data$y_w) & (weekday == i))

            altered_df_list[[i]] <- alter_verification_column(non_inspect_verification_decision, daily_plot_data)
          }
        }

        for (i in days_to_inspect) {

          daily_plot_data <- site_param_df %>%
            filter((y_w %in% df_weekly_data$y_w) & (weekday == i))

          daily_plot_object <- generate_daily_plot(plot_data_arg = daily_plot_data,
                                                   df_list_arg = df_list_arg,
                                                   site_arg = site_arg,
                                                   parameter_arg = parameter_arg)

          # weekly_plot_object

          weekly_plot_object <- generate_supplemental_weekly_plot(daily_plot_data_arg = daily_plot_data,
                                                                  df_list_arg = df_list_arg,
                                                                  site_arg = site_arg,
                                                                  parameter_arg = parameter_arg)

          # print ggarrange daily and weekly plots
          print(ggarrange(daily_plot_object, weekly_plot_object, nrow = 2, ncol = 1))

          # Generate the prompt for the user to answer
          day_dt <- as.character(head(daily_plot_data$DT_round, 1)) # reformat so that this reads as a date
          day_choice_prompt <- paste("Would you like to (pass/fail/skip/inspect) ", day_dt, "?\n(layer~", layer, ")(px/fx/sx/ix): ")
          daily_verification_decision <- get_verification_decision(day_choice_prompt)

          while (layer == "daily" & daily_verification_decision == "QUIT") {
            cat("You cannot quit during this step. Please choose a different option.\n")
            daily_verification_decision <- get_verification_decision(day_choice_prompt)
          }

          if (daily_verification_decision %in% simple_decision_list){
            # browser()
            altered_df_list[[i]] <- alter_verification_column(daily_verification_decision, daily_plot_data) # this should get put in a list
          }

          if (daily_verification_decision %in% inspect_decision_list){
            altered_df_list[[i]] <- verify_subdaily(daily_verification_decision_arg = daily_verification_decision,
                                                    daily_plot_data_arg = daily_plot_data,
                                                    daily_plot_object_arg = daily_plot_object,
                                                    weekly_plot_object_arg = weekly_plot_object
            )
          }
        }
      }

    }


    if (is_tibble(altered_df_list)) {
      altered_df_list <- altered_df_list %>%
        mutate(mean_verified = if_else(grepl(non_removable_flags, flag, ignore.case = TRUE), NA, mean_verified),
               is_verified = TRUE,
               verification_status = ifelse(grepl(non_removable_flags, flag, ignore.case = TRUE), "PASS", verification_status))

    } else if (is.list(altered_df_list)) {
      if (!length(altered_df_list) %in% c(51,49,46)){
        # browser()
        altered_df_list <- map(altered_df_list, function(df) {
          if(!is.null(df)) {
            # browser()
            df %>%
              mutate(mean_verified = if_else(grepl(non_removable_flags, flag, ignore.case = TRUE), NA, mean_verified),
                     is_verified = TRUE,
                     verification_status = ifelse(grepl(non_removable_flags, flag, ignore.case = TRUE), "PASS", verification_status))
          }
        })
      } else {
        # browser()
        altered_df_list <- altered_df_list %>%
          mutate(mean_verified = if_else(grepl(non_removable_flags, flag, ignore.case = TRUE), NA, mean_verified),
                 is_verified = TRUE,
                 verification_status = ifelse(grepl(non_removable_flags, flag, ignore.case = TRUE), "PASS", verification_status))
      }

    }

    return(altered_df_list) # this should not be an altered df, but rather a list of dfs that have been altered.
    # this means that the output of this method will be a list that contains a list of lists

    # If decision in quit (this can only happen before they start verifying a week of data)
    # save the data into the intermediary directory
  } else if (QUIT) {
    df_weekly_name <- ggplot_build(weekly_plot_object)$plot$labels$title
    cat("quitting out of", df_weekly_name, "\n")
    return(NULL)
  }

}


verify_subdaily <- function(daily_verification_decision_arg,
                            daily_plot_data_arg,
                            daily_plot_object_arg,
                            weekly_plot_object_arg){

  # set layer
  layer <- "subdaily"

  # inspect all
  if (daily_verification_decision_arg == "INSPECT ALL") {

    for (j in 1:nrow(daily_plot_data_arg)) {

      # create a data frame with a single row containing the desired DT_round value
      # and add it to the week plot and day plot
      flag_demarcation_data <- tibble(DT_round = daily_plot_data_arg$DT_round[j])

      # ignore flags that you cannot alter
      daily_plot_object_new <- daily_plot_object_arg +
        geom_vline(data = flag_demarcation_data,
                   aes(xintercept = DT_round,
                       color = "data point of interest",
                       # size = 2, # fix this later
                       alpha = 0.05))

      weekly_plot_object_new <- weekly_plot_object_arg +
        geom_vline(data = flag_demarcation_data,
                   aes(xintercept = DT_round,
                       color = "data point of interest",
                       # size = 2, # fix this later
                       alpha = 0.05))

      # print ggarrange daily and weekly plots
      print(ggarrange(daily_plot_object_new, weekly_plot_object_new, nrow = 2, ncol = 1))

      # generate the prompt for the user to answer
      # the user should not be able to inspect again
      day_dt <- daily_plot_data_arg$DT_join[j] # reformat so that this reads as a date
      dt_choice_prompt <- paste("Would you like to (pass/fail/skip) ", day_dt, "?\n(layer~", layer, ")(px/fx/sx)")
      dt_verification_decision <- get_verification_decision(dt_choice_prompt)

      # Handle the case when the user selects an "INSPECT" or "QUIT" option inappropriately
      while (layer == "subdaily" & dt_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME", "QUIT")) {
        if (dt_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME")){
          cat("You cannot inspect at this level. Please choose a different option.\n")
          dt_verification_decision <- get_verification_decision(dt_choice_prompt)
        }
        if (dt_verification_decision == "QUIT"){
          cat("You cannot quit during this step. Please choose a different option.\n")
          dt_verification_decision <- get_verification_decision(dt_choice_prompt)
        }
      }

      # update the specific row in df_daily_data based on the index
      j_index <- which(daily_plot_data_arg$DT_join == day_dt)

      altered_row <- alter_verification_column(dt_verification_decision, daily_plot_data_arg[j_index, ])

      daily_plot_data_arg <- daily_plot_data_arg %>%
        mutate(
          mean_verified = if_else(row_number() == j_index, altered_row$mean_verified, mean_verified),
          is_verified = if_else(row_number() == j_index, altered_row$is_verified, is_verified),
          verification_status = if_else(row_number() == j_index, altered_row$verification_status, verification_status)
        )
    }
    # inform the user that they are done verifying which ever df
    cat("Finished verifying: ", as.character(min(daily_plot_data_arg$DT_round)), " data.\n")

  }

  # inspect valid
  if (daily_verification_decision_arg == "INSPECT VALID") {

    # get the dt of the valid data to inspect to get indexes of that data
    valid_daily_plot_data_dt_list <- daily_plot_data_arg %>%
      filter(is.na(flag)) %>%
      pull(DT_join)

    # get the dt of the flagged data to get the non-inspect decision
    new_dt_list <- pull(daily_plot_data_arg, DT_join)
    non_dts_to_inspect <- new_dt_list[!new_dt_list %in% valid_daily_plot_data_dt_list]

    # decide what to do with the flagged data ----
    # for each non-inspect dt decision:
    if (length(non_dts_to_inspect) > 0) {

      # get non-inspect day decision:
      noninspect_prompt_text = paste("Would you like to (pass/fail/skip) non-inspect dts?\n(layer~", layer, ")(px/fx/sx): ")
      non_inspect_verification_decision <- get_verification_decision(noninspect_prompt_text)

      # Handle the case when the user selects an "INSPECT" or "QUIT" option inappropriately
      while (non_inspect_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME", "QUIT")) {
        if (non_inspect_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME")){
          cat("You cannot inspect non-inspect dts. Please choose a different option.\n")
          non_inspect_verification_decision <- get_verification_decision(noninspect_prompt_text)
        }
        if (non_inspect_verification_decision == "QUIT"){
          cat("You cannot quit during this step. Please choose a different option.\n")
          non_inspect_verification_decision <- get_verification_decision(noninspect_prompt_text)
        }
      }

      for (j in non_dts_to_inspect) {
        # j <- format(as.POSIXct(j, tz = "MST"), "%Y-%m-%d %H:%M:%S")

        j_index <- which((daily_plot_data_arg$DT_join) == j)

        altered_row <- alter_verification_column(non_inspect_verification_decision, daily_plot_data_arg[j_index, ])

        daily_plot_data_arg <- daily_plot_data_arg %>%
          mutate(
            mean_verified = if_else(row_number() == j_index, altered_row$mean_verified, mean_verified),
            is_verified = if_else(row_number() == j_index, altered_row$is_verified, is_verified),
            verification_status = if_else(row_number() == j_index, altered_row$verification_status, verification_status)
          )

        # some sort of print statement to show what happened to the non-target dts
      }
    }

    # update those flagged data in the df

    for (j in valid_daily_plot_data_dt_list) {

      # j <- format(as.POSIXct(j, tz = "MST"), "%Y-%m-%d %H:%M:%S")

      j_index <- which((daily_plot_data_arg$DT_join) == j)

      # create a data frame with a single row containing the desired DT_round value
      # and add it to the week plot and day plot
      flag_demarcation_data <- tibble(DT_round = daily_plot_data_arg$DT_round[j_index])

      daily_plot_object_new <- daily_plot_object_arg +
        geom_vline(data = flag_demarcation_data,
                   aes(xintercept = DT_round,
                       color = "data point of interest",
                       # size = 2, # fix this later
                       alpha = 0.05))

      weekly_plot_object_new <- weekly_plot_object_arg +
        geom_vline(data = flag_demarcation_data,
                   aes(xintercept = DT_round,
                       color = "data point of interest",
                       # size = 2, # fix this later
                       alpha = 0.05))

      # print ggarrange daily and weekly plots
      print(ggarrange(daily_plot_object_new, weekly_plot_object_new, nrow = 2, ncol = 1))

      # generate the prompt for the user to answer
      # the user should not be able to inspect again
      day_dt <- as.character(daily_plot_data_arg$DT_round[j_index]) # reformat so that this reads as a date
      dt_choice_prompt <- paste("Would you like to (pass/fail/skip) ", day_dt, "?\n(layer~", layer, ")(px/fx/sx)")
      dt_verification_decision <- get_verification_decision(dt_choice_prompt)

      # Handle the case when the user selects an "INSPECT" or "QUIT" option inappropriately
      while (layer == "subdaily" & dt_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME", "QUIT")) {
        if (dt_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME")){
          cat("You cannot inspect at this level. Please choose a different option.\n")
          dt_verification_decision <- get_verification_decision(dt_choice_prompt)
        }
        if (dt_verification_decision == "QUIT"){
          cat("You cannot quit during this step. Please choose a different option.\n")
          dt_verification_decision <- get_verification_decision(dt_choice_prompt)
        }
      }

      # update the specific row in df_daily_data based on the index
      # row_index <- which(as.character(daily_plot_data_arg$DT_round) == day_dt)
      altered_row <- alter_verification_column(dt_verification_decision, daily_plot_data_arg[j_index, ])

      daily_plot_data_arg <- daily_plot_data_arg %>%
        mutate(
          mean_verified = if_else(row_number() == j_index, altered_row$mean_verified, mean_verified),
          is_verified = if_else(row_number() == j_index, altered_row$is_verified, is_verified),
          verification_status = if_else(row_number() == j_index, altered_row$verification_status, verification_status)
        )
    }
    # inform the user that they are done verifying which ever df
    cat("Finished verifying: ", as.character(min(daily_plot_data_arg$DT_round)), " data.\n")
    # return the altered df
    # altered_df_list_arg[[i]] <- daily_plot_data_arg # make sure that this is the updated data frame

  }

  # inspect flagged
  if (daily_verification_decision_arg == "INSPECT FLAGGED") {

    # get the dt of the flagged data to inspect to get indexes of that data
    flagged_daily_plot_data_dt_list <- daily_plot_data_arg %>%
      filter(!is.na(flag)) %>%
      pull(DT_join)

    # get the dt of the flagged data to get the non-inspect decision
    new_dt_list <- pull(daily_plot_data_arg, DT_join)
    non_dts_to_inspect <- new_dt_list[!new_dt_list %in% flagged_daily_plot_data_dt_list]

    # decide what to do with the flagged data ----
    # for each non-inspect dt decision:
    if (length(non_dts_to_inspect) > 0) {

      # get non-inspect day decision:
      noninspect_prompt_text = paste("Would you like to (pass/fail/skip) non-inspect dts?\n(layer~", layer, ")(px/fx/sx): ")
      non_inspect_verification_decision <- get_verification_decision(noninspect_prompt_text)

      # Handle the case when the user selects an "INSPECT" or "QUIT" option inappropriately
      while (non_inspect_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME", "QUIT")) {
        if (non_inspect_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME")){
          cat("You cannot inspect non-inspect dts. Please choose a different option.\n")
          non_inspect_verification_decision <- get_verification_decision(noninspect_prompt_text)
        }
        if (non_inspect_verification_decision == "QUIT"){
          cat("You cannot quit during this step. Please choose a different option.\n")
          non_inspect_verification_decision <- get_verification_decision(noninspect_prompt_text)
        }
      }

      for (j in non_dts_to_inspect) {
        # j <- format(as.POSIXct(j, tz = "MST"), "%Y-%m-%d %H:%M:%S") # time zone

        j_index <- which((daily_plot_data_arg$DT_join) == j)

        altered_row <- alter_verification_column(non_inspect_verification_decision, daily_plot_data_arg[j_index, ])

        daily_plot_data_arg <- daily_plot_data_arg %>%
          mutate(
            mean_verified = if_else(row_number() == j_index, altered_row$mean_verified, mean_verified),
            is_verified = if_else(row_number() == j_index, altered_row$is_verified, is_verified),
            verification_status = if_else(row_number() == j_index, altered_row$verification_status, verification_status)
            )
        # some sort of print statement to show what happened to the non-target dts
      }
    }

    # update those flagged data in the df

    for (j in flagged_daily_plot_data_dt_list) {

      # j <- format(as.POSIXct(j, tz = "MST"), "%Y-%m-%d %H:%M:%S")

      j_index <- which((daily_plot_data_arg$DT_join) == j)

      # create a data frame with a single row containing the desired DT_round value
      # and add it to the week plot and day plot
      flag_demarcation_data <- tibble(DT_round = daily_plot_data_arg$DT_round[j_index])

      daily_plot_object_new <- daily_plot_object_arg +
        geom_vline(data = flag_demarcation_data,
                   aes(xintercept = DT_round,
                       color = "data point of interest",
                       # size = 2, # fix this later
                       alpha = 0.05))

      weekly_plot_object_new <- weekly_plot_object_arg +
        geom_vline(data = flag_demarcation_data,
                   aes(xintercept = DT_round,
                       color = "data point of interest",
                       # size = 2, # fix this later
                       alpha = 0.05))


      # print ggarrange daily and weekly plots
      print(ggarrange(daily_plot_object_new, weekly_plot_object_new, nrow = 2, ncol = 1))

      # generate the prompt for the user to answer
      # the user should not be able to inspect again
      day_dt <- as.character(daily_plot_data_arg$DT_round[j_index]) # reformat so that this reads as a date
      dt_choice_prompt <- paste("Would you like to (pass/fail/skip) ", day_dt, "?\n(layer~", layer, ")(px/fx/sx)")
      dt_verification_decision <- get_verification_decision(dt_choice_prompt)

      # Handle the case when the user selects an "INSPECT" or "QUIT" option inappropriately
      while (layer == "subdaily" & dt_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME", "QUIT")) {
        if (dt_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME")){
          cat("You cannot inspect at this level. Please choose a different option.\n")
          dt_verification_decision <- get_verification_decision(dt_choice_prompt)
        }
        if (dt_verification_decision == "QUIT"){
          cat("You cannot quit during this step. Please choose a different option.\n")
          dt_verification_decision <- get_verification_decision(dt_choice_prompt)
        }
      }

      # update the specific row in df_daily_data based on the index
      altered_row <- alter_verification_column(dt_verification_decision, daily_plot_data_arg[j_index, ])

      daily_plot_data_arg <- daily_plot_data_arg %>%
        mutate(
          mean_verified = if_else(row_number() == j_index, altered_row$mean_verified, mean_verified),
          is_verified = if_else(row_number() == j_index, altered_row$is_verified, is_verified),
          verification_status = if_else(row_number() == j_index, altered_row$verification_status, verification_status)
        )

    }
    # inform the user that they are done verifying which ever df
    cat("Finished verifying: ", as.character(min(daily_plot_data_arg$DT_round)), " data.\n")

  }

  # inspect some
  if (daily_verification_decision_arg == "INSPECT SOME") {

    # get which dts to inspect
    dt_intervals_to_inspect <- get_dt_inspection_decisions(daily_plot_data = daily_plot_data_arg)

    # decided what to do with the non-inspect data
    new_dt_list <- pull(daily_plot_data_arg, DT_round)
    non_dts_to_inspect <- new_dt_list[!new_dt_list %within% dt_intervals_to_inspect]

    # decide what to do with the flagged data ----
    # for each non-inspect dt decision:
    if (length(non_dts_to_inspect) > 0) {

      # get non-inspect day decision:
      noninspect_prompt_text = paste("Would you like to (pass/fail/skip) non-inspect dts?\n(layer~", layer, ")(px/fx/sx): ")
      non_inspect_verification_decision <- get_verification_decision(noninspect_prompt_text)

      # Handle the case when the user selects an "INSPECT" or "QUIT" option inappropriately
      while (non_inspect_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME", "QUIT")) {
        if (non_inspect_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME")){
          cat("You cannot inspect non-inspect dts. Please choose a different option.\n")
          non_inspect_verification_decision <- get_verification_decision(noninspect_prompt_text)
        }
        if (non_inspect_verification_decision == "QUIT"){
          cat("You cannot quit during this step. Please choose a different option.\n")
          non_inspect_verification_decision <- get_verification_decision(noninspect_prompt_text)
        }
      }

      for (j in non_dts_to_inspect) {
        j_time <- format(as.POSIXct(j, tz = "MST"), "%H:%M:%S")
        j <- if_else(
          j_time == "00:00:00",
          format(as.POSIXct(j, tz = "MST"), "%Y-%m-%d"),
          format(as.POSIXct(j, tz = "MST"), "%Y-%m-%d %H:%M:%S")
        )

        j_index <- which((daily_plot_data_arg$DT_join) == j)

        # daily_plot_data_arg[j_index, ] <- alter_verification_column(non_inspect_verification_decision, daily_plot_data_arg[j_index, ])

        altered_row <- alter_verification_column(non_inspect_verification_decision, daily_plot_data_arg[j_index, ])

        daily_plot_data_arg <- daily_plot_data_arg %>%
          mutate(
            mean_verified = if_else(row_number() == j_index, altered_row$mean_verified, mean_verified),
            is_verified = if_else(row_number() == j_index, altered_row$is_verified, is_verified),
            verification_status = if_else(row_number() == j_index, altered_row$verification_status, verification_status)
          )
        # some sort of print statement to show what happened to the non-target dts
      }
    }

    # update those flagged data in the df

    for (j in dt_intervals_to_inspect) { # make sure that j can not be 0

      interval_df <- daily_plot_data_arg %>%
        filter(DT_round %within% j)

      for (k in 1:nrow(interval_df)) {

        k_index <- which((daily_plot_data_arg$DT_round) == interval_df$DT_round[k])

        # create a data frame with a single row containing the desired DT_round value
        # and add it to the week plot and day plot
        flag_demarcation_data <- tibble(DT_round = daily_plot_data_arg$DT_round[k_index])

        daily_plot_object_new <- daily_plot_object_arg +
          geom_vline(data = flag_demarcation_data,
                     aes(xintercept = DT_round,
                         color = "data point of interest",
                         # size = 2, # fix this later
                         alpha = 0.05))

        weekly_plot_object_new <- weekly_plot_object_arg +
          geom_vline(data = flag_demarcation_data,
                     aes(xintercept = DT_round,
                         color = "data point of interest",
                         # size = 2, # fix this later
                         alpha = 0.05))


        # print ggarrange daily and weekly plots
        print(ggarrange(daily_plot_object_new, weekly_plot_object_new, nrow = 2, ncol = 1))

        # generate the prompt for the user to answer
        # the user should not be able to inspect again
        day_dt <- as.character(daily_plot_data_arg$DT_round[k_index]) # reformat so that this reads as a date
        dt_choice_prompt <- paste("Would you like to (pass/fail/skip) ", day_dt, "?\n(layer~", layer, ")(px/fx/sx)")
        dt_verification_decision <- get_verification_decision(dt_choice_prompt)

        # Handle the case when the user selects an "INSPECT" or "QUIT" option inappropriately
        while (layer == "subdaily" & dt_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME", "QUIT")) {
          if (dt_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME")){
            cat("You cannot inspect at this level. Please choose a different option.\n")
            dt_verification_decision <- get_verification_decision(dt_choice_prompt)
          }
          if (dt_verification_decision == "QUIT"){
            cat("You cannot quit during this step. Please choose a different option.\n")
            dt_verification_decision <- get_verification_decision(dt_choice_prompt)
          }
        }

        # update the specific row in df_daily_data based on the index
        # row_index <- which(as.character(daily_plot_data_arg$DT_round) == day_dt)
        # daily_plot_data_arg[k_index, ] <- alter_verification_column(dt_verification_decision, daily_plot_data_arg[k_index, ])

        altered_row <- alter_verification_column(dt_verification_decision, daily_plot_data_arg[k_index, ])

        daily_plot_data_arg <- daily_plot_data_arg %>%
          mutate(
            mean_verified = if_else(row_number() == k_index, altered_row$mean_verified, mean_verified),
            is_verified = if_else(row_number() == k_index, altered_row$is_verified, is_verified),
            verification_status = if_else(row_number() == k_index, altered_row$verification_status, verification_status)
          )
      }
    }
    # inform the user that they are done verifying which ever df
    cat("Finished verifying: ", as.character(min(daily_plot_data_arg$DT_round)), " data.\n")
  }

  return(daily_plot_data_arg)

}
