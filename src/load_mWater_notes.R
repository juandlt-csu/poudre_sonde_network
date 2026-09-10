#' @title Load and tidy mWater field notes
#'
#' @description A function that uploads and cleans the field notes submitted to mWater.
#'
#' @param creds A .yml file with necessary credentials for accessing the field notes. Contact Sam Struthers if you need access.
#' @param summarize_interval At what time interval the user would like the data set to be aggregated and rounded to. Default is 15 minutes.
#' @return A dataframe with the field notes.

load_mWater_notes <- function(creds = yaml::read_yaml(here::here("creds", "mWaterCreds.yml")), summarize_interval = "15 minutes"){

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
      #did sonde/vulinks change
      sensor_pulled = as.character(sn_removed),
      vulink_pulled = as.character(sn_removed_vulink),
      sensor_deployed = as.character(sn_deployed),
      vulink_deployed = as.character(sn_deployed_vulink),
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
      DT_round = lubridate::floor_date(start_DT, unit = summarize_interval)) %>%
    # arrange by most recent visit
    dplyr::arrange(DT_round)%>%
    # Remove other columns
    dplyr::select(-c(photos_downloaded_other,visit_type_other, site_other, other_which_sensor_malfunction ))

  return(all_notes_cleaned)

}
