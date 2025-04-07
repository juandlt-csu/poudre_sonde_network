#' @title Load and tidy old field notes
#'
#' @description A function that uploads and cleans the field notes excel file. This function adds datetime
#' columns to the field notes dataframe and filters out field notes where the sensor
#' was not handled.
#' @param filepath A file path to the raw field notes.
#' @param summarize_interval At what time interval the user would like the data set to be aggregated and rounded to. Default is 15 minutes.
#' @return A dataframe with the field notes.
#'
#' @examples
#' clean_old_field_notes(filepath = "data/sensor_field_notes.xlsx")

load_old_field_notes <- function(filepath,  summarize_interval = "15 minutes"){

  raw_field_notes <- readxl::read_excel(filepath) # TODO: convert the old field notes to csv

  field_notes <- raw_field_notes %>%
    dplyr::mutate(start_DT = lubridate::ymd_hm(paste(date, start_time_mst), tz = "MST")) %>%
    dplyr::mutate(
      DT_round = lubridate::floor_date(start_DT, unit = summarize_interval),
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
                                             is.na(sensor_pulled) & is.na(sensor_deployed) ~ 0),
           end_dt  = as.POSIXct(NA, tz = "MST")) %>%
    # remove field dates where sensor was not handled:
    dplyr::filter(grepl("Sensor Cleaning or Check|Sensor Calibration", visit_type, ignore.case = TRUE))

  return(field_notes)

}
