#' @title Summarize site parameter data from the API and field notes data frames.
#'
#' @description A function that summarizes and joins site parameter data from the API with the field notes data frames.
#'
#' @param site_arg A site name.
#' @param parameter_arg A parameter name.
#' @param api_data A dataframe with the munged API data.
#' @param notes The munged field notes
#'
#' @return A dataframe with summary statistics for a given site parameter data frame, plus field notes data
#'
#'  @examples
# summarize_site_param(site_arg = "archery", parameter_arg = "Actual Conductivity", api_data = incoming_data_collated_csvs)
# summarize_site_param(site_arg = "boxelder", parameter_arg = "Temperature", api_data = incoming_data_collated_csvs)

tidy_and_add_field_notes <- function(site_arg, parameter_arg, api_data, notes) {

  # filter deployment records for the full join
  site_field_notes <- notes %>%
    dplyr::filter(grepl(paste(unlist(stringr::str_split(site_arg, " ")), collapse = "|"), site, ignore.case = TRUE))

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
      # pad the dataset so that all user-selected interval time stamps are present
      padr::pad(by = "DT_round", interval = "15 minutes") %>%
      # add a DT_join column to join field notes to
      dplyr::mutate(DT_join = as.character(DT_round),
                    site = site_arg,
                    parameter = parameter_arg,
                    flag = NA) %>%
      # join our tidied data frame with our field notes data:
      dplyr::full_join(., dplyr::select(site_field_notes, sonde_employed,
                                        last_site_visit, DT_join, visit_comments,
                                        sensor_malfunction, cals_performed),
                       by = c('DT_join')) %>%
      arrange((DT_join)) %>%
    # make sure DT_join is still correct:
    dplyr::mutate(DT_round = lubridate::as_datetime(DT_join, tz = "MST")) %>%
      # Use fill() to determine when sonde was in the field, and when the last site visit was
      # Necessary step for FULL dataset only (this step occurs in combine_hist_inc_data.R for auto QAQC)
      dplyr::mutate(sonde_employed = ifelse(is.na(sonde_employed), 0, sonde_employed)) %>%
      tidyr::fill(c(sonde_employed, last_site_visit, sensor_malfunction)) %>%
      # for instances at the top of df's where log was running ahead of deployment:
      dplyr::mutate(sonde_employed = ifelse(is.na(last_site_visit), 1, sonde_employed)) %>%
      dplyr::distinct(.keep_all = TRUE) %>%
      dplyr::filter(!is.na(site))

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
