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

