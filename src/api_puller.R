#' @title Download API data from HydroVu
#'
#' @description A function designed to download raw API data from HydroVu and store in a user-selected file location.
#'
#' @param site List of sites to download data for.
#' @param start_dt Start DT for downloading data.
#' @param end_dt End DT for downloading data. Default is the current time.
#' @param api_token API token for accessing HydroVu. Pulled in via `hv_auth()`.
#' @param dump_dir Directory where raw API data will be saved to.

# To Do:
# - hard code:
#   - sites
#   - network
#   - dump_dir
# - Append strategy: Check that historical data exists, and then

api_puller <- function(site, network = "CSU", start_dt, end_dt = Sys.time(), api_token, dump_dir) {

  locs <- hv_locations_all(hv_token) %>%
    dplyr::filter(!grepl("virridy", name, ignore.case = TRUE))

  # make a list of site data frames
  incoming_data_list <- list()

  for(i in 1:length(site)){

    site_loc <- locs %>%
      dplyr::mutate(name = tolower(name)) %>%
      dplyr::filter(grepl(site[i], name, ignore.case = TRUE))

    site_loc_list <- site_loc$id

    # Fetches data for all site locations. Expect 404 errors for inactive sites.
    # 'alldata' will contain data from active sites for specified timeframe.
    # Note: Processing takes ~2 mins per site per month. Uses UTC timezone.
    # Use most recent download as start date to avoid duplicates.

    utc_start_date <- format(as.POSIXct(start_dt, tz = "UTC") + lubridate::hours(0), format = "%Y-%m-%d %H:%M:%S")
    utc_end_date <- format(as.POSIXct(end_dt, tz = "UTC") + lubridate::hours(0), format = "%Y-%m-%d %H:%M:%S")
    timezone <- "UTC"

    # Map over the location ids
    alldata <- site_loc_list %>% purrr::map(~hv_data_id(.,
                                                        start_time = utc_start_date,
                                                        end_time = utc_end_date,
                                                        token = api_token,
                                                        tz = timezone))

    # grab only locations with data (stored as a data frame) / drop 404 errors
    filtered <- purrr::keep(alldata, is.data.frame) # add is_tibble to this as well

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

        readr::write_csv(one_df,
                         paste0(dump_dir, "/", site[i], "_", stringr::str_replace(stringr::str_replace(substr(end_dt, 1, 16), "[:\\s]", "_"), ":", ""), ".csv"))

        incoming_data_list[[site[i]]] <- one_df
    }
  }
  return(incoming_data_list)
}

#' @title Summarize site parameter data from the API and field notes data frames.
#'
#' @description A function that summarizes and cleans site parameter data from HydroVu API.
#'
#' @param site_arg A site name.
#' @param parameter_arg A parameter name.
#' @param api_data A dataframe with the munged API data.
#' @param summarize_interval At what time interval the user would like the data set to be aggregated and rounded to. Default is 15 minutes
#' @param notes The munged field notes
#'
#' @return A dataframe with summary statistics for a given site parameter data frame
#' that will get saved into the `all_data_tidied_list` object in the linear
#' workflow
#'
#' @examples
# summarize_site_param(site_arg = "archery", parameter_arg = "Actual Conductivity", api_data = incoming_data_collated_csvs)
# summarize_site_param(site_arg = "boxelder", parameter_arg = "Temperature", api_data = incoming_data_collated_csvs)

tidy_data <- function(site_list_arg, parameter_list_arg, incoming_api_data_list_arg){

  # Munge the incoming API data list
  munged_incoming_data <- incoming_api_data_list %>%
    # remove overlapping API-pull data
    dplyr::distinct() %>%
    # remove VuLink data and Virridy sondes
    dplyr::filter(!grepl("vulink|virridy", name, ignore.case = TRUE)) %>%
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

  # filtering the data and generating results
  summary <-
    map2(.x = site_list_arg,
         .y = parameter_list_arg,
         ~{
           tryCatch({
             munged_incoming_data %>%
               # subset to single site-parameter combo:
               dplyr::filter(site == .x & parameter == .y) %>%
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
               padr::pad(by = "DT_round", interval = summarize_interval) %>%
               # add a DT_join column to join field notes to (make DT_round character string, so no
               # funky DT issues occur during the join):
               dplyr::mutate(DT_join = as.character(DT_round),
                             site = site_arg,
                             parameter = parameter_arg,
                             flag = NA) %>% # add "flag" column for future processing
               # join our tidied data frame with our field notes data:
               # dplyr::full_join(., dplyr::select(site_field_notes, sonde_employed,
               #                                   last_site_visit, DT_join, visit_comments,
               #                                   sensor_malfunction, cals_performed),
               #                  by = c('DT_join')) %>%
               # arrange((DT_join)) %>%
               # make sure DT_join is still correct:
               # dplyr::mutate(DT_round = lubridate::as_datetime(DT_join, tz = "MST")) %>%
               # Use fill() to determine when sonde was in the field, and when the last site visit was
               # Necessary step for FULL dataset only (this step occurs in combine_hist_inc_data.R for auto QAQC)
               # dplyr::mutate(sonde_employed = ifelse(is.na(sonde_employed), 0, sonde_employed)) %>%
               # tidyr::fill(c(sonde_employed, last_site_visit, sensor_malfunction)) %>%
               # # for instances at the top of df's where log was running ahead of deployment:
               # dplyr::mutate(sonde_employed = ifelse(is.na(last_site_visit), 1, sonde_employed)) %>%
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
           }) %>%
    # set the names for the dfs in the list
    purrr::set_names(paste0(site_param_combos$sites, "-", site_param_combos$params)) %>%
    # remove NULL from the list (indicating a site-param combo that doesn't actually exist)
    purrr::keep( ~ !is.null(.))


  return(summary)

}

#' @title Combine 24hrs of past auto-qaqc data to the new incoming data
#'
#' @description A function that combines the incoming API data with 24hrs of past API data
#' that has already been auto-qaqc'd so that we can accurately attribute rolling window
#' flags.
#'
#' @param tidy_data_list A list of summarized and cleaned data frames that are coming
#' from the HydroVu API.
#' @param past_data_list A list of data frames that have already been flagged.
#'
#' @return A list of data frames that have a 3-hour chunk of incoming data appended to
#' 24 hours of previously flagged data.
#'
#' @examples
#' combine_hist_inc_data(incoming_data_list = incoming_data_list, past_data_list = past_data_list)

combine_past_and_incoming_data <- function(incoming_data_list, past_data_list) {

  # Get the matching index names
  matching_indexes <- intersect(names(incoming_data_list), names(past_data_list))

  # Get the last 24 hours of the past flagged data based on earliest
  # DT listed in the api pull for that site/data combo

  last_24_hours <- map(past_data_list,
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
      # pad the data in case the latest past DT and the earliest incoming DT are not 15 min apart
      pad(by = "DT_round", interval = "15 mins") %>%
      # make sure new data has info about last site visit, etc. filled out
      fill(c(sonde_employed, last_site_visit, sensor_malfunction))
  }) %>%
    set_names(matching_indexes)

  return(combined_hist_inc_data)
}
