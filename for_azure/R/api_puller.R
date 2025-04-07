#' @title Download API data from HydroVu
#'
#' @description
#' A function designed to download raw API data from HydroVu and store in a user-selected file location.
#'
#' @param site List of sites to download data for.
#' @param start_dt Start DT for downloading data.
#' @param end_dt End DT for downloading data. Default is the current time.
#' @param api_token API token for accessing HydroVu. Pulled in via `hv_auth()`.
#' @param dump_dir Directory where raw API data will be saved to.
#'
#' @return A CSV file for each site containing the raw API data.
#'
#' @examples
#' dontrun
#' api_puller(site = c("archery", "boxelder"), start_dt = "2021-01-01 00:00:00", api_token = hv_auth(), dump_dir = "data/raw")
#' }

api_puller <- function(site, start_dt, end_dt = Sys.time(), api_token, dump_dir) {

  locs <- hv_locations_all(hv_token) %>% # this should be 'api_token' arg
    dplyr::filter(!grepl("virridy", name, ignore.case = TRUE))

  # make a list of site names
  options(scipen = 999)

  for(i in 1:length(site)){

    site_loc <- locs %>%
      dplyr::mutate(name = tolower(name)) %>%
      dplyr::filter(grepl(site[i], name, ignore.case = TRUE))

    site_loc_list <- site_loc$id

    # Get data for each site location. Note this maps over the entire list of locations,
    # many of which are unlikely to be active during the time you specify. It is common
    # to see many '404 Not Found' errors. These are the list of locations
    # that are not active. The data frame 'alldata' should contain data from all applicable
    # sites during the time frame indicated. Note that this will take some time (one month of
    # data for 5 sites will take roughly ~10 minutes).

    utc_start_date <- format(as.POSIXct(start_dt, tz = "UTC") + lubridate::hours(0), format = "%Y-%m-%d %H:%M:%S")
    utc_end_date <- format(as.POSIXct(end_dt, tz = "UTC") + lubridate::hours(0), format = "%Y-%m-%d %H:%M:%S")

    # Map over the location ids
    alldata <- site_loc_list %>% purrr::map(~hv_data_id(.,
                                                        start_time = utc_start_date,
                                                        end_time = utc_end_date,
                                                        token = api_token,
                                                        tz = "UTC"))

    # grab only locations with data (stored as a data frame)/drop 404 errors
    filtered <- purrr::keep(alldata, is.data.frame)

    if(length(filtered) != 0){

      # bind lists together (now that all are dataframes, we can just collate quickly)
      one_df <- dplyr::bind_rows(filtered) %>%
        dplyr::rename(id = Location,
                      parameter = Parameter,
                      units = Units) %>%
        dplyr::left_join(., site_loc, by = "id") %>%
        dplyr::mutate(site = tolower(site[i])) %>%
        dplyr::select(site, id, name, timestamp, parameter, value, units)

      readr::write_csv(one_df,
                       here::here(dump_dir,
                                  paste0(site[i],"-", stringr::str_replace(stringr::str_replace(substr(end_dt, 1, 16), "[:\\s]", "_"), ":", ""), ".csv")))

      # add an appending step here
      ## 1. Check if there is a file present in the dump_dir that also has site[i] in it

    } else {

      print(paste0("No data at ", site[i], " during this time frame"))

    }
  }

}
