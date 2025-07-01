cal_extract_markup_data <- function(html_markup, calibration_site) {

  # Establish sensors we are retrieiving data for
  sensor_list <- c("chlorophyll_a",
                   "conductivity",
                   "fdom",
                   "p_h_orp",
                   "pressure",
                   "rdo",
                   "turbidity")

  # Grab the metadata for the overall calibration
  cal_metadata <- html_markup %>%
    html_elements("table") %>%
    html_table() %>%
    pluck(1) %>%
    pivot_wider(names_from = X1, values_from = X2) %>%
    clean_names() %>%
    mutate(instrument = make_clean_names(instrument)) %>%
    rename(sonde_serial = serial_number)

  # Grab the div metadata for the html_markup file
  # Each div is for one parameter
  html_divs <- html_markup %>%
    html_elements("div")

  html_div_info <- html_divs %>%
    map_dfr(function(div){

      sensor <- div %>%
        html_elements("table") %>%
        html_table() %>%
        pluck(1) %>%
        pivot_wider(names_from = X1, values_from = X2, names_repair = make_clean_names) %>%
        mutate(sensor = make_clean_names(sensor)) %>%
        pull(sensor)

      if (sensor %in% sensor_list){
        div_table_data <- switch(
          EXPR = sensor,
          "conductivity" = cal_extract_conductivity_data(div),
          "rdo" = cal_extract_rdo_data(div),
          "p_h_orp" = cal_extract_ph_orp_data(div),
          "pressure" = cal_extract_pressure_data(div),
          "turbidity" = cal_extract_turbidity_data(div),
          "fdom" = cal_extract_fdom_data(div),
          "chlorophyll_a" = cal_extract_chla_data(div)
        )
      } else {
        return()
      }

      return(div_table_data)

    })

  # Generate the metadata for the whole file
  metadata <- tibble(
    calibration_site,
    cal_metadata,
    html_div_info
  )

  return(metadata)
}





