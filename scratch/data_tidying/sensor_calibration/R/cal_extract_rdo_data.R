# Extract RDO data
cal_extract_rdo_data <- function(div) {
  # Instantiate a named div tables list ----
  div_tables <- div %>%
    html_elements("table") %>%
    discard(\(x) length(html_elements(x, "caption")) == 0) %>% # remove tables with no captions (empty placeholder tables)
    html_table() %>% # extract the tabtle data
    set_names( # use the table captions to name the tables in the div_tables list
      {div %>%
        html_elements("caption") %>%
        html_text() %>%
        make_clean_names()
      }
    )

  # Initial div data availability check ----
  div_check <- cal_div_table_check(
    table_list = div_tables,
    table_name = "sensor",
    col_names = c("sensor", "serial_number", "last_calibrated")
  )

  # If there is only 1 table in the div exit early and return metadata ...
  if (!div_check) {
    div_metadata <- div_tables[["sensor"]] %>%
      pivot_wider(names_from = X1, values_from = X2, names_repair = make_clean_names) %>%
      rename(sensor_serial = serial_number) %>%
      mutate(
        sensor = make_clean_names(sensor),
        calibration_coefs = NULL,
        driftr_input = NULL
      )
    return(div_metadata)
  }
  # Else, continue with the function...

  # RDO div metadata ----
  div_metadata <- div_tables[["sensor"]] %>%
    pivot_wider(names_from = X1, values_from = X2, names_repair = make_clean_names)  %>%
    rename(sensor_serial = serial_number) %>%
    mutate(sensor = make_clean_names(sensor))

  # In-situ calibration coefficients ----
  # Check if the table we expect for the calibration coefficients exists in the structure we expect...
  cal_coef_check <- cal_div_table_check(
    table_list = div_tables,
    table_name = "calibration_details",
    col_names = c("slope", "offset")
  )

  if(cal_coef_check){
    calibration_coefs <- div_tables[["calibration_details"]] %>%
      pivot_wider(names_from = X1, values_from = X2, names_repair = make_clean_names) %>%
      separate_wider_delim(offset, delim = " ", names = c("offset", "units")) %>%
      mutate(across(c(slope, offset), ~as.numeric(.x)))
  } else {
    calibration_coefs <- NULL
  }

  # Driftr calibration input ----
  # Check if the table we expect for the driftr inputs exists in the structure that we expect...
  driftr_data_check <- cal_div_table_check(
    table_list = div_tables,
    table_name = "calibration_point_100_percent",
    col_names = c("pre_measurement", "post_measurement")
  )

  if(driftr_data_check){
    calibration <- div_tables[["calibration_point_100_percent"]] %>%
      pivot_wider(names_from = X1, values_from = X2, names_repair = make_clean_names) %>%
      select(pre_measurement, post_measurement)

    driftr_input <- tibble(point = c(1), bind_rows(calibration)) %>%
      mutate(
        units = word(post_measurement, 2),
        across(c(pre_measurement, post_measurement), ~word(.x, 1)), # get rid of units in measurements
        across(c(pre_measurement, post_measurement), ~gsub(",", "", .x)), # get rid of commas
        across(c(pre_measurement, post_measurement), ~as.numeric(.x)) # convert to numeric
      )
  } else {
    driftr_input <- NULL
  }

  # Return ----
  rdo_cal_info <- tibble(
    div_metadata,
    calibration_coefs = list(calibration_coefs),
    driftr_input = list(driftr_input)
  )

  return(rdo_cal_info)
}
