# Extract conductivity data
cal_extract_conductivity_data <- function(div) {
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

  # Conductivity div metadata ----
  div_metadata <- div_tables[["sensor"]] %>%
    pivot_wider(names_from = X1, values_from = X2, names_repair = make_clean_names)  %>%
    rename(sensor_serial = serial_number) %>%
    mutate(sensor = make_clean_names(sensor))

  # In-situ calibration coefficients ----
  # Check if the table we expect for the calibration coefficients exists in the structure we expect...
  cal_slope_check <- cal_div_table_check(
    table_list = div_tables,
    table_name = "calibration_details",
    col_names = c("cell_constant")
  )

  cal_offset_check <- cal_div_table_check(
    table_list = div_tables,
    table_name = "calibration_details",
    col_names = c("offset")
  )

  if(cal_slope_check && cal_offset_check){
    calibration_coefs <- div_tables[["calibration_details"]] %>%
      pivot_wider(names_from = X1, values_from = X2, names_repair = make_clean_names) %>%
      separate_wider_delim(offset, delim = " ", names = c("offset", "units")) %>%
      select(slope = cell_constant, offset, units) %>%
      mutate(across(c(slope, offset), ~as.numeric(.x)))
  } else if(cal_slope_check && !cal_offset_check) {
    calibration_coefs <- div_tables[["calibration_details"]] %>%
      pivot_wider(names_from = X1, values_from = X2, names_repair = make_clean_names) %>%
      mutate(offset = NA_integer_, units = "µS/cm") %>%
      select(slope = cell_constant, offset, units) %>%
      mutate(across(c(slope, offset), ~as.numeric(.x)))
  } else {
    calibration_coefs <- NULL
  }

  # Driftr calibration input ----
  # Check if the table we expect for the calibration coefficient exists in the structure we expect...
  driftr_data_check_1 <- cal_div_table_check(
    table_list = div_tables,
    table_name = "pre_measurement",
    col_names = c("actual_conductivity", "specific_conductivity")
  )

  driftr_data_check_2 <- cal_div_table_check(
    table_list = div_tables,
    table_name = "post_measurement",
    col_names = c("actual_conductivity", "specific_conductivity")
  )

  if (driftr_data_check_1 && driftr_data_check_2) {

  calibration_1 <- div_tables[["pre_measurement"]] %>%
    pivot_wider(names_from = X1, values_from = X2, names_repair = make_clean_names)

  calibration_2 <- div_tables[["post_measurement"]] %>%
    pivot_wider(names_from = X1, values_from = X2, names_repair = make_clean_names)

  driftr_input <- tibble(point = c(1,2), bind_rows(calibration_1,calibration_2)) %>%
    mutate(
      units = word(specific_conductivity, 2),
      across(c(actual_conductivity, specific_conductivity), ~word(.x, 1)), # get rid of units in measurements
      across(c(actual_conductivity, specific_conductivity), ~gsub(",", "", .x)), # get rid of commas
      across(c(actual_conductivity, specific_conductivity), ~as.numeric(.x)) # convert to numeric
    )
  } else {
    driftr_input <- NULL
  }

  # Return ----
  cond_cal_info <- tibble(
    div_metadata,
    calibration_coefs = list(calibration_coefs),
    driftr_input = list(driftr_input)
  )

  return(cond_cal_info)
}


