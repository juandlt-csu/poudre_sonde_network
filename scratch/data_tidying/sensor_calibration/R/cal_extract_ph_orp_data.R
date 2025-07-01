# Extract pH and ORP data
cal_extract_ph_orp_data <- function(div) {

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
      mutate(sensor = str_split(sensor, "/")) %>%
      unnest(sensor) %>%
      mutate(
        sensor = make_clean_names(sensor),
        calibration_coefs = NULL,
        driftr_input = NULL
      )
    return(div_metadata)
  }
  # Else, continue with the function...

  # pH and ORP metadata ----
  div_metadata <- div_tables[["sensor"]] %>%
    pivot_wider(names_from = X1, values_from = X2, names_repair = make_clean_names)  %>%
    rename(sensor_serial = serial_number) %>%
    mutate(sensor = str_split(sensor, "/")) %>%
    unnest(sensor) %>%
    mutate(sensor = make_clean_names(sensor))

  # pH ----
  ## In-situ calibration coefficients ----
  # Check if the table we expect for the calibration coefficients exists in the structure we expect...
  cal_coef_check_ph_1 <- cal_div_table_check(
    table_list = div_tables,
    table_name = "slope_and_offset_1",
    col_names = c("slope", "offset")
  )

  cal_coef_check_ph_2 <- cal_div_table_check(
    table_list = div_tables,
    table_name = "slope_and_offset_2",
    col_names = c("slope", "offset")
  )

  if(cal_coef_check_ph_1 && cal_coef_check_ph_2){

    calibration_coefs_ph <- div_tables[c("slope_and_offset_1", "slope_and_offset_2")] %>%
      map_dfr(function(ph_coef_df){
        pivot_df <- ph_coef_df %>%
          pivot_wider(names_from = X1, values_from = X2, names_repair = make_clean_names)
      }) %>%
      separate_wider_delim(slope, delim = " ", names = c("slope", "slope_units")) %>%
      separate_wider_delim(offset, delim = " ", names = c("offset", "offset_units")) %>%
      mutate(
        point = c(1,2),
        across(c(slope, offset), ~as.numeric(.x))
      ) %>%
      relocate(point, .before = "slope")

  } else {
    calibration_coefs_ph <- NULL
  }

  ## Driftr calibration input ----
  # Check if the tables we expect for the driftr input exists in the structure we expect...
  calibration_table_name_list <- c("calibration_point_1", "calibration_point_2", "calibration_point_3")
  measurement_table_name_list <- c("pre_measurement", "post_measurement",
                                   "pre_measurement_2", "post_measurement_2",
                                   "pre_measurement_3", "post_measurement_3")

  calibration_column_name_list <- c("p_h_of_buffer", "p_h_m_v", "temperature")
  measurement_column_name_list <- c("p_h", "p_h_m_v")

  calibration_table_check <- calibration_table_name_list %>%
    map(\(table_name){
      cal_div_table_check(
        table_list = div_tables,
        table_name = table_name,
        col_names = calibration_column_name_list
      )
    }) %>%
    unlist()

  measurement_table_check <- measurement_table_name_list %>%
    map(\(table_name){
      cal_div_table_check(
        table_list = div_tables,
        table_name = table_name,
        col_names = measurement_column_name_list
      )
    }) %>%
    unlist()

  if(all(calibration_table_check) && all(measurement_table_check)){

    # Calibration function specifically for pH data
    driftr_ph_cal_extractor <- function(calibration_n, point_table, pre_table, post_table){

      point_df <- point_table %>%
        pivot_wider(names_from = X1, values_from = X2, names_repair = make_clean_names) %>%
        rename(p_h = p_h_of_buffer)

      pre_df <- pre_table %>%
        pivot_wider(names_from = X1, values_from = X2, names_repair = make_clean_names)

      post_df <- post_table %>%
        pivot_wider(names_from = X1, values_from = X2, names_repair = make_clean_names)

      driftr_cal_info <- tibble(
        point = c(calibration_n),
        type = c(0,1,2), # 0 is target, 1 is pre-measurement, 2 is post-measurement
        bind_rows(point_df, pre_df, post_df)
      ) %>%
        tidyr::fill(temperature, .direction = "down") %>%
        separate_wider_delim(p_h, delim = " ", names = c("p_h", "p_h_units")) %>%
        separate_wider_delim(p_h_m_v, delim = " ", names = c("p_h_m_v", "p_h_m_v_units")) %>%
        separate_wider_delim(temperature, delim = " ", names = c("temp", "temp_units"))

      return(driftr_cal_info)
    }

    calibration_1_ph <- driftr_ph_cal_extractor(
      calibration_n = 1,
      point_table = div_tables[["calibration_point_1"]],
      pre_table = div_tables[["pre_measurement"]],
      post_table = div_tables[["post_measurement"]]
    )

    calibration_2_ph <- driftr_ph_cal_extractor(
      calibration_n = 2,
      point_table = div_tables[["calibration_point_2"]],
      pre_table = div_tables[["pre_measurement_2"]],
      post_table = div_tables[["post_measurement_2"]]
    )

    calibration_3_ph <- driftr_ph_cal_extractor(
      calibration_n = 3,
      point_table = div_tables[["calibration_point_3"]],
      pre_table = div_tables[["pre_measurement_3"]],
      post_table = div_tables[["post_measurement_3"]]
    )

    driftr_input_ph <- bind_rows(calibration_1_ph,calibration_2_ph,calibration_3_ph)
  } else {
    driftr_input_ph <- NULL
  }

  # ORP ----
  ## In-situ calibration coefficients ----
  # Check if the table we expect for the calibration coefficients exists in the structure we expect...
  cal_coef_check_orp <- cal_div_table_check(
    table_list = div_tables,
    table_name = "orp",
    col_names = c("orp_solution", "offset", "temperature")
  )

  if(cal_coef_check_orp){
    calibration_coefs_orp <- div_tables[["orp"]] %>%
      pivot_wider(names_from = X1, values_from = X2, names_repair = make_clean_names) %>%
      mutate(
        slope = 1,
        orp_solution = make_clean_names(orp_solution)
      ) %>%
      separate_wider_delim(offset, delim = " ", names = c("offset", "offset_units")) %>%
      separate_wider_delim(temperature, delim = " ", names = c("temperature", "temperature_units")) %>%
      mutate(point = 1) %>%
      select(point, slope, offset, offset_units, temperature, temperature_units, orp_solution)
  } else {
    calibration_coefs_orp <- NULL
  }

  ## Driftr calibration input ----
  # Check if the table we expect for the driftr input exists in the structure we expect...
  driftr_data_check_orp <- cal_div_table_check(
    table_list = div_tables,
    table_name = "orp",
    col_names = c("pre_measurement", "post_measurement")
  )

  if (driftr_data_check_orp){
    calibration_orp <- div_tables[["orp"]] %>%
      pivot_wider(names_from = X1, values_from = X2, names_repair = make_clean_names) %>%
      select(pre_measurement, post_measurement)

    driftr_input_orp <- tibble(point = c(1), bind_rows(calibration_orp)) %>%
      mutate(
        units = word(post_measurement, 2),
        across(c(pre_measurement, post_measurement), ~word(.x, 1)), # get rid of units in measurements
        across(c(pre_measurement, post_measurement), ~gsub(",", "", .x)), # get rid of commas
        across(c(pre_measurement, post_measurement), ~as.numeric(.x)) # convert to numeric
      )
  } else {
    driftr_input_orp <- NULL
  }

  # Return ----
  ph_orp_info <- tibble(
    div_metadata,
    calibration_coefs = list(calibration_coefs_ph, calibration_coefs_orp),
    driftr_input = list(driftr_input_ph, driftr_input_orp)
  )

  return(ph_orp_info)
}
