setup_directories_from_upload <- function(uploaded_file_path, timezone = "MST"){

  testing_path = here("ver_tool_v1")
  data_path = in_progress_path
  all_path = all_data_path
  pre_verification_path_local = pre_verification_path
  intermediary_path_local = intermediary_path
  verified_path_local = verified_path
  raw_data_path_local = raw_data_path

  uploaded_file_type <- tools::file_ext(uploaded_file_path)
  org_file_type <- uploaded_file_type


  #Check if zip file and unzip to main working directory, this should be a data folder that was previously created by the user
  if(uploaded_file_type == "zip"){
    unzip(uploaded_file_path,exdir = testing_path )
    #Remove unneeded macoxs files
    macosx_path <- list.files(pattern = "MACOSX",full.names = TRUE)
    if (dir.exists(macosx_path) & !is_empty(macosx_path)){
      unlink(macosx_path, recursive = TRUE)  # Remove the folder and its contents
    }
    #check/ create all necessary folders
    if(!dir.exists(all_path)){
      dir.create(all_path)
    }
    #check to see if data folder has pre_verification_directory
    if(!dir.exists(pre_verification_path)){
      dir.create(pre_verification_path)
    }
    #check to see if data folder has intermediary_directory
    if(!dir.exists(intermediary_path)){
      dir.create(intermediary_path)
    }

    #check to see if data folder has verified_directory
    if(!dir.exists(verified_path)){
      dir.create(verified_path)
    }
    #create raw data folder
    if(!dir.exists(raw_data_path)){
      dir.create(raw_data_path)
      if(org_file_type != "zip"){
        file.copy(uploaded_file_path, raw_data_path)
      }
    }
    if(!dir.exists(meta_path)){
      dir.create(meta_path)
    }

    #if a user uploads a zip file, the data will be loaded from the zip file. If there is no data in the all path, check the raw data folder for a file to fill the all path
    if(length(list.files(all_path)) == 0){
      #check to see if raw_data folder contains a file
      if(length(list.files(raw_data_path)) > 0){
        #change uploaded data to data in raw_data folder to be processed downstream
        uploaded_file_path <- list.files(raw_data_path, full.names = TRUE)
      }else{
        return("No data found in all_data_directory and no raw_data folder found")
      }
    }else{
      #if no data exists in the pre_data folder but does exist in the all data folder, copy from all data directory
      if(length(list.files(pre_verification_path)) == 0){
        R.utils::copyDirectory(all_path, pre_verification_path)
        return("Data loaded from zip, no files found in pre directory, copied from all_data_directory")
        }
      }
      return("Data is loaded from zip file and contains pre and all_data_directory data.")
    }

  #check to see if data folder/sub folders exist
    if(!dir.exists(data_path)){
      dir.create(data_path)
    }
    #check to see if data folder has all_data_directory
    if(!dir.exists(all_path)){
      dir.create(all_path)
    }
    #check to see if data folder has pre_verification_directory
    if(!dir.exists(pre_verification_path)){
      dir.create(pre_verification_path)
    }
    #check to see if data folder has intermediary_directory
    if(!dir.exists(intermediary_path)){
      dir.create(intermediary_path)
    }
    #check to see if data folder has verified_directory
    if(!dir.exists(verified_path)){
      dir.create(verified_path)
    }
    #create raw data folder
    if(!dir.exists(raw_data_path)){
      dir.create(raw_data_path)
      if(org_file_type != "zip"){
        raw_data_filename <- file.path(raw_data_path, paste0("raw_data.", org_file_type))

        file.copy(uploaded_file_path, raw_data_filename)
      }
    }
    if(!dir.exists(meta_path)){
      dir.create(meta_path)
    }

  #this will read in the file based on its type and check to make sure there is no missing columns
  parse_raw_file <- function(file_path){
    #get type
    file_type <- tools::file_ext(file_path)

    if(file_type == "feather"){
      raw_data <- read_feather(file_path)
    } else if(file_type == "csv"){
      raw_data <- read_csv(file_path)
    } else if(file_type == "rds"){
      raw_data <- read_rds(file_path)
    } else if(file_type == "xlsx"){
      raw_data <- read_xlsx(file_path)
    }else if(file_type == "parquet"){
      raw_data <- read_parquet(file_path)
    } else{
      return(paste0("File type ", file_type, " not supported. Please upload a .feather, .csv, .rds, .xlsx, or .parquet file."))
    }
    #check to make sure site, parameter, DT round and mean all exist as columns
    if(!all(c("site", "parameter", "DT_round", "mean") %in% colnames(raw_data))){
      #find which column is missing
      missing_cols <- paste(setdiff(c("site", "parameter", "DT_round","mean"), colnames(raw_data)), collapse = ", ")
      return(paste0("Columns (", missing_cols, ") must exist in the file: ", file_path))
    }

    raw_data_parsed <- raw_data %>%
      #convert DT_round to date time
      mutate(DT_round = anytime(DT_round, tz = timezone))%>%
      ross.wq.tools::add_column_if_not_exists("flag", NA)

    return(raw_data_parsed)
  }

  #read in all data (either from uploaded file or from raw data folder)
  all_data_parsed <- map(uploaded_file_path, parse_raw_file)%>%
    bind_rows()
  #Taking all_data_parsed and formatting to save into all and pre directories
  site_params <- expand_grid(site = unique(all_data_parsed$site), parameter = unique(all_data_parsed$parameter))
  #create a list were each is named after the site and parameter combination
  flagged_list <- pmap(site_params, ~ all_data_parsed %>% filter(site == ..1, parameter == ..2))
  rm(all_data_parsed) #save some memory here
  names(flagged_list) <- paste(site_params$site, site_params$parameter, sep = "-")
  #remove empty lists
  flagged_list <- keep(flagged_list, ~nrow(.x) > 0)
  #add additional columns needed
  pre_processed_data <- map(.x = flagged_list, \(x)
                            x %>%
                              mutate(
                                mean_verified = mean, # if verification status pass or flag, this is mean. set default to mean so that it can be used in the app
                                is_verified = FALSE, # whether or not the data has been verified
                                verification_status = NA, # can only be pass/fail/skip
                                year = year(DT_round),
                                week = week(DT_round),
                                weekday = lubridate::wday(DT_round, week_start = 7),
                                y_w = paste(year, "-", week),
                                day = yday(DT_round),
                                y_d = paste(year, "-", day),
                                #based on Colorado Hydrology, should be removed for DS version
                                season = case_when(
                                  month(DT_round) %in% c(12, 1, 2, 3, 4) ~ "winter_baseflow",
                                  month(DT_round) %in% c(5, 6) ~ "snowmelt",
                                  month(DT_round) %in% c(7, 8, 9) ~ "monsoon",
                                  month(DT_round) %in% c(10, 11) ~ "fall_baseflow"
                                ),
                                #auto_flag = flag,
                                #merge auto_flag and mal_flag by ; into flag. don't include NAs
                                flag = case_when(
                                #both are flagged then combine by ;
                                !is.na(auto_flag)&!is.na(mal_flag) ~ paste(auto_flag, mal_flag, sep = ";"),
                                #only auto_flag is flagged, only use auto_flag
                                !is.na(auto_flag)&is.na(mal_flag) ~ auto_flag,
                                #only mal_flag is flagged, only use mal_flag
                                !is.na(mal_flag)&is.na(auto_flag) ~ mal_flag
                                         ),
                                user_flag = flag, # Keep Auto flags in flag column, user added flags/alterations live in user flag
                                brush_omit = FALSE, # User Brush Omit instances, default to FALSE
                                user = NA, #User initials
                                final_status = NA, # Final status of the data point after weekly decision
                                week_decision = NA,  #store weekly decision
                                is_finalized = FALSE) #Only TRUE if user submits final decision button
  )

  rm(flagged_list) #save some memory here

  iwalk(pre_processed_data, \(x, idx) {
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    data_hash <- digest::digest(x)
    new_filename <- glue("{idx}_{timestamp}_{data_hash}.parquet")
    write_parquet(x, file.path(all_path, new_filename))  # Ensure you use new_filename, not idx
  })

  R.utils::copyDirectory(all_path, pre_verification_path)
  return("Files saved to all and pre directory")
}
