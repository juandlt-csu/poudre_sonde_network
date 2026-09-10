sync_file_system <- function(){


  pre_dir_path <- pre_verification_path # pre_verification_path
  int_dir_path <- intermediary_path # intermediary_path
  ver_dir_path <- verified_path # verified_path

  pre_dir_names <- list.files(pre_dir_path, pattern = "\\.parquet$")
  int_dir_names <- list.files(int_dir_path, pattern = "\\.parquet$")
  ver_dir_names <- list.files(ver_dir_path, pattern = "\\.parquet$")

  # Duplicate file checks
  # Check if there are files that have duplicates, if there are any, run the file duplicate fixer
  duplicate_dfs_list <- list(pre_dir_names, int_dir_names, ver_dir_names) %>%
    set_names(c("pre", "int", "ver")) %>%
    keep(~length(.x) > 0) %>%
    imap(~{

      directory_path <- case_when(
        .y == "pre" ~ pre_dir_path,
        .y == "int" ~ int_dir_path,
        .y == "ver" ~ ver_dir_path
      )

      file_duplicate_alert_df <- map(.x, split_filename) %>%
        bind_rows()%>%
        mutate(full_file_path = here(directory_path, filename)) %>%
        select(!filename)

      return(file_duplicate_alert_df)
    }) %>%
    # keep those dataframes that have the duplicate_alert column
    keep(~"duplicate_alert" %in% names(.x))

  if (length(duplicate_dfs_list) != 0) {
    iwalk(duplicate_dfs_list, \(data, index){
      fix_duplicate_files(dataframe_with_duplicate_info = data, directory = index)
    })
  } else {
    print("No duplicates in any directories")
  }

  # Directory Checks
  # Check Pre-verification Directory
  # Keep if not in Verified Directory or Intermediary Directory, else if in
  # Verified Directory delete from pre Directory

  if (length(pre_dir_names) > 0) {
    walk(pre_dir_names, check_pre_ver_dir)
  }

  # Intermediary Directory
  # keep if skips in the data OR any(!is_verified), else move to Verified Directory
  if (length(int_dir_names) > 0) {
    walk(int_dir_names, check_int_ver_dir)
  }

  # Verified Directory
  # final destination for data, check that the data that is here should be here
  if (length(ver_dir_names) > 0) {
    walk(ver_dir_names, check_fin_ver_dir)
  }
}

# check pre verification directory (this is redundant, just double checking)
check_pre_ver_dir <- function(pre_file_name) {

  pre_dir_path <- pre_verification_path # pre_verification_path
  int_dir_path <- intermediary_path # intermediary_path
  ver_dir_path <- verified_path # verified_path

  pre_dir_names <- list.files(pre_dir_path, pattern = "\\.parquet$")
  int_dir_names <- list.files(int_dir_path, pattern = "\\.parquet$")
  ver_dir_names <- list.files(ver_dir_path, pattern = "\\.parquet$")

   tryCatch({
    file_info <- split_filename(pre_file_name)
    
    # Check if a file with same site/parameter exists in int or ver directory
    int_files <- int_dir_names[grepl(paste0("^", file_info$site, "-", file_info$parameter), int_dir_names)]
    ver_files <- ver_dir_names[grepl(paste0("^", file_info$site, "-", file_info$parameter), ver_dir_names)]
    
    if (length(int_files) > 0) {
      file.remove(here(pre_dir_path, pre_file_name))
      cat("removed file ", pre_file_name, " from ", pre_dir_path, " because it exists in intermediary\n")
    } else if (length(ver_files) > 0) {
      file.remove(here(pre_dir_path, pre_file_name))
      cat("removed file ", pre_file_name, " from ", pre_dir_path, " because it exists in verified\n")
    }
  }, error = function(e) {
    warning(sprintf("Error processing file %s: %s", pre_file_name, e$message))
  })
}

# Check intermediary verification directory
check_int_ver_dir <- function(int_file_name) {

  pre_dir_path <- pre_verification_path # pre_verification_path
  int_dir_path <- intermediary_path # intermediary_path
  ver_dir_path <- verified_path # verified_path

  pre_dir_names <- list.files(pre_dir_path, pattern = "\\.parquet$")
  int_dir_names <- list.files(int_dir_path, pattern = "\\.parquet$")
  ver_dir_names <- list.files(ver_dir_path, pattern = "\\.parquet$")

  tryCatch({
    #first check if the file is in the verified directory, if it is, remove it from the intermediary directory
    file_info <- split_filename(int_file_name)
    # Check for existing files with same site/parameter in verified dir
    ver_files <- list.files(ver_dir_path, pattern = "\\.parquet$")
    existing_files <- ver_files[grepl(paste0("^", file_info$site, "-", file_info$parameter), ver_files)]

    if (length(existing_files) > 0){
      file.remove(here(int_dir_path, int_file_name))
      cat("Removed ", int_file_name, " from ", int_dir_path, " because it already exists in verified directory\n")
    }

    # Only read in the verification columns
    df <- read_parquet(here(int_dir_path, int_file_name)) %>%
      select(verification_status, is_verified, is_finalized)

    # Check conditions without storing full data
    if (!any(df$verification_status == 'SKIP') & all(df$is_verified) & all(df$is_finalized)) {
      # Parse filename to get site and parameter
      file_info <- split_filename(int_file_name)
      # Check for existing files with same site/parameter in verified dir
      ver_files <- list.files(ver_dir_path, pattern = "\\.parquet$")
      existing_files <- ver_files[grepl(paste0("^", file_info$site, "-", file_info$parameter), ver_files)]


      # If existing files found, compare timestamps
      if (length(existing_files) > 0) {
        existing_info <- map_dfr(existing_files, split_filename)
        existing_times <- as.POSIXct(existing_info$datetime, format="%Y%m%d_%H%M%S")
        new_time <- as.POSIXct(file_info$datetime, format="%Y%m%d_%H%M%S")

        # Only proceed if new file is more recent
        if (new_time > max(existing_times)) {
          source_path <- here(int_dir_path, int_file_name)
          destination_path <- here(ver_dir_path, int_file_name)
          file.copy(source_path, destination_path, overwrite = TRUE)
          file.remove(source_path)
          # Remove older files
          file.remove(here(ver_dir_path, existing_files))
          cat("Updated file ", int_file_name, " in ", ver_dir_path, "\n")
        }
      } else {
        # No existing files, proceed with move
        source_path <- here(int_dir_path, int_file_name)
        destination_path <- here(ver_dir_path, int_file_name)
        file.copy(source_path, destination_path, overwrite = TRUE)
        file.remove(source_path)
        cat("Moved file ", int_file_name, " from ", int_dir_path, " to ", ver_dir_path, "\n")
      }
    }
  }, error = function(e) {
    warning(sprintf("Error processing file %s: %s", int_file_name, e$message))
  })
}

# Check final verified directory
check_fin_ver_dir <- function(ver_file_name) {

  pre_dir_path <- pre_verification_path # pre_verification_path
  int_dir_path <- intermediary_path # intermediary_path
  ver_dir_path <- verified_path # verified_path

  pre_dir_names <- list.files(pre_dir_path, pattern = "\\.parquet$")
  int_dir_names <- list.files(int_dir_path, pattern = "\\.parquet$")
  ver_dir_names <- list.files(ver_dir_path, pattern = "\\.parquet$")

  tryCatch({
    # Only read verification columns
    df <- read_parquet(here(ver_dir_path, ver_file_name)) %>%
      select(verification_status, is_verified)

    if (any(df$verification_status == 'SKIP') | any(!df$is_verified)) {
      msg <- sprintf("File %s contains SKIPS or non-verified data and should not be in verified directory", ver_file_name)
      cat(msg, "\n")
      stop(msg)
    } else {
      #cat("File", ver_file_name, "is correctly verified\n")
    }
  }, error = function(e) {
    warning(sprintf("Error processing file %s: %s", ver_file_name, e$message))
  })
}

# split file name for data synchronization
split_filename <- function(filename) {
  # Extract main parts
  parts <- strsplit(filename, "_\\d{8}_\\d{6}_")[[1]]
  datetime <- regmatches(filename, regexpr("\\d{8}_\\d{6}", filename))[[1]]

  # Check for extra content after hash
  hash_parts <- strsplit(parts[2], "-")[[1]]
  hash <- hash_parts[1]
  duplicate_alert <- if(length(hash_parts) > 1) paste(hash_parts[-1], collapse="-") else NULL

  # Split site/parameter
  if (grepl("Chl-a Fluorescence", parts[1])) {
    site_param <- strsplit(parts[1], "-(?=Chl-a)", perl=TRUE)[[1]]
  } else {
    site_param <- strsplit(parts[1], "-")[[1]]
  }

  #remove _Final if found in site_param[2]
  site_param[2] <- gsub("_FINAL", "", site_param[2])

  list(
    filename = filename,
    site = site_param[1],
    parameter = site_param[2],
    datetime = datetime,
    hash = hash,
    duplicate_alert = duplicate_alert
  )
}

# Function to fix file duplicates in a single directory
fix_duplicate_files <- function(dataframe_with_duplicate_info, directory){

  duplicate_alert_combos <- dataframe_with_duplicate_info %>%
    filter(!is.na(duplicate_alert)) %>%
    select(site, parameter) %>%
    distinct(site, parameter)

  map2(.x = duplicate_alert_combos$site,
       .y = duplicate_alert_combos$parameter,
       ~{
         # filter the data for the relevant files
         relevant_files <- dataframe_with_duplicate_info %>%
           filter(site == .x, parameter == .y)

         # Compare dates to get the most up to date file(s)...
         file <- relevant_files %>%
           mutate(datetime = as.POSIXct(datetime, format="%Y%m%d_%H%M%S")) %>%
           filter(datetime == max(datetime)) %>%
           # Sort by file size (larger files often more complete)
           arrange(desc(file.size(full_file_path))) %>%
           # Then by filename (if sizes equal) (shorter names on top)
           arrange(full_file_path) %>%
           slice(1)

         # if the file with a duplicate alert is the most recent file
         if (!is.na(file$duplicate_alert)) {

           # delete the other files...
           # TODO: archive, don't delete
           incorrect_file_paths <- relevant_files %>%
             filter(full_file_path != file$full_file_path) %>%
             pull(full_file_path)

           file.remove(incorrect_file_paths)

           # and rename the file ...
           dt_string <- format(file$datetime, "%Y%m%d_%H%M%S")
           new_file_name <- paste0(file$site, "-", file$parameter,
                                  "_", dt_string, "_", file$hash, ".parquet")

           file.rename(here(file$full_file_path), here(dirname(file$full_file_path), new_file_name))


         } else {
           # else, just keep the newest one and archive the rest
           incorrect_file_paths <- relevant_files %>%
             filter(full_file_path != file$full_file_path) %>%
             pull(full_file_path)

           file.remove(incorrect_file_paths)
         }
       })
  }

# In progress functions (pseudocode) ----
# TODO: Connect front end decisions to back end functions
# TODO: finish functions to move data around (different from the directory checks)
# TODO: finish function to update intermediary data
# TODO: current split_filename function does not take into account the fact that
  # the files in the pre directory do not have hashes or datetimes in the file names
# TODO: In the Shiny app, after moving a file, the app should reload the datasets
  # to reflect changes. So after these functions are called, the app needs to
  # reload the data from the directories.
# TODO: Need to implement locking mechanisms or checks to prevent data loss.

# Move file from pre to int directory after a user has made the decision to
# use a pre file and has made a weekly decision.
move_file_to_intermediary_directory <- function(pre_to_int_filename, pre_to_int_df) {

  pre_dir_path <- pre_verification_path # pre_verification_path
  int_dir_path <- intermediary_path # intermediary_path
  ver_dir_path <- verified_path # verified_path

  pre_dir_names <- list.files(pre_dir_path, pattern = "\\.parquet$")
  int_dir_names <- list.files(int_dir_path, pattern = "\\.parquet$")
  ver_dir_names <- list.files(ver_dir_path, pattern = "\\.parquet$")


  # Run initial sync to resolve pre-existing duplicates
  sync_file_system()

  # Extract identifiers and validate
  file_meta <- split_filename(pre_to_int_filename)
  site <- file_meta$site
  parameter <- file_meta$parameter

  # Generate conflict-safe filename
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  data_hash <- digest::digest(pre_to_int_df)
  new_filename <- glue("{site}-{parameter}_{timestamp}_{data_hash}.parquet")

  tryCatch({
    # Check for existing intermediary versions
    existing_files <- list.files(int_dir_path, pattern = glue("^{site}-{parameter}"))

    if(length(existing_files) > 0) {
      # Use sync logic to resolve conflicts
      duplicate_info <- split_filename(existing_files) %>%
        bind_rows() %>%
        mutate(full_file_path = file.path(int_dir_path, filename))

      fix_duplicate_files(duplicate_info, "int")
    }

    # Save new version
    write_parquet(pre_to_int_df, file.path(int_dir_path, new_filename))

    # Remove from pre-verification using sync logic
    pre_files <- list.files(pre_dir_path, pattern = glue("^{site}-{parameter}"))
    file.remove(file.path(pre_dir_path, pre_files))

    # Post-sync validation
    sync_file_system()

    return(new_filename)
  }, error = function(e) {
    # Rollback and sync
    if(file.exists(file.path(int_dir_path, new_filename))) {
      file.remove(file.path(int_dir_path, new_filename))
    }
    sync_file_system()
    return(list(success = FALSE, error = e$message))
  })
}

# Update intermediary data. Needs to be able to update the file and rename the file
# with an updated datetime and hash. file name structure is: `site-parameter_DT_hash`.
# the file should be an parquet
update_intermediary_data <- function(int_df_filename, updated_df) {

   pre_dir_path <- pre_verification_path # pre_verification_path
  int_dir_path <- intermediary_path # intermediary_path
  ver_dir_path <- verified_path # verified_path

  pre_dir_names <- list.files(pre_dir_path, pattern = "\\.parquet$")
  int_dir_names <- list.files(int_dir_path, pattern = "\\.parquet$")
  ver_dir_names <- list.files(ver_dir_path, pattern = "\\.parquet$")

  # Get current file metadata
  file_meta <- split_filename(int_df_filename)
  site <- file_meta$site
  parameter <- file_meta$parameter

  # Generate versioned filename
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  data_hash <- digest::digest(updated_df)
  new_filename <- glue("{site}-{parameter}_{timestamp}_{data_hash}.parquet")

  tryCatch({
    # Archive previous versions
    existing_files <- list.files(int_dir_path, pattern = glue("^{site}-{parameter}"), full.names = TRUE)
    file.remove(existing_files)

    # Save new version
    write_parquet(updated_df, file.path(int_dir_path, new_filename))

    # Run sync and duplicate check
    sync_file_system()

    # Validate verification status
    # if(all(updated_df$is_verified)&all(updated_df$is_finalized == TRUE) ) {
    #   check_int_ver_dir(new_filename)
    # }

    return(new_filename)
  }, error = function(e) {
    # Rollback to previous version
    if(file.exists(file.path(int_dir_path, new_filename))) {
      file.remove(file.path(int_dir_path, new_filename))
    }
    sync_file_system()
    return(list(success = FALSE, error = e$message))
  })
}

# Move file from int to final directory after a user has made the final decision
# on an int file.
move_file_to_verified_directory <- function(int_to_fin_filename, int_to_fin_df) {

  pre_dir_path <- pre_verification_path # pre_verification_path
  int_dir_path <- intermediary_path # intermediary_path
  ver_dir_path <- verified_path # verified_path

  pre_dir_names <- list.files(pre_dir_path, pattern = "\\.parquet$")
  int_dir_names <- list.files(int_dir_path, pattern = "\\.parquet$")
  ver_dir_names <- list.files(ver_dir_path, pattern = "\\.parquet$")

  # Pre-move validation
  if(any(int_to_fin_df$verification_status == 'SKIP') ||
     any(!int_to_fin_df$is_verified)) {
    stop("Cannot move unverified data to final directory")
  }

  # Generate final filename
  file_meta <- split_filename(int_to_fin_filename)
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  data_hash <- digest::digest(int_to_fin_df)
  new_filename <- glue("{file_meta$site}-{file_meta$parameter}_FINAL_{timestamp}_{data_hash}.parquet")

  tryCatch({
    # Check for existing final versions
    ver_files <- list.files(ver_dir_path, pattern = glue("^{file_meta$site}-{file_meta$parameter}"))

    if(length(ver_files) > 0) {
      # Use sync conflict resolution
      duplicate_info <- split_filename(ver_files) %>%
        bind_rows() %>%
        mutate(full_file_path = file.path(ver_dir_path, filename))%>%
        mutate(duplicate_alert = TRUE)

      fix_duplicate_files(dataframe_with_duplicate_info = duplicate_info, directory = "ver")
    }

    # Move and validate
    write_parquet(int_to_fin_df, file.path(ver_dir_path, new_filename))
    #remove all old int files that match the site-param
    int_files <- list.files(int_dir_path, pattern = glue("^{file_meta$site}-{file_meta$parameter}"), full.names = T)
    file.remove(int_files)

    # Post-move sync and validation
    sync_file_system()
#Q: Reduncant? check_fin_ver_dir is already called in sync_file_system call. Doesn't result in any errors, just repeats in console.
    check_fin_ver_dir(new_filename)

    return(new_filename)
  }, error = function(e) {
    # Rollback move
    if(file.exists(file.path(ver_dir_path, new_filename))) {
      file.remove(file.path(ver_dir_path, new_filename))
    }
    sync_file_system()
    return(list(success = FALSE, error = e$message))
  })
}


