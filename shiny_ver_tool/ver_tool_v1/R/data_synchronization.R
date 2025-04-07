sync_file_system <- function() {

  pre_dir_path <- here("shiny_ver_tool", "ver_tool_v1", "data", "scratch_data")# pre_verification_path
  int_dir_path <- here("shiny_ver_tool", "ver_tool_v1", "data", "scratch_data")# intermediary_path
  ver_dir_path <- here("shiny_ver_tool", "ver_tool_v1", "data", "scratch_data")# verified_path

  pre_dir_names <- list.files(pre_dir_path)
  int_dir_names <- list.files(int_dir_path)
  ver_dir_names <- list.files(ver_dir_path)

  # Duplicate file checks
  # Check if there are files that have duplicates, if there are any, run the file duplicate fixer
  duplicate_dfs_list <- list(pre_dir_names, int_dir_names, ver_dir_names) %>%
    set_names(c("pre", "int", "ver")) %>%
    imap(~{

      directory_path <- case_when(
        .y == "pre" ~ pre_dir_path,
        .y == "int" ~ int_dir_path,
        .y == "ver" ~ ver_dir_path
      )

      file_duplicate_alert_df <- map(.x, split_filename) %>%
        bind_rows() %>%
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
  tryCatch({
    if (length(int_dir_names) > 0 & pre_file_name %in% int_dir_names) {
      file.remove(here(pre_dir_path, pre_file_name))
      cat("removed file ", i, " from ", int_dir_path, "\n")
    }
    if (length(ver_dir_names) > 0 & pre_file_name %in% ver_dir_names) {
      file.remove(here(pre_dir_path, pre_file_name))
      cat("removed file ", pre_file_name, " from ", ver_dir_path, "\n")
    }
  }, error = function(e) {
    warning(sprintf("Error processing file %s: %s", pre_file_name, e$message))
  })
}

# Check intermediary verification directory
check_int_ver_dir <- function(int_file_name) {
  tryCatch({
    # Only read in the verification columns
    df <- read_rds(here(int_dir_path, int_file_name)) %>%
      select(verification_status, is_verified)

    # Check conditions without storing full data
    if (!any(df$verification_status == 'SKIP') & all(df$is_verified)) {
      # Parse filename to get site and parameter
      file_info <- split_filename(int_file_name)

      # Check for existing files with same site/parameter in verified dir
      ver_files <- list.files(ver_dir_path)
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
  tryCatch({
    # Only read verification columns
    df <- read_rds(here(ver_dir_path, ver_file_name)) %>%
      select(verification_status, is_verified)

    if (any(df$verification_status == 'SKIP') | any(!df$is_verified)) {
      msg <- sprintf("File %s contains SKIPS or non-verified data and should not be in verified directory", ver_file_name)
      cat(msg, "\n")
      stop(msg)
    } else {
      cat("File", ver_file_name, "is correctly verified\n")
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
           new_file_name <- paste(file$site, file$parameter,
                                  dt_string, file$hash, sep="_")

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
move_file_to_intermediary_directory <- function(pre_to_int_df) {
  # Run initial sync to resolve pre-existing duplicates
  sync_file_system()

  # Extract identifiers and validate
  file_meta <- split_filename(pre_to_int_df$source_file[1]) # TODO: this won't work on pre directory file names
  site <- file_meta$site
  parameter <- file_meta$parameter

  # Generate conflict-safe filename
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  data_hash <- digest::digest(pre_to_int_df)
  new_filename <- glue("{site}_{parameter}_{timestamp}_{data_hash}.rds")

  tryCatch({
    # Check for existing intermediary versions
    existing_files <- list.files(int_dir_path, pattern = glue("^{site}_{parameter}"))

    if(length(existing_files) > 0) {
      # Use sync logic to resolve conflicts
      duplicate_info <- split_filename(existing_files) %>%
        bind_rows() %>%
        mutate(full_file_path = file.path(int_dir_path, filename))

      fix_duplicate_files(duplicate_info, "int")
    }

    # Save new version
    saveRDS(pre_to_int_df, file.path(int_dir_path, new_filename))

    # Remove from pre-verification using sync logic
    pre_files <- list.files(pre_dir_path, pattern = glue("^{site}-{parameter}"))
    file.remove(file.path(pre_dir_path, pre_files))

    # Post-sync validation
    sync_file_system()

    return(list(success = TRUE, new_file = new_filename))
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
# with an updated datetime and hash. file name structure is: `site_parameter_DT_hash`.
# the file should be an RDS
update_intermediary_data <- function(updated_df) {
  # Get current file metadata
  file_meta <- split_filename(updated_df$source_file[1])
  site <- file_meta$site
  parameter <- file_meta$parameter

  # Generate versioned filename
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  data_hash <- digest::digest(updated_df)
  new_filename <- glue("{site}_{parameter}_{timestamp}_{data_hash}.rds")

  tryCatch({
    # Archive previous versions but keep last 3
    existing_files <- list.files(int_dir_path, pattern = glue("^{site}_{parameter}"))
    if(length(existing_files) > 3) {
      to_archive <- existing_files[1:(length(existing_files)-3)]
      archive_files(to_archive, "int_archive")
    }

    # Save new version
    saveRDS(updated_df, file.path(int_dir_path, new_filename))

    # Run sync and duplicate check
    sync_file_system()

    # Validate verification status
    if(all(updated_df$is_verified)) {
      check_int_ver_dir(new_filename)
    }

    return(list(success = TRUE, new_version = new_filename))
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
move_file_to_verified_directory <- function(int_to_fin_df) {
  # Pre-move validation
  if(any(int_to_fin_df$verification_status == 'SKIP') ||
     any(!int_to_fin_df$is_verified)) {
    stop("Cannot move unverified data to final directory")
  }

  # Generate final filename
  file_meta <- split_filename(int_to_fin_df$source_file[1])
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  data_hash <- digest::digest(int_to_fin_df)
  new_filename <- glue("{file_meta$site}_{file_meta$parameter}_FINAL_{timestamp}_{data_hash}.rds")

  tryCatch({
    # Check for existing final versions
    ver_files <- list.files(ver_dir_path, pattern = glue("^{file_meta$site}_{file_meta$parameter}"))

    if(length(ver_files) > 0) {
      # Use sync conflict resolution
      duplicate_info <- split_filename(ver_files) %>%
        bind_rows() %>%
        mutate(full_file_path = file.path(ver_dir_path, filename))

      fix_duplicate_files(duplicate_info, "ver")
    }

    # Move and validate
    file.copy(file.path(int_dir_path, int_to_fin_df$source_file[1]),
              file.path(ver_dir_path, new_filename))
    file.remove(file.path(int_dir_path, int_to_fin_df$source_file[1]))

    # Post-move sync and validation
    sync_file_system()
    check_fin_ver_dir(new_filename)

    return(list(success = TRUE, final_file = new_filename))
  }, error = function(e) {
    # Rollback move
    if(file.exists(file.path(ver_dir_path, new_filename))) {
      file.remove(file.path(ver_dir_path, new_filename))
    }
    sync_file_system()
    return(list(success = FALSE, error = e$message))
  })
}


