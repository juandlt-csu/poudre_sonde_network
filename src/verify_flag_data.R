verify_flag_data <- function(df_list_arg,
                             site_arg,
                             parameter_arg, # do we want this to be able to be null eventually?
                             flag_arg = NULL,
                             weekly_plot_object) {

  # set layer
  layer <- "weekly"

  #quit

  if (!QUIT) {
    # verification decision lists ----
    simple_decision_list <- c("PASS ALL", "FAIL NONE",
                              "PASS VALID", "FAIL FLAGGED",
                              "PASS FLAGGED", "FAIL VALID",
                              "PASS NONE", "FAIL ALL",
                              "SKIP", "QUIT")

    inspect_decision_list <- c("INSPECT ALL", "INSPECT VALID",
                               "INSPECT FLAGGED", "INSPECT SOME")

    # initialize site_param_df ----
    site_param_df <- df_list_arg[[paste0(site_arg, "-", parameter_arg)]]

    # get data from initial weekly plot object ----
    df_weekly_data <- ggplot_build(weekly_plot_object)$plot$data %>%
      filter(site == site_arg) %>%
      select(intersect(names(.), names(site_param_df)))

    df_weekly_name <- ggplot_build(weekly_plot_object)$plot$labels$title

    # initialize list to contain updated dfs OUTPUT, VERY IMPORTANT ----
    altered_df_list <- list()

    # initialize week list ----
    weekday_list <- c(1,2,3,4,5,6,7)

    # Make sure that non-removable flags are still in place ----
    non_removable_flags <- paste0(c("sensor biofouling","sonde burial",
                                    "sensor malfunction","frozen","missing data",
                                    "sonde not employed","sonde unsubmerged",
                                    "depth calibration malfunction"),
                                  collapse = "|")

    # Which data are you working on: ----
    cat("Currently verifying: ", df_weekly_name, "data. \n")

    print(weekly_plot_object)

    # Get weekly decision(s) ----
    verification_decision <- get_verification_decision(
      paste("Would you like to (pass/fail/skip/inspect/quit) ", df_weekly_name, "?\n(layer~", layer, ")(px/fx/sx/ix/q): ") # this should be weekly data
    )

    if (verification_decision %in% simple_decision_list){
      altered_df_list <- alter_verification_column(verification_decision, df_weekly_data)
    }

    if (verification_decision %in% inspect_decision_list){
      # set layer
      layer <- "daily"

      # read in the data relevant to the current site-parameter that we are doing
      # this will be used to make the inspect decision plots. pull this data from the working directory
      # site_param_df <<- df_list[[paste0(site_arg, "-", parameter_arg)]] # if parameter arg is null then pull it all in

      # inspect all weekly ----
      if (verification_decision == "INSPECT ALL"){

        for (i in unique(df_weekly_data$weekday)) {

          daily_plot_data <- site_param_df %>%
            filter((y_w %in% df_weekly_data$y_w) & (weekday == i))

          daily_plot_object <- generate_daily_plot(plot_data_arg = daily_plot_data,
                                                   df_list_arg = df_list_arg,
                                                   site_arg = site_arg,
                                                   parameter_arg = parameter_arg)

          # weekly_plot_object
          weekly_plot_object <- generate_supplemental_weekly_plot(daily_plot_data_arg = daily_plot_data,
                                                                  df_list_arg = df_list_arg,
                                                                  site_arg = site_arg,
                                                                  parameter_arg = parameter_arg) # need to make sure that these functions can read in information from outside of them

          # print ggarrange daily and weekly plots
          print(ggarrange(daily_plot_object, weekly_plot_object, nrow = 2, ncol = 1))

          # Generate the prompt for the user to answer
          day_dt <- as.character(head(daily_plot_data$DT_round, 1)) # reformat so that this reads as a date
          day_choice_prompt <- paste("Would you like to (pass/fail/skip/inspect) ", day_dt, "?\n(layer~", layer, ")(px/fx/sx/ix): ")
          daily_verification_decision <- get_verification_decision(day_choice_prompt)

          while (layer == "daily" & daily_verification_decision == "QUIT") {
            cat("You cannot quit during this step. Please choose a different option.\n")
            daily_verification_decision <- get_verification_decision(day_choice_prompt)
          }

          if (daily_verification_decision %in% simple_decision_list){
            altered_df_list[[i]] <- alter_verification_column(daily_verification_decision, daily_plot_data) # this should get put in a list
          }

          if (daily_verification_decision %in% inspect_decision_list){
            altered_df_list[[i]] <- verify_subdaily(daily_verification_decision_arg = daily_verification_decision,
                                                    daily_plot_data_arg = daily_plot_data,
                                                    daily_plot_object_arg = daily_plot_object,
                                                    weekly_plot_object_arg = weekly_plot_object
            )
          }
        }
      }

      # inspect valid weekly ----
      if (verification_decision == "INSPECT VALID") {

        # get which days to inspect
        days_to_inspect <- unique(df_weekly_data %>%
                                    filter(is.na(flag)) %>%
                                    pull(weekday))

        new_weekday_list <- weekday_list[weekday_list >= min(df_weekly_data$weekday) & weekday_list <= max(df_weekly_data$weekday)]
        non_days_to_inspect <- new_weekday_list[!new_weekday_list %in% days_to_inspect]
        # for each non-inspect day decision:
        if (length(non_days_to_inspect) > 0) {

          # get non-inspect day decision:
          noninspect_prompt_text = paste("Would you like to (pass/fail/skip) non-inspect days?\n(layer~", layer, ")(px/fx/sx): ")
          non_inspect_verification_decision <- get_verification_decision(prompt_text = noninspect_prompt_text)

          # Handle the case when the user selects an "INSPECT" or "QUIT" option for non-inspect days
          while (layer == "daily" & non_inspect_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME", "QUIT")) {
            if(non_inspect_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME")) {
              cat("You cannot inspect non-inspect days. Please choose a different option.\n")
              non_inspect_verification_decision <- get_verification_decision(prompt_text = noninspect_prompt_text)
            }
            if(non_inspect_verification_decision == "QUIT") {
              cat("You cannot quit during this step. Please choose a different option.\n")
              non_inspect_verification_decision <- get_verification_decision(prompt_text = noninspect_prompt_text)
            }
          }

          for (i in non_days_to_inspect) {
            daily_plot_data <- site_param_df %>%
              filter((y_w %in% df_weekly_data$y_w) & (weekday == i))

            altered_df_list[[i]] <- alter_verification_column(non_inspect_verification_decision, daily_plot_data)
          }
        }

        for (i in days_to_inspect) {
          daily_plot_data <- site_param_df %>%
            filter((y_w %in% df_weekly_data$y_w) & (weekday == i))

          daily_plot_object <- generate_daily_plot(plot_data_arg = daily_plot_data,
                                                   df_list_arg = df_list_arg,
                                                   site_arg = site_arg,
                                                   parameter_arg = parameter_arg)

          # weekly_plot_object
          weekly_plot_object <- generate_supplemental_weekly_plot(daily_plot_data_arg = daily_plot_data,
                                                                  df_list_arg = df_list_arg,
                                                                  site_arg = site_arg,
                                                                  parameter_arg = parameter_arg) # need to make sure that these functions can read in information from outside of them

          # print ggarrange daily and weekly plots
          print(ggarrange(daily_plot_object, weekly_plot_object, nrow = 2, ncol = 1))

          # Generate the prompt for the user to answer
          day_dt <- as.character(head(daily_plot_data$DT_round, 1)) # reformat so that this reads as a date
          day_choice_prompt <- paste("Would you like to (pass/fail/skip/inspect) ", day_dt, "?\n(layer~", layer, ")(px/fx/sx/ix): ")
          daily_verification_decision <- get_verification_decision(day_choice_prompt)

          while (daily_verification_decision == "QUIT") {
            cat("You cannot quit during this step. Please choose a different option.\n")
            daily_verification_decision <- get_verification_decision(day_choice_prompt)
          }

          if (daily_verification_decision %in% simple_decision_list){
            altered_df_list[[i]] <- alter_verification_column(daily_verification_decision, daily_plot_data) # this should get put in a list
          }

          if (daily_verification_decision %in% inspect_decision_list){
            altered_df_list[[i]] <- verify_subdaily(daily_verification_decision_arg = daily_verification_decision,
                                                    daily_plot_data_arg = daily_plot_data,
                                                    daily_plot_object_arg = daily_plot_object,
                                                    weekly_plot_object_arg = weekly_plot_object
            )
          }
        }
      }

      # inspect flagged weekly ----
      if (verification_decision == "INSPECT FLAGGED") {

        # get which days to inspect
        days_to_inspect <- unique(df_weekly_data %>%
                                    filter(!is.na(flag)) %>%
                                    pull(weekday))

        new_weekday_list <- weekday_list[weekday_list >= min(df_weekly_data$weekday) & weekday_list <= max(df_weekly_data$weekday)]
        non_days_to_inspect <- new_weekday_list[!new_weekday_list %in% days_to_inspect]

        # for each non-inspect day decision:
        if (length(non_days_to_inspect) > 0) {

          # get non-inspect day decision:
          noninspect_prompt_text = paste("Would you like to (pass/fail/skip) non-inspect days?\n(layer~", layer, ")(px/fx/sx): ")
          non_inspect_verification_decision <- get_verification_decision(prompt_text = noninspect_prompt_text)

          # Handle the case when the user selects an "INSPECT" or "QUIT" option for non-inspect days
          while (layer == "daily" & non_inspect_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME", "QUIT")) {
            if(non_inspect_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME")) {
              cat("You cannot inspect non-inspect days. Please choose a different option.\n")
              non_inspect_verification_decision <- get_verification_decision(prompt_text = noninspect_prompt_text)
            }
            if(non_inspect_verification_decision == "QUIT") {
              cat("You cannot quit during this step. Please choose a different option.\n")
              non_inspect_verification_decision <- get_verification_decision(prompt_text = noninspect_prompt_text)
            }
          }

          for (i in non_days_to_inspect) {
            daily_plot_data <- site_param_df %>%
              filter((y_w %in% df_weekly_data$y_w) & (weekday == i))

            altered_df_list[[i]] <- alter_verification_column(non_inspect_verification_decision, daily_plot_data)
          }
        }

        for (i in days_to_inspect) {
          daily_plot_data <- site_param_df %>%
            filter((y_w %in% df_weekly_data$y_w) & (weekday == i))

          daily_plot_object <- generate_daily_plot(plot_data_arg = daily_plot_data,
                                                   df_list_arg = df_list_arg,
                                                   site_arg = site_arg,
                                                   parameter_arg = parameter_arg)


          # weekly_plot_object
          weekly_plot_object <- generate_supplemental_weekly_plot(daily_plot_data_arg = daily_plot_data,
                                                                  df_list_arg = df_list_arg,
                                                                  site_arg = site_arg,
                                                                  parameter_arg = parameter_arg) # need to make sure that these functions can read in information from outside of them

          # print ggarrange daily and weekly plots
          print(ggarrange(daily_plot_object, weekly_plot_object, nrow = 2, ncol = 1))

          # Generate the prompt for the user to answer
          day_dt <- as.character(head(daily_plot_data$DT_round, 1)) # reformat so that this reads as a date
          day_choice_prompt <- paste("Would you like to (pass/fail/skip/inspect) ", day_dt, "?\n(layer~", layer, ")(px/fx/sx/ix): ")
          daily_verification_decision <- get_verification_decision(day_choice_prompt)

          while (layer == "daily" & daily_verification_decision == "QUIT") {
            cat("You cannot quit during this step. Please choose a different option.\n")
            daily_verification_decision <- get_verification_decision(day_choice_prompt)
          }

          if (daily_verification_decision %in% simple_decision_list){
            altered_df_list[[i]] <- alter_verification_column(daily_verification_decision, daily_plot_data) # this should get put in a list
          }

          if (daily_verification_decision %in% inspect_decision_list){
            altered_df_list[[i]] <- verify_subdaily(daily_verification_decision_arg = daily_verification_decision,
                                                    daily_plot_data_arg = daily_plot_data,
                                                    daily_plot_object_arg = daily_plot_object,
                                                    weekly_plot_object_arg = weekly_plot_object
            )
          }
        }
      }

      # inspect some weekly ----
      if (verification_decision == "INSPECT SOME") {

        # get which days to inspect
        days_to_inspect <- get_weekly_inspection_decision(weekly_plot_data = df_weekly_data) # rename to daily ***

        new_weekday_list <- weekday_list[weekday_list >= min(df_weekly_data$weekday) & weekday_list <= max(df_weekly_data$weekday)]
        non_days_to_inspect <- new_weekday_list[!new_weekday_list %in% days_to_inspect]

        # for each non-inspect day decision:
        if (length(non_days_to_inspect)>0) {

          # get non-inspect day decision:
          noninspect_prompt_text = paste("Would you like to (pass/fail/skip) non-inspect days?\n(layer~", layer, ")(px/fx/sx): ")
          non_inspect_verification_decision <- get_verification_decision(prompt_text = noninspect_prompt_text)

          # Handle the case when the user selects an "INSPECT" or "QUIT" option for non-inspect days
          while (layer == "daily" & non_inspect_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME", "QUIT")) {
            if(non_inspect_verification_decision %in% c("INSPECT ALL", "INSPECT FLAGGED", "INSPECT VALID", "INSPECT SOME")) {
              cat("You cannot inspect non-inspect days. Please choose a different option.\n")
              non_inspect_verification_decision <- get_verification_decision(prompt_text = noninspect_prompt_text)
            }
            if(non_inspect_verification_decision == "QUIT") {
              cat("You cannot quit during this step. Please choose a different option.\n")
              non_inspect_verification_decision <- get_verification_decision(prompt_text = noninspect_prompt_text)
            }
          }

          for (i in non_days_to_inspect) {
            daily_plot_data <- site_param_df %>%
              filter((y_w %in% df_weekly_data$y_w) & (weekday == i))

            altered_df_list[[i]] <- alter_verification_column(non_inspect_verification_decision, daily_plot_data)
          }
        }

        for (i in days_to_inspect) {

          daily_plot_data <- site_param_df %>%
            filter((y_w %in% df_weekly_data$y_w) & (weekday == i))
          daily_plot_object <- generate_daily_plot(plot_data_arg = daily_plot_data,
                                                   df_list_arg = df_list_arg,
                                                   site_arg = site_arg,
                                                   parameter_arg = parameter_arg)

          # weekly_plot_object

          weekly_plot_object <- generate_supplemental_weekly_plot(daily_plot_data_arg = daily_plot_data,
                                                                  df_list_arg = df_list_arg,
                                                                  site_arg = site_arg,
                                                                  parameter_arg = parameter_arg)

          # print ggarrange daily and weekly plots
          print(ggarrange(daily_plot_object, weekly_plot_object, nrow = 2, ncol = 1))

          # Generate the prompt for the user to answer
          day_dt <- as.character(head(daily_plot_data$DT_round, 1)) # reformat so that this reads as a date
          day_choice_prompt <- paste("Would you like to (pass/fail/skip/inspect) ", day_dt, "?\n(layer~", layer, ")(px/fx/sx/ix): ")
          daily_verification_decision <- get_verification_decision(day_choice_prompt)

          while (layer == "daily" & daily_verification_decision == "QUIT") {
            cat("You cannot quit during this step. Please choose a different option.\n")
            daily_verification_decision <- get_verification_decision(day_choice_prompt)
          }

          if (daily_verification_decision %in% simple_decision_list){
            #
            altered_df_list[[i]] <- alter_verification_column(daily_verification_decision, daily_plot_data) # this should get put in a list
          }

          if (daily_verification_decision %in% inspect_decision_list){
            altered_df_list[[i]] <- verify_subdaily(daily_verification_decision_arg = daily_verification_decision,
                                                    daily_plot_data_arg = daily_plot_data,
                                                    daily_plot_object_arg = daily_plot_object,
                                                    weekly_plot_object_arg = weekly_plot_object
            )
          }
        }
      }

    }


    if (is_tibble(altered_df_list)|is.data.frame(altered_df_list)) {
      altered_df_list <- altered_df_list %>%
        mutate(mean_verified = if_else(grepl(non_removable_flags, flag, ignore.case = TRUE), NA, mean_verified),
               is_verified = TRUE,
               verification_status = ifelse(grepl(non_removable_flags, flag, ignore.case = TRUE), "PASS", verification_status))

    } else if (is.list(altered_df_list)) {
      if (!length(altered_df_list) %in% c(51,49,46)){
        #
        altered_df_list <- map(altered_df_list, function(df) {
          if(!is.null(df)) {
            #
            df %>%
              mutate(mean_verified = if_else(grepl(non_removable_flags, flag, ignore.case = TRUE), NA, mean_verified),
                     is_verified = TRUE,
                     verification_status = ifelse(grepl(non_removable_flags, flag, ignore.case = TRUE), "PASS", verification_status))
          }
        })
      } else {
        #
        altered_df_list <- altered_df_list %>%
          mutate(mean_verified = if_else(grepl(non_removable_flags, flag, ignore.case = TRUE), NA, mean_verified),
                 is_verified = TRUE,
                 verification_status = ifelse(grepl(non_removable_flags, flag, ignore.case = TRUE), "PASS", verification_status))
      }

    }

    return(altered_df_list) # this should not be an altered df, but rather a list of dfs that have been altered.
    # this means that the output of this method will be a list that contains a list of lists

    # If decision in quit (this can only happen before they start verifying a week of data)
    # save the data into the intermediary directory
  } else if (QUIT) {
    df_weekly_name <- ggplot_build(weekly_plot_object)$plot$labels$title
    cat("quitting out of", df_weekly_name, "\n")
    return(NULL)
  }

}

