generate_initial_weekly_plots <- function(all_df_list, pending_df_list, site_arg, parameter_arg, flag_arg = NULL) {

  site_param <- paste0(site_arg, "-", parameter_arg)

  site_flag_dates <- pending_df_list[[site_param]] %>%
    group_by(y_w) %>%
    filter(any(verification_status == "SKIP") | any(!is_verified)) %>%
    ungroup()

  if (!is.null(site_flag_dates)){

    # establish which sites are relevant to the site of interest
    # the order of the sites are hardcoded in relevant_sonde_selector function
    plot_filter <- relevant_sonde_selector(site_arg = site_arg)

    if (nrow(site_flag_dates >0)) {

      if (is.null(flag_arg)) {
        # This for loop generates an overlayed plot of weekly data for the site of
        # interest sandwiched by the site above and below it for each day that was
        # tagged by a flag of interest
        plot_list <- list()

        grouped_data <- site_flag_dates %>%
          group_by(y_w) %>%
          group_split()

        for(i in 1:length(grouped_data)) {

          group_data <- grouped_data[[i]]

          year_week <- unique(group_data$y_w)

          # filtering dfs of interest for the week of interest
          site_df <- site_flag_dates %>%
            filter(y_w == year_week)

          # Get the relevant sonde data
          relevant_sondes <- map(plot_filter, ~ {

            sonde_name <- paste0(.x, "-", parameter_arg)
            data_source <- NULL
            sonde_df <- NULL

            # Determine which directory to pull data from
            tryCatch({
              data_source <- retrieve_relevant_data_name(sonde_name, year_week)
              # cat("Data for",sonde_name,"will be pulled from",data_source,"\n")
            }, error = function(err) {
              # cat("Data for",sonde_name,"not found.\n")
              return(NULL)  # Return NULL if data source can't be determined
            })

            # Only try to pull in the data if data_source was successfully determined
            if (!is.null(data_source)) {
              tryCatch({
                sonde_df <- get(data_source)[[sonde_name]] %>%
                  filter(y_w == group_data$y_w)

              }, error = function(err) {
                cat("Sonde", sonde_name, "not found.\n")
                return(NULL)  # Return NULL if sonde data can't be retrieved
              })
            }

            # Only return a list if both data_source and sonde_df are available
            if (!is.null(data_source) & !is.null(sonde_df)) {
              return(list(sonde_df = tibble(sonde_df), data_source = data_source))
            } else {
              return(NULL)  # Return NULL if either data_source or sonde_df is NULL
            }

          })

          # Remove any NULL results from the list
          relevant_sondes <- compact(relevant_sondes)

          # append site_df to relevant sonde list, clean list, and bind dfs
          # to find plot info
          relevant_dfs <- map(relevant_sondes, ~.x[[1]])
          week_plot_data <- append(relevant_dfs, list(site_df)) %>%
            keep(~ !is.null(.)) %>%
            keep(~ nrow(.)>0) %>%
            bind_rows() %>%
            arrange(day)

          # Create a sequence of dates for the vertical lines
          start_date <- floor_date(min(week_plot_data$DT_round), "day")
          end_date <- ceiling_date(max(week_plot_data$DT_round), "day")
          vline_dates <- seq(start_date, end_date, by = "day")

          date_seq <- seq(min(vline_dates), max(vline_dates) - days(1), by = "1 day") + hours(12)

          # Use the first day of the group as flag_day
          flag_day <- min(group_data$DT_round)

          week_plot <- ggplot(data = week_plot_data) +
            geom_point(data = filter(week_plot_data, (site == site_arg)),
                       aes(x=DT_round, y=mean, color=flag)) +
            map(relevant_sondes, function(sonde_data) {
              data <- sonde_data[[1]]
              data_source <- sonde_data[[2]]

              y_column <- if (data_source == "all_data") "mean" else "mean_verified"

              geom_line(data = data, aes(x = DT_round, y = .data[[y_column]], color = site))
            }) +
            geom_vline(xintercept = vline_dates, color = "black") +
            ggtitle(paste0(str_to_title(site_arg), " ", parameter_arg, " (", format(flag_day, "%B %d, %Y"), ")")) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            labs(x = "Day",
                 y = "Mean")

          week_plot <- add_threshold_lines(plot = week_plot,
                                           plot_data = week_plot_data,
                                           site_arg = site_arg,
                                           parameter_arg = parameter_arg)

          week_plot <- week_plot +
            theme_bw() +
            scale_x_datetime(date_breaks = "1 day",
                             date_labels = "%b %d",
                             minor_breaks = date_seq,
                             sec.axis = sec_axis(~., breaks = date_seq, labels = unique(week_plot_data$weekday))) +
            theme(legend.position = 'bottom',
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()) +
            guides(color = guide_legend(nrow = 4, byrow = TRUE))

          plot_list[[paste(site_param, as.character(flag_day))]] <- week_plot
          sorted_plot_names <- names(plot_list)[order(names(plot_list))]

          plot_list <- plot_list[sorted_plot_names]

        }
      }
      # ---- if flag arg is not null ----
      if (!is.null(flag_arg)) {
        print("under construction...")
        # # This for loop generates an overlayed plot of weekly data for the site of
        # # interest sandwiched by the site above and below it for each day that was
        # # tagged by a flag of interest
        # plot_list <- list()
        #
        # site_flag_weeks <- site_flag_dates %>%
        #   filter(str_detect(flag, flag_arg)) %>%
        #   group_by(week, year) %>%
        #   slice(1)
        #
        # grouped_data <- site_flag_dates %>%
        #   filter(y_w %in% site_flag_weeks$y_w) %>%
        #   group_by(week, year) %>% # group_by(week, month, year) %>%
        #   group_split()
        #
        #
        # for(i in 1:length(grouped_data)) {
        #
        #   group_data <- grouped_data[[i]]
        #
        #   # flag_title <- site_flag_dates$flag[i] # no flag title ***
        #
        #   # filtering dfs of interest for the weeks where a flag was detected
        #   site_df <- site_flag_dates %>%
        #     filter(y_w == group_data$y_w)
        #
        #   # TryCatch used here to avoid erroring out on the first and last values of
        #   # sites_order object (there is no prior/next record after the first/last record).
        #   # Return df as NULL in case of an error
        #   prev_site_df <- NULL
        #   next_site_df <- NULL
        #
        #   tryCatch({
        #     previous_site <- paste0(sites_order[site_index-1],"-",parameter_arg)
        #     prev_site_df <- all_df_list[[previous_site]] %>%
        #       filter(y_w == group_data$y_w)},
        #     error = function(err) {
        #       cat("No previous site.\n")})
        #
        #   tryCatch({
        #     next_site <- paste0(sites_order[site_index+1],"-",parameter_arg)
        #     next_site_df <- all_df_list[[next_site]] %>%
        #       filter(y_w == group_data$y_w)},
        #     error = function(err) {
        #       cat("No next site.\n")})
        #
        #   # Bind all three dfs
        #   week_plot_data <- list(site_df, prev_site_df, next_site_df) %>%
        #     # remove NULL values from the list
        #     keep(~ !is.null(.)) %>%
        #     bind_rows()
        #
        #   # Create a sequence of dates for the vertical lines
        #   start_date <- floor_date(min(week_plot_data$DT_round), "day")
        #   end_date <- ceiling_date(max(week_plot_data$DT_round), "day")
        #   vline_dates <- seq(start_date, end_date, by = "day")
        #
        #   date_seq <- seq(min(vline_dates), max(vline_dates) - days(1), by = "1 day") + hours(12) ## here ----
        #
        #   # Use the first day of the group as flag_day
        #   flag_day <- min(group_data$DT_round)
        #
        #   week_plot <- ggplot(data = week_plot_data) +
        #     geom_point(data = filter(week_plot_data, (site == site_arg)),
        #                aes(x=DT_round, y=mean, color=flag)) +
        #     geom_line(data = filter(week_plot_data, (site != site_arg)),
        #               aes(x=DT_round, y=mean, color=site)) +
        #     geom_vline(xintercept = vline_dates, color = "black") +
        #     ggtitle(paste0(str_to_title(site_arg), " ", parameter_arg, " (", format(flag_day, "%B %d, %Y"), ")")) +
        #     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        #     labs(x = "Day",
        #          y = "Mean")
        #
        #   week_plot <- add_threshold_lines(plot = week_plot,
        #                                    plot_data = week_plot_data,
        #                                    site_arg = site_arg,
        #                                    parameter_arg = parameter_arg)
        #
        #   week_plot <- week_plot +
        #     theme_bw() +
        #     scale_x_datetime(date_breaks = "1 day",
        #                      date_labels = "%b %d",
        #                      minor_breaks = date_seq) +
        #     theme(legend.position = 'bottom',
        #           panel.grid.major = element_blank(),
        #           panel.grid.minor = element_blank()) +
        #     annotate("text", x = date_seq, y = min(week_plot_data$mean, na.rm = TRUE) - 1, label = 1:length(date_seq), hjust = 0) +
        #     guides(color = guide_legend(nrow = 4, byrow = TRUE))
        #
        #   plot_list[[paste(site_param, as.character(flag_day))]] <- week_plot
        #   sorted_plot_names <- names(plot_list)[order(names(plot_list))]
        #
        #   plot_list <- plot_list[sorted_plot_names]
        # }
      }
      return(plot_list)
    } else {
      return(paste(flag_arg, "not detected.\n"))
    }
  } else {
    return(paste(site_arg, parameter_arg, "combination not available.\n"))
  }

}
