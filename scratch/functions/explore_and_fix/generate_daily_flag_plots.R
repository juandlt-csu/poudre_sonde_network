# # Generate daily plots of flagged data.
# # This function will generate a list of daily plots for site-parameters that have been tagged by a specific flag.
# # @param site_arg A string of the site name.
# # @param parameter_arg A string of the parameter name.
# # @param flag_arg A string of the flag name.
# # @param df_list A list of data frames that have been flagged.
# # @return A list of daily plots for site-parameters that have been tagged by a specific flag.
# # @examples
# # generate_daily_flag_plots(site_arg = "archery", parameter_arg = "Actual Conductivity", flag_arg = "outside of Actual Conductivity sensor specification range", df_list = df_list)
# # generate_daily_flag_plots(site_arg = "boxelder", parameter_arg = "Temperature", flag_arg = "outside of Temperature sensor specification range", df_list = df_list)

generate_daily_flag_plots <- function(site_arg, parameter_arg, flag_arg = NULL, df_list) {

  # Generating df name to pull from df_list list
  site_param <- paste0(site_arg, "-", parameter_arg)

  # filter for all the days that are tagged within the site-param df of interest
  site_flag_dates <- df_list[[site_param]]

  if (!is.null(site_flag_dates)){

    site_flag_dates <- site_flag_dates %>%
      filter(if (!is.null(flag_arg)) str_detect(flag, flag_arg) else !is.na(flag)) %>%
      group_by(day(DT_join), month, year) %>%
      slice(1)

    if (nrow(site_flag_dates > 0)) {
      # for loop to generate plots for everyday that was tagged by a flag
      plot_list <- list()

      for (i in 1:nrow(site_flag_dates)) {

        # do this for dates instead : date = format(DT_round, "%m-%d-%Y")
        flag_title <- site_flag_dates$flag[i]
        flag_year <- site_flag_dates$year[i]
        flag_month <- site_flag_dates$month[i]
        flag_day <- site_flag_dates$DT_round[i]

        plot_data <- df_list[[site_param]] %>%
          filter(year == flag_year,
                 month == flag_month,
                 day(DT_round) == day(flag_day))

        y_min <- site_flag_dates$m_mean01[i]
        y_max <- site_flag_dates$m_mean99[i]

        plot <- ggplot(data = plot_data) +
          geom_point(aes(x=DT_round, y = mean, color = flag)) +
          # exceeding sd visualized
          geom_line(aes(x = DT_round, y = rollavg, color = "mean"), show.legend = TRUE) +
          geom_ribbon(aes(ymin = rollavg - m_sd_0199, ymax = rollavg + m_sd_0199, x = DT_round), alpha = 0.1, color = NA) +
          geom_ribbon(aes(ymin = rollavg - (m_sd_0199*2), ymax = rollavg + (m_sd_0199*2), x = DT_round), alpha = 0.1, color = NA) +
          geom_ribbon(aes(ymin = rollavg - (m_sd_0199*3), ymax = rollavg + (m_sd_0199*3), x = DT_round), alpha = 0.1, color = NA) +
          # exceeding range limits visualized
          # geom_vline(data = (plot_data %>% filter(is.na(mean))), aes(xintercept = DT_round, color = flag)) +
          geom_hline(yintercept = y_min, color = "red") +
          geom_hline(yintercept = y_max, color = "red") +
          theme_bw() +
          theme(legend.position = 'bottom') +
          ggtitle(paste(site_arg, parameter_arg, flag_arg, "on", as.character(flag_day))) +
          labs(x = "Datetime",
               y = "Mean")

        plot_list[[paste(site_param, as.character(flag_day))]] <- plot #c(y_min, y_max)

        sorted_plot_names <- names(plot_list)[order(names(plot_list))]

        plot_list <- plot_list[sorted_plot_names]

      }
      return(plot_list)
    } else {
      return(paste(flag_arg, "not detected"))
    }
  } else {
    return(paste(site_arg, parameter_arg, "combination not available"))
  }
}
