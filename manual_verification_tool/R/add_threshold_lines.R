add_threshold_lines <- function(plot, plot_data, site_arg, parameter_arg) {

  # pull in threshold data
  # real_thresholds <- read_csv(file.path(meta_path, "realistic_thresholds.csv"), show_col_types = FALSE) %>%
  #   filter(parameter == parameter_arg)
  #
  # sensor_thresholds <- yaml::read_yaml(file.path(meta_path, "sensor_spec_thresholds.yml"))[[parameter_arg]] %>%
  #   bind_rows()

  seasonal_thresholds <- read_csv(file.path(meta_path, 'seasonal_thresholds.csv'), show_col_types = FALSE) %>%
    # to do: Check to make sure seasonal thresholds csv is not necessary
    #bind_rows(read_csv(here('data', 'qaqc', 'seasonal_thresholds.csv'), show_col_type = FALSE),
    distinct(site, parameter, season, .keep_all = TRUE) %>%
    filter(parameter == parameter_arg,
           site == site_arg)

  # Determine the unique seasons in plot_data
  unique_seasons <- unique(plot_data$season)

  # Filter seasonal_thresholds for the seasons present in plot_data
  seasonal_thresholds <- seasonal_thresholds %>%
    filter(season %in% unique_seasons)


  if (nrow(seasonal_thresholds) > 1) { # make sure this works

    season_1 <- case_when(
      all(c("winter_baseflow", "snowmelt") %in% unique_seasons) ~ "winter_baseflow",
      all(c("snowmelt", "monsoon") %in% unique_seasons) ~ "snowmelt",
      all(c("monsoon", "fall_baseflow") %in% unique_seasons) ~ "monsoon",
      all(c("fall_baseflow", "winter_baseflow") %in% unique_seasons) ~ "fall_baseflow",
      TRUE ~ NA_character_
    )
    season_2 <- case_when(
      all(c("winter_baseflow", "snowmelt") %in% unique_seasons) ~ "snowmelt",
      all(c("snowmelt", "monsoon") %in% unique_seasons) ~ "monsoon",
      all(c("monsoon", "fall_baseflow") %in% unique_seasons) ~ "fall_baseflow",
      all(c("fall_baseflow", "winter_baseflow") %in% unique_seasons) ~ "winter_baseflow",
      TRUE ~ NA_character_
    )

    seasonal_threshold_s1_quantiles <- unname(quantile(c(seasonal_thresholds %>% filter(season == season_1) %>% pull(t_mean01),
                                                         seasonal_thresholds %>% filter(season == season_1) %>% pull(t_mean99)),
                                                       c(0.1, 0.9)))

    seasonal_threshold_s2_quantiles <- unname(quantile(c(seasonal_thresholds %>% filter(season == season_2) %>% pull(t_mean01),
                                                         seasonal_thresholds %>% filter(season == season_2) %>% pull(t_mean99)),
                                                       c(0.1, 0.9)))

    slice_data <- plot_data %>%
      group_by(season) %>%
      slice(1) %>%
      ungroup()

    site_data <- filter(plot_data, site == site_arg)

    if (!all(is.na(site_data$mean))) {
      # Lower bound
        # Xs

      lower_bounds <- tibble(
        start_x_s1 = min(slice_data$DT_round),

        transition_date = slice_data %>%
          filter(season == season_2) %>%
          arrange(DT_round) %>%
          pull(DT_round),

        end_x_s2  = ceiling_date(max(plot_data$DT_round), "day"),

        #Ys
        y_lower_s1 =  seasonal_thresholds %>%
          filter(season == season_1) %>%
          pull(t_mean01),

        y_lower_s2 =  seasonal_thresholds %>%
          filter(season == season_2) %>%
          pull(t_mean01)
      )
        # #Testing
        # annotate(geom = "segment", x = start_x_s1, xend = transition_date, y = y_lower_s1, yend = y_lower_s1, linetype = 2)+
        #   annotate(geom = "segment", x = transition_date, xend = end_x_s2, y = y_lower_s2, yend = y_lower_s2, linetype = 2)

        plot <- plot +
          # season 1
          ## lower bound
          geom_segment(data = lower_bounds, aes(x = start_x_s1, y = y_lower_s1,
                           xend = transition_date, yend = y_lower_s1,
                           color = "Seasonal Min: 1", linetype = "Seasonal Min 1")) +
          # season 2
          ## lower bound
          geom_segment(data = lower_bounds, aes(x = transition_date, y = y_lower_s2,
                           xend = end_x_s2, yend = y_lower_s2,
                           color = "Seasonal Min: 2", linetype = "Seasonal Min 2"))

      # Upper bound
        # Xs
        upper_bound <- tibble(
        start_x_s1 =  min(slice_data$DT_round),
        transition_date =  slice_data %>%
          filter(season == season_2) %>%
          arrange(DT_round) %>%
          pull(DT_round),
        end_x_s2 = ceiling_date(max(plot_data$DT_round), "day"),
        #Ys
        y_upper_s1 =  seasonal_thresholds %>%
          filter(season == season_1) %>%
          pull(t_mean99),
        y_upper_s2 =  seasonal_thresholds %>%
          filter(season == season_2) %>%
          pull(t_mean99)
        )



        plot <- plot +
          # season 1
          ## upper Bound
          geom_segment(data = upper_bound, aes(x = start_x_s1, y = y_upper_s1,
                           xend = transition_date, yend = y_upper_s1,
                           color = "Seasonal Max: 1", linetype = "Seasonal Max 1")) +
          # season 2
          ## Upper Bound
          geom_segment(data = upper_bound, aes(x = transition_date, y = y_upper_s2,
                           xend = end_x_s2, yend = y_upper_s2,
                           color = "Seasonal Max: 2", linetype = "Seasonal Max 2"))+
          labs(linetype = "Thresholds")
          # plot <- plot +
          # # season 1
          # ## upper Bound
          # annotate(geom = "segment", x = start_x_s1, xend = transition_date, y = y_upper_s1, yend = y_upper_s1, linetype = 1)+
          # annotate(geom = "segment", x = transition_date, xend = end_x_s2, y = y_upper_s2, yend = y_upper_s2, linetype = 1) +
          # scale_linetype_manual(values = c( "Seasonal Max", "Seasonal Min"))+
          # labs(linetype = "Thresholds")

      return(plot)
    }

  } else if (nrow(seasonal_thresholds) == 1){

    site_data <- filter(plot_data, site == site_arg)

    # Filter seasonal_thresholds for the seasons present in site_data
    seasonal_thresholds <- seasonal_thresholds %>%
      filter(season %in% unique(site_data$season))

    if (!all(is.na(site_data$mean))) {
      # seasonal thresholds
      seasonal_thresholds_quantiles <- unname(quantile(c(seasonal_thresholds$t_mean01, seasonal_thresholds$t_mean99), c(0.1, 0.9)))

      #if (min(site_data$mean, na.rm = TRUE) <= seasonal_thresholds_quantiles[1]){
        plot <- plot +
          geom_hline(aes(yintercept = seasonal_thresholds$t_mean01,
                         color = "Seasonal Min",
                         linetype = "Seasonal Min"))
      #}

      #if (max(site_data$mean, na.rm = TRUE) >= seasonal_thresholds_quantiles[2]) {
        plot <- plot +
          geom_hline(aes(yintercept = seasonal_thresholds$t_mean99,
                         color = "Seasonal Max",
                         linetype = "Seasonal Max"))
      #}

      # # real thresholds
      # real_thresholds_quantiles <- unname(quantile(c(real_thresholds$min, real_thresholds$max), c(0.1, 0.9)))
      #
      # if (min(site_data$mean, na.rm = TRUE) <= real_thresholds_quantiles[1]) {
      #   plot <- plot +
      #     geom_hline(aes(yintercept = real_thresholds$min,
      #                    color = "Real Min",
      #                    linetype = "Real"))
      # }
      # if (max(site_data$mean, na.rm = TRUE) >= real_thresholds_quantiles[2]) {
      #   plot <- plot +
      #     geom_hline(aes(yintercept = real_thresholds$min,
      #                    color = "Real Max",
      #                    linetype = "Real"))
      # }
      #
      # # sensor thresholds
      # sensor_thresholds_quantiles <- unname(quantile(c(sensor_thresholds$min, sensor_thresholds$max), c(0.1, 0.9)))
      #
      # # if (min(site_data$mean, na.rm = TRUE) <= sensor_thresholds_quantiles[1]) {
      # #   plot <- plot +
      # #     geom_hline(aes(yintercept = sensor_thresholds$mix,
      # #                color = "Sensor Mine",
      # #                linetype = "Sensor"))
      # # }
      #
      # if (max(site_data$mean, na.rm = TRUE) >= sensor_thresholds_quantiles[2]) {
      #   plot <- plot +
      #     geom_hline(aes(yintercept = sensor_thresholds$max,
      #                    color = "Sensor Max",
      #                    linetype = "Sensor"))
      # }
        plot <- plot+
          labs(linetype = "Thresholds")
      return(plot)
    }
    return(plot)
  }
  return(plot)
}
