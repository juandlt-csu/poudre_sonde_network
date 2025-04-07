# # Generate weekly plots of flagged data.
# # This function will generate a list of weekly plots for site-parameters that have been tagged by a specific flag.
# # @param site_arg A string of the site name.
# # @param parameter_arg A string of the parameter name.
# # @param flag_arg A string of the flag name.
# # @param df_list A list of dataframes that have been flagged.
# # @return A list of weekly plots for site-parameters that have been tagged by a specific flag.
# # @examples
# # generate_weekly_flag_plots(site_arg = "archery", parameter_arg = "Actual Conductivity", flag_arg = "outside of Actual Conductivity sensor specification range", df_list = df_list)
# # generate_weekly_flag_plots(site_arg = "boxelder", parameter_arg = "Temperature", flag_arg = "outside of Temperature sensor specification range", df_list = df_list)
#
# generate_weekly_flag_plots <- function(site_arg, parameter_arg, flag_arg = NULL, df_list) {
#
#   site_param <- paste0(site_arg, "-", parameter_arg)
#
#   site_flag_dates <- df_list[[site_param]]
#
#   if (!is.null(site_flag_dates)){
#     # vector of sites in the order that they are in spatially
#     # some sites have some funkiness going on
#     sites_order <- c("tamasag", # rist
#                      "legacy",
#                      "lincoln",
#                      "timberline",
#                      "prospect",
#                      "boxelder", # elc
#                      "archery",
#                      "river bluffs")
#     # determining the index for the site of interest.
#     site_index <- which(sites_order == site_arg)
#
#     # Generating df name to pull from df_list list
#     site_param <- paste0(site_arg, "-", parameter_arg)
#     # filter for all the days that are tagged within the site-param df of interest
#     site_flag_dates <- site_flag_dates %>%
#       filter(if (is.null(flag_arg)) !is.na(flag) else str_detect(flag, flag_arg)) %>%
#       group_by(day(DT_join), month, year) %>%
#       slice(1)
#
#     if (nrow(site_flag_dates > 0)) {
#       # This for loop generates an overlayed plot of weekly data for the site of
#       # interest sandwiched by the site above and below it for each day that was
#       # tagged by a flag of interest
#       plot_list <- list()
#
#       for(i in 1:nrow(site_flag_dates)) {
#
#         flag_title <- site_flag_dates$flag[i]
#         flag_year <- site_flag_dates$year[i]
#         flag_month <- site_flag_dates$month[i]
#         flag_day <- site_flag_dates$DT_round[i]
#         # Getting the prior and subsequent 3 days to the flag day
#         start_day <- flag_day - days(3)
#         end_day <- flag_day + days(3)
#
#         # filtering dfs of interest for the weeks where a flag was detected
#         site_df <- df_list[[site_param]] %>%
#           filter(year == flag_year,
#                  month == flag_month,
#                  DT_round >= start_day & DT_round <= end_day)
#
#         # TryCatch used here to avoid erroring out on the first and last values of
#         # sites_order object (there is no prior/next record after the first/last record).
#         # Return df as NULL in case of an error
#         prev_site_df <- NULL
#         next_site_df <- NULL
#
#         tryCatch({
#           previous_site <- paste0(sites_order[site_index-1],"-",parameter_arg)
#           prev_site_df <- df_list[[previous_site]] %>%
#             filter(year == flag_year,
#                    month == flag_month,
#                    DT_round >= start_day & DT_round <= end_day)},
#           error = function(err) {
#             cat("No previous site.\n")})
#
#         tryCatch({
#           next_site <- paste0(sites_order[site_index+1],"-",parameter_arg)
#           next_site_df <- df_list[[next_site]] %>%
#             filter(year == flag_year,
#                    month == flag_month,
#                    DT_round >= start_day & DT_round <= end_day)},
#           error = function(err) {
#             cat("No next site.\n")})
#
#         # Bind all three dfs
#         week_plot_data <- list(site_df, prev_site_df, next_site_df) %>%
#           # remove NULL values from the list
#           keep(~ !is.null(.)) %>%
#           bind_rows()
#
#         # Using the data from the day where a flag was detected to generate a window
#         # to easily distinguish the data of interest in comparison with the rest of
#         # the data
#         site_day_data <- df_list[[site_param]] %>%
#           filter(year == flag_year,
#                  month == flag_month,
#                  day(DT_round) == day(flag_day))
#
#         y_min <- site_day_data$m_mean01[i]
#         y_max <- site_day_data$m_mean99[i]
#         # there is a join somewhere that adds two site and parameter cols
#         # need to find it and fix that
#         week_plot <- ggplot(data = week_plot_data, aes(x=DT_round, y=mean, color=site)) +
#           geom_rect(data = site_day_data, aes(xmin = min(DT_round), xmax = max(DT_round),
#                                               ymin = -Inf, ymax = Inf),
#                     fill = "grey",
#                     alpha = 0.01,
#                     color = NA)+
#           geom_line() +
#           geom_hline(yintercept = y_min, color = "red") +
#           geom_hline(yintercept = y_max, color = "red") +
#           theme_bw() +
#           theme(legend.position = 'bottom') +
#           ggtitle(paste(flag_arg,"at", site_arg, "on", as.character(flag_day))) +
#           labs(x = "Datetime",
#                y = "Mean")
#
#         plot_list[[paste(site_param, as.character(flag_day))]] <- week_plot
#
#         sorted_plot_names <- names(plot_list)[order(names(plot_list))]
#
#         plot_list <- plot_list[sorted_plot_names]
#       }
#       return(plot_list)
#     } else {
#       return(paste(flag_arg, "not detected.\n"))
#     }
#   } else {
#     return(paste(site_arg, parameter_arg, "combination not available.\n"))
#   }
# }
