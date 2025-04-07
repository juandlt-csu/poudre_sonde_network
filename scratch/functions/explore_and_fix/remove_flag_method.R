# # Verify flag data ----
# ## This uses the data that is generated from the function generate_daily_flag_plots()
# ## to ask the user whether or not a flagged data point should pass or fail.
#
# # sourcing the functions
# sapply(list.files("src/qaqc/verification_functions", pattern = "*.R", full.names = TRUE), source, .GlobalEnv)
#
# # Verify flag data function and testing ----
# ## Right now this function takes in a ggplot object, extracts its data and then
# ## makes a plot for each point of interest to ask the user if it is valid point
# ## or not. As it is how the user answers has no effect on any generated data.
# ## Instead, new dfs will be generated with the users decisions.
#
# # Examples ----
#
# # Right now this takes in a list of plots, I don't love this but for now its the
# # easiest way to quickly get to the data that has flagged data.
#
# # Generating the list of plots (this example only has one plot that has a repeated value flag)
# # we are doing it this way so that we can get the data from the daily ggplot
# test_weekly <- generate_weekly_flag_plots("archery", "Temperature", "slope violation")
#
# test_daily <- generate_daily_flag_plots("archery", "Temperature", "slope violation")
#
# # Using the list of plots that we generated to verify just one of the dfs in the list
#
# # maybe we will use map2() so that it can take the input from both and then combine them?
# verify_flag_data(test_weekly[[1]], test_daily[[1]]) %>% View()
#
# # Using the function that was made to map over the dfs in the list
# # make sure to put in the weekly information first
# # test_map <- map2(head(test_weekly), head(test_daily), verify_flag_data)
#
# # Need to make a method to track what we have done
