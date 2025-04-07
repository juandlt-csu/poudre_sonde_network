
# Run this for loop to be able to verify the data
# Type R to see the verification options available to you

for (i in weekly_plot_objects) {

  if(QUIT){
    break()
  }

  update <- verify_flag_data(
    df_list_arg = all_data,
    site_arg = site,
    parameter_arg = parameter,
    flag_arg = NULL,
    weekly_plot_object = i
  )

  if (!is.null(update)){

    if(is_tibble(update)|is.data.frame(update)) {
      update <- bind_rows(update)
    } else {
      if (!length(update) %in% c(51,49,46)) { # 46 is because sometimes dfs decompose into lists... idk why -jd
        update <-  update %>%
          # remove NULL values from the list
          keep(~ !is.null(.)) %>%
          # remove empty dfs from the list
          keep(~ nrow(.)>0) %>%
          bind_rows(.) %>%
          distinct()
      } else {
        update <- bind_rows(update) %>%
          distinct()
      }
    }

    # update site param df
    updated_site_param_df <- updated_site_param_df %>%
      mutate(
        mean_verified = as.double(mean_verified),
        verification_status = as.character(verification_status)) %>%
      rows_update(update, by = "DT_join")

    # update the saved data in the intermediary dir
    saveRDS(object = updated_site_param_df, file = here(intermediary_path, site_param_name))

    graphics.off()
    gc()
  } else {
    next
  }
}

# Quit toggle
QUIT <- FALSE

# check data
updated_df <- read_rds(here(intermediary_path, site_param_name)) %>%
  select(DT_join, site, parameter, mean, mean_verified, is_verified, verification_status, flag) %>%
  filter(is_verified) %>%
  View()

# clean the directories
clean_directories()
