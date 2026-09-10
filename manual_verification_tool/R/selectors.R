#copied over from `source`, select site for sub plots
relevant_sonde_selector <- function(site_arg) {
  if (site_arg == "joei") {
    plot_filter <- c("cbri")
  }
  if (site_arg == "cbri") {
    plot_filter <- c("joei", "chd")
  }
  if (site_arg == "chd") {
    plot_filter <- c("cbri", "pfal")
  }
  if (site_arg == "pfal") {
    plot_filter <- c("chd", "pbr")
  }
  if (site_arg == "pbr") {
    plot_filter <- c("pfal", "pman")
  }
  if (site_arg == "sfm") {
    plot_filter <- c("pfal", "pbr")
  }
  if (site_arg == "pman") {
    plot_filter <- c("pbr", "pbd")
  }
  if (site_arg == "penn") {
    plot_filter <- c("sfm")
  }
  if (site_arg == "lbea") {
    plot_filter <- c("sfm")
  }
  if (site_arg == "pbd") {
    plot_filter <- c("pman", "bellvue")
  }
  if (site_arg == "bellvue") {
    plot_filter <- c("pbd", "salyer")
  }
  if (site_arg == "salyer") {
    plot_filter <- c("bellvue", "riverbend")
  }
  if (site_arg == "udall") {
    plot_filter <- c("salyer", "riverbend")
  }
  if (site_arg == "riverbend") {
    plot_filter <- c("salyer", "cottonwood")
  }
  if (site_arg == "riverbend_virridy") {
    plot_filter <- c("udall", "riverbend", "cottonwood")
  }
  if (site_arg == "springcreek") {
    plot_filter <- c("riverbend", "cottonwood")
  }
  if (site_arg == "cottonwood") {
    plot_filter <- c("riverbend", "cottonwood virridy", "elc")
  }
  if (site_arg == "cottonwood_virridy") {
    plot_filter <- c( "cottonwood", "elc")
  }
  if (site_arg == "elc") {
    plot_filter <- c("cottonwood",
                     "archery")
  }
  if (site_arg == "boxcreek") {
    plot_filter <- c("archery", "archery_virridy")
  }
  if (site_arg == "archery") {
    plot_filter <- c("elc", "archery_virridy", "riverbluffs")
  }
  if (site_arg == "archery_virridy") {
    plot_filter <- c("elc", "archery", "riverbluffs")
  }
  if (site_arg == "riverbluffs") {
    plot_filter <- c("archery", "archery_virridy")
  }
  if (site_arg == "mtncampus") {
    plot_filter <- c("sfm")
  }
  return(plot_filter)
}



# # Get parameters available for a site
# get_parameters <- function(datasets, directory, site) {
#   if (directory == "pre_verification") {
#     data_list <- datasets$pre_verification_data
#   } else {
#     data_list <- datasets$intermediary_data
#   }
#
#   names(data_list) %>%
#     keep(str_detect(., paste0("^", site, "-"))) %>%
#     str_remove(paste0(site, "-"))
# }

get_auto_parameters <- function(parameter) {
  tryCatch({
     read_csv(file.path(meta_path, "parameter_autoselections.csv"), show_col_types = F) %>%
      filter(main_parameter == parameter) %>%
      pull(sub_parameters) %>%
      first() %>%
      str_split(",", simplify = TRUE) %>%
      as.character() %>%
      str_trim() %>%
      .[. != ""] # Remove any empty strings
  }, error = function(e) character(0))
}


retrieve_relevant_data_name <- function(df_name_arg, year_week_arg = NULL) {

  if (df_name_arg %in% names(verified_data) & any(year_week_arg %in% verified_data[[df_name_arg]]$y_w)) {
    return("verified_data")
  }
  if (df_name_arg %in% names(intermediary_data) & any(year_week_arg %in% intermediary_data[[df_name_arg]]$y_w)) {
    return("intermediary_data")
  }
  if (df_name_arg %in% names(pre_verification_data) & any(year_week_arg %in% pre_verification_data[[df_name_arg]]$y_w)) {
    return("pre_verification_data")
  }
  if (df_name_arg %in% names(all_data) & any(year_week_arg %in% all_data[[df_name_arg]]$y_w)) {
    return("all_data")
  }

}
