#' Unfinalize a site-parameter dataset
#' This function takes a finalized site-parameter dataset from the verified directory, sets the `is_finalized` column to FALSE, and saves it to the intermediary directory with a new filename that removes the `_FINAL_` suffix.
#' @param site The site name (e.g., "archery")
#' @param parameter The parameter name (e.g., "Temperature")
#' Parameters Available: "DO", "Temperature", "pH","ORP", "Specific Conductivity", "Turbidity", "FDOM Fluorescence", "Chl-a Fluorescence", "Depth
#' @param verified_directory The path to the verified directory (default: "shiny_ver
#' tool/data/verified_directory")
#' @param intermediary_directory The path to the intermediary directory (default: "shiny_ver
#' tool/data/intermediary_directory")
#' @return The site-parameter dataset with `is_finalized` set to FALSE, saved to the intermediary directory. The function returns the path to the new file invisibly.
#'
#'@examples
#' \dontrun{
#' unfinalize_site_parameter(site = "archery", parameter = "Temperature")
#'

unfinalize_site_parameter <- function(site,
                                      parameter,
                                      verified_directory = here(verified_path),
                                      intermediary_directory = here(intermediary_path)) {

  # Find matching file in verified directory
  pattern <- paste0("^", site, "-", parameter, "_FINAL_")
  files <- list.files(verified_directory, pattern = pattern, full.names = TRUE)

  if (length(files) == 0) {
    stop("No finalized file found for site: ", site, " parameter: ", parameter)
  }
  if (length(files) > 1) {
    stop("Multiple finalized files found, using most recent: ", basename(files[length(files)]))
  }

  # Read file
  df <- read_parquet(files)

  # Set is_finalized to FALSE
  df <- df %>% mutate(is_finalized = FALSE)

  # Build new filename without _FINAL_
  new_filename <- basename(files) %>% str_replace("_FINAL_", "_")
  output_path <- file.path(intermediary_directory, new_filename)

  # Save to intermediary directory
  write_parquet(df, output_path)
  file.remove(files)

  message("Unfinalized file saved to: ", output_path)

  invisible(output_path)
}
