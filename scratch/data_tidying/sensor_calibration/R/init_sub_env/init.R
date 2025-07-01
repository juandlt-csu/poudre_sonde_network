# Only edit those things under "====" drop downs

# PACKAGE SET UP () ====

## CRAN PACKAGES (EDIT HERE) ====
### format for name is: "{package_name}"
CRAN_packages <- c(
  "tidyverse", # Data manipulation
  "janitor", # Clean dirty data
  "lubridate", # Date-Time Manipulation
  "rvest", # HTML Retrieval and Manipulation
  "readxl", # Reading excel files
  "here", # Easy, cross platform file referencing
  "ggplot2", # Plotting libraries
  "ggpubr",
  "plotly",
  "devtools", # For downloading GitHub packages
  "remotes",
  "yaml"
)

## GITHUB PACKAGES (EDIT HERE) ====
### format for name is: "{owner/repo}", where the repo is the package name
github_packages <- c(
  "shaughnessyar/driftR", # Drift correction
  "rossyndicate/fcw.qaqc"
)

# PACKAGE LOADER FUNCTIONS ----

## Function to load CRAN packages into environment ----
CRAN_package_loader <- function(x) {
  if (x %in% installed.packages()) {
    library(x, character.only = TRUE, quietly = TRUE)
  } else {
    install.packages(x)
    library(x, character.only = TRUE, quietly = TRUE)
  }
}

## Function to load GitHub packages into environment ----
# If DriftR package is failing, try to download XQuartz to see if that solves the problem
github_package_loader <- function(x) {
  pkg_name <- basename(x)
  if (pkg_name %in% installed.packages()) {
    library(pkg_name, character.only = TRUE, quietly = TRUE)
  } else {
    remotes::install_github(x)
    library(pkg_name, character.only = TRUE, quietly = TRUE)
  }
}

# LOADING PACKAGES ----

message("Loading CRAN packages into environment...")
for (package in CRAN_packages) {
  invisible(suppressMessages(CRAN_package_loader(package)))
  cat("Loaded:", package, "\n")
}

message("Loading GitHub packages into environment...")
for (package in github_packages) {
  invisible(suppressMessages(github_package_loader(package)))
  cat("Loaded:", package, "\n")
}

# SUB-ENVIRONMENT SET UP (EDIT HERE) ====

## EDIT this {here} path to the directory you will be working in (EDIT HERE) ====
wd <- here("scratch", "data_tidying", "sensor_calibration")

## EDIT this {file.path} to where your custom functions are relative to your wd (EDIT HERE) ====
custom_functions_path <- file.path("R")

# ESTABLISHING SUB-ENVIRONMENT function ----
sub_here <- function(...){
  sub_path <- here(wd, ...)
  return(sub_path)
}

# LOADING CUSTOM CALIBRATION FUNCTIONS INTO THE ENVIRONMENT ----
message("Loading custom calibration functions into the environment...")
walk(list.files(sub_here(custom_functions_path), pattern = "*.R", full.names = TRUE), source)
cat("Loaded:", length(list.files(sub_here(custom_functions_path), pattern = "*.R")), "custom functions")
