# =============================================================================
# 00_config.R
# Central configuration for the monthly budget pipeline.
# Sourced by 01_check_hours.Rmd, 02_prep_charges_report_gen.Rmd.
# Shared functions for the monthly budget pipeline.
# =============================================================================


suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(here)
})

`%nin%` <- Negate(`%in%`)

# -----------------------------------------------------------------------------
# 1. SHARED FILE PATHS  --------------------------------------------------------
#    All budget data (rates + timeclock + equipment + travel) lives in ONE
#    workbook with a `project` column, so both projects read the same file.
# -----------------------------------------------------------------------------
budget_dir      <- here("docs", "pwqn", "docs", "monthly_budget")
rates_file      <- file.path(budget_dir, "employee_hourly_rates.csv")
reporting_xlsx  <- file.path(budget_dir, "pwqn_budget_reporting.xlsx")
mwater_creds_file <- here("creds", "mWaterCreds.yml")

idc_rate <- 0.10  # indirect cost rate, applied to (personnel + equipment + travel)

# -----------------------------------------------------------------------------
# 2. HOURLY RATES  ------------------------------------------------------------
#    The hourly rates for each employee, used to calculate the cost of their time.
# -----------------------------------------------------------------------------
projects <- list(

  PWQN = list(
    code           = "PWQN",
    report_title   = "PWQN",
    caption_suffix = "",
    output_dir     = here("docs", "pwqn", "docs", "monthly_budget"),
    category_order = c("Field", "Equipment", "Travel", "Admin", "Annual Report", "QAQC"),
    sites          = c("bellvue", "salyer", "udall", "riverbend", "cottonwood",
                        "elc", "archery", "riverbluffs"),
    site_labels    = c(bellvue = "Bellvue", salyer = "Salyer", udall = "Udall",
                        riverbend = "Riverbend", cottonwood = "Cottonwood",
                        elc = "ELC", archery = "Archery", riverbluffs = "River Bluffs"),
    # pwqn counts a "site visit" as any field-note row that ISN'T a water-sampling row
    field_note_method = "exclude_sampling"
  ),

  UCLP = list(
    code           = "UCLP",
    report_title   = "PDSS",
    caption_suffix = " to City of Fort Collins IGA",
    output_dir     = here("docs", "uclp_dss", "docs", "monthly_budget"),
    category_order = c("Field", "Equipment", "Travel", "Admin", "Annual Report", "QAQC"),
    sites          = c("chd", "pfal", "pbr", "sfm", "pman", "pbd"),
    site_labels    = NULL,  # falls back to toupper(site)
    # uclp combines sensor notes + water-sampling notes and counts unique visit dates
    field_note_method = "unique_dates"
  )
)

# convenience: pull one project's config by name, with a friendly error
get_project_cfg <- function(project_code) {
  cfg <- projects[[project_code]]
  if (is.null(cfg)) {
    stop("Unknown project '", project_code, "'. Options are: ",
         paste(names(projects), collapse = ", "))
  }
  cfg
}

#Functions for loading and processing data for the monthly budget pipeline

# -----------------------------------------------------------------------------
# Hourly rates ------------------------------------------------------------
# -----------------------------------------------------------------------------
load_hourly_rates <- function(rates_file) {
  read_csv(rates_file, show_col_types = FALSE) %>%
    mutate(
      fringe = case_when(
        type == "Salary" ~ 0.286,
        type == "Hourly" ~ 0.012
      ),
      hourly_n_fringe = case_when(
        type == "Salary" ~ rate / 1230,
        type == "Hourly" ~ rate
      ),
      hourly_billed = case_when(
        type == "Salary" ~ (rate + (rate * fringe)) / 1230,
        type == "Hourly" ~ rate + (rate * fringe)
      )
    )
}

# -----------------------------------------------------------------------------
# Timeclock / equipment / travel loaders -----------------------------------
# Each returns ALL projects unless `project` is supplied, so the same call
# works for both the cross-project reconciliation step (01) and the
# per-project prep/gen step (02).
# -----------------------------------------------------------------------------
load_timeclock <- function(reporting_xlsx, month_select, year_select, project = NULL) {
  out <- read_excel(reporting_xlsx, sheet = "staff_timeclock") %>%
    mutate(month = month(date), year = year(date)) %>%
    filter(month == month_select, year == year_select)
  if (!is.null(project)) out <- filter(out, project == !!project)
  out
}

load_equipment <- function(reporting_xlsx, month_select, year_select, project = NULL) {
  out <- read_excel(reporting_xlsx, sheet = "equipment_charges") %>%
    mutate(month = month(date), year = year(date)) %>%
    filter(month == month_select, year == year_select)
  if (!is.null(project)) out <- filter(out, project == !!project)
  out %>%
    mutate(category = "Equipment",
           include = ifelse(startup == TRUE, "No", "Yes")) %>%
    select(project, category, staff = vendor, description = item,
           cost = amount, include) %>%
    mutate(include = as.character(include))
}

load_travel <- function(reporting_xlsx, month_select, year_select, project = NULL) {
  out <- read_excel(reporting_xlsx, sheet = "travel") %>%
    mutate(month = month(date), year = year(date)) %>%
    filter(month == month_select, year == year_select)
  if (!is.null(project)) out <- filter(out, project == !!project)
  out %>%
    mutate(category = "Travel",
           include = ifelse(startup == TRUE, "No", "Yes")) %>%
    select(project, category, description, cost = amount, include) %>%
    mutate(include = as.character(include))
}

# -----------------------------------------------------------------------------
# Personnel cost per staff x project (used by both 01 and 02) ----------------
# -----------------------------------------------------------------------------
summarize_personnel_hours <- function(timeclock) {
  timeclock %>%
    group_by(category, staff, project, Extension) %>%
    summarize(description = paste(unique(description), collapse = ", "),
              hours = sum(hours), .groups = "drop")
}

cost_personnel <- function(hours_summary, hourly_rates, exclude_categories = character(0)) {
  hours_summary %>%
    left_join(select(hourly_rates, person, hourly_billed), by = c("staff" = "person")) %>%
    mutate(
      hourly_billed = round(hourly_billed, 2),
      hours_rate = paste0(hours, " x $", hourly_billed, "/hr"),
      cost = round(hours * hourly_billed, 2),
      include = case_when(
        category %in% exclude_categories ~ "No",
        Extension == TRUE ~ "No",   # staff/intern time paid for by Extension
        TRUE ~ "Yes"
      )
    ) %>%
    select(-hours, -hourly_billed)
}

# -----------------------------------------------------------------------------
# Reconciliation: hours worked (from $ billed / hourly rate) vs. hours logged
# in the timeclock sheet, PER PERSON ACROSS BOTH PROJECTS (single job code).
# -----------------------------------------------------------------------------
reconcile_hours_and_cost <- function(timeclock, hourly_rates, hours_billed_raw) {

  hours_billed <- hours_billed_raw %>%
    left_join(hourly_rates, by = "person") %>%
    mutate(hours_worked = ceiling(monthly_amount / hourly_billed)) %>%
    select(person, monthly_amount, hourly_billed, hours_worked)

  hours_by_project <- timeclock %>%
    group_by(staff, project) %>%
    summarize(hours = sum(hours), .groups = "drop")

  hours_check <- hours_by_project %>%
    pivot_wider(names_from = project, values_from = hours, values_fill = 0) %>%
    mutate(total_hours = rowSums(across(where(is.numeric)))) %>%
    left_join(select(hours_billed, person, hours_worked), by = c("staff" = "person")) %>%
    mutate(remaining_hours = hours_worked - total_hours)

  cost_check <- hours_by_project %>%
    left_join(select(hourly_rates, person, hourly_billed), by = c("staff" = "person")) %>%
    mutate(hourly_billed = round(hourly_billed, 2),
           cost = hours * hourly_billed) %>%
    select(-hourly_billed, -hours) %>%
    pivot_wider(names_from = project, values_from = cost, values_fill = 0) %>%
    mutate(total_amount = rowSums(across(where(is.numeric)))) %>%
    left_join(select(hours_billed, person, monthly_amount), by = c("staff" = "person")) %>%
    mutate(remaining_amount = total_amount - monthly_amount)

  list(hours_check = hours_check, cost_check = cost_check, hours_billed = hours_billed)
}

# -----------------------------------------------------------------------------
# Final totals: personnel + equipment + travel + IDC, per project ------------
# -----------------------------------------------------------------------------
calc_final_totals <- function(timeclock, hourly_rates, equipment, travel, idc_rate = 0.10) {

  personnel_cost_per_project <- timeclock %>%
    summarize_personnel_hours() %>%
    cost_personnel(hourly_rates) %>%
    group_by(project) %>%
    summarise(personnel_cost = sum(cost), .groups = "drop")

  non_personnel <- bind_rows(
    equipment %>% group_by(project) %>% summarise(equipment_materials_cost = sum(cost), .groups = "drop"),
    travel %>% group_by(project) %>% summarise(travel_cost = sum(cost), .groups = "drop"),
    .id = NULL
  ) %>%
    group_by(project) %>%
    summarise(across(everything(), ~ sum(., na.rm = TRUE)), .groups = "drop")

  personnel_cost_per_project %>%
    full_join(non_personnel, by = "project") %>%
    mutate(across(where(is.numeric), ~ replace_na(., 0))) %>%
    mutate(
      orig_total = personnel_cost + travel_cost + equipment_materials_cost,
      idc_amount = orig_total * idc_rate,
      final_total = orig_total + idc_amount
    ) %>%
    mutate(across(where(is.numeric), ~ round(., 2)))
}

# -----------------------------------------------------------------------------
# Field notes (mWater) -------------------------------------------------------
# Generalizes both projects' logic. `method`:
#   "exclude_sampling" (PWQN): a visit = any non-"Water Sampling" field-note row
#   "unique_dates"     (UCLP): combine sensor + sampling notes, count unique
#                              visit dates per site
# -----------------------------------------------------------------------------
grab_field_notes <- function(mWater_creds, cfg, month_select, year_select, full = FALSE) {

  mWater_data <- load_mWater(creds = mWater_creds)

  sensor_notes <- grab_mWater_sensor_notes(mWater_api_data = mWater_data) %>%
    mutate(date = as.Date(date), month = month(date), year = year(date),
           site = tolower(site)) %>%
    filter(site %in% cfg$sites, month == month_select, year == year_select) %>%
    select(date, crew, site, cals_performed, DT_round, visit_type, visit_comments)

  if (cfg$field_note_method == "unique_dates") {
    sampling_notes <- mWater_data %>%
      filter(visit_type == "Water Sampling") %>%
      mutate(date = as.Date(date), month = month(date), year = year(date),
             site = tolower(site)) %>%
      filter(site %in% cfg$sites, month == month_select, year == year_select) %>%
      select(date, crew, site, cals_performed, DT_round, visit_type, visit_comments)

    field_notes <- bind_rows(sensor_notes, sampling_notes)
    monthly_visits <- field_notes %>%
      group_by(site) %>%
      summarise(visits = length(unique(date)), .groups = "drop")
  } else {
    field_notes <- sensor_notes
    monthly_visits <- field_notes %>%
      filter(visit_type != "Water Sampling") %>%
      group_by(site) %>%
      summarise(visits = n(), .groups = "drop")
  }

  cal_dates <- field_notes %>%
    filter(!grepl("none", cals_performed, ignore.case = TRUE), !is.na(cals_performed)) %>%
    select(date, site, cals_performed) %>%
    mutate(date_mmdd = format(date, "%m/%d")) %>%
    group_by(site) %>%
    summarise(cal_dates = paste(date_mmdd, collapse = ", "), .groups = "drop")

  site_labels <- cfg$site_labels
  notes_final <- tibble(site = cfg$sites) %>%
    left_join(monthly_visits, by = "site") %>%
    left_join(cal_dates, by = "site") %>%
    mutate(
      `Site Visits` = ifelse(is.na(visits), 0, visits),
      `Calibration Dates` = ifelse(is.na(cal_dates), "None", cal_dates),
      Site = if (!is.null(site_labels)) unname(site_labels[site]) else toupper(site),
      Site = factor(Site, levels = if (!is.null(site_labels)) unname(site_labels[cfg$sites]) else toupper(cfg$sites))
    ) %>%
    arrange(Site) %>%
    select(Site, `Site Visits`, `Calibration Dates`)

  if (full) {
    return(field_notes %>% arrange(site, DT_round))
  }
  notes_final
}

# -----------------------------------------------------------------------------
# Path helpers -- keep all file-naming conventions in one place -------------
# -----------------------------------------------------------------------------
path_field_summary <- function(cfg, month_select, year_select) {
  file.path(cfg$output_dir, "field_summary",
            paste0("field_summary_", month_select, "_", year_select, ".csv"))
}
path_charges_final <- function(cfg, month_select, year_select) {
  file.path(cfg$output_dir, "charges_final",
            paste0("charges_final_", month_select, "_", year_select, ".csv"))
}

#--------------------------------------------------------------------------
# Monthly site-wise summary for PWQN sites
#-------------------------------------------------------------------------
generate_monthly_boxplot_summary <- function(sensor_data, year, month){

  # Parse the month string as part of a date and extract the month number
  month_number <- mdy(paste(month, "1, 2020")) %>% lubridate::month()
  month_name <- mdy(paste(month, "1, 2020")) %>% lubridate::month(label = TRUE)


  box_data <- sensor_data%>%
    filter(parameter %nin%  c("Turbidity", "Depth", "FDOM Fluorescence", "Chl-a Fluorescence", "ORP"))%>%
    filter(site %in% c("bellvue", "salyer", "udall", "riverbend", "cottonwood", "elc", "archery", "riverbluffs"))%>%
    filter(lubridate::month(DT_round) == month_number & lubridate::year(DT_round) == year)


  # Merge with your data using left_join
  box_data_with_limits <- box_data %>%
    mutate(natural_name = factor(natural_name, levels = site_order ))


  convert_season_to_months <- function(season_range) {
    # Define month order
    months <- c("january", "february", "march", "april", "may", "june",
                "july", "august", "september", "october", "november", "december")

    # Split the season range
    season_parts <- strsplit(season_range, "-")[[1]]
    start_month <- tolower(season_parts[1])
    end_month <- tolower(season_parts[2])

    # Find positions of start and end months
    start_pos <- which(months == start_month)
    end_pos <- which(months == end_month)

    # Handle cases where season crosses year boundary
    if (start_pos <= end_pos) {
      selected_months <- months[start_pos:end_pos]
    } else {
      selected_months <- c(months[start_pos:12], months[1:end_pos])
    }

    return(paste(selected_months, collapse = "|"))
  }

  #read in cdphe standards collected by Sam Struthers from https://cdphe.colorado.gov/water-quality-control-commission-regulations
  #Parameter (DO, pH, Temp) thresholds: Reg 31, Table I
  # Site classifications: Reg 38, ~ Pg 410

  month_standard <- readxl::read_xlsx(path = here("data","raw", "spatial", "metadata", "site_cdphe_classification_2025.xlsx")) %>%
    # Convert season_1 to pipe-separated months
    mutate(season = sapply(season, convert_season_to_months))%>%
    # Split ph column into ph_low and ph_high
    separate(ph, into = c("ph_low", "ph_high"), sep = "-", convert = TRUE) %>%
    filter(grepl(x = season, pattern = month_name, ignore.case = TRUE),
           site %in% unique(box_data_with_limits$natural_name))

  # Get unique parameters from the data
  parameters <- unique(box_data_with_limits$label)

  # Create individual plots for each parameter
  parameter_plots <- list()

  for (param in parameters) {
    # Filter data for current parameter
    param_data <- box_data_with_limits %>%
      filter(label == param)

    standards <- param_data %>%
      select(natural_name) %>%
      distinct() %>%
      left_join(month_standard, by = c("natural_name" = "site"))%>%
      mutate(site_num =  as.numeric(factor(natural_name, levels = site_order))-1)

    # Create base plot for current parameter
    p <- param_data %>%
      ggplot(aes(x = natural_name, y = mean, fill = natural_name)) +
      geom_boxplot() +
      #geom_blank(aes(y = y_min)) +  # Forces inclusion of minimum
      #geom_blank(aes(y = y_max)) +  # Forces inclusion of maximum
      scale_fill_manual(values = bg_colors_full) +
      labs(x = "Site", y = param) +  # Use parameter name as y-axis label
      theme_bw(base_size = 14) +  # Slightly smaller base size for multiple plots
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.x = element_text(face = "bold"),
            axis.title.y = element_text(face = "bold"),
            plot.title = element_text(size = 12))  # Smaller title for individual plots

    # Add reference lines based on parameter type
    if (grepl("pH", param, ignore.case = TRUE)) {
      # Add horizontal lines for pH standards
      p <- p +
        geom_hline(yintercept = 6.5, color = "red", linetype = "dashed", alpha = 0.7) +
        geom_hline(yintercept = 9, color = "red", linetype = "dashed", alpha = 0.7)
    }

    if (grepl("Temperature|Temp", param, ignore.case = TRUE)) {
      # Add chronic temperature lines (orange, dashed)
      if ("temp_chronic" %in% names(standards)) {
        chronic_data <- standards %>%
          filter(!is.na(temp_chronic))

        if (nrow(chronic_data) > 0) {
          p <- p +
            geom_segment(data = chronic_data,
                         aes(x = site_num - 0.5, xend = site_num + 0.5,
                             y = temp_chronic, yend = temp_chronic),
                         color = "orange", linetype = "dashed", linewidth = 1, alpha = 0.8,
                         inherit.aes = FALSE)
        }
      }

      # Add acute temperature lines (red, solid)
      if ("temp_acute" %in% names(standards)) {
        acute_data <- standards %>%
          filter(!is.na(temp_acute))

        if (nrow(acute_data) > 0) {
          p <- p +
            geom_segment(data = acute_data,
                         aes(x = site_num - 0.5, xend = site_num + 0.5,
                             y = temp_acute, yend = temp_acute),
                         color = "red", linetype = "solid", linewidth = 1, alpha = 0.8,
                         inherit.aes = FALSE)
        }
      }
    }
    if (grepl("DO", param, ignore.case = TRUE)) {

      # Add chronic temperature lines (orange, dashed)
      if ("do_chronic" %in% names(standards)) {
        do_chronic_data <- standards %>%
          filter(!is.na(temp_chronic))

        if (nrow(do_chronic_data) > 0) {
          p <- p +
            geom_segment(data = do_chronic_data,
                         aes(x = site_num - 0.4, xend = site_num + 0.4,
                             y = do_chronic, yend = do_chronic),
                         color = "red", linetype = "solid", linewidth = 1, alpha = 0.8,
                         inherit.aes = FALSE)
        }
      }
    }

    parameter_plots[[param]] <- p +
      annotate("text",
               x = -Inf,
               y = Inf,
               label = "PRELIMINARY RESULTS",
               hjust = -0.1,
               vjust = 2,
               size = 4,
               fontface = "bold",
               color = "red")
  }


  parameter_plots[[1]] <- parameter_plots[[1]]+
    theme(axis.title.x = element_blank())
  parameter_plots[[2]] <- parameter_plots[[2]]+
    theme(axis.title.x = element_blank())

  # Combine plots using patchwork
  combined_plot <- wrap_plots(parameter_plots, ncol = 2) +  # Adjust ncol as needed
    plot_annotation(
      title = paste0("Poudre Water Quality Network Data Summary: ", str_to_title(month_name),  " ", year),
      caption = "Preliminary Data, subject to revision.\n Red solid lines indicate state of Colorado chronic thresholds and orange dashed lines indicate acute thresholds for aquatic life.",
      theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
                    plot.caption = element_text(size = 13))
    )

  # Display the combined plot
  return(combined_plot)

}
