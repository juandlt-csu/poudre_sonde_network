

#### Server ####

server <- function(input, output, session) {
  #### Not used in internal version ####
  #
  #
  #   output$uploaded_files <- renderTable({
  #     req(input$upload)
  #
  #     # Get the uploaded file
  #     zip_file <- input$upload$datapath
  #     # Create a unique temporary directory for this upload
  #     temp_dir <- file.path(tempdir(), paste0("upload_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  #     dir.create(temp_dir)
  #
  #     # Unzip the file to the temporary directory
  #     unzip(zip_file, exdir = temp_dir)
  #     files <- list.files(
  #       path = temp_dir,
  #       recursive = TRUE,
  #       full.names = TRUE
  #     )%>%
  #     #remove the temp_dir from the file path
  #     gsub(paste0(temp_dir, "/data/"), "", .)
  #     # Remove metadata folder entries
  #     # Extract folder and site-parameter combinations
  #     file_info <- data.frame(
  #       Folder = dirname(files),
  #       Site_Parameter = basename(files) %>%
  #         tools::file_path_sans_ext()  # Remove .csv extension
  #     ) %>%
  #       # Remove any empty folders or "."
  #       filter(Folder != "." & Folder != "" ) %>%
  #       filter(!grepl("meta", Folder)) %>%
  #       # Sort by folder and site-parameter
  #       arrange(Folder, Site_Parameter)%>%
  #       # Split Site_Parameter into Site and Parameter
  #       tidyr::separate(Site_Parameter,
  #                       into = c("Site", "Parameter"),
  #                       sep = "-",
  #                       remove = TRUE)%>%
  #       # Pivot wider to make folders as columns
  #       pivot_wider(
  #         names_from = Folder,
  #         values_from = Parameter,
  #         values_fn = function(x) paste(unique(x), collapse = ","),
  #         values_fill = ""
  #       ) %>%
  #       arrange(Site)
  #
  #
  # file_info
  #
  #     })
  #### Not used in internal version ####

  #### Reactive values ####
  data <- reactiveVal(NULL)
  current_week <- reactiveVal(NULL) #controlled by next/prev week buttons, and submit weekly decision
  selected_data <- reactiveVal(NULL) # This is essentially site param df
  all_datasets <- reactiveVal(NULL) # List of all datasets and is used in generating sub plots? (why is there a question mark here)
  brush_active <- reactiveVal(FALSE) #internal shiny tracker for brush tool
  selected_data_cur_filename <- reactiveVal(NULL) # Current filename of selected data, to be updated as filename is saved

  auto_refresh <- reactiveTimer(30000) #refresh every 30 sec

  # Check if data folder exists and if data/all_data subfolder has files, if files are available, show table of available files and allow user selection
  output$conditional_data_ui <- renderUI({
    # Check if data folder exists and if data/all_data subfolder has files
    data_folder_exists <- dir.exists(in_progress_path)
    all_data_subfolder_exists <- dir.exists(all_data_path)
    all_data_subfolder_empty <- FALSE
      if(!all_data_subfolder_exists | length(list.files(all_data_path)) == 0) {
        all_data_subfolder_empty <- TRUE
      }

    if(!data_folder_exists|all_data_subfolder_empty){
      # Show file upload and timezone input if conditions are met
      tagList(
        textInput("timezone", "Enter your timezone:", value = "MST"),
        fileInput("data_upload", "Upload your data files:", multiple = FALSE,
                  accept = c(".csv", ".xlsx", ".zip", ".feather", ".rds", ".parquet"))
      )
    } else {
      #Show regular UI if files are present
      tagList(
        DT::dataTableOutput("data_files_table"),
        # Directory selection
        fluidRow(
          column(3,
                 selectInput("directory", "Choose Directory:",
                             choices = c("pre_verification", "intermediary"),
                             selected = "pre_verification")
          ),
          # User Selection
          column(3,
                 selectInput("user", "Select User:",
                             choices = c("SJS", "MNR","DRD", "JDT", "KW", "BS"),
                             selected = "SJS")
          ),
          # Site Selection
          column(3,
                 selectInput("site", "Select Site:", choices = NULL)
          ),
          # Parameter Selection
          column(3,
                 selectInput("parameter", "Select Parameter:", choices = NULL)
          )
        ),
        br(),

        actionButton("load_data", "Load Data", class = "btn-primary"),

      )
    }
  })

  #constantly updating file paths for shiny app
  all_filepaths <- reactive({

    sync_file_system()

    get_filenames()%>%mutate(
      site = map_chr(filename, ~ split_filename(.x)$site),
      parameter = map_chr(filename, ~ split_filename(.x)$parameter),
      datetime = map_chr(filename, ~ split_filename(.x)$datetime))

  })


  #data table for available parameters and which folder they belong to
  output$data_files_table <- renderDataTable({
    files <- all_filepaths() %>%
      filter(directory %in% c("pre_verification", "intermediary", "verified")) %>%
      select(-filename, -datetime) %>%
      mutate(parameter = gsub("_FINAL", "", parameter) )%>% #remove _FINAL from parameter names to make it play nice with final directory
      distinct() %>%
      pivot_wider(
        names_from = parameter,
        values_from = directory,
        values_fill = NA
      ) %>%
      filter(if_any(-site, ~ .x %in% c("pre_verification", "intermediary"))) %>%
      arrange(site)

    DT::datatable(files, options = list(pageLength = 25,
                                        scrollY = "400px")) %>%
      DT::formatStyle(
        columns = names(files)[-1],  # All columns except the first (site)
        valueColumns = names(files)[-1],
        backgroundColor = DT::styleEqual(
          c("pre_verification", "intermediary", "verified"),
          c("#D3D3D3", "#FFA500", "#90EE90")
        )
      )
  })

  # Observer to handle data uploads and timezone input
  observeEvent(input$data_upload, {
    req(input$data_upload, input$timezone)

    # Check if both data_upload has files and timezone has a value
    if (!is.null(input$data_upload) & !is.null(input$timezone) & input$timezone != "") {

      upload_path <- input$data_upload$datapath
      # Call the setup function with the uploaded files and timezone
      # and capture its return value
      result_message <- setup_directories_from_upload(
        uploaded_file_path = upload_path,
        timezone = input$timezone
      )

      # Force auto-refresh to update the UI
      session$reload()

      # Show notification with the returned message
      showNotification(
        result_message,
        type = "message",
        duration = 5
      )
    }
  })

  #### Data Selection functions ####

  # Update site choices when directory changes
  observe({
    req(input$directory, all_filepaths())
    sites <- all_filepaths() %>%
      filter(directory == input$directory) %>%
      pull(site) %>%
      unique()

    updateSelectInput(session, "site",
                      choices = sites) #based on directory
  })

  # Update parameter choices when site changes
  observe({
    req(input$directory, input$site, all_filepaths())

    parameters <- all_filepaths() %>% filter(directory == input$directory & site == input$site)%>%
      pull(parameter) %>%
      unique()

    updateSelectInput(session, "parameter",
                      choices = parameters) #based on site and directory
  })

  # Show/hide and update sub parameters UI based on main parameter selection
  observe({
    req(input$parameter)

    # Get auto-selected parameters for the chosen parameter
    auto_params <- get_auto_parameters(input$parameter)

    # Update sub-parameters selection
    updateSelectInput(session, "sub_parameters",
                      choices = available_parameters,
                      selected = auto_params)
  })
  #Show/hide and update additional sites UI based on site selection
  observe({
    req(input$site)

    # Get auto-selected parameters for the chosen parameter
    auto_sites <- relevant_sonde_selector(site_arg = input$site)

    # Update sub-parameters selection
    updateSelectInput(session, "sub_sites",
                      choices = available_sites,
                      selected = auto_sites)
    updateSelectInput(session, "add_sites",
                      choices = available_sites,
                      selected = auto_sites)

  })


  # Load data when button is clicked
  observeEvent(input$load_data, {
    req(input$directory, input$site, input$parameter)
    # Initialize data directories and load datasets

    sync_file_system()

    # TODO: This loads all data and is probably inefficient, should be updated to only load the data needed (efficiency notes)
    datasets <- load_all_datasets()

    datasets <- map(datasets, function(data_list) {
      # Return immediately if data_list is empty
      if (is_empty(data_list)) {
        return(data_list)
      }

      # Extract filenames from list names and maintain their original order
      file_names <- tibble(list_name = names(data_list)) %>%
        mutate(split_data = map(list_name, split_filename)) %>%
        unnest_wider(split_data) %>%
        mutate(parameter = gsub("_FINAL", "", parameter) )%>% #remove _FINAL from parameter names to make it play nice with final directory
        mutate(site_param = paste(site, parameter, sep = "-")) %>%
        group_by(site_param) %>%
        arrange(desc(datetime)) %>%  # Sort so the latest entry remains unchanged
        mutate(site_param = if_else(row_number() == 1, site_param, paste0(site_param, "_backup"))) %>%
        ungroup()

      # Ensure names are applied in the correct order by matching filenames
      names(data_list) <- file_names$site_param[match(names(data_list), file_names$filename)]

      return(data_list)
    })


    all_datasets(datasets)

    # Get the site-parameter name
    site_param_name <- paste0(input$site, "-", input$parameter)

    # Try to get the specific dataset
    tryCatch({

      if(input$directory == "pre_verification") {
        site_param_df <- datasets$pre_verification_data[[site_param_name]]

        if (is.null(site_param_df)) {
          stop(paste("Dataset", site_param_name, "not found"))
        }

        predata_file_name <- all_filepaths() %>%
          filter(site == input$site & parameter == input$parameter & directory == input$directory) %>%
          pull(filename)

        selected_data_cur_filename(move_file_to_intermediary_directory(pre_to_int_filename = predata_file_name, pre_to_int_df = site_param_df))

        #refresh all_datafiles
      } else {
        site_param_df <- datasets$intermediary_data[[site_param_name]]

        if (is.null(site_param_df)) {
          stop(paste("Dataset", site_param_name, "not found"))
        }

        int_file <- all_filepaths() %>%
          filter(site == input$site & parameter == input$parameter & directory == input$directory) %>%
          #grab the most recent version!
          arrange(desc(datetime))%>%
          slice(1)%>%
          pull(filename)

        selected_data_cur_filename(update_intermediary_data(int_file, site_param_df))
      }

      # Store the processed data
      # selected data is instantiated
      selected_data(site_param_df)

      # Set initial week to earliest week with missing final_status values (unverified) - stores data as final values and automatically sends you to the only data without final status
    first_non_ver_week <-   site_param_df$week[min(which(is.na(site_param_df$final_status)))]

    #if all weeks are verified but data is not "finalized" ^ will be NA and we should move directly to the finalize page
    if(is.na(first_non_ver_week)| is.null(first_non_ver_week)){
      showNotification(
        "All weeks verified! Moving to Finalize Data tab.",
        type = "message"
      )
      #Move to finalize tab
      updateTabsetPanel(session, inputId = "tabs", selected = "Finalize Data")
    }else{
      current_week(first_non_ver_week)
      #Move to next tab
      updateTabsetPanel(session, inputId = "tabs", selected = "Data Verification")
    }



    }, error = function(e) {
      showNotification(
        paste("Error loading data:", e$message),
        type = "error"
      )
    })



  })

  #### Data Verification functions ####
  # Previous Tab
  observeEvent(input$prev_tab, {
    updateNavbarPage(session, "tabs", selected = "Data Selection")

    #Q: Should this update the data files or no?
  })

  ## Week navigation handlers
  observeEvent(input$prev_week, {
    req(selected_data())
    session$resetBrush("plot_brush")
    brushed_areas(list())
    weeks <- unique(selected_data()$week)
    current <- current_week()
    idx <- which(weeks == current)
    if (idx > 1) {
      current_week(weeks[idx - 1])
    }
  })
  # Go to next week
  observeEvent(input$next_week, {
    req(selected_data())
    session$resetBrush("plot_brush")
    brushed_areas(list())
    weeks <- unique(selected_data()$week)
    current <- current_week()
    idx <- which(weeks == current)
    if (idx < length(weeks)) {
      current_week(weeks[idx + 1])
    }else{
      showNotification(
        "No more weeks to verify. Click Final Verification to see unverified weeks",
        type = "warning"
      )
      current <- current_week()
      idx <- which(weeks == current)
    }
  })

  # Reactive to hold USGS flow data
  usgs_flow_data <- reactive({
    req(input$site, current_week(), isolate(selected_data()))

    week_data <- isolate(selected_data()) %>% filter(week == current_week())
    week_min_day <- min(week_data$DT_round, na.rm = TRUE)
    week_max_day <- max(week_data$DT_round, na.rm = TRUE)

    all_flow_sites <- input$site
    flow_data_list <- purrr::map(all_flow_sites, function(site_name) {
      abbrev_vals <- if (site_name %in% c("pbd", "bellvue", "pman", "pbr", "pfal", "cbri", "joei")) {
        "CLAFTCCO"
      } else if (site_name %in% c("salyer", "udall", "riverbend", "riverbend_virridy")) {
        "CLAFORCO"
      } else if (site_name %in% c("cottonwood", "cottonwood_virridy", "elc", "archery", "archery_virridy", "boxcreek", "springcreek", "riverbluffs")) {
        "CLABOXCO"
      } else if (site_name %in% c("sfm")) {
        c("CLASRKCO", "CLAFTCCO")
      } else if (site_name %in% c("chd")) {
        c("JWCCHACO", "CLAFTCCO")
      } else {
        NA_character_
      }
      
      if (all(is.na(abbrev_vals))) return(NULL)

      tryCatch({
        s_date <- week_min_day - days(2)
        e_date <- week_max_day + days(2)

        if (!is.null(global_usgs_flow_data) && nrow(global_usgs_flow_data) > 0) {
          # Use pre-fetched global data
          
          # Safely determine the column name (station_abbrev or abbrev)
          col_name <- if("station_abbrev" %in% names(global_usgs_flow_data)) "station_abbrev" else "abbrev"
          
          for (abbrev_val in abbrev_vals) {
            data <- global_usgs_flow_data %>%
              filter(!!sym(col_name) == abbrev_val,
                     datetime >= s_date,
                     datetime <= e_date)
                     
            if (nrow(data) > 0) {
              data$site <- site_name
              data$abbrev <- abbrev_val
              return(data)
            }
          }
        }
        
        return(NULL)
      }, error = function(e) NULL)
    }) %>% purrr::compact()

    if (length(flow_data_list) > 0) bind_rows(flow_data_list) else NULL
  })

  ## Main plot (main plot starts here)
  output$main_plot <- renderPlot({

    # requirements for this function to work
    req(selected_data(), current_week(), all_datasets(), input$weekly_decision)

    # TODO: See previous note, this should be loaded only once rather than each time the plot updates (efficiency notes)
    pre_verification_data <- all_datasets()[["pre_verification_data"]]
    intermediary_data <- all_datasets()[["intermediary_data"]]
    verified_data <- all_datasets()[["verified_data"]]


    # user selected site parameter combo
    week_data <- selected_data() %>%
      filter(week == current_week())

    week_min_day = min(week_data$DT_round, na.rm = T)
    week_max_day = max(week_data$DT_round, na.rm = T)
    #grab two days on either end of the week
    week_plus_data <- selected_data() %>%
      filter(DT_round >= week_min_day - days(2) & DT_round <= week_max_day + days(2))
    #remove data from the week

    #Setup year week (eg "2024 - 13")
    year_week <- paste0(as.character(year(week_min_day)) ," - ", current_week())
    flag_day <- week_min_day


    #This isn't being sourced correctly in the shiny app
    #TODO: Source this correctly
    #picking the correct name of the data frame
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

    }

    plot_filter <- input$add_sites # remember to check how things interact with UI -JD
    # Get the relevant sonde data THIS IS GOOD, no touch for now -JD
    relevant_sondes <- map(plot_filter, ~ {
      sonde_name <- paste0(.x, "-", input$parameter)
      data_source <- NULL
      sonde_df <- NULL

      # Determine which directory to pull data from
      tryCatch({
        data_source <- retrieve_relevant_data_name(sonde_name, year_week)
        # cat("Data for",sonde_name,"will be pulled from",data_source,"\n")
      }, error = function(err) {
        cat("Data for",sonde_name,"not found.\n")
        return(NULL)  # Return NULL if data source can't be determined
      })

      # Only try to pull in the data if data_source was successfully determined
      if (!is.null(data_source)) {
        tryCatch({
          sonde_df <- get(data_source)[[sonde_name]] %>%
            filter(DT_round >= week_min_day - days(2) & DT_round <= week_max_day + days(2))
        }, error = function(err) {
          #cat("Sonde", sonde_name, "not found.\n")
          return(NULL)  # Return NULL if sonde data can't be retrieved
        })
      }

      # Only return a list if both data_source and sonde_df are available
      if (!is.null(data_source) & !is.null(sonde_df)) {
        return(list(sonde_df = sonde_df, data_source = data_source))
      } else {
        return(NULL)  # Return NULL if either data_source or sonde_df is NULL
      }

    })

    # Remove any NULL results from the list
    relevant_sondes <- compact(relevant_sondes)

    # append week_data to relevant sonde list, clean list, and bind dfs
    relevant_dfs <- map(relevant_sondes, ~.x[[1]])
    week_plot_data <- append(relevant_dfs, list(week_data)) %>%
      append(., list(week_plus_data)) %>%
      keep(~ !is.null(.)) %>%
      keep(~ nrow(.)>0) %>%
      bind_rows() %>%
      arrange(day)


    # Check the decision and create appropriate plot
    # This is where the tweaks will probably happen
    # What is this nested inside of: output$main_plot <- renderPlot
    if (input$weekly_decision != "s") { # make sure these are interacting with the UI correctly
      # Show final decision to user to preview weekly decision
      week_choice_data <- week_data %>%
        mutate(
          final_decision = case_when(
            #AA:Pass all data
            input$weekly_decision == "aa"  ~ "PASS", # how do shiny apps update reactive values. if they are in renderPlot, do they rerun render plot automatically
            #ANO: Accept Non Omit
            input$weekly_decision == "ano" & !brush_omit ~ "PASS", # pass data that is not user select omit
            #KF: Keep Flags, data becomes flagged but kept in dataset (mildly sus)
            input$weekly_decision == "kf" & is.na(user_flag) & !brush_omit ~ "PASS", # pass data that is not user select omit
            input$weekly_decision == "kf" & !is.na(user_flag) & !brush_omit ~ "FLAGGED", # tag data that is flagged
            #OF: Omit Flagged
            input$weekly_decision == "of" & is.na(user_flag) & !brush_omit ~ "PASS", # pass data that is not user select omit
            input$weekly_decision == "of" & !is.na(user_flag) & !brush_omit ~ "OMIT", # omit data that is flagged
            #OA: Omit All
            input$weekly_decision == "oa"  ~ "OMIT",
            # Omit any user selected omit data (assuming AA was not the choice)
            input$weekly_decision != "aa" & brush_omit ~ "OMIT"))

      #Remove omitted data (user or from weekly decision)
      if (("remove_omit" %in% input$plot_options)) {
        week_choice_data <- week_choice_data %>%
          filter(final_decision != "OMIT")

        week_plus_data <- week_plus_data %>%
          filter(!brush_omit)
      }
      #Remove flagged data if user desired
      if(("remove_flag" %in% input$plot_options)){
        week_choice_data <- week_choice_data %>%
          filter(final_decision != "FLAGGED")

        week_plus_data <- week_plus_data %>%
          filter(final_status != "FLAGGED" | is.na(final_status))
      }

      week_min_check <- week_plus_data %>%
        filter(week < current_week())

      week_max_check <- week_plus_data %>%
        filter(week > current_week())

      # Create plot for preview of weekly decision
      # actually building plot
      p <- ggplot(week_choice_data, aes(x = DT_round))

      if (nrow(week_min_check) > 0) {

        week_plus_min <- week_min_check %>%
          summarise(xmin = min(DT_round), xmax = max(DT_round))

        p <- p + annotate(geom = "rect", xmin = week_plus_min$xmin, xmax = week_plus_min$xmax, ymin = -Inf, ymax = Inf, color = "transparent", fill = "grey", alpha = 0.2)
      }
      if (nrow(week_max_check) > 0) {
        week_plus_max <- week_max_check %>%
          summarise(xmin = min(DT_round), xmax = max(DT_round))

        p <- p + annotate(geom = "rect", xmin = week_plus_max$xmin, xmax = week_plus_max$xmax, ymin = -Inf, ymax = Inf, color = "transparent", fill = "grey", alpha = 0.2)
      }

      p <- p +
        # Add relevant sites
        map(relevant_sondes, function(sonde_data) {
          add_data <- sonde_data[[1]]
          data_source <- sonde_data[[2]]

          y_column <- ifelse(data_source %in% c("all_data", "pre_verification_data"), "mean", "mean_verified")

          # Add interpolation for a hour to try to improve isolated point issues
          # If we have missing data for a single point, the line breaks and if there is no point on the other side
          # of it is just disappears
          # Adding in a na interpolation of 3 points (45 min) so that we can still "see" the trends of other sites that might have alot
          # of flags
          add_data_with_interpolation <- add_data %>%
            arrange(site, DT_round) %>%
            mutate(
              !!sym(y_column) := zoo::na.approx(
                x = as.numeric(DT_round),
                object = !!sym(y_column),
                maxgap = 3, # interpolate over 45 min
                na.rm = FALSE
              )
            ) %>%
            ungroup()

          # Return a list of both geom_line
          list(
            geom_line(data = add_data_with_interpolation,
                      aes(x = DT_round, y = .data[[y_column]], color = site),
                      linewidth = 1, na.rm = TRUE)
          )
        }) + #plot other sites
        geom_point(aes(y = mean, fill = final_decision),shape = 21, stroke = 0, size = 2) #plot main site with colors matching final decision

      if(("add_line" %in% input$plot_options)){
        p <- p + geom_line(data = week_choice_data %>% filter(!is.na(mean)), aes(y = mean), color = "grey", linewidth = 1, alpha = 0.4)
      }

      #if incl_ex_days = T, then add in the extra data as points
      if(("incl_ex_days" %in% input$plot_options)){
        p <- p +
          geom_point(data = week_plus_data%>%filter(week != current_week()),
                     aes(y = mean, fill = final_status), shape = 21, stroke = 0, size = 1.5, alpha = 0.5)
      }
      # Do final styling
      p <- p +
        scale_fill_manual(values = final_status_colors, na.value = "grey")+ #set fill colors to weekly status colors
        scale_color_manual(values = setNames(site_color_combo$color, site_color_combo$site))+ #set other colors for add sites
        labs(
          title = paste0(str_to_title(input$site), " ", input$parameter, " (", format(flag_day, "%B %d, %Y"), ")"),
          x = "Date",
          y = input$parameter,
          fill = "Preview of Final Decision",
          color = "Sites" )+
        theme_bw(base_size = 14)

      if(("incl_thresholds" %in% input$plot_options)){
        p <- add_threshold_lines(plot = p,
                                 plot_data = week_plot_data,
                                 site_arg = input$site,
                                 parameter_arg = input$parameter)
      }
      # This is where the xlim is set, a problem -JD
      if(("incl_ex_days" %in% input$plot_options)){
        p <- p +
          scale_x_datetime(
            limits = c(min(week_plus_data$DT_round, na.rm = TRUE),
                       max(week_plus_data$DT_round, na.rm = TRUE)),
            #expand = c(0, 0), # Remove extra white space
            date_labels = "%b %d", # Formats as "Jan 01", "Feb 15", etc.
            date_breaks = "1 day" # Remove extra white space
          )
      }else{
        p <- p +
          scale_x_datetime(
            limits = c(min(week_choice_data$DT_round, na.rm = TRUE),
                       max(week_choice_data$DT_round, na.rm = TRUE)),
            #expand = c(0, 0), # Remove extra white space
            date_labels = "%b %d", # Formats as "Jan 01", "Feb 15", etc.
            date_breaks = "1 day" # Remove extra white space
          )
      }

      if(("plot_log10" %in% input$plot_options)){
        p <- p + scale_y_log10()
      }

      if(!("show_legend" %in% input$plot_options)) {
        p <- p + theme(legend.position = "none")
      }

      p
    } else {
      # TODO: Swap with create weekly plot function call, adding in other sites, etc

      if(("remove_omit" %in% input$plot_options)){
        week_data <- week_data %>%
          filter(!brush_omit) # remove omitted data

        week_plus_data <- week_plus_data %>%
          filter(!brush_omit)
      }
      if(("remove_flag" %in% input$plot_options)){
        week_data <- week_data %>%
          filter(is.na(user_flag) & !brush_omit) # remove flagged data unless it is omitted - why keep omitted data?

        week_plus_data <- week_plus_data %>%
          filter( is.na(user_flag) & !brush_omit)
      }

      p <-ggplot(week_data, aes(x = DT_round))

      #Adding in extra days (+- 2 days on each side)
      if(("incl_ex_days" %in% input$plot_options)){
        # check to make sure there is data from last week
        week_min_check <- week_plus_data %>%
          filter(week < current_week())
        #check to see if there is data from next week
        week_max_check <- week_plus_data %>%
          filter(week > current_week())

        if (nrow(week_min_check) > 0) { #if there is data from last week, add in the grey box for the last X days

          week_plus_min <- week_min_check %>%
            summarise(xmin = min(DT_round), xmax = max(DT_round))

          p <- p + annotate(geom = "rect", xmin = week_plus_min$xmin, xmax = week_plus_min$xmax, ymin = -Inf, ymax = Inf, color = "transparent", fill = "grey", alpha = 0.2)
        }
        if (nrow(week_max_check) > 0) { #if there is data for next week, add in the grey box for the next X days
          week_plus_max <- week_max_check %>%
            summarise(xmin = min(DT_round), xmax = max(DT_round))

          p<- p + annotate(geom = "rect", xmin = week_plus_max$xmin, xmax = week_plus_max$xmax, ymin = -Inf, ymax = Inf, color = "transparent", fill = "grey", alpha = 0.2)
        }


      }

      #TODO: Not sure if this is working how I'd expect
      #set the value for the geom_crossbar, should adjust in size based on range of data
      adjustment_value = sd(week_data$mean, na.rm = TRUE)*.05

      p <- p +
        map(relevant_sondes, function(sonde_data) {
          add_data <- sonde_data[[1]]
          data_source <- sonde_data[[2]]

          # Deciding which column to use based on the verification status
          # Unverified data should use "mean" and verified data should use "mean_verified"
          y_column <- ifelse(data_source %in% c("all_data", "pre_verification_data"), "mean", "mean_verified")

          # Add interpolation for a hour to try to improve isolated point issues
          # If we have missing data for a single point, the line breaks and if there is no point on the other side
          # of it is just disappears
          # Adding in a na interpolation of 3 points (45 min) so that we can still "see" the trends of other sites that might have alot
          # of flags
          add_data_with_interpolation <- add_data %>%
            arrange(site, DT_round) %>%
            mutate(
              !!sym(y_column) := zoo::na.approx(
                x = as.numeric(DT_round),
                object = !!sym(y_column),
                maxgap = 3, # interpolate over 45 min
                na.rm = FALSE
              )
            ) %>%
            ungroup()

          #filter to just single week if incl_ex_days is false
          if(!("incl_ex_days" %in% input$plot_options)){
            add_data_with_interpolation <- add_data_with_interpolation %>%
              filter(week == current_week())
          }

          # Return a list of both geom_line and geom_linerange
          list(
            geom_line(data = add_data_with_interpolation,
                      aes(x = DT_round, y = .data[[y_column]], color = site),
                      linewidth = 1, na.rm = TRUE)
          )
        })

      # add in extra days if incl_ex_days is true
      if(("incl_ex_days" %in% input$plot_options)){
        p <- p + geom_point(data = week_plus_data%>%filter(week != current_week()), aes(y = mean),fill = "black",shape = 21, stroke = 0, size = 1.5, alpha = 0.5) #add two extra days on the side
      }

      #add in primary data
      p <- p+ geom_point(aes(y = mean, fill = user_flag),shape = 21, stroke = 0, size = 2)+ #plot main site with colors matching  user flag column
        #Add Omitted data in red
        geom_point(data = week_data %>%filter(brush_omit == TRUE),aes(y = mean),shape = 21, stroke = 0, size = 2, fill = "#ff1100")+
        scale_color_manual(
          name = "Sites",
          values = setNames(site_color_combo$color, site_color_combo$site)) +
        # Scale for points (flag-based colors)
        scale_fill_viridis_d(
          name = "Flags",
          option = "plasma",
          begin = 0.1,
          end = 0.9,
          na.value = "grey" )+
        scale_x_datetime(
          date_labels = "%b %d", # Formats as "Jan 01"
          date_breaks = "1 day" # Remove extra white space
        )+
        labs(
          title = paste0(str_to_title(input$site), " ", input$parameter, " (", format(flag_day, "%B %d, %Y"), ")"),
          x = "Date",
          y = input$parameter,
          color = "Sites",
          fill = "Flags")+
        theme_bw(base_size = 14)

      if(("add_line" %in% input$plot_options)){
        p <- p + geom_line(data = week_data %>% filter(!is.na(mean)), aes(y = mean), color = "grey", linewidth = 1)
      }


      # Check if there are any brushed areas, THIS IS GOOD -JD
      if(length(brushed_areas()) > 0) {

        # Create a data frame of all brush boundaries
        brush_boxes <- map_dfr(
          brushed_areas()[seq(1, length(brushed_areas()), by = 2)], #for some reason this likes to duplicate values so just grabbing half of them
          ~data.frame(
            xmin_DT = .x$brush_dt_min,
            xmax_DT = .x$brush_dt_max,
            ymin_mean = .x$brush_mean_min,
            ymax_mean = .x$brush_mean_max)
        )
        #   # Add rectangles to plot for all brushed areas
        p <- p +
          geom_rect(data = brush_boxes,
                    aes(xmin = xmin_DT,
                        xmax = xmax_DT,
                        ymin = ymin_mean,
                        ymax = ymax_mean),
                    fill = NA,
                    color = "blue",
                    alpha = 0.3, inherit.aes = F)

      }
      #create plot
      if(("incl_thresholds" %in% input$plot_options)){
        p <- add_threshold_lines(plot = p,
                                 plot_data = week_plot_data,
                                 site_arg = input$site,
                                 parameter_arg = input$parameter)
      }

      if(("plot_log10" %in% input$plot_options)){
        p <- p + scale_y_log10()
      }

      if(!("show_legend" %in% input$plot_options)) {
        p <- p + theme(legend.position = "none")
      }

      # Return the plot from renderPlot
      p
    }
  })
# main plotting code ends here

  ## Sub plots output
  # Anything that is flagged/omitted is removed, we want to keep flags and differentiate those points somehow -JD
  output$sub_plots <- renderPlotly({
    req(all_datasets(), current_week(), input$site)

    pre_verification_data <- all_datasets()[["pre_verification_data"]]
    intermediary_data <- all_datasets()[["intermediary_data"]]
    verified_data <- all_datasets()[["verified_data"]]

    week_data <- selected_data() %>%
      filter(week == current_week())

    week_min_day = min(week_data$DT_round, na.rm = T)
    week_max_day = max(week_data$DT_round, na.rm = T)


    year_week <- paste0(as.character(year(week_min_day)), " - ", current_week())
    all_sub_sites <- c(input$site, input$sub_sites)
    #picking the correct name of the data frame
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
    }

    # Create individual plots for each sub parameter
    plots <- list()
    if (length(input$sub_parameters) > 0) {
      all_sub_plot_data <- map_dfr(input$sub_parameters, function(param) {
      map_dfr(all_sub_sites, function(sub_site) {
        # Get the relevant sonde data
        relevant_sondes <- map(sub_site, ~ {
          sonde_name <- paste0(.x, "-", param)
          data_source <- NULL
          sonde_df <- NULL
          # Determine which directory to pull data from
          tryCatch({
            data_source <- retrieve_relevant_data_name(sonde_name, year_week)
            # cat("Data for",sonde_name,"will be pulled from",data_source,"\n")
          }, error = function(err) {
            #cat("Data for",sonde_name,"not found.\n")
            return(NULL)  # Return NULL if data source can't be determined
          })
          # Only try to pull in the data if data_source was successfully determined
          if (!is.null(data_source)) {
            tryCatch({
              sonde_df <- get(data_source)[[sonde_name]] %>%
                filter(DT_round >= week_min_day - days(2) & DT_round <= week_max_day + days(2))
            }, error = function(err) {
              #cat("Sonde", sonde_name, "not found.\n")
              return(NULL)  # Return NULL if sonde data can't be retrieved
            })
          }
          # Only return a list if both data_source and sonde_df are available
          if (!is.null(data_source) & !is.null(sonde_df)) {
            return(sonde_df)
          } else {
            return(NULL)  # Return NULL if either data_source or sonde_df is NULL
          }
        })
        # Remove any NULL results from the list
        relevant_sondes <- compact(relevant_sondes)
        if(length(relevant_sondes) == 0) {
          return(NULL)
        } else {
          relevant_sondes
        }
      })

    })%>%
      mutate(mean_verified = signif(mean_verified, digits = 3))


    max_param <- all_sub_plot_data %>%
      distinct(parameter, site) %>%  # ensure unique param-site pairs
      group_by(parameter) %>%
      summarise(site_count = n()) %>%
      ungroup()%>%
      #find the max and pull the parameter
      arrange(desc(site_count)) %>%
      slice(1) %>%
      pull(parameter)

    plots <- map(input$sub_parameters, function(param){

      #if(param == max_param){

      # Create plotly plot
      p <- plot_ly()

      # Add main site as grey points
      main_site_data <- all_sub_plot_data %>%
        filter(site == input$site,
               parameter == param)

      if (nrow(main_site_data) > 0) {
        main_site_data <- main_site_data %>%
          mutate(mean_plotting = case_when(
            final_status == "OMIT" ~ NA_real_,
            TRUE ~ mean
          ),
          final_status = case_when(
            is.na(final_status) ~ "PASS",
            TRUE ~ final_status))

          p <- p %>%
            add_markers(
              data = main_site_data,
              x = ~DT_round,
              y = ~mean_plotting,
              # Mapping color to the final_status column
              color = ~final_status,
              # Defining the specific hex or name colors for those levels
              colors = c( "PASS" = "gray", "FLAGGED" = "orange"),
              marker = list(size = 8),
              name = input$site,
              legendgroup = input$site,
              showlegend = (param == max_param)
            # to not have extra legends using plotly, one plot needs to have a legend, pick the one with the most sites (max_param)
            #otherwise it is the exact same as the other param plots
          )
      }

      # Add sub sites as colored lines
      sub_site_data <- all_sub_plot_data %>%
        filter(site %in% input$sub_sites,
               parameter == param)

      if (nrow(sub_site_data) > 0) {
        # Get site colors using joins
        sub_site_data_with_colors <- sub_site_data %>%
          left_join(site_color_combo, by = "site") %>%
          replace_na(list(color = "red"))%>%  # fallback color
          #Adding in interpolation where mean_verified is NA so that we can have more continous lines in the plotsq
            arrange(site, DT_round) %>%
            group_by(site)%>%
            mutate(
              mean_verified = zoo::na.approx(
                x = as.numeric(DT_round),
                object = mean_verified,
                maxgap = 3, # interpolate over 45 min
                na.rm = FALSE
              )
            ) %>%
            ungroup()

        # Add lines for each sub site using walk
        sub_site_data_with_colors %>%
          split(.$site) %>%
          iwalk(~ {
            site_color <- unique(.x$color)[1]  # Get the color for this site

            p <<- p %>%
              add_lines(
                data = .x,
                x = ~DT_round,
                y = ~mean_verified,
                line = list(color = site_color, width = 3),
                name = .y,
                legendgroup = .y,

                showlegend = param == max_param
                #showlegend = TRUE
              )
          })
      }

      # Configure layout
      p <- p %>%
        layout(
          xaxis = list(title = "Date"),
          yaxis = list(title = param),
          showlegend = param == max_param,
          legend = list(
            orientation = "h",
            x = 0.5,
            xanchor = "center",
            y = 1.1,
            title = list(text = "Sites")
          ),
          margin = list(t = 80, b = 40)
        )

      return(p)
    }) %>%
      compact() # remove any null plots
    }

    # Create USGS Flow Plot
    flow_df <- usgs_flow_data()

    if (!is.null(flow_df)) {
      p_flow <- plot_ly()

      # Add main site flow
      main_site_flow <- flow_df %>% filter(site == input$site)
      if (nrow(main_site_flow) > 0) {
        p_flow <- p_flow %>%
          add_lines(
            data = main_site_flow,
            x = ~datetime,
            y = ~meas_value,
            line = list(color = "black", width = 3),
            name = unique(main_site_flow$abbrev)[1],
            legendgroup = input$site,
            showlegend = FALSE
          )
      }

      p_flow <- p_flow %>%
        layout(
          xaxis = list(title = "Date"),
          yaxis = list(title = "USGS Flow (cfs)", type = "log"),
          margin = list(t = 80, b = 40)
        )

      plots <- c(plots, list(p_flow))
    }

    if (length(plots) > 0) {
      # Filter out any NULL or invalid plots using keep
      valid_plots <- plots %>%
        keep(~ !is.null(.x))

      if (length(valid_plots) == 1) {
        # If only one plot, return it directly with proper layout
        all_plot <- valid_plots[[1]] %>%
          layout(
            showlegend = TRUE,
            legend = list(
              orientation = "h",
              x = 0.5,
              xanchor = "center",
              y = 1.02,
              title = list(text = "Sites")
            ),
            margin = list(t = 100)
          )


      } else if (length(valid_plots) > 1) {
        # Create subplot with shared x-axis for multiple plots
        tryCatch({
          all_plot <- subplot(
            valid_plots,
            nrows = length(valid_plots),
            shareX = TRUE,
            titleY = TRUE,
            margin = 0.05
          ) %>%
            layout(
              showlegend = TRUE,
              legend = list(
                orientation = "h",
                xanchor = "center",
                title = list(text = "Sites")
              ),
              margin = list(t = 100)
            )

        }, error = function(e) {
          # Fallback: create a simple combined plot
          warning("Subplot creation failed, creating alternative layout")

          # Use the first valid plot as fallback
          all_plot <- valid_plots[[1]] %>%
            layout(
              showlegend = TRUE,
              legend = list(
                orientation = "h",
                x = 0.5,
                xanchor = "center",
                y = 1.02,
                title = list(text = "Sites")
              ),
              margin = list(t = 100),
              title = "Plot creation partially failed - showing first available plot"
            )
        })
      } else {
        # No valid plots - create empty plot
        all_plot <- plot_ly() %>%
          layout(
            title = "No data available",
            showlegend = FALSE
          )
      }

      all_plot
    }

  })

  output$field_notes_table <- DT::renderDataTable({
    req(current_week(), selected_data(), input$site)

    week_data <- isolate(selected_data()) %>% filter(week == current_week())
    week_min_day <- min(week_data$DT_round, na.rm = TRUE)
    week_max_day <- max(week_data$DT_round, na.rm = TRUE)

    # Use the existing meta folder path for this project
    if(year(week_min_day) <= 2023){
      # Use the existing meta folder path for this project
      field_notes_path <- here(meta_path, "field_notes_21-22.parquet")

    if (file.exists(field_notes_path)) {
      notes <- arrow::read_parquet(field_notes_path)

      # Try filtering by DT_round if it exists, otherwise by date or datetime
      if ("DT_round" %in% names(notes)) {
        notes <- notes %>% filter(DT_round >= week_min_day & DT_round <= week_max_day, site == input$site)
      } else if ("datetime" %in% names(notes)) {
        notes <- notes %>% filter(datetime >= week_min_day & datetime <= week_max_day, site == input$site)
      }

      # Rearrange columns so that visit comments is viewable
      notes <- notes %>%
        select(any_of(c("site", "date", "start_time_mst", "visit_comments", "visit_type")), everything())

      DT::datatable(notes, options = list(pageLength = 10, scrollX = TRUE))
    } else {
      data.frame(Message = paste("Field notes file not found at", field_notes_path))
    }

    } else{
      mWater_creds <- read_yaml(here("creds", "mWaterCreds.yml")) #this should be moved to global.R so that it is only read in once, but for now it is here
      #otherwise use mWater field notes
      mWater_data <- ross.wq.tools::load_mWater(creds = mWater_creds)
      # Only extract sensor visits for our site during the current week
      notes <- ross.wq.tools::grab_mWater_sensor_notes(mWater_api_data = mWater_data)%>%
        filter(DT_round >= week_min_day & DT_round <= week_max_day, site == input$site) %>% # DT round is in MST
        select(any_of(c("site", "date", "start_DT", "visit_comments", "visit_type")), everything())

      if(nrow(notes) > 0) {
        DT::datatable(notes, options = list(pageLength = 10, scrollX = TRUE))
      } else {
        data.frame(Message = "No field notes available for the selected week and site.")
      }
    }

  })

  #### Brush Tools ####

  # Create a reactive value to store multiple brush selections
  brushed_areas <- reactiveVal(list())

  # Add an observer to collect brush selections
  observeEvent(input$plot_brush, {
    req(input$plot_brush)

    # Get current week's data
    week_data <- selected_data() %>%
      filter(week == current_week())

    # Get brushed points
    brushed <- brushedPoints(week_data, input$plot_brush,
                             xvar = "DT_round", yvar = "mean")

    if (nrow(brushed) > 0) {
      # Store brush coordinates
      current_brush <- list(
        brush_dt_max = max(brushed$DT_round, na.rm = TRUE),
        brush_dt_min = min(brushed$DT_round, na.rm = TRUE),
        brush_mean_max = max(brushed$mean, na.rm = TRUE),
        brush_mean_min = min(brushed$mean, na.rm = TRUE)
      )
      #TODO: This creates duplicates for some reason but does work...
      existing_brushes <- brushed_areas()
      brushed_areas(c(existing_brushes, list(current_brush)))
    }

  })


  # Observer to clear brushed areas
  observeEvent(input$clear_brushes, {
    brushed_areas(list())
    session$resetBrush("plot_brush")

  })


  # Helper function to apply brush actions
  apply_brush_action <- function(action) {
    req(brushed_areas(), selected_data())

    user_brush_select <- action
    # Initialize updated data with current data
    updated_data <- selected_data()

    #deal with empty brush areas
    if(is_empty(brushed_areas())){
      # Update the data
      selected_data(updated_data)

      # Clear brushed areas after submission
      brushed_areas(list())
      session$resetBrush("plot_brush")
      #reset input$user_brush_flags to nothing
      updateSelectizeInput(session, "user_brush_flags", selected = "")

      showNotification("No Points Brushed, no changes applied", type = "message")

    }else{

      #go through brushed areas and apply user flags, brush omit or user id as needed
      updated_data <- reduce(
        brushed_areas()[seq(1, length(brushed_areas()), by = 2)], #skip everyother to reduce duplicates
        function(data, brush) {
          data %>%
            mutate(
              user_flag = case_when(

                between(DT_round, brush$brush_dt_min, brush$brush_dt_max) &
                  between(mean, brush$brush_mean_min, brush$brush_mean_max) &
                  user_brush_select == "A" ~ as.character(NA), # remove flags

                between(DT_round, brush$brush_dt_min, brush$brush_dt_max) &
                  between(mean, brush$brush_mean_min, brush$brush_mean_max) &
                  user_brush_select == "F" ~ paste(input$user_brush_flags, collapse = ";\n"), # add user selected flags

                between(DT_round, brush$brush_dt_min, brush$brush_dt_max) &
                  between(mean, brush$brush_mean_min, brush$brush_mean_max) &
                  user_brush_select == "O" ~ paste(input$user_brush_flags, collapse = ";\n"), # add user selected flags
                TRUE ~ user_flag
              ),
              brush_omit = case_when(
                between(DT_round, brush$brush_dt_min, brush$brush_dt_max) &
                  between(mean, brush$brush_mean_min, brush$brush_mean_max) &
                  user_brush_select == "O" ~ TRUE, # set omit to true if brushed with Omit
                between(DT_round, brush$brush_dt_min, brush$brush_dt_max) &
                  between(mean, brush$brush_mean_min, brush$brush_mean_max) &
                  user_brush_select %in% c("A", "F") ~ FALSE, # set omit to false if brushed with Accept or Flag
                TRUE ~ brush_omit
              ),
              user = ifelse(
                between(DT_round, brush$brush_dt_min, brush$brush_dt_max) &
                  between(mean, brush$brush_mean_min, brush$brush_mean_max),
                input$user, # set user to input$user if brushed
                user
              )
            )
        },
        .init = updated_data
      )

      # Update the data
      selected_data(updated_data)

      selected_data_cur_filename(update_intermediary_data(selected_data_cur_filename(), selected_data()))

      # Clear brushed areas after submission
      brushed_areas(list())
      session$resetBrush("plot_brush")
      #reset input$user_brush_flags to nothing
      updateSelectizeInput(session, "user_brush_flags", selected = "")

      showNotification("Brush Changes saved.", type = "message")
    }
  }

  observeEvent(input$btn_accept_brush, {
    apply_brush_action("A")
  })

  observeEvent(input$btn_flag_brush, {
    if (is.null(input$user_brush_flags)) {
      showNotification("Please select at least one flag.", type = "warning")
      return()
    }
    apply_brush_action("F")
  })

  observeEvent(input$btn_omit_brush, {
    if (is.null(input$user_brush_flags)) {
      showNotification("Please select at least one flag.", type = "warning")
      return()
    }
    apply_brush_action("O")
  })


  # Resets data to remove all brush inputs and weekly decisions
  observeEvent(input$reset_week, {

    req(selected_data(), current_week())

    updated_data <- selected_data() %>%
      filter(week == current_week())%>%
      mutate(user_flag = flag,
             brush_omit = FALSE,
             user = NA,
             final_status = NA,
             week_decision = NA,
             is_verified = NA)

    other_data <- selected_data() %>%
      filter(week != current_week())

    selected_data(rbind(updated_data, other_data))
    #update int file and save new filename to selected_data_cur_filename
    selected_data_cur_filename(update_intermediary_data(selected_data_cur_filename(), selected_data()))

  })

  #### Weekly Decision ####

  # UI for weekly decision radio buttons
  output$weekly_decision_radio <- renderUI({
    #filtering data to selected week
    week_data <-  selected_data()%>%
      filter(week == current_week())
    #create a button to show weeekly decision options
    radioButtons(
      "weekly_decision",
      label = "Make Weekly Decision:",
      choices = c("Accept ALL" = "aa",
                  "Accept Non Omit" = "ano",
                  "Keep Flags" = "kf",
                  "Omit Flagged" = "of",
                  "Omit ALL" = "oa",
                  "Skip" = "s"),

      selected = ifelse(all(is.na(week_data$week_decision)), "s", unique(week_data$week_decision)[1]), # If a week has an existing decision made, it will show up first here
      inline = TRUE
    )
  })

  # Submit decision button UI
  # Toggles on (green) if user has correctly selected a decision
  output$submit_decision_ui <- renderUI({
    req(input$weekly_decision)
    can_submit <- FALSE

    if (input$weekly_decision != "s") {
      can_submit <- T
    }

    actionButton(
      "submit_decision",
      "Submit Weekly Decision",
      class = ifelse(can_submit, "btn-success", "btn-secondary"),
      disabled = !can_submit
    )
  })

  # Update selected_data() on backend with submitted decision
  observeEvent(input$submit_decision, {
    req(input$weekly_decision != "s", selected_data())
    #update backend data

    #weekly_decision <- input$weekly_decision
    updated_week_data <- selected_data() %>%
      filter(week == current_week())%>%
      mutate(
        final_status = case_when(
          #AA:Pass all data
          input$weekly_decision  == "aa"  ~ "PASS",
          #ANO: Accept Non Omit
          input$weekly_decision == "ano" & !brush_omit ~ "PASS", # pass data that is not user select omit
          #KF: Keep FLagged, retain flag into final data (sus but on the edge)
          input$weekly_decision == "kf" & is.na(user_flag) & !brush_omit ~ "PASS", # pass data that is not user select omit
          input$weekly_decision == "kf" & !is.na(user_flag) & !brush_omit ~ "FLAGGED", # tag data that is flagged
          #OF: Omit Flagged
          input$weekly_decision == "of" & is.na(user_flag) & !brush_omit ~ "PASS", # pass data that is not user select omit
          input$weekly_decision == "of" & !is.na(user_flag) & !brush_omit ~ "OMIT", # omit data that is flagged
          #OA: Omit All
          input$weekly_decision == "oa"  ~ "OMIT",
          # Omit any user selected omit data (assuming AA was not the choice)
          input$weekly_decision != "aa" & brush_omit ~ "OMIT"),
        # update brush omit to T if needed (OA, OF) and to F as needed
        brush_omit = case_when(
          # AA: Pass All (remove any omits)
          input$weekly_decision  == "aa"  ~ FALSE,
          #OA: Omit All
          input$weekly_decision == "oa"  ~ TRUE,
          #OF: Omit Flagged
          input$weekly_decision == "of" & !is.na(user_flag) & !brush_omit ~ TRUE,
          #Otherwise keep the original brushed decision
          TRUE ~ brush_omit
        ),
        user_flag = case_when(
          #removing all flags if user selects accept all
          input$weekly_decision == "aa" ~ NA,
          TRUE ~ user_flag),
        user = input$user,
        week_decision = input$weekly_decision,
        is_verified = TRUE,
        #TODO: verification status seems to have a different role in the previous version
        verification_status = final_status,
        #Omit flagged and omitted data from mean_verified
        mean_verified = case_when(
          final_status %in%c("PASS", "FLAG") ~ mean,
          TRUE ~ NA_real_)
      )

    other_data <- selected_data() %>%
      filter(week != current_week())

    selected_data(bind_rows(other_data, updated_week_data)%>%arrange(DT_round))

    #update int file and save new filename to selected_data_cur_filename
    selected_data_cur_filename(update_intermediary_data(selected_data_cur_filename(), selected_data()))

    # Get all weeks and current week
    weeks <- unique(selected_data()$week)
    current <- current_week()
    idx <- which(weeks == current)


    # Move to next week if available
    if(all(!is.na(selected_data()$final_status))){

      showNotification("All weeks have been reviewed.", type = "message")
      updateTabsetPanel(session, inputId = "tabs", selected = "Finalize Data")

    }else{

      if(idx == length(weeks)){
        # Find min week where is_verified is NA
        week_min <- selected_data()%>%filter(is.na(is_verified))%>%pull(week)%>%min()
        idx <- which(weeks == week_min)
        current_week(weeks[idx])

        showNotification(
          "No more weeks to verify. Moving to earliest unverified week",
          type = "warning")

      }else{ #move to next week, if this is verified, user should probably just move to final data tab to see what remains
        current_week(weeks[idx + 1])
      }

    }

    # Show notification of submission
    showNotification(
      paste("Decision", toupper(input$weekly_decision), "submitted"),
      type = "message"
    )
    # Reset weekly decision back to "s"
    updateRadioButtons(session, "weekly_decision", selected = "s")
    updateCheckboxInput(session, "remove_omit", value = FALSE)


  })


  ##### Final Verification Tab ####
  # get weeks for final week selection select input
  observe({
    req(selected_data())
    weeks <- selected_data() %>%
      pull(week) %>%
      unique() %>%
      sort()

    updateSelectInput(session, "final_week_selection",
                      choices = weeks)
  })

  # Handle week selection in final tab
  observeEvent(input$goto_final_week, {
    req(input$final_week_selection)
    selected_week <- as.numeric(input$final_week_selection)
    current_week(selected_week)
    updateNavbarPage(session, inputId = "tabs", selected = "Data Verification")
  })

  # Add plotly plot for final overview
  output$final_plot <- renderPlotly({
    req(selected_data())

    final_plot_data <- selected_data()

    if (input$remove_omit_finalplot) {

      final_plot_data <- final_plot_data %>%
        filter(final_status != "OMIT"|is.na(final_status))
    }


    #get year from data, should only be one year but just in case we have multiple years we will take the minimum
    year <- min(selected_data()$year)
    # Find the first day of the minimum week in the data
    start_date <- as.POSIXct(min(selected_data()$DT_round, na.rm = T))
    # Find the last day of the maximum week in the data
    end_date <- as.POSIXct(max(selected_data()$DT_round, na.rm = T))
    # Create vertical lines at the beginning of each week


    vline_dates <- seq(as.POSIXct(paste0(year, "-01-01")),
                       as.POSIXct(paste0(year, "-12-31")),
                       by = "week")%>%
      #filter to the start and end date of the data
      keep(~ .x >= start_date & .x <= end_date)


    #add 3 days to each vertical line to center the week
    week_dates <- vline_dates + days(3)

    final_status_colors <- c("PASS" = "#008a18",
                             "OMIT" = "#ff1100",
                             "FLAGGED" = "#ff8200",
                             "NA" = "grey")

    # Replace NA values in final_status (this will not affect the actual saved data, just for plotting)
    final_plot_data$final_status <- ifelse(is.na(final_plot_data$final_status), "NA", final_plot_data$final_status)
    # This seems to be important for the plotly to work
    final_plot_data$final_status <- as.factor(final_plot_data$final_status)


    # Create plotly object
    p_plotly <- plot_ly(
      data = final_plot_data,
      x = ~DT_round,
      y = ~mean,
      type = 'scatter',
      mode = 'markers',
      color = ~final_status,
      colors = final_status_colors,
      text = ~paste0("Week ", week, "\nStatus: ", final_status),
      hoverinfo = "text"
    ) %>%
      layout(
        title = list(
          text = paste0("Complete Dataset Overview: ", input$site, "-", input$parameter),
          x = 0.5 # Centers the title
        ),
        xaxis = list(
          title = "Date",
          tickformat = "%b %d",
          tickmode = "array",
          tickvals = week_dates,
          ticktext = paste0("Week ", week(week_dates), "\n", format(week_dates, "%b %d")),
          showgrid = TRUE,
          domain = c(0, 1)  # Ensure it spans the full width
        ),
        yaxis = list(
          title = input$parameter
        ),
        xaxis2 = list(
          title = "Week #",
          tickmode = "array",
          tickvals = week_dates,
          ticktext = as.character(week(week_dates)),
          overlaying = "x",
          side = "top",
          showgrid = FALSE,
          zeroline = FALSE,
          ticks = "outside",
          ticklen = 5
        ),
        shapes = lapply(vline_dates, function(date) {
          list(
            type = "line",
            x0 = date, x1 = date, y0 = 0, y1 = 1,
            xref = "x", yref = "paper",
            line = list(color = "black", width = 1)
          )
        })
      )
    if (input$remove_omit_finalplot) {
      p_plotly <- p_plotly %>%
        layout(
          annotations = list(
            list(
              text = "Omitted data removed",
              x = 0.5,
              y = 1.02,
              xref = "paper",
              yref = "paper",
              showarrow = FALSE,
              font = list(size = 12, color = "gray50")
            )
          )
        )
    }

    if(input$log10_finalplot){

      p_plotly <- layout(p_plotly, yaxis = list(type = "log"))
    }


    p_plotly

  })

  output$submit_final_button <- renderUI({
    # Check if all selected data is verified
    all_verified <- FALSE
    if (!is.null(selected_data())) {
      all_verified <- all(!is.na(selected_data()$final_status))
    }

    # Create the button, disabled if not all data is verified
    if(all_verified) {
      actionButton("submit_final", "Submit Finalized Dataset",
                   class = "btn-success w-100")
    } else {
      actionButton("submit_final", "Submit Finalized Dataset",
                   class = "btn-success w-100 disabled",
                   disabled = TRUE)
    }
  })


  # Handle final submission
  observeEvent(input$submit_final, {

    #set is_finalized in selected_data() to true
    update_finalized <- selected_data()%>%
      mutate(is_finalized = TRUE)
    # Move the dataset to the finalized directory and print the file name for users to see
    final_name <- move_file_to_verified_directory(int_to_fin_filename = selected_data_cur_filename(), int_to_fin_df = update_finalized)

    showNotification(paste0(input$site, "-", input$parameter," finalized and saved to ", final_name ), type = "message")
    updateNavbarPage(session, inputId = "tabs", selected = "Data Selection")
    #reload the session to update the data displayed in the data selection tab and reset all reactive elements
    session$reload()

  })

  #### Extras ####
  # Handle quit button
  observeEvent(input$quit_app, {
    #update int file and save new filename to selected_data_cur_filename
    selected_data_cur_filename(update_intermediary_data(selected_data_cur_filename(), selected_data()))
    stopApp()
  })

  # Add this to handle the keyboard shortcut
  observeEvent(input$q_key, {
    if (input$q_key == "q") {
      #update int file and save new filename to selected_data_cur_filename
      selected_data_cur_filename(update_intermediary_data(selected_data_cur_filename(), selected_data()))
      stopApp()
    }
  })

}

