#UI Script
ui <- page_navbar(
  title = "Data Processing Pipeline",
  id = "tabs",
  nav_item(
    input_dark_mode(id = "dark_mode", mode = "light") #Toggle light vs dark mode
  ),
  theme = bs_theme(preset = "bootstrap"),



  #### Tab 1: Data Selection ####
  nav_panel(
    title = "Data Selection",
    card(
      card_header("Select Your Data"),
      card_body(
        uiOutput("conditional_data_ui")
      )
    )
  ),


  #### Tab 2: Data Verification ####
  nav_panel(
    title = "Data Verification",
    layout_sidebar(
      sidebar = sidebar(
        width = 350,

        # Navigation
        div(
          class = "d-flex justify-content-between align-items-center mb-3",
          actionButton("prev_week", "← Prev", class = "btn-secondary btn-sm"),
          actionButton("reset_week", "Reset", class = "btn-danger btn-sm"),
          actionButton("next_week", "Next →", class = "btn-secondary btn-sm")
        ),
        div(
          class = "d-flex justify-content-between align-items-center mb-3",
          keys::useKeys(),
          keys::keysInput("q_key", "q"),
          actionButton("quit_app", "Quit", class = "btn-danger btn-sm w-100")
        ),

        # Weekly Decision
        card(
          card_header(h6("Weekly Decision", class = "m-0")),
          card_body(
            class = "p-2",
            uiOutput("weekly_decision_radio"),
            uiOutput("submit_decision_ui")
          )
        ),

        # Data Brush Tools
        card(
          card_header(h6("Brush Actions", class = "m-0")),
          card_body(
            class = "p-2",
            selectizeInput("user_brush_flags", "Select Flag(s):",
                           choices = available_flags,
                           multiple = TRUE,
                           options = list(plugins = "remove_button")),
            div(
              class = "d-flex flex-column gap-2",
              actionButton("btn_accept_brush", "Accept Selection", class = "btn-success btn-sm"),
              actionButton("btn_flag_brush", "Flag Selection", class = "btn-warning btn-sm"),
              actionButton("btn_omit_brush", "Omit Selection", class = "btn-danger btn-sm"),
              actionButton("clear_brushes", "Clear Brush", class = "btn-secondary btn-sm mt-2")
            )
          )
        ),

        # Additional Sites
        card(
          card_header(h6("Additional Sites", class = "m-0")),
          card_body(
            class = "p-2",
            selectizeInput("add_sites", label = NULL,
                           choices = available_sites,
                           multiple = TRUE,
                           options = list(plugins = "remove_button"),
                           width = "100%")
          )
        ),

        # Plot Options
        card(
          card_header(h6("Plot Options", class = "m-0")),
          card_body(
            class = "p-2",
            checkboxGroupInput("plot_options", label = NULL,
                               choices = c("Remove Omit" = "remove_omit",
                                           "Remove Flag" = "remove_flag",
                                           "Plot Line" = "add_line",
                                           "Thresholds" = "incl_thresholds",
                                           "Log 10" = "plot_log10",
                                           "Extra Data" = "incl_ex_days",
                                           "Show Legend" = "show_legend"),
                               selected = c("incl_ex_days", "show_legend",
                                            "add_line", "remove_omit"))
          )
        )
      ), # end sidebar

      # Main layout for plots
      layout_columns(
        col_widths = c(8, 4),

        # Main plot card
        card(
          full_screen = TRUE,
          card_body(
            plotOutput("main_plot",
                       height = "100%",
                       brush = brushOpts(
                         id = "plot_brush",
                         resetOnNew = TRUE  # Simplified: immediate actions
                       ))
          )
        ),

        # Sub plots card
        navset_card_tab(
          id = "additional_tabs",
          nav_panel("Plots",
            layout_columns(
              col_widths = c(6, 6),
              selectizeInput("sub_parameters", "Select Parameters:",
                             choices = available_parameters,
                             multiple = TRUE,
                             options = list(plugins = "remove_button")),
              selectizeInput("sub_sites", "Select Sites:",
                             choices = available_sites,
                             multiple = TRUE,
                             options = list(plugins = "remove_button"))
            ),
            div(
              class = "flex-fill d-flex flex-column",
              style = "overflow-y: auto; min-height: 600px;",
              plotlyOutput("sub_plots", width = "100%", height = "100%")
            )
          ),
          nav_panel("Field Notes",
            div(
              class = "flex-fill d-flex flex-column",
              style = "overflow-y: auto; min-height: 600px;",
              DT::dataTableOutput("field_notes_table")
            )
          )
        )
      )
    )
  ),

  #### Tab 3: Final Data View ####

  nav_panel(
    title = "Finalize Data",
    layout_columns(
      col_widths = c(10, 2),

      # Main Plot Card
      #To Do: Convert to plotly object for better data vis
      card(
        card_body(
          plotlyOutput("final_plot", height = "100%", width = "100%")
        )
      ),

      # Week Selection and Actions Card
      card(
        card_header("Modify Verification"),
        card_body(
          materialSwitch(
            inputId = "remove_omit_finalplot",
            label = "Remove omitted data from plot",
            value = FALSE,
            width = "200px",
            status = "success"
          ),
          materialSwitch(
            inputId = "log10_finalplot",
            label = "Log Transform",
            value = FALSE,
            width = "200px",
            status = "success"
          ),
          selectInput("final_week_selection", "Select Week:", choices = NULL),
          actionButton("goto_final_week", "Return to Selected Week",
                       class = "btn-primary w-100 mb-3"),
          hr(),
          uiOutput("submit_final_button") # Replaced the direct button with a dynamic UI output
        )
      )
    )
  )

)
