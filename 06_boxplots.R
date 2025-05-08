#library(shiny)
#library(shinyWidgets)  # This is required for switchInput
#library(bslib)
#library(htmlwidgets)

# inputs:
#dmrs_tbl
#annotated_with_betas
#pheno_data <- pData(ratio_geno_Swan_NoSNP)

# functions:
# 1- region_cpgs        <-  extract_cpgs_in_DMR(DMRx, DMR_table, annotation_with_betas, pheno_data)
# 2- long_format_table <-  reshape_to_long_beta(region_cpgs, pheno_data)
# 3- create_boxplot(long_format_table, interactive = FALSE, use_positional_spacing = FALSE)


ui_boxplot <- page_navbar(
  title = "EPIC Array Pipeline",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  
  nav_panel(
    "Boxplots",
    layout_sidebar(
      sidebar = sidebar(
        textInput("dmr", "Enter DMR Name", placeholder = "DMR1"),  # Text input for DMR
        switchInput("interactive", "Interactive Boxplot", value = FALSE),  # Switch button for interactive plot
        switchInput("pos_spacing", "Use Positional Spacing", value = FALSE),  # Switch button for positional spacing
        helpText("Positional spacing will arrange the boxes according to their position on the chromosome, whereas PositionLabel will just use the CpG names."),
        actionButton("create_plot", "Create Boxplot"),
        br(),
        downloadButton("download_plot", "Download Plot")  
      ),
      
      mainPanel(
        # Status card
        card(
          card_header("Boxplot Creating Status"),
          card_body(
            verbatimTextOutput("boxplot_status")  # Output area for status messages
          )
        ),
        uiOutput("plot_area")
      )
    )
  )
)


server_boxplots <- function(input, output, session) {
  
  # Initial status message
  status <- reactiveVal("ℹ️ Waiting to enter DMR number.")
  output$boxplot_status <- renderText({ status() })
  
  observeEvent(input$create_plot, {
    req(input$dmr)
    
    DMRx <- input$dmr
    interactive_flag <- input$interactive
    use_spacing <- input$pos_spacing
    
    # Step 1: Extract CpGs in DMR
    status_lines <- c("🔁 Step 1: Extracting table of CpGs in the DMR")
    
    region_cpgs <- tryCatch(
      {
        extract_cpgs_in_DMR(DMRx, dmrs_df, annotated_with_betas_df, pheno_data)
      },
      error = function(e) {
        error_msg <- paste0("❌ Error in Step 1: ", e$message)
        status_lines <<- c(status_lines, error_msg)
        status(paste(status_lines, collapse = "\n"))
        output$plot_area <- renderUI({ NULL })
        return(NULL)
      }
    )
    
    if (is.null(region_cpgs)) {
      output$plot_area <- renderUI({ NULL })
      return(NULL)
    }
    
    if (is.data.frame(region_cpgs)) {
      status_lines <- c(status_lines, "✅ Boxplot table of CpGs is created...")
    } else {
      status_lines <- c(status_lines, paste0("❌ ", region_cpgs))
      status(paste(status_lines, collapse = "\n"))
      output$plot_area <- renderUI({ NULL })
      return(NULL)
    }
    
    # Step 2: Convert table to long format
    status_lines <- c(status_lines, "🔁 Step 2: Converting the table to long format...")
    
    long_table <- tryCatch(
      {
        reshape_to_long_beta(region_cpgs, pheno_data)
      },
      error = function(e) {
        error_msg <- paste0("❌ Error in Step 2: ", e$message)
        status_lines <<- c(status_lines, error_msg)
        status(paste(status_lines, collapse = "\n"))
        output$plot_area <- renderUI({ NULL })
        return(NULL)
      }
    )
    
    if (is.null(long_table)) {
      output$plot_area <- renderUI({ NULL })
      return(NULL)
    }
    
    if (is.data.frame(long_table)) {
      status_lines <- c(status_lines, "✅ Table is converted.")
    } else {
      status_lines <- c(status_lines, paste0("❌ ", long_table))
      status(paste(status_lines, collapse = "\n"))
      output$plot_area <- renderUI({ NULL })
      return(NULL)
    }
    
    # Step 3: Create the boxplot
    status_lines <- c(status_lines, "🔁 Step 3: Creating the boxplot...")
    
    plot_result <- tryCatch(
      {
        create_boxplot(long_table, interactive = interactive_flag, use_positional_spacing = use_spacing)
      },
      error = function(e) {
        error_msg <- paste0("❌ Error in Step 3: ", e$message)
        status_lines <<- c(status_lines, error_msg)
        status(paste(status_lines, collapse = "\n"))
        output$plot_area <- renderUI({ NULL })
        return(NULL)
      }, 
      warning = function(w) {
        warn_msg <- paste0("⚠️ Warning in Step 3: ", w$message)
        status_lines <<- c(status_lines, warn_msg)
        status(paste(status_lines, collapse = "\n"))
      }
    )
    
    if (is.null(plot_result)) {
      output$plot_area <- renderUI({ NULL })
      return(NULL)
    }
    
    if (inherits(plot_result, "ggplot") || inherits(plot_result, "plotly")) {
      status_lines <- c(status_lines, "✅ Plot successfully created.")
    } else {
      status_lines <- c(status_lines, paste0("❌ ", plot_result))
      status(paste(status_lines, collapse = "\n"))
      output$plot_area <- renderUI({ NULL })
      return(NULL)
    }
    
    output$plot_area <- renderUI({
      if (interactive_flag && inherits(plot_result, "plotly")) {
        plotlyOutput("interactive_boxplot", width = "100%", height = "80vh")
      } else if (!interactive_flag && inherits(plot_result, "ggplot")) {
        plotOutput("boxplot_output", width = "100%", height = "80vh")
      } else {
        NULL
      }
    })
    
    output$boxplot_output <- renderPlot({
      if (!interactive_flag && inherits(plot_result, "ggplot")) plot_result else NULL
    })
    
    output$interactive_boxplot <- renderPlotly({
      if (interactive_flag && inherits(plot_result, "plotly")) plot_result else NULL
    })
    
    
    # Download setup
    output$download_plot <- downloadHandler(
      filename = function() {
        if (interactive_flag) paste0(DMRx, "_plot.html") else paste0(DMRx, "_plot.pdf")
      },
      content = function(file) {
        if (interactive_flag) {
          htmlwidgets::saveWidget(ggplotly(plot_result), file)
        } else {
          ggsave(file, plot_result)
        }
      }
    )
    
    # Final status update
    status(paste(status_lines, collapse = "\n"))
  })
}

  
# Run the UI (server logic will need to be added separately)
shinyApp(ui = ui_boxplot, server_boxplots)

