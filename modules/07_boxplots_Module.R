# 07_Boxplot_Module.R
# Author: Ghazal Sinjar
# Date: 23.06.2025
# Shiny module for generating CpG beta value boxplots for DMRs





# ─────────────────────────────────────
# User Interface (UI) FUNCTION
# ─────────────────────────────────────
boxplotUI <- function(id) {
  ns <- NS(id)
  
  layout_sidebar(
    sidebar = sidebar(
      selectInput(ns("dmr"), "Select DMR ID", choices = NULL),
      switchInput(ns("interactive"), "Interactive Boxplot", value = FALSE),
      switchInput(ns("pos_spacing"), "Use Positional Spacing", value = FALSE),
      helpText("Positional spacing will arrange the boxes according to their position on the chromosome, whereas PositionLabel will just use the CpG names."),
      actionButton(ns("create_plot"), "Create Boxplot"),
      br(),
      downloadButton(ns("download_plot"), "Download Plot")
    ),
    
    div(
      style = "padding-left: 15px; padding-right: 15px;",
      layout_columns(
        col_widths = c(12),
        card(
          card_header("Boxplot Creating Status"),
          card_body(verbatimTextOutput(ns("boxplot_status")))
        ),
        uiOutput(ns("plot_area"))
      )
    )
  )
}



# ─────────────────────────────────────
# SERVER FUNCTION
# ─────────────────────────────────────
boxplotServer <- function(id, dmrs_table, annotated_with_betas_df, pheno_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Update DMR dropdown
    observe({
      updateSelectInput(session, "dmr", choices = sort(unique(dmrs_table$DMR_ID)))
    })
    
    status <- reactiveVal("ℹ️ Waiting to select a DMR ID.")
    output$boxplot_status <- renderText({ status() })
    
    observeEvent(input$create_plot, {
      req(input$dmr)
      
      DMRx <- input$dmr
      interactive_flag <- input$interactive
      use_spacing <- input$pos_spacing
      
      status_lines <- c("🔁 Step 1: Extracting CpGs in the DMR...")
      
      region_cpgs <- tryCatch(
        {
          extract_cpgs_in_DMR(DMRx, dmrs_table, annotated_with_betas_df, pheno_data)
        },
        error = function(e) {
          status_lines <<- c(status_lines, paste0("❌ Error in Step 1: ", e$message))
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
        status_lines <- c(status_lines, "✅ CpG table created.")
      } else {
        status_lines <- c(status_lines, paste0("❌ ", region_cpgs))
        status(paste(status_lines, collapse = "\n"))
        output$plot_area <- renderUI({ NULL })
        return(NULL)
      }
      
      status_lines <- c(status_lines, "🔁 Step 2: Reshaping to long format...")
      
      long_table <- tryCatch(
        {
          reshape_to_long_beta(region_cpgs, pheno_data)
        },
        error = function(e) {
          status_lines <<- c(status_lines, paste0("❌ Error in Step 2: ", e$message))
          status(paste(status_lines, collapse = "\n"))
          output$plot_area <- renderUI({ NULL })
          return(NULL)
        }
      )
      
      if (is.null(long_table)) {
        output$plot_area <- renderUI({ NULL })
        return(NULL)
      }
      
      status_lines <- c(status_lines, "✅ Table reshaped.")
      
      status_lines <- c(status_lines, "🔁 Step 3: Creating boxplot...")
      
      plot_result <- tryCatch(
        {
          create_boxplot(long_table, interactive = interactive_flag, use_positional_spacing = use_spacing)
        },
        error = function(e) {
          status_lines <<- c(status_lines, paste0("❌ Error in Step 3: ", e$message))
          status(paste(status_lines, collapse = "\n"))
          output$plot_area <- renderUI({ NULL })
          return(NULL)
        }
      )
      
      if (is.null(plot_result)) {
        output$plot_area <- renderUI({ NULL })
        return(NULL)
      }
      
      status_lines <- c(status_lines, "✅ Plot successfully created.")
      
      output$plot_area <- renderUI({
        if (interactive_flag && inherits(plot_result, "plotly")) {
          plotlyOutput(ns("interactive_boxplot"), width = "100%", height = "80vh")
        } else if (!interactive_flag && inherits(plot_result, "ggplot")) {
          plotOutput(ns("boxplot_output"), width = "100%", height = "80vh")
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
      
      status(paste(status_lines, collapse = "\n"))
    })
  })
}


'# app.R

# Load required libs
library(shiny)
library(shinyWidgets)
library(bslib)
library(plotly)
library(ggplot2)

# Source utility and module files
source("../utils/dmrs_boxplot_utils.R")
#source("./modules/boxplot_module.R")

# Load data
#dmr_results <- readRDS("../modules/intermediate_data/DMRs_cutoff_neg0.15_to_0.15_B0_2025-06-23.rds")
dmrs_table <- dmr_results$dmr_table
#annotation_results <- readRDS("../modules/intermediate_data/annotated_object_20250623.rds")
#annotated_with_betas_df <- annotation_results$annotated_table
pheno_data <- dmr_results$pheno_data

# App UI
ui <- page_navbar(
  title = "EPIC Array Pipeline",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  nav_panel("Boxplots", boxplotUI("boxplot"))
)

# App server
server <- function(input, output, session) {
  boxplotServer("boxplot", dmrs_table, annotated_with_betas_df, pheno_data)
}

# Run app
shinyApp(ui, server)
'