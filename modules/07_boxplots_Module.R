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
      selectInput(ns("dmr"), "Select DMR ID (chr_gene)", choices = NULL),
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
boxplotServer <- function(id, dmr_output_reactive, annotation_output_reactive) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    dmr_data_container <- reactive({
      req(dmr_output_reactive())
      dmr_output_reactive() 
    })
    
    annotation_data_container <- reactive({
      req(annotation_output_reactive())
      annotation_output_reactive() 
    })
    
    
    dmrs_table_r <- reactive({
      req(dmr_data_container())
      dmr_data_container()$dmr_table() 
    })
    
    pheno_data_r <- reactive({
      req(dmr_data_container())
      dmr_data_container()$pheno() 
    })
    
    annotated_with_betas_df_r <- reactive({
      req(annotation_data_container())
      annotation_data_container()$annotated_table()
    })

    'observe({
      req(dmrs_table_r())
      updateSelectInput(session, "dmr", choices = sort(unique(dmrs_table_r()$DMR_ID)))
    })'
    observe({
      req(dmrs_table_r()) # Ensure dmrs_table_r() (the reactive value of the DMRs table) is available
      
      # Get the DMRs data frame
      dmr_data <- dmrs_table_r()
      
      # Create a named vector for choices
      # The names will be what the user sees, and the values will be the actual DMR_ID to be used internally by your app.
      choices_list <- setNames(
        dmr_data$DMR_ID, # The actual value passed to input$dmr
        # Updated to include the chromosome
        paste0(dmr_data$DMR_ID, " (", dmr_data$chr, "_", dmr_data$overlapped_gene_name, ")") 
      )
      
      updateSelectInput(session, "dmr", choices = choices_list) 
    })
    
    status <- reactiveVal("ℹ️ Waiting to select a DMR ID.")
    output$boxplot_status <- renderText({ status() })
    
    observeEvent(input$create_plot, {
      # All req() calls here are correct as they are already calling the final reactives
      req(input$dmr, dmrs_table_r(), annotated_with_betas_df_r(), pheno_data_r())
      
      DMRx <- input$dmr
      interactive_flag <- input$interactive
      use_spacing <- input$pos_spacing
      
      # Start progress tracking
      withProgress(message = 'Creating boxplot', value = 0, {
        # Step 1: Extracting CpGs
        incProgress(0.1, detail = "Extracting CpGs in the DMR...")
        status_lines <- c("🔁 Step 1: Extracting CpGs in the DMR...")
        
        region_cpgs <- tryCatch(
          {
            # Pass the values (actual data frames) to the utility function
            extract_cpgs_in_DMR(DMRx, dmrs_table_r(), annotated_with_betas_df_r(), pheno_data_r())
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
          return()
        }
        
        if (is.data.frame(region_cpgs)) {
          status_lines <- c(status_lines, "✅ CpG table created.")
          incProgress(0.3, detail = "CpGs extracted successfully")
        } else {
          status_lines <- c(status_lines, paste0("❌ ", region_cpgs))
          status(paste(status_lines, collapse = "\n"))
          output$plot_area <- renderUI({ NULL })
          return()
        }
        
        # Step 2: Reshaping data
        incProgress(0.1, detail = "Reshaping to long format...")
        status_lines <- c(status_lines, "🔁 Step 2: Reshaping to long format...")
        
        long_table <- tryCatch(
          {
            reshape_to_long_beta(region_cpgs, pheno_data_r()) # Pass the value
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
          return()
        }
        
        status_lines <- c(status_lines, "✅ Table reshaped.")
        incProgress(0.3, detail = "Data reshaped successfully")
        
        # Step 3: Creating plot
        incProgress(0.1, detail = "Creating boxplot...")
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
          return()
        }
        
        status_lines <- c(status_lines, "✅ Plot successfully created.")
        incProgress(0.2, detail = "Boxplot created")
        
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
              if (inherits(plot_result, "ggplot")) {
                htmlwidgets::saveWidget(ggplotly(plot_result), file)
              } else {
                htmlwidgets::saveWidget(plot_result, file)
              }
            } else {
              ggsave(file, plot_result)
            }
          }
        )
        
        status(paste(status_lines, collapse = "\n"))
      }) # End of withProgress
    })
  })
}


'#test app
# Load libraries
library(shiny)
library(bslib)
library(shinyWidgets)
library(plotly)
library(ggplot2)

# Load module and utilities
source("../utils/dmrs_boxplot_utils.R")

# Load RDS objects (as static lists)
results_dmr <- readRDS("./intermediate_data/DMRs_cutoff_neg0.15_to_0.15_B0_2025-07-05.rds")
results_anno <- readRDS("./intermediate_data/annotated_object_20250623.rds")

# Wrap each element in reactiveVal(), and return a list of reactive functions
dmr_data <- reactive({
  list(
    dmr_table = reactiveVal(results_dmr$dmr_table),
    pheno = reactiveVal(results_dmr$pheno)
  )
})

anno_data <- reactive({
  list(
    annotated_table = reactiveVal(results_anno$annotated_table)
  )
})

# Create wrapper reactivity so that inside the module,
# we can do dmr_output_reactive()$dmr_table()
dmr_output_reactive <- reactive({
  list(
    dmr_table = dmr_data()$dmr_table,
    pheno = dmr_data()$pheno
  )
})

annotation_output_reactive <- reactive({
  list(
    annotated_table = anno_data()$annotated_table
  )
})

# ---- UI ----
ui <- page_navbar(
  title = "DMR Boxplot Test App",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  nav_panel("Boxplots", boxplotUI("boxplot"))
)

# ---- Server ----
server <- function(input, output, session) {
  boxplotServer(
    id = "boxplot",
    dmr_output_reactive = dmr_output_reactive,
    annotation_output_reactive = annotation_output_reactive
  )
}

# ---- Run the app ----
shinyApp(ui = ui, server = server)
'