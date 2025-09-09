# dmp_module.R

# Load necessary libraries for the module
library(shiny)
library(bslib)
library(DT)
library(minfi)
library(fs)
library(shinyjs)

# ─────────────────────────────────────
# User Interface (UI) FUNCTION
# ─────────────────────────────────────
dmp_UI <- function(id) {
  ns <- NS(id)
  
  tagList(
    shinyjs::useShinyjs(),
    
    page_sidebar(
      sidebar = sidebar(
        width = "300px",
        
        # q-value cutoff numeric input
        numericInput(ns("qvalue_cutoff"), "q-value Cutoff (FDR):",
                     value = 1, min = 0, max = 1, step = 0.00001),
        helpText("It is recommended not to apply a q-value cutoff (i.e., set it to 1). This ensures that p-values for all CpGs are retained in the output, allowing you to filter the table afterward based on your own criteria."),
        br(),
        # Reference group display
        # Reference group display
        div(
          class = "mb-2",  # smaller margin below
          HTML('<label class="control-label fw-semibold fs-7">Reference Group:</label>'),
          textOutput(ns("display_ref_group"))
        ),
        
        helpText(
          "Methylation differences calculated:",
          list(
            div(strong("-"), "Hypomethylation relative to reference"),
            div(strong("+"), "Hypermethylation relative to reference")
          )
        ),
        
        # Action button
        actionButton(ns("run_dmp_finder"), "Run DMP Finder", 
                     class = "btn-primary"),
        
        # Toggle notes section
        div(
          style = "margin-top: 20px;",
          actionLink(ns("toggle_notes"), "📘 Show interpretation guide"),
          shinyjs::hidden(
            div(
              id = ns("note_section"),
              h5(strong("Column Descriptions:")),
              list(
                div(strong("intercept:"), "Reference group baseline methylation"),
                div(strong("f:"), "Between-group vs within-group variability"),
                div(strong("pval:"), "Unadjusted probability of chance difference"),
                div(strong("qval:"), "FDR-adjusted significance value"),
                div(strong("M_mean_difference:"), "Logit-transformed difference"),
                div(strong("Beta_mean_difference:"), "Log2 ratio of proportions")
              )
            )
          ),
          hr()
        ),
        
        # Download UI
        uiOutput(ns("download_ui"))
      ),
      
      # The main content is now passed directly as the second argument to page_sidebar
      div(
        style = "padding-left: 15px; padding-right: 15px;",
        layout_columns(
          col_widths = c(12),
          fill = TRUE,
          card(
            card_header("DMP Identification Status"),
            htmlOutput(ns("status_message"))  # Instead of textOutput 
          )
        ),
        
        # Bottom section: DMP Table
        div(
          style = "margin-top: 20px;",
          h4("Detected DMPs"),
          helpText("The table below shows all CpG sites identified as differentially methylated positions (DMPs) based on your selected q-value cutoff."),
          br(),
          DT::dataTableOutput(ns("dmp_table"))
        )
      )
    )
  )
}
# ─────────────────────────────────────
# SERVER FUNCTION 
# ─────────────────────────────────────
dmp_Server <- function(id,filtered_data , pheno, ref_group, project_output_dir) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Reactive value to store the combined M and Beta mean table
    combined_m_beta_tbl <- reactiveVal(NULL)
    # Reactive value to store the final edited DMP table (exported as module output)
    dmp_final_table <- reactiveVal(NULL)
    # Reactive value to hold the dmpFinder results for debugging
    dmp_results_reactive <- reactiveVal(NULL)
    last_run_qvalue <- reactiveVal(NULL)
    # Inside moduleServer
    status_messages <- reactiveVal("")
    #--------------------------
    reset_analysis <- function() {
      combined_m_beta_tbl(NULL)
      dmp_final_table(NULL)
      dmp_results_reactive(NULL)
      output$dmp_table <- DT::renderDataTable(NULL)
      output$download_ui <- renderUI(NULL)
      status_messages("Ready for new analysis...")  # Reset to initial state
    }
    #----------------------------------
    # Display the reference group
    output$display_ref_group <- renderText({
      req(ref_group())
      ref_group()
    })
  
    #-------------------------------------
    # Toggle notes section
    observeEvent(input$toggle_notes, {
      shinyjs::toggle("note_section")
    })
    #-------------------------------------------
    # Simplified append_message function
    append_message <- function(msg) {
      current <- isolate(status_messages())
      status_messages(paste(current, msg, sep = "<br>"))
    }
    
    # Render the status messages
    output$status_message <- renderUI({
      HTML(status_messages())
    })
    #----------------------------------------------------
    # Function to save DMP results automatically
    save_dmp_results <- function(results_table, output_dir, last_qvalue) {
      dmp_dir <- file.path(output_dir, "DMP_results")
      if (!dir.exists(dmp_dir)) {
        dir.create(dmp_dir)
      }
      
      timestamp <- format(Sys.time(), "%Y%m%d")
      filename <- paste0("DMP_results_", last_qvalue, "qval_",timestamp, ".csv")
      filepath <- file.path(dmp_dir, filename)
      
      write.csv(results_table, filepath, row.names = TRUE)
      return(filepath)
    }
    #--------------------------------------
    # Function to save results as RDS
    save_results_rds <- function(dmp_table, last_qvalue, output_dir) {
      intermediate_dir <- file.path(output_dir, "intermediate_data")
      if (!dir.exists(intermediate_dir)) {
        dir.create(intermediate_dir, recursive = TRUE)
      }
      filename <- paste0("DMP_results_", last_qvalue, "qval.rds")
      filepath <- file.path(intermediate_dir, filename)
      
      saveRDS(list(
        dmp_table = dmp_table,
        last_qvalue = last_qvalue
      ), file = filepath)
      
      return(filepath)
    }
    #-----------------------------------------------------------
    # Observe changes to the q-value cutoff and reset the UI
    observeEvent(input$qvalue_cutoff, {
      message("q-value cutoff changed. Resetting UI.")
      reset_analysis()  # Use the reset function
    }, ignoreInit = TRUE)
    #----------------------------------------------------------
    # Observe the "Run DMP Finder" button click
    observeEvent(input$run_dmp_finder, {
      
      reset_analysis()  # Clear everything first
      status_messages("Starting DMP analysis... ⏳") 
      
      # ----------------- Debugging Step 0: Input Validation -----------------
      append_message("Step 0: Validating inputs...\n")
      req(filtered_data(), pheno(), ref_group(), project_output_dir())
      
      if (!dir.exists(project_output_dir())) {
        append_message(paste("Error: Project Directory", project_output_dir(), "does not exist. ❌ Analysis aborted.\n"))
        return()
      }
      
      withProgress(message = "DMP Analysis in Progress", value = 0, {
        
        # --- Step 1: Access MethylSet ---
        incProgress(0.1, detail = "Accessing normalized data...")
        append_message("Step 1: Accessing filtered data...")
        
        # NEW: More detailed methylset information
        current_GenoRationSet <- filtered_data()
        # debugging
        message(paste0("Methylset class: ", class(current_GenoRationSet)[1]))
        message(paste0("Dimensions: ", paste(dim(current_GenoRationSet), collapse = " x ")))
        message(paste0("First 3 samples: ", paste(sampleNames(current_GenoRationSet)[1:3], collapse = ", ")))
        message(paste0("Number of probes: ", nrow(current_GenoRationSet)))
        
        if (is.null(current_GenoRationSet)) {
          append_message("Error: Filtered data not found. ❌ Analysis aborted.\n")
          return()
        }
        
        # --- Step 2: Calculate Methylation Mean Differences ---
        incProgress(0.2, detail = "Calculating mean differences for Beta and M values...")
        append_message("Step 1: Calculating mean differences of Beta and M values for each CpG...\n")
        
        tryCatch({
          calculated_mean_diffs <- calculateMethylationMeanDifferences(
            GenoRationSet = current_GenoRationSet,
            pheno_table = pheno(),
            ref_group = ref_group(),
            round_digits = 5
          )
          combined_m_beta_tbl(calculated_mean_diffs)
          print(paste("number of rows in m_beta_combine table: ", nrow(calculated_mean_diffs)))
          append_message("✅ completed.\n")
        }, error = function(e) {
          append_message(paste0("Error in Step 1: ", e$message, ". ❌ Analysis aborted.\n"))
          return()
        })
        
        if (is.null(combined_m_beta_tbl())) {
          return()
        }
        
        # --- Step 2: Run DMP Finder ---
        incProgress(0.4, detail = "Running dmpFinder...")
       append_message("Step 2: Running dmpFinder...\n")
        
        print(paste("Data dimensions of dmpFinder input:", paste(dim(getM(current_GenoRationSet)), collapse = " x ")))
        print(paste("Pheno vector length:", length(pheno()$Sample_Group)))
        
        tryCatch({
          # count samples from pheno
          n_samples <- length(pheno()$Sample_Group)
          
          # set shrinkVar dynamically
          shrink_flag <- if (n_samples < 10) TRUE else FALSE
          
          # run dmpFinder with chosen shrinkVar
          dmp_raw <- suppressWarnings(
            minfi::dmpFinder(
              dat = getM(current_GenoRationSet),
              pheno = pheno()$Sample_Group,
              type = "categorical",
              qCutoff = input$qvalue_cutoff, # change to 1 since the table is filterable
              shrinkVar = shrink_flag
            )
          )
          dmp_results_reactive(dmp_raw)
          print(paste("number of rows in dmpfinder original results table: ", nrow(dmp_raw)))
          append_message("✅ completed.\n")
        }, error = function(e) {
         append_message(paste0("Error in Step 2 (dmpFinder): ", e$message, ". ❌ Analysis aborted.\n"))
          return()
        })
        
        if (is.null(dmp_results_reactive())) {
         append_message("Error: dmpFinder did not return any results. ❌Check your input data and parameters. \n")
          return()
        }
        
        if (nrow(dmp_results_reactive()) == 0) {
          append_message({
            paste0("DMP analysis completed. No significant DMPs were found at the specified q-value cutoff (FDR) of ", 
                   input$qvalue_cutoff, ". Please try a less stringent cutoff. ⚠️\n")
          })
          return()
        }
        
        # --- Step 3: Add Mean Differences to DMP table ---
        incProgress(0.2, detail = "Merging mean differences with DMP table...")
        append_message("Step 3: Merging results and saving the table in project directory...\n")
        
        tryCatch({
          final_dmp <- addMeanDifferencesToDMP(
            dmp_table = dmp_results_reactive(),
            combined_Beta_M_mean_table = combined_m_beta_tbl()
          )
          dmp_final_table(final_dmp)
          
          # Check for no results after final merging
          if (is.null(dmp_final_table()) || nrow(dmp_final_table()) == 0) {
           append_message("Error in adding beta and M values to DMP table. ❌  Analysis aborted.\n") 
           return()
          }
          
          saved_path <- save_dmp_results(final_dmp, project_output_dir(), input$qvalue_cutoff)
          append_message(paste0("✅ DMP analysis completed successfully! " , nrow(dmp_final_table()), " significant DMPs found. Results saved to: " , saved_path, "\n"))
        }, error = function(e) {
         append_message(paste0("Error in saving the DMP table: ", e$message, ". ❌ Analysis aborted.\n"))
          return()
        })
    
        incProgress(0.2, detail = "Rendering table and preparing downloads... ")
        output$dmp_table <- DT::renderDT({
          req(dmp_final_table())
          DT::datatable(dmp_final_table(),
                        options = list(
                          pageLength = 10,
                          searching = TRUE,
                          scrollX = TRUE
                        ),
                        rownames = TRUE,
                        filter = "top"
          )
        })
      })
      # At the very end of successful processing:
      if (!is.null(dmp_final_table()) && nrow(dmp_final_table()) > 0) {
        last_run_qvalue(input$qvalue_cutoff)  # Store the q-value used
        # Save the results as RDS
        rds_path <- save_results_rds(
          dmp_table = dmp_final_table(),
          last_qvalue = input$qvalue_cutoff,
          output_dir = project_output_dir()
        )
        
        append_message(paste0("📁 Saved intermediate results to: ", rds_path, "\n"))
      }
    })
    #-----------------------------------------------------
    # NEW: Render the DMP table based on the dmp_final_table reactive value
    observe({
      dmp_final_table()
  
      output$download_ui <- renderUI({
        if (is.null(dmp_final_table())) return()
        tagList(
          h5(strong("Download DMP Table:")),
          radioButtons(ns("download_format"), "Select format:",
                       choices = c("CSV" = "csv", "Excel" = "xlsx"),
                       selected = "csv", inline = TRUE),
          downloadButton(ns("download_dmp"), "Download Table",
                         icon = icon("download"))
        )
      })  
    })
    #-------------------------------------
    # Download handler
    output$download_dmp <- downloadHandler(
      filename = function() {
        paste0("dmp_results_", Sys.Date(), ".", input$download_format)
      },
      content = function(file) {
        if (input$download_format == "csv") {
          write.csv(dmp_final_table(), file, row.names = TRUE)
        } else if (input$download_format == "xlsx") {
          writexl::write_xlsx(as.data.frame(dmp_final_table()), file, row.names = TRUE)
        }
      }
    )
    #-------------------------------------
    # Return the final DMP table as a module output
    # Modify the return statement to always return the q-value, even when no DMPs are found
    return(list(
      DMP_tbl = reactive({ 
        if (is.null(dmp_final_table()) || nrow(dmp_final_table()) == 0) {
          NULL
        } else {
          dmp_final_table() 
        }
      }),
      last_qvalue = reactive({ 
        if (is.null(dmp_final_table()) || nrow(dmp_final_table()) == 0) {
          input$qvalue_cutoff  # Return the cutoff that was used, even if no DMPs found
        } else {
          last_run_qvalue() 
        }
      })
    ))
  })
}


'# test module
# app.R

# Load required libraries
library(shiny)
library(DT)
library(shinyjs)
library(writexl)
library(fs)
library(minfi)
library(bslib)
library(purrr)

# --- Define fixed variables and source files ---
# These are loaded once when the app starts.
dir <- "C:/Users/ghaza/Documents/ghazal/Bioinformatik_Fächer/Masterarbeit_Project/Scripts/R_Scripts/modules/main_app_tests/box_pheno_pass/intermediate_data/filtered_GRset_SWAN_SNPsremoved_SexChrProbes_kept_20250804.rds"
message("Attempting to load data from: ", dir)
tryCatch({
  GenoRationSet <- readRDS(dir)
  message("Data loaded successfully.")
}, error = function(e) {
  stop("Error loading RDS file: ", e$message)
})


# Extract pheno data and set reference group
pheno_table <- pData(GenoRationSet)
ref_group <- "unguided"
project_path <- "../main_app_tests/test"

# Source the DMP module and utilities
#source("11_dmp_module.R")
source("../utils/DMP_utils.R")

# --- Shiny UI ---
ui <- page_navbar(
  title = "EPIC Array Pipeline",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  
  nav_panel("DMP identification",
            dmp_UI("dmp")
  )
)

# --- Shiny Server ---
server <- function(input, output, session) {
  message("Server function has started.")
  
  # Solution 1: List of functions
  filtered_data_r <- reactive({GenoRationSet})
  
  # Other reactives remain the same
  pheno_r <- reactive(pheno_table)
  project_output_dir_r <- reactive(project_path)
  ref_group_r <- reactive(ref_group)

  # Call module
  dmp_Server("dmp",
             filtered_data = filtered_data_r,
             pheno = pheno_r,
             project_output_dir = project_output_dir_r,
             ref_group = ref_group_r
             )
}

# --- Run the application ---
shinyApp(ui, server)'