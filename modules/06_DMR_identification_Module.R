# 06_DMR_identification_Module.R
# Author: Ghazal Sinjar
# Date: 17.06.2025
# Description:
# This Shiny module provides UI and server logic for identifying Differentially Methylated Regions (DMRs).
# It utilizes the bumphunter algorithm to detect DMRs based on specified methylation cutoff values and permutation settings.
# The module outputs an annotated table of detected DMRs, which can be downloaded in CSV or Excel format.

#-------------------------------------------------------------------
# libraries especially for this module
library(DT)        # For interactive data tables
library(openxlsx)
library(GenomicRanges)
library(shiny) # Add shiny explicitly here for module to be self-contained in its library requirements
library(minfi) # Add minfi explicitly here
library(bslib) # Add bslib explicitly here
library(EnsDb.Hsapiens.v86)
# In your app.R or global.R
library(shinyjs)


# ─────────────────────────────────────
# User Interface (UI) FUNCTION
# ─────────────────────────────────────
dmrs_ui <- function(id) {
  ns <- NS(id)
  
  page_sidebar(
    sidebar = sidebar(
      width = "300px",
      
      useShinyjs(), # Important for shinyjs functions to work
      
      # Input: methylation cutoff range
      numericInput(ns("cutoff_from"), "Cutoff from (min -1 to 0):", value = -0.15, min = -1, max = 0, step= 0.01),
      numericInput(ns("cutoff_to"), "Cutoff to (0 to 1):", value = 0.15, min = 0, max = 1, step = 0.01),
      
      # --- CUTOFF HELP IMPLEMENTATION START ---
      actionLink(ns("toggle_cutoff_help"), "Click for help on cutoff."),
      div(
        id = ns("cutoff_help_text_div"),
        style = "display: none;",
        helpText("Cutoff defines how large a difference in methylation must be to consider a region a potential DMR.")
      ),
      # --- CUTOFF HELP IMPLEMENTATION END ---
      
      br(),
      
      # Input: number of permutations
      numericInput(ns("B_val"), "Number of permutations (B):", value = 0, min = 0),
      
      # --- PERMUTATIONS HELP IMPLEMENTATION START ---
      actionLink(ns("toggle_B_help"), "Click for help on permutations."),
      div(
        id = ns("B_help_text_div"),
        style = "display: none;",
        helpText("B controls the number of permutations used to assess significance, reducing false positives. More permutations = higher accuracy.")
      ),
      # --- PERMUTATIONS HELP IMPLEMENTATION END ---
      
      br(), # Add a break for spacing
      
      # Core options radio buttons
      radioButtons(
        ns("core_choice"),
        "Choose number of CPU cores:",
        choices = c(
          "Detected physical cores - 1" = "auto_cores",
          "Choose cores manually" = "manual_cores"
        ),
        selected = "auto_cores"
      ),
      uiOutput(ns("detected_cores_info")), 
      # New: Numeric input for manual core selection (dynamically rendered)
      uiOutput(ns("manual_cores_input")),
      
      # --- CORE HELP IMPLEMENTATION START ---
      actionLink(ns("toggle_cores_help"), "Click for help on core choice."),
      div(
        id = ns("cores_help_text_div"),
        style = "display: none;",
        helpText("This setting controls the number of CPU cores used for calculations. Using more cores can significantly speed up the analysis, especially for larger datasets or when running permutations (B > 0).\n\n",
                 "  * **'Detected cores - 1'**: Recommended for most users. This option automatically uses all but one of your computer's available CPU cores. This leaves one core free for system operations, preventing your computer from becoming unresponsive.\n",
                 "  * **'Choose cores manually'**: Allows you to specify an exact number of cores. Use this if you have specific performance requirements or if 'Detected cores - 1' does not suit your needs (e.g., if you want to use fewer cores to save resources for other applications)."
        )
      ),
      # --- CORE HELP IMPLEMENTATION END ---
      
      
      
      br(), # Add a break for spacing (optional, but can help visual spacing)
      
      # Run button
      actionButton(ns("run_dmr"), "Detect DMRs"),
      helpText("Note: This step can take at least 5 minutes (B=0). Higher B values will increase runtime. For example B=100 takes 1.5 hours using 6 Cores."),
      
      hr(),
      
      # Download options and button (only visible when results available)
      uiOutput(ns("download_ui"))
    ),
    
    # Main output panel
    div(
      style = "padding-left: 15px; padding-right: 15px;", # Add right padding for balance
      layout_columns(
        col_widths = c(12), # This specifies one column that takes up all 12 Bootstrap columns
        fill = TRUE, # Make sure the card fills the available space
        card(
          card_title("DMR Identification Status"),
          verbatimTextOutput(ns("dmr_status"), placeholder = TRUE)
        )
      ),
      
      # Bottom section: DMR Table
      div(
        style = "margin-top: 20px;", # Add some space above the table section
        h4("Detected DMRs"),
        dataTableOutput(ns("dmr_table")) 
      )
    )
  )
}


# ─────────────────────────────────────
# SERVER FUNCTION
# ─────────────────────────────────────
dmrs_server <- function(id, filtered_rgset_reactive, tx_gr_filtered_static) {
  moduleServer(id, function(input, output, session) {
    
    # Reactive value to store status text
    dmr_status_text <- reactiveVal("Ready to identify DMRs.")
    # Reactive value to store final DMR and pheno tables
    dmr_result <- reactiveVal(NULL)
    pheno_result <- reactiveVal(NULL)
    
    # Reactive value to store the detected number of cores for display
    detected_cores_display <- reactiveVal(NULL)
    
    # Display current status
    output$dmr_status <- renderText({dmr_status_text()})
    
    # Reactive expression for the number of cores to be used in calculations
    num_cores_to_use <- reactive({
      if (input$core_choice == "auto_cores") {
        if (!requireNamespace("parallel", quietly = TRUE)) {
          warning("Package 'parallel' is required for automatic core detection. Defaulting to 1 core.")
          return(1)
        }
        max(1, parallel::detectCores(logical = FALSE) - 1)
      } else {
        req(input$manual_cores)
        input$manual_cores
      }
    })
    
    # Dynamic UI for manual core selection
    output$manual_cores_input <- renderUI({
      ns <- session$ns
      if (input$core_choice == "manual_cores") {
        max_detected_cores <- if (requireNamespace("parallel", quietly = TRUE)) {
          parallel::detectCores()
        } else {
          6 # Fallback if parallel package is not available
        }
        tagList(
          numericInput(
            ns("manual_cores"),
            "Number of cores:",
            value = max(1, max_detected_cores - 1),
            min = 1,
            max = max_detected_cores
          ),
          helpText(
            paste0(
              "Note: The default number that appears at first is the number of logical cores detected minus one. ",
              "The detected number of your logical cores is ", max_detected_cores, ", which is the maximum, and it is not recommended to be used."
            )
          )
        )
      }
    })
    
    # Observe core_choice to display detected cores information
    observeEvent(input$core_choice, {
      if (input$core_choice == "auto_cores") {
        current_detected_cores <- if (requireNamespace("parallel", quietly = TRUE)) {
          max(1, parallel::detectCores(logical = FALSE) - 1)
        } else {
          1 # Fallback if 'parallel' package is not available
        }
        detected_cores_display(paste("The number of physical CPUs/cores detected: ",current_detected_cores+1 , " \n.The number of cores to be used for the analysis: ",current_detected_cores, "."))
      } else {
        detected_cores_display(NULL) # Clear the text if manual is chosen
      }
    })
    
    # Render the detected cores information UI
    output$detected_cores_info <- renderUI({
      if (!is.null(detected_cores_display())) {
        helpText(detected_cores_display())
      }
    })
    
    # --- HINT SERVER LOGIC ADDITIONS START ---
    observeEvent(input$toggle_cutoff_help, {
      shinyjs::toggle("cutoff_help_text_div")
    })
    
    observeEvent(input$toggle_B_help, {
      shinyjs::toggle("B_help_text_div")
    })
    
    observeEvent(input$toggle_cores_help, {
      shinyjs::toggle("cores_help_text_div")
    })
    # --- HINT SERVER LOGIC ADDITIONS END ---
    
    # Main observer: Run on "Detect DMRs" button click
    observeEvent(input$run_dmr, {
      req(filtered_rgset_reactive())
      req(tx_gr_filtered_static)
      
      current_num_cores <- num_cores_to_use()
      if (is.null(current_num_cores) || current_num_cores < 1) {
        dmr_status_text("❌ Error: Invalid number of cores selected.")
        return()
      }
      
      message(paste("Using", current_num_cores, "cores for bumphunter analysis."))
      dmr_status_text(paste0("Step 1: Identifying DMRs using ", current_num_cores, " cores..."))
      
      withProgress(message = "Running DMR detection...", value = 0, {
        
        current_rgset <- filtered_rgset_reactive()
        
        if (is.null(current_rgset)) {
          dmr_status_text("❌ Error: Input data 'filtered_rgset' is missing.")
          return()
        }
        
        if (input$cutoff_from >= input$cutoff_to) {
          dmr_status_text("❌ Cutoff 'from' must be less than 'to'.")
          return()
        }
        
        cutoff_vals <- c(input$cutoff_from, input$cutoff_to)
        B_val <- input$B_val
        
        incProgress(0.1, detail = paste("Running bumphunter with", current_num_cores, "cores..."))
        dmrs_step1_result <- tryCatch({
          result <- run_bumphunter_dmrs(
            rgSet = current_rgset,
            cutoff = cutoff_vals,
            B = B_val,
            num_cores = current_num_cores
          )
          
          if (!("DMR_table" %in% names(result))) {
            stop("Missing 'DMR_table' in the result returned by run_bumphunter_dmrs().")
          }
          if (!("pd" %in% names(result))) {
            stop("Missing 'pd' in the result returned by run_bumphunter_dmrs().")
          }
          
          pheno_result(result$pd)
          result
        }, error = function(e) {
          dmr_status_text(paste(dmr_status_text(), "\n❌ Error during DMR detection:\n", e$message))
          return(NULL)
        })
        
        if (is.null(dmrs_step1_result)) {
          dmr_status_text(paste(dmr_status_text(), "\n❌ Step 1 failed ."))
          return()
        }
        if (nrow(dmrs_step1_result$DMR_table) == 0) {
          dmr_status_text(paste(dmr_status_text(), "\n❌ No DMRs found."))
          return()
        }
        
        dmr_status_text(paste0(dmr_status_text(), sprintf("\n✅ Step 1 completed: %d DMRs found.", nrow(dmrs_step1_result$DMR_table))))
        
        dmr_status_text(paste0(dmr_status_text(), "\nStep 2: Converting DMRs to genomic ranges (GRanges)..."))
        incProgress(0.3, detail = "Converting to GRanges...")
        GR_Dmrs <- tryCatch({
          prepare_dmrs_granges(dmrs_step1_result$DMR_table)
        }, error = function(e) {
          dmr_status_text(paste(dmr_status_text(), "\n❌ Error in Step 2:", e$message))
          return(NULL)
        })
        
        if (is.null(GR_Dmrs)) {
          return()
        }
        dmr_status_text(paste0(dmr_status_text(), "\n✅ Step 2 completed!"))
        
        dmr_status_text(paste0(dmr_status_text(), "\nStep 3: Annotating DMRs with gene information..."))
        incProgress(0.6, detail = "Annotating with genes...")
        GR_Dmrs_annotated <- tryCatch({
          annotate_dmrs_with_genes(GR_Dmrs, tx_gr_filtered_static)
        }, error = function(e) {
          dmr_status_text(paste(dmr_status_text(), "\n❌ Error in Step 3 (annotation):", e$message))
          return(NULL)
        })
        
        if (is.null(GR_Dmrs_annotated)) {
          return()
        }
        
        num_annotated <- sum(!is.na(mcols(GR_Dmrs_annotated)$overlapped_gene_name))
        dmr_status_text(paste0(
          dmr_status_text(),
          "\n✅ Step 3 completed: ", num_annotated, " DMRs overlapped with genes."
        ))
        
        dmr_status_text(paste0(dmr_status_text(), "\nStep 4: Preparing final DMR table..."))
        incProgress(0.9, detail = "Finalizing table...")
        
        df <- as.data.frame(GR_Dmrs_annotated)
        if ("seqnames" %in% colnames(df)) {
          colnames(df)[colnames(df) == "seqnames"] <- "chr"
        }
        if ("DMR_ID" %in% colnames(df)) {
          df <- df[, c("DMR_ID", setdiff(colnames(df), "DMR_ID"))]
        }
        dmr_result(df)
        dmr_status_text(paste0(dmr_status_text(), "\n✅ Step 4 completed: DMR table is ready for download."))
        
        base_filename <- get_dmr_base_filename(input$cutoff_from, input$cutoff_to, input$B_val)
        
        if (!dir.exists("intermediate_data")) dir.create("intermediate_data")
        output_path <- file.path("intermediate_data", paste0(base_filename, ".rds"))
        
        saveRDS(list(
          dmr_table = dmr_result(),
          pheno_data = pheno_result()
        ), file = output_path)
        dmr_status_text(paste0(dmr_status_text(), "\n📁 Saved full results automatically to ", output_path))
      })
    })
    
    output$dmr_table <- renderDataTable({
      req(dmr_result())  
      datatable(dmr_result(), options = list(scrollX = TRUE,
                                             pageLength = 10,
                                             autoWidth = TRUE))
    })
    
    output$download_ui <- renderUI({
      req(dmr_result())
      tagList(
        radioButtons(session$ns("download_format"), "Choose format:", choices = c("CSV", "Excel"), inline = TRUE),
        downloadButton(session$ns("download_dmr"), "Download DMR Table")
      )
    })
    
    output$download_dmr <- downloadHandler(
      filename = function() {
        base_name <- get_dmr_base_filename(input$cutoff_from, input$cutoff_to, input$B_val)
        ext <- if (input$download_format == "Excel") "xlsx" else "csv"
        paste0(base_name, "_", Sys.Date(), ".", ext)
      },
      content = function(file) {
        df <- dmr_result()
        if (input$download_format == "Excel") {
          if (!requireNamespace("openxlsx", quietly = TRUE)) {
            stop("Package 'openxlsx' is required for Excel export. Please install it.")
          }
          openxlsx::write.xlsx(df, file)
        } else {
          write.csv(df, file, row.names = FALSE)
        }
      }
    )
    
    return(list(
      dmr_table = dmr_result,
      pheno = pheno_result
    )
    )
  })
}


'# test module
## libraries in main app
library(shiny)
library(minfi)     # For methylation analysis
library(bslib)     # For Bootstrap 5 theming

## load input data
filtered_rgset <- readRDS("../intermediate_data/filtered_GRset_SWAN_SNPsremoved_SexChrProbes_kept_20250608.rds")

# Load custom utility functions for DMR processing
# In a module, you might consider if these utils should always loaded by the main app. For now, well keep the source here.
# Adjust path for sourcing relative to the module file itself
source("../utils/dmrs_utils.R")

## creating the gene granges for the input
library(EnsDb.Hsapiens.v86)
edb <- EnsDb.Hsapiens.v86
options(ucscChromosomeNames = TRUE)
#To extract specific data, such as transcript or gene information, from an EnsDb object, you can use functions like transcripts(), genes(), and exons().
######################### Get gene ranges
tx_gr <- genes(edb)
#length(tx_gr) # 63970
head(tx_gr)
# Filter to standard chromosomes only????
tx_gr_filtered <- keepSeqlevels(tx_gr, standardChromosomes(tx_gr), pruning.mode = "coarse")
#length(tx_gr_filtered) #58650
#tx_gr_filtered
seqlevelsStyle(tx_gr_filtered) <- "UCSC"
#tx_gr_filtered




# UI
ui <- page_navbar(
  title = "EPIC Array Pipeline",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  
  nav_panel("DMR identification",
            dmrs_ui("dmrs")
  )
)
# Server
server <- function(input, output, session) {
  # Wrap your object as a reactive expression
  filtered_rgset_reactive <- reactive({ filtered_rgset })
  
  # Call the module
  annotated_table <- dmrs_server("dmrs", filtered_rgset_reactive,tx_gr_filtered)
}

shinyApp(ui, server)
'