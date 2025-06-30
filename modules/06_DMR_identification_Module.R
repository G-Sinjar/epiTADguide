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


# ─────────────────────────────────────
# User Interface (UI) FUNCTION
# ─────────────────────────────────────
dmrs_ui <- function(id) {
  ns <- NS(id) 
  
  page_sidebar(
    sidebar = sidebar(
      width = "300px",
      
      # Input: methylation cutoff range
      numericInput(ns("cutoff_from"), "Cutoff from (min -1 to 0):", value = -0.15, min = -1, max = 0),
      numericInput(ns("cutoff_to"), "Cutoff to (0 to 1):", value = 0.15, min = 0, max = 1),
      helpText("Cutoff defines how large a difference in methylation must be to consider a region a potential DMR."),
      
      br(),
      
      # Input: number of permutations
      numericInput(ns("B_val"), "Number of permutations (B):", value = 0, min = 0),
      helpText("B controls the number of permutations used to assess significance, reducing false positives. More permutations = higher accuracy."),
      
      # Run button
      actionButton(ns("run_dmr"), "Detect DMRs"),
      helpText("Note: This step can take at least 5 minutes (B=0). Higher B values will increase runtime. For example B=100 takes 1.5 hours."),
      
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
      # This div will stretch to fill the available width
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

    # Display current status
    output$dmr_status <- renderText({dmr_status_text()})
    
    # Main observer: Run on "Detect DMRs" button click
    observeEvent(input$run_dmr, {
      # Ensure filtered_rgset_reactive has a value before proceeding
      req(filtered_rgset_reactive())
      # Also ensure tx_gr_filtered_static is not NULL, as it's directly used
      req(tx_gr_filtered_static)
      
      
      # Step 1 to 4: all inside progress
      withProgress(message = "Running DMR detection...", value = 0, {
        
        # Access the reactive value
        current_rgset <- filtered_rgset_reactive()
        
        if (is.null(current_rgset)) {
          dmr_status_text("❌ Error: Input data 'filtered_rgset' is missing.")
          return()
        }
        
        dmr_status_text("Step 1: Identifying DMRs...")
        if (input$cutoff_from >= input$cutoff_to) {
          dmr_status_text("❌ Cutoff 'from' must be less than 'to'.")
          return()
        }
        
        cutoff_vals <- c(input$cutoff_from, input$cutoff_to)
        B_val <- input$B_val
        
        incProgress(0.1, detail = "Running bumphunter...")
        dmrs_step1_result <- tryCatch({
          result <- run_bumphunter_dmrs(rgSet = current_rgset, cutoff = cutoff_vals, B = B_val)
          
          # ✅ Sanity check to avoid NULL or malformed return structure
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
        
        # Step 2 converting DMRs to genomic ranges
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
        
        # Step 3
        dmr_status_text(paste0(dmr_status_text(), "\nStep 3: Annotating DMRs with gene information..."))
        incProgress(0.6, detail = "Annotating with genes...")
        GR_Dmrs_annotated <- tryCatch({
          # tx_gr_filtered_static is correctly passed as a static (non-reactive) object
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
        
        # Step 4 Preparing final DMR table
        dmr_status_text(paste0(dmr_status_text(), "\nStep 4: Preparing final DMR table..."))
        incProgress(0.9, detail = "Finalizing table...") # Added progress for step 4
        
        df <- as.data.frame(GR_Dmrs_annotated)
        if ("seqnames" %in% colnames(df)) {
          colnames(df)[colnames(df) == "seqnames"] <- "chr"
        }
        if ("DMR_ID" %in% colnames(df)) {
          df <- df[, c("DMR_ID", setdiff(colnames(df), "DMR_ID"))]
        }
        dmr_result(df) # Update the reactiveVal
        dmr_status_text(paste0(dmr_status_text(), "\n✅ Step 4 completed: DMR table is ready for download."))
        # Save both DMR and phenotype results to RDS
        base_filename <- get_dmr_base_filename(input$cutoff_from, input$cutoff_to, input$B_val)
        
        'if (!dir.exists("intermediate_data")) dir.create("intermediate_data")
        output_path <- file.path("intermediate_data", paste0(base_filename, ".rds"))
        
        saveRDS(list(
          dmr_table = dmr_result(),
          pheno_data = pheno_result()
        ), file = output_path)
        dmr_status_text(paste0(dmr_status_text(), "\n📁 Saved full results automatically to ", output_path))'
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