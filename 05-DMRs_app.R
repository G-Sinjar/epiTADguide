#test App
#filtered_rgset <- readRDS("./intermediate_data/filtered_GRset_SWAN_SNPsremoved_SexChrProbes_kept_20250608.rds")



# Load required libraries
library(shiny)
library(minfi)   # For methylation analysis
library(bslib)   # For Bootstrap 5 theming
library(DT)      # For interactive data tables
library(openxlsx)

# Load custom utility functions for DMR processing
source("./utils/dmrs_utils.R")

# Define the UI for DMR detection
ui_dmrs <- page_sidebar(
  sidebar = sidebar(
    width = "300px",
    
    # Input: methylation cutoff range
    numericInput("cutoff_from", "Cutoff from (min -1 to 0):", value = -0.15, min = -1, max = 0),
    numericInput("cutoff_to", "Cutoff to (0 to 1):", value = 0.15, min = 0, max = 1),
    helpText("Cutoff defines how large a difference in methylation must be to consider a region a potential DMR."),
    
    br(),
    
    # Input: number of permutations
    numericInput("B_val", "Number of permutations (B):", value = 0, min = 0),
    helpText("B controls the number of permutations used to assess significance, reducing false positives. More permutations = higher accuracy."),
    
    # Run button
    actionButton("run_dmr", "Detect DMRs"),
    helpText("Note: This step can take at least 5 minutes (B=0). Higher B values will increase runtime. For example B=100 takes 1.5 hours."),
    
    hr(),
    
    # Download options and button (only visible when results available)
    uiOutput("download_ui")
  ),
  
  # Main output panel
  mainPanel( fillable = TRUE,
    # Status display box
    card(fill = TRUE,
      card_title("DMR Identification Status"),
      style = "margin-bottom: 20px;",
      verbatimTextOutput("dmr_status", placeholder = TRUE)
    ),
    
    # Annotated DMR table
    div( # Optional div wrapper, gives you a handle for CSS if needed
      h4("Detected DMRs"),
      dataTableOutput(ns("dmr_table"), width = "100%")
      # You might want to add some styling to this div if you removed the card
      # e.g., style = "padding: 15px; margin-top: 20px;" to mimic some card spacing
    )
  )
)


# Define server logic for DMR detection and annotation
server_dmrs <- function(input, output, session) {
  
  # Reactive value to store status text
  dmr_status_text <- reactiveVal("Ready to identify DMRs.")
  
  # Reactive value to store final DMR table (after annotation)
  dmr_result <- reactiveVal(NULL)
  
  # Display current status
  output$dmr_status <- renderText({
    dmr_status_text()
  })
  
  # Main observer: Run on "Detect DMRs" button click
  observeEvent(input$run_dmr, {
    # Step 1 to 3: inside progress
    withProgress(message = "Running DMR detection...", value = 0, {
      
      if (is.null(filtered_rgset)) {
        dmr_status_text("❌ Error: Input data 'filtered_rgset' is missing.")
        return()
      }
      
      dmr_status_text("Step 1: Identifying DMRs...")
      dmr_result(NULL)
      
      cutoff_vals <- c(input$cutoff_from, input$cutoff_to)
      B_val <- input$B_val
      
      incProgress(0.1)
      dmrs_step1_result <- tryCatch({
        run_bumphunter_dmrs(rgSet = filtered_rgset, cutoff = cutoff_vals, B = B_val)
      }, error = function(e) {
        dmr_status_text(paste(dmr_status_text(), "\n❌ Error during DMR detection:\n", e$message))
        return(NULL)
      })
      
      if (is.null(dmrs_step1_result) || nrow(dmrs_step1_result) == 0) {
        dmr_status_text(paste(dmr_status_text(), "\n❌ Step 1 failed or no DMRs found."))
        return()
      }
      
      dmr_status_text(paste0(dmr_status_text(), sprintf("\n✅ Step 1 completed: %d DMRs found.", nrow(dmrs_step1_result))))
      
      # Step 2
      dmr_status_text(paste0(dmr_status_text(), "\nStep 2: Converting DMRs to genomic ranges (GRanges)..."))
      incProgress(0.3)
      GR_Dmrs <- tryCatch({
        prepare_dmrs_granges(dmrs_step1_result)
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
      incProgress(0.6)
      GR_Dmrs_annotated <- tryCatch({
        annotate_dmrs_with_genes(GR_Dmrs, tx_gr_filtered)
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
    })  # <-- End of withProgress here
    
    # Step 4: Now outside withProgress, safe for reactive updates
    dmr_status_text(paste0(dmr_status_text(), "\nStep 4: Preparing final DMR table..."))
    df <- as.data.frame(GR_Dmrs_annotated)
    if ("seqnames" %in% colnames(df)) {
      colnames(df)[colnames(df) == "seqnames"] <- "chr"
    }
    if ("DMR_ID" %in% colnames(df)) {
      df <- df[, c("DMR_ID", setdiff(colnames(df), "DMR_ID"))]
    }
    dmr_result(df)
    dmr_status_text(paste0(dmr_status_text(), "\n✅ Step 4 completed: DMR table is ready for download."))
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
      radioButtons("download_format", "Choose format:", choices = c("CSV", "Excel"), inline = TRUE),
      downloadButton("download_dmr", "Download DMR Table")
    )
  })
  output$download_dmr <- downloadHandler(
    filename = function() {
      ext <- if (input$download_format == "Excel") "xlsx" else "csv"
      paste0("DMRs_", Sys.Date(), ".", ext)
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
}

# Run the Shiny app
shinyApp(ui = ui_dmrs, server = server_dmrs)


