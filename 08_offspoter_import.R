# app.R
# Author: Ghazal Sinjar
# Date: 29.06.2025
# Description: Shiny app for importing and processing off-target guide files for EPIC array data.
# Users input guide file names and directory, view processed results, and download outputs.

library(shiny)
library(bslib)
library(readr)
library(dplyr)
library(stringr)
library(tools)
library(GenomicRanges)
library(IRanges)
library(S4Vectors)
library(DT)
library(writexl)

# Source external utility functions for processing
#source("./utils/offTargets_processing_utils.R")

# ------------------- UI -------------------

ui <- page_navbar(
  title = "EPIC Array Pipeline",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  
  nav_panel(
    "Offtargets import",
    layout_sidebar(
      fillable = TRUE,
      sidebar = sidebar(
        h4("Input Data"),
        
        # Input: Directory path for guide files
        textInput(
          "directory_path",
          "Directory Path:",
          placeholder = "e.g. C:\\path\\to\\data"
        ),
        
        # Input: List of guide filenames (1 per line)
        textAreaInput(
          "guide_file_names",
          "Guide File Names (one per line):",
          placeholder = "e.g.,\nguide1.txt\nguide2a.txt\nguide2b.txt",
          rows = 5
        ),
        
        # Helpful tip for users
        p(tags$small(
          "Note: Ensure your guide files are tab-separated, located in the specified directory, and named using only the guide ID (e.g., 'guide1.txt', 'guide2a.txt')."
        )),
        
        # Action: Trigger file reading
        actionButton("run_processing", "Process Guide Files into one file", class = "btn-primary"),
        hr(),
        
        # Conditional: Show download options once data is ready
        conditionalPanel(
          condition = "output.tableReady === true",
          h4("Download Results"),
          selectInput("download_format", "Select Format:", choices = c("CSV", "Excel")),
          downloadButton("download_data", "Download", class = "btn-success")
        )
      ),
      
      # Main panel: Status + results table
      div(
        style = "padding-left: 15px; padding-right: 15px;",
        layout_columns(
          col_widths = c(12),
          
          # Card: Display messages/warnings/errors
          card(
            card_header("Off-targets File Processing Status"),
            card_body(uiOutput("status_output_card_body"))
          ),
          
          # Table: Display combined guide data
          div(
            h5("Combined all guides offtargets"),
            uiOutput("row_count_text"),
            DTOutput("offtargets_table", width = "100%")
          )
        )
      )
    )
  )
)

# ------------------- SERVER -------------------

server <- function(input, output, session) {
  
  output_messages <- reactiveVal(list())  # Stores messages for UI
  processed_data <- reactiveVal(NULL)     # Stores final data
  
  # DEBUG: main processing observer
  observeEvent(input$run_processing, {
    
    output_messages(list())       # reset messages
    processed_data(NULL)          # reset old table
    
    msgs <- list()
    all_guides <- list()
    
    dir_path <- gsub("\\\\", "/", input$directory_path)
    msgs <- append(msgs, paste0("🛠️ Processing started... Directory: ", dir_path))  # DEBUG
    
    if (!dir.exists(dir_path)) {
      msgs <- append(msgs, "❌ Error: The path you provided doesn't exist. Please check and try again.")
      output_messages(msgs)
      return()
    } 
    
    # DEBUG: parse guide file names
    file_names <- strsplit(input$guide_file_names, "\\s+")[[1]]
    
    col_names <- c("Chrom", "Strand", "Start", "End", "Given query", "Actual genomic hit",
                   "Number_of_mismatches", "Pre-mRNA (Unspliced)", "mRNA (5UTR)",
                   "mRNA (CDS)", "mRAN (3UTR)", "lincRNA (Unspliced)",
                   "lincRNA (Spliced)", "GC content")
    
    for (file_name in file_names) {
      file_path <- file.path(dir_path, file_name)
      msgs <- append(msgs, paste0("🔍 Reading file: ", file_path))  # DEBUG
      
      guide_data <- tryCatch({
        read_delim(file_path, delim = "\t", col_names = FALSE, show_col_types = FALSE)
      }, error = function(e) {
        msgs <<- append(msgs, paste0("❌ Step 1 - Reading failed: ", e$message))
        return(NULL)
      })
      if (is.null(guide_data)) next
      msgs <- append(msgs, "✅ Step 1 - File read successfully")
      
      # Step 2: Check columns
      if (ncol(guide_data) != length(col_names)) {
        msgs <- append(msgs, paste0("❌ Step 2 - Column mismatch: expected ", length(col_names), ", got ", ncol(guide_data)))
        next
      }
      msgs <- append(msgs, "✅ Step 2 - Column count verified")
      
      # Step 3: Processing
      processed <- tryCatch({
        colnames(guide_data) <- col_names
        guide_data <- guide_data[, !colnames(guide_data) %in% "GC content"]
        guide_data$guide <- tools::file_path_sans_ext(file_name)
        guide_data$Start <- as.integer(guide_data$Start)
        guide_data$End <- as.integer(guide_data$End)
        guide_data$Start_minus_70 <- guide_data$Start - 70
        guide_data$End_plus_70 <- guide_data$End + 70
        guide_data$Number_of_mismatches <- as.integer(guide_data$Number_of_mismatches)
        guide_data
      }, error = function(e) {
        msgs <<- append(msgs, paste0("❌ Step 3 - Processing failed: ", e$message))
        return(NULL)
      })
      if (is.null(processed)) next
      if (!is.null(processed)) {
        msgs <- append(msgs, "✅ Step 3 - Data processed")
      }
      
      all_guides[[file_name]] <- processed
      msgs <- append(msgs, "✅ Step 4 - Data in the guide file stored")
    }
    
    
    # Step 5: Combine all
    msgs <- append(msgs, "📦 Step 5 - All off-targets in one file...")
    combined <- tryCatch({
      bind_rows(all_guides)
    }, error = function(e) {
      msgs <<- append(msgs, paste0("❌ Step 5a - Combining all off-targets of successfully processed guide files : ", e$message))
      return(NULL)
    })
    if (is.null(combined)) {
      output_messages(msgs)
      return()
    }
    msgs <- append(msgs, "✅ Step 5a- all off-targets of successfully processed guide files are combined")
    
    # Step 5b: Filter
    final <- tryCatch({
      filtered <- combined %>%
        dplyr::filter(Number_of_mismatches != 0) %>%
        dplyr::mutate(
          guide_num = stringr::str_match(guide, "(?i)guide([[:alnum:]]+)")[, 2],
          ID = paste0("G", guide_num, ".M", Number_of_mismatches, ".chr", Chrom, "_", Start)
        ) %>%
        dplyr::select(-guide_num) %>%
        dplyr::select(ID, guide, Chrom, Start, End, Number_of_mismatches,
                      Start_minus_70, End_plus_70, dplyr::everything())
      
      filtered
    }, error = function(e) {
      msgs <<- append(msgs, paste0("❌ Step 5b - Filtering out on-targets failed: ", e$message))
      return(NULL)
    })
    
    if (is.null(final)) {
      output_messages(msgs)
      return()
    }
    if (!is.null(final)) {
      msgs <- append(msgs, "✅ Step 5b - filtering out on-targets")
    }
    
    processed_data(final)
    output_messages(msgs)
  })
  
  output$status_output_card_body <- renderUI({
    tags$ul(
      style = "padding-left: 1.2em; margin: 0;",
      lapply(output_messages(), function(msg) {
        tags$li(style = "margin-bottom: 4px;", msg)
      })
    )
  })
  
  
  output$offtargets_table <- renderDT({
    req(processed_data())  
    processed_data()
  }, options = list(pageLength = 25))
  
  
  output$row_count_text <- renderUI({
    req(processed_data())
    n <- nrow(processed_data())
    HTML(paste0("<p><strong>", n, "</strong> off-targets found</p>"))
  })
}


# Run the Shiny application
shinyApp(ui = ui, server = server)
