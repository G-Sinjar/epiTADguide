# 01_loadData_app.R
# Author: Ghazal Sinjar
# Date: 30.05.2025

# Load required libraries
library(shiny)  # For building interactive web applications
library(minfi)  # For processing Illumina methylation arrays
library(IlluminaHumanMethylationEPICv2manifest)  # Manifest for EPIC v2 arrays
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)  # Annotation package for EPIC v2 arrays
library(bslib)  # For Bootstrap 5 theming in Shiny apps

# User interface
ui <- page_navbar(
  title = "EPIC Array Pipeline",
  theme = bs_theme(version = 5, bootswatch = "flatly"),  # Use "flatly" theme
  
  nav_panel("Load Raw Data",
            layout_sidebar(
              sidebar = sidebar(
                strong("Choose the path to raw data from sequencer"),
                
                # Text box for entering raw data folder path
                textInput("manual_path", "Enter folder path manually:",
                          value = "",
                          placeholder = "E.g. C:\\Users\\yourname\\Documents\\sequencer_output"),
                
                # Button to start loading process
                actionButton("load_btn", "Load Data"),
                
                # Help text explaining required folder structure
                helpText("Please select a folder that includes: (1) a CSV file with slide and sample information, and (2) a subfolder named after the slide, which contains the green and red IDAT files.")
              ),
              layout_columns(
                # Card showing the status of loading process
                card(
                  card_header("Load Status"),
                  verbatimTextOutput("status")
                ),
                # Card showing the manifest details
                card(
                  card_header("Array Manifest"),
                  verbatimTextOutput("manifest_output")
                )
              )
            )
  )
)

# Server logic
server <- function(input, output, session) {
  
  # Holds and updates status messages
  status_log <- reactiveVal("Waiting for user input...")
  
  # Helper function to append new messages to the status
  append_status <- function(new_message) {
    current <- status_log()
    updated <- paste(current, new_message, sep = "\n")
    status_log(updated)
  }
  
  # Render the status log in the UI
  output$status <- renderText({ status_log() })
  
  # Observe when user clicks "Load Data"
  observeEvent(input$load_btn, {
    # Step 1: Validate path
    status_log("Step 1: Checking path...")
    req(input$manual_path)
    
    # Normalize the path (replace backslashes with forward slashes)
    RawDataDir <- gsub("\\\\", "/", input$manual_path)
    
    # Check if the directory exists
    if (!dir.exists(RawDataDir)) {
      append_status(paste0("❌ Path does not exist: ", RawDataDir, "\nPlease check the path and try again."))
      return()
    }
    append_status("✅ Step 1: completed.")
    
    # Step 2: Look for CSV sample sheet
    append_status("Step 2: Reading Sample Sheet...")
    csv_files <- list.files(RawDataDir, pattern = "\\.csv$", full.names = TRUE)
    if (length(csv_files) == 0) {
      append_status("❌ No CSV Sample Sheet found in the selected directory.")
      return()
    }
    
    # Try reading the sample sheet using minfi
    targets <- tryCatch({
      read.metharray.sheet(RawDataDir)
    }, error = function(e) {
      append_status(paste("❌ Error reading Sample Sheet:", e$message))
      return(NULL)
    })
    if (is.null(targets)) return()
    append_status("✅ Step 2: completed.")
    
    # Step 3: Load raw IDAT files
    append_status("Step 3: Reading IDAT files...")
    RGset <- tryCatch({
      read.metharray.exp(targets = targets)
    }, error = function(e) {
      append_status(paste("❌ Error reading IDATs:", e$message))
      return(NULL)
    })
    if (is.null(RGset)) return()
    append_status("✅ Step 3: completed.")
    
    # Step 4: Retrieve the array manifest
    append_status("Step 4: Getting Manifest...")
    manifest <- tryCatch({
      getManifest(RGset)
    }, error = function(e) {
      append_status(paste("❌ Error getting manifest:", e$message))
      return(NULL)
    })
    if (is.null(manifest)) return()
    append_status("✅ Step 4: completed")
    
    # Output manifest details
    output$manifest_output <- renderPrint({ manifest })
    
    # Step 5: Perform raw preprocessing
    append_status("Step 5: Converting raw signals to locus level...")
    raw_normalised <- tryCatch({
      preprocessRaw(RGset)
    }, error = function(e) {
      append_status(paste("❌ Error preprocessing:", e$message))
      return(NULL)
    })
    if (is.null(raw_normalised)) return()
    append_status("✅ Step 5: completed.")
    
    # Step 6: Display sample and CpG statistics
    append_status("Step 6: Extracting sample info and CpG count...")
    
    sample_names <- sampleNames(RGset)
    num_samples <- length(sample_names)
    num_cpgs <- nrow(raw_normalised)
    
    append_status(paste0(
      "✅ Sample and CpG info retrieved:\n",
      "  - Number of CpGs: ", num_cpgs, "\n",
      "  - Number of Samples: ", num_samples, "\n",
      "  - Sample Names:\n", paste(sample_names, collapse = "\n")
    ))
    
    # Optional: Save processed data
    saveRDS(list(RGset = RGset, raw_normalised = raw_normalised, targets = targets), "./intermediate_data/preprocessed_data.rds")
    append_status("✅ Step 6: Data saved to 'preprocessed_data.rds'")
  })
}

# Run the Shiny application
shinyApp(ui, server)
