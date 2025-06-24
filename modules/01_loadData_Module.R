# loadDataModule.R
# Author: Ghazal Sinjar
# Date: 30.05.2025
# Description: This Shiny module allows users to load EPIC array data
# by specifying a path to a directory containing a sample sheet (.csv) 
# and corresponding IDAT files. The data is preprocessed and summarized.

# ─────────────────────────────────────
# User Interface (UI) FUNCTION
# ─────────────────────────────────────

loadDataUI <- function(id) {
  ns <- NS(id)   # 'id' allows unique namespacing of inputs and outputs in this module
  
  # Return only the content inside the tab 
  layout_sidebar(
    sidebar = sidebar(
      strong("Choose the path to raw data from sequencer"),
      
      # Text box input for directory path from user
      textInput(
        inputId = ns("manual_path"),
        label = "Enter folder path manually:",
        value = "",
        placeholder = "E.g. C:\\Users\\yourname\\Documents\\sequencer_output"
      ),
      
      # Button to trigger loading of data
      actionButton(ns("load_btn"), "Load Data"),
      
      # Help text with folder requirements
      helpText(
        "Please select a folder that includes: (1) a CSV file with slide and sample information, ",
        "and (2) a subfolder named after the slide, which contains the green and red IDAT files."
      )
    ),
    
    layout_columns(
      # Card showing the status of loading process
      card(
        card_header("Load Status"),
        verbatimTextOutput(ns("status")) 
      ),
      
      # Card showing the manifest details
      card(
        card_header("Array Manifest"),
        verbatimTextOutput(ns("manifest_output")) # Output for manifest
      )
    )
  )
}


# ─────────────────────────────────────
# SERVER FUNCTION
# ─────────────────────────────────────

loadDataServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    # Reactive string to store and update the step-by-step status messages
    status_log <- reactiveVal("Waiting for user input...")
    
    # Append new message to current status text
    append_status <- function(new_message) {
      current <- status_log()
      updated <- paste(current, new_message, sep = "\n")
      status_log(updated)
    }
    
    # Render status text in UI
    output$status <- renderText({ status_log() })
    
    # Reactive values for results
    RGset <- reactiveVal(NULL)
    raw_normalised <- reactiveVal(NULL)
    targets <- reactiveVal(NULL)
    
    # Observe when "Load Data" button is clicked
    observeEvent(input$load_btn, {
      # Reset reactive values to NULL on each new load attempt
      RGset(NULL)
      raw_normalised(NULL)
      targets(NULL)
      status_log("Waiting for user input...")
      # Clear manifest output on new load attempt
      output$manifest_output <- renderPrint({"Manifest will appear here after successful loading." })
      
      req(input$manual_path)
      RawDataDir <- gsub("\\\\", "/", input$manual_path)
      
      withProgress(message = 'Loading data...', value = 0, {
        
        # Step 1: Validate Path
        if (!dir.exists(RawDataDir)) {
          append_status(paste0("❌ Path does not exist: ", RawDataDir))
          incProgress(amount = 0, detail = "Path error") # Indicate immediate failure, no progress
          return()
        }
        append_status("✅ Step 1: Validating path completed.")
        incProgress(1/6, detail = "Path validated") # Progress: 1/6
        
        
        # Step 2: Reading Sample Sheet
        append_status("Step 2: Reading Sample Sheet...")
        incProgress(1/6, detail = "Reading Sample Sheet") # Progress: 2/6
        targets_val <- tryCatch({
          read.metharray.sheet(RawDataDir)
        }, error = function(e) {
          append_status(paste("❌ Error reading Sample Sheet:", e$message))
          incProgress(amount = 0, detail = "Sample Sheet error") # Indicate failure, no progress
          return(NULL)
        })
        if (is.null(targets_val)) return()
        append_status("✅ Step 2: Sample Sheet read.")
        
        
        # Step 3: Reading IDAT files
        append_status("Step 3: Reading IDAT files...")
        incProgress(1/6, detail = "Reading IDAT files") # Progress: 3/6
        RGset_val <- tryCatch({
          read.metharray.exp(targets = targets_val)
        }, error = function(e) {
          append_status(paste("❌ Error reading IDATs:", e$message))
          incProgress(amount = 0, detail = "IDAT error") # Indicate failure, no progress
          return(NULL)
        })
        if (is.null(RGset_val)) return()
        append_status("✅ Step 3: IDAT files read.")
        
        
        # Step 4: Getting Manifest
        append_status("Step 4: Getting Manifest...")
        incProgress(1/6, detail = "Getting Manifest") # Progress: 4/6
        manifest <- tryCatch({
          getManifest(RGset_val)
        }, error = function(e) {
          append_status(paste("❌ Error getting manifest:", e$message))
          incProgress(amount = 0, detail = "Manifest error") # Indicate failure, no progress
          return(NULL)
        })
        if (is.null(manifest)) return()
        append_status("✅ Step 4: Manifest retrieved.")
        output$manifest_output <- renderPrint({ manifest })
        
        
        # Step 5: Converting raw signals to locus level (preprocessRaw)
        append_status("Step 5: Converting raw signals to locus level...")
        incProgress(1/6, detail = "Preprocessing raw data") # Progress: 5/6
        raw_norm_val <- tryCatch({
          preprocessRaw(RGset_val)
        }, error = function(e) {
          append_status(paste("❌ Error preprocessing:", e$message))
          incProgress(amount = 0, detail = "Preprocessing error") # Indicate failure, no progress
          return(NULL)
        })
        if (is.null(raw_norm_val)) return()
        append_status("✅ Step 5: Preprocessing completed.")
        
        
        # Step 6: Extracting sample info and CpG count
        append_status("Step 6: Extracting sample info and CpG count...")
        incProgress(1/6, detail = "Extracting summary info") # Progress: 6/6 (or 100%)
        sample_names <- sampleNames(RGset_val)
        num_samples <- length(sample_names)
        num_cpgs <- nrow(raw_norm_val)
        
        append_status(paste0(
          "✅ Sample and CpG info retrieved:\n",
          "  - Number of CpGs: ", num_cpgs, "\n",
          "  - Number of Samples: ", num_samples, "\n",
          "  - Sample Names:\n", paste(sample_names, collapse = "\n")
        ))
        
        # Update reactive values for loaded data
        RGset(RGset_val)
        raw_normalised(raw_norm_val)
        targets(targets_val)
        
        # Final progress update: Ensures it hits 100% and then fades.
        incProgress(0, detail = "Done!")
      }) # Closes withProgress block
      
      # Optional: Save data to disk
      # saveRDS(list(RGset = RGset, raw_normalised = raw_normalised, targets = targets), "preprocessed_data.rds")
      #append_status("✅ Step 7: Data saved to 'preprocessed_data.rds'")
      
    }) # Closes observeEvent(input$load_btn, { ... })
    
    # Return a list of reactive values
    return(list(
      RGset = RGset,
      raw_normalised = raw_normalised,
      targets = targets
    ))
  }) # Closes moduleServer(id, function(input, output, session) { ... })
} # Closes loadDataServer <- function(id) { ... }



'# test_loadDataApp.R

library(shiny)
library(bslib)
library(minfi)

ui <- page_navbar(
  title = "EPIC Array Pipeline",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  
  nav_panel("Load Raw Data", loadDataUI("loader")),
  #nav_panel("Another Tab", ...)
)



server <- function(input, output, session) {
  loadDataServer("loader") 
  # other server code or modules
}

shinyApp(ui, server)'