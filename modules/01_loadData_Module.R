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
  ns <- NS(id)
  
  layout_sidebar(
    sidebar = sidebar(
      # ADDED: Project Name Input
      strong("1. Set Project Name and Output Directory"),
      textInput(
        inputId = ns("project_name"),
        label = "Enter Project Name:",
        value = "",
        placeholder = "E.g. My_EPIC_Study"
      ),
      textInput(
        inputId = ns("project_manual_path"),
        label = "Enter Base Folder Path for Outputs:",
        value = "",
        placeholder = "E.g. C:\\Users\\yourname\\Documents\\Shiny_Projects"
      ),
      actionButton(ns("check_set_project_path_btn"), "Check and Set Project Path"),
      helpText("All generated reports and plots will be saved in a subfolder named after your project within this base path."),
      hr(),
      
      # These elements will be conditionally shown
      uiOutput(ns("conditional_raw_data_input"))
    ),
    
    layout_columns(
      # Card showing the status of loading process
      card(
        card_header("Load Status"),
        verbatimTextOutput(ns("status"))
      ),
      
      # Card showing the array infos details
      card(
        card_header("Array information"),
        uiOutput(ns("arry_info_ui")) # Changed to uiOutput to allow for conditional display of a table
      )
    )
  )
}

# ─────────────────────────────────────
# SERVER FUNCTION
# ─────────────────────────────────────
loadDataServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns 
    
    # Reactive string to store and update the step-by-step status messages
    status_log <- reactiveVal("Waiting for user to set project name and output path...")
    
    # Append new message to current status text
    append_status <- function(new_message) {
      current <- status_log()
      updated <- paste(current, new_message, sep = "\n")
      status_log(updated)
    }
    
    # Render status text in UI
    output$status <- renderText({ status_log() })
    
    #-----------------------------------------------
    # Reactive values for results
    RGset <- reactiveVal(NULL)
    raw_normalised <- reactiveVal(NULL)
    targets <- reactiveVal(NULL)
    
    # Reactive value for the FINAL project output directory (base + project_name)
    project_output_dir <- reactiveVal(NULL)
    
    #--------------------------------------------------------
    # Reactive value to control visibility of other sidebar elements
    raw_data_input_visible <- reactiveVal(FALSE)
    
    # Render conditional UI elements
    output$conditional_raw_data_input <- renderUI({
      if (raw_data_input_visible()) {
        tagList(
          strong("2. Choose the path to raw data from sequencer"),
          textInput(
            inputId = ns("manual_path"), 
            label = "Enter folder path manually:",
            value = "",
            placeholder = "E.g. C:\\Users\\yourname\\Documents\\sequencer_output"
          ),
          actionButton(ns("load_btn"), "Load Data"),
          helpText(
            "Please select a folder that includes: (1) a CSV file with slide and sample information, ",
            "and (2) a subfolder named after the slide, which contains the green and red IDAT files."
          )
        )
      }
    })
    
    # Observe "Check and Set Path" button for project output directory
    observeEvent(input$check_set_project_path_btn, {
      if (is.null(input$project_name) || nchar(input$project_name) == 0 ||
          is.null(input$project_manual_path) || nchar(input$project_manual_path) == 0) {
        showNotification(
          "Please enter both a Project Name and a Base Folder Path for Outputs before trying again.",
          type = "error",
          duration = 5
        )
        return()
      }
      
      base_path <- gsub("\\\\", "/", input$project_manual_path)
      project_name <- make.names(input$project_name)
      project_name <- gsub("\\.", "_", project_name)
      project_name <- gsub(" ", "_", project_name)
      project_name <- gsub("[\\\\/]", "_", project_name)
      
      if (nchar(project_name) == 0) {
        append_status("❌ Project name cannot be empty or invalid after sanitization.")
        project_output_dir(NULL)
        raw_data_input_visible(FALSE)
        return()
      }
      
      status_log(paste0("Checking and setting project output path for '", project_name, "'..."))
      raw_data_input_visible(FALSE)
      
      if (!dir.exists(base_path)) {
        append_status(paste0("❌ Error: The specified Base Folder Path for Outputs does not exist: ", base_path))
        showNotification(
          paste0("Error: The base folder path '", base_path, "' does not exist. Please enter a valid path."),
          type = "error",
          duration = 8
        )
        project_output_dir(NULL)
        raw_data_input_visible(FALSE)
        return()
      } else {
        append_status(paste0("✅ Base Folder Path exists: ", base_path))
      }
      
      final_project_path <- file.path(base_path, project_name)
      
      if (!dir.exists(final_project_path)) {
        tryCatch({
          dir.create(final_project_path, recursive = TRUE, showWarnings = TRUE)
          append_status("✅ Successfully created project directory")
          project_output_dir(final_project_path)
          append_status(paste0("✅ Project output path set successfully.", final_project_path))
          raw_data_input_visible(TRUE)
        }, error = function(e) {
          append_status(paste("❌ Error creating project directory:", e$message))
          showNotification(
            paste("Error creating project directory:", e$message),
            type = "error",
            duration = 8
          )
          project_output_dir(NULL)
          raw_data_input_visible(FALSE)
        }, warning = function(w) {
          append_status(paste("⚠️ Warning creating project directory:", w$message))
          showNotification(
            paste("Warning creating project directory:", w$message),
            type = "warning",
            duration = 8
          )
          project_output_dir(final_project_path)
          append_status("✅ Project output path set successfully (with warning).")
          raw_data_input_visible(TRUE)
        })
      } else {
        append_status(paste0("✅ Project directory already exists: ", final_project_path))
        project_output_dir(final_project_path)
        append_status("✅ Project output path set successfully.\nNow please enter the path to your raw data from sequencer in the side bar bellow.")
        raw_data_input_visible(TRUE)
      }
    })
    
    # Observe when "Load Data" button is clicked (for raw data)
    observeEvent(input$load_btn, {
      
      RGset(NULL)
      raw_normalised(NULL)
      targets(NULL)
      
      # Reset and start new status log
      status_log(paste0("Step 1: Validating raw data path: ", input$manual_path, "..."))
      
      req(input$manual_path, project_output_dir())
      
      RawDataDir <- gsub("\\\\", "/", input$manual_path)
      
      withProgress(message = 'Loading data...', value = 0, {
        
        # Step 1: Validate Raw Data Path
        if (!dir.exists(RawDataDir)) {
          append_status(paste0("❌ Path does not exist: ", RawDataDir))
          incProgress(amount = 0, detail = "Path error")
          showNotification(
            paste0("Error: The raw data path '", RawDataDir, "' does not exist. Please enter a valid path."),
            type = "error",
            duration = 8
          )
          return()
        }
        append_status("✅ Step 1: completed. ")
        
        # Step 2: Reading Sample Sheet
        append_status("Step 2: Reading Sample Sheet...")
        incProgress(1/6, detail = "Reading Sample Sheet")
        targets_val <- tryCatch({
          read.metharray.sheet(RawDataDir)
        }, error = function(e) {
          append_status(paste("❌ Error reading Sample Sheet:", e$message))
          incProgress(amount = 0, detail = "Sample Sheet error")
          showNotification(
            paste("Error reading Sample Sheet:", e$message),
            type = "error",
            duration = 8
          )
          return(NULL)
        })
        if (is.null(targets_val)) return()
        append_status("✅ Step 2: completed.")
        
        # Step 3: Reading IDAT files
        append_status("Step 3: Reading IDAT files...")
        incProgress(1/6, detail = "Reading IDAT files")
        RGset_val <- tryCatch({
          read.metharray.exp(targets = targets_val)
        }, error = function(e) {
          append_status(paste("❌ Error reading IDATs:", e$message))
          incProgress(amount = 0, detail = "IDAT error")
          showNotification(
            paste("Error reading IDATs:", e$message),
            type = "error",
            duration = 8
          )
          return(NULL)
        })
        if (is.null(RGset_val)) return()
        append_status("✅ Step 3: completed.")
        
        # Step 4: Getting Array Information
        append_status("Step 4: Getting Array Information...")
        incProgress(1/6, detail = "Getting array infos")
        arry_info <- tryCatch({
          annotation(RGset_val)
        }, error = function(e) {
          append_status(paste("❌ Error getting array informations:", e$message))
          incProgress(amount = 0, detail = "Array information error")
          showNotification(
            paste("Error getting array information:", e$message),
            type = "error",
            duration = 8
          )
          return(NULL)
        })
        if (is.null(arry_info)) return()
        
        # Format the array information for display
        array_type <- ifelse("array" %in% names(arry_info), arry_info["array"], "Unknown")
        array_annotation <- ifelse("annotation" %in% names(arry_info), arry_info["annotation"], "Unknown")
        
        append_status(paste0("✅ Step 4: completed."))
        
        # Step 5: Converting raw signals to locus level (preprocessRaw)
        append_status("Step 5: Converting raw signals to locus level...")
        incProgress(1/6, detail = "Preprocessing raw data")
        raw_norm_val <- tryCatch({
          preprocessRaw(RGset_val)
        }, error = function(e) {
          append_status(paste("❌ Error preprocessing:", e$message))
          incProgress(amount = 0, detail = "Preprocessing error")
          showNotification(
            paste("Error preprocessing:", e$message),
            type = "error",
            duration = 8
          )
          return(NULL)
        })
        if (is.null(raw_norm_val)) return()
        append_status("✅ Step 5: completed.")
        
        # Step 6: Extracting sample info and CpG count
        append_status("Step 6: Extracting sample and CpG info...")
        incProgress(1/7, detail = "Extracting summary info")
        
        sample_slide_array <- sampleNames(RGset_val)
        num_samples <- length(sample_slide_array)
        num_cpgs <- nrow(raw_norm_val)
        sample_name <- targets_val$Sample_Name
        sample_group <- targets_val$Sample_Group
        
        # Create a data frame for display in the "Array information" card
        sample_info_df <- data.frame(
          Slide_Array = sample_slide_array,
          Sample_Name = sample_name,
          Sample_Group = sample_group
        )
        
        # Update output$arry_info_ui to display the table
        output$arry_info_ui <- renderUI({
          tagList(
            renderTable(sample_info_df)
          )
        })
        
        # Final, combined status summary message
        final_summary <- paste0(
          "\n✅ All steps completed! Data is ready.\n",
          "  - Array Type: ", array_type, "\n",
          "  - Annotation: ", array_annotation, "\n",
          "  - Number of CpGs: ", num_cpgs, "\n",
          "  - Number of Samples: ", num_samples,
          "\n\n⚠️ Important Note:\n",
          "If the sample names and sample groups aren't correct, please check your sample sheet:\n",
          "1. Ensure sample names are in a column named exactly 'Sample_Name'\n",
          "2. Ensure sample groups are in a column named exactly 'Sample_Group'\n",
          "3. Avoid editing the sample sheet in Excel as it may alter the file structure\n",
          "4. Make corrections and restart the analysis"
        )
        append_status(final_summary)
        
        RGset(RGset_val)
        raw_normalised(raw_norm_val)
        targets(targets_val)
        
        incProgress(0, detail = "Done!")
        # New Step: Save intermediate data
        append_status("Step 7: Saving intermediate data...")
        incProgress(1/7, detail = "Saving data")
        
        intermediate_dir <- file.path(project_output_dir(), "intermediate_data")
        if (!dir.exists(intermediate_dir)) {
          dir.create(intermediate_dir, recursive = TRUE)
        }
        
        # Define file paths
        rgset_filepath <- file.path(intermediate_dir, "RGset.rds")

        # Save the objects
        tryCatch({
          saveRDS(RGset_val, file = rgset_filepath)
          append_status(paste0("✅ Intermediate data saved successfully to ", intermediate_dir, "."))
        }, error = function(e) {
          append_status(paste0("❌ Error saving intermediate data: ", e$message))
          showNotification("Error saving intermediate data.", type = "error", duration = 8)
        })
      }) # Closes withProgress block
      
    }) # Closes observeEvent(input$load_btn, { ... })
    
    # Return a list of reactive values, including the new project_output_dir
    return(list(
      RGset = RGset,
      raw_normalised = raw_normalised,
      targets = targets,
      project_dir = project_output_dir
    ))
  })
}



'# test_loadDataApp.R
library(shiny)
library(bslib)
library(minfi)

ui <- page_navbar(
  title = "EPIC Array Pipeline",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  
  nav_panel("Load Raw Data", loadDataUI("loader")),
)

server <- function(input, output, session) {
  loadDataServer("loader") 
}

shinyApp(ui, server)'