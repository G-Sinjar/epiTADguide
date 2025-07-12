# loadDataModule.R
# Author: Ghazal Sinjar
# Date: 30.05.2025
# Description: This Shiny module allows users to load EPIC array data
# by specifying a path to a directory containing a sample sheet (.csv)
# and corresponding IDAT files. The data is preprocessed and summarized.

library(shiny)
library(bslib)
library(minfi)

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
      helpText("This name will be used to create a subfolder for all project outputs."),
      hr(),
      
      # Project Output Base Folder Input (now for the *base* directory)
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
# loadDataModule.R
# Author: Ghazal Sinjar
# Date: 30.05.2025
# Description: This Shiny module allows users to load EPIC array data
# by specifying a path to a directory containing a sample sheet (.csv)
# and corresponding IDAT files. The data is preprocessed and summarized.

library(shiny)
library(bslib)
library(minfi)

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
      #helpText("This name will be used to create a subfolder for all project outputs."),
      #hr(),
      
      # Project Output Base Folder Input (now for the *base* directory)
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
    ns <- session$ns # <--- Keep this line
    
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
    
    # Reactive values for results
    RGset <- reactiveVal(NULL)
    raw_normalised <- reactiveVal(NULL)
    targets <- reactiveVal(NULL)
    
    # Reactive value for the FINAL project output directory (base + project_name)
    project_output_dir <- reactiveVal(NULL) # NULL initially, only set on successful check
    
    # Reactive value to control visibility of other sidebar elements
    raw_data_input_visible <- reactiveVal(FALSE)
    
    # Render conditional UI elements
    output$conditional_raw_data_input <- renderUI({
      if (raw_data_input_visible()) {
        tagList(
          strong("2. Choose the path to raw data from sequencer"),
          
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
        )
      }
    })
    
    
    # Observe "Check and Set Path" button for project output directory
    observeEvent(input$check_set_project_path_btn, {
      # Check for null fields and show error message
      if (is.null(input$project_name) || nchar(input$project_name) == 0 ||
          is.null(input$project_manual_path) || nchar(input$project_manual_path) == 0) {
        showNotification(
          "Please enter both a Project Name and a Base Folder Path for Outputs before trying again.",
          type = "error",
          duration = 5
        )
        return() # Stop the execution of the observeEvent
      }
      
      base_path <- gsub("\\\\", "/", input$project_manual_path)
      project_name <- make.names(input$project_name) # Sanitize project name for folder
      project_name <- gsub("\\.", "_", project_name) # Replace dots with underscores (make.names uses dots)
      project_name <- gsub(" ", "_", project_name) # Replace spaces with underscores
      # Replace backslashes and forward slashes with underscores
      project_name <- gsub("[\\\\/]", "_", project_name)
      
      if (nchar(project_name) == 0) {
        append_status("❌ Project name cannot be empty or invalid after sanitization.")
        project_output_dir(NULL)
        raw_data_input_visible(FALSE)
        return()
      }
      
      status_log(paste0("Checking and setting project output path for '", project_name, "'..."))
      raw_data_input_visible(FALSE) # Hide raw data input until path is valid
      
      # --- START MODIFIED LOGIC ---
      
      # 1. Check if the BASE PATH exists first
      if (!dir.exists(base_path)) {
        append_status(paste0("❌ Error: The specified Base Folder Path for Outputs does not exist: ", base_path))
        showNotification(
          paste0("Error: The base folder path '", base_path, "' does not exist. Please enter a valid path."),
          type = "error",
          duration = 8
        )
        project_output_dir(NULL)
        raw_data_input_visible(FALSE)
        return() # Stop execution
      } else {
        append_status(paste0("✅ Base Folder Path exists: ", base_path))
      }
      
      # 2. Construct the full, final project output path (base + project_name)
      final_project_path <- file.path(base_path, project_name)
      
      # 3. Now, handle the project-specific subfolder
      if (!dir.exists(final_project_path)) {
        append_status(paste0("Attempting to create project directory: ", final_project_path))
        tryCatch({
          dir.create(final_project_path, recursive = TRUE, showWarnings = TRUE) # recursive=TRUE is still good practice
          append_status("✅ Successfully created project directory")
          project_output_dir(final_project_path)
          append_status("✅ Project output path set successfully.")
          raw_data_input_visible(TRUE) # Show raw data input now
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
          raw_data_input_visible(TRUE) # Show raw data input now
        })
      } else {
        append_status(paste0("✅ Project directory already exists: ", final_project_path))
        project_output_dir(final_project_path)
        append_status("✅ Project output path set successfully.\nNow please enter the path to your raw data from sequencer in the side bar bellow.")
        raw_data_input_visible(TRUE) # Show raw data input now
      }
      # --- END MODIFIED LOGIC ---
    })
    
    
    # Observe when "Load Data" button is clicked (for raw data)
    observeEvent(input$load_btn, {
      
      # Reset reactive values to NULL on each new load attempt
      RGset(NULL)
      raw_normalised(NULL)
      targets(NULL)
      
      status_log(paste0("Loading raw data... (Project output: ", project_output_dir(), ")"))
      # Require both paths to be set
      req(input$manual_path, project_output_dir())
      
      RawDataDir <- gsub("\\\\", "/", input$manual_path)
      
      withProgress(message = 'Loading data...', value = 0, {
        
        # Step 1: Validate Raw Data Path
        if (!dir.exists(RawDataDir)) {
          append_status(paste0("❌ Path does not exist: ", RawDataDir))
          incProgress(amount = 0, detail = "Path error") # Indicate immediate failure, no progress
          showNotification(
            paste0("Error: The raw data path '", RawDataDir, "' does not exist. Please enter a valid path."),
            type = "error",
            duration = 8
          )
          return()
        }
        append_status("✅ Step 1: Validating raw data path completed.")
        #incProgress(1/6, detail = "Path validated") # Progress: 1/6
        
        
        # Step 2: Reading Sample Sheet
        append_status("Step 2: Reading Sample Sheet...")
        incProgress(1/6, detail = "Reading Sample Sheet") # Progress: 2/6
        targets_val <- tryCatch({
          read.metharray.sheet(RawDataDir)
        }, error = function(e) {
          append_status(paste("❌ Error reading Sample Sheet:", e$message))
          incProgress(amount = 0, detail = "Sample Sheet error") # Indicate failure, no progress
          showNotification(
            paste("Error reading Sample Sheet:", e$message),
            type = "error",
            duration = 8
          )
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
          showNotification(
            paste("Error reading IDATs:", e$message),
            type = "error",
            duration = 8
          )
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
          showNotification(
            paste("Error getting manifest:", e$message),
            type = "error",
            duration = 8
          )
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
          showNotification(
            paste("Error preprocessing:", e$message),
            type = "error",
            duration = 8
          )
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
'# ADDED: Project Output Folder Input
strong("Choose Project Output Directory:"),
shinyDirButton(ns("projectPathDir"), "Select Output Folder", "Please select a folder for project outputs"),
verbatimTextOutput(ns("selected_project_path")), # To display the selected path to the user
helpText("All generated reports and plots will be saved here."),
hr(), # Separator'



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


