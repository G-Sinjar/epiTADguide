# modules/9_tadcalling_module.R

# This module provides a user interface for performing TAD (Topologically Associating Domain)
# calling on Hi-C contact maps. It integrates a Python script for initial data processing
# from .mcool files and then executes the deDoc2 Java tool for TAD and subTAD identification.
# The module handles input validation, file management (including copying .mcool files
# to a project directory), execution of external tools, and displays the processed TAD
# and subTAD tables within the Shiny application.

# Author: Ghazal Sinjar
# Date: 29.07.2025

# Load necessary libraries (these should also be in the main app, but good to have here for module testing)
library(shiny)
library(bslib)
library(DT)
library(reticulate)
library(processx)
library(readr)
library(openxlsx)

# Module UI
tadcalling_ui <- function(id) {
  ns <- NS(id) # Create a namespace function for this module's inputs/outputs
  
  tagList(
    layout_sidebar(
      sidebar = sidebar(
        textInput(
          inputId = ns("tissue"), # Namespace input ID
          label = "Tissue Type",
          placeholder = "CAKI2"
        ),
        helpText("Tissue of your mcool file/ The source tissue of the hic-map."),
        
        textInput(
          inputId = ns("chromosome"), # Namespace input ID
          label = "Chromosome",
          placeholder = "ChrX, chrx or x"
        ),
        
        numericInput(
          inputId = ns("resolution"), # Namespace input ID
          label = "Hi-C Map Resolution (in kb)",
          value = 25,
          min = 1,
          step = 1
        ),
        
        textInput(
          inputId = ns("mcool_path"), # Namespace input ID
          label = "Path to .mcool file",
          placeholder = "C:\\path\\to\\your\\file.mcool",
          width = "100%"
        ),
        
        selectInput(
          inputId = ns("java_memory"), # Namespace input ID
          label = "Java Memory Allocation",
          choices = c("4g", "8g", "16g"),
          selected = "8g"
        ),
        helpText("💡 Don't choose more than the limit of your RAM. Recommended: use the maximum available RAM for better performance, as this process is demanding."),
        
        actionButton(ns("start_tadcalling"), "🚀 Start TADcalling", class = "btn-primary"), # Namespace button ID
        
        hr(),
        
        radioButtons(
          inputId = ns("table_choice"), # Namespace input ID
          label = "Table to Display",
          choices = c("TADs table", "SubTADs table"),
          selected = "TADs table"
        )
      ),
      
      div(class = "col-lg-9 col-md-9 col-sm-12",
          uiOutput(ns("status_card")), # Namespace output ID
          br(),
          DT::dataTableOutput(ns("tad_table")) # Namespace output ID
      )
    )
  )
}

# Module Server
# ─────────────────────────────────────
# SERVER FUNCTION
# ─────────────────────────────────────
#' @param id Module ID.
#' @param project_output_dir A reactive expression for the project's main output directory path.
#' @return A list of reactive expressions containing the processed TADs and SubTADs data frames.
tadcalling_server <- function(id,project_output_dir = NULL) {
  moduleServer(id, function(input, output, session) {
    
    # Set up Python environment (assuming this path is globally valid or relative to app root)
    use_virtualenv("shiny_py", required = TRUE)
    source_python("./python_scripts/1-Hi-c_maps_exploring.py")
    
    # Reactive values to store status messages and generated dataframes
    status_msg <- reactiveVal(
      "📥 Enter the required inputs on the left and click 'Start TADcalling' to begin.<br>The resolution refers to the bin size used for the contact map, e.g. 25 kb means each matrix bin spans 25,000 base pairs."
    )
    # Reactive values to store the TAD and SubTAD data frames
    # These will be returned by the module
    tads_df_rv <- reactiveVal(NULL)
    subtads_df_rv <- reactiveVal(NULL)
    
    # Reactive values to persist paths and parameters after a successful run
    results_root_folder_rv <- reactiveVal(NULL)
    processed_tads_output_dir_rv <- reactiveVal(NULL)
    java_raw_output_dir_rv <- reactiveVal(NULL)
    current_tissue_rv <- reactiveVal(NULL)
    current_chrom_rv <- reactiveVal(NULL)
    current_resolution_rv <- reactiveVal(NULL)
    current_java_memory_rv <- reactiveVal(NULL) 
    new_mcool_path_rv <- reactiveVal(NULL) 
    
    output$status_card <- renderUI({
      msg <- status_msg()
      card_color <- if (grepl("❌", msg) || grepl("Error", msg)) {
        "bg-danger text-white"
      } else {
        "bg-light"
      }
      div(
        class = paste("card shadow-sm mb-4 p-3", card_color),
        h4("🧬 Analysis Status"),
        HTML(gsub("\n", "<br>", msg))
      )
    })
    
    selected_table <- reactiveVal("TAD")
    observeEvent(input$table_choice, {
      selected_table(ifelse(input$table_choice == "SubTADs table", "SubTAD", "TAD"))
    })
    
    # --- Main Analysis Logic ---
    observeEvent(input$start_tadcalling, {
      # Wrap the entire analysis process in withProgress
      withProgress(message = 'TAD Calling in Progress', value = 0, {
        
        # Capture current input values immediately
        current_tissue <- input$tissue
        current_chromosome <- input$chromosome
        current_resolution <- input$resolution
        current_mcool_path <- input$mcool_path
        current_java_memory <- input$java_memory
        
        
        # Input Validation - use the captured current_ values
        missing_inputs <- c()
        if (is.null(current_tissue) || current_tissue == "") {
          missing_inputs <- c(missing_inputs, "Tissue Type")
        }
        if (is.null(current_chromosome) || current_chromosome == "") {
          missing_inputs <- c(missing_inputs, "Chromosome")
        }
        if (is.null(current_mcool_path) || current_mcool_path == "") {
          missing_inputs <- c(missing_inputs, "Path to .mcool file")
        }
        
        if (length(missing_inputs) > 0) {
          status_msg(paste0("❌ Error: Please enter values for the following required inputs: ",
                            paste(missing_inputs, collapse = ", "),
                            ". Then click 'Start TADcalling' again."))
          return(NULL)
        }
        
        tissue_val <- current_tissue
        chromosome_val <- current_chromosome
        resolution_val <- current_resolution
        mcool_path_val <- current_mcool_path
        java_memory_val <- current_java_memory
        print(mcool_path_val)
        # 1. Normalize chromosome input
        chrom_temp <- gsub("^chr", "", chromosome_val, ignore.case = TRUE)
        
        # Now, handle specific cases (X, Y, M) and general numeric chromosomes
        if (toupper(chrom_temp) %in% c("X", "Y", "M")) {
          chrom <- paste0("chr", toupper(chrom_temp)) # Keep X, Y, M uppercase: chrX, chrY, chrM
        } else {
          chrom <- paste0("chr", tolower(chrom_temp)) # For numeric: chr1, chr2, etc. (ensures consistent lowercase numbers)
        }
        
        
        # 2. Normalize mcool path
        base_project_dir_val <- if (!is.null(project_output_dir)) project_output_dir() else NULL
        
        # Ensure project_output_dir is an absolute path for consistency
        if (!is.null(base_project_dir_val) && base_project_dir_val != "") {
          base_project_dir_val <- normalizePath(base_project_dir_val, winslash = "/", mustWork = FALSE)
        }
        
        # Normalize mcool_path_val: Ensure it has a .mcool extension
        if (!grepl("\\.mcool$", mcool_path_val, ignore.case = TRUE)) {
          mcool_path_val <- paste0(mcool_path_val, ".mcool")
          status_msg(paste0(status_msg(),"<br>ℹ️ Appended '.mcool' extension to input path: "))
        }

        # Determine the final mcool_path that the Python script will use
        mcool_path <- "" 
        
        if (!is.null(base_project_dir_val) && base_project_dir_val != "") {
          # Define the base directory for all TADcaller results within the project output
          tadcaller_base_output_dir <- file.path(base_project_dir_val, "TADcaller_results")
          
          # Extract filename and filename without extension
          mcool_filename <- basename(mcool_path_val)
          print(mcool_filename) # Debug print
          mcool_filename_no_ext <- tools::file_path_sans_ext(mcool_filename)
          print(mcool_filename_no_ext) # Debug print
          
          # Create the specific subfolder named after the mcool file
          specific_mcool_output_dir <- file.path(tadcaller_base_output_dir, mcool_filename_no_ext)
          print(specific_mcool_output_dir) # Debug print
          
          # Ensure this new specific directory exists
          if (!dir.exists(specific_mcool_output_dir)) {
            dir.create(specific_mcool_output_dir, recursive = TRUE)
            status_msg(paste0(status_msg(),"<br>📁 Created directory for .mcool file: ", specific_mcool_output_dir))
            print(paste(specific_mcool_output_dir, "doesnt excists so it is created")) # Debug print
          }
          
          # Define the full new path for the copied mcool file
          new_mcool_path_for_copy <- file.path(specific_mcool_output_dir, mcool_filename) # Corrected variable name here
          print(new_mcool_path_for_copy) # Debug print
          
          # --- DEBUGGING FILE.COPY ---
          print(paste("DEBUG: Source mcool_path_val (after extension check):", mcool_path_val))
          print(paste("DEBUG: Does source file exist?", file.exists(mcool_path_val)))
          print(paste("DEBUG: Destination for copy:", new_mcool_path_for_copy))
          print(paste("DEBUG: Does destination directory exist?", dir.exists(dirname(new_mcool_path_for_copy))))
          
          copy_success <- FALSE
          if (file.exists(mcool_path_val)) {
            # Ensure 'to' argument is correct here. It should be new_mcool_path_for_copy
            copy_success <- file.copy(from = mcool_path_val, to = new_mcool_path_for_copy, overwrite = TRUE)
          } else {
            status_msg(paste0(status_msg(),"<br>❌ Error: Original .mcool file not found at: ", mcool_path_val))
            return(NULL) # Stop execution if source file is missing
          }
          
          if (copy_success) {
            status_msg(paste0(status_msg(),"<br>🔄 Copied .mcool file to project output directory: ", new_mcool_path_for_copy))
            mcool_path <- gsub("\\\\", "/", new_mcool_path_for_copy) # Update mcool_path to the copied file
            new_mcool_path_rv(mcool_path)
          } else {
            status_msg(paste0(status_msg(),"<br>❌ Error: Failed to copy .mcool file from '", mcool_path_val, "' to '", new_mcool_path_for_copy, "'. Check permissions and paths."))
            return(NULL) # Stop execution if copy failed
          }
          # --- END DEBUGGING FILE.COPY ---
          
        } else {
          # If project_output_dir is null, use the original mcool_path and normalize slashes
          mcool_path <- gsub("\\\\", "/", mcool_path_val)
          status_msg(paste0(status_msg(),"<br>ℹ️ Using original .mcool file path: ", mcool_path))
          new_mcool_path_rv(mcool_path)
        }
        
        # Print the final mcool_path that will be passed to Python
        print(paste("DEBUG: Final mcool_path passed to Python:", mcool_path))
        print(paste("DEBUG: Does final mcool_path exist?", file.exists(mcool_path)))
        
        incProgress(0.1, detail = "Starting Step 1: Processing .mcool file")
        status_msg("Step 1: Processing Hi-c map (.mcool file) to extract contact matrix of the chosen chr and resolution...")
        # Step 1: Process .mcool file using Python
        python_step_success <- FALSE
        result <- NULL
        try({
          result <- process_contact_matrix(
            tissue = tissue_val,
            chrom = chrom,
            resolution_kb = resolution_val,
            path_mcool_file = mcool_path
          )
          python_step_success <- result$success
        }, silent = TRUE)
        
        if (is.null(result) || !python_step_success) {
          error_message <- if (!is.null(result) && !is.null(result$error) && !is.null(result$message)) {
            paste0("❌ Error in Step 1: ", result$error, "<br>", result$message)
          } else {
            "❌ Unknown error occurred during Python Step 1. Check console for Python errors."
          }
          status_msg(paste0(status_msg(),"<br>", error_message))
          return(NULL)
        }
        
        status_msg(paste0(status_msg(),"<br>✅ Step 1 completed.<br>",
                          "Total contacts: in whole map ", format(result$total_contacts, big.mark = ",")))
        
        # Define results folders and update reactive values
        base_dir_of_mcool <- dirname(mcool_path)
        if (is.null(base_dir_of_mcool) || base_dir_of_mcool == "") {
          base_dir_of_mcool <- getwd()
        }
        results_root_folder <- file.path(base_dir_of_mcool, "TADcaller_Results", paste0("TADs_", tissue_val))
        java_raw_output_dir <- file.path(results_root_folder, "java_raw_output")
        processed_tads_output_dir <- file.path(results_root_folder, "processed_tads")
        
        # Update reactive values for paths and current parameters
        results_root_folder_rv(results_root_folder)
        java_raw_output_dir_rv(java_raw_output_dir)
        processed_tads_output_dir_rv(processed_tads_output_dir)
        current_tissue_rv(tissue_val)
        current_chrom_rv(chrom)
        current_resolution_rv(resolution_val)
        current_java_memory_rv(java_memory_val) 
        
        if (!dir.exists(results_root_folder)) dir.create(results_root_folder, recursive = TRUE)
        if (!dir.exists(java_raw_output_dir)) dir.create(java_raw_output_dir, recursive = TRUE)
        if (!dir.exists(processed_tads_output_dir)) dir.create(processed_tads_output_dir, recursive = TRUE)
        
        incProgress(0.4, detail = "Step 2: Running TADcaller deDoc2")
        
        # Step 2: Run TADcaller (Java deDoc2)
        status_msg(paste0(status_msg(),"<br>  Step 2: Running TADcaller deDoc2..."))
        
        input_matrix <- result$output_file
        if (is.null(input_matrix) || !file.exists(input_matrix)) {
          status_msg(paste0(status_msg(),"<br>❌ Input matrix file not found. Ensure Step 1 completed correctly and generated the file."))
          return(NULL)
        }
        
        jar_path <- "../PythonProject/deDoc2-main/deDoc2.jar" # Relative to app's working directory
        if (!file.exists(jar_path)) {
          status_msg(paste0(status_msg(),"<br>❌ .jar file not found at: ", normalizePath(jar_path, mustWork = FALSE), ". Please ensure the JAR file is in the correct location and accessible."))
          return(NULL)
        }
        
        output_base <- file.path(java_raw_output_dir, paste0(tissue_val, "_", chrom, "_", resolution_val, "kb_TADlines"))
        
        java_args <- c(
          paste0("-Xms", java_memory_val),
          "-jar", jar_path,
          "-inputfile", input_matrix,
          "-binsize", resolution_val,
          "-outputfile", output_base
        )
        
        result_java <- tryCatch({
          processx::run("java", args = java_args, echo = TRUE, error_on_status = FALSE)
        }, error = function(e) {
          list(status = 1, stdout = "", stderr = as.character(e))
        })
        
        if (result_java$status != 0) {
          err_msg <- result_java$stderr
          if (grepl("Unable to access jarfile", err_msg)) {
            status_msg(paste0(status_msg(),"<br>❌ Java error: Unable to access .jar file. Check the path and file permissions."))
          } else {
            status_msg(paste0(status_msg(),"<br>❌ Java TADcaller error:\n", err_msg))
          }
          return(NULL)
        }
        
        status_msg(paste0(status_msg(),"<br>✅ Step 2 completed"))
        
        incProgress(0.7, detail = "Step 3: Processing TAD results")
        
        # Step 3: Post-process TAD results
        status_msg(paste0(status_msg(),"<br>Step 3: Processing TAD results..."))
        
        tryCatch({
          processed_tads_df <- process_tad_results( # This function is sourced from TADcalling_utils.R
            chrom = chrom,
            tissue = tissue_val,
            binsize = resolution_val,
            input_dir = java_raw_output_dir,
            output_dir = processed_tads_output_dir
          )
          
          if (is.null(processed_tads_df) || !is.data.frame(processed_tads_df) || nrow(processed_tads_df) == 0) {
            status_msg(paste0(status_msg(),"<br>❌ Error: No TADs were processed or `process_tad_results` returned an empty data frame. Check the Java raw output files."))
            tads_df_rv(NULL) # Set reactive value to NULL if no data
            subtads_df_rv(NULL) # Also clear subtads if main TADs failed
            return(NULL)
          }
          
          # Update the reactive value for TADs
          tads_df_rv(processed_tads_df)
          
          processed_subtads_df <- process_subtad_results( # This function is sourced from TADcalling_utils.R
            chrom = chrom,
            tissue = tissue_val,
            binsize = resolution_val,
            input_dir = java_raw_output_dir,
            output_dir = processed_tads_output_dir,
            main_tads_df_input = processed_tads_df # Pass the actual TADs dataframe
          )
          
          # Update the reactive value for SubTADs
          subtads_df_rv(processed_subtads_df) # Will be NULL if subtad processing failed
          
          status_msg(paste0(status_msg(),"<br>✅ Step 3 completed"))
        }, error = function(e) {
          status_msg(paste0(status_msg(),"<br>❌ Error in Step 3 (Post-processing): ", as.character(e)))
          tads_df_rv(NULL) # Clear reactive value on error
          subtads_df_rv(NULL) # Clear reactive value on error
          return(NULL)
        })
        
        status_msg(paste0(status_msg(), "<br>📁 All results are saved under: '", normalizePath(results_root_folder, winslash = "/"), "'"))
        
        incProgress(1, detail = "Analysis Complete!")
      }) # End withProgress
      
      # Display the table based on current selection after successful run
      output$tad_table <- DT::renderDataTable({
        req(processed_tads_output_dir_rv()) # Ensure path is set
        req(current_tissue_rv(), current_chrom_rv(), current_resolution_rv()) # Ensure parameters are set
        
        table_path <- if (selected_table() == "TAD") {
          file.path(processed_tads_output_dir_rv(), sprintf("%s_%s_%skb_TADs.txt", current_tissue_rv(), current_chrom_rv(), current_resolution_rv()))
        } else {
          file.path(processed_tads_output_dir_rv(), sprintf("%s_%s_%skb_SubTADs.txt", current_tissue_rv(), current_chrom_rv(), current_resolution_rv()))
        }
        
        print(paste("R (Debug): Attempting to read table from:", table_path))
        
        if (!file.exists(table_path)) {
          current_selection <- if (selected_table() == "TAD") "TADs" else "SubTADs"
          status_msg(paste0(status_msg(), "<br>⚠️ Warning: The requested '", current_selection, "' table file not found: '", table_path, "'"))
          return(NULL)
        }
        
        tryCatch({
          data <- readr::read_delim(table_path, delim = "\t", show_col_types = FALSE) # Use readr::read_delim for tsv
          print(paste("R (Debug): Successfully read data with dimensions:", paste(dim(data), collapse = "x")))
          DT::datatable(data, options = list(pageLength = 10)) # Added pageLength option
        }, error = function(e) {
          status_msg(paste0(status_msg(), "<br>❌ Error displaying table: Could not read or process the data from '", basename(table_path), "'. Error: ", as.character(e)))
          return(NULL)
        })
      })
    })
    
    # Observer to update the table when the radio button choice changes
    observeEvent(input$table_choice, {
      req(processed_tads_output_dir_rv()) # Only proceed if an analysis has successfully run
      req(current_tissue_rv(), current_chrom_rv(), current_resolution_rv())
      
      output$tad_table <- DT::renderDataTable({
        table_path <- if (selected_table() == "TAD") {
          file.path(processed_tads_output_dir_rv(), sprintf("%s_%s_%skb_TADs.txt", current_tissue_rv(), current_chrom_rv(), current_resolution_rv()))
        } else {
          file.path(processed_tads_output_dir_rv(), sprintf("%s_%s_%skb_SubTADs.txt", current_tissue_rv(), current_chrom_rv(), current_resolution_rv()))
        }
        
        print(paste("R (Debug): (table_choice) Attempting to read table from:", table_path))
        
        if (!file.exists(table_path)) {
          current_selection <- if (selected_table() == "TAD") "TADs" else "SubTADs"
          status_msg(paste0(status_msg(), "<br>⚠️ Warning: The requested '", current_selection, "' table file not found: '", table_path, "'. Perhaps the analysis hasn't run yet or failed."))
          return(NULL)
        }
        tryCatch({
          data <- readr::read_delim(table_path, delim = "\t", show_col_types = FALSE) # Use readr::read_delim for tsv
          print(paste("R (Debug): (table_choice) Successfully read data with dimensions:", paste(dim(data), collapse = "x")))
          DT::datatable(data, options = list(pageLength = 10))
        }, error = function(e) {
          status_msg(paste0(status_msg(), "<br>❌ Error displaying table from table_choice: ", as.character(e)))
          return(NULL)
        })
      })
    }, ignoreInit = TRUE)
    
    # Return a list of reactive values that the parent app can access
    return(list(
      tads_table = tads_df_rv,
      subtads_table = subtads_df_rv,
      current_tissue = current_tissue_rv,
      current_chrom = current_chrom_rv,
      current_resolution = current_resolution_rv,
      processed_output_dir = processed_tads_output_dir_rv, 
      java_memory_used = current_java_memory_rv,             
      copied_mcool_path = new_mcool_path_rv
      #status_message = status_msg # Optionally return status if parent needs to monitor
    ))
  })
}

#the test code of this module is in the R_scrips called 09-test.R