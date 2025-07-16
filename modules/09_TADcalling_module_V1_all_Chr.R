# modules/09_TADcalling_module_V1_all_Chr.R
# Date: 16-07-2025
# Author: Ghazal Sinjar

# Load necessary libraries (these should also be in the main app, but good to have here for module testing)
library(shiny)
library(bslib)
library(DT)
library(reticulate)
library(processx)
library(readr)
library(openxlsx)

# ─────────────────────────────────────
# User Interface (UI) FUNCTION for the Module
# ─────────────────────────────────────
tadcalling_ui <- function(id) {
  ns <- NS(id) # Create a namespace function for this module's inputs/outputs
  
  tagList(
    layout_sidebar(
      sidebar = sidebar(
        textInput(
          inputId = ns("tissue"),
          label = "Tissue Type",
          placeholder = "CAKI2"
        ),
        helpText("Tissue of your mcool file/ The source tissue of the hic-map."),
        
        selectInput(
          inputId = ns("resolution"),
          label = "Hi-C Map Resolution (in kb)",
          choices = c(5, 10, 25, 50, 100, 250, 500, 1000), # Hardcoded choices
          selected = 25,
          width = "100%"
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
        
        actionButton(ns("start_tadcalling"), "🚀 Start TADcalling for ALL Chromosomes", class = "btn-primary"),
        
        hr(), # Separator
        
        uiOutput(ns("chromosome_selector_ui")), # This should be fine here
        
        radioButtons(
          inputId = ns("table_choice"),
          label = "Table to Display for selected Chr:",
          choices = c("TADs table", "SubTADs table"),
          selected = "TADs table"
        ),
        helpText("Select one of the processed chromosomes and its type of table (TADs or SubTADs) to display.")
        
      ), # End sidebar
      
      # Main content area - This div MUST wrap ALL content in the main area
      div( 
        style = "padding-left: 15px; padding-right: 15px;", # Apply padding to the entire main content block
        
        # Status card - now correctly contained within its own layout_columns
        layout_columns(
          col_widths = c(12), # This specifies one column that takes up all 12 Bootstrap columns
          fill = TRUE,
          card(
            uiOutput(ns("status_card"))
          )
        ),
        # The br() (line break) was commented out, so I'll keep it out.
        # Table Section - Now a direct sibling to the status card's layout_columns
        div(
          style = "margin-top: 20px;", # Add some space above the table section
          uiOutput(ns("table_title")), # Now correctly separated by a comma
          DT::dataTableOutput(ns("tad_table"))
        ) 
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
tadcalling_server <- function(id, project_output_dir = NULL) {
  moduleServer(id, function(input, output, session) {
    
    # Set up Python environment
    use_virtualenv("shiny_py", required = TRUE)
    source_python("./python_scripts/process_contact_matrix_all_chroms_function.py")
    
    # Reactive values to store status messages and generated dataframes
    status_msg <- reactiveVal(
      "📥 Enter the path to your .mcool file and other details. Then, click 'Start TADcalling' to begin.<br>This will process all chromosomes in the .mcool file."
    )
    
    # Reactive values to persist paths and parameters after a successful run
    results_root_folder_rv <- reactiveVal(NULL)
    processed_tads_output_dir_rv <- reactiveVal(NULL)
    java_raw_output_dir_rv <- reactiveVal(NULL)
    current_tissue_rv <- reactiveVal(NULL)
    current_resolution_rv <- reactiveVal(NULL)
    current_java_memory_rv <- reactiveVal(NULL)

    # Store a list of successfully processed chromosomes for display validation
    processed_chroms_list_rv <- reactiveVal(NULL) # Stores all successfully processed chroms
    
    #--------------------------------------------
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
    #--------------------------------------------------
    selected_table <- reactiveVal("TAD")
    observeEvent(input$table_choice, {
      selected_table(ifelse(input$table_choice == "SubTADs table", "SubTAD", "TAD"))
    })
    #-------------------------------------------
    # Reactive for the dynamic table title
    output$table_title <- renderUI({
      # Require necessary reactives to be available
      req(input$selected_chrom, current_tissue_rv(), current_resolution_rv(), selected_table())
      
      chrom <- input$selected_chrom
      table_type <- selected_table() # "TAD" or "SubTAD"
      tissue <- current_tissue_rv()
      resolution <- current_resolution_rv()
      
      # Construct the title
      title_text <- sprintf("Processed %s Data for %s (%s, %skb Resolution)",
                            table_type, chrom, tissue, resolution)
      
      h4(title_text)
    })
    #---------------------------------------
    # --- Main Analysis Logic ---
    observeEvent(input$start_tadcalling, {
      # Wrap the entire analysis process in withProgress
      withProgress(message = 'TAD Calling in Progress', value = 0, {
        
        #-------------------------------------------
        # step1: check if there is missing values in input
        # Capture current input values immediately
        current_tissue <- input$tissue
        current_resolution <- as.numeric(input$resolution)
        current_mcool_path <- input$mcool_path
        current_java_memory <- input$java_memory
        
        # Input Validation - use the captured current_ values
        missing_inputs <- c()
        if (is.null(current_tissue) || current_tissue == "") {
          missing_inputs <- c(missing_inputs, "Tissue Type")
        }
        if (is.null(current_resolution) || is.na(current_resolution)) {
          missing_inputs <- c(missing_inputs, "Hi-C Map Resolution")
        }
        if (is.null(current_mcool_path) || current_mcool_path == "") {
          missing_inputs <- c(missing_inputs, "Path to .mcool file")
        }
        
        if (length(missing_inputs) > 0) {
          status_msg(paste0("❌ Error: Please enter/select values for the following required inputs: ",
                            paste(missing_inputs, collapse = ", "),
                            ". Then click 'Start TADcalling' again."))
          return(NULL)
        }
        tissue_val <- current_tissue
        resolution_val <- current_resolution
        mcool_path_val <- current_mcool_path
        java_memory_val <- current_java_memory
        
        
        #--------------------------------------
        # Step 2: Normalize mcool path and copy it to project directory
        base_project_dir_val <- if (!is.null(project_output_dir)) project_output_dir() else NULL
        
        if (!is.null(base_project_dir_val) && base_project_dir_val != "") {
          base_project_dir_val <- normalizePath(base_project_dir_val, winslash = "/", mustWork = FALSE)
        }
        
        mcool_path_for_python <- gsub("\\\\", "/", mcool_path_val)
        if (!grepl("\\.mcool$", mcool_path_for_python, ignore.case = TRUE)) {
          mcool_path_for_python <- paste0(mcool_path_for_python, ".mcool")
        }
        
        if (!is.null(base_project_dir_val) && base_project_dir_val != "") {
          tadcaller_base_output_dir <- file.path(base_project_dir_val, "TADcaller_results")
          mcool_filename <- basename(mcool_path_for_python)
          mcool_filename_no_ext <- tools::file_path_sans_ext(mcool_filename)
          specific_mcool_output_dir <- file.path(tadcaller_base_output_dir, mcool_filename_no_ext)
          
          if (!dir.exists(specific_mcool_output_dir)) {
            dir.create(specific_mcool_output_dir, recursive = TRUE)
            status_msg(paste0(status_msg(),"<br>📁 Created directory for .mcool file specific results: ", specific_mcool_output_dir))
          }
          
          new_mcool_path_for_copy <- file.path(specific_mcool_output_dir, mcool_filename)
          
          copy_success <- FALSE
          if (file.exists(mcool_path_for_python)) { # Check original path existence
            copy_success <- file.copy(from = mcool_path_for_python, to = new_mcool_path_for_copy, overwrite = TRUE)
          } else {
            status_msg(paste0(status_msg(),"<br>❌ Error: Original .mcool file not found at: ", mcool_path_for_python))
            return(NULL)
          }
          
          if (copy_success) {
            status_msg(paste0(status_msg(),"<br>🔄 Copied .mcool file to project output directory: ", new_mcool_path_for_copy))
            mcool_path_for_python <- gsub("\\\\", "/", new_mcool_path_for_copy)
          } else {
            status_msg(paste0(status_msg(),"<br>❌ Error: Failed to copy .mcool file from '", mcool_path_for_python, "' to '", new_mcool_path_for_copy, "'. Check permissions and paths."))
            return(NULL)
          }
        } else {
          # If project_output_dir is null, use the original mcool_path and normalize slashes
          status_msg(paste0(status_msg(),"<br>ℹ️ Using original .mcool file path: ", mcool_path_for_python))
        }
        
        # Print the final mcool_path that will be passed to Python
        print(paste("DEBUG: Final mcool_path passed to Python:", mcool_path_for_python))
        print(paste("DEBUG: Does final mcool_path exist?", file.exists(mcool_path_for_python)))
        
        
        #---------------------------------------------
        # Step 3. extract contact matrixes of all chromosomes
        incProgress(0.1, detail = "Step 1: Extracting contact matrix for all chromosomes from .mcool file ")
        status_msg("Step 1: Processing Hi-C map to extract contact matrices for ALL chromosomes at chosen resolution...")
        
        # Step A: Process .mcool file for ALL chromosomes using Python
        python_overall_success <- FALSE
        all_chrom_results <- NULL
        try({
          all_chrom_results <- process_contact_matrix_all_chroms( # Call the Python function
            tissue = tissue_val,
            resolution_kb = resolution_val,
            path_mcool_file = mcool_path_for_python
          )
          python_overall_success <- all_chrom_results$success
        }, silent = TRUE)
        
        if (is.null(all_chrom_results) || !python_overall_success) {
          error_message <- if (!is.null(all_chrom_results$overall_message)) {
            paste0("❌ Error in Step 1 (Overall Python process): ", all_chrom_results$overall_message)
          } else {
            "❌ Unknown error occurred during Python Step 1 (extracting contact matrices of all chromosomes). Check console for Python errors."
          }
          status_msg(paste0(status_msg(),"<br>", error_message))
          return(NULL)
        }
        
        status_msg(paste0(status_msg(),"<br>✅ Step 1 completed. Found ",
                          length(all_chrom_results$chromosomes_processed),
                          " chromosomes with valid data at ", resolution_val, "kb resolution."))
        if ('total_contacts_overall_mcool' %in% names(all_chrom_results)) {
          status_msg(paste0(status_msg(), "<br>Total contacts in .mcool file: ", format(all_chrom_results$total_contacts_overall_mcool, big.mark = ",")))
        }
        
        # Step b: Define results folders (consistent for all chroms from this mcool file)
        # These paths are where the Python script *already* saved the contact matrices
        # and where the Java and R post-processing will save their outputs.
        base_dir_of_mcool <- dirname(mcool_path_for_python)
        if (is.null(base_dir_of_mcool) || base_dir_of_mcool == "") {
          base_dir_of_mcool <- getwd()
        }
        results_root_folder <- file.path(base_dir_of_mcool, "TADcaller_Results", paste0("TADs_", tissue_val))
        java_raw_output_dir_base <- file.path(results_root_folder, "java_raw_output")
        processed_tads_output_dir <- file.path(results_root_folder, "processed_tads")
        
        # Update reactive values for paths and current parameters
        results_root_folder_rv(results_root_folder)
        java_raw_output_dir_rv(java_raw_output_dir_base) # Base directory for java raw output
        processed_tads_output_dir_rv(processed_tads_output_dir)
        current_tissue_rv(tissue_val)
        current_resolution_rv(resolution_val)
        current_java_memory_rv(java_memory_val)
        
        # Ensure directories exist for Java and R outputs
        if (!dir.exists(java_raw_output_dir_base)) dir.create(java_raw_output_dir_base, recursive = TRUE)
        if (!dir.exists(processed_tads_output_dir)) dir.create(processed_tads_output_dir, recursive = TRUE)
        
        #-------------------------------------------------
        # Step4:  Iterate through each chromosome and run TADcaller + post-processing ---
        processed_chroms <- c()
        num_successful_chroms <- 0
        total_chroms_from_python <- length(all_chrom_results$chromosomes_processed)
        status_msg(paste0(status_msg(),"<br>  Step 2+3: (2)Running TADcaller deDoc2 for all chromosomes provided from Step1 while directly (3)extracting their TAD and SubTAD tables..."))
        
        for (i in seq_along(all_chrom_results$chromosomes_processed)) {
          chrom_info <- all_chrom_results$chromosomes_processed[[i]]
          current_chrom_normalized <- chrom_info$chrom # e.g., "chr1", "chrX"
          input_matrix_path <- chrom_info$output_file # Path to the .tsv generated by Python
          
          if (!chrom_info$success || is.null(input_matrix_path) || !file.exists(input_matrix_path)) {
            status_msg(paste0(status_msg(),"<br>⚠️ Skipping ", current_chrom_normalized, ": Python step failed or input matrix not found for this chromosome. Message: ", chrom_info$message, " Error: ", chrom_info$error))
            next # Skip to the next chromosome
          }
          incProgress(0.1 + (0.6 * (i / total_chroms_from_python)),
                      detail = paste0("Step 2: Running TADcaller for ", current_chrom_normalized, " (", i, "/", total_chroms_from_python, ")"))
          
          jar_path <- "../PythonProject/deDoc2-main/deDoc2.jar"
          if (!file.exists(jar_path)) {
            status_msg(paste0(status_msg(),"<br>❌ .jar file not found at: ", normalizePath(jar_path, mustWork = FALSE), ". Please ensure the JAR file is in the correct location and accessible."))
            return(NULL) # Stop entire process if JAR is missing
          }
          
          # Specific output paths for each chromosome's raw java output
          java_raw_output_file_base <- file.path(java_raw_output_dir_base, paste0(tissue_val, "_", current_chrom_normalized, "_", resolution_val, "kb_TADlines"))
          
          java_args <- c(
            paste0("-Xms", java_memory_val),
            "-jar", jar_path,
            "-inputfile", input_matrix_path,
            "-binsize", resolution_val,
            "-outputfile", java_raw_output_file_base
          )
          
          result_java <- tryCatch({
            processx::run("java", args = java_args, echo = TRUE, error_on_status = FALSE)
          }, error = function(e) {
            list(status = 1, stdout = "", stderr = as.character(e))
          })
          
          if (result_java$status != 0) {
            err_msg <- result_java$stderr
            status_msg(paste0(status_msg(),"<br>❌ Java TADcaller error for ", current_chrom_normalized, ":\n", err_msg))
            next # Skip to next chromosome on Java error for this one
          }
          

          incProgress(0.1 + (0.8 * (i / total_chroms_from_python)),
                      detail = paste0("Step 3: Processing TAD results for ", current_chrom_normalized))
          tryCatch({
            processed_tads_df <- process_tad_results(
              chrom = current_chrom_normalized, # Pass normalized chrom for consistent internal naming
              tissue = tissue_val,
              binsize = resolution_val,
              input_dir = java_raw_output_dir_base, # Pass the base directory for raw java outputs
              output_dir = processed_tads_output_dir # Pass the processed output directory
            )
            
            if (is.null(processed_tads_df) || !is.data.frame(processed_tads_df) || nrow(processed_tads_df) == 0) {
              status_msg(paste0(status_msg(),"<br>❌ Error: No TADs were processed for ", current_chrom_normalized, ". Check the Java raw output files."))
              next # Skip to next chromosome
            }
            
            processed_subtads_df <- process_subtad_results(
              chrom = current_chrom_normalized, # Pass normalized chrom
              tissue = tissue_val,
              binsize = resolution_val,
              input_dir = java_raw_output_dir_base, # Pass the base directory for raw java outputs
              output_dir = processed_tads_output_dir, # Pass the processed output directory
              main_tads_df_input = processed_tads_df
            )
            
            processed_chroms <- c(processed_chroms, current_chrom_normalized) # Add to list of successfully processed chroms
            num_successful_chroms <- num_successful_chroms + 1
          
          }, error = function(e) {
            status_msg(paste0(status_msg(),"<br>❌ Error in Step 3 (Post-processing) for ", current_chrom_normalized, ": ", as.character(e)))
            next # Skip to next chromosome on error
          })
        } # End for loop for chromosomes
        
        print(paste0("processed_chroms: ", paste(processed_chroms, collapse = ", ")))
        processed_chroms_list_rv(processed_chroms) # Update reactive value with list of processed chroms
        print(paste0("chr list: " ,processed_chroms_list_rv(),collapse = ", "))
        #----------------------------------------------------------------
        
        status_msg(paste0(status_msg(), "<br>🎉 TAD calling analysis complete for the following chromosomes:<br>", paste(processed_chroms_list_rv(), collapse = ", ")))
        status_msg(paste0(status_msg(), "<br>📁 All results are saved under: '", normalizePath(results_root_folder, winslash = "/"), "'"))
        incProgress(1, detail = "Analysis Complete!")
      }) 
    }) 
    
    #----------------------------------
    output$chromosome_selector_ui <- renderUI({
      # Debugging print statement (already present, good!)
      print(paste("R (Debug): Rendering chromosome_selector_ui. processed_chroms_list_rv() has", length(processed_chroms_list_rv()), "elements."))
      if (length(processed_chroms_list_rv()) == 0) {
        return(
          selectInput(
            inputId = session$ns("selected_chrom"), # <--- Corrected: use session$ns()
            label = "Select Chromosome:",
            choices = character(0), # Empty choices
            selected = NULL,
            width = "100%",
            multiple = FALSE
          )
        )
      } else {
        selectInput(
          inputId = session$ns("selected_chrom"), # <--- Corrected: use session$ns()
          label = "Select Chromosome:",
          choices = processed_chroms_list_rv(),
          selected = processed_chroms_list_rv()[1],
          width = "100%",
          multiple = FALSE
        )
      }
    })
    #---------------------------------------------------------------------------------------------------
    # --- Reactive for displaying the selected table (now defaults to first processed chrom) ---
    observeEvent({
      input$table_choice
      input$selected_chrom # Correctly included as a dependency
      processed_tads_output_dir_rv()    # Trigger if analysis just finished
    }, {

      # --- PREVENT RACE CONDITIONS ---
      req(processed_tads_output_dir_rv())       # Ensure analysis ran
      req(current_tissue_rv(), current_resolution_rv())  # Ensure parameters set
      req(processed_chroms_list_rv())           # Ensure list is ready
      print(paste("R (Debug): Checking reqs - input$selected_chrom is", !is.null(input$selected_chrom)))
      req(input$selected_chrom)                 # Ensure user has selected chrom
      
      if (is.null(input$selected_chrom) || input$selected_chrom == "") {
        print("R (Debug): input$selected_chrom is NULL or empty, returning.")
        output$tad_table <- DT::renderDataTable(NULL) # Clear table if no chrom is selected
        return(NULL)
      }
      print(paste("R (Debug): input$selected_chrom is OK:", input$selected_chrom))
      
      # --- MAIN LOGIC ---
      chrom_to_display <- input$selected_chrom
      tissue_val <- current_tissue_rv()
      resolution_val <- current_resolution_rv()
      table_type_val <- selected_table() # This will be "TAD" or "SubTAD"
      
      print(paste("R (Debug): Variables for table path - chrom:", chrom_to_display,
                  "tissue:", tissue_val, "resolution:", resolution_val, "table_type:", table_type_val))
      
      table_file_name <- if (table_type_val == "TAD") {
        sprintf("%s_%s_%skb_TADs.txt", tissue_val, chrom_to_display, resolution_val)
      } else {
        sprintf("%s_%s_%skb_SubTADs.txt", tissue_val, chrom_to_display, resolution_val)
      }
      
      table_path <- file.path(processed_tads_output_dir_rv(), table_file_name)
      
      print(paste("R (Debug): Attempting to read table for display from:", table_path))
      
      if (!file.exists(table_path)) {
        current_selection_label <- if (table_type_val == "TAD") "TADs" else "SubTADs"
        status_msg(paste0(status_msg(), "<br>⚠️ Warning: The requested '", current_selection_label,
                          "' table file for ", chrom_to_display,
                          " not found: '", table_path, "'"))
        output$tad_table <- DT::renderDataTable(NULL)
        print(paste("R (Debug): Table file not found:", table_path))
        return(NULL)
      }
      
      tryCatch({
        data <- readr::read_delim(table_path, delim = "\t", show_col_types = FALSE)
        print(paste("R (Debug): Successfully read data with dimensions:", paste(dim(data), collapse = "x")))
        if (nrow(data) == 0) {
          status_msg(paste0(status_msg(), "<br>⚠️ Warning: Table for ", chrom_to_display, " (", table_type_val, ") is empty. No data to display."))
          output$tad_table <- DT::renderDataTable(NULL)
          print("R (Debug): Data table is empty.")
        } else {
          output$tad_table <- DT::renderDataTable(data, options = list(pageLength = 10))
          print("R (Debug): renderDataTable called with data.")
        }
      }, error = function(e) {
        status_msg(paste0(status_msg(), "<br>❌ Error displaying table for ", chrom_to_display,
                          ": Could not read or process the data from '", basename(table_path),
                          "'. Error: ", as.character(e)))
        output$tad_table <- DT::renderDataTable(NULL)
        print(paste("R (Debug): Error reading/rendering table:", as.character(e)))
      })
    }, ignoreInit = TRUE)
    
    
    
    # Return a list of reactive values that the parent app can access
    return(list(
      current_tissue = current_tissue_rv,
      current_resolution = current_resolution_rv,
      processed_chroms_list_rv = processed_chroms_list_rv,
      processed_output_dir = processed_tads_output_dir_rv # Path to folder with all TAD/SubTAD files
    ))
  }) # moduleServer close
} # function close


# test code of chromosome is in the r.scripts in file 09_test