# utils/tadcalling_processing_utils.R 
# Date: 14.07.2025

# ─────────────────────────────────────
# SERVER FUNCTION (PURE PROCESSING VERSION)
# ─────────────────────────────────────
#' @param id Module ID.
#' @param project_output_dir A reactive expression for the project's main output directory path.
#' @param trigger_button A reactive expression (e.g., input$button) to trigger the analysis.
#' @param tissue_r A reactive expression for the tissue type.
#' @param chromosome_r A reactive expression for the chromosome.
#' @param resolution_r A reactive expression for the resolution.
#' @param mcool_path_r A reactive expression for the path to the .mcool file.
#' @param java_memory_r A reactive expression for Java memory settings.
#' @return A list of reactive expressions containing the processed TADs and SubTADs data frames,
#'         the output directory for processed TADs, and the Java memory used.
tadcalling_processing_server <- function(id,
                                         project_output_dir = reactive(NULL),
                                         trigger_button = reactive(NULL),
                                         tissue_r = reactive(NULL),
                                         chromosome_r = reactive(NULL),
                                         resolution_r = reactive(NULL),
                                         mcool_path_r = reactive(NULL),
                                         java_memory_r = reactive(NULL)) {
  
  moduleServer(id, function(input, output, session) {
    
    # Set up Python environment (assuming this path is globally valid or relative to app root)
    use_virtualenv("shiny_py", required = TRUE)
    source_python("../python_scripts/1-Hi-c_maps_exploring.py")
    
    # Reactive values to store status messages (for internal debugging or logging, not UI)
    status_msg <- reactiveVal("Processing module initialized.")
    
    # Reactive values to store the TAD and SubTAD data frames (these are the primary outputs)
    tads_df_rv <- reactiveVal(NULL)
    subtads_df_rv <- reactiveVal(NULL)
    
    
    # Reactive values to persist paths and parameters after a successful run (for return list)
    results_root_folder_rv <- reactiveVal(NULL)
    processed_tads_output_dir_rv <- reactiveVal(NULL) # Added to return list
    java_raw_output_dir_rv <- reactiveVal(NULL)
    current_tissue_rv <- reactiveVal(NULL)
    current_chrom_rv <- reactiveVal(NULL)
    current_resolution_rv <- reactiveVal(NULL)
    current_java_memory_rv <- reactiveVal(NULL) # Added to return list
    
    
    # --- Main Analysis Logic ---
    # Observe the 'trigger_button' reactive argument instead of an internal input$start_tadcalling
    observeEvent(trigger_button(), {
      # Use req() to ensure all necessary reactive inputs are available before proceeding
      req(tissue_r(), chromosome_r(), resolution_r(), mcool_path_r(), java_memory_r())
      
      # Capture current input values from the *passed reactive arguments*
      current_tissue <- tissue_r()
      current_chromosome <- chromosome_r()
      current_resolution <- resolution_r()
      current_mcool_path <- mcool_path_r()
      current_java_memory <- java_memory_r()
      
      # Wrap the entire analysis process in withProgress
      withProgress(message = 'TAD Calling in Progress', value = 0, {
        
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
          return(NULL) # Stop execution if validation fails
        }
        
        tissue_val <- current_tissue
        chromosome_val <- current_chromosome
        resolution_val <- current_resolution
        mcool_path_val <- current_mcool_path
        java_memory_val <- current_java_memory
        print(mcool_path_val) # Debug print
        
        # 1. Normalize chromosome input
        chrom_temp <- gsub("^chr", "", chromosome_val, ignore.case = TRUE)
        if (toupper(chrom_temp) %in% c("X", "Y", "M")) {
          chrom <- paste0("chr", toupper(chrom_temp))
        } else {
          chrom <- paste0("chr", tolower(chrom_temp))
        }
        
        
        # 2. Normalize mcool path
        base_project_dir_val <- if (!is.null(project_output_dir())) project_output_dir() else NULL # Access reactive with ()
        
        if (!is.null(base_project_dir_val) && base_project_dir_val != "") {
          base_project_dir_val <- normalizePath(base_project_dir_val, winslash = "/", mustWork = FALSE)
        }
        
        if (!grepl("\\.mcool$", mcool_path_val, ignore.case = TRUE)) {
          mcool_path_val <- paste0(mcool_path_val, ".mcool")
          status_msg(paste0(status_msg(),"<br>ℹ️ Appended '.mcool' extension to input path: "))
        }
        
        mcool_path <- ""
        if (!is.null(base_project_dir_val) && base_project_dir_val != "") {
          tadcaller_base_output_dir <- file.path(base_project_dir_val, "TADcaller_results")
          mcool_filename <- basename(mcool_path_val)
          mcool_filename_no_ext <- tools::file_path_sans_ext(mcool_filename)
          specific_mcool_output_dir <- file.path(tadcaller_base_output_dir, mcool_filename_no_ext)
          
          if (!dir.exists(specific_mcool_output_dir)) {
            dir.create(specific_mcool_output_dir, recursive = TRUE)
            status_msg(paste0(status_msg(),"<br>📁 Created directory for .mcool file: ", specific_mcool_output_dir))
          }
          
          new_mcool_path_for_copy <- file.path(specific_mcool_output_dir, mcool_filename)
          
          copy_success <- FALSE
          if (file.exists(mcool_path_val)) {
            copy_success <- file.copy(from = mcool_path_val, to = new_mcool_path_for_copy, overwrite = TRUE)
          } else {
            status_msg(paste0(status_msg(),"<br>❌ Error: Original .mcool file not found at: ", mcool_path_val))
            return(NULL)
          }
          
          if (copy_success) {
            status_msg(paste0(status_msg(),"<br>🔄 Copied .mcool file to project output directory: ", new_mcool_path_for_copy))
            mcool_path <- gsub("\\\\", "/", new_mcool_path_for_copy)
          } else {
            status_msg(paste0(status_msg(),"<br>❌ Error: Failed to copy .mcool file from '", mcool_path_val, "' to '", new_mcool_path_for_copy, "'. Check permissions and paths."))
            return(NULL)
          }
        } else {
          mcool_path <- gsub("\\\\", "/", mcool_path_val)
          status_msg(paste0(status_msg(),"<br>ℹ️ Using original .mcool file path: ", mcool_path))
        }
        
        print(paste("DEBUG: Final mcool_path passed to Python:", mcool_path))
        print(paste("DEBUG: Does final mcool_path exist?", file.exists(mcool_path)))
        
        incProgress(0.1, detail = "Starting Step 1: Processing .mcool file")
        status_msg(paste0(status_msg(), "<br>Step 1: Processing Hi-c map (.mcool file) to extract contact matrix..."))
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
        
        base_dir_of_mcool <- dirname(mcool_path)
        if (is.null(base_dir_of_mcool) || base_dir_of_mcool == "") {
          base_dir_of_mcool <- getwd()
        }
        results_root_folder <- file.path(base_dir_of_mcool, "TADcaller_Results", paste0("TADs_", tissue_val))
        java_raw_output_dir <- file.path(results_root_folder, "java_raw_output")
        processed_tads_output_dir <- file.path(results_root_folder, "processed_tads")
        
        results_root_folder_rv(results_root_folder)
        java_raw_output_dir_rv(java_raw_output_dir)
        processed_tads_output_dir_rv(processed_tads_output_dir) # Update this reactive value
        current_tissue_rv(tissue_val)
        current_chrom_rv(chrom)
        current_resolution_rv(resolution_val)
        current_java_memory_rv(java_memory_val) # Update this reactive value
        
        if (!dir.exists(results_root_folder)) dir.create(results_root_folder, recursive = TRUE)
        if (!dir.exists(java_raw_output_dir)) dir.create(java_raw_output_dir, recursive = TRUE)
        if (!dir.exists(processed_tads_output_dir)) dir.create(processed_tads_output_dir, recursive = TRUE)
        
        incProgress(0.4, detail = "Step 2: Running TADcaller deDoc2")
        status_msg(paste0(status_msg(),"<br>  Step 2: Running TADcaller deDoc2..."))
        
        input_matrix <- result$output_file
        if (is.null(input_matrix) || !file.exists(input_matrix)) {
          status_msg(paste0(status_msg(),"<br>❌ Input matrix file not found. Ensure Step 1 completed correctly and generated the file."))
          return(NULL)
        }
        
        jar_path <- "../PythonProject/deDoc2-main/deDoc2.jar"
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
        status_msg(paste0(status_msg(),"<br>Step 3: Processing TAD results..."))
        
        tryCatch({
          processed_tads_df <- process_tad_results(
            chrom = chrom,
            tissue = tissue_val,
            binsize = resolution_val,
            input_dir = java_raw_output_dir,
            output_dir = processed_tads_output_dir
          )
          
          if (is.null(processed_tads_df) || !is.data.frame(processed_tads_df) || nrow(processed_tads_df) == 0) {
            status_msg(paste0(status_msg(),"<br>❌ Error: No TADs were processed or `process_tad_results` returned an empty data frame. Check the Java raw output files."))
            tads_df_rv(NULL)
            subtads_df_rv(NULL)
            return(NULL)
          }
          
          tads_df_rv(processed_tads_df)
          
          processed_subtads_df <- process_subtad_results(
            chrom = chrom,
            tissue = tissue_val,
            binsize = resolution_val,
            input_dir = java_raw_output_dir,
            output_dir = processed_tads_output_dir,
            main_tads_df_input = processed_tads_df
          )
          
          subtads_df_rv(processed_subtads_df)
          
          status_msg(paste0(status_msg(),"<br>✅ Step 3 completed"))
        }, error = function(e) {
          status_msg(paste0(status_msg(),"<br>❌ Error in Step 3 (Post-processing): ", as.character(e)))
          tads_df_rv(NULL)
          subtads_df_rv(NULL)
          return(NULL)
        })
        
        status_msg(paste0(status_msg(), "<br>📁 All results are saved under: '", normalizePath(results_root_folder, winslash = "/"), "'"))
        
        incProgress(1, detail = "Analysis Complete!")
      }) # End withProgress
      
      # THIS MODULE DOES NOT RENDER TABLES. It only processes and returns data.
    }) # End observeEvent(trigger_button())
    
    # Return a list of reactive values that the parent app/module can access
    return(list(
      tads_table = tads_df_rv,
      subtads_table = subtads_df_rv,
      status_message = status_msg # Optionally return status for parent to monitor
    ))
  })
}