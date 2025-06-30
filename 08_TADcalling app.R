library(shiny)
library(bslib)
library(reticulate)
library(processx)
library(readr)
library(openxlsx)
library(DT) 

source("utils/TADcalling_utils.R")

ui_TADcalling <- page_navbar(
  title = "EPIC Array Pipeline",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  
  nav_panel(
    "TADcalling",
    layout_sidebar(
      sidebar = sidebar(
        textInput(
          inputId = "tissue",
          label = "Tissue Type",
          placeholder = "CAKI2"
        ),
        helpText("Tissue of your mcool file/ The source tissue of the hic-map."),
        
        textInput(
          inputId = "chromosome",
          label = "Chromosome",
          placeholder = "ChrX, chrx or x"
        ),
        
        numericInput(
          inputId = "resolution",
          label = "Hi-C Map Resolution (in kb)",
          value = 25,
          min = 1,
          step = 1
        ),
        
        textInput(
          inputId = "mcool_path",
          label = "Path to .mcool file",
          placeholder = "C:\\path\\to\\your\\file.mcool",
          width = "100%"
        ),
        
        selectInput(
          inputId = "java_memory",
          label = "Java Memory Allocation",
          choices = c("4g", "8g", "16g"),
          selected = "8g"
        ),
        helpText("💡 Don't choose more than the limit of your RAM. Recommended: use the maximum available RAM for better performance, as this process is demanding."),
        
        actionButton("start_tadcalling", "🚀 Start TADcalling", class = "btn-primary"),
        
        hr(),
        
        radioButtons(
          inputId = "table_choice",
          label = "Table to Display",
          choices = c("TADs table", "SubTADs table"),
          selected = "TADs table"
        )
      ),
      
      div(class = "col-lg-9 col-md-9 col-sm-12", # Use col-lg-9 for large, col-md-9 for medium, col-sm-12 for small screens
          uiOutput("status_card"),
          br(),
          DT::dataTableOutput("tad_table")
      )
      
    )
  )
)

server_TADcalling <- function(input, output, session) {
  # Set up environment
  use_virtualenv("shiny_py", required = TRUE)
  source_python("../PythonProject/1-Hi-c_maps_exploring.py") # Ensure this path is correct for your project
  
  # Reactive value to store status messages
  status_msg <- reactiveVal(
    "📥 Enter the required inputs on the left and click 'Start TADcalling' to begin.<br>The resolution refers to the bin size used for the contact map, e.g. 25 kb means each matrix bin spans 25,000 base pairs."
  )
  
  output$status_card <- renderUI({
    msg <- status_msg()
    
    # Determine background color class based on message content
    card_color <- if (grepl("❌", msg) || grepl("Error", msg)) {
      "bg-danger text-white"
    } else {
      "bg-light"
    }
    
    # Render the card UI
    div(
      class = paste("card shadow-sm mb-4 p-3", card_color),
      h4("🧬 Analysis Status"),
      HTML(gsub("\n", "<br>", msg))  # Converts \n into proper line breaks
    )
  })
  
  
  # Reactive to track which table to show
  selected_table <- reactiveVal("TAD")
  
  observeEvent(input$table_choice, {
    selected_table(ifelse(input$table_choice == "SubTADs table", "SubTAD", "TAD"))
  })
  
  
  # --- Main Analysis Logic ---
  observeEvent(input$start_tadcalling, {
    # --- ALWAYS RESET STATUS MESSAGE AT THE VERY BEGINNING OF THE ANALYSIS ---
    status_msg("Step 1: Processing .mcool file...")
    
    # Input Validation: Check for missing inputs and provide specific messages
    missing_inputs <- c()
    if (is.null(input$tissue) || input$tissue == "") {
      missing_inputs <- c(missing_inputs, "Tissue Type")
    }
    if (is.null(input$chromosome) || input$chromosome == "") {
      missing_inputs <- c(missing_inputs, "Chromosome")
    }
    # numericInput for resolution always has a value, but you can check if it's reasonable
    # if (is.null(input$resolution) || !is.numeric(input$resolution)) {
    #   missing_inputs <- c(missing_inputs, "Hi-C Map Resolution")
    # }
    if (is.null(input$mcool_path) || input$mcool_path == "") {
      missing_inputs <- c(missing_inputs, "Path to .mcool file")
    }
    # selectInput for java_memory always has a value
    
    if (length(missing_inputs) > 0) {
      status_msg(paste0("❌ Error: Please enter values for the following required inputs: ",
                        paste(missing_inputs, collapse = ", "),
                        ". Then click 'Start TADcalling' again."))
      return(NULL) # Stop execution if inputs are missing
    }
    
    # Fetch inputs directly from input$ (since req() ensures they exist or we've returned early)
    tissue_val <- input$tissue
    chromosome_val <- input$chromosome
    resolution_val <- input$resolution
    mcool_path_val <- input$mcool_path
    java_memory_val <- input$java_memory
    
    
    #1.. Normalize chromosome input (e.g., "2" becomes "chr2")
    chrom <- chromosome_val
    if (!grepl("^chr", chrom, ignore.case = TRUE)) {
      chrom <- paste0("chr", chrom)
    }
    # Ensure the entire string is lowercase (e.g., "Chr3" becomes "chr3")
    chrom <- tolower(chrom)
    
    #2.. Normalize mcool path to use forward slashes (important for Python consistency)
    mcool_path <- gsub("\\\\", "/", mcool_path_val)
    print(paste("R (Debug): Path being passed to Python:", mcool_path))
    
    
    # Step 1: Process .mcool file using Python
    python_step_success <- FALSE
    result <- NULL # Initialize result to NULL
    try({
      result <- process_contact_matrix(
        tissue = tissue_val,
        chrom = chrom,
        resolution_kb = resolution_val,
        path_mcool_file = mcool_path # Use the normalized mcool_path here
      )
      python_step_success <- result$success
    }, silent = TRUE)
    
    # Debugging print
    print("Python process_contact_matrix result:")
    print(result)
    
    # Check for Python script's success
    if (is.null(result) || !python_step_success) { # Simplified check
      error_message <- if (!is.null(result) && !is.null(result$error) && !is.null(result$message)) {
        paste0("❌ Error in Step 1: ", result$error, "<br>", result$message)
      } else {
        "❌ Unknown error occurred during Python Step 1. Check console for Python errors."
      }
      status_msg(paste0(status_msg(),"<br>", error_message))
      return(NULL) # Exit the observeEvent if Step 1 failed
    }
    
    # Append success message for Step 1
    status_msg(paste0(status_msg(),"<br>✅ Step 1 completed.<br>",
                      "Total contacts: in whole map ", format(result$total_contacts, big.mark = ",")))
    
    # --- Define the base root for all results for this run (NEW HIERARCHICAL STRUCTURE) ---
    # This will be like: /path/to/mcool_dir/TADcaller_Results/TADs_YourTissue/
    base_dir_of_mcool <- dirname(mcool_path)
    if (is.null(base_dir_of_mcool) || base_dir_of_mcool == "") {
      base_dir_of_mcool <- getwd()
    }
    # This is the overall folder for all results related to this specific run (tissue, chrom, resolution)
    results_root_folder <- file.path(base_dir_of_mcool, "TADcaller_Results", paste0("TADs_", tissue_val))
    
    # Create the root folder if it doesn't exist
    if (!dir.exists(results_root_folder)) dir.create(results_root_folder, recursive = TRUE)
    
    
    # Define specific subfolders for Java raw output and R processed output
    java_raw_output_dir <- file.path(results_root_folder, "java_raw_output")
    processed_tads_output_dir <- file.path(results_root_folder, "processed_tads")
    
    # Ensure these specific subfolders exist
    if (!dir.exists(java_raw_output_dir)) dir.create(java_raw_output_dir, recursive = TRUE)
    if (!dir.exists(processed_tads_output_dir)) dir.create(processed_tads_output_dir, recursive = TRUE)
    
    
    # Step 2: Run TADcaller (Java deDoc2)
    status_msg(paste0(status_msg(),"<br>  Step 2: Running TADcaller deDoc2..."))
    
    # input_matrix is the path to the contact matrix generated by Python, which is now in its dedicated subfolder
    input_matrix <- result$output_file
    
    print(paste("R (Debug): Attempting to use input_matrix for Java:", input_matrix))
    if (is.null(input_matrix) || !file.exists(input_matrix)) {
      print(paste("R (Debug): input_matrix is NULL or file doesn't exist:", input_matrix))
      status_msg(paste0(status_msg(),"<br>❌ Input matrix file not found. Ensure Step 1 completed correctly and generated the file."))
      return(NULL) # Exit if matrix file is missing
    }
    
    binsize_kb <- resolution_val
    tissue <- tissue_val
    chrom_for_java <- chrom # Use the normalized chrom
    java_memory <- java_memory_val
    jar_path <- "../PythonProject/deDoc2-main/deDoc2.jar" # Path relative to your R app's working directory
    
    if (!file.exists(jar_path)) {
      status_msg(paste0(status_msg(),"<br>❌ .jar file not found at: ", normalizePath(jar_path, mustWork = FALSE), ". Please ensure the JAR file is in the correct location and accessible."))
      return(NULL) # Exit if JAR is missing
    }
    
    # Java will output its raw files (.TAD, .window.TAD) to the dedicated java_raw_output_dir
    output_base <- file.path(java_raw_output_dir, paste0(tissue, "_", chrom_for_java, "_", binsize_kb, "kb_TADlines"))
    
    java_args <- c(
      paste0("-Xms", java_memory),
      "-jar", jar_path,
      "-inputfile", input_matrix,
      "-binsize", binsize_kb,
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
      return(NULL) # Exit if Java execution failed
    }
    
    status_msg(paste0(status_msg(),"<br>✅ Step 2 completed"))
    print(paste("tADcalling finished see folder: ",output_base ))
    
    # Step 3: Post-process TAD results
    status_msg(paste0(status_msg(),"<br>🧬 Step 3: Processing TAD results..."))
    
    tryCatch({
      # Call process_tad_results and capture its return value (the TADs data frame)
      processed_tads_df <- process_tad_results(
        chrom = chrom_for_java,
        tissue = tissue,
        binsize = binsize_kb,
        input_dir = java_raw_output_dir,
        output_dir = processed_tads_output_dir
      )
      print(paste("process_tad_results finished see folder: ", processed_tads_output_dir))
      
      # Check if processed_tads_df is valid before proceeding
      if (is.null(processed_tads_df) || !is.data.frame(processed_tads_df) || nrow(processed_tads_df) == 0) {
        status_msg(paste0(status_msg(),"<br>❌ Error: No TADs were processed or `process_tad_results` returned an empty data frame. Check the Java raw output files."))
        return(NULL)
      }
      
      # Pass the processed_tads_df directly to process_subtad_results
      process_subtad_results(
        chrom = chrom_for_java,
        tissue = tissue,
        binsize = binsize_kb,
        input_dir = java_raw_output_dir,
        output_dir = processed_tads_output_dir,
        main_tads_df_input = processed_tads_df 
      )
      print(paste("process_Subtad_results finished see folder: ", processed_tads_output_dir))
      
      status_msg(paste0(status_msg(),"<br>✅ Step 3 completed"))
    }, error = function(e) {
      status_msg(paste0(status_msg(),"<br>❌ Error in Step 3 (Post-processing): ", as.character(e)))
      return(NULL)
    })
    
    
    status_msg(paste0(status_msg(), "<br>📁 All results are saved under: '", normalizePath(results_root_folder, winslash = "/"), "'"))
    
    
    
    # Render result table from the processed_tads_output_dir
    output$tad_table <- DT::renderDataTable({ 
      table_path <- if (selected_table() == "TAD") {
        file.path(processed_tads_output_dir, sprintf("%s_%s_%skb_TADs.txt", tissue, chrom_for_java, binsize_kb))
      } else {
        file.path(processed_tads_output_dir, sprintf("%s_%s_%skb_SubTADs.txt", tissue, chrom_for_java, binsize_kb))      }
      
      print(paste("R (Debug): Attempting to read table from:", table_path)) 
      
      if (!file.exists(table_path)) {
        # Update status message to indicate which table it tried to load
        current_selection <- if (selected_table() == "TAD") "TADs" else "SubTADs"
        status_msg(paste0(status_msg(), "<br>⚠️ Warning: The requested '", current_selection, "' table file not found: '", table_path, "'"))
        return(NULL)
      }
      
      tryCatch({
        data <- read.delim(table_path)
        print(paste("R (Debug): Successfully read data with dimensions:", paste(dim(data), collapse = "x"))) # Added debug print
        DT::datatable(data) # Wrapped data with DT::datatable
      }, error = function(e) {
        status_msg(paste0(status_msg(), "<br>❌ Error displaying table: Could not read or process the data from '", basename(table_path), "'. Error: ", as.character(e)))
        return(NULL)
      })
    })
  })
  
  # Ensure the table updates when the radio button changes AFTER the initial run
  observeEvent(input$table_choice, {
    # This observeEvent will trigger output$tad_table to re-render if the relevant
    # variables (processed_tads_output_dir, tissue, chrom_for_java, binsize_kb)
    # have been set by the prior execution of observeEvent(input$start_tadcalling, ...)
    # If a run hasn't occurred yet, this will just try to render from potentially uninitialized paths,
    # leading to the "file not found" warning, which is acceptable initial behavior.
    if (exists("processed_tads_output_dir")) { # Check if path variable exists from a previous run
      # Manually trigger the renderDataTable if it's already set up
      output$tad_table <- DT::renderDataTable({ # Recalculate the table for the new selection
        table_path <- if (selected_table() == "TAD") {
          file.path(processed_tads_output_dir, sprintf("%s_%s_%skb_TADs.txt", tissue_val, chrom_for_java, resolution_val)) # Use global vars
        } else {
          file.path(processed_tads_output_dir, sprintf("%s_%s_%skb_SubTADs.txt", tissue_val, chrom_for_java, resolution_val))
        }
        
        print(paste("R (Debug): (table_choice) Attempting to read table from:", table_path))
        
        if (!file.exists(table_path)) {
          current_selection <- if (selected_table() == "TAD") "TADs" else "SubTADs"
          status_msg(paste0(status_msg(), "<br>⚠️ Warning: The requested '", current_selection, "' table file not found: '", table_path, "'. Perhaps the analysis hasn't run yet or failed."))
          return(NULL)
        }
        tryCatch({
          data <- read.delim(table_path)
          print(paste("R (Debug): (table_choice) Successfully read data with dimensions:", paste(dim(data), collapse = "x")))
          DT::datatable(data)
        }, error = function(e) {
          status_msg(paste0(status_msg(), "<br>❌ Error displaying table from table_choice: ", as.character(e)))
          return(NULL)
        })
      })
    } else {
      status_msg(paste0(status_msg(), "<br>Please run the TADcalling analysis first to generate tables."))
    }
  }, ignoreInit = TRUE) # ignoreInit prevents it from running when the app first starts
}

shinyApp(ui = ui_TADcalling, server = server_TADcalling)