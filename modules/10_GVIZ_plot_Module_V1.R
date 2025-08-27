# 10_Gviz_plot_Module_v1.R
# Author: Ghazal Sinjar
# Date: 02.07.2025
# Description:
# This Shiny module provides UI and server logic for visualizing genomic regions using the Gviz package.
# It enables dynamic plotting of differentially methylated regions (DMRs) along with gene annotations and other genomic features.
# The module integrates with upstream DMR identification results and uses EnsDb or TxDb objects for gene tracks.
# The visualization includes gene models, ideograms, axis tracks, and optional custom annotation tracks.

# ─────────────────────────────────────
# User Interface (UI) FUNCTION for the Module
# ─────────────────────────────────────
GvizPlotUI <- function(id) {
  ns <- NS(id)
  
  layout_sidebar(
    sidebar = sidebar(
      width = 300,
      
      selectInput(ns("region_type"), "Choose Region Type:",
                  choices = c("DMRs", "Off-targets", "Desired targeted region"),
                  selected = "DMRs"),
      
      uiOutput(ns("region_selector")),
      uiOutput(ns("chromosome_input")),
      numericInput(ns("from"), "From:", value = 1, min = 1, max = 1),
      numericInput(ns("to"), "To:", value = 1000, min = 1, max = 1),
      checkboxInput(ns("include_tads"), "Include TAD/SubTAD tracks", value = FALSE),
      br(),
      # Static display for Binsize
      div(
        class = "mb-1", 
        HTML('<label class="control-label">TAD/SubTAD Bin Size (kb):</label>'),
        uiOutput(ns("display_binsize"))
      ),
      
      # Static display for Tissue
      div(
        class = "mb-1", 
        HTML('<label class="control-label">TAD/SubTAD Tissue:</label>'),
        uiOutput(ns("display_tissue")) # This UI output will display the static value
      ),
      # Static display for qval
      div(
        class = "mb-1", 
        HTML('<label class="control-label">qval cutoff of DMPs:</label>'),
        uiOutput(ns("qval"))
      ),
      helpText("CpGs in green have q-values ≤ the value above."),
      
      actionButton(ns("create_plot"), "Create Plot", class = "btn btn-info"),
      helpText("Note: The boundaries of TADs and sub-TADs are approximate and depend on the resolution of the TAD-calling algorithm."),
      
      div(
        actionButton(ns("zoom_in"), "🔍 Zoom In", class = "btn btn-success btn-sm"),
        actionButton(ns("zoom_out"), "🔎 Zoom Out", class = "btn btn-warning btn-sm"),
        br(), br(),
        actionButton(ns("go_left"), "⬅️ Go Left", class = "btn btn-secondary btn-sm"),
        actionButton(ns("go_right"), "➡️ Go Right", class = "btn btn-secondary btn-sm")
      ),
      
      br(),
      downloadButton(ns("downloadPlot"), "Save Plot as PDF")
    ),
    plotOutput(ns("gvizPlot"), height = "1000px")
  )
}


# ─────────────────────────────────────
# SERVER FUNCTION for the Module
# ─────────────────────────────────────
GvizPlotServer <- function(id,
                           dmr_results,
                           dmp_results,
                           boxplot_results,
                           offtarget_table,
                           tadcalling_results,
                           chr_size_df_global,
                           tx_gr_filtered_global,
                           gr_cpgIslands_global,
                           genome = "hg38") {
  
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    print(paste0("Module Server (ID: ", id, ") started."))
    # --- Reactive Values for Plotting Range ---
    processed_chroms_list_rv <- reactiveVal(NULL)
    # These now specifically hold the range/chr for the *currently displayed plot*
    selectedRange <- reactiveVal(c(NA, NA))
    selectedChr <- reactiveVal(NA_character_)
    all_gviz_tracks_rv <- reactiveVal(NULL)
    # these hold the granges of tads and subtads for zooming and paning since they are on the same chromosome of the plot created
    gr_TADs <- reactiveVal(NULL)
    gr_SUBTADs <- reactiveVal(NULL)
    #--------------------------------------------------------------------
    # --- Static Binsize and Tissue Display ---
    output$display_binsize <- renderUI({
      # First, check if the reactive itself exists and is not NULL.
      # This prevents the "could not find function" error.
      if (!is.null(tadcalling_results)) {
        # Now that we know it's a valid object (likely a reactive), we can check its value.
        # req() is also safe to use here because it will stop execution if tadcalling_results() returns NULL.
        req(tadcalling_results())
        
        # The rest of your logic is sound.
        print(paste("DEBUG: display_binsize: class(tadcalling_results()) is", class(tadcalling_results())))
        print(paste("DEBUG: display_binsize: names(tadcalling_results()) are", paste(names(tadcalling_results()), collapse = ", ")))
        
        req(tadcalling_results()$current_resolution())
        print(paste("DEBUG: display_binsize: class(tadcalling_results()$current_resolution) value is", class(tadcalling_results()$current_resolution())))
        
        p(paste0(tadcalling_results()$current_resolution(), " kb"), class = "form-control-plaintext")
      } else {
        p("Not available, No TADs are called.", class = "form-control-plaintext")
      }
    })
    
    #--------------------------
    output$display_tissue <- renderUI({
      # First, check if the reactive expression itself is not NULL
      if (!is.null(tadcalling_results)) {
        # If it exists, call it and then use req() to ensure its value is not NULL
        req(tadcalling_results())
        
        # DEBUG: Check what tadcalling_results() returns
        print(paste("DEBUG: display_tissue: class(tadcalling_results()) is", class(tadcalling_results())))
        print(paste("DEBUG: display_tissue: names(tadcalling_results()) are", paste(names(tadcalling_results()), collapse = ", ")))
        
        req(tadcalling_results()$current_tissue()) # Get the reactive value from tadcalling_results
        # DEBUG: Check the class of the reactive object itself
        print(paste("DEBUG: display_tissue: class(tadcalling_results()$current_tissue()) is", class(tadcalling_results()$current_tissue())))
        
        p(tadcalling_results()$current_tissue(), class = "form-control-plaintext")
      } else {
        p("Not available, No TADs are called.", class = "form-control-plaintext")
      }
    })
    # --- End Static Binsize and Tissue Display ---
    #------------------------------------
    output$qval <- renderUI({
      
      req(dmp_results()$last_qvalue()) 
      # DEBUG: Check the class of the reactive object itself
      print(paste("DEBUG: qval: class(dmp_results()$last_qvalue()) is", class(dmp_results()$last_qvalue())))
      
      p(dmp_results()$last_qvalue(), class = "form-control-plaintext")
    })
    
    #-------------------------------------------
    # get chr list from tadcalling module
    observe({
      # First, check if the reactive exists and is not NULL before trying to call it.
      if (!is.null(tadcalling_results) && !is.null(tadcalling_results())) {
        # Now that we know it's a valid reactive, we can use req()
        req(tadcalling_results())
        processed_chroms_list_rv(tadcalling_results()$processed_chroms_list_rv())
      } else {
        # You might want to handle the case where tadcalling_results() is NULL
        # For example, by setting the list to an empty value or showing a notification.
        processed_chroms_list_rv(NULL)
        showNotification("TAD calling results are not available.", type = "warning", duration = 8)
      }
    })
    #--------------------------------------------------------
    # --- Reactive GRanges Objects from input tables ---
    gr_dmrs <- reactive({
      req(dmr_results())
      req(dmr_results()$dmr_table())
      print("Reactive: gr_dmrs is calculating...")
      tryCatch({
        gr <- create_gr_dmrs(dmr_results()$dmr_table())
        print(paste("Reactive: gr_dmrs created with", length(gr), "DMRs."))
        gr
      }, error = function(e) {
        message("Error in create_gr_dmrs(): ", e$message)
        NULL
      })
    })
    
    num_samples <- reactive({
      req(dmr_results())
      req(dmr_results()$pheno)
      req(dmr_results()$pheno()$Sample_Group)
      print("Reactive: num_samples is calculating...")
      tryCatch({
        n_samples <- length(dmr_results()$pheno()$Sample_Group)
        print(paste("Reactive: num_samples calculated as", n_samples))
        n_samples
      }, error = function(e) {
        message("Error calculating num_samples: ", e$message)
        0
      })
    })
    
    gr_cpgs <- reactive({
      req(boxplot_results())
      
      print("Reactive: gr_cpgs is calculating...")
      tryCatch({
        gr <- create_gr_cpgs(boxplot_results())
        print(paste("Reactive: gr_cpgs created with", length(gr), "CpGs."))
        gr
      }, error = function(e) {
        message("Error in create_gr_cpgs(): ", e$message)
        NULL
      })
    })
    
    
    gr_offtargets <- reactive({
      # First, req() the outer reactiveVal to ensure it's not NULL
      req(offtarget_table()) 
      
      # Now, get the INNER reactiveVal's value. 
      # This means calling the result of the first call, like this:
      df_offtargets <- offtarget_table()() 
      
      req(df_offtargets) # Ensure the data.frame is not NULL
      
      print("Reactive: gr_offtargets is calculating...")
      tryCatch({
        gr <- create_gr_offtargets(df_offtargets) # Pass the actual data.frame
        print(paste("Reactive: gr_offtargets created with", length(gr), "off-targets."))
        gr
      }, error = function(e) {
        message("Error in gr_offtargets reactive: ", e$message) 
        NULL
      })
    })
    #--------------------------------------------------------------------
    # --- Chromosome Length Management ---
    # This reactive now reflects the chromosome selected in the UI inputs,
    # *not* necessarily the one currently plotted.
    selected_chr_length <- reactive({
      # The chromosome name comes from input$chromosome_input if "Desired targeted region"
      # or from selectedChr() if "DMRs" or "Off-targets"
      # To prevent triggering on selectedChr() before create_plot,
      # let's make this depend on the UI inputs for the *proposed* chromosome.
      chr_name <- if (input$region_type == "Desired targeted region") {
        input$chromosome_input
      } else {
        # When changing region_choice, this logic determines the chromosome for validation
        if (!is.null(input$region_choice) && input$region_choice != "") {
          if (input$region_type == "DMRs" && !is.null(gr_dmrs())) {
            selected_dmr_id <- sub("^(\\S+)\\s.*$", "\\1", input$region_choice)
            dmr_gr <- gr_dmrs()[gr_dmrs()$DMR_ID == selected_dmr_id]
            if (length(dmr_gr) > 0) as.character(seqnames(dmr_gr[1])) else NA_character_
          } else if (input$region_type == "Off-targets" && !is.null(gr_offtargets())) {
            offtarget_gr <- gr_offtargets()[gr_offtargets()$ID == input$region_choice]
            if (length(offtarget_gr) > 0) as.character(seqnames(offtarget_gr[1])) else NA_character_
          } else {
            NA_character_
          }
        } else {
          NA_character_ # Fallback if nothing selected yet
        }
      }
      
      # Crucial: Ensure chr_name is valid before proceeding
      req(chr_name) # Ensure chr_name is not NULL or empty string
      validate(
        need(!is.na(chr_name) && chr_name != "" && chr_name %in% chr_size_df_global$Chromosome,
             paste0("Chromosome '", chr_name, "' not found in chromosome size data."))
      )
      
      print(paste("Reactive: selected_chr_length calculating for", chr_name))
      
      len <- chr_size_df_global %>% dplyr::filter(Chromosome == chr_name) %>% pull(Length)
      print(paste("Reactive: selected_chr_length for", chr_name, "is", len))
      len
    })
    
    
    #----------------------------------------------------------------------
    # observe block that updates *max values* of 'from' and 'to' numeric inputs
    observe({
      chr_len <- selected_chr_length() # This depends on the *proposed* chromosome from UI
      if (!is.na(chr_len) && !is.null(chr_len) && length(chr_len) > 0) { # Add more robust check
        print(paste("Observe: Updating 'from' and 'to' max values to", chr_len))
        updateNumericInput(session, "from", max = chr_len)
        updateNumericInput(session, "to", max = chr_len)
      }
    })
    
    # ---------------------------------------------------------------------------------------
    # --- Dynamic UI for Chromosome based on region_type ---
    output$chromosome_input <- renderUI({
      ns <- session$ns
      chrom_choices <- processed_chroms_list_rv()
      if (input$region_type == "Desired targeted region") {
        selectInput(session$ns("chromosome_input"), 
                    label = "Select Chromosome:",
                    choices = chrom_choices,
                    selected = chrom_choices[1])
      } else {
        print("renderUI: Hiding chromosome_input.")
        NULL
      }
    })
    
    #----------------------------------------------------------------------------------------
    # Handle region_type change: (DMR, Off-targets, disered): 
    #this changes the first region choices by getting it from the table. this triggers the observeEvent(input$region_choice)
    observeEvent(input$region_type, {
      print(paste("ObserveEvent: region_type changed to", input$region_type))
      
      # Always reset inputs (and implicitly the proposed range) on type switch
      updateNumericInput(session, "from", value = 1)
      updateNumericInput(session, "to", value = 10000)
      
      # Reset the *displayed* plot's parameters and tracks
      selectedRange(c(NA, NA))
      selectedChr(NA_character_)
      all_gviz_tracks_rv(NULL) # CRITICAL: Clear the old plot when changing region type
      
      if (input$region_type == "Desired targeted region") {
        updateTextInput(session, "chromosome", value = "")
        print("ObserveEvent: Switched to 'Desired targeted region'. Manual inputs visible, values reset.")
      } else {
        # Ensure the manual inputs are cleared internally even if they become hidden
        session$sendInputMessage("chromosome", list(value = ""))
        
        # If switching to DMRs/Off-targets, programmatically select the first available region.
        # This will then trigger the input$region_choice observer to update 'from'/'to' inputs.
        if (input$region_type == "DMRs") {
          dmrs_data <- isolate(gr_dmrs())
          if (!is.null(dmrs_data) && length(dmrs_data) > 0) {
            first_dmr <- dmrs_data[1]
            
            # Access the chromosome from the object. This might be 'seqnames' or 'chr'
            # Example assuming the columns are named 'DMR_ID', 'seqnames', and 'first_overlapped_gene'
            first_dmr_id_display <- paste0(
              first_dmr$DMR_ID, 
              " (", 
              as.character(first_dmr@seqnames), # or first_dmr$chr
              "_", 
              first_dmr$first_overlapped_gene, 
              ")"
            )
            
            # This update will trigger observeEvent(input$region_choice)
            updateSelectizeInput(session, "region_choice", selected = first_dmr_id_display, server = TRUE)
            print("ObserveEvent: Switched to DMRs, first DMR set as selected_choice (will trigger update).")
          } else {
            updateSelectizeInput(session, "region_choice", choices = character(0), selected = NULL, server = TRUE)
            showNotification("No DMRs found. Please load DMR data.", type = "warning")
            print("ObserveEvent: Switched to DMRs, but no data available to pre-select.")
          }
        }else if (input$region_type == "Off-targets") {
          offtargets_data <- isolate(gr_offtargets())
          if (!is.null(offtargets_data) && length(offtargets_data) > 0) {
            first_offtarget <- offtargets_data[1]
            first_offtarget_id_display <- first_offtarget$ID
            # This update will trigger observeEvent(input$region_choice)
            updateSelectizeInput(session, "region_choice", selected = first_offtarget_id_display, server = TRUE)
            print("ObserveEvent: Switched to Off-targets, first Off-target set as selected_choice (will trigger update).")
          } else {
            updateSelectizeInput(session, "region_choice", choices = character(0), selected = NULL, server = TRUE)
            showNotification("No Off-targets found. Please load off-target data.", type = "warning")
            print("ObserveEvent: Switched to Off-targets, but no data available to pre-select.")
          }
        }
      }
    })
    
    #------------------------------------------------------------------------
    # --- Region Selection Logic (Dropdown) ---
    # handels choices appearing and choosing one of them   
    output$region_selector <- renderUI({
      req(input$region_type)
      
      current_choices <- character(0)
      selected_choice <- NULL # Default to NULL, let the observeEvent handle specific selection
      
      if (input$region_type == "DMRs") {
        print("renderUI: region_selector for type DMRs")
        if (!is.null(gr_dmrs()) && length(gr_dmrs()) > 0) {
          dmrs <- gr_dmrs()
          
          # Construct the gene part of the string conditionally
          gene_part <- ifelse(
            is.na(dmrs$first_overlapped_gene),
            "",  # If gene is NA, add an empty string
            paste0("_", dmrs$first_overlapped_gene) # If gene is not NA, add "_GeneName"
          )
          
          # Combine everything to create the full display string
          current_choices <- paste0(dmrs$DMR_ID, " (", as.character(seqnames(dmrs)), gene_part, ")")
          
          names(current_choices) <- current_choices
          
          # The rest of your code remains the same
          if (!is.null(input$region_choice) && input$region_choice %in% current_choices) {
            selected_choice <- input$region_choice
          } else if (length(current_choices) > 0) {
            selected_choice <- current_choices[1] # Fallback to first if no previous selection
          }
          print(paste("renderUI: DMRs dropdown choices generated. Selected_choice set to:", selected_choice))
          
        } else {
          print("renderUI: No DMRs data available for dropdown.")
        }
      } else if (input$region_type == "Off-targets") {
        print("renderUI: region_selector for type Off-targets")
        if (!is.null(gr_offtargets()) && length(gr_offtargets()) > 0) {
          offtargets <- gr_offtargets()
          current_choices <- offtargets$ID
          names(current_choices) <- current_choices
          # If the input$region_choice already has a value, try to keep it.
          if (!is.null(input$region_choice) && input$region_choice %in% current_choices) {
            selected_choice <- input$region_choice
          } else if (length(current_choices) > 0) {
            selected_choice <- current_choices[1] # Fallback to first if no previous selection
          }
          print(paste("renderUI: Off-targets dropdown choices generated. Selected_choice set to:", selected_choice))
        } else {
          print("renderUI: No Off-targets data available for dropdown.")
        }
      } else if (input$region_type == "Desired targeted region") {
        print("renderUI: Desired targeted region selected, region_selector is NULL.")
        return(NULL) # No dropdown needed for this type
      }
      
      selectInput(session$ns("region_choice"),
                  label =paste("Choose", input$region_type, ":"),
                  choices = current_choices,
                  selected = selected_choice, # This will be the initially selected value
                  multiple = FALSE
      )
    })
    #-------------------------------------------------------------------------------------------------------
    # Handles region_choice changes (by updating 'from' and 'to' inputs)
    observeEvent(input$region_choice, {
      req(input$region_choice)
      # This observer should ONLY run if region_type is DMRs or Off-targets
      req(input$region_type %in% c("DMRs", "Off-targets"))
      
      print(paste("ObserveEvent: region_choice changed to", input$region_choice, "for type", input$region_type))
      
      selected_gr <- NULL
      current_chr_name <- NULL
      if (input$region_type == "DMRs") {
        req(gr_dmrs())
        selected_dmr_id <- str_extract(input$region_choice, "^\\S+")
        selected_gr <- gr_dmrs()[gr_dmrs()$DMR_ID == selected_dmr_id]
        print(paste("ObserveEvent: Selected DMR ID for update:", selected_dmr_id))
      } else if (input$region_type == "Off-targets") {
        req(gr_offtargets())
        selected_offtarget <- input$region_choice
        selected_gr <- gr_offtargets()[gr_offtargets()$ID == selected_offtarget]
        print(paste("ObserveEvent: Selected Off-target ID for update:", selected_offtarget))
      }
      
      if (!is.null(selected_gr) && length(selected_gr) == 1) {
        current_chr_name <- as.character(seqnames(selected_gr))
        pad <- 2000
        # Use chr_size_df_global as it's passed as an argument/global
        chr_len <- chr_size_df_global %>% dplyr::filter(Chromosome == current_chr_name) %>% pull(Length)
        if (length(chr_len) == 0) {
          message("Error: Chromosome length not found for ", current_chr_name)
          showNotification(paste("Chromosome length not found for", current_chr_name), type = "warning")
          return()
        }
        
        new_from <- max(1, start(selected_gr) - pad)
        new_to <- min(chr_len, end(selected_gr) + pad)
        
        # Update the UI numeric inputs directly.
        # These now act as the *proposed* new range/chr, not the *active* plot range.
        updateNumericInput(session, "from", value = new_from)
        updateNumericInput(session, "to", value = new_to)
        # We don't update selectedChr/selectedRange here anymore.
        # They will only be updated when "Create Plot" or zoom/pan is clicked.
        
        print(paste("ObserveEvent: UI 'from' and 'to' updated to", new_from, "-", new_to, "for Chr:", current_chr_name))
        # No longer printing that selectedChr/selectedRange are updated here.
      } else {
        print("ObserveEvent: No matching region found or selected_gr is invalid. Clearing UI inputs.")
        # Clear UI inputs
        updateNumericInput(session, "from", value = 1)
        updateNumericInput(session, "to", value = 10000)
        # No need to reset selectedChr/selectedRange here, as they only reflect the *current plot*
      }
    })
    
    #--------------------------------------------------------------------
    # --- Handle 'Create Plot' button click ---
    observeEvent(input$create_plot, {
      print("ObserveEvent: 'Create Plot' button clicked. Starting track creation.")
      observe({
        if (!is.null(tadcalling_results)) {
          updateCheckboxInput(session, "include_tads", value = TRUE)
        } else {
          updateCheckboxInput(session, "include_tads", value = FALSE)
        }
      })
      # Early validation
      if (isTRUE(input$include_tads)) {
        if (is.null(tadcalling_results)) {
          showNotification("Cannot include TADs - no TAD calling results available",
                           type = "warning")
          updateCheckboxInput(session, "include_tads", value = FALSE)
          return()
        }
      }
      withProgress(message = 'Setting up plot region and creating tracks...', value = 0.1, {
        # step1: set the from, to, chr 
        local_from <- input$from
        local_to <- input$to
        local_chr <- if (input$region_type == "Desired targeted region") {
          # ADDED CHECK if chr is entered by the user or not
          if (is.null(input$chromosome_input) || input$chromosome_input == "") {
            showNotification("Please enter a chromosome name for 'Desired targeted region'.", type = "error", duration = 5)
            return(NULL) # Return NULL to stop further execution in this observeEvent/reactive
          }
          input$chromosome_input
        }else {
          # Get the chromosome from the currently selected region_choice
          if (input$region_type == "DMRs") {
            req(gr_dmrs())
            selected_dmr_id <- sub("^(\\S+)\\s.*$", "\\1", input$region_choice)
            dmr_gr <- gr_dmrs()[gr_dmrs()$DMR_ID == selected_dmr_id]
            req(length(dmr_gr) > 0)
            as.character(seqnames(dmr_gr[1]))
          } else if (input$region_type == "Off-targets") {
            req(gr_offtargets())
            offtarget_gr <- gr_offtargets()[gr_offtargets()$ID == input$region_choice]
            req(length(offtarget_gr) > 0)
            as.character(seqnames(offtarget_gr[1]))
          } else {
            NA_character_ # Should not happen with req above
          }
        }
        req(local_chr) # Ensure we have a chromosome before proceeding
        
        
        #--------------------------------------------
        # steps 2: validate them, get chr length
        # Validate values (keep this important validation)
        if (is.null(local_chr) || is.na(local_chr) || local_chr == "") {
          showNotification("Please specify a chromosome for the plot.", type = "error")
          return()
        }
        if (is.na(local_from) || is.na(local_to) || local_from >= local_to) {
          showNotification("Invalid 'From' or 'To' coordinates, or 'From' must be less than 'To'.", type = "error")
          return()
        }
        
        # Validate against chromosome length
        # Re-fetch chromosome length based on local_chr for robustness
        current_chr_len_for_plot <- chr_size_df_global %>% dplyr::filter(Chromosome == local_chr) %>% pull(Length)
        if (length(current_chr_len_for_plot) == 0 || is.na(current_chr_len_for_plot) || local_from < 1 || local_to > current_chr_len_for_plot) {
          showModal(modalDialog(
            title = "Input Error: Range Out of Bounds",
            paste0("The specified genomic range (", local_from, " - ", local_to, ") is out of bounds for chromosome '", local_chr, "'. ",
                   "The maximum length for this chromosome is: ", current_chr_len_for_plot, "."),
            footer = modalButton("Dismiss"), # A button to close the modal
            easyClose = TRUE # Allows closing by clicking outside the modal or pressing Esc
          ))
          return()
        }
        
        #--------------------------------------
        # Step 3: Load TAD and SubTAD data into local variables into the tables
        incProgress(0.1, detail = "Checking for TAD/SubTAD data")
        # Initialize local variables to NULL
        gr_TADs_local <- NULL
        gr_SUBTADs_local <- NULL
        if (isTRUE(input$include_tads)) {
          # First, check if the reactive expression object itself is not NULL
          if (!is.null(tadcalling_results)) {
            # Now that we know the object exists, we can safely call it.
            if (is.null(tadcalling_results())) {
              gr_TADs(NULL)
              gr_SUBTADs(NULL)
              showNotification("Creating without TAD and SubTAD tracks because they are not available.", type = "warning", duration = 8)
            } else {
              # step a: valisate values
              message("include_tads is TRUE, loading TAD data.")
              current_tissue <- tadcalling_results()$current_tissue()
              current_resolution <- tadcalling_results()$current_resolution()
              current_processed_data_path <- tadcalling_results()$processed_output_dir()
              # step b: check if tad file for the chosen chr is there
              check_results <- find_existing_tads_subtads_chrs_with_check(
                tissue = current_tissue,
                resolution = current_resolution,
                chr = local_chr,
                processed_data_path = current_processed_data_path
              )
              print(paste("find_existing_tads_subtads_chrs_with_check result for", local_chr, ": chr_in_list =", check_results$chr_in_list))
              
              if (check_results$chr_in_list) {
                tad_file_path <- file.path(current_processed_data_path, check_results$matched_tad_filename)
                subtad_file_path <- file.path(current_processed_data_path, check_results$matched_subtad_filename)
                
                tryCatch({
                  tad_tbl <- read_delim(file.path(current_processed_data_path, check_results$matched_tad_filename), delim = "\t", col_names = TRUE, show_col_types = FALSE)
                  subtad_tbl <- read_delim(file.path(current_processed_data_path, check_results$matched_subtad_filename), delim = "\t", col_names = TRUE, show_col_types = FALSE)
                  
                  gr_TADs_local <- create_gr_TADs_SUBTADs(local_chr, tad_tbl)
                  gr_SUBTADs_local <- create_gr_TADs_SUBTADs(local_chr, subtad_tbl)
                  
                  gr_SUBTADs(gr_SUBTADs_local)
                  gr_TADs(gr_TADs_local)
                  
                  showNotification(paste("Loaded existing TAD/SubTAD data for", local_chr), type = "message")
                }, error = function(e) {
                  showNotification(paste("Error reading TAD/SubTAD files:", e$message), type = "error", duration = 8)
                })
              } else {
                showNotification(paste("TAD/SubTAD data not found for", local_chr, ". Not plotting these tracks."),
                                 type = "warning", duration = 8)
                gr_TADs(NULL)
                gr_SUBTADs(NULL)
              }
            }
          } else {
            # This else block is for the check on the tadcalling_results object itself
            gr_TADs(NULL)
            gr_SUBTADs(NULL)
            showNotification("Creating without TAD and SubTAD tracks because the data source is missing.", type = "warning", duration = 8)
          }
        } else {
          # This is the else block for the outer if statement (if (isTRUE(input$include_tads)))
          message("include_Tads is off.")
          gr_TADs(NULL)
          gr_SUBTADs(NULL)
        }
        #------------------------------------------
        #Step4: create the tracks
        incProgress(0.3, detail = "Creating Gviz tracks")
        # This is CRITICAL: Update the reactives that *define the currently plotted region*
        selectedRange(c(local_from, local_to))
        selectedChr(local_chr)
        
        
        # Now, create the tracks and store them in the reactiveVal
        # Pass the locally created GRanges objects directly to the function
        tracks_list <- create_tracks(
          genome = genome,
          chr = local_chr, # Use the validated local_chr
          gr_cpgs = gr_cpgs(),
          gr_cpgIslands = gr_cpgIslands_global,
          dmrs_gr = gr_dmrs(),
          gr_offtargets = gr_offtargets(),
          tx_gr_filtered = tx_gr_filtered_global,
          gr_SUBTADs = gr_SUBTADs_local,
          gr_TADs =gr_TADs_local ,
          binsize = input$binsize,
          num_samples = num_samples(),
          tissue = input$tissue,
          pheno_data = dmr_results()$pheno()
        )
        all_gviz_tracks_rv(tracks_list) # Update the reactiveVal here
        
        incProgress(1, detail = "Ready to plot")
      }) # End of withProgress
    })
    
    
    #---------------------------------------------------------
    # --- Plot Rendering ---
    output$gvizPlot <- renderPlot({
      # This req now waits for all_gviz_tracks_rv to be non-NULL (i.e., updated by create_plot or zoom/pan)
      req(all_gviz_tracks_rv())
      req(selectedChr(), selectedRange()[1], selectedRange()[2]) # Ensure these are also set (should be from create_plot)
      
      # Debug print
      'print(paste("ref_group value:", dmr_results()$ref_group()))
      print(paste("ref_group class:", class(dmr_results()$ref_group())))'
      
      withProgress(message = 'Rendering Gviz Plot...', value = 0, {
        
        incProgress(0.1, detail = "Retrieving plot data")
        tracks_to_plot <- all_gviz_tracks_rv()
        
        # The range and chr should now *always* be consistent with what's in all_gviz_tracks_rv
        from <- selectedRange()[1]
        to <- selectedRange()[2]
        chr <- selectedChr()
        
        # Re-validate the range against the chromosome length for the *current* plot.
        # This should theoretically pass if selectedRange and selectedChr were set correctly by create_plot.
        current_chr_len_for_plot <- chr_size_df_global %>% dplyr::filter(Chromosome == chr) %>% pull(Length)
        
        validate(
          need(from < to, "'From' must be less than 'To'. This indicates an internal error if 'Create Plot' passed."),
          need(!is.na(current_chr_len_for_plot), "Chromosome length not found for current plot. Internal error."),
          need(from >= 1 && to <= current_chr_len_for_plot,
               paste0("Plot range out of bounds for chromosome '", chr, "'. Max length: ", current_chr_len_for_plot, ". Internal error."))
        )
        
        incProgress(0.8, detail = "Drawing plot")
        print(paste0("renderPlot: Calling plotGvizTracks for Chr:", chr, " From:", from, " To:", to))
        plotGvizTracks(tracks = tracks_to_plot, from = from, to = to, pheno = dmr_results()$pheno(), gr_cpgs= gr_cpgs(), ref_group= dmr_results()$ref_group())
        print("renderPlot: Plotting complete.")
        
        incProgress(1, detail = "Done")
      })
    })
    
    
    # --- Plot Download ---
    output$downloadPlot <- downloadHandler(
      filename = function() {
        # It's good practice to ensure selectedChr and selectedRange are valid
        # even though they should be if all_gviz_tracks_rv() is not NULL
        req(selectedChr(), selectedRange()[1], selectedRange()[2])
        paste0("GVizPlot_", selectedChr(), "_", selectedRange()[1], "-", selectedRange()[2], ".pdf")
      },
      content = function(file) {
        print(paste("Download: Initiating PDF download to", file))
        
        # REQUIRE that the tracks have been generated
        # This means the user must have clicked "Create Plot" at least once,
        # or used the zoom/pan buttons.
        req(all_gviz_tracks_rv(), dmr_results()$pheno(), gr_cpgs())
        
        # Get the already generated tracks
        tracks_for_download <- all_gviz_tracks_rv()
        
        # Get the range and chromosome that these tracks were generated for
        current_chr <- selectedChr()
        current_from <- selectedRange()[1]
        current_to <- selectedRange()[2]
        
        # Add a sanity check, though selectedChr/Range should be valid by req(all_gviz_tracks_rv())
        validate(
          need(!is.null(tracks_for_download), "No plot data available for download. Please create a plot first."),
          need(!is.na(current_chr) && !is.null(current_chr), "Chromosome not specified for download."),
          need(!is.na(current_from) && !is.na(current_to) && current_from < current_to, "Invalid plot range for download.")
        )
        
        print("Download: Using existing tracks for PDF.")
        pdf(file = file, width = 12, height = 8) # Specify output device and dimensions
        plotGvizTracks(tracks_for_download, from = current_from, to = current_to, pheno = dmr_results()$pheno(), gr_cpgs=gr_cpgs())
        dev.off()
        print("Download: PDF creation complete.")
      }
    )
    
    #------------------------------------------------------------------------
    # --- Zoom and Pan controls ---
    zoom_factor <- 0.2
    min_range_width <- 1000
    
    # Helper function to update plot and UI after range/chr change
    update_plot_and_ui <- function(new_from, new_to, new_chr) {
      # Update UI inputs
      updateNumericInput(session, "from", value = new_from)
      updateNumericInput(session, "to", value = new_to)
      if (input$region_type == "Desired targeted region") {
        updateTextInput(session, "chromosome", value = new_chr)
      }
      
      # Update the reactives that define the *currently plotted* region
      selectedRange(c(new_from, new_to))
      selectedChr(new_chr)
      
      # Re-create and update the tracks for the new range/chr
      tracks_list <- create_tracks(
        genome = genome,
        chr = new_chr,
        gr_cpgs = gr_cpgs(),
        gr_cpgIslands = gr_cpgIslands_global,
        dmrs_gr = gr_dmrs(),
        gr_offtargets = gr_offtargets(),
        tx_gr_filtered = tx_gr_filtered_global,
        gr_SUBTADs = gr_SUBTADs(),
        gr_TADs = gr_TADs(),
        binsize = input$binsize,
        num_samples = num_samples(),
        tissue = input$tissue
      )
      all_gviz_tracks_rv(tracks_list)
    }
    
    
    observeEvent(input$zoom_in, {
      print("ObserveEvent: 'Zoom In' button clicked.")
      req(selectedRange()[1], selectedRange()[2], selectedChr()) # Ensure a plot is already rendered
      current_chr_len <- chr_size_df_global %>% dplyr::filter(Chromosome == selectedChr()) %>% pull(Length)
      req(current_chr_len) # Ensure chromosome length is available
      
      current_from <- selectedRange()[1]
      current_to <- selectedRange()[2]
      current_width <- current_to - current_from
      mid_point <- (current_from + current_to) / 2
      
      new_width <- max(min_range_width, current_width * (1 - zoom_factor))
      new_from <- floor(mid_point - new_width / 2)
      new_to <- ceiling(mid_point + new_width / 2)
      
      # Apply Chromosome Boundary Constraints
      new_from <- max(1, new_from)
      new_to <- min(current_chr_len, new_to)
      
      if (new_from >= new_to) { # Handle case where zoom-in makes range too small
        new_to <- new_from + min_range_width
        if (new_to > current_chr_len) {
          new_to <- current_chr_len
          new_from <- max(1, new_to - min_range_width)
        }
      }
      print(paste("Zoom In: New range calculated as", new_from, "-", new_to))
      update_plot_and_ui(new_from, new_to, selectedChr())
    })
    
    observeEvent(input$zoom_out, {
      print("ObserveEvent: 'Zoom Out' button clicked.")
      req(selectedRange()[1], selectedRange()[2], selectedChr())
      current_chr_len <- chr_size_df_global %>% dplyr::filter(Chromosome == selectedChr()) %>% pull(Length)
      req(current_chr_len)
      
      current_from <- selectedRange()[1]
      current_to <- selectedRange()[2]
      current_width <- current_to - current_from
      mid_point <- (current_from + current_to) / 2
      
      new_width <- current_width * (1 + zoom_factor)
      new_from <- floor(mid_point - new_width / 2)
      new_to <- ceiling(mid_point + new_width / 2)
      
      new_from <- max(1, new_from)
      new_to <- min(current_chr_len, new_to)
      
      print(paste("Zoom Out: New range calculated as", new_from, "-", new_to))
      update_plot_and_ui(new_from, new_to, selectedChr())
    })
    
    observeEvent(input$go_left, {
      print("ObserveEvent: 'Go Left' button clicked.")
      req(selectedRange()[1], selectedRange()[2], selectedChr())
      current_chr_len <- chr_size_df_global %>% dplyr::filter(Chromosome == selectedChr()) %>% pull(Length)
      req(current_chr_len)
      
      current_from <- selectedRange()[1]
      current_to <- selectedRange()[2]
      current_width <- current_to - current_from
      
      shift_amount <- current_width * zoom_factor
      
      new_from <- max(1, current_from - shift_amount)
      new_to <- new_from + current_width
      
      if (new_to > current_chr_len) { # Prevent going past end of chromosome
        new_to <- current_chr_len
        new_from <- max(1, new_to - current_width)
      }
      print(paste("Go Left: New range calculated as", new_from, "-", new_to))
      update_plot_and_ui(new_from, new_to, selectedChr())
    })
    
    observeEvent(input$go_right, {
      print("ObserveEvent: 'Go Right' button clicked.")
      req(selectedRange()[1], selectedRange()[2], selectedChr())
      current_chr_len <- chr_size_df_global %>% dplyr::filter(Chromosome == selectedChr()) %>% pull(Length)
      req(current_chr_len)
      
      current_from <- selectedRange()[1]
      current_to <- selectedRange()[2]
      current_width <- current_to - current_from
      
      shift_amount <- current_width * zoom_factor
      
      new_to <- min(current_chr_len, current_to + shift_amount)
      new_from <- new_to - current_width
      
      if (new_from < 1) { # Prevent going past start of chromosome
        new_from <- 1
        new_to <- min(current_chr_len, new_from + current_width)
      }
      print(paste("Go Right: New range calculated as", new_from, "-", new_to))
      update_plot_and_ui(new_from, new_to, selectedChr())
    })
    
    observe({
      req(input$region_type)
      if (input$region_type == "Desired targeted region") {
        req(input$from, input$to) # Only require these if inputs are visible
        if (input$from > input$to) {
          showNotification("Start coordinate ('From') cannot be greater than end coordinate ('To').", type = "error", duration = 5)
          print("Validation: 'From' > 'To' detected for 'Desired targeted region'.")
        }
      }
    })
  })
  
} 




