# 10_Gviz_Visualization_Module.R
# Author: Ghazal Sinjar
# Date: 02.07.2025
# Description:
# This Shiny module provides UI and server logic for visualizing genomic regions using the Gviz package.
# It enables dynamic plotting of differentially methylated regions (DMRs) along with gene annotations and other genomic features.
# The module integrates with upstream DMR identification results and uses EnsDb or TxDb objects for gene tracks.
# The visualization includes gene models, ideograms, axis tracks, and optional custom annotation tracks.

# Load Required Libraries for the entire app
library(shiny)
library(shinyWidgets) # For selectizeInput
library(bslib)
library(Gviz)
library(GenomicRanges)
library(IRanges)
library(S4Vectors)
library(GenomeInfoDb)
library(BiocGenerics)
library(DT) # If you plan to use data tables
library(EnsDb.Hsapiens.v86)
library(readr) # For read_delim
library(dplyr) # For data manipulation




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
      
      # Static display for Binsize
      div(
        class = "mb-3", # Add some margin bottom for spacing
        tags$label("TAD/SubTAD Bin Size (kb):", class = "form-label"),
        uiOutput(ns("display_binsize")) # This UI output will display the static value
      ),
      
      # Static display for Tissue
      div(
        class = "mb-3", # Add some margin bottom for spacing
        tags$label("TAD/SubTAD Tissue:", class = "form-label"),
        uiOutput(ns("display_tissue")) # This UI output will display the static value
      ),
      
      helpText("Note: The boundaries of TADs and sub-TADs are approximate and depend on the resolution of the TAD-calling algorithm."),
      
      actionButton(ns("create_plot"), "Create Plot", class = "btn btn-info"),
      
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
                           annotation_results,
                           offtarget_table,
                           tadcalling_results, 
                           chr_size_df_global,
                           tx_gr_filtered_global,
                           gr_cpgIslands_global,
                           project_output_dir,
                           genome = "hg38") {
  
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    print(paste0("Module Server (ID: ", id, ") started."))
    Tadcalling_TAD_tbl <- reactiveVal(NULL)
    Tadcalling_SubTAD_tbl <- reactiveVal(NULL)
    # --- Reactive Values for Plotting Range ---
    selectedRange <- reactiveVal(c(NA, NA))
    selectedChr <- reactiveVal(NA_character_)
    all_gviz_tracks_rv <- reactiveVal(NULL)
    
    
    tad_calling_module_output <- reactiveVal(NULL)
    # Reactive trigger for the TAD calling module
    # We will increment this to trigger the TAD calling when needed
    trigger_tad_processing_module <- reactiveVal(0)
    
    # 2. Observer for TAD calling module's output (Place it here)
    observeEvent(list(tad_calling_module_output()$tads_table(),
                      tad_calling_module_output()$subtads_table()), {
                        req(tad_calling_module_output()$tads_table())
                        req(tad_calling_module_output()$subtads_table())
                        
                        Tadcalling_TAD_tbl(tad_calling_module_output()$tads_table())
                        Tadcalling_SubTAD_tbl(tad_calling_module_output()$subtads_table())
                        
                        showNotification(paste("TAD calling completed for", selectedChr(),
                                               "Status:", tad_calling_module_output()$status_message()),
                                         type = "message")
                        print(paste("Received TADs (rows):", nrow(Tadcalling_TAD_tbl())))
                        print(paste("Received SubTADs (rows):", nrow(Tadcalling_SubTAD_tbl())))
                        
                      }, ignoreNULL = TRUE, ignoreInit = TRUE)
    
    
    # --- Static Binsize and Tissue Display ---
    # These will now react to the current_resolution and current_tissue from tadcalling_results
    output$display_binsize <- renderUI({
      req(tadcalling_results$current_resolution()) # Get the reactive value from tadcalling_results
      tags$p(paste0(tadcalling_results$current_resolution(), " kb"), class = "form-control-plaintext")
    })
    
    output$display_tissue <- renderUI({
      req(tadcalling_results$current_tissue()) # Get the reactive value from tadcalling_results
      tags$p(tadcalling_results$current_tissue(), class = "form-control-plaintext")
    })
    # --- End Static Binsize and Tissue Display ---
    
    # --- Reactive GRanges Objects from input tables ---
    gr_dmrs <- reactive({
      req(dmr_results())
      req(dmr_results()$dmr_table)
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
      req(annotation_results())
      req(annotation_results()$annotated_table)
      print("Reactive: gr_cpgs is calculating...")
      tryCatch({
        gr <- create_gr_cpgs(annotation_results()$annotated_table())
        print(paste("Reactive: gr_cpgs created with", length(gr), "CpGs."))
        gr
      }, error = function(e) {
        message("Error in create_gr_cpgs(): ", e$message)
        NULL
      })
    })
    
    gr_offtargets <- reactive({
      req(offtarget_table())
      
      df_offtargets <- offtarget_table()#()
      
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
    
    # ADDED: Reactive GRanges for TADs and SubTADs based on Tadcalling_TAD_tbl and Tadcalling_SubTAD_tbl
    gr_TADs <- reactive({
      req(Tadcalling_TAD_tbl()) # Require the data frame to be available
      print("Reactive: gr_TADs is calculating...")
      tryCatch({
        gr <- create_gr_TADs(Tadcalling_TAD_tbl())
        print(paste("Reactive: gr_TADs created with", length(gr), "TADs."))
        gr
      }, error = function(e) {
        message("Error in create_gr_TADs(): ", e$message)
        NULL
      })
    })
    
    gr_SUBTADs <- reactive({
      req(Tadcalling_SubTAD_tbl()) # Require the data frame to be available
      print("Reactive: gr_SUBTADs is calculating...")
      tryCatch({
        gr <- create_gr_SUBTADs(Tadcalling_SubTAD_tbl())
        print(paste("Reactive: gr_SUBTADs created with", length(gr), "SubTADs."))
        gr
      }, error = function(e) {
        message("Error in create_gr_SUBTADs(): ", e$message)
        NULL
      })
    })
    
    #--------------------------------------------------------------------
    # --- Chromosome Length Management ---
    selected_chr_length <- reactive({
      chr_name <- if (input$region_type == "Desired targeted region") {
        input$chromosome
      } else {
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
      req(chr_name)
      validate(
        need(!is.na(chr_name) && chr_name != "" && chr_name %in% chr_size_df_global$Chromosome,
             paste0("Chromosome '", chr_name, "' not found in chromosome size data."))
      )
      
      print(paste("Reactive: selected_chr_length calculating for", chr_name))
      
      len <- chr_size_df_global %>% dplyr::filter(Chromosome == chr_name) %>% pull(Length)
      print(paste("Reactive: selected_chr_length for", chr_name, "is", len))
      len
    })
    
    observe({
      chr_len <- selected_chr_length()
      if (!is.na(chr_len) && !is.null(chr_len) && length(chr_len) > 0) {
        print(paste("Observe: Updating 'from' and 'to' max values to", chr_len))
        updateNumericInput(session, "from", max = chr_len)
        updateNumericInput(session, "to", max = chr_len)
      }
    })
    # ---------------------------------------------------------------------------------------
    # --- Dynamic UI for Chromosome based on region_type ---
    output$chromosome_input <- renderUI({
      ns <- session$ns
      if (input$region_type == "Desired targeted region") {
        print("renderUI: Showing chromosome_input.")
        textInput(ns("chromosome"), "Chromosome:", value = "")
      } else {
        print("renderUI: Hiding chromosome_input.")
        NULL
      }
    })
    #----------------------------------------------------------------------------------------
    # Handle region_type change:
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
        if (input$region_type == "DMRs") {
          dmrs_data <- isolate(gr_dmrs())
          if (!is.null(dmrs_data) && length(dmrs_data) > 0) {
            first_dmr <- dmrs_data[1]
            first_dmr_id_display <- paste0(first_dmr$DMR_ID, " (", first_dmr$overlapped_gene_name, ")")
            updateSelectizeInput(session, "region_choice", selected = first_dmr_id_display, server = TRUE)
            print("ObserveEvent: Switched to DMRs, first DMR set as selected_choice (will trigger update).")
          } else {
            updateSelectizeInput(session, "region_choice", choices = character(0), selected = NULL, server = TRUE)
            showNotification("No DMRs found. Please load DMR data.", type = "warning")
            print("ObserveEvent: Switched to DMRs, but no data available to pre-select.")
          }
        } else if (input$region_type == "Off-targets") {
          offtargets_data <- isolate(gr_offtargets())
          if (!is.null(offtargets_data) && length(offtargets_data) > 0) {
            first_offtarget <- offtargets_data[1]
            first_offtarget_id_display <- first_offtarget$ID
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
    output$region_selector <- renderUI({
      req(input$region_type)
      
      current_choices <- character(0)
      selected_choice <- NULL # Default to NULL, let the observeEvent handle specific selection
      
      if (input$region_type == "DMRs") {
        print("renderUI: region_selector for type DMRs")
        if (!is.null(gr_dmrs()) && length(gr_dmrs()) > 0) {
          dmrs <- gr_dmrs()
          current_choices <- paste0(dmrs$DMR_ID, " (", dmrs$overlapped_gene_name, ")")
          names(current_choices) <- current_choices
          
          if (!is.null(input$region_choice) && input$region_choice %in% current_choices) {
            selected_choice <- input$region_choice
          } else if (length(current_choices) > 0) {
            selected_choice <- current_choices[1]
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
          if (!is.null(input$region_choice) && input$region_choice %in% current_choices) {
            selected_choice <- input$region_choice
          } else if (length(current_choices) > 0) {
            selected_choice <- current_choices[1]
          }
          print(paste("renderUI: Off-targets dropdown choices generated. Selected_choice set to:", selected_choice))
        } else {
          print("renderUI: No Off-targets data available for dropdown.")
        }
      } else if (input$region_type == "Desired targeted region") {
        print("renderUI: Desired targeted region selected, region_selector is NULL.")
        return(NULL)
      }
      
      selectizeInput(
        ns("region_choice"),
        paste("Choose", input$region_type, ":"),
        choices = current_choices,
        selected = selected_choice,
        multiple = FALSE,
        options = list(placeholder = paste("Select a", input$region_type, "ID"))
      )
    })
    #-------------------------------------------------------------------------------------------------------
    # Handles region_choice changes (by updating 'from' and 'to' inputs)
    observeEvent(input$region_choice, {
      req(input$region_choice)
      req(input$region_type %in% c("DMRs", "Off-targets"))
      
      print(paste("ObserveEvent: region_choice changed to", input$region_choice, "for type", input$region_type))
      
      selected_gr <- NULL
      current_chr_name <- NULL
      if (input$region_type == "DMRs") {
        req(gr_dmrs())
        selected_dmr_id <- sub("^(\\S+)\\s.*$", "\\1", input$region_choice)
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
        chr_len <- chr_size_df_global %>% dplyr::filter(Chromosome == current_chr_name) %>% pull(Length)
        if (length(chr_len) == 0) {
          message("Error: Chromosome length not found for ", current_chr_name)
          showNotification(paste("Chromosome length not found for", current_chr_name), type = "warning")
          return()
        }
        
        new_from <- max(1, start(selected_gr) - pad)
        new_to <- min(chr_len, end(selected_gr) + pad)
        
        updateNumericInput(session, "from", value = new_from)
        updateNumericInput(session, "to", value = new_to)
        
        print(paste("ObserveEvent: UI 'from' and 'to' updated to", new_from, "-", new_to, "for Chr:", current_chr_name))
      } else {
        print("ObserveEvent: No matching region found or selected_gr is invalid. Clearing UI inputs.")
        updateNumericInput(session, "from", value = 1)
        updateNumericInput(session, "to", value = 10000)
      }
    })
    
    #--------------------------------------------------------------------
    # --- Handle 'Create Plot' button click ---
    observeEvent(input$create_plot, {
      print("ObserveEvent: 'Create Plot' button clicked. Starting track creation.")
      
      withProgress(message = 'Setting up plot region and creating tracks...', value = 0.1, {
        local_from <- input$from
        local_to <- input$to
        local_chr <- if (input$region_type == "Desired targeted region") {
          if (is.null(input$chromosome) || input$chromosome == "") {
            showNotification("Please enter a chromosome name for 'Desired targeted region'.", type = "error", duration = 5)
            return(NULL)
          }
          input$chromosome
        }else {
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
            NA_character_
          }
        }
        req(local_chr)
        
        if (is.null(local_chr) || is.na(local_chr) || local_chr == "") {
          showNotification("Please specify a chromosome for the plot.", type = "error")
          return()
        }
        if (is.na(local_from) || is.na(local_to) || local_from >= local_to) {
          showNotification("Invalid 'From' or 'To' coordinates, or 'From' must be less than 'To'.", type = "error")
          return()
        }
        
        current_chr_len_for_plot <- chr_size_df_global %>% dplyr::filter(Chromosome == local_chr) %>% pull(Length)
        if (length(current_chr_len_for_plot) == 0 || is.na(current_chr_len_for_plot) || local_from < 1 || local_to > current_chr_len_for_plot) {
          showNotification(paste0("Range out of bounds for chromosome '", local_chr, "'. Max length: ", current_chr_len_for_plot), type = "error")
          return()
        }
        
      
        # START OF NEW TAD/SUBTAD LOGIC
        incProgress(0.3, detail = "Checking for existing TAD/SubTAD data...")
        
        current_tissue <- tadcalling_results$current_tissue()()
        current_resolution <- tadcalling_results$current_resolution()()
        current_processed_data_path <- tadcalling_results$processed_output_dir()
        
        # Step 1: Check for existing TAD/SubTAD files
        check_results <- find_existing_tads_subtads_chrs_with_check(
          tissue = current_tissue,
          resolution = current_resolution,
          chr = local_chr, # The currently selected chromosome for plotting
          processed_data_path = current_processed_data_path
        )
        print(paste("find_existing_tads_subtads_chrs_with_check result for", local_chr, ": chr_in_list =", check_results$chr_in_list))
        
        if (check_results$chr_in_list) {
          incProgress(0.5, detail = "Loading existing TAD/SubTAD tables...")
          # If files exist, read them directly
          tad_file_path <- file.path(current_processed_data_path, check_results$matched_tad_filename)
          subtad_file_path <- file.path(current_processed_data_path, check_results$matched_subtad_filename)
          
          tryCatch({
            df_tad <- read_delim(tad_file_path, delim = "\t", col_names = TRUE, show_col_types = FALSE)
            df_subtad <- read_delim(subtad_file_path, delim = "\t", col_names = TRUE, show_col_types = FALSE)
            
            Tadcalling_TAD_tbl(df_tad)
            Tadcalling_SubTAD_tbl(df_subtad)
            print(paste("Loaded TADs from:", tad_file_path, "Rows:", nrow(df_tad)))
            print(paste("Loaded SubTADs from:", subtad_file_path, "Rows:", nrow(df_subtad)))
            showNotification(paste("Loaded existing TAD/SubTAD data for", local_chr), type = "message")
          }, error = function(e) {
            showNotification(paste("Error reading TAD/SubTAD files:", e$message), type = "error", duration = 8)
            Tadcalling_TAD_tbl(NULL)
            Tadcalling_SubTAD_tbl(NULL)
            print(paste("Error loading files:", e$message))
          })
          
        } else {
          incProgress(0.5, detail = "Running TAD calling module for current chromosome...")
          # If files do not exist, trigger the TAD calling module
          print(paste("Calling tadcalling_processing_server for:", local_chr, current_tissue, current_resolution))
          
          trigger_tad_processing_module(trigger_tad_processing_module() + 1)
          module_output <- callModule(
            tadcalling_processing_server,
            "dynamic_tad_caller_gviz",
            project_output_dir = reactive(project_output_dir()),
            trigger_button = reactive(trigger_tad_processing_module()),
            tissue_r = reactive(current_tissue),
            chromosome_r = reactive(local_chr),
            resolution_r = reactive(current_resolution),
            mcool_path_r = reactive(current_mcool_path),
            java_memory_r = reactive(current_java_memory)
          )
          tad_calling_module_output(module_output) # Store the returned reactive list
        }
          
        incProgress(0.7, detail = "Creating Gviz tracks")
        selectedRange(c(local_from, local_to))
        selectedChr(local_chr)
        
        
        
        
        # Pass the reactive values from tadcalling_results directly to create_tracks
        tracks_list <- create_tracks(
          genome = genome,
          chr = local_chr,
          gr_cpgs = gr_cpgs(),
          gr_cpgIslands = gr_cpgIslands_global,
          dmrs_gr = gr_dmrs(),
          gr_offtargets = gr_offtargets(),
          tx_gr_filtered = tx_gr_filtered_global,
          gr_SUBTADs = gr_SUBTADs(), # Use the reactive that processes the new subtad table
          gr_TADs = gr_TADs(),       # Use the reactive that processes the new tad table
          binsize = current_resolution,
          num_samples = num_samples(),
          tissue = current_tissue
        )
        all_gviz_tracks_rv(tracks_list)
        
        incProgress(1, detail = "Ready to plot")
      })
    })
    
    
    #---------------------------------------------------------
    # --- Plot Rendering ---
    output$gvizPlot <- renderPlot({
      req(all_gviz_tracks_rv())
      req(selectedChr(), selectedRange()[1], selectedRange()[2])
      
      
      withProgress(message = 'Rendering Gviz Plot...', value = 0, {
        
        incProgress(0.1, detail = "Retrieving plot data")
        tracks_to_plot <- all_gviz_tracks_rv()
        
        from <- selectedRange()[1]
        to <- selectedRange()[2]
        chr <- selectedChr()
        
        current_chr_len_for_plot <- chr_size_df_global %>% dplyr::filter(Chromosome == chr) %>% pull(Length)
        
        validate(
          need(from < to, "'From' must be less than 'To'. This indicates an internal error if 'Create Plot' passed."),
          need(!is.na(current_chr_len_for_plot), "Chromosome length not found for current plot. Internal error."),
          need(from >= 1 && to <= current_chr_len_for_plot,
               paste0("Plot range out of bounds for chromosome '", chr, "'. Max length: ", current_chr_len_for_plot, ". Internal error."))
        )
        
        incProgress(0.8, detail = "Drawing plot")
        print(paste0("renderPlot: Calling plotGvizTracks for Chr:", chr, " From:", from, " To:", to))
        plotGvizTracks(tracks = tracks_to_plot, from = from, to = to)
        print("renderPlot: Plotting complete.")
        
        incProgress(1, detail = "Done")
      })
    })
    
    #---------------------------------------------------------------------------------
    # --- Plot Download ---
    output$downloadPlot <- downloadHandler(
      filename = function() {
        req(selectedChr(), selectedRange()[1], selectedRange()[2])
        paste0("GVizPlot_", selectedChr(), "_", selectedRange()[1], "-", selectedRange()[2], ".pdf")
      },
      content = function(file) {
        print(paste("Download: Initiating PDF download to", file))
        
        req(all_gviz_tracks_rv())
        
        tracks_for_download <- all_gviz_tracks_rv()
        
        current_chr <- selectedChr()
        current_from <- selectedRange()[1]
        current_to <- selectedRange()[2]
        
        validate(
          need(!is.null(tracks_for_download), "No plot data available for download. Please create a plot first."),
          need(!is.na(current_chr) && !is.null(current_chr), "Chromosome not specified for download."),
          need(!is.na(current_from) && !is.na(current_to) && current_from < current_to, "Invalid plot range for download.")
        )
        
        print("Download: Using existing tracks for PDF.")
        pdf(file = file, width = 12, height = 8)
        plotGvizTracks(tracks_for_download, from = current_from, to = current_to)
        dev.off()
        print("Download: PDF creation complete.")
      }
    )
    
    #------------------------------------------------------------------------
    # --- Zoom and Pan controls ---
    zoom_factor <- 0.2
    min_range_width <- 1000
    
    update_plot_and_ui <- function(new_from, new_to, new_chr) {
      updateNumericInput(session, "from", value = new_from)
      updateNumericInput(session, "to", value = new_to)
      if (input$region_type == "Desired targeted region") {
        updateTextInput(session, "chromosome", value = new_chr)
      }
      
      selectedRange(c(new_from, new_to))
      selectedChr(new_chr)
      
      tracks_list <- create_tracks(
        genome = genome,
        chr = new_chr,
        gr_cpgs = gr_cpgs(),
        gr_cpgIslands = gr_cpgIslands_global,
        dmrs_gr = gr_dmrs(),
        gr_offtargets = gr_offtargets(),
        tx_gr_filtered = tx_gr_filtered_global,
        gr_SUBTADs = gr_SUBTADs(), # Use the reactive that processes the new subtad table
        gr_TADs = gr_TADs(),       # Use the reactive that processes the new tad table
        binsize = tadcalling_results$current_resolution(), # Get from tadcalling_results
        num_samples = num_samples(),
        tissue = tadcalling_results$current_tissue()       # Get from tadcalling_results
      )
      all_gviz_tracks_rv(tracks_list)
    }
    
    
    observeEvent(input$zoom_in, {
      print("ObserveEvent: 'Zoom In' button clicked.")
      req(selectedRange()[1], selectedRange()[2], selectedChr())
      current_chr_len <- chr_size_df_global %>% dplyr::filter(Chromosome == selectedChr()) %>% pull(Length)
      req(current_chr_len)
      
      current_from <- selectedRange()[1]
      current_to <- selectedRange()[2]
      current_width <- current_to - current_from
      mid_point <- (current_from + current_to) / 2
      
      new_width <- max(min_range_width, current_width * (1 - zoom_factor))
      new_from <- floor(mid_point - new_width / 2)
      new_to <- ceiling(mid_point + new_width / 2)
      
      new_from <- max(1, new_from)
      new_to <- min(current_chr_len, new_to)
      
      if (new_from >= new_to) {
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
      
      if (new_to > current_chr_len) {
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
      
      if (new_from < 1) {
        new_from <- 1
        new_to <- min(current_chr_len, new_from + current_width)
      }
      print(paste("Go Right: New range calculated as", new_from, "-", new_to))
      update_plot_and_ui(new_from, new_to, selectedChr())
    })
    
    observe({
      req(input$region_type)
      if (input$region_type == "Desired targeted region") {
        req(input$from, input$to)
        if (input$from > input$to) {
          showNotification("Start coordinate ('From') cannot be greater than end coordinate ('To').", type = "error", duration = 5)
          print("Validation: 'From' > 'To' detected for 'Desired targeted region'.")
        }
      }
    })
    
  })
}



# Load necessary libraries
library(shiny)
library(GenomicRanges) # For GRanges objects
library(rtracklayer) # For importing BED/GTF/GFF
library(Gviz) # For Gviz plotting
library(EnsDb.Hsapiens.v86) # For gene annotations
library(BSgenome.Hsapiens.UCSC.hg38) # For chromosome lengths
library(dplyr) # For data manipulation
library(readr) # For read_delim
library(bslib) # For layout_sidebar

# Source your UI and utility functions
source("../utils/GVIZ_plot_utils.R")
source("../utils/tadcalling_processing_utils.R")

# This is a static object that will be passed to the DMR module
message("Loading/Preparing tx_gr_filtered for annotation...")
edb <- EnsDb.Hsapiens.v86
options(ucscChromosomeNames = TRUE)
tx_gr <- genes(edb)
tx_gr_filtered_global <- keepSeqlevels(tx_gr, standardChromosomes(tx_gr), pruning.mode = "coarse")
seqlevelsStyle(tx_gr_filtered_global) <- "UCSC"
message("tx_gr_filtered prepared.")
#-------------------------------------------------------------------------
# Chromosome Lengths Table from BSgenome
library(BSgenome.Hsapiens.UCSC.hg38)
hg38 <- BSgenome.Hsapiens.UCSC.hg38
chr_lengths <- seqlengths(hg38)
chr_size_df_global <- data.frame(
  Chromosome = names(chr_lengths),
  Length = as.numeric(chr_lengths),
  stringsAsFactors = FALSE
)
print("Global: chr_size_df_global loaded.")
#------------------------------------------------------------------
# CpG Islands - Loaded once globally
# Uses the loadCpGIslands_gr function from utils/GVIZ_plot_utils.R
gr_cpgIslands_global <- loadCpGIslands_gr(destfile = "cpgIslandExt_hg38.txt.gz")
print(paste("Global: gr_cpgIslands_global loaded. Number of CpG Islands:", length(gr_cpgIslands_global)))
#------------------------------

dmr_list <- reactiveVal(NULL)
annotated <- reactiveVal(NULL)
offtargets_combined <- reactiveVal(NULL)
# SUBTADs <- reactiveVal(NULL) # No longer needed as direct input
# tads <- reactiveVal(NULL)     # No longer needed as direct input
print("Main App: Initializing data loading...")

# Replace with your actual file paths
dmr_list_path <- "C:/Users/ghaza/Documents/ghazal/Bioinformatik_Fächer/Masterarbeit_Project/Scripts/R_Scripts/intermediate_data/DMRs_cutoff_neg0.15_to_0.15_B0_2025-06-24.rds"
annotated_path <- "C:/Users/ghaza/Documents/ghazal/Bioinformatik_Fächer/Masterarbeit_Project/Scripts/R_Scripts/intermediate_data/annotated_object_20250627.rds"
offtargets_path <- "C:/Users/ghaza/Documents/ghazal/Bioinformatik_Fächer/Masterarbeit_Project/Scripts/R_Scripts/modules/intermediate_data/guide1_guide2_guide3_guide4_guide5_guide6.rds"
subtad_path <- "C:/Users/ghaza/Documents/ghazal/Bioinformatik_Fächer/Masterarbeit_Project/Scripts/PythonProject/TADs_CAKI2/3_TADs_as_txt_and_Excel/CAKI2_chr13_25kb_SubTADs_noDuplicates.txt"
tad_path <- "C:/Users/ghaza/Documents/ghazal/Bioinformatik_Fächer/Masterarbeit_Project/Scripts/PythonProject/TADs_CAKI2/3_TADs_as_txt_and_Excel/CAKI2_chr13_25kb_TADs.txt"
processed_output_dir_path <- "C:/Users/ghaza/Documents/ghazal/Bioinformatik_Fächer/Masterarbeit_Project/Scripts/PythonProject/test/TADcaller_Results/TADs_CAKI2/processed_tads"
project_dir = "C:/Users/ghaza/Documents/ghazal/Bioinformatik_Fächer/Masterarbeit_Project/Scripts/R_Scripts/my_epic_test"
# Load the raw data for TADs and SubTADs
raw_subtad_df <- NULL
raw_tad_df <- NULL

if (file.exists(dmr_list_path)) {
  results_dmr <- readRDS(dmr_list_path)
  dmr_list <- reactive({
    list(dmr_table = reactiveVal(results_dmr$dmr_table),
         pheno = reactiveVal(results_dmr$pheno))
  })
  print(paste("Main App: DMRs file loaded from:", dmr_list_path))
} else {
  print(paste("Main App: ERROR - DMRs file not found:", dmr_list_path))
  dmr_list <- reactive({
    list(dmr_table = reactiveVal(data.frame()), pheno = reactiveVal(data.frame()))
  })
}


if (file.exists(annotated_path)) {
  results_anno <- readRDS(annotated_path)
  annotated <- reactive({
    list(annotated_table = reactiveVal(results_anno$annotated_table))
  })
  print(paste("Main App: Annotation file loaded from:", annotated_path))
} else {
  print(paste("Main App: ERROR - Annotation file not found:", annotated_path))
  annotated <- reactive({
    list(annotated_table = reactiveVal(data.frame()))
  })
}


if (file.exists(offtargets_path)) {
  offtargets_combined(readRDS(offtargets_path))
  print(paste("Main App: Off-targets file loaded from:", offtargets_path))
} else {
  print(paste("Main App: ERROR - Off-targets file not found:", offtargets_path))
  offtargets_combined(data.frame()) # Initialize with an empty data frame
}


if (file.exists(subtad_path)) {
  raw_subtad_df <- read_delim(subtad_path, delim = "\t", show_col_types = FALSE)
  print(paste("Main App: SubTADs file loaded from:", subtad_path))
} else {
  print(paste("Main App: ERROR - SubTADs file not found:", subtad_path))
}

if (file.exists(tad_path)) {
  raw_tad_df <- read_delim(tad_path, delim = "\t", show_col_types = FALSE)
  print(paste("Main App: TADs file loaded from:", tad_path))
} else {
  print(paste("Main App: ERROR - TADs file not found:", tad_path))
}

print("Main App: Data loading complete.")

# --- Extract tissue and chromosome from filenames ---
# For subtad_path: "CAKI2_chr13_25kb_SubTADs_noDuplicates.txt"
# For tad_path: "CAKI2_chr13_25kb_TADs.txt"

extract_info_from_path <- function(path) {
  if (!file.exists(path)) {
    return(list(tissue = "Unknown", chrom = "Unknown", resolution = "Unknown"))
  }
  filename <- basename(path)
  parts <- strsplit(filename, "_")[[1]]
  
  tissue <- if (length(parts) > 0) parts[1] else "Unknown"
  chrom <- if (length(parts) > 1) parts[2] else "Unknown"
  
  # Try to extract resolution (e.g., "25kb")
  resolution_match <- regmatches(filename, regexpr("\\d+kb", filename))
  resolution <- if (length(resolution_match) > 0) gsub("kb", "", resolution_match[1]) else "Unknown"
  
  list(tissue = tissue, chrom = chrom, resolution = as.numeric(resolution))
}

subtad_info <- extract_info_from_path(subtad_path)
tad_info <- extract_info_from_path(tad_path)

# You might want to ensure consistency if both files are expected to have the same info
# For this test, we'll take from one (e.g., subtad_info) or ensure they match
current_tissue_val <- subtad_info$tissue
current_chrom_val <- subtad_info$chrom
current_resolution_val <- subtad_info$resolution


# --- Create the dummy tadcalling_results reactive list ---
tadcalling_results_dummy <- list(
  tads_table = reactiveVal(raw_tad_df),
  subtads_table = reactiveVal(raw_subtad_df),
  current_tissue = reactiveVal(current_tissue_val),
  current_chrom = reactiveVal(current_chrom_val),
  current_resolution = reactiveVal(current_resolution_val),
  processed_output_dir = reactiveVal(processed_output_dir_path), # Static string now wrapped in reactiveVal
  java_memory_used = reactiveVal(16), # Static numeric now wrapped in reactiveVal
  status_message = reactiveVal(NULL) # Static NULL now wrapped in reactiveVal
)

ui <- page_navbar(
  title = "EPIC Array Pipeline",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  nav_panel("Final Plot", GvizPlotUI("gviz_test_id"))
)

server <- function(input, output, session) {
  print("Main Server started.")
  
  GvizPlotServer(
    id = "gviz_test_id",
    dmr_results = dmr_list,
    annotation_results = annotated,
    offtarget_table = offtargets_combined,
    tadcalling_results = tadcalling_results_dummy, # Pass the dummy results here
    chr_size_df_global = chr_size_df_global,
    tx_gr_filtered_global = tx_gr_filtered_global,
    gr_cpgIslands_global = gr_cpgIslands_global,
    project_output_dir = reactive({project_dir})
  )
  print("Main Server: GvizPlotServer module called.")
}

shinyApp(ui, server)