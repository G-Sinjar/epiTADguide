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

# Source the utility file (plotting and GRanges creation functions)
# IMPORTANT: Ensure this path is correct relative to where 10_GVIZ_plot_Module.R is run.
# Also, ensure you have manually commented out the 'groups' line in plotGvizTracks in this file.
source("../utils/GVIZ_plot_utils.R")

# --- Global Data Loading (Static and large objects, loaded once) ---

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

# Gene annotations (EnsDb) - Loaded once globally
edb <- EnsDb.Hsapiens.v86
options(ucscChromosomeNames = TRUE)
tx_gr <- genes(edb)
tx_gr_filtered_global <- keepSeqlevels(tx_gr, standardChromosomes(tx_gr), pruning.mode = "coarse")
seqlevelsStyle(tx_gr_filtered_global) <- "UCSC"
print(paste("Global: tx_gr_filtered_global loaded. Number of genes:", length(tx_gr_filtered_global)))

# CpG Islands - Loaded once globally
# Uses the loadCpGIslands_gr function from utils/GVIZ_plot_utils.R
# IMPORTANT: This will download the file to the 'data' directory if it doesn't exist.
gr_cpgIslands_global <- loadCpGIslands_gr(destfile = "cpgIslandExt_hg38.txt.gz")
print(paste("Global: gr_cpgIslands_global loaded. Number of CpG Islands:", length(gr_cpgIslands_global)))


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
      
      # Optional inputs for TADs/Beta values if you want to make them dynamic
      #numericInput(ns("binsize"), "TAD/SubTAD Bin Size (kb):", value = 25, min = 1),
      #textInput(ns("tissue"), "TAD/SubTAD Tissue:", value = "CAKI2"),
      # num_samples is derived from dmr_results, no direct UI input needed for it here,
      # but it's passed to create_tracks.
      
      tags$div(style = "font-size: 15px; color: #555;",
               helpText("Note: The boundaries of TADs and sub-TADs are approximate and depend on the resolution of the TAD-calling algorithm.")
      ),
      
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
                           dmr_results_table,      # Reactive: dmr_list$dmr_table
                           annotation_table,       # Reactive: annotated$annotated_table
                           offtarget_table,        # Reactive: offtargets_combined
                           tad_table,              # Reactive: tads
                           subtad_table,           # Reactive: SUBTADs
                           chr_size_df_global,     # Global: chr_size_df_global (renamed from chr_size_df for clarity of scope)
                           tx_gr_filtered_global,  # Global: tx_gr_filtered_global
                           gr_cpgIslands_global,   # Global: gr_cpgIslands_global
                           genome = "hg38") {
  
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    print(paste0("Module Server (ID: ", id, ") started."))
    
    # --- Reactive Values for Plotting Range ---
    # pendingRange removed
    selectedRange <- reactiveVal(c(NA, NA))
    selectedChr <- reactiveVal(NA_character_) # Reactive value to hold the currently selected chromosome
    print("ReactiveVals: selectedRange and selectedChr initialized.")
    
    
    # --- Reactive GRanges Objects from input tables ---
    gr_dmrs <- reactive({
      req(dmr_results_table())
      req(dmr_results_table()$dmr_table)
      print("Reactive: gr_dmrs is calculating...")
      tryCatch({
        gr <- create_gr_dmrs(dmr_results_table()$dmr_table)
        print(paste("Reactive: gr_dmrs created with", length(gr), "DMRs."))
        gr
      }, error = function(e) {
        message("Error in create_gr_dmrs(): ", e$message)
        NULL
      })
    })
    
    num_samples <- reactive({
      req(dmr_results_table())
      req(dmr_results_table()$pheno)
      req(dmr_results_table()$pheno$Sample_Group)
      print("Reactive: num_samples is calculating...")
      tryCatch({
        n_samples <- length(dmr_results_table()$pheno$Sample_Group)
        print(paste("Reactive: num_samples calculated as", n_samples))
        n_samples
      }, error = function(e) {
        message("Error calculating num_samples: ", e$message)
        0
      })
    })
    
    gr_cpgs <- reactive({
      req(annotation_table())
      req(annotation_table()$annotated_table)
      print("Reactive: gr_cpgs is calculating...")
      tryCatch({
        gr <- create_gr_cpgs(annotation_table()$annotated_table)
        print(paste("Reactive: gr_cpgs created with", length(gr), "CpGs."))
        gr
      }, error = function(e) {
        message("Error in create_gr_cpgs(): ", e$message)
        NULL
      })
    })
    
    gr_offtargets <- reactive({
      req(offtarget_table())
      print("Reactive: gr_offtargets is calculating...")
      tryCatch({
        gr <- create_gr_offtargets(offtarget_table())
        print(paste("Reactive: gr_offtargets created with", length(gr), "off-targets."))
        gr
      }, error = function(e) {
        message("Error in create_gr_offtargets(): ", e$message)
        NULL
      })
    })
    
    # Corrected gr_TADs reactive to use dynamic chromosome
    gr_TADs <- reactive({
      req(tad_table())
      
      current_chr_for_tad <- if (input$region_type == "Desired targeted region") {
        req(input$chromosome) # Ensure input$chromosome has a value
        input$chromosome
      } else {
        req(selectedChr()) # Ensure selectedChr has a value
        selectedChr()
      }
      
      print(paste("Reactive: gr_TADs is calculating for chromosome", current_chr_for_tad, "..."))
      tryCatch({
        gr <- create_gr_TADs_SUBTADs(tad_table(), chr = current_chr_for_tad) # CORRECTED: dynamic chr
        print(paste("Reactive: gr_TADs created with", length(gr), "TADs."))
        gr
      }, error = function(e) {
        message("Error in create_gr_TADs_SUBTADs() for TADs: ", e$message)
        NULL
      })
    })
    
    # Corrected gr_SUBTADs reactive to use dynamic chromosome
    gr_SUBTADs <- reactive({
      req(subtad_table())
      
      current_chr_for_subtad <- if (input$region_type == "Desired targeted region") {
        req(input$chromosome)
        input$chromosome
      } else {
        req(selectedChr())
        selectedChr()
      }
      
      print(paste("Reactive: gr_SUBTADs is calculating for chromosome", current_chr_for_subtad, "..."))
      tryCatch({
        gr <- create_gr_TADs_SUBTADs(subtad_table(), chr = current_chr_for_subtad) # CORRECTED: dynamic chr
        print(paste("Reactive: gr_SUBTADs created with", length(gr), "SubTADs."))
        gr
      }, error = function(e) {
        message("Error in create_gr_TADs_SUBTADs() for SUBTADs: ", e$message)
        NULL
      })
    })
    
    
    # --- Chromosome Length Management ---
    selected_chr_length <- reactive({
      # The chromosome name comes from input$chromosome if "Desired targeted region"
      # or from selectedChr() if "DMRs" or "Off-targets"
      chr_name <- if (input$region_type == "Desired targeted region") {
        input$chromosome
      } else {
        selectedChr()
      }
      
      # Crucial: Ensure chr_name is valid before proceeding
      req(chr_name) # Ensure chr_name is not NULL or empty string
      validate(
        need(!is.na(chr_name) && chr_name != "" && chr_name %in% chr_size_df_global$Chromosome,
             paste0("Chromosome '", chr_name, "' not found in chromosome size data."))
      )
      
      print(paste("Reactive: selected_chr_length calculating for", chr_name))
      
      len <- chr_size_df_global %>% filter(Chromosome == chr_name) %>% pull(Length)
      print(paste("Reactive: selected_chr_length for", chr_name, "is", len))
      len
    })
    
    # observe block that updates max values of 'from' and 'to' numeric inputs
    observe({
      chr_len <- selected_chr_length()
      if (!is.na(chr_len)) {
        print(paste("Observe: Updating 'from' and 'to' max values to", chr_len))
        updateNumericInput(session, "from", max = chr_len)
        updateNumericInput(session, "to", max = chr_len)
      }
    })
    
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
      
      # Always reset reactive ranges on type switch
      # pendingRange removed
      selectedRange(c(NA, NA)) 
      selectedChr(NA_character_) # Reset here
      
      if (input$region_type == "Desired targeted region") {
        updateTextInput(session, "chromosome", value = "")
        updateNumericInput(session, "from", value = 1)
        updateNumericInput(session, "to", value = 1000)
        print("ObserveEvent: Switched to 'Desired targeted region'. Manual inputs visible, values reset.")
        # No need to set selectedChr/selectedRange here, as it will come from user input$chromosome
      } else {
        # Ensure the manual inputs are cleared internally even if they become hidden
        session$sendInputMessage("chromosome", list(value = ""))
        
        # If switching to DMRs/Off-targets, programmatically select the first available region.
        # This will then trigger the input$region_choice observer to update 'from'/'to' and reactives.
        if (input$region_type == "DMRs") {
          # Use isolate() to prevent this section from reacting to gr_dmrs() changing
          dmrs_data <- isolate(gr_dmrs()) 
          if (!is.null(dmrs_data) && length(dmrs_data) > 0) {
            first_dmr <- dmrs_data[1]
            first_dmr_id_display <- paste0(first_dmr$DMR_ID, " (", first_dmr$overlapped_gene_name, ")")
            
            # This update will trigger observeEvent(input$region_choice)
            updateSelectizeInput(session, "region_choice", selected = first_dmr_id_display, server = TRUE)
            print("ObserveEvent: Switched to DMRs, first DMR set as selected_choice (will trigger update).")
          } else {
            updateSelectizeInput(session, "region_choice", choices = character(0), selected = NULL, server = TRUE)
            updateNumericInput(session, "from", value = 1)
            updateNumericInput(session, "to", value = 1000)
            showNotification("No DMRs found. Please load DMR data.", type = "warning")
            print("ObserveEvent: Switched to DMRs, but no data available to pre-select.")
          }
        } else if (input$region_type == "Off-targets") {
          offtargets_data <- isolate(gr_offtargets()) 
          if (!is.null(offtargets_data) && length(offtargets_data) > 0) {
            first_offtarget <- offtargets_data[1]
            first_offtarget_id_display <- first_offtarget$ID
            
            # This update will trigger observeEvent(input$region_choice)
            updateSelectizeInput(session, "region_choice", selected = first_offtarget_id_display, server = TRUE)
            print("ObserveEvent: Switched to Off-targets, first Off-target set as selected_choice (will trigger update).")
          } else {
            updateSelectizeInput(session, "region_choice", choices = character(0), selected = NULL, server = TRUE)
            updateNumericInput(session, "from", value = 1)
            updateNumericInput(session, "to", value = 1000)
            showNotification("No Off-targets found. Please load off-target data.", type = "warning")
            print("ObserveEvent: Switched to Off-targets, but no data available to pre-select.")
          }
        }
      }
    })
    
    
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
          
          # If the input$region_choice already has a value, try to keep it.
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
      
      selectizeInput(
        ns("region_choice"),
        paste("Choose", input$region_type, ":"),
        choices = current_choices,
        selected = selected_choice, # This will be the initially selected value
        multiple = FALSE,
        options = list(placeholder = paste("Select a", input$region_type, "ID"))
      )
    })
    #-------------------------------------------------------------------------------------------------------
    # Handles region_choice changes (by updating 'from' and 'to' inputs, and selectedChr/selectedRange)
    observeEvent(input$region_choice, {
      req(input$region_choice)
      # This observer should ONLY run if region_type is DMRs or Off-targets
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
        # Use chr_size_df_global as it's passed as an argument/global
        chr_len <- chr_size_df_global %>% filter(Chromosome == current_chr_name) %>% pull(Length)
        if (length(chr_len) == 0) {
          message("Error: Chromosome length not found for ", current_chr_name)
          showNotification(paste("Chromosome length not found for", current_chr_name), type = "warning")
          return()
        }
        
        new_from <- max(1, start(selected_gr) - pad)
        new_to <- min(chr_len, end(selected_gr) + pad)
        
        # Update the UI numeric inputs directly
        updateNumericInput(session, "from", value = new_from)
        updateNumericInput(session, "to", value = new_to)
        
        # Also update the reactive values that trigger the plot
        selectedChr(current_chr_name)
        selectedRange(c(new_from, new_to))
        
        print(paste("ObserveEvent: UI 'from' and 'to' updated to", new_from, "-", new_to, "for Chr:", current_chr_name))
        print(paste("ObserveEvent: selectedChr and selectedRange also updated."))
      } else {
        print("ObserveEvent: No matching region found or selected_gr is invalid. Clearing UI inputs and reactives.")
        # Clear UI inputs
        updateNumericInput(session, "from", value = 1)
        updateNumericInput(session, "to", value = 1000)
        # Reset reactives
        selectedChr(NA_character_)
        selectedRange(c(NA, NA))
      }
    })
    
    
    # --- Handle 'Create Plot' button click ---
    observeEvent(input$create_plot, {
      print("ObserveEvent: 'Create Plot' button clicked.")
      
      withProgress(message = 'Setting up plot region...', value = 0.1, {
        local_from <- input$from # Always read from current UI input
        local_to <- input$to     # Always read from current UI input
        local_chr <- if (input$region_type == "Desired targeted region") {
          input$chromosome
        } else {
          selectedChr() # For DMRs/Off-targets, this is already updated by region_choice observer
        }
        
        # Validate derived values
        if (is.null(local_chr) || is.na(local_chr) || local_chr == "") {
          showNotification("Please specify a chromosome for the plot.", type = "error")
          return()
        }
        if (is.na(local_from) || is.na(local_to) || local_from >= local_to) {
          showNotification("Invalid 'From' or 'To' coordinates, or 'From' must be less than 'To'.", type = "error")
          return()
        }
        
        current_chr_len <- selected_chr_length() # This reactive now has proper req and validation
        if (is.na(current_chr_len) || local_from < 1 || local_to > current_chr_len) {
          showNotification(paste0("Range out of bounds for chromosome '", local_chr, "'. Max length: ", current_chr_len), type = "error")
          return()
        }
        
        incProgress(0.5, detail = "Updating plot region")
        
        # These updates will trigger renderPlot and all_gviz_tracks
        selectedRange(c(local_from, local_to))
        selectedChr(local_chr)
        
        incProgress(1, detail = "Ready to plot")
      }) # End of withProgress for create_plot observer
    })
    
    
    # --- Reactive Track Creation ---
    all_gviz_tracks <- reactive({
      req(selectedChr(), selectedRange()[1], selectedRange()[2]) # Now depend on selectedChr reactiveVal
      
      current_chr <- selectedChr()
      current_from <- selectedRange()[1]
      current_to <- selectedRange()[2]
      
      validate(
        need(current_from < current_to, "Invalid plot range detected during track creation (From >= To)."),
        need(!is.na(current_chr), "Chromosome is not set for track creation."),
        need(current_chr %in% chr_size_df_global$Chromosome, paste0("Chromosome '", current_chr, "' not recognized for track creation."))
      )
      
      print(paste0("Reactive: all_gviz_tracks is creating tracks for Chr:", current_chr, " From:", current_from, " To:", current_to))
      
      tracks_list <- create_tracks(
        genome = genome,
        chr = current_chr, # Use current_chr from selectedChr()
        gr_cpgs = gr_cpgs(),
        gr_cpgIslands = gr_cpgIslands_global, # Use global argument
        dmrs_gr = gr_dmrs(),
        gr_offtargets = gr_offtargets(),
        tx_gr_filtered = tx_gr_filtered_global, # Use global argument
        gr_SUBTADs = gr_SUBTADs(),
        gr_TADs = gr_TADs(),
        binsize = input$binsize,
        num_samples = num_samples(),
        tissue = input$tissue
      )
      print(paste("Reactive: all_gviz_tracks created", length(tracks_list), "tracks."))
      tracks_list
    })
    
    
    # --- Plot Rendering ---
    output$gvizPlot <- renderPlot({
      # This req is essential to ensure values are ready before plotting.
      # Without it, the plot might try to render before selectedRange/selectedChr are properly set.
      req(selectedChr(), selectedRange()[1], selectedRange()[2])
      
      # IMPORTANT: Place withProgress directly inside renderPlot
      withProgress(message = 'Rendering Gviz Plot...', value = 0, {
        
        incProgress(0.1, detail = "Retrieving plot range")
        from <- selectedRange()[1]
        to <- selectedRange()[2]
        chr <- selectedChr()
        
        validate(
          need(from < to, "'From' must be less than 'To'. Please correct the range and click Create Plot."),
          need(!is.na(selected_chr_length()), "Please select a valid chromosome."),
          need(from >= 1 && to <= selected_chr_length(),
               paste0("Range out of bounds for chromosome '", chr, "'. Max length: ", selected_chr_length()))
        )
        
        incProgress(0.4, detail = "Creating Gviz tracks")
        # Call the reactive that creates the tracks.
        # This reactive will only re-run if its dependencies change.
        tracks_to_plot <- all_gviz_tracks()
        
        # It's good practice to add a req here as well, in case all_gviz_tracks() somehow returns NULL
        req(tracks_to_plot)
        
        incProgress(0.8, detail = "Drawing plot")
        print(paste0("renderPlot: Calling plotGvizTracks for Chr:", chr, " From:", from, " To:", to))
        plotGvizTracks(tracks = tracks_to_plot, from = from, to = to)
        print("renderPlot: Plotting complete.")
        
        incProgress(1, detail = "Done") # Mark as complete
      }) # End of withProgress
    })
    
    
    # --- Plot Download ---
    output$downloadPlot <- downloadHandler(
      filename = function() {
        paste0("GVizPlot_", selectedChr(), "_", selectedRange()[1], "-", selectedRange()[2], ".pdf")
      },
      content = function(file) {
        print(paste("Download: Initiating PDF download to", file))
        req(selectedChr(), selectedRange()[1], selectedRange()[2])
        
        tracks_for_download <- create_tracks(
          genome = genome,
          chr = selectedChr(), # Use selectedChr() for download
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
        print("Download: Tracks re-created for PDF.")
        plotGvizTracks(tracks_for_download, from = selectedRange()[1], to = selectedRange()[2])
        dev.off()
        print("Download: PDF creation complete.")
      }
    )
    
    # --- Zoom and Pan controls ---
    zoom_factor <- 0.2
    min_range_width <- 1000
    
    observeEvent(input$zoom_in, {
      print("ObserveEvent: 'Zoom In' button clicked.")
      req(selectedRange()[1], selectedRange()[2], selected_chr_length())
      
      current_from <- selectedRange()[1]
      current_to <- selectedRange()[2]
      current_width <- current_to - current_from
      mid_point <- (current_from + current_to) / 2
      
      new_width <- max(min_range_width, current_width * (1 - zoom_factor))
      new_from <- floor(mid_point - new_width / 2)
      new_to <- ceiling(mid_point + new_width / 2)
      
      new_from <- max(1, new_from)
      new_to <- min(selected_chr_length(), new_to)
      
      if (new_from >= new_to) {
        new_to <- new_from + min_range_width
        if (new_to > selected_chr_length()) {
          new_to <- selected_chr_length()
          new_from <- max(1, new_to - min_range_width)
        }
      }
      
      # Update UI inputs always
      updateNumericInput(session, "from", value = new_from)
      updateNumericInput(session, "to", value = new_to)
      selectedRange(c(new_from, new_to))
      print(paste("Zoom In: New range set to", new_from, "-", new_to))
    })
    
    observeEvent(input$zoom_out, {
      print("ObserveEvent: 'Zoom Out' button clicked.")
      req(selectedRange()[1], selectedRange()[2], selected_chr_length())
      
      current_from <- selectedRange()[1]
      current_to <- selectedRange()[2]
      current_width <- current_to - current_from
      mid_point <- (current_from + current_to) / 2
      
      new_width <- current_width * (1 + zoom_factor)
      new_from <- floor(mid_point - new_width / 2)
      new_to <- ceiling(mid_point + new_width / 2)
      
      new_from <- max(1, new_from)
      new_to <- min(selected_chr_length(), new_to)
      
      updateNumericInput(session, "from", value = new_from)
      updateNumericInput(session, "to", value = new_to)
      
      selectedRange(c(new_from, new_to))
      print(paste("Zoom Out: New range set to", new_from, "-", new_to))
    })
    
    observeEvent(input$go_left, {
      print("ObserveEvent: 'Go Left' button clicked.")
      req(selectedRange()[1], selectedRange()[2], selected_chr_length())
      
      current_from <- selectedRange()[1]
      current_to <- selectedRange()[2]
      current_width <- current_to - current_from
      
      shift_amount <- current_width * zoom_factor
      
      new_from <- max(1, current_from - shift_amount)
      new_to <- new_from + current_width
      
      if (new_to > selected_chr_length()) {
        new_to <- selected_chr_length()
        new_from <- max(1, new_to - current_width)
      }
      
      updateNumericInput(session, "from", value = new_from)
      updateNumericInput(session, "to", value = new_to)
      
      selectedRange(c(new_from, new_to))
      print(paste("Go Left: New range set to", new_from, "-", new_to))
    })
    
    observeEvent(input$go_right, {
      print("ObserveEvent: 'Go Right' button clicked.")
      req(selectedRange()[1], selectedRange()[2], selected_chr_length())
      
      current_from <- selectedRange()[1]
      current_to <- selectedRange()[2]
      current_width <- current_to - current_from
      
      shift_amount <- current_width * zoom_factor
      
      new_to <- min(selected_chr_length(), current_to + shift_amount)
      new_from <- new_to - current_width
      
      if (new_from < 1) {
        new_from <- 1
        new_to <- min(selected_chr_length(), new_from + current_width)
      }
      
      updateNumericInput(session, "from", value = new_from)
      updateNumericInput(session, "to", value = new_to)
      
      selectedRange(c(new_from, new_to))
      print(paste("Go Right: New range set to", new_from, "-", new_to))
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


# ─────────────────────────────────────
# Main Shiny App (app.R equivalent)
# ─────────────────────────────────────

dmr_list <- reactiveVal(NULL)
annotated <- reactiveVal(NULL)
offtargets_combined <- reactiveVal(NULL)
SUBTADs <- reactiveVal(NULL)
tads <- reactiveVal(NULL)
print("Main App: Initializing data loading...")
  
# Replace with your actual file paths
dmr_list_path <- "C:/Users/ghaza/Documents/ghazal/Bioinformatik_Fächer/Masterarbeit_Project/Scripts/R_Scripts/intermediate_data/DMRs_cutoff_neg0.15_to_0.15_B0_2025-06-24.rds"
annotated_path <- "C:/Users/ghaza/Documents/ghazal/Bioinformatik_Fächer/Masterarbeit_Project/Scripts/R_Scripts/intermediate_data/annotated_object_20250627.rds"
offtargets_path <- "C:/Users/ghaza/Documents/ghazal/Bioinformatik_Fächer/Masterarbeit_Project/Scripts/R_Scripts/modules/intermediate_data/guide1_guide2_guide3_guide4_guide5_guide6.rds"
subtad_path <- "C:/Users/ghaza/Documents/ghazal/Bioinformatik_Fächer/Masterarbeit_Project/Scripts/PythonProject/TADs_CAKI2/3_TADs_as_txt_and_Excel/CAKI2_chr13_25kb_TADs.txt"
tad_path <- "C:/Users/ghaza/Documents/ghazal/Bioinformatik_Fächer/Masterarbeit_Project/Scripts/PythonProject/TADs_CAKI2/3_TADs_as_txt_and_Excel/CAKI2_chr13_25kb_SubTADs_noDuplicates.txt"
  
  if (file.exists(dmr_list_path)) {
    dmr_list(readRDS(dmr_list_path))
    print(paste("Main App: DMRs file loaded from:", dmr_list_path))
  } else {
    message(paste("DMRs file not found:", dmr_list_path))
    print(paste("Main App: ERROR - DMRs file not found:", dmr_list_path))
  }
  
  if (file.exists(annotated_path)) {
    annotated(readRDS(annotated_path))
    print(paste("Main App: Annotated CpGs file loaded from:", annotated_path))
  } else {
    message(paste("Annotated CpGs file not found:", annotated_path))
    print(paste("Main App: ERROR - Annotated CpGs file not found:", annotated_path))
  }
  
  if (file.exists(offtargets_path)) {
    offtargets_combined(readRDS(offtargets_path))
    print(paste("Main App: Off-targets file loaded from:", offtargets_path))
  } else {
    message(paste("Off-targets file not found:", offtargets_path))
    print(paste("Main App: ERROR - Off-targets file not found:", offtargets_path))
  }
  
  if (file.exists(subtad_path)) {
    SUBTADs(read_delim(subtad_path, delim = "\t", show_col_types = FALSE))
    print(paste("Main App: SubTADs file loaded from:", subtad_path))
  } else {
    message(paste("SubTADs file not found:", subtad_path))
    print(paste("Main App: ERROR - SubTADs file not found:", subtad_path))
  }
  
  if (file.exists(tad_path)) {
    tads(read_delim(tad_path, delim = "\t", show_col_types = FALSE))
    print(paste("Main App: TADs file loaded from:", tad_path))
  } else {
    message(paste("TADs file not found:", tad_path))
    print(paste("Main App: ERROR - TADs file not found:", tad_path))
  }
  print("Main App: Data loading complete.")



ui <- page_navbar(
  title = "EPIC Array Pipeline",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  nav_panel("Final Plot", GvizPlotUI("gviz_test_id"))
)

server <- function(input, output, session) {
  print("Main Server started.")
  
  GvizPlotServer(
    id = "gviz_test_id",
    dmr_results_table = dmr_list,
    annotation_table = annotated,
    offtarget_table = offtargets_combined,
    tad_table = tads,
    subtad_table = SUBTADs,
    chr_size_df = chr_size_df_global,
    tx_gr_filtered = tx_gr_filtered_global,
    gr_cpgIslands = gr_cpgIslands_global
  )
  print("Main Server: GvizPlotServer module called.")
}

shinyApp(ui, server)