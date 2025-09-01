# 06_DMR_identification_Module_v1.R
# Author: Ghazal Sinjar
# Date: 04.08.2025
# Description:
# This Shiny module provides UI and server logic for identifying Differentially Methylated Regions (DMRs).
# It utilizes the bumphunter algorithm to detect DMRs based on specified methylation cutoff values and permutation settings.
# The module outputs an annotated table of detected DMRs, which can be downloaded in CSV or Excel format.

#-------------------------------------------------------------------
# ─────────────────────────────────────
# User Interface (UI) FUNCTION
# ─────────────────────────────────────
dmrs_ui <- function(id) {
  ns <- NS(id)
  
  page_sidebar(
    sidebar = sidebar(
      width = "300px",
      
      useShinyjs(), 
      
      # --- NEW UI: Dynamic dropdowns for group selection ---
      uiOutput(ns("ref_group_select")),
      uiOutput(ns("tested_group_select")),
      
      # Input: methylation cutoff range
      numericInput(ns("cutoff_from"), "Cutoff from (min -1 to 0):", value = -0.15, min = -1, max = 0, step= 0.01),
      numericInput(ns("cutoff_to"), "Cutoff to (0 to 1):", value = 0.15, min = 0, max = 1, step = 0.01),
      
      # --- CUTOFF HELP IMPLEMENTATION ---
      actionLink(ns("toggle_cutoff_help"), "Click for help on cutoff."),
      div(
        id = ns("cutoff_help_text_div"),
        style = "display: none;",
        helpText("The cutoff defines how large a difference in Beta-values (methylation ratio) is needed for a region to be considered a candidate. For example, a cutoff of 0.2 means regions with Beta-value differences greater than ±0.2 (i.e., 20%) between groups will be selected. You can specify one value (symmetric around zero) or two values (asymmetric cutoff).")
      ),
      
      #----------------------------------------------------------------
      # Input: number of permutations
      numericInput(ns("B_val"), "Number of permutations (B):", value = 0, min = 0),
      
      # --- PERMUTATIONS HELP IMPLEMENTATION ---
      actionLink(ns("toggle_B_help"), "Click for help on permutations."),
      div(
        id = ns("B_help_text_div"),
        style = "display: none;",
        helpText("B controls the number of permutations used to assess significance, reducing false positives. More permutations = higher accuracy. If the number of samples is large this can be set to a large number, such as 1000. Note that this will take longer.")
      ),
      
      #-------------------------------------------------------
      # Core options radio buttons
      radioButtons(
        ns("core_choice"),
        "Choose number of CPU cores:",
        choices = c(
          "Detected physical cores - 1" = "auto_cores",
          "Choose cores manually" = "manual_cores"
        ),
        selected = "auto_cores"
      ),
      uiOutput(ns("detected_cores_info")),
      # New: Numeric input for manual core selection (dynamically rendered)
      uiOutput(ns("manual_cores_input")),
      
      # --- CORE HELP IMPLEMENTATION START ---
      actionLink(ns("toggle_cores_help"), "Click for help on core choice."),
      div(
        id = ns("cores_help_text_div"),
        style = "display: none;",
        helpText("This setting controls the number of CPU cores used for calculations. Using more cores can significantly speed up the analysis, especially for larger datasets or when running permutations (B > 0).\n\n",
                 "  * **'Detected cores - 1'**: Recommended for most users. This option automatically uses all but one of your computer's available CPU cores. This leaves one core free for system operations, preventing your computer from becoming unresponsive.\n",
                 "  * **'Choose cores manually'**: Allows you to specify an exact number of cores. Use this if you have specific performance requirements or if 'Detected cores - 1' does not suit your needs (e.g., if you want to use fewer cores to save resources for other applications)."
        )
      ),
      
      #-----------------------------------------------------
      # Run button
      actionButton(ns("run_dmr"), "Detect DMRs"),
      helpText("Note: This step can take at least 5 minutes (B=0). Higher B values will increase runtime. For example B=100 takes 1.5 hours using 6 Cores."),
      
      # --- NEW: TABLE NOTES TOGGLEABLE HELP TEXT ---
      actionLink(ns("toggle_table_notes"), "Click for notes on reading the table."),
      div(
        id = ns("table_notes_div"),
        style = "display: none; font-size: 0.9em;",
        h5("Notes on reading the table"),
        p(strong("Width:"), " The length of the DMR in Bp."),
        p(strong("CpGs.inDMR:"), " The number of probes (CpGs) contained within the identified bump."),
        p(strong("CpGs.inCluster:"), " The number of probes (CpGs) around the DMR region, including the CpGs in the DMR."),
        p(strong("first_overlapped_gene:"), " First gene hit in the EnsDb.Hsapiens.v86 (Ensembl based annotation package)."),
        p(strong("All_overlapped_genes:"), " All gene hits in the EnsDb.Hsapiens.v86 (Ensembl based annotation package)."),
        p(strong("*p.value:"), " The unadjusted p-value for the DMR.  A small p-value indicates that a DMR with such a large peak height is unlikely to occur by random chance."),
        p(strong("*fwer:"), " The Family-Wise Error Rate (FWER) adjusted p-value.  It adjusts the p-value to control the probability of making even a single false discovery. A low FWER value (e.g., < 0.05) is strong evidence that the bump is a true finding."),
        p(strong("*p.valueArea:"), "The unadjusted p-value for the bump's area. It answers the question: How likely is it to find a bump with this large of an area by chance? "),
        p(strong("*fwerArea:"), "The FWER-adjusted p-value for the bump's area. It is often a more robust and recommended measure of significance than fwer because it considers both the magnitude and the length (number of CpGs) of the differential methylation. A region with a modest methylation change over many probes might be more biologically significant than a very large change in a single probe, and the area metric captures this."),
        p("(*) for columns which only appear if the number of permutations (B) is not 0.")
      ),
      
      hr(),
      
      # Download options and button (only visible when results available)
      uiOutput(ns("download_ui"))
    ),
    
    # Main output panel
    div(
      style = "padding-left: 15px; padding-right: 15px;",
      layout_columns(
        col_widths = c(12),
        fill = TRUE,
        card(
          card_title("DMR Identification Status"),
          verbatimTextOutput(ns("dmr_status"), placeholder = TRUE)
        )
      ),
      
      # Bottom section: DMR Table
      div(
        style = "margin-top: 20px;",
        h4("Detected DMRs"),
        helpText("Note: A Differentially Methylated Region (DMR) containing just one CpG site is equivalent to a Differentially Methylated Position (DMP)."),
        br(),
        DT::dataTableOutput(ns("dmr_table"))
      )
    )
  )
}


# ─────────────────────────────────────
# SERVER FUNCTION
# ─────────────────────────────────────
#' @param id Module ID.
#' @param filtered_rgset_reactive A reactive expression holding the filtered RGChannelSet or GenomicRatioSet object.
#' @param tx_gr_filtered_static A static (non-reactive) GRanges object of filtered gene transcripts for annotation.
#' @param project_output_dir A reactive expression for the project's output directory.
#' @return A list of reactive values: `dmr_table` (the detected DMRs table) and `pheno` (phenotype data).
dmrs_server <- function(id, filtered_rgset_reactive, tx_gr_filtered_static,project_output_dir) {
  moduleServer(id, function(input, output, session) {
    
    # Reactive value to store final DMR and pheno tables
    dmr_result <- reactiveVal(NULL)
    pheno_result <- reactiveVal(NULL)
    
    # Reactive value to store status text
    dmr_status_text <- reactiveVal("Ready to identify DMRs.")
    
    # Reactive value to store the detected number of cores for display
    detected_cores_display <- reactiveVal(NULL)
    
    #-----------------------------------------------------------------
    # Reactive expression for phenotype data (pData)
    pd_reactive <- reactive({
      req(filtered_rgset_reactive())
      minfi::pData(filtered_rgset_reactive())
    })
    
    # Reactive expression for a design matrix to get column names for the UI
    design_matrix_colnames_reactive <- reactive({
      req(pd_reactive())
      req(input$ref_group)
      
      pd <- pd_reactive()
      sample_group <- pd$Sample_Group
      if (!is.factor(sample_group)) sample_group <- factor(sample_group)
      
      # Check if the reference group exists before releveling
      if (!input$ref_group %in% levels(sample_group)) {
        return(NULL)
      }
      
      sample_group <- relevel(sample_group, ref = input$ref_group)
      designMatrix <- model.matrix(~ sample_group)
      
      # Return all column names except the intercept
      colnames(designMatrix)[-1]
    })
    
    # Dynamic UI for reference group selection
    output$ref_group_select <- renderUI({
      ns <- session$ns
      pd <- pd_reactive()
      req(pd)
      
      sample_groups <- unique(pd$Sample_Group)
      selectInput(ns("ref_group"),
                  "Choose reference sample group:",
                  choices = sample_groups)
    })
    
    # Dynamic UI for tested group selection
    output$tested_group_select <- renderUI({
      ns <- session$ns
      design_colnames <- design_matrix_colnames_reactive()
      req(design_colnames)
      
      selectInput(ns("tested_group"),
                  "Choose group to test against the reference:",
                  choices = design_colnames)
    })
    
    #---------------------------------------------------
    # Display current status
    output$dmr_status <- renderText({dmr_status_text()})
    
    #------------------------------------------------
    # Reactive expression for the number of cores to be used in calculations
    num_cores_to_use <- reactive({
      if (input$core_choice == "auto_cores") {
        if (!requireNamespace("parallel", quietly = TRUE)) {
          warning("Package 'parallel' is required for automatic core detection. Defaulting to 1 core.")
          return(1)
        }
        max(1, parallel::detectCores(logical = FALSE) - 1)
      } else {
        req(input$manual_cores)
        input$manual_cores
      }
    })
    
    #-------------------------------------------------
    # Dynamic UI for manual core selection
    output$manual_cores_input <- renderUI({
      ns <- session$ns
      if (input$core_choice == "manual_cores") {
        max_detected_cores <- if (requireNamespace("parallel", quietly = TRUE)) {
          parallel::detectCores()
        } else {
          4
        }
        tagList(
          numericInput(
            ns("manual_cores"),
            "Number of cores:",
            value = max(1, max_detected_cores - 1),
            min = 1,
            max = max_detected_cores
          ),
          helpText(
            paste0(
              "Note: The default number that appears at first is the number of logical cores detected minus one. ",
              "The detected number of your logical cores is ", max_detected_cores, ", which is the maximum, and it is not recommended to be used."
            )
          )
        )
      }
    })
    
    #---------------------------------------------------
    # Observe core_choice to display detected cores information
    observeEvent(input$core_choice, {
      if (input$core_choice == "auto_cores") {
        current_detected_cores <- if (requireNamespace("parallel", quietly = TRUE)) {
          max(1, parallel::detectCores(logical = FALSE) - 1)
        } else {
          1
        }
        detected_cores_display(paste("The number of physical CPUs/cores detected: ",current_detected_cores+1 , " \n.The number of cores to be used for the analysis: ",current_detected_cores, "."))
      } else {
        detected_cores_display(NULL)
      }
    })
    
    # Render the detected cores information UI
    output$detected_cores_info <- renderUI({
      if (!is.null(detected_cores_display())) {
        helpText(detected_cores_display())
      }
    })
    
    #---------------------------------------------------
    # --- HINT SERVER LOGIC ADDITIONS START ---
    observeEvent(input$toggle_cutoff_help, { shinyjs::toggle("cutoff_help_text_div") })
    observeEvent(input$toggle_B_help, { shinyjs::toggle("B_help_text_div") })
    observeEvent(input$toggle_cores_help, { shinyjs::toggle("cores_help_text_div") })
    observeEvent(input$toggle_table_notes, { shinyjs::toggle("table_notes_div") })
    # --- HINT SERVER LOGIC ADDITIONS END ---
    
    #---------------------------------------------------
    # Helper function to generate a consistent filename base
    get_dmr_base_filename <- function(cutoff_from, cutoff_to, B_val, ref_group, tested_group) {
      paste0("DMRs_" ,ref_group, "_vs_", tested_group,
             "_cutoff_", format(cutoff_from, nsmall = 2), "_",
             format(cutoff_to, nsmall = 2), "_B", B_val)
    }
    
    #------------------------------------------------------------
    # Main observer: Run on "Detect DMRs" button click
    observeEvent(input$run_dmr, {
      req(filtered_rgset_reactive(), tx_gr_filtered_static, input$ref_group, input$tested_group)
      
      current_num_cores <- num_cores_to_use()
      if (is.null(current_num_cores) || current_num_cores < 1) {
        dmr_status_text("❌ Error: Invalid number of cores selected.")
        return()
      }
      
      message(paste("Using", current_num_cores, "cores for bumphunter analysis."))
      dmr_status_text(paste0("Step 1: Identifying DMRs using ", current_num_cores, " cores..."))
      
      withProgress(message = "Running DMR detection...", value = 0, {
        
        # Step 1: get inputs and validate them
        ## current_rgset
        current_rgset <- filtered_rgset_reactive()
        if (is.null(current_rgset)) {
          dmr_status_text("❌ Error: Input data 'filtered_rgset' is missing.")
          return()
        }
        
        ## both cutoff values
        if (input$cutoff_from >= input$cutoff_to) {
          dmr_status_text("❌ Cutoff 'from' must be less than 'to'.")
          return()
        }
        cutoff_vals <- c(input$cutoff_from, input$cutoff_to)
        
        ## B value
        B_val <- input$B_val
        
        dmr_status_text(paste0(
          dmr_status_text(), 
          "\nComparison: '", input$tested_group, "' vs. reference '", input$ref_group, "'."
        ))
        
        # Step 2: run bumphunter
        incProgress(0.1, detail = "Running bumphunter...")
        dmrs_step1_result <- tryCatch({
          run_bumphunter_dmrs(
            rgSet = current_rgset,
            ref_sample_group = input$ref_group,
            tested_group_colname = input$tested_group,
            cutoff = cutoff_vals,
            B = B_val,
            num_cores = current_num_cores
          )
        }, error = function(e) {
          msg <- e$message
          if (grepl("Schreibfehler in Verbindung", msg)) {
            msg <- "A connection write error during the parallelized bumphunter.\nLikely causes: Out of memory (RAM) or worker crash during parallel execution."
          }
          dmr_status_text(paste(dmr_status_text(), "\n❌ Error during DMR detection:\n", msg))
          return(NULL)
        })
        
        ## check bumphunter results
        if (is.null(dmrs_step1_result)) {
          dmr_status_text(paste(dmr_status_text(), "\n❌ Step 1 failed."))
          return()
        }
        if (nrow(dmrs_step1_result) == 0) {
          dmr_status_text(paste(dmr_status_text(), "\n❌ No DMRs found."))
          return()
        }
        dmr_status_text(paste0(dmr_status_text(), sprintf("\n✅ Step 1 completed: %d DMRs found.", nrow(dmrs_step1_result))))
        
        # Step 3: convert DMRs to genomic ranges
        dmr_status_text(paste0(dmr_status_text(), "\nStep 2: Converting DMRs to genomic ranges (GRanges)..."))
        incProgress(0.3, detail = "Converting to GRanges...")
        GR_Dmrs <- tryCatch({
          prepare_dmrs_granges(dmrs_step1_result)
        }, error = function(e) {
          dmr_status_text(paste(dmr_status_text(), "\n❌ Error in Step 2:", e$message))
          return(NULL)
        })
        
        if (is.null(GR_Dmrs)) {
          return()
        }
        dmr_status_text(paste0(dmr_status_text(), "\n✅ Step 2 completed!"))
        
        # Step 4: Annotate DMRs
        dmr_status_text(paste0(dmr_status_text(), "\nStep 3: Annotating DMRs with gene information..."))
        incProgress(0.6, detail = "Annotating with genes...")
        GR_Dmrs_annotated <- tryCatch({
          annotate_dmrs_with_genes(GR_Dmrs, tx_gr_filtered_static)
        }, error = function(e) {
          dmr_status_text(paste(dmr_status_text(), "\n❌ Error in Step 3 (annotation):", e$message))
          return(NULL)
        })
        
        if (is.null(GR_Dmrs_annotated)) {
          return()
        }
        
        num_annotated <- sum(!is.na(mcols(GR_Dmrs_annotated)$first_overlapped_gene))
        dmr_status_text(paste0(
          dmr_status_text(),
          "\n✅ Step 3 completed: ", num_annotated, " DMRs overlapped with genes."
        ))
        
        # Step 5: adjust DMR table to output to user
        dmr_status_text(paste0(dmr_status_text(), "\nStep 4: Preparing final DMR table..."))
        incProgress(0.9, detail = "Finalizing table...")
        
        df <- as.data.frame(GR_Dmrs_annotated)
        if ("seqnames" %in% colnames(df)) {
          colnames(df)[colnames(df) == "seqnames"] <- "chr"
        }
        if ("DMR_ID" %in% colnames(df)) {
          df <- df[, c("DMR_ID", base::setdiff(colnames(df), "DMR_ID"))]
        }
        dmr_result(df)
        dmr_status_text(paste0(dmr_status_text(), "\n✅ Step 4 completed: DMR table is ready for download."))
        
        # Step 6: Save full results as rds automatically
        base_filename <- get_dmr_base_filename(input$cutoff_from, input$cutoff_to, input$B_val, input$ref_group, input$tested_group)
        output_dir_full <- file.path(project_output_dir(), "DMR_results")
        if (!dir.exists(output_dir_full)) dir.create(output_dir_full, recursive = TRUE)
        output_path <- file.path(output_dir_full, paste0(base_filename, ".rds"))
        
        # Get the phenotype data for saving
        pheno_data <- pd_reactive()
        
        saveRDS(list(
          dmr_table = dmr_result(),
          pheno_data = pheno_data
        ), file = output_path)
        dmr_status_text(paste0(dmr_status_text(), "\n📁 Saved full results automatically to ", output_path))
      })
    })
    
    #-------------------------------------------------
    # render the table in the main panel
    output$dmr_table <- DT::renderDataTable({
      req(dmr_result())
      datatable(dmr_result(), options = list(scrollX = TRUE,
                                             pageLength = 10,
                                             autoWidth = TRUE))
    })
    
    #-------------------------------------------------
    # view the download button as soon as DMR table is ready
    output$download_ui <- renderUI({
      req(dmr_result())
      tagList(
        radioButtons(session$ns("download_format"), "Choose format:", choices = c("CSV", "Excel"), inline = TRUE),
        downloadButton(session$ns("download_dmr"), "Download DMR Table")
      )
    })
    
    #-------------------------------------------------
    # Handel clicking on Download button
    output$download_dmr <- downloadHandler(
      filename = function() {
        base_name <- get_dmr_base_filename(input$cutoff_from, input$cutoff_to, input$B_val, input$ref_group, input$tested_group)
        ext <- if (input$download_format == "Excel") "xlsx" else "csv"
        paste0(base_name, "_", Sys.Date(), ".", ext)
      },
      content = function(file) {
        df <- dmr_result()
        if (input$download_format == "Excel") {
          if (!requireNamespace("openxlsx", quietly = TRUE)) {
            stop("Package 'openxlsx' is required for Excel export. Please install it.")
          }
          openxlsx::write.xlsx(df, file)
        } else {
          write.csv(df, file, row.names = FALSE)
        }
      }
    )
    
    #------------------------------------
    # The server returns the final table and pheno data, ref_group and tested_group for the main app
    return(list(
      dmr_table = dmr_result,
      pheno = pd_reactive,
      ref_group = reactive({ input$ref_group })
    ))
  })
}

# test module in module/test_dmrs.R
