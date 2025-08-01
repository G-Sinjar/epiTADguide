#04_Filtering_Module.R
#Author: Ghazal Sinjar
#Date: 05.06.2025
#Description:
# This Shiny module provides UI and server logic for filtering Illumina EPIC methylation array data.
# Filtering includes detection P-value thresholding, removal of SNP-affected probes, and exclusion of sex chromosome probes.
# It outputs filtered GenomicRatioSet along with beta and M values and saves them automatically.



# ─────────────────────────────────────
# User Interface (UI) FUNCTION
# Creates a sidebar + main panel layout for filtering configuration and output
# ─────────────────────────────────────
#' @param id Module ID.
#' @return A Shiny UI page_sidebar object.
filter_data_ui <- function(id) {
  ns <- NS(id) 
  
  page_sidebar(
    sidebar = sidebar(
      width = "300px",
      
      # Dropdown: normalization method
      selectInput(
        inputId = ns("norm_method"), 
        label = "Choose normalization method:",
        choices = c(
          "Raw" = "raw_normalised",
          "SWAN" = "SWAN",
          "Statified Quantile" = "Quantile",
          "Funnorm" = "Funnorm",
          "Noob" = "Noob",
          "Noob_Swan" = "Noob_Swan"
        ),
        selected = "SWAN"
      ),
      helpText("Select a normalization method to apply before filtering."),
      
      # Checkbox: remove SNPs
      checkboxInput(
        inputId = ns("remove_snps"), 
        label = "Remove SNPs",
        value = TRUE
      ),
      helpText("Probes with known SNPs either at the targeted CpG site or at the Single Base Extension (SBE) site — the base immediately following the CpG site — are flagged for removal."),
      
      br(),
      
      # Sex chromosomes to remove
      checkboxGroupInput(
        inputId = ns("remove_sex_chr"),
        label = "Remove probes on sex chromosomes:",
        choices = c("Remove probes on chrX" = "chrX", "Remove probes on chrY" = "chrY"),
        selected = NULL
      ),
      
      # Action button to start filtering
      actionButton(ns("run_filtering"), "Run Filtering", class = "btn-primary"), 
      
      hr(),
      
      # Output controls: View and download filtered values
      selectInput(ns("value_type"), "Choose value table to view:", 
                  choices = c("Beta values", "M values"),
                  selected = "Beta values"),
      helpText(" Both β-value and M-value have been used as metrics to measure methylation levels. "),
      br(),
      
      # Format selection (shared for both tables)
      radioButtons(ns("table_format"), "Select Download Format:", 
                   choices = c("CSV" = "csv", "Excel" = "xlsx"), 
                   inline = TRUE),
      
      # Download buttons for each table
      downloadButton(ns("download_beta"), "Download Beta Values Table"),
      downloadButton(ns("download_m"), "Download M Values Table")
    ),
    
    # Main content area
    div(
      style = "padding-left: 15px;",
      layout_columns(
        col_widths = c(6, 6),
        card(
          card_title("Filtering Status"),
          verbatimTextOutput(ns("filter_status")) 
        ),
        card(
          card_title("Summary Statistics"),
          tableOutput(ns("filter_stats")) 
        )
      ),
      
      div(
        h5(textOutput(ns("value_table_title"))),
        helpText(textOutput(ns("value_table_help"))), 
        br(),
        DTOutput(ns("value_table")) 
      )
    )
  )
}

# ─────────────────────────────────────
# SERVER FUNCTION 
# Handles filtering steps, updates UI, saves filtered data, provides download handlers
# ─────────────────────────────────────
#' @param id Module ID.
#' @param RGset The RGChannelSet object from minfi.
#' @param raw_normalised The reactive expression holding the raw normalized data.
#' @param normalized_all_methods A reactive expression holding a list containing all normalized methods data.
#' @param project_output_dir A reactive expression for the project's output directory. 
#' @return A list of reactive values: `filtered_data` (the final GenomicRatioSet), `beta_vals`, and `m_vals`.

filter_data_server <- function(id, RGset, raw_normalised, normalized_all_methods, project_output_dir) {
  moduleServer(id, function(input, output, session) {
    shinyjs::useShinyjs()
    
    
    # Reactive containers for status, data tables, filtered object
    status <- reactiveVal("Waiting to start...")
    stats_data <- reactiveVal(NULL)
    beta_table <- reactiveVal(NULL)
    mval_table <- reactiveVal(NULL)
    # Reactive value to store the final filtered GenomicRatioSet
    filtered_grset <- reactiveVal(NULL)
    save_rds_status_msg <- reactiveVal("")
    
    # Disable the save button by default using shinyjs
    shinyjs::disable("save_grset_rds")
    
    # Display current filtering status
    output$filter_status <- renderText({
      status()
    })
    
    output$save_rds_status <- renderText({
      save_rds_status_msg()
    })
    
    # Update save button state based on object availability
    observe({
      if (!is.null(filtered_grset())) {
        shinyjs::enable("save_grset_rds")
      } else {
        shinyjs::disable("save_grset_rds")
      }
    })
    
    #--------------------------------------------------
    # Main filtering procedure
    observeEvent(input$run_filtering, {
      save_rds_status_msg("") 
      filtered_grset(NULL) 
      
      withProgress(message = "Running Filtering Steps...", value = 0, {
        status("Step 1: Filtering CpGs reliable in all samples (p-value > 0.01)")
        incProgress(0.1, detail = "Step 1: Detection P-Values")
        
        # Clear previous results
        stats_data(NULL)
        beta_table(NULL)
        mval_table(NULL)
        
        # Get selected normalized data
        current_normalized_data_value <- if (input$norm_method == "raw_normalised") {
          req(raw_normalised()) 
          raw_normalised()
        } else {
          req(normalized_all_methods()) 
          normalized_all_methods()[[input$norm_method]]()
        }
        
        
        raw_n_initial <- nrow(raw_normalised()) 
        normalized_data_to_filter <- current_normalized_data_value
        
        
        # Step 1: Detection P-value filter
        filtered_by_detectionP <- tryCatch({
          filter_by_detectionP(normalized_data_to_filter, RGset = RGset()) 
        }, error = function(e) {
          status(paste(status(), "\n❌ Error in Step 1:", e$message))
          return(NULL)
        })
        if (is.null(filtered_by_detectionP)) return()
        status(paste(status(), "\n✅ Step 1 is done."))
        reliable_n <- nrow(filtered_by_detectionP)
        
        
        # Step 2: Map to genome
        incProgress(0.2, detail = "Step 2: Mapping to Genome")
        status(paste(status(), "\nStep 2: Mapping to genome"))
        if (inherits(normalized_data_to_filter, "GenomicRatioSet")) {
          mapped <- filtered_by_detectionP
          mapped_n <- nrow(filtered_by_detectionP)
        } else {
          mapped <- tryCatch({
            mapToGenome(filtered_by_detectionP)
          }, error = function(e) {
            status(paste(status(), "\n❌ Error in Step 2:", e$message))
            return(NULL)
          })
          if (is.null(mapped)) return()
          mapped_n <- nrow(mapped)
        }
        status(paste(status(), "\n✅ Step 2 is done."))
        
        
        # Step 3: SNP Removal
        incProgress(0.3, detail = "Step 3: SNP Removal")
        if (isTRUE(input$remove_snps)) {
          status(paste(status(), "\nStep 3: Removing SNPs"))
          ratio_geno <- tryCatch({
            temp <- dropLociWithSnps(mapped, snps = c("SBE", "CpG"), maf = 0)
            temp <- addSnpInfo(temp)
            temp
          }, error = function(e) {
            status(paste(status(), "\n❌ Error in Step 3:", e$message))
            return(NULL)
          })
          if (is.null(ratio_geno)) return()
        } else {
          status(paste(status(), "\nStep 3: Skipping SNP removal"))
          ratio_geno <- addSnpInfo(mapped)
        }
        after_snp_n <- nrow(ratio_geno)
        status(paste(status(), "\n✅ Step 3 is done."))
        
        
        # Step 4: Ratio Conversion
        incProgress(0.5, detail = "Step 4: Ratio Conversion")
        status(paste(status(), "\nStep 4: Getting Beta and M values"))
        if (!inherits(ratio_geno, "GenomicRatioSet")) {
          ratio_geno <- tryCatch({
            ratioConvert(ratio_geno, what = "both", keepCN = TRUE)
          }, error = function(e) {
            status(paste(status(), "\n❌ Error in Step 4:", e$message))
            return(NULL)
          })
          if (is.null(ratio_geno)) return()
        }
        status(paste(status(), "\n✅ Step 4 is done."))
        
        
        # Step 5: Sex Prediction
        incProgress(0.7, detail = "Step 5: Predicting Sex")
        status(paste(status(), "\nStep 5: Predicting sample sex"))
        predictedSex <- suppressWarnings(getSex(ratio_geno, cutoff = -2))
        ratio_geno <- addSex(ratio_geno, sex = predictedSex)
        status(paste(status(), "\n✅ Step 5 is done."))
        
        sex_table <- data.frame(Sample = rownames(predictedSex), Predicted_Sex = predictedSex$predictedSex)
        status(paste(status(), "\n📌 Predicted Sex:\n", paste(capture.output(print(sex_table, row.names = FALSE)), collapse = "\n")))
        
        
        # Step 6: Remove sex chromosome probes
        incProgress(0.9, detail = "Step 6: Removing Sex Chromosome Probes")
        remove_chrs <- input$remove_sex_chr
        if (length(remove_chrs) > 0) {
          status(paste(status(), "\nStep 6: Removing probes on:", paste(remove_chrs, collapse = ", ")))
          gr <-  SummarizedExperiment::rowRanges(ratio_geno)
          ratio_geno <- ratio_geno[!(seqnames(gr) %in% remove_chrs), ]
          status(paste(status(), "\n✅ Removed:", paste(remove_chrs, collapse = ", ")))
        } else {
          status(paste(status(), "\nℹ️ No sex chromosome probes removed"))
        }
        after_sexchr_n <- nrow(ratio_geno)
        
        incProgress(1, detail = "Finalizing tables and object")
        status(paste(status(), "\n✅ Filtering is done!"))
        
        stats <- data.frame(
          "Number of CpGs" = c(raw_n_initial, reliable_n, mapped_n, after_snp_n, after_sexchr_n),
          "Percentage from raw" = round(c(
            100,
            reliable_n / raw_n_initial * 100,
            mapped_n / raw_n_initial * 100,
            after_snp_n / raw_n_initial * 100,
            after_sexchr_n / raw_n_initial * 100
          ), 2)
        )
        rownames(stats) <- c("Raw", "Reliable", "Mapped", "Post-SNP", "Post-SexChr")
        stats_data(stats)
        
        
        
        # Store beta/mval and object
        beta_table(round(getBeta(ratio_geno), 5))
        mval_table(round(getM(ratio_geno), 5))
        filtered_grset(ratio_geno) # Store the final filtered object
        
        
        # Save automatically with descriptive name
        # 1. Normalization method
        norm_method_label <- switch(input$norm_method, 
                                    "raw_normalised" = "Raw",
                                    "SWAN" = "SWAN",
                                    "Quantile" = "Quantile",
                                    "Funnorm" = "Funnorm", 
                                    "Noob" = "Noob",
                                    "Noob_Swan" = "NoobSwan",
                                    "UnknownNorm" # Fallback
        )
        print(paste("Selected norm_method:", input$norm_method, "-> norm_method_label:", norm_method_label))
        
        # 2. SNP removal status
        snp_status_label <- if (isTRUE(input$remove_snps)) "SNPsremoved" else "SNPsKept"

        # 3. Sex chromosome removal status
        sex_chr_status <- if (length(input$remove_sex_chr) == 0) {
          "SexChrProbes_kept"
        } else {
          paste0(paste(input$remove_sex_chr, collapse = "_"), "_removed")
        }
        
        # 4. Date
        current_date <- format(Sys.Date(), "%Y%m%d")
        
        file_name <- paste0(
          "filtered_GRset_",
          norm_method_label, 
          "_", snp_status_label, 
          "_", sex_chr_status,
          "_", current_date,
          ".rds"
        )
        
        output_dir_full <- file.path(project_output_dir(), "intermediate_data") 
        if (!dir.exists(output_dir_full)) {
          dir.create(output_dir_full, recursive = TRUE)
        }
        file_path <- file.path(output_dir_full, file_name)
        
        tryCatch({
          saveRDS(filtered_grset(), file = file_path)
          status(paste0(status(), "\n✅ Filtered GenomicRatioSet automatically saved to: ", file_path))
        }, error = function(e) {
          status(paste0(status(), "\n❌ Error automatically saving RDS: ", e$message))
        })
        # --- End of Automatic Saving ---
      })
    })
    
    # Reactive container for the data table to be displayed
    displayed_table_data <- reactive({
      if (input$value_type == "Beta values") {
        req(beta_table()) # Ensure beta_table has data
        beta_table()
      } else { # M values
        req(mval_table()) # Ensure mval_table has data
        mval_table()
      }
    })
    
    output$value_table <- DT::renderDT({
      req(displayed_table_data()) # Ensure there's data to display
      datatable(
        displayed_table_data(),
        options = list(
          scrollX = TRUE, # Allows horizontal scrolling for wide tables
          pageLength = 10, # Number of rows to display per page
          lengthMenu = c(10, 25, 50, 100), # Options for number of rows
          search = list(regex = FALSE, smart = TRUE), 
          language = list(search = "Search CpGs:")
        ),
        filter = 'top', # Add search filters at the top of each column
        rownames = TRUE # Keep row names (e.g., CpG IDs)
      )
    })
    
    # Render UI elements
    output$filter_stats <- renderTable({
      req(stats_data())
      stats_data()
    }, rownames = TRUE)
    
    output$value_table_title <- renderText({
      req(input$value_type)
      input$value_type
    })
    
    output$value_table_help <- renderText({
      req(input$value_type)
      
      if (input$value_type == "Beta values") {
        "Beta value (β) represents the proportion of methylation at a given CpG site and is widely used due to its biological interpretability. It is calculated as β = Meth / (Meth + Unmeth), where 'Meth' and 'Unmeth' are the intensities measured by methylated and unmethylated probes, respectively. (Xie et al., 2019)"
      } else {
        "M-value offers improved statistical validity for differential methylation analysis. It is calculated as M = log2(β / (1 − β)), transforming the bounded β-values into an unbounded scale more suitable for statistical modeling. (Xie et al., 2019)"
      }
    })
    
    
    # Download Handlers
    output$download_beta <- downloadHandler(
      filename = function() {
        paste0("beta_values.", input$table_format)
      },
      content = function(file) {
        if (input$table_format == "csv") {
          write.csv(beta_table(), file, row.names = FALSE)  
        } else {
          writexl::write_xlsx(beta_table(), path = file)    
        }
      }
    )
    
    output$download_m <- downloadHandler(
      filename = function() {
        paste0("m_values.", input$table_format)
      },
      content = function(file) {
        if (input$table_format == "csv") {
          write.csv(mval_table(), file, row.names = TRUE) 
        } else {
          writexl::write_xlsx(mval_table(), path = file)   
        }
      }
    )
    
    # Return filtered output reactives
    return(list(
      filtered_data = filtered_grset,
      beta_vals = beta_table,
      m_vals = mval_table
    ))
  })
}


'# test module
# ───────────────────────────────────────────────────────────────────────
# LIBRARIES
library(shiny)
library(bslib)
library(writexl)
library(openxlsx)
library(DT)
library(minfi)
library(GenomicRanges)
library(shinyjs)
# ───────────────────────────────────────────────────────────────────────
# Sources:
source("../utils/preprocessing_utils.R")
# ───────────────────────────────────────────────────────────────────────
ui <- page_navbar(
  title = "EPIC Array Pipeline",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  nav_panel("Filtering",
            filter_data_ui("myFilterModule"))
)


# Corrected test module server
server <- function(input, output, session) {
  # Load your data first
  # Ensure these paths are correct for your environment
  all_normalized_methods_data <- readRDS("C:/Users/ghaza/Documents/ghazal/Bioinformatik_Fächer/Masterarbeit_Project/Scripts/R_Scripts/intermediate_data/normalised_all_methods.rds")
  preprocessed_data_object <- readRDS("C:/Users/ghaza/Documents/ghazal/Bioinformatik_Fächer/Masterarbeit_Project/Scripts/R_Scripts/intermediate_data/preprocessed_data.rds")
  project_base_path <- "./main_app_tests"
  
  
  # Create reactive expressions for each required input of the module
  reactive_RGset <- reactive({ preprocessed_data_object$RGset })
  reactive_raw_normalised <- reactive({ preprocessed_data_object$raw_normalised })
  reactive_normalized_all_methods_list <- reactive({ all_normalized_methods_data })
  reactive_project_output_dir <- reactive({ project_base_path }) 
  
  
  
  filtered_output <- filter_data_server(
    "myFilterModule",
    RGset = reactive_RGset,
    raw_normalised = reactive_raw_normalised,
    normalized_all_methods = reactive_normalized_all_methods_list,
    project_output_dir = reactive_project_output_dir
  )
  
}


shinyApp(ui = ui, server = server)'