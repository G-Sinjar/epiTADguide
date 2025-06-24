# LIBRARIES
library(shiny)
library(bslib)
library(writexl)
library(openxlsx)
library(DT)
library(minfi)
library(GenomicRanges) 

#nputs
# --- GLOBAL DATA LOADING (RUNS ONCE WHEN APP STARTS) ---
normalized_all_methods <- readRDS("C:/Users/ghaza/Documents/ghazal/Bioinformatik_Fächer/Masterarbeit_Project/Scripts/R_Scripts/intermediate_data/normalised_all_methods.rds")
preprocessed <- readRDS("C:/Users/ghaza/Documents/ghazal/Bioinformatik_Fächer/Masterarbeit_Project/Scripts/R_Scripts/intermediate_data/preprocessed_data.rds")

RGset <- preprocessed$RGset
raw_normalised <- preprocessed$raw_normalised

#------------------------------------------------------------------
# function to filter both methylset and GenomicRatioSet
filter_by_detectionP <- function(normalized_object, RGset) {
  
  detectp <- detectionP(RGset)
  unreliable_probes_mask <- rowSums(detectp > 0.01) > 0 
  unreliable_probe_names <- rownames(detectp)[unreliable_probes_mask]
  normalized_object_probe_names <- rownames(normalized_object)
  
  # Check if the normalized object is a minfi S4 object that commonly alters probe sets/order.
  is_minfi_S4_object <- inherits(normalized_object, "GenomicRatioSet") ||
    inherits(normalized_object, "MethylSet") ||
    inherits(normalized_object, "MSet") ||
    inherits(normalized_object, "RGChannelSet")
  
  
  can_do_direct_subset <- !is_minfi_S4_object && isTRUE(all.equal(normalized_object_probe_names, rownames(detectp)))
  if (can_do_direct_subset) {
    message("Applying detection P-value filtering using direct logical subsetting (row names match).")
    # Direct subsetting, assuming the order and presence of probes matches detectp
    filtered_object <- normalized_object[!unreliable_probes_mask, ]
  } else {
    message("Applying detection P-value filtering using probe name intersection (likely S4 object or mismatch).")
    # This is the more robust path, using setdiff for name-based filtering
    probes_to_keep <- setdiff(normalized_object_probe_names, unreliable_probe_names)
    filtered_object <- normalized_object[probes_to_keep, ]
  }
  
  message(paste0("Original probes: ", nrow(normalized_object), 
                 ". Probes after detection P-value filtering: ", nrow(filtered_object), "."))
  
  return(filtered_object)
}
#---------------------------------------------------------------------------------------




ui_filter <- page_navbar(
  title = "EPIC Array Pipeline",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  
  nav_panel(
    "Filtering",
    layout_sidebar(
      sidebar = sidebar(
        width = "300px",
        
        # Normalization method dropdown
        selectInput(
          inputId = "norm_method",
          label = "Choose normalization method:",
          choices = c(
            "Raw" = "raw_normalised", # This refers to preprocessed$raw_normalised
            "SWAN" = "SWAN",
            "Quantile" = "Quantile",
            "Functional" = "Funnorm",
            "Noob" = "Noob",
            "Noob_Swan" = "Noob_Swan"
          ),
          selected = "SWAN"
        ),
        helpText("Select a normalization method to apply before filtering.")
        ,
        
        
        # SNPs removal toggle
        checkboxInput(
          inputId = "remove_snps",
          label = "Remove SNPs",
          value = TRUE
        ),
        helpText("Toggle on to remove SNP-affected probes."),
        
        br(),
        
        # Sex chromosomes removal checkboxes
        checkboxGroupInput(
          inputId = "remove_sex_chr",
          label = "Remove probes on sex chromosomes:",
          choices = c("Remove probes on chrX" = "chrX", "Remove probes on chrY" = "chrY"),
          selected = NULL
        ),
        
        
        # Run filtering button
        actionButton("run_filtering", "Run Filtering", class = "btn-primary"),
        
        hr(),
        
        # Value type selection
        selectInput("value_type", "Choose value table to view:",
                    choices = c("Beta values", "M values"),
                    selected = "Beta values"),
        br(),
        
        # Beta values download
        selectInput("beta_format", "Choose file format to download:", choices = c("CSV" = "csv", "Excel" = "xlsx")),
        downloadButton("download_beta", "Download Beta Values table"),
        
        
        # M values download
        selectInput("mval_format", "Choose file format to download:", choices = c("CSV" = "csv", "Excel" = "xlsx")),
        downloadButton("download_m", "Download M Values table")
      ),
      
      # Removed scroll styling here — let page scroll naturally
      div(
        style = "padding-left: 15px;",
        layout_columns(
          col_widths = c(6, 6),
          card(
            card_title("Filtering Status"),
            verbatimTextOutput("filter_status")
          ),
          card(
            card_title("Summary Statistics"),
            tableOutput("filter_stats")
          )
        ),
        
        div(
          h5(textOutput("value_table_title")),
          tags$small(em("Use | for searching multiple CpGs. e.g: cg00000029|cg00000108")),
          br(), br(),
          DTOutput("value_table")
        )
        
      )
    )
  )
)



server_filter <- function(input, output, session) {
  
  status <- reactiveVal("Waiting to start...")
  stats_data <- reactiveVal(NULL)
  beta_table <- reactiveVal(NULL)
  mval_table <- reactiveVal(NULL)
  
  output$filter_status <- renderText({
    status()
  })
  
  observeEvent(input$run_filtering, {
    withProgress(message = "Running Filtering Steps...", value = 0, {
      status("Step 1: Filtering CpGs reliable in all samples (p-value > 0.01)")
      incProgress(0.1, detail = "Step 1: Detection P-Values")
      
      # Step 1: Detection P
      stats_data(NULL)
      beta_table(NULL)
      mval_table(NULL)
      
      current_normalized_data <- reactive({
        if (input$norm_method == "raw_normalised") {
          preprocessed$raw_normalised
        } else {
          normalized_all_methods[[input$norm_method]]
        }
      })
      
      raw_n_initial <- nrow(raw_normalised)
      normalized_data_to_filter <- current_normalized_data()
      
      filtered_by_detectionP <- tryCatch({
        filter_by_detectionP(normalized_data_to_filter, RGset = RGset)
      }, error = function(e) {
        status(paste(status(), "\n❌ Error in Step 1:", e$message))
        return(NULL)
      })
      if (is.null(filtered_by_detectionP)) return()
      status(paste(status(), "\n✅ Step 1 is done."))
      reliable_n <- nrow(filtered_by_detectionP)
      
      incProgress(0.2, detail = "Step 2: Mapping to Genome")
      status(paste(status(), "\nStep 2: Mapping to genome"))
      if (inherits(normalized_data_to_filter, "GenomicRatioSet")) {
        mapped <- filtered_by_detectionP
        mapped_n <- nrow(filtered_by_detectionP)
      } else {
        mapped <- tryCatch({
          mapToGenome(filtered_by_detectionP, mergeManifest = TRUE)
        }, error = function(e) {
          status(paste(status(), "\n❌ Error in Step 2:", e$message))
          return(NULL)
        })
        if (is.null(mapped)) return()
        mapped_n <- nrow(mapped)
      }
      status(paste(status(), "\n✅ Step 2 is done."))
      
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
      
      incProgress(0.7, detail = "Step 5: Predicting Sex")
      status(paste(status(), "\nStep 5: Predicting sample sex"))
      predictedSex <- suppressWarnings(getSex(ratio_geno, cutoff = -2))
      ratio_geno <- addSex(ratio_geno, sex = predictedSex)
      status(paste(status(), "\n✅ Step 5 is done."))
      
      sex_table <- data.frame(Sample = rownames(predictedSex), Predicted_Sex = predictedSex$predictedSex)
      status(paste(status(), "\n📌 Predicted Sex:\n", paste(capture.output(print(sex_table, row.names = FALSE)), collapse = "\n")))
      
      incProgress(0.9, detail = "Step 6: Removing Sex Chromosome Probes")
      remove_chrs <- input$remove_sex_chr
      if (length(remove_chrs) > 0) {
        status(paste(status(), "\nStep 6: Removing probes on:", paste(remove_chrs, collapse = ", ")))
        gr <- rowRanges(ratio_geno)
        ratio_geno <- ratio_geno[!(seqnames(gr) %in% remove_chrs), ]
        status(paste(status(), "\n✅ Removed:", paste(remove_chrs, collapse = ", ")))
      } else {
        status(paste(status(), "\nℹ️ No sex chromosome probes removed"))
      }
      after_sexchr_n <- nrow(ratio_geno)
      
      incProgress(1, detail = "Finalizing tables")
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
      
      beta_table(round(getBeta(ratio_geno), 5))
      mval_table(round(getM(ratio_geno), 5))
    })
  })
  
  
  # Status and output rendering
  output$filter_stats <- renderTable({
    req(stats_data())
    stats_data()
  }, rownames = TRUE)
  
  output$value_table_title <- renderText({
    req(input$value_type)
    input$value_type
  })
  
  output$value_table <- renderDT({
    req(input$value_type)
    table <- if (input$value_type == "Beta values") {
      beta_table()
    } else {
      mval_table()
    }
    req(table)  # <- ensures table is not NULL
    datatable(table,
              options = list(pageLength = 10, scrollX = TRUE),
              rownames = TRUE)
  })
  
  
  # Download handlers
  output$download_stats <- downloadHandler(
    filename = function() {
      paste0("filtering_stats_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(stats_data(), file, row.names = TRUE)
    }
  )
  
  output$download_beta <- downloadHandler(
    filename = function() {
      ext <- input$beta_format
      paste0("beta_values_", Sys.Date(), ".", ext)
    },
    content = function(file) {
      table <- beta_table()
      if (input$beta_format == "csv") {
        write.csv(table, file, row.names = TRUE)
      } else if (input$beta_format == "xlsx") {
        openxlsx::write.xlsx(table, file, rowNames = TRUE)
      }
    }
  )
  
  output$download_m <- downloadHandler(
    filename = function() {
      ext <- input$mval_format
      paste0("m_values_", Sys.Date(), ".", ext)
    },
    content = function(file) {
      table <- mval_table()
      if (input$mval_format == "csv") {
        write.csv(table, file, row.names = TRUE)
      } else if (input$mval_format == "xlsx") {
        openxlsx::write.xlsx(table, file, rowNames = TRUE)
      }
    }
  )
}


shinyApp(ui = ui_filter, server = server_filter)