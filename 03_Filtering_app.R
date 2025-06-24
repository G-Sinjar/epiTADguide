# 1st filtering: CPGs reliable in all samples
#detectp <- detectionP(RGset)
#dim(detectp) #936990      8

#if (all.equal(row.names(SWAN_normalised), row.names(detectp))){
#  absent <- apply(detectp, 1, function(x) sum(x>0.01)) 
  #table(absent) #-> 911931 -_>  97.32559%
#}
#SWAN_filtered <- SWAN_normalised[absent == 0] # 911931 


# 2nd filtering: Map to the Genome
#Geno_SWAN_filtered <- mapToGenome(SWAN_filtered, mergeManifest = TRUE)
#dim(Geno_SWAN_filtered) #905075

# 3rd filtering: Convert GenoMethylSet to a GenoRatioSet (Beta + M): 
#ratio_geno_Swan <- ratioConvert(Geno_SWAN_filtered,what = "both", keepCN= TRUE)
#dim(ratio_geno_Swan) #905075      8

# 4th filtering: removing SNPs 
#snps <- getSnpInfo(ratio_geno_Swan)
#ratio_geno_Swan <- addSnpInfo(ratio_geno_Swan)
#ratio_geno_Swan_NoSNP <- dropLociWithSnps(ratio_geno_Swan, snps = c("SBE", "CpG"), maf = 0)
#dim(ratio_geno_Swan_NoSNP)#891492    8

# 5th filtering: # Get annotation data
#annotation <- getAnnotation(ratio_geno_Swan_NoSNP)
#dim(annotation) #891492     42 -> 100%

# filtering stats
#beta vallues
# m vallues

library(shiny)
library(bslib)
library(writexl)


ui_filter <- page_navbar(
  title = "EPIC Array Pipeline",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  
  nav_panel(
    "Filtering",
    layout_sidebar(
      sidebar = sidebar(
        width = "300px",
        actionButton("run_filtering", "Run Filtering"),
        br(), br(),
        uiOutput("download_stats_ui"),
        br(),br(),
        selectInput("value_type", "Choose value table:", 
                    choices = c("Beta values", "M values"),
                    selected = "Beta values"),
        downloadButton("download_beta", "Download Beta Values Table"),
        downloadButton("download_m", "Download M Values Table")
      ),
      layout_columns(
        col_widths = c(6, 6),
        
        # Left: Filtering Status
        card(
          card_title("Filtering Status"),
          style = "min-height: 400px; max-height: 400px; overflow-y: auto;",
          verbatimTextOutput("filter_status")
        ),
        
        # Right: Summary
        card(
          card_title("Summary Statistics"),
          style = "min-height: 400px; max-height: 400px; overflow-y: auto;",
          tableOutput("filter_stats")
        )
      ),
      br(),
      # Below both: Beta/M values table
      card(
        card_title(textOutput("value_table_title")),
        tableOutput("value_table"),
        helpText("Displaying the first 10 rows only.")
        
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
    # Reset status
    status("Step 1: Filtering CpGs reliable in all samples (p-value>0.01)")
    
    detectp <- tryCatch({
      detectionP(RGset)
    }, error = function(e) {
      status(paste(status(), "\n❌ Error in Step 1:", e$message))
      return(NULL)
    })
    if (is.null(detectp)) return()
    
    if (!all.equal(row.names(SWAN_normalised), row.names(detectp))) {
      status(paste(status(), "\n❌ Row names do not match in Step 1."))
      return()
    }
    
    absent <- apply(detectp, 1, function(x) sum(x > 0.01))
    SWAN_filtered <- SWAN_normalised[absent == 0, ]
    status(paste(status(), "\n✅ Step 1 is done.", sep = ""))
    
    
    status(paste(status(), "\nStep 2: Mapping to genome"))
    Geno_SWAN_filtered <- tryCatch({
      mapToGenome(SWAN_filtered, mergeManifest = TRUE)
    }, error = function(e) {
      status(paste(status(), "\n❌ Error in Step 2:", e$message))
      return(NULL)
    })
    if (is.null(Geno_SWAN_filtered)) return()
    status(paste(status(), "\n✅ Step 2 is done.", sep = ""))
    
    
    status(paste(status(), "\nStep 3: Getting Beta and M values"))
    ratio_geno_Swan <- tryCatch({
      ratioConvert(Geno_SWAN_filtered, what = "both", keepCN = TRUE)
    }, error = function(e) {
      status(paste(status(), "\n❌ Error in Step 3:", e$message))
      return(NULL)
    })
    if (is.null(ratio_geno_Swan)) return()
    status(paste(status(), "\n✅ Step 3 is done.", sep = ""))
    
    
    status(paste(status(), "\nStep 4: Removing SNPs"))
    ratio_geno_Swan_NoSNP <- tryCatch({
      snps <- getSnpInfo(ratio_geno_Swan)
      ratio_geno_Swan <- addSnpInfo(ratio_geno_Swan)
      dropLociWithSnps(ratio_geno_Swan, snps = c("SBE", "CpG"), maf = 0)
    }, error = function(e) {
      status(paste(status(), "\n❌ Error in Step 4:", e$message))
      return(NULL)
    })
    if (is.null(ratio_geno_Swan_NoSNP)) return()
    status(paste(status(), "\n✅ Step 4 is done.", sep = ""))
    
    status(paste(status(), "\n🎉 Filtering complete!", sep = ""))
    
    # ✅ Create and assign stats table here
    raw_n <- nrow(raw_normalised)
    reliable_n <- nrow(SWAN_filtered)
    mapped_n <- nrow(Geno_SWAN_filtered)
    no_snp_n <- nrow(ratio_geno_Swan_NoSNP)
    snp_n <- nrow(ratio_geno_Swan) - no_snp_n
    
    stats <- data.frame(
      "Number of CpGs" = c(raw_n, reliable_n, mapped_n, snp_n, no_snp_n),
      "Percentage from raw" = round(c(
        100,
        reliable_n / raw_n * 100,
        mapped_n / raw_n * 100,
        snp_n / raw_n * 100,
        no_snp_n / raw_n * 100
      ), 2)
    )
    
    rownames(stats) <- c(
      "Raw data",
      "Reliable CpGs",
      "Mapped CpGs",
      "CpGs representing SNPs",
      "CpGs without SNPs"
    )
    
    stats_data(stats)  # Now this is only called after filtering is done
    
    # Create Beta and M value tables
    beta_values_table <- getBeta(ratio_geno_Swan_NoSNP)
    mval_values_table <- getM(ratio_geno_Swan_NoSNP)
    
    # Assign them to reactive values for the tables
    beta_table(beta_values_table)
    mval_table(mval_values_table)
    
    # Enable download buttons once the tables are ready
    output$download_stats_ui <- renderUI({
      req(stats_data())  # Ensure stats table is available before showing the download button
      downloadButton("download_stats", "Download Stats")
    })
    
    output$download_beta <- renderUI({
      req(beta_table())  # Ensure the Beta table is ready
      downloadButton("download_beta", "Download Beta Values")
    })
    
    output$download_m <- renderUI({
      req(mval_table())  # Ensure the M values table is ready
      downloadButton("download_m", "Download M Values")
    })
  })
  
  # Show stats table only when it's ready
  output$filter_stats <- renderTable({
    req(stats_data())  # only display after it's created
    stats_data()
  }, rownames = TRUE)
  
  # Show Beta or M values table depending on the selection
  output$value_table_title <- renderText({
    req(input$value_type)  # Ensure selection exists
    input$value_type
  })
  
  output$value_table <- renderTable({
    req(input$value_type)
    if (input$value_type == "Beta values") {
      head(beta_table(), 10)  # Show first 10 rows of Beta table
    } else {
      head(mval_table(), 10)  # Show first 10 rows of M values table
    }
  }, rownames = TRUE)
  
  # Download handler for Stats table
  output$download_stats <- downloadHandler(
    filename = function() {
      paste0("filtering_stats_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(stats_data(), file, row.names = TRUE)
    }
  )
  
  # Download handler for Beta values
  output$download_beta <- downloadHandler(
    filename = function() {
      paste0("beta_values_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(beta_table(), file, row.names = TRUE)
    }
  )
  
  # Download handler for M values
  output$download_m <- downloadHandler(
    filename = function() {
      paste0("m_values_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(mval_table(), file, row.names = TRUE)
    }
  )
}

shinyApp(ui = ui_filter, server = server_filter)
