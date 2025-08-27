# 02_QC_Module.R
# Author: Ghazal Sinjar
# Date: 30.05.2025
# Description: Shiny module to display quality control plots for EPIC array data.
# The module loads preprocessed data and displays different QC plots based on user selection.
# It also allows the user to download a full QC report as a PDF.


# ─────────────────────────────────────
# User Interface (UI) FUNCTION
# ─────────────────────────────────────
qcUI <- function(id) {
  ns <- NS(id)
  
  page_sidebar(
    # Sidebar with plot selection and download option
    sidebar = sidebar(
      width = 300,
      radioButtons(ns("qc_plot"), "Select QC Plot:",
                   choices = list(
                     "Plot 1: Channel Intensity" = "plot1",
                     "Plot 2: Medien Intensity (Sample-specific)" = "plot2",
                     "Plot 3: Multi-dimensional scaling (MDS) plot" = "plot3",
                     "Plot 4: Beta value densities of the samples" = "plot4",
                     "Plot 5: Bean plot of Beta values densities" = "plot5"
                   ),
                   selected = "plot1"),
      radioButtons(ns("plot_format"), "Select Format:", choices = c("PDF" = "pdf", "PNG" = "png"), inline = TRUE),
      downloadButton(ns("download_current_plot"), "Download Current Plot"),
      hr(),
      downloadButton(ns("download_qc_report"), "Download QC Report"),
      helpText("This button generates a comprehensive PDF report that compiles multiple standard QC plots into a single document. Note that some plots—such as control probe plots—are included only in the report and are not displayed interactively within the app. We recommend reviewing the report for a complete quality assessment.")
    ),
    
    # Main content panel for plots
    uiOutput(ns("qc_plot_ui"))
  )
}




# ─────────────────────────────────────
# SERVER FUNCTION
# ─────────────────────────────────────
qcServer <- function(id, RGset, raw_normalised, targets, project_output_dir) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Render the appropriate plot UI based on selection
    output$qc_plot_ui <- renderUI({
      switch(input$qc_plot,
             
             # Plot 1: Channel intensity
             "plot1" = tagList(
               h3("Channel Intensity"),
               helpText("The Channel Intensity plot displays the density distribution of raw red and green channel fluorescence intensities across all samples. This visualization is crucial for rapidly identifying samples with aberrant overall signal levels, which could indicate issues such as poor hybridization, insufficient DNA input, or errors during array scanning."),
               plotOutput(ns("plot1"), height = "100vh")
             ),
             
             # Plot 2: sample-specific QC
             "plot2" = tagList(
               h3("Medien Intensity (Sample-specific)"),
               helpText("	This is a simple quality control plot that uses the log median intensity in both the methylated (M) and unmethylated (U) channels. When plotting these two medians against each other, it has been observed that good samples cluster together, while failed samples tend to separate and have lower median intensities."),
               plotOutput(ns("plot2"), height = "100vh")
             ),
             
             # Plot 3: MDS plot
             "plot3" = tagList(
               h3("Multi-dimensional scaling (MDS) plot"),
               helpText("Multi-dimensional scaling plot gives an overview of similarities and differences between samples using Euclidean distance."),
               plotOutput(ns("plot3"), height = "100vh")
             ),
             
             # Plot 4: Beta value densities of the samples
             "plot4" = tagList(
               h3("Beta values densities of the samples"),
               helpText("A Distribution of beta values for each sample. This plot helps assessing the overall methylation profile of a sample and identify potential problems.\nThis distribution is expected to be bimodal with the 2 peaks (around 0 and 1 beta values) representing methylated and unmethylated signals. Any center peaks should be further investigated (e.g. in Plot5) for problems."),
               plotOutput(ns("plot4"), height = "100vh")
             ),
             
             # Plot 5: Bean plot of Beta values densities
             "plot5" = tagList(
               h3("Bean plot of Beta values densities"),
               helpText("Displays the density of beta values for each individual sample, with samples colored according to their respective experimental groups.\nthis granular view allows for the precise identification of individual samples exhibiting unusual or defective beta value distributions. Such deviations can indicate quality control issues at the single-sample level, necessitating further investigation or potential exclusion of affected samples from downstream analyses."),
               plotOutput(ns("plot5"), height = "100vh")
             )
      )
    })
    
    #-----------------------------------------
    # Plot 1: Channel intensities
    output$plot1 <- renderPlot({
      req(RGset())
      # Wrap the plotting code in withProgress
      withProgress(message = 'Generating Channel Intensity Plot...', value = 0, {
        incProgress(0.1, detail = "Processing data for plot")
        plot(density(as.vector(assay(RGset(), "Red"))), main = "Channel Intensities", col= "red" ,lwd = 2)
        lines(density(as.vector(assay(RGset(), "Green"))), col = "green", lwd = 2)
        legend("topright", legend = c("Red Channel", "Green Channel"), col = c("red", "green"), lwd = 2)
        incProgress(1, detail = "Plot ready")
      })
    })
    
    # Plot 2: sample-specific QC: median-intensity QC
    output$plot2 <- renderPlot({
      req(raw_normalised())
      withProgress(message = 'Generating Sample Reliability Plot...', value = 0, {
        incProgress(0.1, detail = "Calculating QC metrics")
        qc <- getQC(raw_normalised())
        raw_normalised_qc <- addQC(raw_normalised(), qc)
        incProgress(0.5, detail = "Plotting QC results")
        plotQC(qc)
        incProgress(1, detail = "Plot ready")
      })
    })
    
    # Plot 3: Multi-dimensional scaling (MDS) plot
    output$plot3 <- renderPlot({
      req(RGset(), targets())
      withProgress(message = 'Generating MDS Plot...', value = 0, {
        incProgress(0.1, detail = "Performing MDS analysis")
        mdsPlot(RGset(), sampNames = targets()$Sample_Name, sampGroups = targets()$Sample_Group,
                main = "Raw Beta MDS", legendNCol = 1, legendPos = "topright")
        incProgress(1, detail = "Plot ready")
      })
    })
    
    # Plot 4: Beta distribution
    output$plot4 <- renderPlot({
      req(RGset(), targets())
      withProgress(message = 'Generating Beta Distribution Plot...', value = 0, {
        incProgress(0.1, detail = "Calculating density")
        densityPlot(RGset(), sampGroups = targets()$Sample_Group,
                    main = "Beta values distribution of raw data", xlab = "Beta values")
        incProgress(1, detail = "Plot ready")
      })
    })
    
    # Plot 5: Bean plot
    output$plot5 <- renderPlot({
      req(RGset(), targets())
      withProgress(message = 'Generating Bean Plot...', value = 0, {
        incProgress(0.1, detail = "Preparing bean plot data")
        par(mar = c(5, 6, 4, 2))
        densityBeanPlot(RGset(), sampGroups = targets()$Sample_Group, sampNames = targets()$Sample_Name)
        incProgress(1, detail = "Plot ready")
      })
    })
    
    #---------------------------------------------
    output$download_current_plot <- downloadHandler(
      filename = function() {
        
        paste0("QC_", input$qc_plot, "_", Sys.Date(), ".", input$plot_format)    
        },
      content = function(file) {
        req(input$qc_plot)
        
        plot_fun <- switch(input$qc_plot,
                           "plot1" = function() {
                             plot(density(as.vector(assay(RGset(), "Red"))), main = "Channel Intensities", lwd = 2)
                             lines(density(as.vector(assay(RGset(), "Green"))), col = "green", lwd = 2)
                             legend("topright", legend = c("Red Channel", "Green Channel"), col = c("black", "green"), lwd = 2)
                           },
                           "plot2" = function() {
                             qc <- getQC(raw_normalised())
                             plotQC(qc)
                           },
                           "plot3" = function() {
                             mdsPlot(RGset(), sampNames = targets()$Sample_Name, sampGroups = targets()$Sample_Group,
                                     main = "Raw Beta MDS", legendNCol = 1, legendPos = "topright")
                           },
                           "plot4" = function() {
                             densityPlot(RGset(), sampGroups = targets()$Sample_Group,
                                         main = "Beta values distribution of raw data", xlab = "Beta values")
                           },
                           "plot5" = function() {
                             par(mar = c(5, 6, 4, 2))
                             densityBeanPlot(RGset(), sampGroups = targets()$Sample_Group, sampNames = targets()$Sample_Name)
                           })
        
        if (input$plot_format == "pdf") {
          pdf(file,width = 12, height = 7)
          plot_fun()
          dev.off()
        } else if (input$plot_format == "png") {
          png(file, width = 1700, height = 1000, res = 150)
          plot_fun()
          dev.off()
        }
      }
    )
    #-------------------------------------------
    # Download Handler for QC report
    output$download_qc_report <- downloadHandler(
      filename = function() {
        paste0("qcReport_", Sys.Date(), ".pdf")
      },
      content = function(file) {
        req(RGset(), targets())
        withProgress(message = 'Generating QC Report PDF...', value = 0, {
          incProgress(0.1, detail = "Starting report generation")
          tryCatch({
            qcReport(
              RGset(),
              sampNames = targets()$Sample_Name,
              sampGroups = targets()$Sample_Group0,
              pdf = file
            )
            incProgress(1, detail = "Report generated successfully")
          }, error = function(e) {
            pdf(file)
            plot.new()
            text(0.5, 0.5, paste("Error generating QC report:", e$message), cex = 1.5)
            dev.off()
            incProgress(0, detail = paste("Error:", e$message))
            showNotification(paste("Error generating QC report:", e$message), type = "error", duration = 10)
          })
        })
      }
    )
  })
}

'library(shiny)
library(bslib)
library(minfi) 

rds_path <- "../intermediate_data/preprocessed_data.rds"
initial_path <- "./test_project_dir"
if (!dir.exists(initial_path)) {
  dir.create(initial_path, recursive = TRUE)
  message("Created dummy project directory: ", initial_path)
}
# Check for the RDS file BEFORE launching the app
if (!file.exists(rds_path)) {
  stop("Error: File not found: preprocessed_data.rds. Please run the data loading module first.")
}

# Load the data once
data <- readRDS(rds_path)

ui <- page_navbar(
  title = "EPIC Array Pipeline",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  nav_panel("QC Plots", qcUI("qc")) # Assuming qcUI is defined somewhere accessible
)

server <- function(input, output, session) {
  reactive_RGset <- reactiveVal(data$RGset)
  reactive_raw_normalised <- reactiveVal(data$raw_normalised)
  reactive_targets <- reactiveVal(data$targets)
  simulated_project_path <- reactiveVal(initial_path)
  qcServer(
    "qc",
    RGset = reactive_RGset,
    raw_normalised = reactive_raw_normalised,
    targets = reactive_targets,
    project_output_dir =simulated_project_path
  )
}

shinyApp(ui, server)'