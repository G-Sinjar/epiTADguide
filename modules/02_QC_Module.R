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
    title = "QC Plots",
    
    # Sidebar with plot selection and download option
    sidebar = sidebar(
      width = 300,
      radioButtons(ns("qc_plot"), "Select QC Plot:",
                   choices = list(
                     "Plot 1: Channel Intensity" = "plot1",
                     "Plot 2: Sample Reliability" = "plot2",
                     "Plot 3: MDS plot (Sample PCA)" = "plot3",
                     "Plot 4: Beta distribution (Badewanne)" = "plot4",
                     "Plot 5: Bean Plot" = "plot5"
                   ),
                   selected = "plot1"),
      br(),
      downloadButton(ns("download_qc_report"), "Download QC Report"),
      helpText("Generating the report may take a moment. Once ready, a dialog will prompt you to choose the download location.")
    ),
    
    # Main content panel for plots
    uiOutput(ns("qc_plot_ui"))
  )
}




# ─────────────────────────────────────
# SERVER FUNCTION
# ─────────────────────────────────────
qcServer <- function(id, RGset, raw_normalised, targets) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Render the appropriate plot UI based on selection
    output$qc_plot_ui <- renderUI({
      switch(input$qc_plot,
             
             # Plot 1: Channel intensity
             "plot1" = tagList(
               h3("Channel Intensity"),
               helpText("This plot shows the intensity distribution of the red and green channels across all samples."),
               plotOutput(ns("plot1"), height = "100vh")
             ),
             
             # Plot 2: Sample reliability
             "plot2" = tagList(
               h3("Sample Reliability"),
               helpText("Visualizes detection p-values to assess the reliability of each sample."),
               plotOutput(ns("plot2"), height = "100vh")
             ),
             
             # Plot 3: MDS plot
             "plot3" = tagList(
               h3("MDS Plot (Sample PCA)"),
               helpText("Multidimensional scaling plot of raw beta values. Used to detect sample outliers and batch effects."),
               plotOutput(ns("plot3"), height = "100vh")
             ),
             
             # Plot 4: Beta distribution
             "plot4" = tagList(
               h3("Beta Distribution (Badewanne)"),
               helpText("Distribution of beta values for each sample, typically bimodal in shape."),
               plotOutput(ns("plot4"), height = "100vh")
             ),
             
             # Plot 5: Bean plot
             "plot5" = tagList(
               h3("Bean Plot"),
               helpText("Displays the distribution and density of methylation beta values per sample group."),
               plotOutput(ns("plot5"), height = "100vh")
             )
      )
    })
    
    # Plot 1: Channel intensities
    output$plot1 <- renderPlot({
      req(RGset())
      # Wrap the plotting code in withProgress
      withProgress(message = 'Generating Channel Intensity Plot...', value = 0, {
        incProgress(0.1, detail = "Processing data for plot")
        plot(density(as.vector(assay(RGset(), "Red"))), main = "Channel Intensities", lwd = 2)
        lines(density(as.vector(assay(RGset(), "Green"))), col = "green", lwd = 2)
        legend("topright", legend = c("Red Channel", "Green Channel"), col = c("black", "green"), lwd = 2)
        incProgress(1, detail = "Plot ready")
      })
    })
    
    # Plot 2: QC plot (Sample reliability)
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
    
    # Plot 3: MDS plot
    output$plot3 <- renderPlot({
      req(RGset(), targets())
      withProgress(message = 'Generating MDS Plot...', value = 0, {
        incProgress(0.1, detail = "Performing MDS analysis")
        mdsPlot(RGset(), sampNames = targets()$Array, sampGroups = targets()$Sample_Group,
                main = "Raw Beta MDS", legendNCol = 1, legendPos = "topright")
        incProgress(1, detail = "Plot ready")
      })
    })
    
    # Plot 4: Beta distribution
    output$plot4 <- renderPlot({
      req(RGset(), targets())
      withProgress(message = 'Generating Beta Distribution Plot...', value = 0, {
        incProgress(0.1, detail = "Calculating density")
        densityPlot(RGset(), sampGroups = targets()$Sample_Label,
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
        densityBeanPlot(RGset(), sampGroups = targets()$Sample_Label, sampNames = targets()$Sample_Name)
        incProgress(1, detail = "Plot ready")
      })
    })
    
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
  
  qcServer(
    "qc",
    RGset = reactive_RGset,
    raw_normalised = reactive_raw_normalised,
    targets = reactive_targets
  )
}

shinyApp(ui, server)'