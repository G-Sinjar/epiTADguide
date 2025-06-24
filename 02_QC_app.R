# 02_QC_app.R
# Author: Ghazal Sinjar
# Date: 30.05.2025
# Description: Shiny app to display quality control plots for EPIC array data.
# The app loads preprocessed data and displays different QC plots based on user selection.
# It also allows the user to download a full QC report as a PDF.

library(shiny)
library(minfi)
library(bslib)

# UI Definition
ui <- page_navbar(
  title = "EPIC Array Pipeline",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  
  nav_panel("QC Plots",
            page_sidebar(
              title = "QC Plots",
              
              # Sidebar with plot selection and download option
              sidebar = sidebar(
                width = 300,
                radioButtons("qc_plot", "Select QC Plot:",
                             choices = list(
                               "Plot 1: Channel Intensity" = "plot1",
                               "Plot 2: Sample Reliability" = "plot2",
                               "Plot 3: MDS plot (Sample PCA)" = "plot3",
                               "Plot 4: Beta distribution (Badewanne)" = "plot4",
                               "Plot 5: Bean Plot" = "plot5"
                             ),
                             selected = "plot1"),
                br(),
                downloadButton("download_qc_report", "Download QC Report"),
                helpText("Generating the report may take a moment. Once ready, a dialog will prompt you to choose the download location.")
              ),
              
              # Main content panel for plots
              uiOutput("qc_plot_ui")
            )
  )
)

# Server Logic
server <- function(input, output, session) {
  # Load preprocessed data
  data <- readRDS("preprocessed_data.rds")
  RGset <- data$RGset
  raw_normalised <- data$raw_normalised
  targets <- data$targets
  
  # Render the appropriate plot UI based on selection
  output$qc_plot_ui <- renderUI({
    switch(input$qc_plot,
           
           # Plot 1: Channel intensity
           "plot1" = tagList(
             h3("Channel Intensity"),
             helpText("This plot shows the intensity distribution of the red and green channels across all samples."),
             plotOutput("plot1", height = "100vh")
           ),
           
           # Plot 2: Sample reliability
           "plot2" = tagList(
             h3("Sample Reliability"),
             helpText("Visualizes detection p-values to assess the reliability of each sample."),
             plotOutput("plot2", height = "100vh")
           ),
           
           # Plot 3: MDS plot
           "plot3" = tagList(
             h3("MDS Plot (Sample PCA)"),
             helpText("Multidimensional scaling plot of raw beta values. Used to detect sample outliers and batch effects."),
             plotOutput("plot3", height = "100vh")
           ),
           
           # Plot 4: Beta distribution
           "plot4" = tagList(
             h3("Beta Distribution (Badewanne)"),
             helpText("Distribution of beta values for each sample, typically bimodal in shape."),
             plotOutput("plot4", height = "100vh")
           ),
           
           # Plot 5: Bean plot
           "plot5" = tagList(
             h3("Bean Plot"),
             helpText("Displays the distribution and density of methylation beta values per sample group."),
             plotOutput("plot5", height = "100vh")
           )
    )
  })
  
  # Plot 1: Channel intensities
  output$plot1 <- renderPlot({
    plot(density(as.vector(assay(RGset, "Red"))), main = "Channel Intensities", lwd = 2)
    lines(density(as.vector(assay(RGset, "Green"))), col = "green", lwd = 2)
    legend("topright", legend = c("Red Channel", "Green Channel"), col = c("black", "green"), lwd = 2)
  })
  
  # Plot 2: QC plot
  output$plot2 <- renderPlot({
    qc <- getQC(raw_normalised)
    raw_normalised <- addQC(raw_normalised, qc)
    plotQC(qc)
  })
  
  # Plot 3: MDS plot
  output$plot3 <- renderPlot({
    mdsPlot(RGset, sampNames = targets$Array, sampGroups = targets$Sample_Group ,
            main = "Raw Beta MDS", legendNCol = 1, legendPos = "topright")
  })
  
  # Plot 4: Beta distribution
  output$plot4 <- renderPlot({
    densityPlot(RGset, sampGroups = targets$Sample_Label, 
                main ="Beta values distribution of raw data", xlab = "Beta values")
  })
  
  # Plot 5: Bean plot
  output$plot5 <- renderPlot({
    par(mar = c(5, 6, 4, 2))
    densityBeanPlot(RGset, sampGroups = targets$Sample_Label, sampNames = targets$Sample_Name)
  })
  
  # Download Handler for QC report
  output$download_qc_report <- downloadHandler(
    filename = function() {
      paste0("qcReport_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      tryCatch({
        qcReport(
          RGset,
          sampNames = targets$Sample_Name,
          sampGroups = targets$Sample_Group0,
          pdf = file
        )
      }, error = function(e) {
        pdf(file)
        plot.new()
        text(0.5, 0.5, paste("Error generating QC report:", e$message), cex = 1.5)
        dev.off()
      })
    }
  )
}

# Run the Shiny App
shinyApp(ui, server)
