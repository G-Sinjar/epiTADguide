# libraries
#library(minfi)
#library("IlluminaHumanMethylationEPICv2manifest")
#library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)

# 1- load data
#RawDataDir <- "C:\\Users\\ghaza\\Documents\\ghazal\\Bioinformatik_Fächer\\Masterarbeit_Project\\Data\\CLDN10_RawData\\epic_data\\2024_071_ILL_METUVG_N_8\\METUVG"
#targets <- read.metharray.sheet(RawDataDir)
# Load methydata from idat files 
#RGset <- read.metharray.exp(targets = targets)
#manifest <- getManifest(RGset)
#manifest

# 2- create a data set on the locus level
#raw_normalised <- preprocessRaw(RGset)
#dim(raw_normalised)

# 3- QC
## raw density plot
## sample specific plot
## MDS plot for arrays
## beta plot
## bean plot
## note QC repot is saved in the results folder created in the raw data folder you p


#add this to first card
sample_names <- sampleNames(RGset)
num_samples <- length(sample_names)
print(num_samples)
print(sample_names)


#library(shiny)
#library(ggplot2)
#library(minfi)
#library(bslib)
#library(shinyjs)

ui <- page_navbar(
  title = "EPIC Array Pipeline",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  
  # -- Summary tab with sidebar input
  nav_panel("Load Raw Data",
            layout_sidebar(
              sidebar = sidebar(
                textInput("dir_path", "Enter Path of Raw Data folder:", placeholder = "e.g. C:/path/to/data"),
                actionButton("load_btn", "Load Data")
              ),
              layout_columns(
                card(
                  card_header("SampleSheet and IDAT Load Status"),
                  textOutput("load_status")
                ),
                card(
                  card_header("Infos about the Array"),
                  verbatimTextOutput("manifest_output")
                )
              )
            )
  ),
  
  # -- QC Plots tab with sidebar pill navigation
  nav_panel("QC Plots",
            page_sidebar(
              title = "QC Plots",
              sidebar = sidebar(
                width = 300,
                
                # Single radio input for plot selection
                radioButtons("qc_plot", "Select QC Plot:",
                             choices = list(
                               "Plot 1: Signal Intensity" = "plot1",
                               "Plot 2: Sample Reliability" = "plot2",
                               "Plot 3: MDS plot (Sample PCA)" = "plot3",
                               "Plot 4: Beta plot (Badenwanne)" = "plot4",
                               "Plot 5: Bean Plot" = "plot5"
                             ),
                             selected = "plot1")
              ),
              
              # Main panel changes based on selected plot
              uiOutput("qc_plot_ui")
            )
  )
)

server <- function(input, output, session) {
  observeEvent(input$load_btn, {
    req(input$dir_path)
    
    try({
      RawDataDir <- input$dir_path
      
      message_output <- capture.output({
        targets <- read.metharray.sheet(RawDataDir)
      }, type = "message")
      
      output_lines <- capture.output({
        targets <- read.metharray.sheet(RawDataDir)
      }, type = "output")
      
      all_output <- c(message_output, output_lines)
      RGset <- read.metharray.exp(targets = targets)
      manifest <- getManifest(RGset)
      
      output$load_status <- renderText({
        paste(all_output, collapse = "\n")
      })
      
      output$manifest_output <- renderPrint({
        manifest
      })
    }, silent = TRUE)
  })
  
  
  output$qc_plot_ui <- renderUI({
    switch(input$qc_plot,
           "plot1" = tagList(
             h3("Signal Intensity"),
             helpText("..."),
             plotOutput("plot1", height = "100vh", width = "100%")
           ),
           "plot2" = tagList(
             h3("Detection P-Values"),
             helpText("..."),
             plotOutput("plot2", height = "100vh", width = "100%")
           ),
           "plot3" = tagList(
             h3("Control Probes Heatmap"),
             helpText("..."),
             plotOutput("plot3", height = "100vh", width = "100%")
           ),
           "plot4" = tagList(
             h3("Sample Clustering"),
             helpText("..."),
             plotOutput("plot4", height = "100vh", width = "100%")
           ),
           "plot5" = tagList(
             h3("Beta Value Density"),
             helpText("..."),
             plotOutput("plot5", height = "100vh", width = "100%")
           )
    )
  })
  
  # Sample plots
  output$plot1 <- renderPlot({ plot(density(as.vector(assay(RGset, "Red"))), main = "Channel Intensities", lwd = 2) 
    lines(density(as.vector(assay(RGset, "Green"))), col = "green", lwd = 2)
    legend("topright", legend = c("Red Channel", "Green Channel"), col = c("black", "green"), lwd = 2) })
  output$plot2 <- renderPlot({ plotQC(qc) })
  output$plot3 <- renderPlot({ mdsPlot(RGset, sampNames = targets$Array, sampGroups = targets$Sample_Group ,main = "Raw Beta MDS", legendNCol = 1, legendPos = "topright")
 })
  output$plot4 <- renderPlot({ densityPlot(RGset, sampGroups = targets$Sample_Label, main ="Beta values distribution of raw data", xlab = "Beta values")
 })
  output$plot5 <- renderPlot({ par(mar = c(5, 12, 4, 2))
    densityBeanPlot(RGset, sampGroups = targets$Sample_Label, sampNames = targets$Sample_Name) })
}


shinyApp(ui, server)
