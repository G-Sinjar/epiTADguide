

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
#sample_names <- sampleNames(RGset)
#num_samples <- length(sample_names)
#print(num_samples)
#print(sample_names)

#libraries
library(minfi)
library("IlluminaHumanMethylationEPICv2manifest")
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(shiny)
library(bslib)
library(shinyjs)
#library(ggplot2)

ui <- page_navbar(
  title = "EPIC Array Pipeline",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  
  # -- Summary tab with sidebar input
  nav_panel("Load Raw Data",
            layout_sidebar(
              sidebar = sidebar(
                textInput("dir_path", "Enter Path of Raw Data folder:", placeholder = r"(e.g. C:\path\to\data)"),
                helpText("The specified directory must contain a completed Sample Sheet as well as a subfolder named after the slide ID, which includes all corresponding IDAT files."),
                actionButton("load_btn", "Load Data")
              ),
              layout_columns(
                card(
                  card_header("SampleSheet and IDAT Load Status"),
                  uiOutput("load_status")
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
                               "Plot 1: Channel Intensity" = "plot1",
                               "Plot 2: Sample Reliability" = "plot2",
                               "Plot 3: MDS plot (Sample PCA)" = "plot3",
                               "Plot 4: Beta distribution (Badenwanne)" = "plot4",
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
  # Initial status message
  status <- reactiveVal("ℹ️ Waiting to enter the path to data from sequencer.")
  
  output$load_status <- renderUI({
    HTML(gsub("\n", "<br>", status()))
  })
  
  observeEvent(input$load_btn, {
    req(input$dir_path)

    # Reset status
    status("🔄 Starting data load...")

    output$load_status <- renderUI({
      HTML(gsub("\n", "<br>", status()))
    })

    # Normalize path
    RawDataDir <- gsub("\\\\", "/", input$dir_path)

    # Step 1: Read Sample Sheet
    status(paste(status(), "\nStep 1: Reading Sample Sheet..."))
    targets <- tryCatch({
      read.metharray.sheet(RawDataDir)
    }, error = function(e) {
      status(paste(status(), "\n❌ Error in Step 1:", e$message))
      return(NULL)
    })
    if (is.null(targets)) return()
    status(paste(status(), "\n✅ Sample Sheet loaded."))

    # Step 2: Load IDAT files
    status(paste(status(), "\nStep 2: Loading IDAT files..."))
    RGset <<- tryCatch({
      read.metharray.exp(targets = targets)
    }, error = function(e) {
      status(paste(status(), "\n❌ Error in Step 2:", e$message))
      return(NULL)
    })
    if (is.null(RGset)) return()
    status(paste(status(), "\n✅ IDAT files loaded."))

    # Step 3: Get Manifest
    status(paste(status(), "\nStep 3: Getting array manifest info..."))
    manifest <- tryCatch({
      getManifest(RGset)
    }, error = function(e) {
      status(paste(status(), "\n❌ Error in Step 3:", e$message))
      return(NULL)
    })
    if (is.null(manifest)) return()
    output$manifest_output <- renderPrint({ manifest })
    status(paste(status(), "\n✅ Manifest info retrieved."))

    # Step 4: Preprocess
    status(paste(status(), "\nStep 4: Converting to locus-level data..."))
    raw_normalised <<- tryCatch({
      preprocessRaw(RGset)
    }, error = function(e) {
      status(paste(status(), "\n❌ Error in Step 4:", e$message))
      return(NULL)
    })
    if (is.null(raw_normalised)) return()
    dims <- dim(raw_normalised)
    status(paste(status(), "\n✅ Locus-level data ready. Dimensions: ", paste(dims, collapse = " x ")))

    # Step 5: Sample info
    status(paste(status(), "\nStep 5: Extracting sample info and CpG count..."))
    sample_names <- sampleNames(RGset)
    num_samples <- length(sample_names)
    num_cpgs <- nrow(raw_normalised)
    status(paste(
      status(),
      "\n✅ Sample and CpG info retrieved:",
      "\n  - Number of CpGs: ", num_cpgs,
      "\n  - Number of Samples: ", num_samples,
      "\n  - Sample Names: \n", paste(sample_names, collapse = "\n")
    ))
  })


  
  output$qc_plot_ui <- renderUI({
    switch(input$qc_plot,
           
           "plot1" = tagList(
             h3("Channel Intensity"),
             helpText("Density plot of signal intensities from Red and Green channels across arrays. Helps assess signal quality."),
             plotOutput("plot1", height = "100vh")
           ),
           
           "plot2" = tagList(
             h3("Sample Reliability"),
             helpText("Shows quality metrics per sample (e.g., signal noise or background), helping detect unreliable arrays."),
             plotOutput("plot2", height = "100vh")
           ),
           
           "plot3" = tagList(
             h3("MDS Plot (Sample PCA)"),
             helpText("Visualizes clustering of samples based on methylation profiles, using multidimensional scaling."),
             plotOutput("plot3", height = "100vh"),
           ),
           
           "plot4" = tagList(
             h3("Beta Distribution (Badewanne)"),
             helpText("Shows the distribution of beta values across all probes and samples. Useful to detect technical issues."),
             plotOutput("plot4", height = "100vh"),
           ),
           
           "plot5" = tagList(
             h3("Bean Plot"),
             helpText("Displays beta-value densities per sample using bean plots. Useful to detect batch effects or outliers."),
             plotOutput("plot5", height = "100vh"),
           )
    )
  })
}  




shinyApp(ui, server)
