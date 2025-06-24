# Inputs
#load -> RGset, raw_normalised, targets

'# Values to store the normalized data
SWAN_normalised <- NULL
Quantile_normalised <- NULL
funnorm_normalised <- NULL
noob_normalised <- NULL
noob.swan_normalised <- NULL

if (!is.null(RGset) && !is.null(raw_normalised)) {
  # Perform normalizations with tryCatch
  SWAN_normalised <<- tryCatch(preprocessSWAN(RGset, raw_normalised), error = function(e) {
    showNotification(paste("Error in SWAN normalization:", e$message), type = "error")
    NULL
  })
  Quantile_normalised <<- tryCatch(preprocessQuantile(RGset, fixOutliers = TRUE, removeBadSamples = TRUE, badSampleCutoff = 10.5, quantileNormalize = TRUE, stratified = TRUE, mergeManifest = FALSE, sex = NULL), error = function(e) {
    showNotification(paste("Error in Quantile normalization:", e$message), type = "error")
    NULL
  })
  funnorm_normalised <<- tryCatch(preprocessFunnorm(RGset), error = function(e) {
    showNotification(paste("Error in Funnorm normalization:", e$message), type = "error")
    NULL
  })
  noob_normalised <<- tryCatch(preprocessNoob(RGset), error = function(e) {
    showNotification(paste("Error in Noob normalization:", e$message), type = "error")
    NULL
  })
  noob.swan_normalised <<- tryCatch(preprocessSWAN(RGset, mSet = preprocessNoob(RGset), verbose = TRUE), error = function(e) {
    showNotification(paste("Error in Noob-SWAN normalization:", e$message), type = "error")
    NULL
  })
}

#add qc report button 

'

ui_norm <- page_navbar(
  title = "EPIC Array Pipeline",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  

  nav_panel("Normalisation",
            page_sidebar(
                sidebar = sidebar(
                width = 300,
                
                # Single radio input for plot selection
                radioButtons("norm_methode", "Select a normalisation methode:",
                             choices = list(
                               "SWAN" = "Swan",
                               "NOOB" = "Noob",
                               "NOOB.SWAN" = "Noob.Swan",
                               "FUNNORM" = "Funnorm",
                               "Quantile" = "Quantile"
                             ),
                             selected = "Swan")
              ),
              
              # Main panel changes based on selected plot
              uiOutput("beta_plot")
            )
  )
)




server_norm <- function(input, output, session) {
  # Access the globally assigned raw_normalised and targets
  output$beta_plot <- renderUI({
    tagList(
      h3("Raw Normalised"),
      helpText("Beta values before applying the selected normalisation method."),
      fluidRow(
        column(2),
        column(8, plotOutput("Raw", height = "300px", width = "100%")),
        column(2)
      ),
      switch(input$norm_methode,
             "Swan" = tagList(
               h3("SWAN (Subset-quantile Within Array Normalization)"),
               helpText("..."),
               fluidRow(
                 column(2),
                 column(8, plotOutput("Swan", height = "300px", width = "100%")),
                 column(2)
               )
             ),
             "Noob" = tagList(
               h3("Noob (Normal-exponential Out-of-Band Correction)"),
               helpText("..."),
               fluidRow(
                 column(2),
                 column(8, plotOutput("Noob", height = "300px", width = "100%")),
                 column(2)
               )
             ),
             "Noob.Swan" = tagList(
               h3("Noob then SWAN normalisation"),
               helpText("..."),
               fluidRow(
                 column(2),
                 column(8, plotOutput("Noob.Swan", height = "300px", width = "100%")),
                 column(2)
               )
             ),
             "Funnorm" = tagList(
               h3("Funnorm (Functional Normalization)"),
               helpText("..."),
               fluidRow(
                 column(2),
                 column(8, plotOutput("Funnorm", height = "300px", width = "100%")),
                 column(2)
               )
             ),
             "Quantile" = tagList(
               h3("Quantile normalisation"),
               helpText("..."),
               fluidRow(
                 column(2),
                 column(8, plotOutput("Quantile", height = "300px", width = "100%")),
                 column(2)
               )
             )
      )
    )
  })
  
  # Sample plots - Access the globally assigned normalized data and targets
  output$Raw <- renderPlot({
    req(raw_normalised)
    req(exists("targets", envir = .GlobalEnv)) # Check if 'targets' exists globally
    densityPlot(raw_normalised, sampGroups = targets$Sample_Label, main = "Beta Values distribution of raw data", pal = rainbow(length(unique(targets$Sample_Label))))
  })
  output$Swan <- renderPlot({
    req(SWAN_normalised)
    req(exists("targets", envir = .GlobalEnv))
    if (!is.null(SWAN_normalised)) {
      densityPlot(SWAN_normalised, sampGroups = targets$Sample_Label, main = "Beta values distribution after SWAN normalisation", pal = rainbow(length(unique(targets$Sample_Label))))
    } else {
      plot(1, type = "n", main = "SWAN Normalization Not Available")
    }
  })
  output$Noob <- renderPlot({
    req(noob_normalised)
    req(exists("targets", envir = .GlobalEnv))
    if (!is.null(noob_normalised)) {
      densityPlot(noob_normalised, sampGroups = targets$Sample_Label, main = "Beta values distribution after Noob", pal = rainbow(length(unique(targets$Sample_Label))))
    } else {
      plot(1, type = "n", main = "Noob Normalization Not Available")
    }
  })
  output$Noob.Swan <- renderPlot({
    req(noob.swan_normalised)
    req(exists("targets", envir = .GlobalEnv))
    if (!is.null(noob.swan_normalised)) {
      densityPlot(noob.swan_normalised, sampGroups = targets$Sample_Label, main = "Beta values distribution after Noob-SWAN", pal = rainbow(length(unique(targets$Sample_Label))))
    } else {
      plot(1, type = "n", main = "Noob-Swan Normalization Not Available")
    }
  })
  output$Funnorm <- renderPlot({
    req(funnorm_normalised)
    req(exists("targets", envir = .GlobalEnv))
    if (!is.null(funnorm_normalised)) {
      densityPlot(getBeta(funnorm_normalised), sampGroups = targets$Sample_Label, main = "Beta Values distribution after Funnorm normalisation", pal = rainbow(length(unique(targets$Sample_Label))))
    } else {
      plot(1, type = "n", main = "Funnorm Normalization Not Available")
    }
  })
  output$Quantile <- renderPlot({
    req(Quantile_normalised)
    req(exists("targets", envir = .GlobalEnv))
    if (!is.null(Quantile_normalised)) {
      densityPlot(getBeta(Quantile_normalised), sampGroups = targets$Sample_Label, main = "Beta Values distribution after Quantile normalisation", pal = rainbow(length(unique(targets$Sample_Label))))
    } else {
      plot(1, type = "n", main = "Quantile Normalization Not Available")
    }
  })
}

shinyApp(ui_norm, server_norm)
