
set.seed(123)
SWAN_normalised <- preprocessSWAN(RGset, raw_normalised)

Quantile_normalised <- preprocessQuantile(RGset, fixOutliers = TRUE,removeBadSamples = TRUE, badSampleCutoff = 10.5, quantileNormalize = TRUE, stratified = TRUE, mergeManifest = FALSE, sex = NULL)

funnorm_normalised <- preprocessFunnorm(RGset)

noob_normalised <- preprocessNoob(RGset)

set.seed(123)
noob.swan_normalised <- preprocessSWAN(RGset, mSet = noob_normalised, verbose = TRUE)

#add qc report button 



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
  output$beta_plot <- renderUI({
    tagList(
      h3("Raw Normalised"),
      helpText("Beta values before applying the selected normalisation method."),
      fluidRow(
        column(2),  # left spacer (15%)
        column(8, plotOutput("Raw", height = "300px", width = "100%")),  # center 70%
        column(2)   # right spacer (15%)
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
  
  
  
  # Sample plots
  output$Raw <- renderPlot({ densityPlot(raw_normalised,sampGroups = targets$Sample_Label,main = "Beta Values distribution of raw data", pal = rainbow(length(unique(targets$Sample_Label))))})
  output$Swan <- renderPlot({ densityPlot(SWAN_normalised, sampGroups = targets$Sample_Label , main = "Beta values distribution after SWAN normalisation", pal = rainbow(length(unique(targets$Sample_Label))))})
  output$Noob <- renderPlot({ densityPlot(noob_normalised, sampGroups = targets$Sample_Label ,main = "Beta values distribution after Noob", pal = rainbow(length(unique(targets$Sample_Label)))) })
  output$Noob.Swan <- renderPlot({ densityPlot(noob.swan_normalised,sampGroups = targets$Sample_Label, main = "Beta values distribution after Noob-SWAN", pal = rainbow(length(unique(targets$Sample_Label))))})
  output$Funnorm <- renderPlot({densityPlot(getBeta(funnorm_normalised),sampGroups = targets$Sample_Label,main = "Beta Values distribution after Funnorm normalisation", pal = rainbow(length(unique(targets$Sample_Label)))) })
  output$Quantile <- renderPlot({ densityPlot(getBeta(Quantile_normalised),sampGroups = targets$Sample_Label,main = "Beta Values distribution after Quantile normalisation", pal = rainbow(length(unique(targets$Sample_Label))))})
}

shinyApp(ui_norm, server_norm)
