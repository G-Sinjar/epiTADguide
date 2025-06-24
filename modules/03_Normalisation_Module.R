# 03_Normalisation_Module.R
# Author: Ghazal Sinjar
# Date: 02.06.2025
# Description:
# This Shiny module provides an interactive UI and server logic
# for visualising the effect of various normalisation methods
# applied to Illumina EPIC methylation array data.

# ─────────────────────────────────────
# User Interface (UI) FUNCTION
# ─────────────────────────────────────
norm_ui <- function(id) {
  ns <- NS(id)
  
  # Layout with sidebar to choose normalization method and main panel for plots
  page_sidebar(
    sidebar = sidebar(
      width = 300,
      
      # Radio buttons for selecting normalization method
      radioButtons(ns("norm_methode"), "Select a normalisation methode:",
                   choices = list(
                     "SWAN" = "Swan",
                     "NOOB" = "Noob",
                     "NOOB.SWAN" = "Noob.Swan",
                     "FUNNORM" = "Funnorm",
                     "Quantile" = "Quantile"
                   ),
                   selected = "Swan"
      )
    ),
    # Placeholder for dynamic UI based on method selected
    uiOutput(ns("beta_plot_container"))
  )
}



# ─────────────────────────────────────
# SERVER FUNCTION (with progress and robust error handling)
# ─────────────────────────────────────

norm_server <- function(id, RGset, raw_normalised, targets) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Reactive values to hold the results of NORMALIZATION methods
    SWAN_normalised_res <- reactiveVal(NULL)
    Quantile_normalised_res <- reactiveVal(NULL)
    funnorm_normalised_res <- reactiveVal(NULL)
    noob_normalised_res <- reactiveVal(NULL)
    noob_swan_normalised_res <- reactiveVal(NULL)
    
    # Use an observeEvent for the RGset to trigger normalizations
    observeEvent(RGset(), {
      req(RGset(), raw_normalised(), targets())
      
      # Reset results when new data is loaded
      SWAN_normalised_res(NULL)
      Quantile_normalised_res(NULL)
      funnorm_normalised_res(NULL)
      noob_normalised_res(NULL)
      noob_swan_normalised_res(NULL)
      
      # --- Start withProgress for the entire normalization process ---
      withProgress(message = 'Performing Normalization...', value = 0, {
        
        # Increment for each normalization method
        num_methods <- 5 
        inc_amount <- 1 / num_methods
        
# SWAN Normalization
        incProgress(0, detail = "Starting SWAN normalization...")
        # Set seed for reproducibility before calling preprocessSWAN
        set.seed(123)
        SWAN_normalised_res(
          tryCatch({
            minfi::preprocessSWAN(RGset(), mSet = raw_normalised())
          }, error = function(e) {
            showNotification(paste("Error in SWAN normalization:", e$message), type = "error", duration = 10)
            message("SWAN normalization failed: ", e$message) # Log to console
            return(NULL) 
          })
        )
        incProgress(inc_amount, detail = "SWAN normalization complete.")
        
        
        # Quantile Normalization
        incProgress(0, detail = "Starting Quantile normalization...")
        Quantile_normalised_res(
          tryCatch({
            minfi::preprocessQuantile(RGset(), fixOutliers = TRUE, removeBadSamples = FALSE,
                                      badSampleCutoff = 10.5, quantileNormalize = TRUE,
                                      stratified = TRUE, mergeManifest = FALSE, sex = NULL)
          }, error = function(e) {
            showNotification(paste("Error in Quantile normalization:", e$message), type = "error", duration = 10)
            message("Quantile normalization failed: ", e$message) # Log to console
            return(NULL) 
          })
        )
        incProgress(inc_amount, detail = "Quantile normalization complete.")
        
        
        # FUNNORM Normalization
        incProgress(0, detail = "Starting Funnorm normalization...")
        funnorm_normalised_res(
          tryCatch({
            minfi::preprocessFunnorm(RGset())
          }, error = function(e) {
            showNotification(paste("Error in Funnorm normalization:", e$message), type = "error", duration = 10)
            message("Funnorm normalization failed: ", e$message) # Log to console
            return(NULL) # Return NULL if an error occurs
          })
        )
        incProgress(inc_amount, detail = "Funnorm normalization complete.")
        
        
        # NOOB Normalization
        incProgress(0, detail = "Starting Noob normalization...")
        noob_normalised_res(
          tryCatch({
            minfi::preprocessNoob(RGset())
          }, error = function(e) {
            showNotification(paste("Error in Noob normalization:", e$message), type = "error", duration = 10)
            message("Noob normalization failed: ", e$message) # Log to console
            return(NULL) # Return NULL if an error occurs
          })
        )
        incProgress(inc_amount, detail = "Noob normalization complete.")
        
        
        # NOOB.SWAN Normalization (depends on noob_normalised_res)
        incProgress(0, detail = "Starting Noob.Swan normalization...")
        # Check if Noob result is available and not NULL due to previous error
        if (is.null(noob_normalised_res())) {
          showNotification("Noob result is missing or failed, cannot perform Noob.Swan.", type = "warning", duration = 10)
          noob_swan_normalised_res(NULL) # Explicitly set to NULL if dependency is missing
        } else {
          set.seed(123)
          noob_swan_normalised_res(
            tryCatch({
              minfi::preprocessSWAN(RGset(), mSet = noob_normalised_res(), verbose = TRUE)
            }, error = function(e) {
              showNotification(paste("Error in Noob.Swan normalization:", e$message), type = "error", duration = 10)
              message("Noob.Swan normalization failed: ", e$message) # Log to console
              return(NULL) # Return NULL if an error occurs
            })
          )
        }
        incProgress(inc_amount, detail = "Noob.Swan normalization complete.")
        
        
        # Final progress update
        incProgress(0, detail = "All normalizations finished!")
        
      }) # --- End withProgress block ---
      
      # Save all normalized datasets only if all were successful
      # We'll check for NULL results before saving
      if (!is.null(SWAN_normalised_res()) && !is.null(Quantile_normalised_res()) &&
          !is.null(funnorm_normalised_res()) && !is.null(noob_normalised_res()) &&
          !is.null(noob_swan_normalised_res())) {
        normalised_list <- list(
          SWAN = SWAN_normalised_res(),
          Quantile = Quantile_normalised_res(),
          Funnorm = funnorm_normalised_res(),
          Noob = noob_normalised_res(),
          Noob_Swan = noob_swan_normalised_res()
        )
        saveRDS(normalised_list, "./intermediate_data/normalised_all_methods.rds")
        message("✔️ Saved all normalised datasets to ./intermediate_data/normalised_all_methods.rds")
        showNotification("All normalisation methods completed and results saved!", type = "message", duration = 5)
      } else {
        showNotification("Some normalisation methods failed. Results not saved.", type = "warning", duration = 10)
      }
      
      
    }, ignoreNULL = TRUE, ignoreInit = FALSE, once = TRUE)
    
    
    # --- Dynamic UI for Normalisation Plot Display ---
    output$beta_plot_container <- renderUI({
      tagList(
        h3("Raw Normalised"),
        helpText("Beta values before applying the selected normalisation method."),
        fluidRow(column(2), column(8, plotOutput(ns("Raw"), height = "300px")), column(2)),
        
        # Always render plot containers for selected method
        switch(input$norm_methode,
               "Swan" = tagList(h3("SWAN"),
                                helpText("This method adjusts probe intensities within each array using quantile normalization."),
                                fluidRow(column(2), column(8, plotOutput(ns("Swan"))), column(2))),
               "Noob" = tagList(h3("Noob"),
                                helpText("Background correction using Noob method."),
                                fluidRow(column(2), column(8, plotOutput(ns("Noob"))), column(2))),
               "Noob.Swan" = tagList(h3("Noob then SWAN"),
                                     helpText("Background correction with Noob, followed by SWAN normalization."),
                                     fluidRow(column(2), column(8, plotOutput(ns("Noob.Swan"))), column(2))),
               "Funnorm" = tagList(h3("Funnorm"),
                                   helpText("Functional normalization using control probes."),
                                   fluidRow(column(2), column(8, plotOutput(ns("Funnorm"))), column(2))),
               "Quantile" = tagList(h3("Quantile"),
                                    helpText("Quantile normalization across samples."),
                                    fluidRow(column(2), column(8, plotOutput(ns("Quantile"))), column(2)))
        )
      )
    })
    
    # --- Density Plots (each wrapped in withProgress) ---
    output$Raw <- renderPlot({
      req(raw_normalised(), targets())
      withProgress(message = 'Generating Raw Data Density Plot...', value = 0, {
        incProgress(0.5, detail = "Plotting raw data")
        densityPlot(raw_normalised(), sampGroups = targets()$Sample_Label,
                    main = "Beta Values distribution of raw data",
                    pal = rainbow(length(unique(targets()$Sample_Label))))
        incProgress(1, detail = "Plot ready")
      })
    })
    
    output$Swan <- renderPlot({
      req(SWAN_normalised_res(), targets()) # Now requires SWAN_normalised_res to be non-NULL
      withProgress(message = 'Generating SWAN Density Plot...', value = 0, {
        incProgress(0.5, detail = "Plotting SWAN normalized data")
        densityPlot(SWAN_normalised_res(), sampGroups = targets()$Sample_Label,
                    main = "Beta values distribution after SWAN",
                    pal = rainbow(length(unique(targets()$Sample_Label))))
        incProgress(1, detail = "Plot ready")
      })
    })
    
    output$Noob <- renderPlot({
      req(noob_normalised_res(), targets()) # Now requires noob_normalised_res to be non-NULL
      withProgress(message = 'Generating Noob Density Plot...', value = 0, {
        incProgress(0.5, detail = "Plotting Noob normalized data")
        densityPlot(noob_normalised_res(), sampGroups = targets()$Sample_Label,
                    main = "Beta values distribution after Noob",
                    pal = rainbow(length(unique(targets()$Sample_Label))))
        incProgress(1, detail = "Plot ready")
      })
    })
    
    output$Noob.Swan <- renderPlot({
      req(noob_swan_normalised_res(), targets()) # Now requires noob_swan_normalised_res to be non-NULL
      withProgress(message = 'Generating Noob.Swan Density Plot...', value = 0, {
        incProgress(0.5, detail = "Plotting Noob.Swan normalized data")
        densityPlot(noob_swan_normalised_res(), sampGroups = targets()$Sample_Label,
                    main = "Beta values distribution after Noob-SWAN",
                    pal = rainbow(length(unique(targets()$Sample_Label))))
        incProgress(1, detail = "Plot ready")
      })
    })
    
    output$Funnorm <- renderPlot({
      req(funnorm_normalised_res(), targets()) # Now requires funnorm_normalised_res to be non-NULL
      withProgress(message = 'Generating Funnorm Density Plot...', value = 0, {
        incProgress(0.5, detail = "Plotting Funnorm normalized data")
        densityPlot(getBeta(funnorm_normalised_res()), sampGroups = targets()$Sample_Label,
                    main = "Beta values distribution after Funnorm",
                    pal = rainbow(length(unique(targets()$Sample_Label))))
        incProgress(1, detail = "Plot ready")
      })
    })
    
    output$Quantile <- renderPlot({
      req(Quantile_normalised_res(), targets()) 
      withProgress(message = 'Generating Quantile Density Plot...', value = 0, {
        incProgress(0.5, detail = "Plotting Quantile normalized data")
        densityPlot(getBeta(Quantile_normalised_res()), sampGroups = targets()$Sample_Label,
                    main = "Beta values distribution after Quantile",
                    pal = rainbow(length(unique(targets()$Sample_Label))))
        incProgress(1, detail = "Plot ready")
      })
    })
    
    # Return a list of reactive values for downstream modules
    return(list(
      SWAN = SWAN_normalised_res,
      Quantile = Quantile_normalised_res,
      Funnorm = funnorm_normalised_res,
      Noob = noob_normalised_res,
      Noob_Swan = noob_swan_normalised_res
    ))
  })
}


# 1) Libraries
library(shiny)  
library(bslib)  
library(minfi)  
library(IlluminaHumanMethylationEPICv2manifest)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)

#source("./03_Normalisation_Module.R")
path <- "C:/Users/ghaza/Documents/ghazal/Bioinformatik_Fächer/Masterarbeit_Project/Scripts/R_Scripts/intermediate_data/preprocessed_data.rds"
file.exists("./intermediate_data/preprocessed_data.rds")


preprocessed <- readRDS(path)
RGset <- preprocessed$RGset
targets <- preprocessed$targets
raw_normalised <- preprocessed$raw_normalised

# 3) Build the top‐level UI with three tabs (Load, QC, Normalization)
ui <- page_navbar(
  title = "EPIC Array Pipeline",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  
  nav_panel("Normalisation",
            norm_ui("norm")
  )
)

# 4) Server: wire up all three modules
server <- function(input, output, session) {
  norm_server("norm",
              RGset = reactive({ RGset }), 
              raw_normalised = reactive({ raw_normalised }), 
              targets = reactive({ targets })
  )
}


shinyApp(ui, server)