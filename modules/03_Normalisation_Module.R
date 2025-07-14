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
                     "SWAN (Subset-quantile within array normalization)" = "Swan",
                     "NOOB (normal-exponential out-of-band normalisation)" = "Noob",
                     "NOOB.SWAN" = "Noob.Swan",
                     "FunNorm (Functional normalisation)" = "Funnorm",
                     "Stratified quantile normalization" = "Quantile"
                   ),
                   selected = "Swan"
      ),
      hr(),
      radioButtons(ns("plot_format"), "Select Format:", choices = c("PDF" = "pdf", "PNG" = "png"), inline = TRUE),
      downloadButton(ns("download_plot"), "Current normalised Plot")),
    # Placeholder for dynamic UI based on method selected
    uiOutput(ns("beta_plot_container"))
  )
}



# ─────────────────────────────────────
# SERVER FUNCTION 
# ─────────────────────────────────────

norm_server <- function(id, RGset, raw_normalised, targets, project_output_dir) {
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
        
        # Step1: SWAN Normalization
        incProgress(0, detail = "Starting SWAN normalization...")
        
        # Set seed for reproducibility
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
        
        
        # Step2: NOOB Normalization
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
        
        
        # Step3: NOOB.SWAN Normalization (depends on noob_normalised_res)
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
        
        
        # Step4: FUNNORM Normalization
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
        
        
        # Step5: Quantile Normalization
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
        
        
        # Final progress update
        incProgress(0, detail = "All normalizations finished!")
        
      }) # --- End withProgress block ---
      
      
      # saving the normalised data as .rds file
      if (!is.null(SWAN_normalised_res()) && !is.null(Quantile_normalised_res()) &&
          !is.null(funnorm_normalised_res()) && !is.null(noob_normalised_res()) &&
          !is.null(noob_swan_normalised_res())) {
        
        withProgress(message = "Saving results", value = 0, {
          incProgress(0.2, detail = "Preparing data...")
          
          # Wrap everything in tryCatch to handle errors during save
          tryCatch({
            normalised_list <- list(
              SWAN = SWAN_normalised_res(),
              Quantile = Quantile_normalised_res(),
              Funnorm = funnorm_normalised_res(),
              Noob = noob_normalised_res(),
              Noob_Swan = noob_swan_normalised_res()
            )
            
            incProgress(0.4, detail = "saving as .rds file...")
            
            # Set output path
            outputDir <-file.path(project_output_dir(), "intermediate_data")
            # Make sure the directory exists
            if (!dir.exists(outputDir)) {
              dir.create(outputDir, recursive = TRUE)
            }
            # Save the file
            saveRDS(normalised_list, file = file.path(outputDir, "Five_objects_normalised_data.rds"))
            showNotification("✅ All normalizations saved successfully!", type = "message", duration = 10)
            incProgress(1, detail = "All normalizations saved successfully!")
            
          }, error = function(e) {
            showNotification(paste("❌ Error saving normalized data:", e$message), type = "error", duration = 15)            
            message("Save failed: ", e$message)
          })
        })
      }
    }, ignoreNULL = TRUE, ignoreInit = FALSE, once = TRUE)
    
    
    # --- Dynamic UI for Normalisation Plot Display ---
    output$beta_plot_container <- renderUI({
      ns <- session$ns
      
      # Dynamic title and description
      header_text <- switch(input$norm_methode,
                            "Swan" = "SWAN (Subset-quantile within array normalization)",
                            "Noob" = "Noob (normal-exponential out-of-band normalisation)",
                            "Noob.Swan" = "Noob then SWAN",
                            "Funnorm" = "Funnorm (functional normalization)",
                            "Quantile" = "Stratified quantile normalization"
      )
      
      description <- switch(input$norm_methode,
                            "Swan" = "SWAN is a within-array normalization method designed to correct for technical differences between Infinium Type I and Type II probe designs. A common source of such variability is the differing intensity distributions observed between these two probe types.",
                            "Noob" = "Noob is a background correction method with dye-bias normalization for Illumina Infinium methylation arrays. It accounts for technical variation in background fluorescence signal. It uses the Infinium I design bead types to measure non-specific fluorescence in the colour channel opposite of their design (Cy3/Cy5).",
                            "Noob.Swan" = "Noob background correction followed by SWAN to address both background and probe-type differences.",
                            "Funnorm" = "The functional normalization algorithm is an extension to quantile normalization that removes unwanted technical variation using control probes. It is particularly useful for studies comparing conditions with known large-scale differences, such as cancer/normal studies, or between-tissue studies. By default, the function applies the Noob function as a first step for background substraction, and uses the first two principal components of the control probes to infer the unwanted variation. This method outperforms all existing normalization methods with respect to replication of results between experiments, and yields robust results even in the presence of batch effects.",
                            "Quantile" = "This normalization procedure is applied to the Meth and Unmeth intensities separately. The distribution of type I and type II signals is forced to be the same by first quantile normalizing the type II probes across samples and then interpolating a reference distribution to which the type I probes are normalised. Since probe types and probe regions are confounded and it is known that DNAm (methylation) distributions vary across regions, the probes are stratified by region  (CpG island, shore, etc.) before applying this interpolation. the functions fixes outliers of both the methylated and unmethylated channels when small intensities are close to zero. It also removes bad samples using the QC criterion discussed in Plot2.\nNote that this algorithm is not recommended for cases where global changes are expected such as in cancer-normal comparisons."
      )
      
      # Choose plot output ID
      plot_output_id <- switch(input$norm_methode,
                               "Swan" = "Swan",
                               "Noob" = "Noob",
                               "Noob.Swan" = "Noob.Swan",
                               "Funnorm" = "Funnorm",
                               "Quantile" = "Quantile"
      )
      
      tagList(
        h3("Raw Beta values distribution"),
        helpText("No normalization is applied here. This shows the initial distribution of methylation signal intensities."),
        fluidRow(column(2), column(8, plotOutput(ns("Raw"))), column(2)),
        hr(),
        h3(header_text),
        helpText(description),
        fluidRow(column(2), column(8, plotOutput(ns(plot_output_id))), column(2))
      )
    })
    
    # --- Density Plots (each wrapped in withProgress) ---
    output$Raw <- renderPlot({
      req(raw_normalised(), targets())
      withProgress(message = 'Generating Raw Data Density Plot...', value = 0, {
        incProgress(0.5, detail = "Plotting raw data")
        densityPlot(raw_normalised(), sampGroups = targets()$Sample_Label,
                    #main = "Beta Values distribution of raw data",
                    pal = rainbow(length(unique(targets()$Sample_Label))))
        incProgress(1, detail = "Plot ready")
      })
    })
    
    output$Swan <- renderPlot({
      req(input$norm_methode == "Swan", SWAN_normalised_res(), targets()) 
      withProgress(message = 'Generating SWAN Density Plot...', value = 0, {
        incProgress(0.5, detail = "Plotting SWAN normalized data")
        densityPlot(SWAN_normalised_res(), sampGroups = targets()$Sample_Label,
                    #main = "Beta values distribution after SWAN",
                    pal = rainbow(length(unique(targets()$Sample_Label))))
        incProgress(1, detail = "Plot ready")
      })
    })
    
    output$Noob <- renderPlot({
      req(input$norm_methode == "Noob", noob_normalised_res(), targets()) 
      withProgress(message = 'Generating Noob Density Plot...', value = 0, {
        incProgress(0.5, detail = "Plotting Noob normalized data")
        densityPlot(noob_normalised_res(), sampGroups = targets()$Sample_Label,
                    #main = "Beta values distribution after Noob",
                    pal = rainbow(length(unique(targets()$Sample_Label))))
        incProgress(1, detail = "Plot ready")
      })
    })
    
    output$Noob.Swan <- renderPlot({
      req(input$norm_methode == "Noob.Swan", noob_swan_normalised_res(), targets()) 
      withProgress(message = 'Generating Noob.Swan Density Plot...', value = 0, {
        incProgress(0.5, detail = "Plotting Noob.Swan normalized data")
        densityPlot(noob_swan_normalised_res(), sampGroups = targets()$Sample_Label,
                    #main = "Beta values distribution after Noob-SWAN",
                    pal = rainbow(length(unique(targets()$Sample_Label))))
        incProgress(1, detail = "Plot ready")
      })
    })
    
    output$Funnorm <- renderPlot({
      req(input$norm_methode == "Funnorm", funnorm_normalised_res(), targets()) 
      withProgress(message = 'Generating Funnorm Density Plot...', value = 0, {
        incProgress(0.5, detail = "Plotting Funnorm normalized data")
        densityPlot(getBeta(funnorm_normalised_res()), sampGroups = targets()$Sample_Label,
                    #main = "Beta values distribution after Funnorm",
                    pal = rainbow(length(unique(targets()$Sample_Label))))
        incProgress(1, detail = "Plot ready")
      })
    })
    
    output$Quantile <- renderPlot({
      req(input$norm_methode == "Quantile", Quantile_normalised_res(), targets()) 
      withProgress(message = 'Generating Quantile Density Plot...', value = 0, {
        incProgress(0.5, detail = "Plotting Quantile normalized data")
        densityPlot(getBeta(Quantile_normalised_res()), sampGroups = targets()$Sample_Label,
                    #main = "Beta values distribution after Quantile",
                    pal = rainbow(length(unique(targets()$Sample_Label))))
        incProgress(1, detail = "Plot ready")
      })
    })
    
    output$download_plot <- downloadHandler(
      filename = function() {
        timestamp <- format(Sys.time(), "%Y%m%d")
        paste0("BetaPlot_", input$norm_methode, "_", timestamp, ".", input$plot_format)
      },
      content = function(file) {
        # Choose graphics device based on file type
        if (input$plot_format == "png") {
          png(file, width = 1700, height = 1000)
        } else {
          pdf(file, width = 12, height = 7)
        }
        
        switch(input$norm_methode,
               "Swan" = {
                 req(SWAN_normalised_res(), targets())
                 densityPlot(SWAN_normalised_res(), sampGroups = targets()$Sample_Label,
                             main = "Beta values distribution after SWAN",
                             pal = rainbow(length(unique(targets()$Sample_Label))))
               },
               "Noob" = {
                 req(noob_normalised_res(), targets())
                 densityPlot(noob_normalised_res(), sampGroups = targets()$Sample_Label,
                             main = "Beta values distribution after Noob",
                             pal = rainbow(length(unique(targets()$Sample_Label))))
               },
               "Noob.Swan" = {
                 req(noob_swan_normalised_res(), targets())
                 densityPlot(noob_swan_normalised_res(), sampGroups = targets()$Sample_Label,
                             main = "Beta values distribution after Noob-SWAN",
                             pal = rainbow(length(unique(targets()$Sample_Label))))
               },
               "Funnorm" = {
                 req(funnorm_normalised_res(), targets())
                 densityPlot(getBeta(funnorm_normalised_res()), sampGroups = targets()$Sample_Label,
                             main = "Beta values distribution after Funnorm",
                             pal = rainbow(length(unique(targets()$Sample_Label))))
               },
               "Quantile" = {
                 req(Quantile_normalised_res(), targets())
                 densityPlot(getBeta(Quantile_normalised_res()), sampGroups = targets()$Sample_Label,
                             main = "Beta values distribution after Quantile",
                             pal = rainbow(length(unique(targets()$Sample_Label))))
               }
        )
        
        dev.off()
      }
    )
    
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


'# 1) Libraries
library(shiny)  
library(bslib)  
library(minfi)  
library(IlluminaHumanMethylationEPICv2manifest)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)

#source("./03_Normalisation_Module.R")
path <- "C:/Users/ghaza/Documents/ghazal/Bioinformatik_Fächer/Masterarbeit_Project/Scripts/R_Scripts/intermediate_data/preprocessed_data.rds"
file.exists("../intermediate_data/preprocessed_data.rds")


preprocessed <- readRDS(path)
RGset <- preprocessed$RGset
targets <- preprocessed$targets
raw_normalised <- preprocessed$raw_normalised
path <- "./epic-test"

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
              targets = reactive({ targets }),
              project_output_dir = reactive({path})
  )
}


shinyApp(ui, server)'