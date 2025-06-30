# 1) Libraries
library(shiny)
library(bslib)
library(minfi)
library(shinyjs)
library(DT)
library(writexl)
library(openxlsx)
library(GenomicRanges)
library(shinyWidgets)

# 2) Source modules
source("modules/01_loadData_Module.R")
source("modules/02_QC_Module.R")
source("modules/03_Normalisation_Module.R")
source("modules/04_filtering_Module.R") 
source("modules/05_Annotation_Module.R")
source("modules/06_DMR_identification_Module.R")
source("utils/dmrs_utils.R")
source("utils/preprocessing_utils.R")
source("modules/07_boxplots_Module.R")
source("utils/dmrs_boxplot_utils.R")

# 3) Prepare static annotation data
message("Loading/Preparing tx_gr_filtered for annotation...")
edb <- EnsDb.Hsapiens.v86
options(ucscChromosomeNames = TRUE)
tx_gr <- genes(edb)
tx_gr_filtered_global <- keepSeqlevels(tx_gr, standardChromosomes(tx_gr), pruning.mode = "coarse")
seqlevelsStyle(tx_gr_filtered_global) <- "UCSC"
message("tx_gr_filtered prepared.")

# 4) UI 
ui <- navbarPage(
  id = "main_tabs",
  title = "EPIC Array Pipeline",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  
  tabPanel("Load Raw Data",
           useShinyjs(),
           loadDataUI("loader"),
           actionButton("to_qc", "Next → QC")
  ),
  
  tabPanel("QC Plots",
           qcUI("qc"),
           actionButton("to_norm", "Next → Normalisation")
  ),
  
  tabPanel("Normalisation",
           norm_ui("norm"),
           actionButton("to_filter", "Next → Filtering") 
  ),
  
  tabPanel("Filtering", 
           filter_data_ui("filter_module"),
           actionButton("to_annot", "Next → Annotation")
  ),
  
  tabPanel("Annotation",
           annotationUI("annot_module"),
           actionButton("to_dmrs", "Next → DMR Identification")
  ),
  
  tabPanel("DMR Identification", 
           dmrs_ui("dmrs_module_id"),
           actionButton("to_boxplots", "Next → DMR Boxplots")
  ),
  
  tabPanel("DMR Boxplots",
           boxplotUI("boxplot_module_id")
  )
)

# 5) Server 
server <- function(input, output, session) {
  # --- Centralized App State ---
  app_data <- reactiveValues(
    loaded_data = NULL,
    qc_initialized = FALSE,
    norm_initialized = FALSE,
    filter_initialized = FALSE,
    annot_initialized = FALSE,
    dmr_initialized = FALSE,
    boxplot_initialized = FALSE,
    
    # Reactive storage
    normalized_output = NULL,
    filter_results = NULL,
    annotation_results = NULL,
    dmr_results = NULL
  )
  
  # --- Helper Functions ---
  disable_tab <- function(tab_name) {
    runjs(sprintf("$('li a[data-value=\"%s\"]').addClass('disabled');", tab_name))
  }
  
  enable_tab <- function(tab_name) {
    runjs(sprintf("$('li a[data-value=\"%s\"]').removeClass('disabled');", tab_name))
  }
  
  # --- Initialize Only the First Module ---
  app_data$loaded_data <- loadDataServer("loader")
  
  # --- Navigation Logic ---
  
  # 1. Load Data -> QC
  observeEvent(input$to_qc, {
    req(app_data$loaded_data$RGset())  # Basic requirement
    
    tryCatch({
      if (!app_data$qc_initialized) {
        qcServer("qc",
                 RGset = reactive({ app_data$loaded_data$RGset() }),
                 raw_normalised = reactive({ app_data$loaded_data$raw_normalised() }),
                 targets = reactive({ app_data$loaded_data$targets() }))
        app_data$qc_initialized <- TRUE
      }
      
      enable_tab("QC Plots")
      updateNavbarPage(session, "main_tabs", selected = "QC Plots")
    }, error = function(e) {
      showNotification(paste("Failed to initialize QC:", e$message), type = "error")
    })
  })
  
  # 2. QC -> Normalization
  observeEvent(input$to_norm, {
    req(app_data$loaded_data$RGset())  # Basic requirement
    
    tryCatch({
      if (!app_data$norm_initialized) {
        norm_mod <- norm_server("norm",
                                RGset = reactive({ app_data$loaded_data$RGset() }),
                                raw_normalised = reactive({ app_data$loaded_data$raw_normalised() }),
                                targets = reactive({ app_data$loaded_data$targets() }))
        
        app_data$normalized_output <- reactive({
          req(norm_mod())
          norm_mod()
        })
        
        app_data$norm_initialized <- TRUE
      }
      
      enable_tab("Normalisation")
      updateNavbarPage(session, "main_tabs", selected = "Normalisation")
    }, error = function(e) {
      showNotification(paste("Normalization Error:", e$message), type = "error")
    })
  })
  
  # 3. Normalization -> Filtering
  observeEvent(input$to_filter, {
    req(app_data$normalized_output())  # Must have normalization results
    
    tryCatch({
      if (!app_data$filter_initialized) {
        filter_mod <- filter_data_server(
          "filter_module",
          RGset = reactive({ app_data$loaded_data$RGset() }),
          raw_normalised = reactive({ app_data$loaded_data$raw_normalised() }),
          normalized_all_methods = reactive({ app_data$normalized_output() })
        )
        
        app_data$filter_results <- reactive({
          req(filter_mod())
          filter_mod()
        })
        
        app_data$filter_initialized <- TRUE
      }
      
      enable_tab("Filtering")
      updateNavbarPage(session, "main_tabs", selected = "Filtering")
    }, error = function(e) {
      showNotification(paste("Filtering Error:", e$message), type = "error")
    })
  })
  
  # 4. Filtering -> Annotation
  observeEvent(input$to_annot, {
    req(app_data$filter_results())
    
    tryCatch({
      if (!app_data$annot_initialized) {
        annot_mod <- annotationServer(
          "annot_module",
          grset_reactive = reactive({ 
            req(app_data$filter_results())
            app_data$filter_results()$filtered_data 
          })
        )
        
        app_data$annotation_results <- reactive({
          req(annot_mod())
          annot_mod()
        })
        
        app_data$annot_initialized <- TRUE
      }
      
      enable_tab("Annotation")
      updateNavbarPage(session, "main_tabs", selected = "Annotation")
    }, error = function(e) {
      showNotification(paste("Annotation Error:", e$message), type = "error")
    })
  })
  
  # 5. Annotation -> DMR Identification
  observeEvent(input$to_dmrs, {
    req(app_data$annotation_results())
    
    tryCatch({
      if (!app_data$dmr_initialized) {
        dmr_mod <- dmrs_server(
          "dmrs_module_id",
          filtered_rgset_reactive = reactive({
            req(app_data$filter_results())
            app_data$filter_results()$filtered_data
          }),
          tx_gr_filtered_static = tx_gr_filtered_global
        )
        
        app_data$dmr_results <- reactive({
          req(dmr_mod())
          dmr_mod()
        })
        
        app_data$dmr_initialized <- TRUE
      }
      
      enable_tab("DMR Identification")
      updateNavbarPage(session, "main_tabs", selected = "DMR Identification")
    }, error = function(e) {
      showNotification(paste("DMR Identification Error:", e$message), type = "error")
    })
  })
  
  # 6. DMR Identification -> DMR Boxplots
  observeEvent(input$to_boxplots, {
    req(app_data$dmr_results())
    
    tryCatch({
      if (!app_data$boxplot_initialized) {
        boxplotServer(
          "boxplot_module_id",
          dmr_output_reactive = reactive({
            req(app_data$dmr_results())
            app_data$dmr_results()
          }),
          annotation_output_reactive = reactive({
            req(app_data$annotation_results())
            app_data$annotation_results()
          })
        )
        app_data$boxplot_initialized <- TRUE
      }
      
      enable_tab("DMR Boxplots")
      updateNavbarPage(session, "main_tabs", selected = "DMR Boxplots")
    }, error = function(e) {
      showNotification(paste("Boxplot Error:", e$message), type = "error")
    })
  })
  
  # --- Initial Setup ---
  observe({
    disable("main_tabs")
    lapply(c("QC Plots", "Normalisation", "Filtering", 
             "Annotation", "DMR Identification", "DMR Boxplots"), disable_tab)
  })
  
  # --- Debugging Outputs ---
  output$debug_state <- renderPrint({
    list(
      loaded_data = !is.null(app_data$loaded_data),
      normalized_output = !is.null(try(app_data$normalized_output(), silent = TRUE)),
      filter_results = !is.null(try(app_data$filter_results(), silent = TRUE)),
      annotation_results = !is.null(try(app_data$annotation_results(), silent = TRUE)),
      dmr_results = !is.null(try(app_data$dmr_results(), silent = TRUE))
    )
  })
}

shinyApp(ui, server) 