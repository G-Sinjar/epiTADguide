# 1) Libraries
library(shiny)
library(bslib)
library(minfi)
library(shinyjs) # for enabling/disabling tabs
# both libraries are automaticaly loaded when they are needed -> to Do test app without thoes
#library(IlluminaHumanMethylationEPICv2manifest) 
#library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(DT) # DT is needed for the filtering module's table output
library(writexl) # For excel downloads in filtering module
library(openxlsx) # For excel downloads in filtering module
library(GenomicRanges) # For GenomicRatioSet and related operations in filtering module

# 2) Source your existing modules
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



# This is a static object that will be passed to the DMR module
message("Loading/Preparing tx_gr_filtered for annotation...")
edb <- EnsDb.Hsapiens.v86
options(ucscChromosomeNames = TRUE)
tx_gr <- genes(edb)
tx_gr_filtered_global <- keepSeqlevels(tx_gr, standardChromosomes(tx_gr), pruning.mode = "coarse")
seqlevelsStyle(tx_gr_filtered_global) <- "UCSC"
message("tx_gr_filtered prepared.")




# 3) UI
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
           dmrs_ui("dmrs_module_id") ,
           actionButton("to_boxplots", "Next → DMR Boxplots")
           
  ),
  tabPanel("DMR Boxplots",
           boxplotUI("boxplot_module_id")
  )
)


# 4) server
server <- function(input, output, session) {
  # --- Disable tab navigation initially ---
  disable("main_tabs") # disables all tab switching via clicks
  
  # Allow only first tab to be selected initially and disable subsequent tabs
  observe({
    runjs("
      $('li a[data-value=\"QC Plots\"]').addClass('disabled');
      $('li a[data-value=\"Normalisation\"]').addClass('disabled');
      $('li a[data-value=\"Filtering\"]').addClass('disabled'); 
      $('li a[data-value=\"Annotation\"]').addClass('disabled');
      $('li a[data-value=\"DMR Identification\"]').addClass('disabled');
      $('li a[data-value=\"DMR Boxplots\"]').addClass('disabled');
    ")
  })
  
  # --- Module Server Calls and Data Flow ---
  
  # Load Data Module
  # This module returns a list of reactive values: RGset, raw_normalised, targets
  loaded_data <- loadDataServer("loader")
  
  normalized_output <- reactiveVal(NULL)
  filter_results <- reactiveVal(NULL)
  annotation_results <- reactiveVal(NULL)
  dmr_results <- reactiveVal(NULL)
  
  
  
  # --- Navigation Logic ---
  
  # Navigate to QC tab
  observeEvent(input$to_qc, {
    # 1. Enable tab
    runjs("$('li a[data-value=\"QC Plots\"]').removeClass('disabled');")
    
    # 2. Switch to the QC tab
    updateNavbarPage(session, "main_tabs", selected = "QC Plots")
    
    # 3. Load QC server logic (deferred execution)
    qcServer("qc",
             RGset = reactive({ loaded_data$RGset() }),
             raw_normalised = reactive({ loaded_data$raw_normalised() }),
             targets = reactive({ loaded_data$targets() })
    )
  })
  
  # Navigate to Normalisation tab
  observeEvent(input$to_norm, {
    # 1. Enable the Normalisation tab
    runjs("$('li a[data-value=\"Normalisation\"]').removeClass('disabled');")
    
    # 2. Switch to the Normalisation tab
    updateNavbarPage(session, "main_tabs", selected = "Normalisation")
    
    # 3. Load the Normalisation server logic
    # This module will return a reactive list containing all normalized methods.
    norm_results <- norm_server("norm",
                                RGset = reactive({ loaded_data$RGset() }),
                                raw_normalised = reactive({ loaded_data$raw_normalised() }),
                                targets = reactive({ loaded_data$targets() })
    )
    normalized_output(norm_results) 
  })
  
  
  
  # Navigate to Filtering tab and run filtering module
  observeEvent(input$to_filter, {
    req(normalized_output())
    
    runjs("$('li a[data-value=\"Filtering\"]').removeClass('disabled');")
    updateNavbarPage(session, "main_tabs", selected = "Filtering")
    
    res <- filter_data_server(
      "filter_module",
      RGset = reactive({ loaded_data$RGset() }),
      raw_normalised = reactive({ loaded_data$raw_normalised() }),
      normalized_all_methods = reactive({ normalized_output() })
    )
    
    filter_results(res)
  })
  
  # Navigate to Annotation tab when clicking "Next → Annotation"
  observeEvent(input$to_annot, {
    req(filter_results())
      
    
    runjs("$('li a[data-value=\"Annotation\"]').removeClass('disabled');")
    updateNavbarPage(session, "main_tabs", selected = "Annotation")
    
    # Pass the reactive filtered GenomicRatioSet to annotation module
    annotation_results(
      annotationServer("annot_module", grset_reactive = filter_results()$filtered_data)
    )
  })
  
'  # If you want to access them:
  annotated_table <- annotation_results$annotated_table
  annotation_object <- annotation_results$annotation_object'
  
  # --- Navigate to DMR Identification tab ---
  observeEvent(input$to_dmrs, {
    req(filter_results()) 
    # No explicit req for tx_gr_filtered_global needed here as it's defined globally
    
    runjs("$('li a[data-value=\"DMR Identification\"]').removeClass('disabled');")
    updateNavbarPage(session, "main_tabs", selected = "DMR Identification")
    
    # Call the DMR module server
    # Pass the filtered GenomicRatioSet from the filtering module
    # Pass the globally prepared tx_gr_filtered object
    dmr_results(
      dmrs_server(
        id = "dmrs_module_id",
        filtered_rgset_reactive = filter_results()$filtered_data,
        tx_gr_filtered_static = tx_gr_filtered_global
      )
    )
  })
  
  
  # --- Navigate to DMR Boxplots tab ---
  observeEvent(input$to_boxplots, {
    runjs("$('li a[data-value=\"DMR Boxplots\"]').removeClass('disabled');")
    updateNavbarPage(session, "main_tabs", selected = "DMR Boxplots")
    boxplotServer(
      id = "boxplot_module_id",
      dmrs_table = dmr_results()$dmr_table(),
      annotated_with_betas_df = annotation_results()$annotated_table(),
      pheno_data = dmr_results()$pheno_data()
    )
  })
  
}

shinyApp(ui, server)