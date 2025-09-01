# Author: Ghazal Sinjar
# Date: 30.07.2025

# 1) Libraries
library(shiny)
library(bslib)
library(minfi)
library(shinyjs) # for enabling/disabling tabs
library(shinyWidgets)
library(DT) # DT is needed for the filtering module's table output
library(writexl) # For excel downloads in filtering module
library(openxlsx) # For excel downloads in filtering module
library(readr)
library(dplyr) # For data manipulation (filter, select, mutate, bind_rows)
library(stringr)
library(reticulate) # ADDED: Needed for Python integration in tadcalling_module
library(processx) # ADDED: Needed for running external processes (Java TADcaller)
# both libraries are automaticaly loaded when they are needed -> to Do test app without thoes
#library(IlluminaHumanMethylationEPICv2manifest)
#library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(GenomicRanges) # For GenomicRatioSet and related operations in filtering module
library(EnsDb.Hsapiens.v86)
library(fs)
library(purrr)

# Libraries only for GVIZ module
library(Gviz)
library(IRanges)
library(S4Vectors)
library(GenomeInfoDb)
library(BiocGenerics)


# 2) Source your existing modules
source("modules/01_loadData_Module.R")
source("modules/02_QC_Module.R")
source("modules/03_Normalisation_Module.R")
source("modules/04_filtering_Module.R")
source("modules/05_Annotation_Module.R")
source("modules/06_DMR_identification_Module_v1.R")
source("utils/dmrs_utils.R")
source("utils/preprocessing_utils.R")
source("modules/07_boxplots_Module.R")
source("utils/dmrs_boxplot_utils.R")
source("modules/08_offspotter_results_processing_module.R")
source("modules/09_TADcalling_module_V1_all_Chr.R")
source("utils/TADcalling_utils.R")
source("utils/GVIZ_plot_utils.R")
source("modules/10_GVIZ_plot_Module_V1.R")
source("modules/11_dmp_module.R")
source("utils/DMP_utils.R")


#3) This is a static object that will be passed to the DMR module
message("Loading/Preparing tx_gr_filtered for annotation...")
edb <- EnsDb.Hsapiens.v86
options(ucscChromosomeNames = TRUE)
tx_gr <- genes(edb)
tx_gr_filtered_global <- keepSeqlevels(tx_gr, standardChromosomes(tx_gr), pruning.mode = "coarse")
seqlevelsStyle(tx_gr_filtered_global) <- "UCSC"
message("Global: tx_gr_filtered prepared.")
#-------------------------------------------------------------------------
# 4) Chromosome Lengths Table from BSgenome
library(BSgenome.Hsapiens.UCSC.hg38)
hg38 <- BSgenome.Hsapiens.UCSC.hg38
chr_lengths <- seqlengths(hg38)
chr_size_df_global <- data.frame(
  Chromosome = names(chr_lengths),
  Length = as.numeric(chr_lengths),
  stringsAsFactors = FALSE
)
print("Global: chr_size_df_global loaded.")
#------------------------------------------------------------------
# 5) CpG Islands - Loaded once globally
# Uses the loadCpGIslands_gr function from utils/GVIZ_plot_utils.R
gr_cpgIslands_global <- loadCpGIslands_gr(destfile = "cpgIslandExt_hg38.txt.gz")
print(paste("Global: gr_cpgIslands_global loaded. Number of CpG Islands:", length(gr_cpgIslands_global)))
#--------------------------------------------------------------------

# 6) UI
ui <- navbarPage(
  id = "main_tabs",
  title = "epiTADGuide",
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
           actionButton("to_dmp", "Next → DMP Identification")
  ),
  
  tabPanel("DMP Identification",
           dmp_UI("dmp_module_id") ,
           actionButton("to_boxplots", "Next → DMR Boxplots")
  ),
  
  tabPanel("DMR Boxplots",
           boxplotUI("boxplot_module_id"),
           actionButton("to_offtargets", "Next → Offtargets Import")
  ),
  tabPanel("Offtargets Import",
           offtargetsUI("myOfftargetModule"),
           actionButton("to_tadcalling", "Next → TAD Calling")
  ),
  tabPanel("TAD Calling",
           tadcalling_ui("my_tadcalling_module"),
           actionButton("to_gvizplot", "Next → final plot")
  ),
  tabPanel("Visulazation",
           GvizPlotUI("myGvizPlot")
  )
)

#-------------------------------------------------------------------------
# 7) server
server <- function(input, output, session) {
  
  # --- Helper Functions ---
  disable_tab <- function(tab_name) {
    runjs(sprintf("$('li a[data-value=\"%s\"]').addClass('disabled');", tab_name))
  }
  
  enable_tab <- function(tab_name) {
    runjs(sprintf("$('li a[data-value=\"%s\"]').removeClass('disabled');", tab_name))
  }
  
  #------------------------------------------------
  # --- Disable tab navigation initially ---
  disable("main_tabs") 
  #-----------------------------------------
 
   # disabeling tabs and next bottuns initially
  observe({
    print("DEBUG: Main App: Initializing UI and disabling tabs/buttons.")
    # Disable all except Offtargets Import, tadcalling and Load Raw Data tab
    tabs_to_disable <- c("QC Plots", "Normalisation", "Filtering", "Annotation", "DMR Identification", "DMP Identification","DMR Boxplots", "Final plot")
    lapply(tabs_to_disable, disable_tab)
    
    # Offtargets Import remains enabled always, so no disable call here
    
    # Disable all "Next" buttons initially
    shinyjs::disable("to_qc")
    # deleted to_norm since this step is quick and not obligatory for the normalisation
    shinyjs::disable("to_filter")
    shinyjs::disable("to_annot")
    shinyjs::disable("to_dmrs")
    shinyjs::disable("to_dmp")
    shinyjs::disable("to_boxplots")
    shinyjs::disable("to_offtargets")
    #shinyjs::disable("to_tadcalling") # commented cause it can work without off-targets -> they wont be ploted if not available but every thing else will be
    shinyjs::disable("to_gvizplot")   
  })
  #-------------------------------------------
  # defining reactive values
  normalized_output <- reactiveVal(NULL)
  filter_results <- reactiveVal(NULL)
  annotation_results <- reactiveVal(NULL)
  dmr_results <- reactiveVal(NULL)
  offtargets_results <- reactiveVal(NULL)
  tadcalling_results <- reactiveVal(NULL)
  dmp_results <- reactiveVal(NULL)
  boxplot_results <- reactiveVal(NULL)
  #--------------------------------
  
  # --- Module Server Calls and Data Flow ---
  # 1- Load Data Module
  # This module returns a list of reactive values: RGset, raw_normalised, targets, project_dir
  print("DEBUG: Main App: Calling loadDataServer.")
  loaded_data <- loadDataServer("loader")
  
  # OBSERVE targets from loaded_data and enable "Next -> QC" when it's created
  observe({
    if (!is.null(loaded_data$targets())) {
      shinyjs::enable("to_qc")
    } else {
      shinyjs::disable("to_qc")
    }
  })
  
  # ---  Navigate to QC tab ---
  observeEvent(input$to_qc, {
    print("DEBUG: Main App: 'to_qc' button clicked.")
    req(loaded_data$targets())
    enable_tab("QC Plots")
    updateNavbarPage(session, "main_tabs", selected = "QC Plots")
    
    # Load QC server logic
    qcServer("qc",
             RGset = reactive({ loaded_data$RGset() }),
             raw_normalised = reactive({ loaded_data$raw_normalised() }),
             targets = reactive({ loaded_data$targets() }),
             project_output_dir = reactive({ loaded_data$project_dir() })
    )
  })
  
  #-----------------------------------------------
  # Navigate to Normalisation tab
  observeEvent(input$to_norm, {
    req(loaded_data$RGset(), loaded_data$raw_normalised(), loaded_data$targets()) 
    enable_tab("Normalisation")
    updateNavbarPage(session, "main_tabs", selected = "Normalisation")
    
    # Load the Normalisation server logic
    norm_results_local <- norm_server("norm",
                                      RGset = reactive({ loaded_data$RGset() }),
                                      raw_normalised = reactive({ loaded_data$raw_normalised() }),
                                      targets = reactive({ loaded_data$targets() }),
                                      project_output_dir = reactive({ loaded_data$project_dir() })
    )
    normalized_output(norm_results_local) 
  
    # Enable "Next -> Filtering" button only when normalized_output has data
    observe({
      if (!is.null(normalized_output())) {
        shinyjs::enable("to_filter")
      } else {
        shinyjs::disable("to_filter")
      }
    })
  })
  
  #-----------------------------------------------------------
  # Navigate to Filtering tab and run filtering module
  observeEvent(input$to_filter, {
    req(normalized_output()) # Ensure normalized data is available
    
    enable_tab("Filtering")
    updateNavbarPage(session, "main_tabs", selected = "Filtering")
    
    res <- filter_data_server(
      "filter_module",
      RGset = reactive({ loaded_data$RGset() }),
      raw_normalised = reactive({ loaded_data$raw_normalised() }),
      normalized_all_methods = reactive({ normalized_output() }),
      project_output_dir = reactive({ loaded_data$project_dir() })
    )
    
    filter_results(res)
    
    # Enable "Next -> Annotation" button only when filter_results has data
    observe({
      if (!is.null(filter_results()$filtered_data())) { 
        shinyjs::enable("to_annot")
      } else {
        shinyjs::disable("to_annot")
      }
    })
  })
  #---------------------------------------------------------------
  # Navigate to Annotation tab when clicking "Next → Annotation"
  observeEvent(input$to_annot, {
    req(filter_results()) # Ensure filtered data is available
    
    enable_tab("Annotation")
    updateNavbarPage(session, "main_tabs", selected = "Annotation")
    
    # Pass the reactive filtered GenomicRatioSet to annotation module
    annotation_results_local <- annotationServer("annot_module", 
                                                 grset_reactive = filter_results()$filtered_data,
                                                 project_output_dir=  reactive({ loaded_data$project_dir() }))
    annotation_results(annotation_results_local) # Store results in reactiveVal
    
    # Enable "Next -> DMR Identification" button only when annotation_results has data
    # Removed `once = TRUE`
    observe({
      if (!is.null(annotation_results()$annotated_table())) { # Assuming annotation_results() will be non-NULL upon successful annotation
        shinyjs::enable("to_dmrs")
      } else {
        shinyjs::disable("to_dmrs")
      }
    })
  })
  
  #---------------------------------------------------
  # --- Navigate to DMR Identification tab ---
  observeEvent(input$to_dmrs, {
    req(filter_results(), annotation_results()) 
    
    enable_tab("DMR Identification")
    updateNavbarPage(session, "main_tabs", selected = "DMR Identification")
    dmr_results_local <- dmrs_server(
      id = "dmrs_module_id",
      filtered_rgset_reactive = filter_results()$filtered_data,
      tx_gr_filtered_static = tx_gr_filtered_global,
      project_output_dir = reactive({ loaded_data$project_dir() })
    )
    dmr_results(dmr_results_local)
    
    observe({
      if (!is.null(dmr_results()$dmr_table())) { 
        shinyjs::enable("to_dmp")
      } else {
        shinyjs::disable("to_dmp")
      }
    })
  })
  #---------------------------------------------------------
  # --- Navigate to DMP Identification tab ---
  observeEvent(input$to_dmp, {
    req(dmr_results(), filter_results()) 
    
    enable_tab("DMP Identification")
    updateNavbarPage(session, "main_tabs", selected = "DMP Identification")
    dmp_results_local <- dmp_Server(
      id = "dmp_module_id",
      filtered_data = filter_results()$filtered_data,
      pheno = dmr_results()$pheno,
      ref_group = dmr_results()$ref_group,
      project_output_dir = reactive({ loaded_data$project_dir() })
    )
    dmp_results(dmp_results_local)
    observe({
      if (!is.null(dmp_results())) { 
        shinyjs::enable("to_boxplots")
      } else {
        shinyjs::disable("to_boxplots")
      }
    })
  })
  #--------------------------------------------------------------
  # --- Navigate to DMR Boxplots tab ---
  observeEvent(input$to_boxplots, {
    req(dmr_results(), annotation_results())
    enable_tab("DMR Boxplots")
    updateNavbarPage(session, "main_tabs", selected = "DMR Boxplots")
    boxplot_res_local <- boxplotServer(
      id = "boxplot_module_id",
      dmr_output_reactive = dmr_results,
      annotation_tbl_reactive = annotation_results,
      dmp_results = dmp_results,
      project_output_dir =  reactive({ loaded_data$project_dir() })
    )
    boxplot_results(boxplot_res_local)
    # Add the debug observer:
    observe({
      tryCatch({
        print("DEBUG: boxplot_results structure:")
        print(str(boxplot_results()))
        print("DEBUG: First few rows:")
        print(head(boxplot_results()))
      }, error = function(e) {
        message("Debug observer error: ", e$message)
      })
    })
    shinyjs::enable("to_offtargets")
  })
  #----------------------------------------------------------
  # # --- Navigate to off-targets tab ---
  observeEvent(input$to_offtargets, {
    enable_tab("Offtargets Import") # just in case, enable tab
    updateNavbarPage(session, "main_tabs", selected = "Offtargets Import")
  })
  project_dir_for_offtargets <- reactive({ # Renamed from tadcalling to be clear
    if (!is.null(loaded_data$project_dir())) {
      loaded_data$project_dir()
    } else {
      NULL # Pass NULL if not available, which triggers the module's fallback
    }
  })  # Initialize the offtargetsServer module and store its reactive outputs
  offtargets_results_returned <- offtargetsServer(
      id = "myOfftargetModule",
      project_output_dir=  project_dir_for_offtargets
  )
  
  # Assign the returned list of reactives to the reactiveVal
  offtargets_results(offtargets_results_returned)
  shinyjs::enable("to_tadcalling")


  
  # --- Navigate to TADcalling tab ---
  observeEvent(input$to_tadcalling, {
    enable_tab("TAD Calling") # Enable the TAD Calling tab
    updateNavbarPage(session, "main_tabs", selected = "TAD Calling")
  })
  project_dir_for_tadcalling <- reactive({
    if (!is.null(loaded_data$project_dir())) {
      loaded_data$project_dir()
    } else {
      NULL # Pass NULL if not available, which triggers the module's fallback
    }
  })
  tadcalling_results_local <- tadcalling_server(
    "my_tadcalling_module",
    project_output_dir = project_dir_for_tadcalling
  )
  tadcalling_results(tadcalling_results_local)
  
  observe({
    if (!is.null(tadcalling_results()$processed_chroms_list_rv)) {
      shinyjs::enable("to_gvizplot")
    } else {
      shinyjs::disable("to_gvizplot")
    }
  })
  
  
  # --- Navigate to final GVIZ plot tab ---
  observeEvent(input$to_gvizplot, {
    enable_tab("Final plot")
    updateNavbarPage(session, "main_tabs", selected = "Final plot")
  })
  
  GvizPlotServer(
    id = "myGvizPlot",
    dmr_results = dmr_results,
    dmp_results = dmp_results,
    boxplot_results = boxplot_results,
    offtarget_table = offtargets_results,
    tadcalling_results = tadcalling_results,
    chr_size_df_global = chr_size_df_global,
    tx_gr_filtered_global = tx_gr_filtered_global,
    gr_cpgIslands_global = gr_cpgIslands_global
  )
}

shinyApp(ui, server)

