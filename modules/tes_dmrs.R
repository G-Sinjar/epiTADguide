'# test 06-DMR module
# libraries especially for this module
library(DT)        # For interactive data tables
library(openxlsx)
library(GenomicRanges)
library(shiny) # Add shiny explicitly here for module to be self-contained in its library requirements
library(minfi) # Add minfi explicitly here
library(bslib) # Add bslib explicitly here
library(EnsDb.Hsapiens.v86)
library(shinyjs)


## load input data
filtered_rgset <- readRDS("../intermediate_data/filtered_GRset_SWAN_SNPsremoved_SexChrProbes_kept_20250608.rds")

# Load custom utility functions for DMR processing
# In a module, you might consider if these utils should always loaded by the main app. For now, well keep the source here.
# Adjust path for sourcing relative to the module file itself
source("../utils/dmrs_utils.R")
source("./06_DMR_identification_Module_v1.R")

## creating the gene granges for the input
library(EnsDb.Hsapiens.v86)
edb <- EnsDb.Hsapiens.v86
options(ucscChromosomeNames = TRUE)
#To extract specific data, such as transcript or gene information, from an EnsDb object, you can use functions like transcripts(), genes(), and exons().
######################### Get gene ranges
tx_gr <- genes(edb)
#length(tx_gr) # 63970
head(tx_gr)
# Filter to standard chromosomes only????
tx_gr_filtered <- keepSeqlevels(tx_gr, standardChromosomes(tx_gr), pruning.mode = "coarse")
#length(tx_gr_filtered) #58650
#tx_gr_filtered
seqlevelsStyle(tx_gr_filtered) <- "UCSC"
#tx_gr_filtered

# Check if the directory exists
dir_path <- "../main_app_tests/epic-test"

if (dir.exists(dir_path)) {
  print(paste("The directory exists:", dir_path))
} else {
  print(paste("The directory does not exist:", dir_path))
}


# UI
ui <- page_navbar(
  title = "EPIC Array Pipeline",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  
  nav_panel("DMR identification",
            dmrs_ui("dmrs")
  )
)
# Server
server <- function(input, output, session) {
  # Wrap your object as a reactive expression
  filtered_rgset_reactive <- reactive({ filtered_rgset })
  
  # Call the module
  annotated_table <- dmrs_server("dmrs", 
                                 filtered_rgset_reactive,
                                 tx_gr_filtered,
                                 project_output_dir = reactive({dir_path}))
}
shinyApp(ui, server)'




#--------------------------------------------------
# test both  06-07 modules
# libraries especially for this module
library(DT)        # For interactive data tables
library(openxlsx)
library(GenomicRanges)
library(shiny) # Add shiny explicitly here for module to be self-contained in its library requirements
library(minfi) # Add minfi explicitly here
library(bslib) # Add bslib explicitly here
library(EnsDb.Hsapiens.v86)
library(shinyjs)
library(shinyWidgets)       # For enhanced widgets like switchInput

## load input data
filtered_rgset <- readRDS("../intermediate_data/filtered_GRset_SWAN_SNPsremoved_SexChrProbes_kept_20250608.rds")
results_anno <- readRDS("./intermediate_data/annotated_object_20250623.rds")

# Load custom utility functions
source("../utils/dmrs_utils.R")
source("../utils/dmrs_boxplot_utils.R")  # For boxplot module
source("./06_DMR_identification_Module_v1.R")
## creating the gene granges for the input
edb <- EnsDb.Hsapiens.v86
options(ucscChromosomeNames = TRUE)
tx_gr <- genes(edb)
tx_gr_filtered <- keepSeqlevels(tx_gr, standardChromosomes(tx_gr), pruning.mode = "coarse")
seqlevelsStyle(tx_gr_filtered) <- "UCSC"

# Check if the directory exists
dir_path <- "../main_app_tests/epic-test"

# UI
ui <- page_navbar(
  id = "main_nav",
  title = "EPIC Array Pipeline",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  
  nav_panel(
    "DMR Identification",
    dmrs_ui("dmrs_module_id"),
    actionButton("to_boxplots", "Next → DMR Boxplots", class = "btn-primary")
  ),
  
  nav_panel(
    "DMR Boxplots",
    boxplotUI("boxplot_module_id")
  )
)

# ui and data loading remain the same

# Server with corrected reactive handling
server <- function(input, output, session) {
  # Initialize shinyjs
  shinyjs::useShinyjs()
  
  # Disable boxplots button initially
  shinyjs::disable("to_boxplots")
  
  # Create a reactive expression for the input data
  filtered_rgset_reactive <- reactive({ filtered_rgset })
  
  # Call the DMR module server function ONCE at the top level
  # This returns a reactive expression that holds the DMR results
  dmr_results_reactive <- dmrs_server(
    "dmrs_module_id",
    filtered_rgset_reactive = filtered_rgset_reactive,
    tx_gr_filtered_static = tx_gr_filtered,
    project_output_dir = reactive({dir_path})
  )
  
  # Create a reactive for the annotation data that will be passed to the boxplot module
  annotation_output_reactive <- reactive({
    list(annotated_table = results_anno$annotated_table)
  })
  
  # Call the Boxplot module server function ONCE at the top level
  # Pass it the reactive expressions it needs
  # The output of the dmr module is passed as an argument
  # to the boxplot module.
  boxplotServer(
    id = "boxplot_module_id",
    dmr_output_reactive = dmr_results_reactive,
    annotation_output_reactive = annotation_output_reactive
  )
  
  # Enable the "Next" button when the DMR results are ready
  observe({
    req(dmr_results_reactive())
    shinyjs::enable("to_boxplots")
  })
  
  # Navigate to the boxplots tab when the button is clicked
  observeEvent(input$to_boxplots, {
    updateNavbarPage(session, "main_nav", selected = "DMR Boxplots")
  })
}

shinyApp(ui, server)