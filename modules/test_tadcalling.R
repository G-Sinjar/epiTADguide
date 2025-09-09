'# app.R (Main Shiny Application - Minimal Host for TADcalling Module)
library(shiny)
library(bslib)
# they should be loaded within tadcalling_module.R or a global setup script.

# Source the utility functions
source("./utils/TADcalling_utils.R")
source("./modules/09_TADcalling_module.R")

ui <- page_navbar(
  title = "EPIC Array Pipeline",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  
  nav_panel(
    "TAD Calling Analysis",
    tadcalling_ui("my_tad_analysis_module")
  )
)

server <- function(input, output, session) {
  #project_dr <- reactive({"./epic-test"})
  tad_results <- tadcalling_server("my_tad_analysis_module")    #, project_output_dir= project_dr)
}
shinyApp(ui, server)'

#---------------------------------------------
# test for V1 all chr app
library(shiny)          # Core Shiny functionality
library(shinyWidgets)   # For advanced UI components (if used)
library(shinyjs)        # Optional, for enabling/disabling UI elements dynamically
library(readr)          # read_delim, write_tsv
library(tools)          # file_path_sans_ext
library(dplyr)          # Optional, if you do further data manipulations
library(stringr)        # Optional, string handling
library(tidyr)          # Optional, for data reshaping
library(fs)             # Optional, for cross-platform file handling
library(processx)       # Running external commands (Java deDoc2)
library(reticulate)     # use_virtualenv, source_python
library(DT)             # DT::dataTableOutput and renderDataTable
library(bslib)          # If using Bootstrap layouts
library(openxlsx)       # createWorkbook, addWorksheet, writeData, saveWorkbook

# app.R or global.R
options(shiny.maxRequestSize = 3 * 1024^3)  # 3 GB


# Source the utility functions
source("../utils/TADcalling_utils.R")
source("./09_TADcalling_module_ya.R")

ui <- page_navbar(
  title = "EPIC Array Pipeline",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  
  nav_panel(
    "TAD Calling Analysis",
    tadcalling_ui("my_tad_analysis_module")
  )
)

server <- function(input, output, session) {
  project_dr <- reactive({"./epic-test"})
  tad_results <- tadcalling_server("my_tad_analysis_module",
                                   project_output_dir= project_dr)
}
shinyApp(ui, server)