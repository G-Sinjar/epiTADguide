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


# test for V1 all chr app
# app.R (Main Shiny Application - Minimal Host for TADcalling Module)
library(shiny)
library(bslib)
# they should be loaded within tadcalling_module.R or a global setup script.

# Source the utility functions
source("./utils/TADcalling_utils.R")
source("./modules/09_TADcalling_module_V1_all_Chr.R")

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