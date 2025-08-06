# dmp_app.R
# inputs



# Load necessary libraries
library(shiny)
library(bslib)
library(DT)
library(minfi) # For dmpFinder and getM functions

# Source the utility functions
source("../utils/DMP_utils.R")

# --- UI for the DMP App ---
ui <- page_navbar(
  title = "DMP Finder App",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  
  nav_panel("DMP Identification",
            layout_sidebar(
              sidebar = sidebar(
                width = 300,
                
                # q-value cutoff numeric input
                numericInput("qvalue_cutoff", "q-value Cutoff (FDR):",
                             value = 0.05, min = 0, max = 1, step = 0.01
                ),
                
                # Reference group display (static)
                h5(strong("Reference Group:")),
                textOutput("display_ref_group"),
                
                # Normalization method display (static)
                h5(strong("Normalization Method:")),
                textOutput("display_norm_method"),
                
                helpText("Differences in methylation values are compared to methylation values of this reference group. (-) means hypomethylated, (+) means hypermethylated against this group."),
                
                # Action button to run the analysis
                actionButton("run_dmp_finder", "Run DMP Finder", class = "btn-primary",
                             icon = icon("play-circle")),
                
                div(style = "margin-top: 20px;"),
                # Action link for notes on DMP table columns
                actionLink("notes_link", "Notes on reading the DMP table"),
                hr(),
                
                # Placeholder for download buttons (conditionally rendered)
                uiOutput("download_ui")
              ),
              main = card(
                card_header("Process Status"),
                textOutput("status_message")
              ),
              # DMP table goes under the card
              DT::dataTableOutput("dmp_table")
            )
  )
)

# --- Server for the DMP App ---
server <- function(input, output, session) {
  # Reactive values to store data (to be set by the user or parent app)
  normalisation_results <- reactiveVal(NULL) # List of methylsets
  pheno <- reactiveVal(NULL) # Phenotype data
  ref_group <- reactiveVal(NULL) # Reference group
  norm_method <- reactiveVal(NULL) # Normalization method
  
  # Display the reference group
  output$display_ref_group <- renderText({
    req(ref_group())
    ref_group()
  })
  
  # Display the normalization method
  output$display_norm_method <- renderText({
    req(norm_method())
    norm_method()
  })
  
  # Reactive value to store the combined M and Beta mean table
  combined_m_beta_tbl <- reactiveVal(NULL)
  # Reactive value to store the final edited DMP table
  dmp_final_table <- reactiveVal(NULL)
  
  # Observe the "Run DMP Finder" button click
  observeEvent(input$run_dmp_finder, {
    # Reset previous results and status messages
    combined_m_beta_tbl(NULL)
    dmp_final_table(NULL)
    output$status_message <- renderText("Starting DMP analysis...")
    output$download_ui <- renderUI(NULL)
    
    # Validate required inputs
    req(normalisation_results(), pheno(), ref_group(), norm_method())
    
    # Use withProgress to show analysis progress
    withProgress(message = "DMP Analysis in Progress", value = 0, {
      
      # --- Step 1: Calculate Methylation Mean Differences ---
      incProgress(0.2, detail = "Calculating mean differences for Beta and M values...")
      output$status_message <- renderText("Step 1/3: Calculating mean differences...")
      
      # Select the correct methylset based on the normalization method
      current_methylset <- normalisation_results()[[norm_method()]]
      if (is.null(current_methylset)) {
        output$status_message <- renderText(paste0("Error: Methylset for '", norm_method(), "' not found."))
        return()
      }
      
      tryCatch({
        calculated_mean_diffs <- calculateMethylationMeanDifferences(
          normalised_methylset = current_methylset,
          pheno_table = pheno(),
          ref_group = ref_group(),
          round_digits = 5
        )
        combined_m_beta_tbl(calculated_mean_diffs)
        output$status_message <- renderText("Step 1/3: Mean differences calculated successfully.")
      }, error = function(e) {
        output$status_message <- renderText(paste0("Error in Step 1: ", e$message))
        return()
      })
      
      req(combined_m_beta_tbl())
      
      # --- Step 2: Run DMP Finder ---
      incProgress(0.4, detail = "Running dmpFinder...")
      output$status_message <- renderText("Step 2/3: Running dmpFinder...")
      
      tryCatch({
        dmp_raw <- minfi::dmpFinder(
          dat = getM(current_methylset),
          pheno = pheno()$Sample_Group,
          type = "categorical",
          qCutoff = input$qvalue_cutoff,
          shrinkVar = FALSE
        )
        output$status_message <- renderText("Step 2/3: dmpFinder completed.")
      }, error = function(e) {
        output$status_message <- renderText(paste0("Error in Step 2: ", e$message))
        return()
      })
      
      req(dmp_raw)
      
      # --- Step 3: Add Mean Differences to DMP table ---
      incProgress(0.2, detail = "Merging mean differences with DMP table...")
      output$status_message <- renderText("Step 3/3: Merging results...")
      
      tryCatch({
        final_dmp <- addMeanDifferencesToDMP(
          dmp_table = dmp_raw,
          combined_Beta_M_mean_table = combined_m_beta_tbl()
        )
        dmp_final_table(final_dmp)
        output$status_message <- renderText("DMP analysis completed successfully!")
      }, error = function(e) {
        output$status_message <- renderText(paste0("Error in Step 3: ", e$message))
        return()
      })
      
      incProgress(0.2, detail = "Rendering table and preparing downloads...")
    })
  })
  
  # Render the DMP table
  output$dmp_table <- DT::renderDataTable({
    req(dmp_final_table())
    DT::datatable(dmp_final_table(),
                  options = list(
                    pageLength = 10,
                    searching = TRUE,
                    scrollX = TRUE
                  ),
                  rownames = TRUE
    )
  })
  
  # Download UI
  output$download_ui <- renderUI({
    req(dmp_final_table())
    tagList(
      h5(strong("Download DMP Table:")),
      radioButtons("download_format", "Select format:",
                   choices = c("CSV" = "csv", "Excel" = "xlsx"),
                   selected = "csv", inline = TRUE),
      downloadButton("download_dmp", "Download Table",
                     icon = icon("download"))
    )
  })
  
  # Download handler
  output$download_dmp <- downloadHandler(
    filename = function() {
      paste0("dmp_results_", Sys.Date(), ".", input$download_format)
    },
    content = function(file) {
      if (input$download_format == "csv") {
        write.csv(dmp_final_table(), file, row.names = TRUE)
      } else if (input$download_format == "xlsx") {
        if (!requireNamespace("writexl", quietly = TRUE)) {
          stop("Package 'writexl' is required for Excel download.")
        }
        writexl::write_xlsx(as.data.frame(dmp_final_table()), file, row.names = TRUE)
      }
    }
  )
  
  # Notes modal
  observeEvent(input$notes_link, {
    showModal(modalDialog(
      title = "Notes on Reading the DMP Table Columns",
      HTML(
        "<ul>
          <li><strong>intercept:</strong> Baseline methylation level for the reference group.</li>
          <li><strong>f:</strong> F-statistic value (between-group vs within-group variability).</li>
          <li><strong>pval:</strong> Raw p-value (probability of difference by chance).</li>
          <li><strong>qval:</strong> Adjusted p-value (FDR corrected).</li>
          <li><strong>M_mean_difference:</strong> Difference in mean M-values.</li>
          <li><strong>Beta_mean_difference:</strong> Log2 ratio of mean Beta values.</li>
        </ul>"
      ),
      footer = modalButton("Close")
    ))
  })
}

# Run the application
shinyApp(ui = ui, server = server)