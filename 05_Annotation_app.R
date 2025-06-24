library(shiny)
library(bslib)
library(DT)
library(openxlsx)
library(minfi)


#ratio_geno_Swan_NoSNP <- readRDS("./intermediate_data/filtered_GRset_SWAN_SNP_removed_SexChr_kept_20250605.rds")


# UI
ui <- page_navbar(
  title = "EPIC Array Pipeline",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  
  nav_panel("Annotation",
            layout_sidebar(
              sidebar = sidebar(
                actionButton("run_annot", "Run Annotation"),
                hr(),
                downloadButton("download_data", "Download Table"),
                radioButtons("download_format", "Choose format:", choices = c("CSV", "Excel"), inline = TRUE)
              ),
              DTOutput("annot_table")
            )
  )
)

# Server
server <- function(input, output, session) {
  annotated_data <- reactiveVal(NULL)
  
  observeEvent(input$run_annot, {
    withProgress(message = "Running annotation steps...", value = 0, {
      
      incProgress(0.3, detail = "Step 1: Getting annotation...")
      annotation <- getAnnotation(ratio_geno_Swan_NoSNP)
      annotation_df <- as.data.frame(annotation)
      
      incProgress(0.6, detail = "Step 2: Extracting beta values...")
      beta <- getBeta(ratio_geno_Swan_NoSNP)
      
      incProgress(0.9, detail = "Combining annotation with beta values...")
      if (identical(rownames(annotation_df), rownames(beta))) {
        annotated_with_betas <- cbind(annotation_df, beta)
        annotated_data(annotated_with_betas)
      } else {
        showNotification("Row names do not match between annotation and beta values.", type = "error")
      }
      
      incProgress(1, detail = "Done.")
    })
  })
  
  
  output$annot_table <- renderDT({
    req(annotated_data())
    datatable(annotated_data(), options = list(pageLength = 10, scrollX = TRUE))
  })
  
  output$download_data <- downloadHandler(
    filename = function() {
      paste0("annotation.", ifelse(input$download_format == "CSV", "csv", "xlsx"))
    },
    content = function(file) {
      data <- annotated_data()
      req(data)
      if (input$download_format == "CSV") {
        write.csv(data, file, row.names = FALSE)
      } else {
        write.xlsx(data, file)
      }
    }
  )
}

# Run App
shinyApp(ui, server)
