#05_Annotation_Module.R
#Author: Ghazal Sinjar
#Date: 07.06.2025
#Description:
#This Shiny module provides an interactive UI and server logic
#for annotating Illumina EPIC methylation array data. It combines
#array probe annotation with beta values and offers an option
#to download the annotated table in CSV or Excel format.


library(shiny)
library(DT)
library(minfi)
library(writexl)

# ─────────────────────────────────────
# User Interface (UI) FUNCTION
# ─────────────────────────────────────
annotationUI <- function(id) {
  ns <- NS(id)
  page_sidebar(
    sidebar = sidebar(
      actionButton(ns("run_annot"), "Run Annotation on filtered data"),
      hr(),
      radioButtons(ns("download_format"), "Choose format:", choices = c("CSV", "Excel"), inline = TRUE),
      downloadButton(ns("download_data"), "Download Table")
    ),
    DTOutput(ns("annot_table"))
  )
}



# ─────────────────────────────────────
# SERVER FUNCTION 
# ─────────────────────────────────────
#' @param id Module ID.
#' @param grset_reactive A reactive expression holding the filtered GenomicRatioSet object.
#' @param project_output_dir A reactive expression for the project's output directory.
#' @return A list of reactive values: `annotation_object` (the minfi annotation object),
#'         and `annotated_table` (the combined annotation and beta values table).
annotationServer <- function(id, grset_reactive, project_output_dir) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    annotation_obj <- reactiveVal(NULL)
    annot_df <- reactiveVal(NULL)
    temp_file <- reactiveVal(NULL)
    
    observeEvent(input$run_annot, {
      req(grset_reactive())
      #req(project_output_dir())
      
      current_grset <- grset_reactive()
      print(paste("Class of object passed to getAnnotation:", class(current_grset)))
      
      if (is.null(current_grset) || !any(class(current_grset) %in% c("GenomicRatioSet", "MethylSet", "RGChannelSet"))) {
        showNotification("Error: Data is not a valid minfi object (GenomicRatioSet, MethylSet, or RGChannelSet).", type = "error", duration = NULL)
        return(NULL)
      }
      
      withProgress(message = "Running annotation steps...", value = 0, {
        incProgress(0.2, detail = "Step 1: Getting annotation...")
        annotation <- getAnnotation(current_grset)
        annotation_obj(annotation)
        annotation_df <- as.data.frame(annotation)
        
        incProgress(0.4, detail = "Step 2: Extracting beta values...")
        beta <- getBeta(current_grset)
        
        incProgress(0.6, detail = "Step 3: Combining data...")
        if (identical(rownames(annotation_df), rownames(beta))) {
          combined_df <- cbind(annotation_df, beta)
          annot_df(combined_df)
          
          tmp_file <- tempfile(fileext = ifelse(input$download_format == "CSV", ".csv", ".xlsx"))
          if (input$download_format == "CSV") {
            write.csv(combined_df, tmp_file, row.names = FALSE)
          } else {
            write_xlsx(combined_df, path = tmp_file)
          }
          temp_file(tmp_file)
          
          # Step 4: Save annotation object to RDS
          incProgress(0.8, detail = "Step 4: Saving annotation object...")
          output_dir_full <- file.path(project_output_dir(), "intermediate_data")
          if (!dir.exists(output_dir_full)) {
            dir.create(output_dir_full, recursive = TRUE)
          }
          saved_path <- file.path(output_dir_full, paste0("annotated_object_", format(Sys.Date(), "%Y%m%d"), ".rds"))
          
          # Save full annotation results as list
          annotation_result_list <- list(
            annotation_object = annotation,
            annotated_table = combined_df
          )
          saveRDS(annotation_result_list, file = saved_path)
          
          
          showNotification(paste("✅ Annotation object automatically saved to:", saved_path),
                           type = "message", duration = 10)
        } else {
          showNotification("Row names do not match between annotation and beta values.", type = "error")
          annot_df(NULL)
          temp_file(NULL)
        }
        
        incProgress(1)
      })
    })
    
    output$annot_table <- renderDT({
      req(annot_df())
      datatable(annot_df(), options = list(pageLength = 25, scrollX = TRUE))
    })
    
    output$download_data <- downloadHandler(
      filename = function() {
        paste0("annotation.", ifelse(input$download_format == "CSV", "csv", "xlsx"))
      },
      content = function(file) {
        req(temp_file())
        file.copy(temp_file(), file)
      }
    )
    # Return only the annotation object
    return(list(
        annotation_object = annotation_obj,
        annotated_table = annot_df
      ))
  })
}


'# test
library(shiny)
library(bslib)

# Load your preprocessed GRanges object
ratio_geno_Swan_NoSNP <- readRDS("../intermediate_data/filtered_GRset_SWAN_SNP_removed_SexChr_kept_20250605.rds")

# UI
ui <- page_navbar(
  title = "EPIC Array Pipeline",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  
  nav_panel("Annotation",
            annotationUI("annot")
  )
)

# Server
server <- function(input, output, session) {
  # Wrap your object as a reactive expression
  grset_reactive <- reactive({ ratio_geno_Swan_NoSNP })
  path_rc <- reactive({"./epic-test"})
  
  # Call the module
  annotated_result <- annotationServer("annot", 
                                       grset_reactive= grset_reactive,
                                       project_output_dir=path_rc)
}

shinyApp(ui, server)'