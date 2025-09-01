#05_Annotation_Module.R
#Author: Ghazal Sinjar
#Date: 07.06.2025
#Description:
#This Shiny module provides an interactive UI and server logic
#for annotating Illumina EPIC methylation array data. It combines
#array probe annotation with beta values and offers an option
#to download the annotated table in CSV or Excel format.


# ─────────────────────────────────────
# User Interface (UI) FUNCTION
# ─────────────────────────────────────
annotationUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    shinyjs::useShinyjs(),
    
    page_sidebar(
      sidebar = sidebar(
        width = "300px",
        actionButton(ns("run_annot"), "Run Annotation on filtered data"),
        
        # Toggle notes section
        div(
          style = "margin-top: 20px;",
          actionLink(ns("toggle_notes"), "📘 Show notes for reading the table"),
          div(
            id = ns("note_section"),
            style = "display: none; margin-top: 10px; font-size: 0.9em;",
            
            h5("Notes for reading the table:"),
            p(strong("ProbeSeqA, ProbeSeqB:"), " the sequence of the probe; B is for Infinium I beadtypes which has 2 probes instead of one."),
            p(strong("Type:"), " Infinium Design Type; Infinium I (2 probes/locus) or Infinium II (1 probe/locus)."),
            p(strong("Next_Base:"), " For Infinium I probes, the nucleotide immediately following the CpG. Blank for Infinium II."),
            p(strong("Color:"), " For Infinium I probes, the color channel of the 'Next_Base' signal."),
            p(strong("Strand_CO:"), " Refers to whether the probe targets the converted or opposite strand."),
            p(strong("Probe_rs:"), " SNP(s) within the probe sequence itself."),
            p(strong("CpG_rs:"), " SNP(s) at the CpG site. Blank if removed."),
            p(strong("SBE_rs:"), " SNP(s) at the Single Base Extension (SBE) site. Blank if removed."),
            p(strong("...maf:"), " Minor Allele Frequency of the SNP in the population."),
            p(strong("Island_Name:"), " CpG Island coordinates from UCSC."),
            p(strong("Relation_to_Island:"), " CpG location relative to island:"),
            HTML("<ul>
              <li>Shore = 0-2 kb from island.</li>
              <li>Shelf = 2-4 kb from island.</li>
              <li>N = upstream (5') of island.</li>
              <li>S = downstream (3') of island.</li>
            </ul>"),
            p(strong("Probe_Type:"), " Probe type: cg=CpG, nv=variant, rs=SNP, ch=Cp<nonG base>"),
            p(strong("UCSC_RefGene:"), " Gene annotations from UCSC."),
            p(strong("Gencode:"), " Gene annotations from GENCODE."),
            p(strong("Regulatory_Feature_Group:"), " Predicted regulatory elements."),
            p(strong("Methyl450K_Enhancer:"), " Marked 'True' for 450K enhancer annotations."),
            p(strong("Methyl...._Loci:"), " CpGs from older platforms."),
            p("Final columns show ", strong("beta values"), " per probe/sample.")
          )
        ),
        
        hr(),
        radioButtons(ns("download_format"), "Choose format:", choices = c("CSV", "Excel"), inline = TRUE),
        downloadButton(ns("download_data"), "Download Table")
      ),
      DTOutput(ns("annot_table"))
    )
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
    #------------------------------------------------------
    observeEvent(input$toggle_notes, {
      shinyjs::toggle(id = "note_section")  # Toggle visibility of the note section
    })
    
    #-----------------------------------------
    observeEvent(input$run_annot, {
      req(grset_reactive())

      current_grset <- grset_reactive()
      print(paste("Class of object passed to getAnnotation:", class(current_grset)))
      
      if (is.null(current_grset) || !any(class(current_grset) %in% c("GenomicRatioSet", "MethylSet", "RGChannelSet"))) {
        showNotification("Error: Data is not a valid minfi object (GenomicRatioSet, MethylSet, or RGChannelSet).", type = "error", duration = NULL)
        return(NULL)
      }
      
      #--------------------------------
      withProgress(message = "Running annotation steps...", value = 0, {
        incProgress(0.2, detail = "Step 1: Getting annotation...")
        annotation <- getAnnotation(current_grset)
        annotation_obj(annotation)
        annotation_df <- as.data.frame(annotation)
        # Specify columns to remove
        cols_to_remove <- c(
          "AddressA", "AddressB", "Strand_FR", "Strand_TB", "col", "Infinium_Design", "Regulatory_Feature_Group","Rep_Num", 
          "Manifest_probe_match", "Phantom5_Enhancers", "HMM_Island", "DMR")
        
        # Drop them if present
        annotation_df <- annotation_df[, !(colnames(annotation_df) %in% cols_to_remove)]
        #--------------------------------------------------------------
        incProgress(0.4, detail = "Step 2: Extracting beta values...")
        beta <- round(getBeta(current_grset), 5)
        
        #------------------------------------------
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
    #------------------
    #rendering the table to UI
    output$annot_table <- renderDT({
      req(annot_df())
      datatable(annot_df(), options = list(pageLength = 25, scrollX = TRUE))
    })
    #---------------
    # download button handler
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
library(DT)
library(minfi)
library(writexl)
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
  path_rc <- reactive({"./main_app_tests"})
  
  # Call the module
  annotated_result <- annotationServer("annot", 
                                       grset_reactive= grset_reactive,
                                       project_output_dir=path_rc)
}

shinyApp(ui, server)'