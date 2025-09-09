# 07_Boxplot_Module.R
# Author: Ghazal Sinjar
# Date: 23.06.2025
# Shiny module for generating CpG beta value boxplots for DMRs

# ─────────────────────────────────────
# User Interface (UI) FUNCTION
# ─────────────────────────────────────
boxplotUI <- function(id) {
  ns <- NS(id)
  
  layout_sidebar(
    sidebar = sidebar( width = "300px",
      selectInput(ns("dmr"), "Select DMR ID (chr_gene)", choices = NULL),
      switchInput(ns("interactive"), "Interactive Boxplot", value = FALSE),
      switchInput(ns("pos_spacing"), "Use Positional Spacing", value = FALSE),
      helpText("Positional spacing will arrange the boxes according to their position on the chromosome, whereas PositionLabel will just use the CpG names."),
      br(),
      h5("Q-value Cutoff:"),
      verbatimTextOutput(ns("qval_cutoff")),
      helpText("All CpG sites with a q-value less than or equal to the threshold q-value indicated above, their box's border will be highlighted in green."),
      actionButton(ns("create_plot"), "Create Boxplot"),
      hr(),
      downloadButton(ns("download_plot"), "Download Plot")
    ),
    
    div(
      style = "padding-left: 15px; padding-right: 15px;",
      layout_columns(
        col_widths = c(12),
        card(
          card_header("Boxplot Creating Status"),
          card_body(verbatimTextOutput(ns("boxplot_status")))
        ),
        uiOutput(ns("plot_area"))
      )
    )
  )
}

# ─────────────────────────────────────
# SERVER FUNCTION
# ─────────────────────────────────────
boxplotServer <- function(id, dmr_output_reactive, annotation_tbl_reactive, dmp_results, project_output_dir) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # --- Access reactive values ---
    dmr_table_r <- reactive({
      req(dmr_output_reactive()$dmr_table())
      dmr_output_reactive()$dmr_table()
    })
    
    pheno_data_r <- reactive({
      req(dmr_output_reactive()$pheno())
      dmr_output_reactive()$pheno()
    })
    
    ref_group_r <- reactive({
      req(dmr_output_reactive()$ref_group())
      dmr_output_reactive()$ref_group()
    })
    
    annotated_with_betas_df_r <- reactive({
      req(annotation_tbl_reactive()$annotated_table())
      annotation_tbl_reactive()$annotated_table()
    })
    dmp_table_r <- reactive({
      req(dmp_results()$DMP_tbl())
      dmp_results()$DMP_tbl()
    })
    
    last_qvalue_r <- reactive({
      req(dmp_results()$last_qvalue())
      dmp_results()$last_qvalue()
    })
    #----------------------------------------------------------------
    # --- RENDER Q-VALUE CUTOFF TEXT ---
    output$qval_cutoff <- renderText({
      req(last_qvalue_r())
      # Format the q-value with exactly 4 decimal places
      paste0(
        "q <= ",
        sprintf("%.5f", as.numeric(last_qvalue_r()))
      )
    })
    #---------------------------------------------
    # update choices names with chr and gene if available
    observe({
      req(dmr_table_r())
      dmr_data <- dmr_table_r()
      display_names <- ifelse(
        is.na(dmr_data$first_overlapped_gene),
        paste0(dmr_data$DMR_ID, " (", dmr_data$chr, ")"),
        paste0(dmr_data$DMR_ID, " (", dmr_data$chr, "_", dmr_data$first_overlapped_gene, ")")
      )
      choices_list <- setNames(dmr_data$DMR_ID, display_names)
      updateSelectInput(session, "dmr", choices = choices_list)
    })
    
    #----------------------------------
    status_message <- reactiveVal("ℹ️ Waiting to select a DMR ID.")
    output$boxplot_status <- renderText({ status_message() })
    #------------------------------------------------------
    dmp_table_r <- reactive({
      # More robust access to DMP table
      dmp_data <- dmp_results()
      if (!is.null(dmp_data$DMP_tbl)) {
        tbl <- dmp_data$DMP_tbl()
        if (is.null(tbl) || nrow(tbl) == 0) {
          return(NULL)
        }
        tbl
      } else {
        NULL
      }
    })
    #-----------------------------------------------------------
    annotated_with_qval_r <- reactive({
      req(annotated_with_betas_df_r())
      
      # Debug 1: Print initial state
      print("DEBUG: Starting annotated_with_qval_r reactive")
      print(paste("DEBUG: annotated_with_betas_df_r() class:", class(annotated_with_betas_df_r())))
      print(paste("DEBUG: annotated_with_betas_df_r() dim:", 
                  if(is.data.frame(annotated_with_betas_df_r())) 
                    paste(dim(annotated_with_betas_df_r()), collapse="x") 
                  else "Not a data frame"))
      
      status_lines <- c("Merging q-values with annotation table...")
      status_message(paste(status_lines, collapse = "\n"))
      
      tryCatch({
        # Debug 2: Print inputs
        annot_table <- annotated_with_betas_df_r()
        dmp_table <- dmp_table_r()
        last_q_val <- last_qvalue_r()
        
        print("DEBUG: Input values:")
        print(paste("dmp_table_r() class:", class(dmp_table)))
        print(paste("dmp_table_r() length:", length(dmp_table)))
        print(paste("last_qvalue_r():", last_q_val))
        
        # Input validation
        validate(
          need(is.data.frame(annot_table), "Annotation table is not a data frame"),
          need(!is.null(last_q_val), "Q-value threshold is NULL"),
          need(!is.na(as.numeric(last_q_val)), "Q-value threshold is not valid")
        )
        
        if (is.null(dmp_table) || nrow(dmp_table) == 0) {
          # Debug 3: No DMPs case
          print("DEBUG: No DMP table or empty DMP table")
          
          status_lines <- c(
            status_lines,
            paste0("⚠️ No significant DMPs found at q <= ", sprintf("%.5f", as.numeric(last_q_val))),
            "✅ Proceeding with annotation data only",
            "🎯 App ready to create boxplots",
            "👉 Please choose a DMR and click 'Create Boxplot' button to proceed"
          )
          status_message(paste(status_lines, collapse = "\n"))
          
          # Debug 4: Print return value
          print("DEBUG: Returning annotation table without q-values")
          print(str(annot_table))
          
          return(annot_table) 
        } else {
          # Debug 5: DMPs available case
          print("DEBUG: Processing DMP table")
          print(paste("DMP table dimensions:", nrow(dmp_table), "x", ncol(dmp_table)))
          print("First few rows of DMP table:")
          print(head(dmp_table))
          
          validate(need(is.data.frame(dmp_table), "DMP table is not a data frame"))
          
          last_q_val_num <- as.numeric(last_q_val)
          qval_vector <- dmp_table$qval
          names(qval_vector) <- rownames(dmp_table)
          
          # Debug 6: Check q-value vector
          print("DEBUG: Q-value vector summary:")
          print(summary(qval_vector))
          
          annot_table$qval <- qval_vector[rownames(annot_table)]
          annot_table$significance_last_qvalue <- ifelse(
            is.na(annot_table$qval) | annot_table$qval > last_q_val_num,
            "Insig",
            "Sig"
          )
          
          # Debug 7: Check merged result
          print("DEBUG: Merged table summary:")
          print(summary(annot_table$qval))
          print(paste("Significant probes:", sum(annot_table$significance_last_qvalue == "Sig")))
          
          status_lines <- c(
            status_lines,
            "✅ Successfully merged q-values",
            "🎯 App ready to create boxplots",
            "👉 Please choose a DMR and click 'Create Boxplot' button to proceed"
          )
          status_message(paste(status_lines, collapse = "\n"))
          
          # Debug 8: Final output
          print("DEBUG: Final merged table structure:")
          print(str(annot_table))
          
          return(annot_table)
        }
        
      }, error = function(e) {
        # Debug 9: Error case
        print(paste("DEBUG: Error in annotated_with_qval_r:", e$message))
        error_msg <- paste0("❌ Merging failed: ", e$message)
        status_lines <<- c(status_lines, error_msg)
        status_message(paste(status_lines, collapse = "\n"))
        return(NULL)
      })
    })
    #------------------------------------------
    # optional: saving the annotated_with_qval_r table as rds
    observe({
      req(annotated_with_qval_r()) 
      
      withProgress(message = 'Saving data', value = 0, {
        # Get the current annotated data
        incProgress(0.2, detail = "Preparing data...")
        annotated_data <- annotated_with_qval_r()
        
        # Only proceed if we have valid data
        if (!is.null(annotated_data)) {
          # --- ADDED CONDITION ---
          # Check for required columns before saving
          required_cols <- c("qval", "significance_last_qvalue")
          if (!all(required_cols %in% names(annotated_data))) {
            warning("Skipping save: no significant CpGs from the DMP analysis run previously")
            showNotification(
              "Save qval skipped: no significant CpGs from the DMP analysis run previously",
              type = "warning",
              duration = 5
            )
            return() # Exit the observe block
          }
          # --- END ADDED CONDITION ---
          
          tryCatch({
            incProgress(0.3, detail = "Creating directory...")
            intermediate_dir <- file.path(project_output_dir(), "intermediate_data")
            if (!dir.exists(intermediate_dir)) {
              dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)
            }
            
            incProgress(0.6, detail = "Saving annotation table with qvalue file...Please wait to get the DMRs!")
            save_path <- file.path(intermediate_dir, "annotated_with_qval.rds")
            saveRDS(annotated_data, save_path)
            
            incProgress(0.9, detail = "Finalizing...")
            # Show message in console
            message(paste("💾 Saved annotated data to:", save_path))
            
          }, error = function(e) {
            warning(paste("Failed to save annotated data:", e$message))
            # Show error notification
            showNotification(
              paste("Save failed:", e$message),
              type = "error",
              duration = NULL  # Persistent until dismissed
            )
          })
        }
        incProgress(1, detail = "Done!")
      })
    })
    #-----------------------------------------------------------------------
    observeEvent(input$create_plot, {
      req(input$dmr, dmr_table_r(), annotated_with_qval_r(), pheno_data_r())
      
      output$plot_area <- renderUI({ NULL })
      
      DMRx <- input$dmr
      interactive_flag <- input$interactive
      use_spacing <- input$pos_spacing
      ref_group <- ref_group_r()
      
      withProgress(message = 'Creating boxplot', value = 0, {
        #--------------------------------------
        # Step 1: Extracting CpGs in the DMR
        incProgress(0.1, detail = "Extracting CpGs in the DMR...")
        status_lines <- c("🔁 Step 1: Extracting CpGs in the DMR...")
        
        region_cpgs <- tryCatch(
          {
            extract_cpgs_in_DMR(DMRx, dmr_table_r(), annotated_with_qval_r(), pheno_data_r())
          },
          error = function(e) {
            status_lines <<- c(status_lines, paste0("❌ Error in Step 1: ", e$message))
            status_message(paste(status_lines, collapse = "\n"))
            output$plot_area <- renderUI({ NULL })
            return(NULL)
          }
        )
        
        if (is.null(region_cpgs)) {
          output$plot_area <- renderUI({ NULL })
          return()
        }
        
        if (is.data.frame(region_cpgs)) {
          status_lines <- c(status_lines, "✅ CpG table created.")
          incProgress(0.3, detail = "CpGs extracted successfully")
        } else {
          status_lines <- c(status_lines, paste0("❌ ", region_cpgs))
          status_message(paste(status_lines, collapse = "\n"))
          output$plot_area <- renderUI({ NULL })
          return()
        }
        
        #-------------------------------------------
        # Step 2: Reshaping data
        incProgress(0.1, detail = "Reshaping to long format...")
        status_lines <- c(status_lines, "🔁 Step 2: Reshaping to long format...")
        
        long_table <- tryCatch(
          {
            reshape_to_long_beta(region_cpgs, pheno_data_r()) 
          },
          error = function(e) {
            status_lines <<- c(status_lines, paste0("❌ Error in Step 2: ", e$message))
            status_message(paste(status_lines, collapse = "\n"))
            output$plot_area <- renderUI({ NULL })
            return(NULL)
          }
        )
        
        if (is.null(long_table)) {
          output$plot_area <- renderUI({ NULL })
          return()
        }
        
        status_lines <- c(status_lines, "✅ Table reshaped.")
        incProgress(0.3, detail = "Data reshaped successfully")
        
        # Step 3: Creating plot
        incProgress(0.1, detail = "Creating boxplot...")
        status_lines <- c(status_lines, "🔁 Step 3: Creating boxplot...")
        
        plot_result <- tryCatch(
          {
            create_boxplot(long_table, interactive = interactive_flag, use_positional_spacing = use_spacing, ref_group = ref_group)
          },
          error = function(e) {
            status_lines <<- c(status_lines, paste0("❌ Error in Step 3: ", e$message))
            status_message(paste(status_lines, collapse = "\n"))
            output$plot_area <- renderUI({ NULL })
            return(NULL)
          }
        )
        
        if (is.null(plot_result)) {
          output$plot_area <- renderUI({ NULL })
          return()
        }
        
        status_lines <- c(status_lines, "✅ Plot successfully created.")
        incProgress(0.2, detail = "Boxplot created")
        
        output$plot_area <- renderUI({
          if (interactive_flag && inherits(plot_result, "plotly")) {
            plotlyOutput(ns("interactive_boxplot"), width = "100%", height = "80vh")
          } else if (!interactive_flag && inherits(plot_result, "ggplot")) {
            plotOutput(ns("boxplot_output"), width = "100%", height = "80vh")
          } else {
            NULL
          }
        })
        #-----------------------
        output$boxplot_output <- renderPlot({
          if (!interactive_flag && inherits(plot_result, "ggplot")) plot_result else NULL
        })
        #-----------------
        output$interactive_boxplot <- renderPlotly({
          if (interactive_flag && inherits(plot_result, "plotly")) plot_result else NULL
        })
        #---------------------------
        output$download_plot <- downloadHandler(
          filename = function() {
            if (interactive_flag) paste0(DMRx, "_boxplot.html") else paste0(DMRx, "_boxplot.pdf")
          },
          content = function(file) {
            if (interactive_flag) {
              if (inherits(plot_result, "ggplot")) {
                htmlwidgets::saveWidget(ggplotly(plot_result), file)
              } else {
                htmlwidgets::saveWidget(plot_result, file)
              }
            } else {
              ggsave(file, plot_result , width = 12, height = 8)
            }
          }
        )
        
        status_message(paste(status_lines, collapse = "\n"))
      }) # End of withProgress
    })
    #return(annotated_with_qval_r())
  })
}

'#test app dont work cause of the error in the return statement eventhough this works in the main app. this test code works only wheen the return done have paranthese return(annotated_with_qval_r) which doesnt work in main app
# Load libraries
library(shiny)
library(bslib)
library(shinyWidgets)
library(plotly)
library(ggplot2)

# Load module and utilities
source("../utils/dmrs_boxplot_utils.R")

# Load RDS objects (as static lists)
results_dmr <- readRDS("../main_app_tests/test/DMR_results/DMRs_unguided_vs_sample_groupguided_cutoff_-0.15_0.15_B99.rds")
results_anno <- readRDS("../main_app_tests/test/intermediate_data/annotated_object_20250831.rds")
dmp_results <- readRDS("./main_app_tests/dmp_qval_test/intermediate_data/DMP_results_SWAN_0.89999qval.rds")

# ---- UI ----
ui <- page_navbar(
  title = "DMR Boxplot Test App",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  nav_panel("Boxplots", boxplotUI("boxplot"))
)

# ---- Server ----
server <- function(input, output, session) {
  # Wrap each element in reactiveVal(), and return a list of reactive functions
  dmr_output_reactive <- reactive({
    list(
      dmr_table = reactiveVal(results_dmr$dmr_table),
      pheno = reactiveVal(results_dmr$pheno),
      ref_group = reactiveVal("unguided")
    )
  })
  
  annotation_tbl_reactive <- reactive({
    list(
      annotated_table = reactiveVal(results_anno$annotated_table)
    )
  })
  
  dmp_results_reactive <- reactive({
    list(
      DMP_tbl = reactiveVal(dmp_results$dmp_table),
      last_qvalue = reactiveVal(dmp_results$last_qvalue)
    )
  })
  #case 2: null
  #dmp_results_reactive <- reactive({
  #  list(
  #    DMP_tbl = reactiveVal(NULL),
  #    last_qvalue = reactiveVal(dmp_results$last_qvalue)
  #  )
  #})
  boxplotServer(
    id = "boxplot",
    dmr_output_reactive = dmr_output_reactive,
    annotation_tbl_reactive = annotation_tbl_reactive,
    dmp_results = dmp_results_reactive,
    project_output_dir = reactive({"../main_app_tests/test"})
  )
}
shinyApp(ui = ui, server = server)'