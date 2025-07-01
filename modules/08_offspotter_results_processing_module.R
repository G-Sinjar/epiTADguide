# offtargets_import_module.R


'library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("mutate", "dplyr")
'

offtargetsUI <- function(id) {
  ns <- NS(id)
  
  layout_sidebar(
    fillable = TRUE,
    sidebar = sidebar(
      h4("Input Data"),
      textInput(ns("directory_path"), "Directory Path:", placeholder = "e.g. C:\\path\\to\\data"),
      textAreaInput(ns("guide_file_names"), "Guide File Names (one per line):",
                    placeholder = "e.g.,\nguide1.txt\nguide2a.txt", rows = 5),
      p(tags$small("Note: Ensure your guide files are tab-separated and properly named.")),
      actionButton(ns("run_processing"), "Process Guide Files into one file", class = "btn-primary"),
      hr(),
      conditionalPanel(
        condition = paste0("output['", ns("tableReady"), "'] === true"),
        h4("Download Results"),
        selectInput(ns("download_format"), "Select Format:", choices = c("CSV", "Excel")),
        downloadButton(ns("download_data"), "Download", class = "btn-success")
      )
    ),
    
    div(
      style = "padding-left: 15px; padding-right: 15px;",
      layout_columns(
        col_widths = c(12),
        card(
          card_header("Off-targets File Processing Status"),
          card_body(uiOutput(ns("status_output_card_body")))
        ),
        div(
          h5("Combined all guides offtargets"),
          uiOutput(ns("row_count_text")),
          DTOutput(ns("offtargets_table"), width = "100%")
        )
      )
    )
  )
}

offtargetsServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    output_messages <- reactiveVal(list())
    processed_data <- reactiveVal(NULL)
    
    observeEvent(input$run_processing, {
      output_messages(list())
      processed_data(NULL)
      
      msgs <- list()
      all_guides <- list()
      
      dir_path <- gsub("\\\\", "/", input$directory_path)
      msgs <- append(msgs, paste0("🛠️ Processing started... Directory: ", dir_path))
      
      if (!dir.exists(dir_path)) {
        msgs <- append(msgs, "❌ Error: The path you provided doesn't exist.")
        output_messages(msgs)
        return()
      }
      
      file_names <- strsplit(input$guide_file_names, "\\s+")[[1]]
      
      col_names <- c("Chrom", "Strand", "Start", "End", "Given query", "Actual genomic hit",
                     "Number_of_mismatches", "Pre-mRNA (Unspliced)", "mRNA (5UTR)",
                     "mRNA (CDS)", "mRAN (3UTR)", "lincRNA (Unspliced)",
                     "lincRNA (Spliced)", "GC content")
      
      for (file_name in file_names) {
        file_path <- file.path(dir_path, file_name)
        msgs <- append(msgs, paste0("🔍 Reading file: ", file_path))
        
        guide_data <- tryCatch({
          read_delim(file_path, delim = "\t", col_names = FALSE, show_col_types = FALSE)
        }, error = function(e) {
          msgs <<- append(msgs, paste0("❌ Step 1 - Reading failed: ", e$message))
          return(NULL)
        })
        if (is.null(guide_data)) next
        msgs <- append(msgs, "✅ Step 1 - File read successfully")
        
        if (ncol(guide_data) != length(col_names)) {
          msgs <- append(msgs, paste0("❌ Step 2 - Column mismatch: expected ", length(col_names), ", got ", ncol(guide_data)))
          next
        }
        msgs <- append(msgs, "✅ Step 2 - Column count verified")
        
        processed <- tryCatch({
          colnames(guide_data) <- col_names
          guide_data <- guide_data[, !colnames(guide_data) %in% "GC content"]
          guide_data$guide <- tools::file_path_sans_ext(file_name)
          guide_data$Start <- as.integer(guide_data$Start)
          guide_data$End <- as.integer(guide_data$End)
          guide_data$Start_minus_70 <- guide_data$Start - 70
          guide_data$End_plus_70 <- guide_data$End + 70
          guide_data$Number_of_mismatches <- as.integer(guide_data$Number_of_mismatches)
          guide_data
        }, error = function(e) {
          msgs <<- append(msgs, paste0("❌ Step 3 - Processing failed: ", e$message))
          return(NULL)
        })
        if (is.null(processed)) next
        msgs <- append(msgs, "✅ Step 3 - Data processed")
        
        all_guides[[file_name]] <- processed
        msgs <- append(msgs, "✅ Step 4 - Data in the guide file stored")
      }
      
      msgs <- append(msgs, "📦 Step 5 - All off-targets in one file...")
      combined <- tryCatch({
        bind_rows(all_guides)
      }, error = function(e) {
        msgs <<- append(msgs, paste0("❌ Step 5a - Combining failed: ", e$message))
        return(NULL)
      })
      if (is.null(combined)) {
        output_messages(msgs)
        return()
      }
      msgs <- append(msgs, "✅ Step 5a - Guides combined")
      
      final <- tryCatch({
        combined %>%
          dplyr::filter(Number_of_mismatches != 0) %>%
          dplyr::mutate(
            guide_num = stringr::str_match(guide, "(?i)guide([[:alnum:]]+)")[, 2],
            ID = paste0("G", guide_num, ".M", Number_of_mismatches, ".chr", Chrom, "_", Start)
          ) %>%
          dplyr::select(-guide_num) %>%
          dplyr::select(ID, guide, Chrom, Start, End, Number_of_mismatches,
                        Start_minus_70, End_plus_70, dplyr::everything())
      }, error = function(e) {
        msgs <<- append(msgs, paste0("❌ Step 5b - Filtering failed: ", e$message))
        return(NULL)
      })
      
      
      if (is.null(final)) {
        output_messages(msgs)
        return()
      }
      msgs <- append(msgs, "✅ Step 5b - Filtered out on-targets")
      
      # ✅ Step 6: Save as RDS in ./intermediate_data ------------------------
      output_dir <- "./intermediate_data"
      if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
      
      # Construct safe filename from guide names
      safe_filename <- gsub("[^A-Za-z0-9_]", "_", paste0(tools::file_path_sans_ext(file_names), collapse = "_"))
      output_path <- file.path(output_dir, paste0(safe_filename, ".rds"))
      
      tryCatch({
        saveRDS(final, output_path)
        msgs <- append(msgs, paste0("💾 Step 6 - Saved processed data as RDS at: ", output_path))
      }, error = function(e) {
        msgs <<- append(msgs, paste0("❌ Step 6 - Saving RDS failed: ", e$message))
      })
      
      processed_data(final)
      output_messages(msgs)
    })
    
    output$status_output_card_body <- renderUI({
      tags$ul(
        style = "padding-left: 1.2em; margin: 0;",
        lapply(output_messages(), function(msg) {
          tags$li(style = "margin-bottom: 4px;", msg)
        })
      )
    })
    
    output$offtargets_table <- renderDT({
      req(processed_data())
      processed_data()
    }, options = list(pageLength = 25))
    
    output$row_count_text <- renderUI({
      req(processed_data())
      n <- nrow(processed_data())
      HTML(paste0("<p><strong>", n, "</strong> off-targets found</p>"))
    })

    return(processed_data)
  })
}

# app.R

library(shiny)
library(DT)
library(readr)
library(dplyr)
library(stringr)


ui <- page_navbar(
  title = "EPIC Array Pipeline",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  nav_panel("Offtargets Import Test",
            offtargetsUI("myOfftargetModule"))
)

server <- function(input, output, session) {
  # Call the module
  processed_data <- offtargetsServer("myOfftargetModule")
  
  # Optional: observe or debug what’s returned from the module
  observe({
    req(processed_data())
    cat("✅ Test App: Processed data returned with", nrow(processed_data()), "rows\n")
  })
}

shinyApp(ui = ui, server = server)
