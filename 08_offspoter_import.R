# app.R

library(shiny)
library(bslib)
library(readr)
library(dplyr)
library(stringr)
library(tools)
library(GenomicRanges)
library(IRanges)
library(S4Vectors)
library(DT)


read_guide_files <- function(directory_path, guide_file_names) {
  captured_output <- list(messages = character(), warnings = character(), errors = character())
  append_output <- function(type, msg) {
    captured_output[[type]] <<- c(captured_output[[type]], msg)
  }
  
  col_names <- c("Chrom", "Strand", "Start", "End", "Given query", "Actual genomic hit",
                 "Number_of_mismatches", "Pre-mRNA (Unspliced)", "mRNA (5UTR)",
                 "mRNA (CDS)", "mRAN (3UTR)", "lincRNA (Unspliced)",
                 "lincRNA (Spliced)", "GC content")
  
  all_guides <- list()
  
  for (file_name in guide_file_names) {
    full_file_path <- file.path(directory_path, file_name)
    append_output("messages", paste("Processing file:", full_file_path))
    
    guide_data <- tryCatch(
      {
        read_delim(full_file_path, delim = "\t", col_names = FALSE, show_col_types = FALSE)
      },
      error = function(e) {
        err_msg <- paste("Critical Error in reading a file: '", full_file_path, "'. Reason: ", e$message, sep = "")
        append_output("errors", err_msg)
        return(NULL)
      }
    )
    if (is.null(guide_data)) next
    
    if (ncol(guide_data) != length(col_names)) {
      err_msg <- paste("Critical Error: File '", full_file_path, "' has ", ncol(guide_data),
                       " columns, but ", length(col_names), " columns were expected.", sep = "")
      append_output("errors", err_msg)
      next
    }
    
    processed_guide_data <- tryCatch(
      {
        colnames(guide_data) <- col_names
        
        # Remove unwanted column
        guide_data <- guide_data[, !(colnames(guide_data) %in% "GC content")]
        
        # Guide name
        guide_name_for_col <- tools::file_path_sans_ext(file_name)
        guide_data$guide <- guide_name_for_col
        
        # Compute and clean numeric columns
        guide_data$Start <- as.integer(guide_data$Start)
        guide_data$End <- as.integer(guide_data$End)
        guide_data$Start_minus_70 <- guide_data$Start - 70
        guide_data$End_plus_70 <- guide_data$End + 70
        guide_data$Start_minus_70 <- as.integer(guide_data$Start_minus_70)
        guide_data$End_plus_70 <- as.integer(guide_data$End_plus_70)
        guide_data$Number_of_mismatches <- as.integer(guide_data$Number_of_mismatches)
        
        guide_data
      },
      error = function(e) {
        err_msg <- paste("Critical Error: Problem with initial data processing for file '", full_file_path,
                         "'. Check column names or data structure. Reason: ", e$message, sep = "")
        append_output("errors", err_msg)
        return(NULL)
      }
    )
    if (is.null(processed_guide_data)) next
    
    all_guides[[file_name]] <- processed_guide_data
  }
  
  if (length(all_guides) == 0) {
    append_output("errors", "No guide files were successfully processed. Please check paths, file names and contents.")
    return(list(data = NULL, output = captured_output))
  }
  
  combined_guides <- tryCatch(
    {
      bind_rows(all_guides)
    },
    error = function(e) {
      err_msg <- paste("Critical Error: Failed to combine all guide tables. Reason: ", e$message, sep = "")
      append_output("errors", err_msg)
      return(NULL)
    }
  )
  if (is.null(combined_guides)) {
    return(list(data = NULL, output = captured_output))
  }
  
  just_offtargets <- tryCatch(
    {
      combined_guides %>%
        dplyr::filter(`Number_of_mismatches` != 0) %>%
        mutate(
          guide_num = str_extract(guide, "(?i)guide([[:alnum:]]+)", group = 1),
          ID = paste0("G", guide_num, ".M", `Number_of_mismatches`, ".chr", Chrom, "_", Start)
        ) %>%
        dplyr::select(-guide_num)
    },
    
    error = function(e) {
      err_msg <- paste("Critical Error: Failed during off-target filtering or ID generation. Reason: ", e$message, sep = "")
      append_output("errors", err_msg)
      return(NULL)
    }
  )
  if (is.null(just_offtargets)) {
    return(list(data = NULL, output = captured_output))
  }
  
  final_data <- tryCatch(
    {
      just_offtargets %>%
        dplyr::select(
          ID,
          guide,
          Chrom,
          Start,
          End,
          Number_of_mismatches,
          Start_minus_70,
          End_plus_70,
          dplyr::everything()
        )
    },
    
    error = function(e) {
      err_msg <- paste("Critical Error: Failed during final column reordering. Reason: ", e$message, sep = "")
      append_output("errors", err_msg)
      return(NULL)
    }
  )
  if (is.null(final_data)) {
    return(list(data = NULL, output = captured_output))
  }
  
  append_output("messages", "Finished processing all guide files successfully.")
  return(list(data = final_data, output = captured_output))
}








# Shiny UI
ui <- page_navbar(
  title = "EPIC Array Pipeline",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  
  nav_panel(
    "Offtargets import",
    layout_sidebar( fillable = TRUE,
      sidebar = sidebar(
        h4("Input Data"),
        textInput(
          "directory_path",
          "Directory Path:",
          placeholder = "e.g. C:\\path\\to\\data" # Corrected placeholder
        ),
        textAreaInput(
          "guide_file_names",
          "Guide File Names (one per line):",
          placeholder = "e.g.,\nguide1.txt\nguide2a.txt\nguide2b.txt",
          rows = 5
        ),
        p(tags$small(
          "Note: Ensure your guide files are tab-separated, located in the specified directory, and named using only the guide ID (e.g., 'guide1.txt', 'guide2a.txt')."
        )
        ),
        actionButton("read_files_btn", "Read Guide Files", class = "btn-primary"),
        hr(),
        conditionalPanel(
          condition = "output.tableReady === true",
          h4("Download Results"),
          selectInput("download_format", "Select Format:", choices = c("CSV", "Excel")),
          downloadButton("download_data", "Download", class = "btn-success")
        )
        ),
      mainPanel(
          uiOutput("status_output"),
          br(),
          h4("Off-targets Table"),
          textOutput("offtargets_row_count"),
          DTOutput("offtargets_table", width = "100%")
      )
    )
  )
)






# Shiny Server
# Inside your server function in app.R

server <- function(input, output, session) {
  
  process_results_reactive <- eventReactive(input$read_files_btn, {
    # >>> CHANGE START: Early Input Validation for Directory Path <<<
    if (is.null(input$directory_path) || input$directory_path == "") {
      return(list(data = NULL, output = list(
        messages = character(),
        warnings = character(),
        errors = "No directory path provided. Please enter a directory path in the sidebar."
      )))
    }
    # >>> CHANGE END: Early Input Validation for Directory Path <<<
    
    file_names_vector <- unlist(str_split(input$guide_file_names, "\n"))
    file_names_vector <- str_trim(file_names_vector)
    file_names_vector <- file_names_vector[file_names_vector != ""]
    
    
    normalized_directory_path <- str_replace_all(input$directory_path, "\\\\", "/")
    
    captured_messages_from_handlers <- character()
    captured_warnings_from_handlers <- character()
    
    result <- NULL
    
    withCallingHandlers(
      expr = {
        result <- read_guide_files(normalized_directory_path, file_names_vector)
      },
      message = function(m) {
        captured_messages_from_handlers <<- c(captured_messages_from_handlers, m$message)
        invokeRestart("muffleMessage")
      },
      warning = function(w) {
        captured_warnings_from_handlers <<- c(captured_warnings_from_handlers, w$message)
        invokeRestart("muffleWarning")
      }
    )
    
    func_messages <- if (!is.null(result) && !is.null(result$output$messages)) result$output$messages else character()
    func_warnings <- if (!is.null(result) && !is.null(result$output$warnings)) result$output$warnings else character()
    func_errors <- if (!is.null(result) && !is.null(result$output$errors)) result$output$errors else character()
    
    final_output <- list(
      messages = c(captured_messages_from_handlers, func_messages),
      warnings = c(captured_warnings_from_handlers, func_warnings),
      errors = func_errors
    )
    
    list(data = if (!is.null(result)) result$data else NULL, output = final_output)
  }, ignoreNULL = TRUE)
  
  status_card_content <- reactive({
    if (is.null(input$read_files_btn) || input$read_files_btn == 0) {
      list(
        title = "Importing Status",
        color = NULL,
        body = p("Enter the **directory path** and **guide file names** in the sidebar, then click 'Read Guide Files' to start the analysis. The process status (messages, warnings, errors) will appear here.")
      )
    } else {
      results <- process_results_reactive()

      all_output_lines <- c(
        if (length(results$output$messages) > 0) paste(results$output$messages, collapse = "\n"),
        if (length(results$output$warnings) > 0) paste("WARNINGS:\n", paste(results$output$warnings, collapse = "\n")),
        if (length(results$output$errors) > 0) paste("ERRORS:\n", paste(results$output$errors, collapse = "\n"))
      )
      
      status_text <- paste(all_output_lines, collapse = "\n\n")
      if (status_text == "") {
        status_text <- "All steps completed successfully. No messages, warnings, or errors."
      }
      
      # The condition here is already correct to catch errors from process_results_reactive
      if (length(results$output$errors) > 0 || is.null(results$data) || (is.data.frame(results$data) && nrow(results$data) == 0)) {
        card_color <- "bg-danger"
        card_title <- "Processing Failed"
      } else if (length(results$output$warnings) > 0) {
        card_color <- "bg-warning"
        card_title <- "Processing Completed with Warnings"
      } else {
        card_color <- "bg-success"
        card_title <- "Processing Successful"
      }
      
      list(
        title = card_title,
        color = card_color,
        body = pre(status_text, style = "color: white; overflow-x: auto; white-space: pre-wrap; word-wrap: break-word;")
      )
    }
  })
  
  output$status_output <- renderUI({
    content <- status_card_content()
    
    card(
      style = "width: 100%;",
      class = content$color,
      card_header(content$title, class = if (!is.null(content$color)) "text-white" else NULL),
      card_body(content$body)
    )
  })
  
  output$offtargets_row_count <- renderText({
    results <- process_results_reactive()
    
    if (!is.null(results$data) && nrow(results$data) > 0) {
      paste("Number of off-targets found in all files:", nrow(results$data))
    } else {
      # This will display nothing or an empty string if no data,
      # as the table message handles the "no data" scenario.
      "" 
    }
  })
  
  output$offtargets_table <- renderDT({
    if (is.null(input$read_files_btn) || input$read_files_btn == 0) {
      return(NULL)
    }
    
    results <- process_results_reactive()
    
    if (!is.null(results$data) && nrow(results$data) > 0) {
      datatable(
        results$data,
        options = list(
          pageLength = 10,
          autoWidth = TRUE
        ),
        rownames = FALSE,
        filter = "top",
        class = "stripe hover compact"
      )
    } else {
      datatable(
        data.frame(Message = "No off-targets data found or processed. Check the status messages for details."),
        options = list(dom = 't'), # hide search box
        rownames = FALSE
      )
    }
  })
  
  
  output$download_data <- downloadHandler(
    filename = function() {
      ext <- if (input$download_format == "Excel") "xlsx" else "csv"
      paste0("offtargets_", Sys.Date(), ".", ext)
    },
    content = function(file) {
      results <- process_results_reactive()
      if (!is.null(results$data) && nrow(results$data) > 0) {
        if (input$download_format == "Excel") {
          writexl::write_xlsx(results$data, path = file)
        } else {
          write.csv(results$data, file, row.names = FALSE)
        }
      } else {
        writeLines("No data to export.", file)
      }
    }
  )
  
  
  output$tableReady <- reactive({
    results <- process_results_reactive()
    !is.null(results$data) && nrow(results$data) > 0
  })
  outputOptions(output, "tableReady", suspendWhenHidden = FALSE)
  
}



shinyApp(ui = ui, server = server)