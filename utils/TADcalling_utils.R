# TADcalling_utils.R
# 2 functions


process_tad_results <- function(chrom, tissue, binsize, input_dir, output_dir) {
  # ✅ Ensure input_dir and output_dir exist
  if (!dir.exists(input_dir)) {
    cat("❌ The specified input directory for TAD results does not exist:\n ➤", normalizePath(input_dir, winslash = "/"), "\n")
    return(invisible(NULL))
  }
  if (!dir.exists(output_dir)) { # Ensure output directory exists before writing
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Define file paths
  input_tad <- file.path(input_dir, sprintf("%s_%s_%skb_TADlines.TAD", tissue, chrom, binsize)) # Use input_dir for reading
  # CHANGE HERE: output_txt now points to input_dir (java_raw_output_dir)
  output_txt <- file.path(input_dir, sprintf("%s_%s_%skb_TADlines_filtered.txt", tissue, chrom, binsize))
  output_tad <- file.path(output_dir, sprintf("%s_%s_%skb_TADs.txt", tissue, chrom, binsize)) # This remains in output_dir
  output_excel <- file.path(output_dir, sprintf("%s_%s_%skb_TADs.xlsx", tissue, chrom, binsize)) # This remains in output_dir
  
  # Check input file
  if (!file.exists(input_tad)) {
    cat("❌ Input TAD file not found:\n ➤", normalizePath(input_tad, winslash = "/"), "\n")
    return(invisible(NULL))
  }
  
  # Initialize output vectors
  tad_lines <- readLines(input_tad)
  filtered_lines <- character()
  start_pos <- c()
  end_pos <- c()
  tad_length <- c()
  
  for (line in tad_lines) {
    bins_str <- strsplit(trimws(line), "\\s+")[[1]]
    bins <- suppressWarnings(as.integer(as.numeric(bins_str)))
    
    if (any(is.na(bins)) || length(bins) <= 1) next
    
    # Store filtered line
    filtered_lines <- c(filtered_lines, line)
    
    # Calculate TAD coordinates
    s <- (bins[1] - 1) * binsize * 1000
    e <- bins[length(bins)] * binsize * 1000
    
    # FIX: Append to start_pos, end_pos, and tad_length vectors
    start_pos <- c(start_pos, s)
    end_pos <- c(end_pos, e)
    tad_length <- c(tad_length, e - s) # Correctly calculate and append tad_length
  }
  
  # Write filtered raw lines to input_dir (java_raw_output_dir)
  writeLines(filtered_lines, output_txt)
  
  # Create TAD coordinate table
  tad_df <- data.frame(
    StartPos = start_pos,
    EndPos = end_pos,
    TadLength = tad_length
  )
  
  # --- NEW LOGIC FOR TAD_ID ---
  # 1. Order the dataframe by StartPos
  tad_df <- tad_df[order(tad_df$StartPos), ]
  
  # 2. Add TAD_ID column at the beginning
  tad_df$TAD_ID <- paste0("TAD", seq_len(nrow(tad_df)))
  
  # 3. Reorder columns to put TAD_ID first
  tad_df <- tad_df[, c("TAD_ID", "StartPos", "EndPos", "TadLength")]
  # --- END NEW LOGIC ---
  
  # Write TAD coordinate table to output_dir (processed_tads_output_dir)
  # Ensure you have 'readr' package loaded for write_tsv, or use write.table
  # library(readr)
  write_tsv(tad_df, output_tad) # Assuming write_tsv is from readr, if not, use write.table as below
  
  cat("📄 Filtered raw lines saved to:\n -", normalizePath(output_txt, winslash = "/"), "\n")
  cat("📄 Processed TADs table saved to:\n -", normalizePath(output_tad, winslash = "/"), "\n")
  
  
  # Write Excel to output_dir (processed_tads_output_dir)
  # Ensure you have 'openxlsx' package loaded for createWorkbook, addWorksheet, writeData, saveWorkbook
  # library(openxlsx)
  wb <- createWorkbook()
  addWorksheet(wb, "Sheet1")
  
  # Header
  writeData(wb, "Sheet1", data.frame(chrom = chrom, binsize = binsize), startRow = 1, colNames = TRUE)
  # Write the sorted tad_df with TAD_ID
  writeData(wb, "Sheet1", tad_df, startRow = 3, colNames = TRUE)
  
  saveWorkbook(wb, output_excel, overwrite = TRUE)
  cat("✅ Excel file created:\n -", normalizePath(output_excel, winslash = "/"), "\n")
  
  return(tad_df)
}



# 4- process SubTAD results 
process_subtad_results <- function(chrom, tissue, binsize, input_dir, output_dir, main_tads_df_input) {
  # ✅ Ensure input_dir and output_dir exist
  if (!dir.exists(input_dir)) {
    cat("❌ The specified input directory for SubTAD results does not exist:\n ➤", normalizePath(input_dir, winslash = "/"), "\n")
    return(invisible(NULL))
  }
  if (!dir.exists(output_dir)) { # Ensure output directory exists before writing
    dir.create(output_dir, recursive = TRUE)
  }
  
  # File paths
  input_window_tad <- file.path(input_dir, sprintf("%s_%s_%skb_TADlines.window.TAD", tissue, chrom, binsize))
  output_window_txt <- file.path(input_dir, sprintf("%s_%s_%skb_TADlines_filtered.window.txt", tissue, chrom, binsize)) # Filtered raw lines go to input_dir
  
  output_subtad_txt <- file.path(output_dir, sprintf("%s_%s_%skb_SubTADs.txt", tissue, chrom, binsize)) # Processed table goes to output_dir
  output_subtad_excel <- file.path(output_dir, sprintf("%s_%s_%skb_SubTADs.xlsx", tissue, chrom, binsize))
  
  
  # ❌ Custom error for missing .window.TAD
  if (!file.exists(input_window_tad)) {
    cat("❌ The expected .window.TAD file was not found.\n")
    cat("📌 Based on your inputs, the missing file is:\n")
    cat("   ➤", normalizePath(input_window_tad, winslash = "/"), "\n")
    cat("📋 Please ensure that the TADcalling step completed successfully for:\n")
    cat("   ➤ Tissue:", tissue, "\n")
    cat("   ➤ Chromosome:", chrom, "\n")
    cat("   ➤ Bin size:", binsize, "kb\n")
    cat("📁 And that the raw Java result is located in the provided input path:\n")
    cat("   ➤", normalizePath(input_dir, winslash = "/"), "\n")
    return(invisible(NULL))
  }
  
  # Use the passed main_tads_df_input directly
  # No need to check for output_tad_main file existence or read it here
  main_tads_df <- main_tads_df_input
  # Ensure main_tads_df is sorted by StartPos for consistent lookup (it should already be from process_tad_results, but good to be sure)
  main_tads_df <- main_tads_df[order(main_tads_df$StartPos), ]
  
  
  # Read input window TADs
  window_lines <- readLines(input_window_tad)
  filtered_lines <- c()
  subtad_start_pos <- c()
  subtad_end_pos <- c()
  subtad_length <- c()
  
  parent_tad_ids <- c()
  sub_tad_sequential_ids <- c()
  
  subtad_counter_per_tad <- list()
  
  
  for (line in window_lines) {
    bins_str <- strsplit(trimws(line), "\\s+")[[1]]
    bins <- suppressWarnings(as.numeric(bins_str))
    bins <- bins[!is.na(bins)]
    
    if (length(bins) <= 1) next
    
    filtered_lines <- c(filtered_lines, line)
    
    s <- (bins[1] - 1) * binsize * 1000
    e <- bins[length(bins)] * binsize * 1000
    tag <- paste(s, e, sep = "_")
    
    # Check if this SubTAD is a main TAD (and thus should be excluded as a true SubTAD)
    is_main_tad <- FALSE
    # Check against the passed main_tads_df_input
    for (i in 1:nrow(main_tads_df)) {
      if (s == main_tads_df$StartPos[i] && e == main_tads_df$EndPos[i]) {
        is_main_tad <- TRUE
        break
      }
    }
    
    if (!is_main_tad) { # Only process if it's a true SubTAD (not a duplicate of a main TAD)
      subtad_start_pos <- c(subtad_start_pos, s)
      subtad_end_pos <- c(subtad_end_pos, e)
      subtad_length <- c(subtad_length, e - s)
      
      # Find parent TAD and assign hierarchical ID
      parent_tad_id_found <- NA
      for (i in 1:nrow(main_tads_df)) {
        tad_start <- main_tads_df$StartPos[i]
        tad_end <- main_tads_df$EndPos[i]
        tad_id <- main_tads_df$TAD_ID[i]
        
        if (s >= tad_start && e <= tad_end) {
          parent_tad_id_found <- tad_id
          break
        }
      }
      
      parent_tad_ids <- c(parent_tad_ids, parent_tad_id_found)
      
      if (!is.na(parent_tad_id_found)) {
        if (is.null(subtad_counter_per_tad[[parent_tad_id_found]])) {
          subtad_counter_per_tad[[parent_tad_id_found]] <- 1
        } else {
          subtad_counter_per_tad[[parent_tad_id_found]] <- subtad_counter_per_tad[[parent_tad_id_found]] + 1
        }
        sub_tad_sequential_ids <- c(sub_tad_sequential_ids, paste0("Sub", subtad_counter_per_tad[[parent_tad_id_found]]))
      } else {
        sub_tad_sequential_ids <- c(sub_tad_sequential_ids, "Unassigned")
      }
    }
  }
  
  # Write filtered raw lines to input_dir (java_raw_output_dir)
  writeLines(filtered_lines, output_window_txt)
  cat("✅ Filtered window lines saved to:\n -", normalizePath(output_window_txt, winslash = "/"), "\n")
  
  
  # Create SubTADs data frame
  sub_tad_df <- data.frame(
    Parent_TAD_ID = parent_tad_ids,
    SubTAD_ID = sub_tad_sequential_ids,
    StartPos = subtad_start_pos,
    EndPos = subtad_end_pos,
    TadLength = subtad_length
  )
  
  # Sort the SubTADs by Parent_TAD_ID and then by StartPos for consistent ordering
  sub_tad_df <- sub_tad_df[order(sub_tad_df$Parent_TAD_ID, sub_tad_df$StartPos), ]
  
  # Write SubTADs to output_dir (processed_tads_output_dir)
  write_tsv(sub_tad_df, output_subtad_txt)
  
  cat("✅ Processed SubTADs table saved to:\n -", normalizePath(output_subtad_txt, winslash = "/"), "\n")
  
  # Write Excel for SubTADs
  wb_subtad <- createWorkbook()
  addWorksheet(wb_subtad, "SubTADs")
  
  # Header
  writeData(wb_subtad, "SubTADs", data.frame(chrom = chrom, binsize = binsize), startRow = 1, colNames = TRUE)
  writeData(wb_subtad, "SubTADs", sub_tad_df, startRow = 3, colNames = TRUE)
  
  saveWorkbook(wb_subtad, output_subtad_excel, overwrite = TRUE)
  cat("✅ Excel file for SubTADs created:\n -", normalizePath(output_subtad_excel, winslash = "/"), "\n")
}
