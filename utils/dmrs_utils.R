# dmrs_utils.R
# Author: Ghazal Sinjar
# Date: 16.06.2025
# Utility functions for DMR detection, processing, and annotation in the EPIC methylation array pipeline.

#libraries
library(doParallel)
library(conflicted) # Make sure you have this package installed
conflicts_prefer(base::setdiff)

# ----------------------------------------------------------------------
# 1. Run bumphunter to detect DMRs
# ----------------------------------------------------------------------
#' Run Bumphunter to detect Differentially Methylated Regions (DMRs)
#'
#' This function applies the bumphunter algorithm to identify differentially
#' methylated regions (DMRs) by comparing a specified group against a
#' defined reference group.
#'
#' @param rgSet A GenomicRatioSet or RGChannelSet object (from the minfi package).
#' @param ref_sample_group Character. The name of the group to be used as the
#'   reference level for the bumphunter comparison. Must be a valid level in
#'   the 'Sample_Group' column.
#' @param tested_group_colname Character. The exact name of the column in the
#'   design matrix that corresponds to the group to be tested against the
#'   reference.
#' @param cutoff Numeric. Methylation cutoff threshold for bump detection.
#' @param B Integer. Number of permutations to use for significance estimation (0 = none).
#' @param num_cores Integer. Number of cores to use for parallel processing.
#'
#' @return A data.frame with the identified DMRs. Returns an empty data.frame
#'   if no DMRs are found or NULL on error.
#'
#' @importFrom minfi pData
#' @importFrom bumphunter bumphunter
#' @importFrom doParallel registerDoParallel
#'
#' @export
run_bumphunter_dmrs <- function(rgSet, ref_sample_group, tested_group_colname,
                                cutoff = 0.15, B = 0, num_cores = 1) {
  
  # --- Helper: cross-platform RAM detection ---
  get_total_ram_gb <- function() {
    sysname <- Sys.info()[["sysname"]]
    
    if (sysname == "Windows") {
      mem_mb <- suppressWarnings(utils::memory.limit()) # returns MB
      if (!is.na(mem_mb) && mem_mb > 0) return(mem_mb / 1024)
    } else if (sysname == "Linux") {
      mem_kb <- suppressWarnings(as.numeric(
        system("grep MemTotal /proc/meminfo | awk '{print $2}'", intern = TRUE)
      ))
      if (!is.na(mem_kb)) return(mem_kb / 1024 / 1024)
    } else if (sysname == "Darwin") {
      mem_bytes <- suppressWarnings(as.numeric(
        system("sysctl -n hw.memsize", intern = TRUE)
      ))
      if (!is.na(mem_bytes)) return(mem_bytes / 1024^3)
    }
    
    return(NA_real_)  # unknown RAM
  }
  
  # --- Memory check ---
  dataset_size_mb <- as.numeric(object.size(rgSet)) / (1024^2)
  total_ram_gb <- get_total_ram_gb()
  
  if (!is.na(total_ram_gb)) {
    est_usage_gb <- (dataset_size_mb * num_cores) / 1024
    if (est_usage_gb > total_ram_gb * 0.8) {
      message(sprintf(
        "⚠ Estimated memory usage %.2f GB exceeds 80%% of total RAM (%.2f GB). Reducing cores.",
        est_usage_gb, total_ram_gb
      ))
      num_cores <- max(1, floor((total_ram_gb * 0.8) / (dataset_size_mb / 1024)))
    }
  }
  
  message(sprintf("Using %d cores for bumphunter analysis.", num_cores))
  
  # --- Extract phenotype data ---
  pd <- tryCatch({
    pData(rgSet)
  }, error = function(e) {
    stop(paste("Error extracting phenotype data:", e$message))
  })
  
  if (!"Sample_Group" %in% colnames(pd)) {
    stop("Column 'Sample_Group' not found in phenotype data.")
  }
  
  sample_group <- pd$Sample_Group
  if (!is.factor(sample_group)) sample_group <- factor(sample_group)
  
  if (!ref_sample_group %in% levels(sample_group)) {
    stop(sprintf("Reference group '%s' not found in 'Sample_Group' levels.", ref_sample_group))
  }
  
  sample_group <- relevel(sample_group, ref = ref_sample_group)
  designMatrix <- model.matrix(~ sample_group)
  
  group_col_index <- which(colnames(designMatrix) == tested_group_colname)
  if (length(group_col_index) == 0) {
    stop(sprintf(
      "The tested group column '%s' was not found in the design matrix. Available columns: %s",
      tested_group_colname, paste(colnames(designMatrix), collapse = ", ")
    ))
  }
  
  # --- Parallel setup with PSOCK ---
  if (!requireNamespace("doParallel", quietly = TRUE)) {
    stop("Package 'doParallel' is required for parallel execution.")
  }
  options(snowFTTimeout = 300)  # 5 min timeout
  cl <- parallel::makeCluster(num_cores, type = "PSOCK")
  doParallel::registerDoParallel(cl)
  
  # --- Run bumphunter with special error handling ---
  dmrs <- tryCatch({
    bumphunter(
      rgSet,
      design = designMatrix,
      coef = group_col_index,
      cutoff = cutoff,
      B = B,
      type = "Beta",
      smooth = TRUE
    )
  }, error = function(e) {
    if (grepl("Schreibfehler in Verbindung", e$message)) {
      warning("❌ Error during DMR detection: a connection write error during the parallelized bumphunter. Likely causes: Out of memory (RAM).")
    } else {
      warning(paste("Error during bumphunter execution:", e$message))
    }
    return(NULL)
  })
  
  # Always stop cluster
  parallel::stopCluster(cl)
  
  if (is.null(dmrs)) return(NULL)
  
  dmrs_tbl <- as.data.frame(dmrs$table)
  rownames(dmrs_tbl) <- NULL
  return(dmrs_tbl)
}




# ----------------------------------------------------------------------
# 2. Convert DMRs to GRanges and process metadata
# ----------------------------------------------------------------------

library(GenomicRanges)
library(dplyr)
library(tibble)

#' Prepare DMRs as a GRanges object with cleaned metadata
#'
#' @param dmrs_df Data.frame of DMRs as returned by `run_bumphunter_dmrs`.
#'
#' @return A sorted GRanges object with renamed and rounded metadata.
prepare_dmrs_granges <- function(dmrs_df) {
  # Drop unwanted columns if they exist
  cols_to_remove <- c("indexStart", "indexEnd", "area", "cluster")
  dmrs_df <- dmrs_df[, !(names(dmrs_df) %in% cols_to_remove)]
  
  gr <- makeGRangesFromDataFrame(
    df = dmrs_df,
    seqnames.field = "chr",
    start.field = "start",
    end.field = "end",
    keep.extra.columns = TRUE
  )
  
  # Rename common fields if present
  if ("value" %in% names(mcols(gr))) names(mcols(gr))[names(mcols(gr)) == "value"] <- "mean.methy.difference"
  if ("L" %in% names(mcols(gr))) names(mcols(gr))[names(mcols(gr)) == "L"] <- "CpGs.inDMR"
  if ("clusterL" %in% names(mcols(gr))) names(mcols(gr))[names(mcols(gr)) == "clusterL"] <- "CpGs.inCluster"
  
  # Add width and sort
  mcols(gr)$width <- width(gr)
  if ("CpGs.inDMR" %in% names(mcols(gr))) {
    gr <- gr[order(gr$CpGs.inDMR, decreasing = TRUE)]
  }
  
  # Reorder and round metadata
  desired_cols <- c("width", "mean.methy.difference", "CpGs.inDMR", "p.value", "CpGs.inCluster")
  all_cols <- names(mcols(gr))
  cols_to_select <- desired_cols[desired_cols %in% all_cols]
  
  reordered_mcols <- mcols(gr) %>%
    as_tibble() %>%
    dplyr::select(any_of(cols_to_select), everything())
  
  is_roundable <- sapply(reordered_mcols, function(x) {
    is.numeric(x) && !all(x == floor(x), na.rm = TRUE)
  })
  reordered_mcols[is_roundable] <- lapply(reordered_mcols[is_roundable], function(x) round(x, 5))
  
  mcols(gr) <- as(reordered_mcols, "DataFrame")
  return(gr)
}


# ----------------------------------------------------------------------
# 3. Annotate DMRs with overlapping gene names
# ----------------------------------------------------------------------

#' Annotate DMR GRanges with overlapping gene names
#'
#' @param GR_Dmrs GRanges object of DMRs.
#' @param tx_gr_filtered GRanges object of transcript or gene annotation.
#'
#' @return GRanges object with additional columns: first_overlapped_gene, All_overlapped_genes, and unique DMR_ID.
annotate_dmrs_with_genes <- function(GR_Dmrs, tx_gr_filtered) {
  
  # -------------------- STEP 1: Get first overlapping gene for each DMR --------------------
  # Finds the first gene (by order) that overlaps with each DMR, ignoring strand direction
  hits_first <- findOverlaps(GR_Dmrs, tx_gr_filtered, ignore.strand = TRUE, select = "first")
  
  # Add the name of the first overlapping gene to the DMR metadata
  mcols(GR_Dmrs)$first_overlapped_gene <- tx_gr_filtered$gene_name[hits_first]
  
  # -------------------- STEP 2: Get all overlapping gene names concatenated --------------------
  # Find all overlaps between DMRs and genes
  hits_all <- findOverlaps(GR_Dmrs, tx_gr_filtered, ignore.strand = TRUE)
  
  # Create a named list: for each DMR index, list of gene names that overlap
  overlapping_gene_names_list <- split(tx_gr_filtered$gene_name[subjectHits(hits_all)], queryHits(hits_all))
  
  # Concatenate gene names per DMR into a single string (e.g. "Gene1;Gene2")
  All_overlapped_genes_vector <- sapply(overlapping_gene_names_list, paste, collapse = ";")
  
  # Create a full-length character vector, fill with NA initially
  all_concat_names <- rep(NA_character_, length(GR_Dmrs))
  
  # Insert concatenated names at the correct positions
  all_concat_names[as.integer(names(All_overlapped_genes_vector))] <- All_overlapped_genes_vector
  
  # Add concatenated gene names to metadata
  mcols(GR_Dmrs)$All_overlapped_genes <- all_concat_names
  
  # -------------------- STEP 3: Create unique DMR IDs based on gene names --------------------
  dmr_numbers <- seq_len(length(GR_Dmrs))
  
  # Paste "DMR" in front of each number to create the desired IDs
  mcols(GR_Dmrs)$DMR_ID <- paste0("DMR", dmr_numbers)
  
  # Return annotated GRanges
  return(GR_Dmrs)
  
  # -------------------- STEP 5: Assign final DMR_IDs back to GRanges --------------------
  # Create vector of same length as GR_Dmrs to hold DMR_IDs
  DMR_IDs <- rep(NA_character_, length(GR_Dmrs))
  
  # Fill in the DMR IDs at their respective positions
  DMR_IDs[dmr_gene_df$DMR_idx] <- dmr_gene_df$DMR_ID
  
  # Add to GRanges metadata
  mcols(GR_Dmrs)$DMR_ID <- DMR_IDs
  
  # Return annotated GRanges
  return(GR_Dmrs)
}


# ----------------------------------------------------------------------
# 4. a helper function to generate the filename
# ----------------------------------------------------------------------

get_dmr_base_filename <- function(cutoff_from, cutoff_to, B_val) {
  cutoff_from_str <- sub("-", "neg", sprintf("%.2f", cutoff_from)) 
  cutoff_to_str <- sprintf("%.2f", cutoff_to)
  b_str <- as.character(B_val)
  paste0("DMRs_cutoff_", cutoff_from_str, "_to_", cutoff_to_str, "_B", b_str, "_", Sys.Date())
}
