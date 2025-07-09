# dmrs_utils.R
# Author: Ghazal Sinjar
# Date: 16.06.2025
# Utility functions for DMR detection, processing, and annotation in the EPIC methylation array pipeline.

#libraries
library(doParallel)


# ----------------------------------------------------------------------
# 1. Run bumphunter to detect DMRs
# ----------------------------------------------------------------------

#' Run Bumphunter to detect Differentially Methylated Regions (DMRs)
#'
#' This function applies the bumphunter algorithm to identify differentially
#' methylated regions (DMRs) using a GenomicRatioSet or RGChannelSet object.
#' It includes internal error handling for both phenotype extraction and
#' the bumphunter process.
#'
#' @param rgSet A GenomicRatioSet or RGChannelSet object (from the minfi package).
#' @param cutoff Numeric. Methylation cutoff threshold for bump detection.
#' @param B Integer. Number of permutations to use for significance estimation (0 = none).
#' @param num_cores 
#'
#' @return A list containing:
#' \describe{
#'   \item{DMR_table}{A data.frame with the identified DMRs. Returns an empty data.frame if no DMRs are found or NULL on error.}
#'   \item{pd}{Phenotype data.frame extracted from rgSet.}
#' }
#' Returns NULL if bumphunter fails to run.
#'
#' @importFrom minfi pData
#' @importFrom bumphunter bumphunter
#' @importFrom doParallel registerDoParallel
#'
#' @export
run_bumphunter_dmrs <- function(rgSet, cutoff = 0.15, B = 0, num_cores = 1) {
  # Extract phenotype data
  pd <- tryCatch({
    pData(rgSet)
  }, error = function(e) {
    stop(paste("Error extracting phenotype data:", e$message))
  })
  
  if (!"Sample_Group" %in% colnames(pd)) {
    stop("Column 'Sample_Group' not found in phenotype data.")
  }
  
  # Prepare phenotype factor
  sample_group <- pd$Sample_Group
  if (!is.factor(sample_group)) sample_group <- factor(sample_group)
  sample_group <- relevel(sample_group, ref = "unguided")
  
  # Create design matrix
  designMatrix <- model.matrix(~ sample_group)
  
  # Set up parallel backend
  if (!requireNamespace("doParallel", quietly = TRUE)) {
    stop("Package 'doParallel' is required for parallel execution.")
  }
  doParallel::registerDoParallel(cores = num_cores) 
  
  # Run bumphunter with error handling
  dmrs <- tryCatch({
    bumphunter(
      rgSet,
      design = designMatrix,
      coef = 2,
      cutoff = cutoff,
      B = B,
      type = "Beta",
      smooth = TRUE
    )
  }, error = function(e) {
    warning(paste("Error during bumphunter execution:", e$message))
    return(NULL)
  })
  
  # If bumphunter failed
  if (is.null(dmrs)) {
    return(NULL)
  }
  
  # Convert to data frame
  dmrs_tbl <- as.data.frame(dmrs$table)
  rownames(dmrs_tbl) <- NULL
  
  return(list(
    DMR_table = dmrs_tbl,
    pd = pd
  ))
}

# test function
#rgset <- readRDS("C:/Users/ghaza/Documents/ghazal/Bioinformatik_Fächer/Masterarbeit_Project/Scripts/R_Scripts/intermediate_data/filtered_GRset_SWAN_SNP_removed_SexChr_kept_20250605.rds")
#dmr <- run_bumphunter_dmrs(rgSet =  rgset, num_cores = 6 )

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
  desired_cols <- c("width", "mean.methy.difference", "CpGs.inDMR", "p.value", "cluster", "CpGs.inCluster")
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
#' @return GRanges object with additional columns: overlapped_gene_name, concatenated_gene_names, and unique DMR_ID.
annotate_dmrs_with_genes <- function(GR_Dmrs, tx_gr_filtered) {
  
  # -------------------- STEP 1: Get first overlapping gene for each DMR --------------------
  # Finds the first gene (by order) that overlaps with each DMR, ignoring strand direction
  hits_first <- findOverlaps(GR_Dmrs, tx_gr_filtered, ignore.strand = TRUE, select = "first")
  
  # Add the name of the first overlapping gene to the DMR metadata
  mcols(GR_Dmrs)$overlapped_gene_name <- tx_gr_filtered$gene_name[hits_first]
  
  # -------------------- STEP 2: Get all overlapping gene names concatenated --------------------
  # Find all overlaps between DMRs and genes
  hits_all <- findOverlaps(GR_Dmrs, tx_gr_filtered, ignore.strand = TRUE)
  
  # Create a named list: for each DMR index, list of gene names that overlap
  overlapping_gene_names_list <- split(tx_gr_filtered$gene_name[subjectHits(hits_all)], queryHits(hits_all))
  
  # Concatenate gene names per DMR into a single string (e.g. "Gene1;Gene2")
  concatenated_gene_names_vector <- sapply(overlapping_gene_names_list, paste, collapse = ";")
  
  # Create a full-length character vector, fill with NA initially
  all_concat_names <- rep(NA_character_, length(GR_Dmrs))
  
  # Insert concatenated names at the correct positions
  all_concat_names[as.integer(names(concatenated_gene_names_vector))] <- concatenated_gene_names_vector
  
  # Add concatenated gene names to metadata
  mcols(GR_Dmrs)$concatenated_gene_names <- all_concat_names
  
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
