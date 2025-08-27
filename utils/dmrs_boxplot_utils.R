# dmrs_boxplot_utils.R
# Author: Ghazal Sinjar
# Date: 23.06.2025
# Utility functions to support DMR-based CpG extraction, data reshaping, and boxplot visualization
# for methylation analysis using EPIC arrays in Shiny applications.

# ----------------------------------------------------------------------
# Libraries
# ----------------------------------------------------------------------
library(tidyr)
library(dplyr)
library(ggplot2)
library(plotly)


# ----------------------------------------------------------------------
# 1. Extract CpGs in a selected DMR region
# ----------------------------------------------------------------------

#' Extract CpGs within a specified DMR
#'
#' This function filters annotated CpG beta values within the genomic region
#' defined by a selected DMR and returns the relevant CpGs for visualization.
#'
#' @param DMRx Character. ID of the DMR to extract CpGs from.
#' @param DMR_table Data.frame. Table of identified DMRs, must contain columns: DMR_ID, chr, start, end.
#' @param annotated_table Data.frame. Annotated methylation beta matrix, must include `chr`, `pos`, and beta values.
#' @param pheno Data.frame. Phenotype information with rownames matching sample IDs and column `Sample_Group`.
#'
#' @return A data.frame containing CpGs in the DMR with beta values per sample.
#' Returns an error message string if inputs are missing or invalid.
#'
#' @export
extract_cpgs_in_DMR <- function(DMRx, DMR_table, annotated_table, pheno) {
  # Step 1: Check input
  if (missing(DMRx)) return("Error: 'DMRx' is missing.")
  if (is.null(annotated_table)) {
    return("Error: annotated_table is missing.")
  }
  if (is.null(pheno)) {
    return("Error: phenotype data (pheno) is missing.")
  }
  
  # Step 2: Check if DMRx exists
  if (!(DMRx %in% DMR_table$DMR_ID)) {
    return(paste("Error:", DMRx, "not found in DMR table"))
  }
  
  # Step 3: Extract DMR region info
  DMR_row <- DMR_table[DMR_table$DMR_ID == DMRx, ]
  chr <- DMR_row$chr
  start <- DMR_row$start
  end <- DMR_row$end
  
  # Step 4: Get beta columns using pheno
  m <- length(pheno)
  beta_cols <- tail(colnames(annotated_table), m)
  
  # Step 5: Subset annotation table for CpGs in the DMR
  region_cpgs <- annotated_table[
    annotated_table$chr == chr &
      annotated_table$pos >= start &
      annotated_table$pos <= end,
    c("chr", "pos", beta_cols)
  ]
  
  if (nrow(region_cpgs) == 0) {
    return("No CpGs found in the specified DMR region.")
  }
  
  return(region_cpgs)
}

'# test function
annotated_tbl_withqval <- readRDS("../modules/main_app_tests/box_pheno_pass/intermediate_data/annotated_with_qval.rds")
dim(annotated_tbl_withqval) #891492     40
sum(!is.na(annotated_tbl_withqval$qval)) # 60
sum(annotated_tbl_withqval$significance_last_qvalue == "Sig")# 60
results_dmr <- readRDS("./main_app_tests/box_pheno_pass/DMR_results/DMRs_unguided_vs_sample_groupguided_cutoff_-0.15_0.15_B0_20250804.rds")
head(results_dmr$dmr_table)
head(results_dmr$pheno_data)
cpgs_inDMR1 <- extract_cpgs_in_DMR(DMRx = "DMR1", DMR_table =results_dmr$dmr_table , annotated_table = annotated_tbl_withqval, pheno =results_dmr$pheno_data )
View(cpgs_inDMR1) # 10 rows

# case two the table is not merged with qval
annotated_tbl <- readRDS("../modules/main_app_tests/dmp_qval_test/intermediate_data/annotated_object_20250808.rds")
dim(annotated_tbl$annotated_table) #891492     38
cpgs_inDMR1_noq <- extract_cpgs_in_DMR(DMRx = "DMR1", DMR_table =results_dmr$dmr_table , annotated_table = annotated_tbl$annotated_table, pheno =results_dmr$pheno_data )
View(cpgs_inDMR1_noq)'
# ----------------------------------------------------------------------
# 2. Reshape CpG beta table to long format for plotting
# ----------------------------------------------------------------------

#' Reshape CpG beta matrix to long format
#'
#' Transforms a wide-format beta value matrix of CpGs into a long format
#' suitable for ggplot or plotly-based boxplots.
#'
#' @param region_cpgs Data.frame. Table returned by `extract_cpgs_in_DMR()`.
#' @param pheno_data Data.frame. Phenotype data with rownames as sample IDs and column `Sample_Group`.
#'
#' @return A long-format data.frame with CpG, SampleID, Group, BetaValue, etc.
#'
#' @export
#' Reshape CpG beta matrix to long format
#'
#' Transforms a wide-format beta value matrix of CpGs into a long format
#' suitable for ggplot or plotly-based boxplots.
#'
#' @param region_cpgs Data.frame. Table returned by `extract_cpgs_in_DMR()`.
#' @param pheno_data Data.frame. Phenotype data with rownames as sample IDs and column `Sample_Group`.
#'
#' @return A long-format data.frame with CpG, SampleID, Group, BetaValue, etc.
#'
#' @export
reshape_to_long_beta <- function(region_cpgs, pheno_data) {
  # Input validation checks
  if (missing(region_cpgs)) stop("Error: 'region_cpgs' table is missing")
  if (missing(pheno_data)) stop("Error: 'pheno_data' table is missing")
  
  # Convert to data frames
  region_cpgs_df <- as.data.frame(region_cpgs)
  pheno_df <- as.data.frame(pheno_data)
  
  # Check for the presence of the required 'chr' and 'pos' columns
  # This part is kept as per the original function's stop condition
  if (!all(c("chr", "pos") %in% colnames(region_cpgs_df))) {
    stop("Missing required columns in region_cpgs: 'chr' and/or 'pos'")
  }
  
  # Check for the required 'Sample_Group' column
  if (!"Sample_Group" %in% colnames(pheno_df)) {
    stop("Missing 'Sample_Group' column in pheno_data")
  }
  
  # --- START MODIFIED LOGIC ---
  # Determine if the optional 'qval' and 'significance_last_qvalue' columns exist
  has_qval_sig <- all(c("qval", "significance_last_qvalue") %in% colnames(region_cpgs_df))
  
  # Sample handling
  m <- nrow(pheno_df)
  
  # Dynamically determine beta columns based on the presence of qval/significance columns
  if (has_qval_sig) {
    # If qval and significance columns are present, assume they are the last two columns
    # Beta columns are the columns before them
    beta_cols <- colnames(region_cpgs_df)[(ncol(region_cpgs_df) - m - 1):(ncol(region_cpgs_df) - 2)]
  } else {
    # If qval and significance columns are missing, assume the last m columns are the beta values
    beta_cols <- tail(colnames(region_cpgs_df), m)
  }
  # --- END MODIFIED LOGIC ---
  
  # Additional validation checks
  if (length(beta_cols) != m) {
    stop("Number of beta columns doesn't match number of samples")
  }
  if (!all(beta_cols %in% rownames(pheno_df))) {
    stop("Some beta columns don't match pheno_data sample names")
  }
  
  # Validation of beta values
  beta_values <- region_cpgs_df[, beta_cols, drop = FALSE]
  if (!all(vapply(beta_values, is.numeric, logical(1)))) {
    stop("Non-numeric values found in beta columns")
  }
  if (any(beta_values < 0 | beta_values > 1, na.rm = TRUE)) {
    stop("Beta values outside [0,1] range detected")
  }
  
  # Create CpG identifier if not already present
  if (!"CpGName" %in% colnames(region_cpgs_df)) {
    region_cpgs_df$CpGName <- rownames(region_cpgs_df)
  }
  
  # Core processing pipeline
  result <- region_cpgs_df %>%
    dplyr::arrange(.data$pos) %>%
    tidyr::pivot_longer(
      cols = all_of(beta_cols),
      names_to = "SampleID",
      values_to = "BetaValue"
    ) %>%
    dplyr::mutate(
      SampleID = gsub("^X", "", SampleID),
      Group = pheno_df$Sample_Group[match(SampleID, rownames(pheno_df))],
      PositionLabel = paste(.data$CpGName, "\npos=", .data$pos)
    )
  
  # Conditionally select and mutate based on optional columns
  if (has_qval_sig) {
    result <- result %>%
      dplyr::select(
        "chr", "pos", "CpGName", "SampleID", "Group", "BetaValue", "PositionLabel", "qval", "significance_last_qvalue"
      ) %>%
      dplyr::mutate(
        qval = as.numeric(.data$qval),
        significance_last_qvalue = factor(.data$significance_last_qvalue, levels = c("Sig", "Insig"))
      )
  } else {
    result <- result %>%
      dplyr::select(
        "chr", "pos", "CpGName", "SampleID", "Group", "BetaValue", "PositionLabel"
      )
  }
  
  return(result)
}


'#test function
# input from privious function
long_dmr1 <-  reshape_to_long_beta(region_cpgs= cpgs_inDMR1, pheno_data =results_dmr$pheno_data) 
dim(long_dmr1) #80  7
View(long_dmr1)
# case 2: no qval
#long_DMR1_noq <- reshape_to_long_beta(region_cpgs= cpgs_inDMR1_noq, pheno_data =results_dmr$pheno_data)
#View(long_DMR1_noq)'
# ----------------------------------------------------------------------
# 3. Create boxplot of beta values across groups and CpGs
# ----------------------------------------------------------------------
# This function creates a boxplot for Beta values grouped by 'Group' and 'CpGName'.
# It can handle both static and interactive plots. The boxplot border colors
# are determined by a 'significance_last_qvalue' column.

create_boxplot <- function(
    long_format_table, 
    interactive = FALSE, 
    use_positional_spacing = FALSE,
    ref_group = NULL) {
  
  # Step 1: Input and data checks
  if (is.null(long_format_table) || nrow(long_format_table) == 0) {
    stop("Error: Input data table is empty or NULL.")
  }
  
  if (length(unique(long_format_table$chr)) != 1) {
    stop("Error: Not all CpGs are on one chromosome.")
  }
  
  # Check for significance data
  has_significance_data <- all(c("qval", "significance_last_qvalue") %in% names(long_format_table))
  
  if (has_significance_data) {
    long_format_table <- long_format_table %>%
      mutate(
        significance_last_qvalue = factor(significance_last_qvalue, levels = c("Sig", "Insig"))
      )
  }
  
  if (!use_positional_spacing) {
    long_format_table <- long_format_table %>%
      mutate(PositionLabel = factor(PositionLabel, levels = unique(long_format_table[order(long_format_table$pos), ]$PositionLabel)))
  }
  
  # Step 3: Dynamic color palette generation
  group_levels <- sort(unique(as.character(long_format_table$Group)))
  n_groups <- length(group_levels)
  
  if (n_groups == 2) {
    if (!is.null(ref_group)) {
      ref_group <- as.character(ref_group)
      if (!ref_group %in% group_levels) {
        warning("Specified ref_group '", ref_group, "' not found in data groups. Using default ordering.")
        ref_group <- NULL
      }
    }
    
    if (!is.null(ref_group)) {
      other_group <- setdiff(group_levels, ref_group)
      group_levels <- c(ref_group, other_group)
      color_palette <- c("blue", "orange")
    } else {
      color_palette <- c("blue", "orange")
    }
    names(color_palette) <- group_levels
  } else {
    color_palette <- rainbow(n_groups)
    names(color_palette) <- group_levels
  }
  
  # Step 4: Plot setup
  x_col <- if (use_positional_spacing) "pos" else "PositionLabel"
  
  # Define text aesthetics for tooltips
  long_format_table$text <- if (interactive) {
    if (use_positional_spacing) {
      paste0(long_format_table$SampleID, "\n", long_format_table$PositionLabel)
    } else {
      long_format_table$SampleID
    }
  } else {
    long_format_table$SampleID
  }
  
  # Calculate mean Beta values
  mean_values <- long_format_table %>%
    group_by(!!sym(x_col), Group) %>%
    summarise(mean_beta = mean(BetaValue, na.rm = TRUE), .groups = "drop")
  
  # Create grouping variable for proper dodging
  long_format_table$group_pos <- interaction(long_format_table[[x_col]], long_format_table$Group)
  mean_values$group_pos <- interaction(mean_values[[x_col]], mean_values$Group)
  
  # Set up position adjustments
  pos_dodge <- if (use_positional_spacing && !interactive) {
    position_dodge(width = 0.75)
  } else {
    position_dodge2(width = 0.75, preserve = "single")
  }
  
  pos_jit <- position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75)
  
  # Step 5: Create the base ggplot object
  p <- ggplot(long_format_table, aes_string(x = x_col, y = "BetaValue", fill = "Group")) +
    theme_minimal() +
    labs(
      x = if (use_positional_spacing) paste("Position on", unique(long_format_table$chr)) else paste("CPG + Position on", unique(long_format_table$chr)),
      y = "Beta Value",
      fill = "Group"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12))
  
  if (length(color_palette) > 0) {
    p <- p + scale_fill_manual(values = color_palette)
  }
  
  # Add boxplot with proper dodging
  if (use_positional_spacing && !interactive) {
    if (has_significance_data) {
      p <- p + geom_boxplot(
        aes(group = group_pos, color = significance_last_qvalue),
        outlier.shape = NA,
        lwd = 0.55,
        position = pos_dodge
      ) +
        scale_color_manual(
          name = "Significance",
          values = c("Sig" = "red", "Insig" = "black"),
          labels = c("Sig" = "Significant", "Insig" = "Non-significant")
        )
    } else {
      p <- p + geom_boxplot(
        aes(group = group_pos),
        color = "black",
        outlier.shape = NA,
        lwd = 0.5,
        position = pos_dodge
      )
    }
  } else {
    if (has_significance_data) {
      p <- p + geom_boxplot(
        aes(color = significance_last_qvalue),
        outlier.shape = NA,
        lwd = 0.55,
        position = pos_dodge
      ) +
        scale_color_manual(
          name = "Significance",
          values = c("Sig" = "red", "Insig" = "black"),
          labels = c("Sig" = "Significant", "Insig" = "Non-significant")
        )
    } else {
      p <- p + geom_boxplot(
        color = "black",
        outlier.shape = NA,
        lwd = 0.5,
        position = pos_dodge
      )
    }
  }
  
  # Add mean diamonds
  p <- p + geom_point(
    data = mean_values,
    aes(x = !!sym(x_col), y = mean_beta, group = group_pos),
    shape = 23, size = 3, color = "black", fill = "white",
    position = pos_dodge
  )
  
  # Add jittered points
  if (interactive) {
    p <- suppressWarnings(
      p + geom_jitter(
        aes(fill = Group, text = text), 
        position = pos_jit, 
        shape = 21,
        color = "black", 
        size = 2.5,
        stroke = 0.5
      )
    )
    return(plotly::ggplotly(p, tooltip = "text"))
  } else {
    p <- suppressWarnings(
      p + geom_jitter(
        aes(fill = Group),
        position = pos_jit, 
        shape = 21,
        color = "black",
        size = 2.5,
        stroke = 0.5
      )
    )
    return(p)
  }
}


'# test function

# make a cpg significant to color in green
#long_dmr1_manipulated <- long_dmr1 %>%
#  mutate(significance_last_qvalue = ifelse(row_number() <= 8, "Sig", "Insig"))

#View(long_dmr1_manipulated)
create_boxplot(
  long_format_table =long_dmr1_manipulated, 
  interactive =FALSE, 
  use_positional_spacing = FALSE,
  ref_group = "unguided")
# case 2 : no qval column
#create_boxplot(
#  long_format_table =long_DMR1_noq, 
#  interactive =FALSE, 
#  use_positional_spacing = TRUE,
#  ref_group = "unguided")
'