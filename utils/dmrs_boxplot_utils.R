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
reshape_to_long_beta <- function(region_cpgs, pheno_data) {
  # Input validation checks
  if (missing(region_cpgs)) stop("Error: 'region_cpgs' table is missing")
  if (missing(pheno_data)) stop("Error: 'pheno_data' table is missing")
  
  # Convert to data frames
  region_cpgs_df <- as.data.frame(region_cpgs)
  pheno_df <- as.data.frame(pheno_data)
  
  # Column validation
  required_cols <- c("chr", "pos")
  if (!all(required_cols %in% colnames(region_cpgs_df))) {
    stop("Missing required columns in region_cpgs: 'chr' and/or 'pos'")
  }
  if (!"Sample_Group" %in% colnames(pheno_df)) {
    stop("Missing 'Sample_Group' column in pheno_data")
  }
  
  # Sample handling
  m <- nrow(pheno_df)
  beta_cols <- tail(colnames(region_cpgs_df), m)
  beta_values <- region_cpgs_df[, beta_cols, drop = FALSE]
  
  if (!all(vapply(beta_values, is.numeric, logical(1)))) {
    stop("Non-numeric values found in beta columns")
  }
  if (any(beta_values < 0 | beta_values > 1, na.rm = TRUE)) {
    stop("Beta values outside [0,1] range detected")
  }
  
  # Create CpG identifier
  if (!"CpGName" %in% colnames(region_cpgs_df)) {
    region_cpgs_df$CpGName <- rownames(region_cpgs_df)
  }
  
  # Core processing pipeline
  region_cpgs_df %>%
    dplyr::arrange(.data$pos) %>%
    tidyr::pivot_longer(
      cols = all_of(beta_cols),
      names_to = "SampleID",
      values_to = "BetaValue"
    ) %>%
    dplyr::mutate(
      SampleID = gsub("^X", "", SampleID),  # X am Anfang entfernen!
      Group = pheno_df$Sample_Group[match(SampleID, rownames(pheno_df))],
      PositionLabel = paste(.data$CpGName, "\npos=", .data$pos)
    ) %>%
    dplyr::select(
      "chr", "pos", "CpGName", "SampleID", "Group", "BetaValue", "PositionLabel"
    )
}

# ----------------------------------------------------------------------
# 3. Create boxplot of beta values across groups and CpGs
# ----------------------------------------------------------------------

#' Create boxplot from reshaped CpG beta data
#'
#' Visualizes beta values of CpGs within a DMR across sample groups using
#' boxplots. Optionally supports interactive plots with plotly.
#'
#' @param long_format_table Data.frame. Output of `reshape_to_long_beta()`.
#' @param interactive Logical. If TRUE, generates a plotly plot.
#' @param use_positional_spacing Logical. If TRUE, x-axis is genomic position instead of CpG labels.
#'
#' @return A ggplot2 or plotly boxplot object.
#'
#' @export
'create_boxplot <- function(long_format_table, interactive = FALSE, use_positional_spacing = FALSE) {
  if (length(unique(long_format_table$chr)) != 1) {
    stop("Error: Not all CpGs are on one chromosome.")
  }
  
  # Correct x-axis variable name
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
  
  # Calculate mean values for each group at each position
  mean_values <- long_format_table %>%
    group_by(!!sym(x_col), Group) %>%
    summarise(mean_beta = mean(BetaValue, na.rm = TRUE), .groups = "drop")
  
  # Create a grouping variable that combines position and group
  long_format_table$group_pos <- interaction(long_format_table[[x_col]], long_format_table$Group)
  mean_values$group_pos <- interaction(mean_values[[x_col]], mean_values$Group)
  
  # Set up position adjustments
  pos_dodge <- if (use_positional_spacing && !interactive) {
    position_dodge(width = 0.75)
  } else {
    position_dodge2(width = 0.75, preserve = "single")
  }
  
  pos_jit <- position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75)
  
  # Create base plot
  p <- ggplot(long_format_table, aes_string(x = x_col, y = "BetaValue", fill = "Group")) +
    theme_minimal() +
    scale_fill_manual(values = c("guided" = "orange", "unguided" = "blue")) +
    scale_color_manual(values = c("guided" = "orange", "unguided" = "blue")) +
    labs(
      x = if (use_positional_spacing) paste("Position on", unique(long_format_table$chr)) else paste("CPG + Position on", unique(long_format_table$chr)),
      y = "Beta Value",
      fill = "Group"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Add boxplots with proper grouping
  if (use_positional_spacing && !interactive) {
    p <- p + geom_boxplot(
      aes(group = group_pos),
      outlier.shape = NA,
      color = "black",
      lwd = 0.5,
      position = pos_dodge
    )
  } else {
    p <- p + geom_boxplot(
      outlier.shape = NA,
      color = "black",
      lwd = 0.5,
      position = pos_dodge
    )
  }
  
  # Add mean diamonds
  p <- p + geom_point(
    data = mean_values,
    aes(x = !!sym(x_col), y = mean_beta, group = group_pos),
    shape = 23, size = 3, color = "black", fill = "white",
    position = pos_dodge
  )
  
  if (interactive) {
    p <- suppressWarnings(
      p + geom_jitter(aes_string(color = "Group", text = "text"), 
                      position = pos_jit, 
                      shape = 16, size = 1.5)
    )
    return(plotly::ggplotly(p, tooltip = "text"))
  } else {
    p <- suppressWarnings(
      p + geom_jitter(aes_string(color = "Group"), 
                      position = pos_jit, 
                      shape = 16, size = 1.5)
    )
    return(p)
  }
}'
create_boxplot <- function(long_format_table, interactive = FALSE, use_positional_spacing = FALSE) {
  if (length(unique(long_format_table$chr)) != 1) {
    stop("Error: Not all CpGs are on one chromosome.")
  }
  
  # Correct x-axis variable name
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
  
  # Calculate mean values for each group at each position
  mean_values <- long_format_table %>%
    group_by(!!sym(x_col), Group) %>%
    summarise(mean_beta = mean(BetaValue, na.rm = TRUE), .groups = "drop")
  
  # Create a grouping variable that combines position and group
  long_format_table$group_pos <- interaction(long_format_table[[x_col]], long_format_table$Group)
  mean_values$group_pos <- interaction(mean_values[[x_col]], mean_values$Group)
  
  # Set up position adjustments
  pos_dodge <- if (use_positional_spacing && !interactive) {
    position_dodge(width = 0.75)
  } else {
    position_dodge2(width = 0.75, preserve = "single")
  }
  
  pos_jit <- position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75)
  
  # Create base plot
  p <- ggplot(long_format_table, aes_string(x = x_col, y = "BetaValue", fill = "Group")) +
    theme_minimal() +
    scale_fill_manual(values = c("guided" = "orange", "unguided" = "blue")) +
    scale_color_manual(values = c("guided" = "orange", "unguided" = "blue")) +
    labs(
      x = if (use_positional_spacing) paste("Position on", unique(long_format_table$chr)) else paste("CPG + Position on", unique(long_format_table$chr)),
      y = "Beta Value",
      fill = "Group"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Add boxplots with proper grouping
  if (use_positional_spacing && !interactive) {
    p <- p + geom_boxplot(
      aes(group = group_pos),
      outlier.shape = NA,
      color = "black",
      lwd = 0.5,
      position = pos_dodge
    )
  } else {
    p <- p + geom_boxplot(
      outlier.shape = NA,
      color = "black",
      lwd = 0.5,
      position = pos_dodge
    )
  }
  
  # Add mean diamonds
  p <- p + geom_point(
    data = mean_values,
    aes(x = !!sym(x_col), y = mean_beta, group = group_pos),
    shape = 23, size = 3, color = "black", fill = "white",
    position = pos_dodge
  )
  
  if (interactive) {
    p <- suppressWarnings(
      p + geom_jitter(
        aes_string(color = "Group", text = "text"), 
        position = pos_jit, 
        shape = 21,  # Use shape 21-25 for fillable points with borders
        size = 2.5,  # Slightly larger to see the border
        stroke = 0.5  # Border thickness
      )
    )
    return(plotly::ggplotly(p, tooltip = "text"))
  } else {
    p <- suppressWarnings(
      p + geom_jitter(
        aes_string(fill = "Group"),  # Use fill instead of color for shape 21
        position = pos_jit, 
        shape = 21,  # Shape that supports both fill and color
        color = "black",  # Edge color
        size = 2.5,      # Slightly larger size
        stroke = 0.5      # Border thickness
      ) +
        # Need to add fill scale since we're using fill aesthetic
        scale_fill_manual(values = c("guided" = "orange", "unguided" = "blue"))
    )
    return(p)
  }
}

