# DMP_utils.R
# Author: Ghazal Sinjar

# Funtion 1 -------------------------------------------------------------------------
#' Calculate mean differences for Beta and M values between two groups.
#' # PS. M_mean_difference is the fold change that can be used to volcano-plot
calculateMethylationMeanDifferences <- function(GenoRationSet, pheno_table, ref_group, round_digits = 5) {
  
  # 1. Get Beta and M values using tryCatch for robust error handling
  beta_values <- tryCatch({
    # Attempt to get Beta values. This assumes 'getBeta' is available
    # and compatible with 'GenoRationSet' (e.g., from 'minfi').
    getBeta(GenoRationSet)
  }, error = function(e) {
    message("Error getting Beta values: ", e$message)
    return(NULL) # Return NULL on error
  })
  
  M_values <- tryCatch({
    # Attempt to get M values. This assumes 'getM' is available
    # and compatible with 'GenoRationSet' (e.g., from 'minfi').
    getM(GenoRationSet)
  }, error = function(e) {
    message("Error getting M values: ", e$message)
    return(NULL) # Return NULL on error
  })
  
  # Check if Beta or M values could not be retrieved
  if (is.null(beta_values) || is.null(M_values)) {
    stop("Could not retrieve Beta or M values. Please check 'GenoRationSet' object.")
  }
  
  # Ensure column names of beta/M matrices match row names of pheno_table.
  # This is crucial for correctly mapping samples to their groups.
  if (!all(colnames(beta_values) %in% rownames(pheno_table))) {
    stop("Column names of beta values do not fully match row names of pheno table.
          Please ensure sample IDs are consistent.")
  }
  if (!all(colnames(M_values) %in% rownames(pheno_table))) {
    stop("Column names of M values do not fully match row names of pheno table.
          Please ensure sample IDs are consistent.")
  }
  
  # Order the pheno_table to match the column order of the methylation data.
  # This ensures correct group assignment for each sample column.
  pheno_ordered <- pheno_table[colnames(beta_values), , drop = FALSE]
  
  # 2. Identify and validate groups
  if (!"Sample_Group" %in% colnames(pheno_ordered)) {
    stop("Phenotype table must contain a column named 'Sample_Group'.")
  }
  
  all_groups <- unique(pheno_ordered$Sample_Group)
  if (!ref_group %in% all_groups) {
    stop(paste0("Reference group '", ref_group, "' not found in the 'Sample_Group' column.
                Available groups are: ", paste(all_groups, collapse = ", ")))
  }
  
  # Identify the 'other' group. This function is designed for exactly two groups.
  other_groups <- setdiff(all_groups, ref_group)
  if (length(other_groups) != 1) {
    stop(paste0("This function is designed for exactly two groups (reference and one other).
                Found groups: ", paste(all_groups, collapse = ", ")))
  }
  other_group_name <- other_groups[1]
  
  # Get the sample IDs for each group based on the ordered phenotype table.
  ref_sample_names <- rownames(pheno_ordered)[pheno_ordered$Sample_Group == ref_group]
  other_sample_names <- rownames(pheno_ordered)[pheno_ordered$Sample_Group == other_group_name]
  
  # 3. Subset Beta values by group
  # 'drop = FALSE' ensures that the result remains a matrix even if only one sample is present.
  beta_ref <- beta_values[, ref_sample_names, drop = FALSE]
  beta_other <- beta_values[, other_sample_names, drop = FALSE]
  
  # 4. Compute row-wise mean for each group (Beta values)
  mean_beta_ref <- rowMeans(beta_ref, na.rm = TRUE) # na.rm=TRUE handles potential NA values
  mean_beta_other <- rowMeans(beta_other, na.rm = TRUE)
  
  # 5. Compute Beta mean difference: log2(other_group) - log2(ref_group).
  # Add a small epsilon to avoid log2(0) or log2 of very small numbers, which can result in -Inf.
  epsilon <- 1e-9 # one of a billion
  Beta_mean_difference <- log2(mean_beta_other + epsilon) - log2(mean_beta_ref + epsilon)
  
  # 6. Subset M values by group
  M_ref <- M_values[, ref_sample_names, drop = FALSE]
  M_other <- M_values[, other_sample_names, drop = FALSE]
  
  # 7. Compute row-wise mean for each group (M values)
  mean_M_ref <- rowMeans(M_ref, na.rm = TRUE)
  mean_M_other <- rowMeans(M_other, na.rm = TRUE)
  
  # 8. Compute M mean difference: other_group - ref_group
  M_mean_difference_FC <- mean_M_other - mean_M_ref
  
  # 9. Round the final mean differences before combining
  M_mean_difference_rounded <- round(M_mean_difference_FC, round_digits)
  Beta_mean_difference_rounded <- round(Beta_mean_difference, round_digits)
  
  # 10. Combine rounded mean differences into a single data frame
  # Row names will be the CpG IDs from the original methylation data.
  combined_results <- data.frame(
    M_mean_difference_FC = M_mean_difference_rounded,
    Beta_mean_difference = Beta_mean_difference_rounded,
    row.names = rownames(beta_values) # CpG IDs
  )
  
  return(combined_results)
}

'# test fucntion 1:
dir <- "C:/Users/ghaza/Documents/ghazal/Bioinformatik_Fächer/Masterarbeit_Project/Scripts/R_Scripts/modules/main_app_tests/box_pheno_pass/intermediate_data/filtered_GRset_SWAN_SNPsremoved_SexChrProbes_kept_20250804.rds"
GenoRationSet <- readRDS(dir)
pheno_table <- pData(GenoRationSet)
ref_group <- "unguided"
mean_beta_M_tbl <- calculateMethylationMeanDifferences(GenoRationSet = GenoRationSet, pheno_table = pheno_table, ref_group = ref_group)
head(mean_beta_M_tbl)
dim(mean_beta_M_tbl)'

# Function 2 -------------------------------------------------------------------
#Add M and Beta mean difference columns to a DMP (Differentially Methylated Positions) table.
addMeanDifferencesToDMP <- function(dmp_table, combined_Beta_M_mean_table) {
  tryCatch({
    # Validate inputs
    if (is.null(rownames(dmp_table))) {
      stop("'dmp_table' must have row names (CpG IDs)")
    }
    if (is.null(rownames(combined_Beta_M_mean_table))) {
      stop("'combined_Beta_M_mean_table' must have row names")
    }
    
    # Convert to data frames
    dmp_df <- as.data.frame(dmp_table)
    mean_df <- as.data.frame(combined_Beta_M_mean_table)
    
    # Find matching CpGs
    common_cpgs <- base::intersect(rownames(dmp_df), rownames(mean_df))
    if (length(common_cpgs) == 0) {
      stop("No matching CpG IDs between tables")
    }
    
    message(paste("Successfully matched", length(common_cpgs), 
                  "CpGs between DMP and mean difference tables"))
    
    # Subset to common CpGs
    dmp_df <- dmp_df[common_cpgs, , drop = FALSE]
    mean_df <- mean_df[common_cpgs, , drop = FALSE]
    
    # Verify alignment before merging
    if (!identical(rownames(dmp_df), rownames(mean_df))) {
      stop("Row name alignment failed after subsetting")
    }
    
    # Add mean difference columns
    merged_df <- cbind(
      dmp_df,
      M_mean_difference_FC = mean_df$M_mean_difference_FC,
      Beta_mean_difference = mean_df$Beta_mean_difference
    )
    
    # Round numeric columns
    numeric_cols <- sapply(merged_df, is.numeric)
    if (any(numeric_cols)) {
      merged_df[numeric_cols] <- lapply(merged_df[numeric_cols], function(x) {
        if (is.numeric(x)) signif(x, digits = 5) else x
      })
    }
    # Remove 'intercept' and 'f' columns
    merged_df <- merged_df[, !names(merged_df) %in% c("intercept", "f"), drop = FALSE]
    
    return(merged_df)
    
  }, error = function(e) {
    # Enhanced error reporting
    message("\nERROR DETAILS:")
    message("Input dimensions:")
    message("- DMP table: ", paste(dim(dmp_table), collapse = " x "))
    message("- Mean table: ", paste(dim(combined_Beta_M_mean_table), collapse = " x "))
    if (exists("common_cpgs")) {
      message("Matching CpGs: ", length(common_cpgs))
    }
    message("Error message: ", e$message)
    return(NULL)
  })
}

'#test function
pheno_table <- pData(GenoRationSet$filtered_data)
pheno <- pheno_table$Sample_Group
dmp <- dmpFinder(GenoRationSet$filtered_data, pheno = pheno, type = "categorical", qCutoff = 0.6, shrinkVar = FALSE) 
dim(dmp)
head(dmp)
edited_dmp <- addMeanDifferencesToDMP(dmp_table = dmp, combined_Beta_M_mean_table = mean_beta_M_tbl)
head(edited_dmp)
dim(edited_dmp)
# check if it worked
# Correct way 1:
edited_dmp["cg12521566_BC21", ]$Beta_mean_difference  
mean_beta_M_tbl["cg12521566_BC21", "Beta_mean_difference"]     
# Correct way 2:
mean_beta_M_tbl["cg15620155_TC21",]$Beta_mean_difference      
mean_beta_M_tbl$Beta_mean_difference[rownames(mean_beta_M_tbl) == "cg15620155_TC21"] '
#--------------------------------------------------------