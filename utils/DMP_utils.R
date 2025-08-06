# DMP_utils.R


# Funtion 1 -------------------------------------------------------------------------
#' Calculate mean differences for Beta and M values between two groups.
#' # PS. M_mean_difference is the fold change that can be used to volcano-plot
calculateMethylationMeanDifferences <- function(normalised_methylset, pheno_table, ref_group, round_digits = 5) {
  
  # 1. Get Beta and M values using tryCatch for robust error handling
  beta_values <- tryCatch({
    # Attempt to get Beta values. This assumes 'getBeta' is available
    # and compatible with 'normalised_methylset' (e.g., from 'minfi').
    getBeta(normalised_methylset)
  }, error = function(e) {
    message("Error getting Beta values: ", e$message)
    return(NULL) # Return NULL on error
  })
  
  M_values <- tryCatch({
    # Attempt to get M values. This assumes 'getM' is available
    # and compatible with 'normalised_methylset' (e.g., from 'minfi').
    getM(normalised_methylset)
  }, error = function(e) {
    message("Error getting M values: ", e$message)
    return(NULL) # Return NULL on error
  })
  
  # Check if Beta or M values could not be retrieved
  if (is.null(beta_values) || is.null(M_values)) {
    stop("Could not retrieve Beta or M values. Please check 'normalised_methylset' object and 'minfi' package availability.")
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
#dir <- "C:/Users/ghaza/Documents/ghazal/Bioinformatik_Fächer/Masterarbeit_Project/Scripts/R_Scripts/modules/main_app_tests/box_pheno_pass/intermediate_data/Five_objects_normalised_data.rds"
#methylsets <- readRDS(dir)
#pheno_table <- pData(methylsets$SWAN)
#ref_group <- "unguided"
mean_beta_M_tbl <- calculateMethylationMeanDifferences(normalised_methylset = methylsets$SWAN, pheno_table = pheno_table, ref_group = ref_goup)
head(mean_beta_M_tbl)'

# Function 2 -------------------------------------------------------------------
#Add M and Beta mean difference columns to a DMP (Differentially Methylated Positions) table.
addMeanDifferencesToDMP <- function(dmp_table, combined_Beta_M_mean_table) {
  tryCatch({
    # Ensure both tables have row names
    if (is.null(rownames(dmp_table))) {
      stop("'dmp_table' must have row names (CpG IDs).")
    }
    if (is.null(rownames(combined_Beta_M_mean_table))) {
      stop("'combined_Beta_M_mean_table' must have row names (CpG IDs).")
    }
    
    # Convert to data.frame if they aren't already
    dmp_df <- as.data.frame(dmp_table)
    mean_df <- as.data.frame(combined_Beta_M_mean_table)
    
    # Add the mean differences directly by row name matching
    dmp_df$M_mean_difference_FC <- mean_df[rownames(dmp_df), "M_mean_difference_FC"]
    dmp_df$Beta_mean_difference <- mean_df[rownames(dmp_df), "Beta_mean_difference"]
    # Round all numeric columns to 5 digits
    numeric_cols <- sapply(dmp_df, is.numeric)
    dmp_df[numeric_cols] <- lapply(dmp_df[numeric_cols], function(x) {
      if (is.numeric(x)) {
        # First round to 5 digits
        x <- round(x, digits = 5)
        # Then format to prevent scientific notation
        format(x, scientific = FALSE, trim = TRUE)
      } else {
        x
      }
    })
    
    # Convert formatted numbers back to numeric (without scientific notation)
    dmp_df[numeric_cols] <- lapply(dmp_df[numeric_cols], function(x) {
      if (is.character(x)) {
        as.numeric(x)
      } else {
        x
      }
    })
    return(dmp_df)
    
  }, error = function(e) {
    message("Error adding mean differences to DMP table: ", e$message)
    return(NULL)
  })
}

'#test function
pheno_table <- pData(methylsets$SWAN)
pheno <- pheno_table$Sample_Group
dmp <- dmpFinder(methylsets$SWAN, pheno = pheno, type = "categorical", qCutoff = 1, shrinkVar = FALSE) 
edited_dmp <- addMeanDifferencesToDMP(dmp_table = dmp, combined_Beta_M_mean_table = mean_beta_M_tbl)
# check if it worked
# Correct way 1:
edited_dmp["cg24426691_TC11", ]$Beta_mean_difference  
mean_beta_M_tbl["cg24426691_TC11", "Beta_mean_difference"]     
# Correct way 2:
mean_beta_M_tbl["cg06428163_BC21",]$Beta_mean_difference      
mean_beta_M_tbl$Beta_mean_difference[rownames(mean_beta_M_tbl) == "cg06428163_BC21"] '
#--------------------------------------------------------