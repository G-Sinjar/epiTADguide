# preprocessing_utils.R
#Author: Ghazal Sinjar
#Date: 16.06.2025


#' Filter probes based on detection p-values
#'
#' @param normalized_object A matrix or minfi object (e.g., GenomicRatioSet)
#' @param RGset The raw RGChannelSet used to compute detection p-values
#' @return A filtered version of `normalized_object` with low-quality probes removed
#' @export



# ───────────────────────────────────────────────────────────────────────
# Filtering function: Detection P-Value Filtering
# Filters out probes with detection p-value > 0.01 in any sample
# This function handles both standard matrices and minfi S4 objects
filter_by_detectionP <- function(normalized_object, RGset) {
  print(paste("Class of RGset in filter_by_detectionP:", class(RGset)[1])) 
  detectp <- detectionP(RGset)
  unreliable_probes_mask <- rowSums(detectp > 0.01) > 0
  unreliable_probe_names <- rownames(detectp)[unreliable_probes_mask]
  normalized_object_probe_names <- rownames(normalized_object)
  
  # Check if normalized object is a minfi object (like GenomicRatioSet)
  is_minfi_S4_object <- inherits(normalized_object, "GenomicRatioSet") ||
    inherits(normalized_object, "MethylSet") ||
    inherits(normalized_object, "MSet") ||
    inherits(normalized_object, "RGChannelSet")
  
  # Use direct subsetting only if probe names match exactly
  can_do_direct_subset <- !is_minfi_S4_object && isTRUE(all.equal(normalized_object_probe_names, rownames(detectp)))
  
  if (can_do_direct_subset) {
    message("Applying detection P-value filtering using direct logical subsetting (row names match).")
    # Direct subsetting, assuming the order and presence of probes matches detectp
    filtered_object <- normalized_object[!unreliable_probes_mask, ]
  } else {
    message("Applying detection P-value filtering using probe name intersection (likely S4 object or mismatch).")
    message("This might take a moment if the object is large or probe names mismatch extensively.")
    # This is the more robust path, using setdiff for name-based filtering
    probes_to_keep <- setdiff(normalized_object_probe_names, unreliable_probe_names)
    filtered_object <- normalized_object[probes_to_keep, ]
  }
  
  message(paste0("Original probes: ", nrow(normalized_object),
                 ". Probes after detection P-value filtering: ", nrow(filtered_object), "."))
  
  return(filtered_object)
}
# ───────────────────────────────────────────────────────────────────────

