# GVIZ_p
library(GenomicRanges)
library(IRanges)  # IRanges is usually loaded with GenomicRanges but can be loaded explicitly
library(readxl)
library(dplyr)
library(GenomicRanges)
library(Gviz)
library(shiny)
library(bslib)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)



# 5 function to convert the 5 tables to granges
#-----------------------------------------------------------------------------
# function 1: create granges of DMR table
create_gr_dmrs <- function(DMR_table) {
  strand_col <- if ("strand" %in% colnames(DMR_table)) "strand" else NULL
  
  gr <- makeGRangesFromDataFrame(DMR_table, seqnames.field = "chr",
                                 start.field = "start", end.field = "end",
                                 strand.field = strand_col,
                                 keep.extra.columns = TRUE)
  
  mcols(gr)$p.value <- if ("p.value" %in% colnames(mcols(gr))) 
    as.numeric(as.character(mcols(gr)$p.value)) else NA_real_
  
  gr
}
# test function
#dmr_list <- readRDS("../modules/main_app_tests/epic-test/DMR_results/DMRs_cutoff_-0.15_0.15_B1_20250801.rds")
#dmrs_gr <- create_gr_dmrs(dmr_list$dmr_table)

#--------------------------------------------------------------------------------------------------------------
# function 2: create granges of annotation table
create_gr_cpgs <- function(cpg_table) {
  gr <- GRanges(
    seqnames = cpg_table$chr,
    ranges = IRanges(start = cpg_table$pos, end = cpg_table$pos)
  )
  mcols(gr) <- cpg_table[, 4:ncol(cpg_table), drop = FALSE]
  return(gr)
}
## test function 2
#annotated <- readRDS("C:/Users/ghaza/Documents/ghazal/Bioinformatik_Fächer/Masterarbeit_Project/Scripts/R_Scripts/intermediate_data/annotated_object_20250627.rds")
#cpgs_gr <- create_gr_cpgs(annotated$annotated_table)

# case 2: with qval
#annotated_tbl_withqval <- readRDS("../modules/main_app_tests/box_pheno_pass/intermediate_data/annotated_with_qval.rds")
#cpgs_gr_qval <- create_gr_cpgs(annotated_tbl_withqval)
# manipulate the data 
#cpgs_gr_qval$significance_last_qvalue[cpgs_gr_qval$Name== "cg16556145_BC11"] <- "Sig"

#--------------------------------------------------------------------------------------------
# function 3: create granges of offtargets table

create_gr_offtargets <- function(offtarget_table) {
  gr <- GRanges(
    seqnames = offtarget_table$Chrom,
    ranges = IRanges(
      start = offtarget_table$Start,
      end = offtarget_table$End
    )
  )
  
  mcols(gr) <- offtarget_table[, !colnames(offtarget_table) %in% 
                                 c("Chrom", "Strand", "Start", "End"), drop = FALSE]
  
  seqlevelsStyle(gr) <- "UCSC"
  
  return(gr)
}

'#test function
offtargets_combined <- readRDS("C:/Users/ghaza/Documents/ghazal/Bioinformatik_Fächer/Masterarbeit_Project/Scripts/R_Scripts/modules/intermediate_data/guide1_guide2_guide3_guide4_guide5_guide6.rds")
gr_offtargets <- create_gr_offtargets(offtargets_combined)
'
#-----------------------------------------------------------------------------------
# function 4: create granges of tads/subtads table
create_gr_TADs_SUBTADs <- function(chr, SUBTADs) {
  gr <- GRanges(
    seqnames = chr,
    ranges = IRanges(start = SUBTADs$StartPos, end = SUBTADs$EndPos)
  )
  
  # Exclude StartPos and EndPos, keep other columns as metadata
  meta_cols <- SUBTADs[, !colnames(SUBTADs) %in% c("StartPos", "EndPos"), drop = FALSE]
  
  mcols(gr) <- meta_cols
  
  return(gr)
}

'#setwd("C:/Users/ghaza/Documents/ghazal/Bioinformatik_Fächer/Masterarbeit_Project/Scripts/R_Scripts/utils")
#test funcgion
library(readr)
chromosome_name <- "chr13"
subtad_file_dir <- "../main_app_tests/epic-test/TADcaller_results/4DNFIIH3SM5N/TADcaller_Results/TADs_CAKI2/processed_tads/CAKI2_chr13_25kb_SubTADs.txt"
SUBTADs <- read_delim(subtad_file_dir, delim = "\t", show_col_types = FALSE)
gr_SUBTADs <- create_gr_TADs_SUBTADs(chromosome_name, SUBTADs)

tad_file_dir <- "../main_app_tests/epic-test/TADcaller_results/4DNFIIH3SM5N/TADcaller_Results/TADs_CAKI2/processed_tads/CAKI2_chr13_25kb_TADs.txt"
tads <- read_delim(tad_file_dir, delim = "\t", show_col_types = FALSE)
gr_tads <- create_gr_TADs_SUBTADs(chromosome_name, tads)'
#----------------------------------------------------------------------------------------
# function 5: create granges of cpgIslands
loadCpGIslands_gr <- function(destfile = "cpgIslandExt_hg38.txt.gz") {
  if (!file.exists(destfile)) {
    url <- "http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cpgIslandExt.txt.gz"
    download.file(url, destfile)
  }
  cpg_Islands <- read.table(gzfile(destfile), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(cpg_Islands) <- c("bin", "chrom","chromStart","chromEnd","name","length","cpgNum","gcNum","perCpg","perGc", "obsExp")
  
  gr <- GRanges(
    seqnames = cpg_Islands$chrom,
    ranges = IRanges(start = cpg_Islands$chromStart + 1, end = cpg_Islands$chromEnd)
  )
  mcols(gr) <- cbind(feature = cpg_Islands$name, cpg_Islands[, 7:ncol(cpg_Islands)])
  return(gr)
}

# test function
# Then somewhere globally (e.g. global.R or your app initialization)
#gr_cpgIslands <- loadCpGIslands_gr()



#----------------------------------------------------------------------------------------
# create_tracks function
# so reactive values are chr() from() to() the range() between them
# after creating the granges  a function that creates the tracks of only the granges that are available so if offtargets arent there so no track
create_tracks <- function(genome, chr,
                          gr_cpgs = NULL,
                          gr_cpgIslands = NULL,
                          dmrs_gr = NULL,
                          gr_offtargets = NULL,
                          tx_gr_filtered = NULL,
                          gr_SUBTADs = NULL,
                          gr_TADs = NULL,
                          binsize = NULL,
                          num_samples = NULL,
                          tissue = NULL,
                          pheno_data = NULL) {
  
  if (!is.null(gr_cpgs)) {
    if (!inherits(gr_cpgs, "GRanges")) {
      message("gr_cpgs is not a GRanges object")
      gr_cpgs <- NULL
    }
  }  # Initialize
  itrack <- gtrack <- sondentrack <- cpgIslandTrack <- DmrTrack <- NULL
  offtargetTrack <- betaTrack <- geneTrack <- SUBTadTrack <- TadTrack <- NULL
  
  # Base tracks
  tryCatch({
    itrack <- IdeogramTrack(genome = genome, chromosome = chr)
    gtrack <- GenomeAxisTrack()
  }, error = function(e) {
    message("Error creating base tracks: ", e$message)
  })
  
  #-------------------------------------
  #message(paste0("class of gr_cpgs", class(gr_cpgs)))
  #message(paste0("class of gr_cpgs()",class(gr_cpgs())))
  # CpG Probes  # V2 Modified CpG Probes 
  tryCatch({
    if (!is.null(gr_cpgs)) {
      # Create basic track
      sondentrack <- AnnotationTrack(
        gr_cpgs,
        chromosome = chr,
        genome = genome,
        col = NA,
        name = "CpG Probes",
        id = gr_cpgs$CpGs,
        showFeatureId = FALSE,
        stacking = "dense"
      )
      
      # Apply coloring using explicit package reference
      if ("significance_last_qvalue" %in% colnames(mcols(gr_cpgs))) {
        # Use Gviz::feature to avoid namespace conflicts
        Gviz::feature(sondentrack) <- ifelse(
          gr_cpgs$significance_last_qvalue == "Sig", 
          "Sig", 
          "Non-Sig"
        )
        
        # Set display parameters
        displayPars(sondentrack) <- list(
          "Sig" = "green",
          "Non-Sig" = "blue"
        )
      } else {
        # No significance column - all blue
        displayPars(sondentrack) <- list(fill = "blue")
      }
    }
  }, error = function(e) {
    message("Error creating CpG Probes track: ", e$message)
  })
  #------------------------------------------------
  # CpG Islands
  tryCatch({
    if (!is.null(gr_cpgIslands)) {
      cpgIslandTrack <- AnnotationTrack(gr_cpgIslands, chromosome = chr, genome = genome,
                                        name = "CpG Islands", fill = "darkgreen",
                                        feature = mcols(gr_cpgIslands)$feature)
    }
  }, error = function(e) {
    message("Error creating CpG Islands track: ", e$message)
  })
  #----------------------
  # DMRs
  tryCatch({
    if (!is.null(dmrs_gr)) {
      DmrTrack <- AnnotationTrack(dmrs_gr, genome = genome, chromosome = chr,
                                  name = "DMRs", col = "orange", fill = "orange",
                                  feature = paste(dmrs_gr$DMR_ID, " P.val=", round(dmrs_gr$p.value, 5)))
    }
  }, error = function(e) {
    message("Error creating DMRs track: ", e$message)
  })
  #----------------------------
  # Off-targets
  tryCatch({
    if (!is.null(gr_offtargets)) {
      offtargetTrack <- AnnotationTrack(gr_offtargets, genome = genome, chromosome = chr,
                                        name = "potential Off-targets",
                                        feature = gr_offtargets$ID, col = NA, fill = "red")
    }
  }, error = function(e) {
    message("Error creating Off-targets track: ", e$message)
  })
  #--------------------------------------
  # Beta values # V2:
  # Beta values - modified version with dynamic column selection
  tryCatch({
    if (!is.null(gr_cpgs) && !is.null(mcols(gr_cpgs)) && !is.null(num_samples)) {
      # Get all metadata column names
      all_cols <- colnames(mcols(gr_cpgs))
      
      # Check if q-value and significance columns exist
      has_qval_sig <- all(c("qval", "significance_last_qvalue") %in% all_cols)
      
      # Determine beta value columns
      if (has_qval_sig) {
        # When qval/sig columns exist, beta values are the columns before them
        beta_cols <- all_cols[(length(all_cols) - num_samples - 1):(length(all_cols) - 2)]
      } else {
        # Otherwise just take the last num_samples columns
        beta_cols <- tail(all_cols, num_samples)
      }
      
      # Extract beta values
      beta_values <- as.data.frame(mcols(gr_cpgs)[, beta_cols, drop = FALSE])
      
      # Create DataTrack
      betaTrack <- DataTrack(gr_cpgs, 
                             chromosome = chr, 
                             genome = genome,
                             name = "Beta values", 
                             data = beta_values)
    }
  }, error = function(e) {
    message("Error creating Beta values track: ", e$message)
  })
  #-------------------------------
  # Genes
  tryCatch({
    if (!is.null(tx_gr_filtered)) {
      geneTrack <- AnnotationTrack(tx_gr_filtered, chromosome = chr,
                                   name = "Gene Ens", id = tx_gr_filtered$gene_name)
    }
  }, error = function(e) {
    message("Error creating Gene track: ", e$message)
  })
  #----------------------------
  # SUBTADs
  tryCatch({
    if (!is.null(gr_SUBTADs)) {
      subtad_name <- paste0("SUBTADs\n", tissue, if (!is.null(binsize)) paste0("\n", binsize, "kb") else "")
      SUBTadTrack <- AnnotationTrack(gr_SUBTADs, genome = genome, chromosome = chr,
                                     name = subtad_name,
                                     #feature = paste0(gr_SUBTADs$Parent_TAD_ID, "_", gr_SUBTADs$SubTAD_ID)
                                     col = "black", fill = "pink")
    }
  }, error = function(e) {
    message("Error creating SUBTADs track: ", e$message)
  })
  #-------------------------------
  # TADs
  tryCatch({
    if (!is.null(gr_TADs)) {
      tad_name <- paste0("TADs\n", tissue, if (!is.null(binsize)) paste0("\n", binsize, "kb") else "")
      TadTrack <- AnnotationTrack(gr_TADs, genome = genome, chromosome = chr,
                                  name = tad_name, 
                                  #feature = gr_TADs$TAD_ID, 
                                  col = "black", fill = "purple")
    }
  }, error = function(e) {
    message("Error creating TADs track: ", e$message)
  })
  
  # Return all tracks as a named list
  return(list(
    itrack = itrack,
    gtrack = gtrack,
    cpgsTrack = sondentrack,
    cpgIslandTrack = cpgIslandTrack,
    DmrTrack = DmrTrack,
    offtargetTrack = offtargetTrack,
    betaTrack = betaTrack,
    geneTrack = geneTrack,
    SUBTadTrack = SUBTadTrack,
    TadTrack = TadTrack
  ))
}

################## test function
# gene granges which is global
'library(EnsDb.Hsapiens.v86)
edb <- EnsDb.Hsapiens.v86
options(ucscChromosomeNames = TRUE)
tx_gr <- genes(edb)
tx_gr_filtered_global <- keepSeqlevels(tx_gr, standardChromosomes(tx_gr), pruning.mode = "coarse")
seqlevelsStyle(tx_gr_filtered_global) <- "UCSC"
# ------------------
tracks <- create_tracks(genome = "hg38", chr = "chr13",
                        gr_cpgs = cpgs_gr,   # or cpgs_gr
                        gr_cpgIslands = gr_cpgIslands,
                        dmrs_gr = dmrs_gr,
                        gr_offtargets = gr_offtargets,
                        tx_gr_filtered = tx_gr_filtered_global,
                        gr_SUBTADs = gr_SUBTADs,
                        gr_TADs = gr_tads,
                        binsize = 25,
                        num_samples = 8,
                        tissue = "CAKI2")'

#---------------------------------------------------------------------------------------------------
# plotGvizTracks function
plotGvizTracks <- function(tracks, from, to, pheno = NULL, gr_cpgs = NULL, ref_group = NULL) {
  
  # Input validation
  if (!is.null(ref_group)) {
    stopifnot(is.character(ref_group) || is.factor(ref_group))
    ref_group <- as.character(ref_group)
  }
  
  # Extract tracks from named list
  itrack <- tracks$itrack
  gtrack <- tracks$gtrack
  cpgsTrack <- tracks$cpgsTrack
  cpgIslandTrack <- tracks$cpgIslandTrack
  DmrTrack <- tracks$DmrTrack
  offtargetTrack <- tracks$offtargetTrack
  betaTrack <- tracks$betaTrack
  geneTrack <- tracks$geneTrack
  SUBTadTrack <- tracks$SUBTadTrack
  TadTrack <- tracks$TadTrack
  
  # Set display parameters for all tracks
  if (!is.null(itrack)) displayPars(itrack) <- list(col = "black", fontcolor = "black")
  if (!is.null(gtrack)) displayPars(gtrack) <- list(col = "black", fontcolor = "black")
  
  if (!is.null(cpgsTrack)) {
    if (inherits(cpgsTrack, "AnnotationTrack")) {
      track_values <- values(cpgsTrack)
      sig_features <- which(track_values$feature == "Sig")
    } else {
      sig_features <- which(Gviz::feature(cpgsTrack) == "Sig")
    }
    
    if (length(sig_features) > 0) {
      displayPars(cpgsTrack) <- list(
        "Sig" = "green",
        "Non-Sig" = "blue",
        fill = ifelse(track_values$feature == "Sig", "green", "blue"),
        showFeatureId = FALSE,
        stacking = "dense"
      )
    }
  }
  if (!is.null(cpgIslandTrack)) displayPars(cpgIslandTrack) <- list( groupAnnotation = "feature")
  if (!is.null(DmrTrack)) displayPars(DmrTrack) <- list(groupAnnotation = "feature")
  if (!is.null(offtargetTrack)) displayPars(offtargetTrack) <- list(groupAnnotation = "feature", fontcolor.group = "black")
  if (!is.null(geneTrack)) displayPars(geneTrack) <- list(col = "black", fontcolor = "black", featureAnnotation = "id", fontcolor.feature = "black", stacking = "full", cex = 0.7)
  
  group_labels <- NULL
  group_colors <- NULL
  
  # --- CORRECTED LOGIC FOR SAMPLE AND GROUP MAPPING ---
  if (!is.null(pheno) && !is.null(gr_cpgs)) {
    if ("Sample_Group" %in% colnames(pheno)) {
      
      # 1. Get the sample IDs from the pheno data (this is the correct order)
      pheno_sample_ids <- rownames(pheno)
      
      # 2. Get all metadata column names from gr_cpgs
      mcols_names <- colnames(mcols(gr_cpgs))
      
      # 3. Filter mcols_names to only include columns that are also in pheno_sample_ids
      beta_colnames <- mcols_names[mcols_names %in% pheno_sample_ids]
      
      # 4. Reorder the beta_colnames to match the order of pheno_sample_ids
      beta_colnames_reordered <- pheno_sample_ids[match(beta_colnames, pheno_sample_ids)]
      
      # 5. Remove any NA values that may have resulted from the match
      #beta_colnames_reordered <- beta_colnames_reordered[!is.na(beta_colnames_reordered)]
      
      # 6. Verify that the reordering worked as expected (for debugging)
      if (length(beta_colnames_reordered) > 0) {
        message("✅ Sample order verification passed after reordering.")
      } else {
        warning("No matching beta value columns found for the given pheno data.")
      }
      
      # 7. Assign group labels based on the corrected order
      if (length(beta_colnames_reordered) > 0) {
        group_labels <- pheno[beta_colnames_reordered, "Sample_Group", drop = TRUE]
        
        if (!is.factor(group_labels)) {
          group_labels <- factor(group_labels)
        }
        
        '# --- Robust color assignment ---
        unique_groups <- levels(group_labels)
        if (length(unique_groups) == 2) {
          group_colors <- setNames(rep("", length(unique_groups)), unique_groups)
          if (!is.null(ref_group) && ref_group %in% unique_groups) {
            group_colors[as.character(ref_group)] <- "blue"
            group_colors[setdiff(unique_groups, ref_group)] <- "orange"
          } else {
            warning("Reference group not found or not provided. Using default blue/orange.")
            group_colors[] <- c("blue", "orange")[seq_along(unique_groups)]
          }
        } else if (length(unique_groups) > 2) {
          group_colors <- setNames(rainbow(length(unique_groups)), unique_groups)
        } else if (length(unique_groups) == 1) {
          group_colors <- setNames("darkblue", unique_groups)
        }
        message("Group levels: ", paste(unique_groups, collapse = ", "))
        print(data.frame(Group = names(group_colors), Color = group_colors))'
      }
      
    } else {
      warning("Column 'Sample_Group' not found in pheno.")
    }
  }
  # --- END OF CORRECTED LOGIC ---
  
  # Prepare tracks to plot
  to_plot <- Filter(Negate(is.null), list(
    itrack, gtrack, cpgsTrack, cpgIslandTrack, DmrTrack,
    offtargetTrack, betaTrack, geneTrack, SUBTadTrack, TadTrack
  ))
  
  if (length(to_plot) == 0) {
    stop("No valid tracks to plot")
  }
  
  message("Plot range: ", from, "-", to)
  
  tryCatch({
    plotTracks(
      to_plot,
      from = from, to = to,
      showBandId = TRUE,
      cex.bands = 0.7,
      groups = group_labels,
      #col = group_colors,
      type = c("a", "p"),
      transcriptAnnotation = "transcript",
      cex.group = 0.9,
      cex.legend = 0.9,
      cex.title = 0.75,
      col.title = "black"
    )
  }, error = function(e) {
    message("Plotting failed with error: ", e$message)
    message("Track details:")
    print(lapply(to_plot, function(x) c(class(x), length(x))))
    stop(e)
  })
}


#pheno <- dmr_list$pheno
#cpgs_gr
#from <- 95433599
#to <- 95579959
## test function
#plotGvizTracks(tracks = tracks, from = 95433599, to = 95579959, pheno = pheno, gr_cpgs=cpgs_gr, ref_group= "unguided")
#----------------------------------------------------------------------------------------

#' Find Chromosomes with TAD and SubTAD Tables in a Given Path and Check for a Specific Chromosome
#'
#' This function scans a specified directory for TAD and SubTAD files
#' based on a predefined naming convention and returns a unique list
#' of chromosomes for which these files exist. It also checks if a
#' provided chromosome name (expected to be already normalized, e.g., "chr1", "chrX")
#' is present in this list. If the chromosome is found, it also returns the
#' basenames (filenames only) of the corresponding TAD and SubTAD files.
#'
#' The expected file naming convention is:
#' "TISSUE_chrCHROMOSOME_RESOLUTIONkb_TADs.txt"
#' "TISSUE_chrCHROMOSOME_RESOLUTIONkb_SubTADs.txt"
#'
#' @param tissue Character string. The tissue name (e.g., "CAKI2").
#' @param resolution Numeric. The resolution in kilobases (e.g., 25).
#' @param chr Character string. The chromosome name to check for.
#'   This argument is expected to be already normalized (e.g., "chr1", "chrX", "chrY", "chrM").
#' @param processed_data_path Character string. The full path to the directory
#'   containing the processed TAD/SubTAD files.
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{available_chrs}: A character vector of unique chromosome names
#'       (e.g., "chr1", "chr13") for which both TAD and SubTAD files exist with
#'       the specified tissue and resolution. Returns an empty vector if no
#'       matching files are found.
#'     \item \code{chr_in_list}: A logical value (TRUE/FALSE) indicating whether
#'       the provided \code{chr} argument is present in the \code{available_chrs} list.
#'     \item \code{matched_tad_filename}: Character string. The basename (filename only)
#'       of the matching TAD file for the given \code{chr}, or \code{NULL} if not found.
#'     \item \code{matched_subtad_filename}: Character string. The basename (filename only)
#'       of the matching SubTAD file for the given \code{chr}, or \code{NULL} if not found.
#'   }
#'
#' @examples
#' # Create a dummy directory and files for demonstration
#' temp_dir <- tempdir()
#' dummy_path <- file.path(temp_dir, "test_tads_check_basenames")
#' dir.create(dummy_path, showWarnings = FALSE)
#'
#' # Create dummy files
#' file.create(file.path(dummy_path, "CAKI2_chr1_25kb_TADs.txt"))
#' file.create(file.path(dummy_path, "CAKI2_chr1_25kb_SubTADs.txt")) # Note: noDuplicates removed
#' file.create(file.path(dummy_path, "CAKI2_chr13_25kb_TADs.txt"))
#' file.create(file.path(dummy_path, "CAKI2_chr13_25kb_SubTADs.txt")) # Note: noDuplicates removed
#'
#' # Run the function with a normalized chromosome input
#' results1 <- find_existing_tads_subtads_chrs_with_check(
#'   tissue = "CAKI2",
#'   resolution = 25,
#'   chr = "chr1", # Already normalized
#'   processed_data_path = dummy_path
#' )
#' print(results1)
#'
#' results2 <- find_existing_tads_subtads_chrs_with_check(
#'   tissue = "CAKI2",
#'   resolution = 25,
#'   chr = "chrX", # Already normalized
#'   processed_data_path = dummy_path
#' )
#' print(results2)
#'
#' # Clean up dummy files
#' unlink(dummy_path, recursive = TRUE)
#'
find_existing_tads_subtads_chrs_with_check <- function(tissue, resolution, chr, processed_data_path) {
  
  # Validate processed_data_path
  if (!dir.exists(processed_data_path)) {
    stop(paste("Directory not found:", processed_data_path, "Please provide a valid path."))
  }
  
  # The 'chr' argument is assumed to be already normalized (e.g., "chr1", "chrX")
  normalized_input_chr <- chr
  
  # Construct patterns for TAD and SubTAD files
  tad_pattern <- paste0("^", tissue, "_chr(.+)_", resolution, "kb_TADs\\.txt$")
  # SubTAD pattern updated: removed '_noDuplicates' as per user confirmation
  subtad_pattern <- paste0("^", tissue, "_chr(.+)_", resolution, "kb_SubTADs\\.txt$")
  
  # List files matching the TAD pattern, getting full names so we can extract base names later
  tad_files_full_paths <- list.files(
    path = processed_data_path,
    pattern = tad_pattern,
    full.names = TRUE
  )
  tad_files_base_names <- basename(tad_files_full_paths)
  
  # List files matching the SubTAD pattern, getting full names
  subtad_files_full_paths <- list.files(
    path = processed_data_path,
    pattern = subtad_pattern,
    full.names = TRUE
  )
  subtad_files_base_names <- basename(subtad_files_full_paths)
  
  # Function to extract chromosome from filename (basename) and prepend "chr"
  extract_chr_from_filename <- function(filename, pattern) {
    chr_match <- sub(pattern, "\\1", filename)
    if (chr_match != filename && chr_match != "") {
      if (toupper(chr_match) %in% c("X", "Y", "M")) {
        return(paste0("chr", toupper(chr_match)))
      } else {
        return(paste0("chr", chr_match))
      }
    } else {
      return(NA_character_)
    }
  }
  
  # Create a named vector mapping extracted chromosome (name) to full path (value) for TADs
  chrs_from_tads_mapping <- setNames(
    tad_files_full_paths,
    sapply(tad_files_base_names, extract_chr_from_filename, pattern = tad_pattern, USE.NAMES = FALSE)
  )
  chrs_from_tads_mapping <- chrs_from_tads_mapping[!is.na(names(chrs_from_tads_mapping))]
  chrs_from_tads <- unique(names(chrs_from_tads_mapping))
  
  # Create a named vector mapping extracted chromosome (name) to full path (value) for SubTADs
  chrs_from_subtads_mapping <- setNames(
    subtad_files_full_paths,
    sapply(subtad_files_base_names, extract_chr_from_filename, pattern = subtad_pattern, USE.NAMES = FALSE)
  )
  chrs_from_subtads_mapping <- chrs_from_subtads_mapping[!is.na(names(chrs_from_subtads_mapping))]
  chrs_from_subtads <- unique(names(chrs_from_subtads_mapping))
  
  # Find chromosomes that have *both* TAD and SubTAD files
  common_chrs <- base::intersect(chrs_from_tads, chrs_from_subtads)
  
  # Sort the chromosomes for consistent output
  if (requireNamespace("gtools", quietly = TRUE)) {
    sorted_chrs <- gtools::mixedsort(common_chrs)
  } else {
    sorted_chrs <- sort(common_chrs)
  }
  
  # Check if the (already normalized) input chromosome is in the list of available chromosomes
  chr_is_in_list <- normalized_input_chr %in% sorted_chrs
  
  # Initialize file basenames to NULL
  matched_tad_filename <- NULL
  matched_subtad_filename <- NULL
  
  # If the input chromosome is found, get its specific file basenames
  if (chr_is_in_list) {
    # Find the specific TAD file full path for the normalized_input_chr
    matching_tad_full_paths <- chrs_from_tads_mapping[names(chrs_from_tads_mapping) == normalized_input_chr]
    if (length(matching_tad_full_paths) > 0) {
      matched_tad_filename <- basename(unname(matching_tad_full_paths[1]))
    }
    
    # Find the specific SubTAD file full path for the normalized_input_chr
    matching_subtad_full_paths <- chrs_from_subtads_mapping[names(chrs_from_subtads_mapping) == normalized_input_chr]
    if (length(matching_subtad_full_paths) > 0) {
      matched_subtad_filename <- basename(unname(matching_subtad_full_paths[1]))
    }
  }
  
  return(list(
    available_chrs = sorted_chrs,
    chr_in_list = chr_is_in_list,
    matched_tad_filename = matched_tad_filename,
    matched_subtad_filename = matched_subtad_filename
  ))
}
'
# Example Usage (assuming you have your files in the specified directory)
processed_tads_output_dir_rv <- "C:/Users/ghaza/Documents/ghazal/Bioinformatik_Fächer/Masterarbeit_Project/Scripts/PythonProject/test/TADcaller_Results/TADs_CAKI2/processed_tads"

tissue_name <- "CAKI2"
resolution_value <- 25 # Assuming resolution is numeric

# CHANGE THIS LINE:
available_chrs <- find_existing_tads_subtads_chrs_with_check( 
  tissue = tissue_name,
  resolution = resolution_value,
  processed_data_path = processed_tads_output_dir_rv,
  chr = "chr13"
)
print(available_chrs)'