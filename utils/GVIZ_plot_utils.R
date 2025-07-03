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

# cpgisland granges should be done once at the top of app server



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
dmr_list <- readRDS("C:/Users/ghaza/Documents/ghazal/Bioinformatik_Fächer/Masterarbeit_Project/Scripts/R_Scripts/intermediate_data/DMRs_cutoff_neg0.15_to_0.15_B0_2025-06-24.rds")
dmrs_gr <- create_gr_dmrs(dmr_list$dmr_table)

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

#test function
offtargets_combined <- readRDS("C:/Users/ghaza/Documents/ghazal/Bioinformatik_Fächer/Masterarbeit_Project/Scripts/R_Scripts/modules/intermediate_data/guide1_guide2_guide3_guide4_guide5_guide6.rds")
gr_offtargets <- create_gr_offtargets(offtargets_combined)

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

#test funcgion
#library(readr)
#chromosome_name <- "chr13"
#subtad_file_dir <- "C:/Users/ghaza/Documents/ghazal/Bioinformatik_Fächer/Masterarbeit_Project/Scripts/PythonProject/TADs_CAKI2/3_TADs_as_txt_and_Excel/CAKI2_chr13_25kb_TADs.txt"
#SUBTADs <- read_delim(subtad_file_dir, delim = "\t", show_col_types = FALSE)
#gr_SUBTADs <- create_gr_SUBTADs(chromosome_name, SUBTADs)

#tad_file_dir <- "C:/Users/ghaza/Documents/ghazal/Bioinformatik_Fächer/Masterarbeit_Project/Scripts/PythonProject/TADs_CAKI2/3_TADs_as_txt_and_Excel/CAKI2_chr13_25kb_SubTADs_noDuplicates.txt"
#tads <- read_delim(tad_file_dir, delim = "\t", show_col_types = FALSE)
#gr_tads <- create_gr_SUBTADs(chromosome_name, tads)
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
                          tissue = NULL) {
  
  # Initialize
  itrack <- gtrack <- sondentrack <- cpgIslandTrack <- DmrTrack <- NULL
  offtargetTrack <- betaTrack <- geneTrack <- SUBTadTrack <- TadTrack <- NULL
  
  # Base tracks
  tryCatch({
    itrack <- IdeogramTrack(genome = genome, chromosome = chr)
    gtrack <- GenomeAxisTrack()
  }, error = function(e) {
    message("Error creating base tracks: ", e$message)
  })
  
  # CpG Probes
  tryCatch({
    if (!is.null(gr_cpgs)) {
      sondentrack <- AnnotationTrack(gr_cpgs, chromosome = chr, genome = genome,
                                     name = "CpG Probes", col = NA, fill = "blue",
                                     id = gr_cpgs$CpGs, showFeatureId = FALSE)
    }
  }, error = function(e) {
    message("Error creating CpG Probes track: ", e$message)
  })
  
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
  
  # Beta values
  tryCatch({
    if (!is.null(gr_cpgs) && !is.null(mcols(gr_cpgs)) && !is.null(num_samples)) {
      beta_values <- as.data.frame(mcols(gr_cpgs)[, tail(seq_along(mcols(gr_cpgs)), num_samples)])
      betaTrack <- DataTrack(gr_cpgs, chromosome = chr, genome = genome,
                             name = "Beta values", data = beta_values)
    }
  }, error = function(e) {
    message("Error creating Beta values track: ", e$message)
  })
  
  # Genes
  tryCatch({
    if (!is.null(tx_gr_filtered)) {
      geneTrack <- AnnotationTrack(tx_gr_filtered, chromosome = chr,
                                   name = "Gene Ens", id = tx_gr_filtered$gene_name)
    }
  }, error = function(e) {
    message("Error creating Gene track: ", e$message)
  })
  
  
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

# -----------------
# gene granges which is global
library(EnsDb.Hsapiens.v86)
edb <- EnsDb.Hsapiens.v86
options(ucscChromosomeNames = TRUE)
tx_gr <- genes(edb)
tx_gr_filtered_global <- keepSeqlevels(tx_gr, standardChromosomes(tx_gr), pruning.mode = "coarse")
seqlevelsStyle(tx_gr_filtered_global) <- "UCSC"
# ------------------
'tracks <- create_tracks(genome = "hg38", chr = "chr13",
                        gr_cpgs = cpgs_gr,
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

# then here plot the available tracks
plotGvizTracks <- function(tracks, from, to) {
  # Extract tracks from named list, or NULL if not present
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
  
  # Apply display parameters safely if tracks exist
  if (!is.null(DmrTrack)) displayPars(DmrTrack) <- list(groupAnnotation = "feature")
  if (!is.null(cpgIslandTrack)) displayPars(cpgIslandTrack) <- list(groupAnnotation = "feature")
  if (!is.null(cpgsTrack)) displayPars(cpgsTrack) <- list(stacking = "dense", showFeatureId = FALSE)
  if (!is.null(geneTrack)) displayPars(geneTrack) <- list(featureAnnotation = "id", fontcolor.feature = "black", stacking = "full", cex = 0.7)
  if (!is.null(offtargetTrack)) displayPars(offtargetTrack) <- list(groupAnnotation = "feature", fontcolor.group = "black")
  #if (!is.null(SUBTadTrack)) displayPars(cpgsTrack) <- list(stacking = "dense", showFeatureId = FALSE) # and add thicker borders of boxes if stacked
  #if (!is.null(TadTrack)) displayPars(cpgsTrack) <- list(stacking = "dense", showFeatureId = FALSE)
  
  
  
  # Prepare list of tracks to plot (exclude NULLs)
  to_plot <- Filter(Negate(is.null), list(
    itrack, gtrack, cpgsTrack, cpgIslandTrack, DmrTrack,
    offtargetTrack, betaTrack, geneTrack, SUBTadTrack, TadTrack
  ))
  
  # Plot tracks with given genomic range
  plotTracks(
    to_plot,
    from = from, to = to,
    showBandId = TRUE, cex.bands = 0.7,
    groups = rep(c("unguided", "guided"), each = 4),
    type = c("a", "p"),
    transcriptAnnotation = "transcript",
    cex.group = 0.9,
    cex.legend = 0.9,
    cex.title = 0.75,
    col.title = "black"
  )
}


# test function
#plotGvizTracks(tracks = tracks, from = 95551415, to = 95554041)
#----------------------------------------------------------------------------------------

# utils_genome_helpers.R

# Load once at top of the app or where utilities are loaded
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)

# Generate chromosome size table for hg38
hg38 <- BSgenome.Hsapiens.UCSC.hg38
chr_lengths <- seqlengths(hg38)

chr_size_df <- data.frame(
  Chromosome = names(chr_lengths),
  Length = as.numeric(chr_lengths),
  stringsAsFactors = FALSE
)

# Helper function to get chromosome size
get_chr_max <- function(chr_name) {
  if (chr_name %in% chr_size_df$Chromosome) {
    chr_size_df$Length[chr_size_df$Chromosome == chr_name]
  } else {
    Inf  # fallback if chromosome name is invalid
  }
}
