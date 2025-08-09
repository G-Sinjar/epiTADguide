library(readr)
library(stringr)

# ─────────────────────────────────────
# Main Shiny App (app.R equivalent)
# ─────────────────────────────────────
source("../utils/GVIZ_plot_utils.R")
source("./10_Gviz_plot_Module_v1.R")

# This is a static object that will be passed to the DMR module
message("Loading/Preparing tx_gr_filtered for annotation...")
edb <- EnsDb.Hsapiens.v86
options(ucscChromosomeNames = TRUE)
tx_gr <- genes(edb)
tx_gr_filtered_global <- keepSeqlevels(tx_gr, standardChromosomes(tx_gr), pruning.mode = "coarse")
seqlevelsStyle(tx_gr_filtered_global) <- "UCSC"
message("tx_gr_filtered prepared.")
#-------------------------------------------------------------------------
# Chromosome Lengths Table from BSgenome
library(BSgenome.Hsapiens.UCSC.hg38)
hg38 <- BSgenome.Hsapiens.UCSC.hg38
chr_lengths <- seqlengths(hg38)
chr_size_df_global <- data.frame(
  Chromosome = names(chr_lengths),
  Length = as.numeric(chr_lengths),
  stringsAsFactors = FALSE
)
print("Global: chr_size_df_global loaded.")
#------------------------------------------------------------------
# CpG Islands - Loaded once globally
# Uses the loadCpGIslands_gr function from utils/GVIZ_plot_utils.R
gr_cpgIslands_global <- loadCpGIslands_gr(destfile = "cpgIslandExt_hg38.txt.gz")
print(paste("Global: gr_cpgIslands_global loaded. Number of CpG Islands:", length(gr_cpgIslands_global)))
#------------------------------

dmr_list <- reactiveVal(NULL)
annotated <- reactiveVal(NULL)
offtargets_combined <- reactiveVal(NULL)
tadcalling_results_dummy <- reactiveVal(NULL)
print("Main App: Initializing data loading...")

# Replace with your actual file paths
dmr_list_path <- "../modules/main_app_tests/df/DMR_results/DMRs_unguided_vs_sample_groupguided_cutoff_-0.15_0.15_B0_20250804.rds"
annotated_path <- "../intermediate_data/annotated_object_20250627.rds"
annotated_tbl_withqval <- "../main_app_tests/dmr_ram_detect/intermediate_data/annotated_with_qval.rds"
offtargets_path <- "./intermediate_data/guide1_guide2_guide3_guide4_guide5_guide6.rds"
tad_subtad_path <- "../main_app_tests/dmr_ram_detect/TADcaller_results/4DNFIIH3SM5N/TADcaller_Results/TADs_CAKI2/processed_tads"


# Create a list of autosomal chromosomes (1 to 22)
autosomal_chrs <- paste0("chr", 1:22)
sex_chrs <- c("chrX", "chrY")
all_chrs <- c(autosomal_chrs, sex_chrs)
#print(all_chrs)



results_dmr <- readRDS(dmr_list_path)
results_anno <- readRDS(annotated_path)
results_ano_qval <- readRDS(annotated_tbl_withqval)

dmr_list <- reactive({
  list(dmr_table = reactiveVal(results_dmr$dmr_table),
       pheno = reactiveVal(results_dmr$pheno),
       ref_group = reactiveVal("unguided")
       )})

annotated <- reactive({ list(
  annotated_table_with_qval = reactiveVal(results_anno$annotated_table)
)
})
annotated_qval <- reactive({ list(
  annotated_table_with_qval = reactiveVal(results_ano_qval)
) })

if (file.exists(offtargets_path)) {
  offtargets_combined(reactiveVal(readRDS(offtargets_path)))
  print(paste("Main App: Off-targets file loaded from:", offtargets_path))
} else {
  print(paste("Main App: ERROR - Off-targets file not found:", offtargets_path))
}

if (file.exists(tad_subtad_path)) {
  print("Main App: Tads_SubTADs path excists")
} else {
  print("Main App: ERROR - STads_SubTADs path doesnt excists")
}
print("Main App: Data loading complete.")

tadcalling_results_dummy(list(
  current_tissue = reactiveVal(toString("CAKI2")),
  current_resolution = reactiveVal(25),
  processed_output_dir= reactiveVal(tad_subtad_path),
  processed_chroms_list_rv = reactiveVal(all_chrs)
))

ui <- page_navbar(
  title = "EPIC Array Pipeline",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  nav_panel("Final Plot", GvizPlotUI("gviz_test_id"))
)

server <- function(input, output, session) {
  print("Main Server started.")
  
  GvizPlotServer(
    id = "gviz_test_id",
    dmr_results = dmr_list,
    annotation_results = annotated,
    offtarget_table = offtargets_combined,
    tadcalling_results = tadcalling_results_dummy,
    chr_size_df = chr_size_df_global,
    tx_gr_filtered = tx_gr_filtered_global,
    gr_cpgIslands = gr_cpgIslands_global
  )
  print("Main Server: GvizPlotServer module called.")
}

shinyApp(ui, server)