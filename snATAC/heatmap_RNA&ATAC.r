#' @title Transcription Factor ChromVAR Activity & RNA Expression Heatmap Pipeline
#' @description Extracts specified TF motifs and gene matrices across cohorts, 
#' normalizes row scales, and generates multi-panel publication-ready analytical heatmaps.
#' @author Chen Ao
#' @date 2026-06-28

library(Seurat)
library(tidyverse)
library(pheatmap)
library(ClusterGVis)
library(scRNAtoolVis)

# ==========================================
# 1. Configured Target Profiles & Coordinates
# ==========================================

# 1.1 Mammary Gland Target Dictionary
MAMMARY_MOTIFS <- c(
  "MA0519.1", "MA1624.1", "MA0703.2", "MA0808.1", "MA0148.4", 
  "MA0157.2", "MA0687.1", "MA0081.2", "MA1123.2", "MA1638.1", 
  "MA1563.1", "MA1120.1", "MA0624.1", "MA0606.1", "MA1604.1", 
  "MA0773.1", "MA0098.3", "MA0750.2"
)

MAMMARY_TFS <- c(
  "STAT5A", "STAT5B", "LMX1B", "TEAD3", "FOXA1", "FOXO3", "SPIC", "SPIB", 
  "TWIST1", "HAND2", "SOX18", "SOX13", "NFATC1", "NFAT5", "EBF2", "MEF2D", 
  "ETS1", "ZBTB7A"
)

# 1.2 Liver Target Dictionary
LIVER_MOTIFS <- c(
  "MA0114.4", "MA0065.2", "MA0259.1", "MA0482.2", "MA0620.3", 
  "MA0871.2", "MA0002.2", "MA0511.2", "MA0684.2", "MA0475.2", 
  "MA0700.2", "MA0090.3", "MA0507.1", "MA0014.3", "MA0645.1", "MA0652.1"
)

LIVER_TFS <- c(
  "HNF4A", "PPARGC1A", "HIF1A", "GATA4", "MITF", "TFEC", "RUNX1", "RUNX2", 
  "RUNX3", "FLI1", "LHX2", "TEAD1", "POU2F2", "PAX5", "ETV6", "IRF8"
)

# Active switch variables for runtime target routing
current_motifs <- LIVER_MOTIFS  # Swap to MAMMARY_MOTIFS depending on tissue batch
current_tfs    <- LIVER_TFS     # Swap to MAMMARY_TFS depending on tissue batch

# Color palette matrices configurations
col_palette  <- c("#0099CC", "white", "#CC0033")
anno_colors  <- c("#0099CC", "white", "#CC0033") # Custom cohort vector (e.g. col / col1)

# ==========================================
# 2. Standard Base ChromVAR Heatmap
# ==========================================
message("--> Generating base pheatmap for motif chromVAR deviations...")

# Calculate aggregated average profiling values across standard active clusters
aver_chromvar   <- AverageExpression(sc, assays = "chromvar", features = current_motifs, return.seurat = FALSE)
chromvar_matrix <- aver_chromvar$chromvar

# Safe execution block using explicit file connection mapping instead of direct ggsave() wrappers
png("motif_activity_pheatmap.png", width = 5, height = 6, units = "in", res = 300)
pheatmap(
  mat           = chromvar_matrix,
  scale         = "row", 
  cluster_cols  = FALSE,
  cluster_rows  = FALSE,
  show_rownames = TRUE,
  color         = colorRampPalette(col_palette)(100)
)
dev.off()

# ==========================================
# 3. High-Level ClusterGVis Analytics Plots
# ==========================================

# 3.1 Motif Accessibility Profiling Heatmap Track
message("--> Exporting high-resolution TF ChromVAR Activity Track...")
pdf("TF_chromvar_activity_advanced.pdf", width = 5, height = 10)
averageHeatmap(
  object           = sc,               
  markerGene       = current_motifs,  
  group.by         = "ident",          
  assays           = "chromvar",          
  slot             = "data",               
  htCol            = col_palette,  
  htRange          = c(-2, 0, 2),       
  showRowNames     = TRUE,         
  row_names_side   = "right",     
  fontsize         = 10,               
  column_names_rot = 45,       
  cluster_columns  = FALSE,     
  cluster_rows     = FALSE,         
  annoCol          = TRUE,              
  myanCol          = anno_colors
)
dev.off()

# 3.2 Dynamic Intercept Map for Parallel Single-Cell RNA Core Gene Expressions
if (exists("scrna")) {
  message("--> Exporting corresponding single-cell Transcriptional Expression Track...")
  pdf("TF_rna_expression_advanced.pdf", width = 5, height = 10)
  averageHeatmap(
    object           = scrna,               
    markerGene       = current_tfs,  
    group.by         = "ident",          
    assays           = "RNA",          
    slot             = "data",               
    htCol            = col_palette,  
    htRange          = c(-2, 0, 2),       
    showRowNames     = TRUE,         
    row_names_side   = "right",     
    fontsize         = 10,               
    column_names_rot = 45,       
    cluster_columns  = FALSE,     
    cluster_rows     = FALSE,         
    annoCol          = TRUE,              
    myanCol          = anno_colors
  )
  dev.off()
} else {
  warning("Skipping step 3.2: Reference scRNA-seq object ('scrna') is absent from current workspace workspace.")
}

# 3.3 Alternative scRNAtoolVis Export Layout Wrapper
message("--> Compiling secondary validation framework grids via scRNAtoolVis...")
p_rna_tool <- AverageHeatmap(
  object     = sc,
  markerGene = current_motifs
)
ggsave("RNA_activity_scRNAtoolVis.png", plot = p_rna_tool, width = 5, height = 6, units = "in", dpi = 300)

message("--> Complete visualization data panels compiled successfully.")