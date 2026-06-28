#' @title Advanced TF Analysis Pipeline: ChromVAR, Pseudotime Cor, and Footprinting
#' @description Computes ChromVAR deviations, identifies pseudotime-associated motifs,
#' resolves Bos taurus chromosome nomenclature mismatches, and plots genomic footprints.
#' @author Chen Ao
#' @date 2026-06-28

library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(tidyverse)
library(BSgenome.Btaurus.UCSC.bosTau9)
library(TFBSTools)

# ==========================================
# 1. Motif Matrix & ChromVAR Initialization
# ==========================================

DefaultAssay(sc) <- "peaks"

# Load local JASPAR2024 PFM repository mapping
message("--> Loading custom JASPAR2024 position frequency matrices...")
pfm <- readJASPARMatrix(
  file        = "./jaspar2024.txt", 
  matrixClass = "PFM"
)

# Append Motif assignments to the Seurat Object
sc <- AddMotifs(
  object = sc,
  genome = BSgenome.Btaurus.UCSC.bosTau9,
  pfm    = pfm
)

# Compute per-cell motif accessibility deviations via chromVAR
message("--> Running ChromVAR background model and deviation scoring...")
sc <- RunChromVAR(
  object = sc,
  genome = BSgenome.Btaurus.UCSC.bosTau9
)

# ==========================================
# 2. Pseudotime & ChromVAR Dynamic Modeling
# ==========================================

DefaultAssay(sc) <- "chromvar"

# Extract latent pseudotime dynamic vector
cell_pseudotime <- sc@meta.data$pseudotime

message("--> Calculating Spearman correlations between TF scores and trajectory...")
chromvar_matrix <- GetAssayData(sc, assay = "chromvar", layer = "data")

# Vectorized correlation analysis to identify differentiation-driven TFs
motif_correlations <- apply(chromvar_matrix, 1, function(motif_row) {
  if (all(motif_row == 0) || sd(motif_row) == 0) return(0) # Prevent variance zero errors
  cor(motif_row, cell_pseudotime, method = "spearman", use = "complete.obs")
})

motif_correlations_df <- data.frame(
  motif_id    = names(motif_correlations),
  correlation = as.numeric(motif_correlations),
  stringsAsFactors = FALSE
)

# Extract dynamic driving motifs (Abs value > 0.3)
significant_motifs <- motif_correlations_df %>% 
  filter(abs(correlation) > 0.3) %>%
  arrange(desc(abs(correlation)))

write.csv(significant_motifs, "pseudotime_correlated_motifs.csv", row.names = FALSE)

# ==========================================
# 3. Structural Refactoring of Chromosome Styles
# ==========================================
message("--> Harmonizing genomic sequence levels between Object and BSgenome reference...")

# NOTE: BSgenome is a read-only locked object. We must modify the Seurat object's 
# seqlevelsStyle instead of attempting to overwrite BSgenome.Btaurus.UCSC.bosTau9 slots.

DefaultAssay(sc) <- "peaks"

# Force Seurat object coordinates to drop "chr" prefixes and use NCBI style (1:29, X, MT)
seqlevelsStyle(sc) <- "NCBI"

# If minor scaffolds or unplaced contigs (e.g., Un/Random) remain in your matrix, 
# stringently subset to standard Bos taurus karyotype elements to prevent Footprint failures
standard_bovine_chroms <- c(as.character(1:29), "X", "MT")
current_sc_levels      <- seqlevels(sc)
valid_sc_levels        <- current_sc_levels[current_sc_levels %in% standard_bovine_chroms]

sc <- keepSeqlevels(sc, value = valid_sc_levels, pruning.mode = "coarse")

# ==========================================
# 4. Transcription Factor Footprinting & Plotting
# ==========================================
message("--> Initiating high-resolution TF genomic footprinting calculations...")

# Run footprinting analysis for prioritized lineages (Runx1 and others)
# Ensure the motif names exactly match the identifiers within names(pfm)
target_footprints <- c("Runx1", "SOX4", "Stat5a::Stat5b", "NR3C1")

# Verify matching motifs exist in the database before processing
available_motifs <- names(pfm)
valid_footprints <- target_footprints[target_footprints %in% available_motifs || 
                                      HasMotif(sc, motif = target_footprints)]

sc <- Footprint(
  object     = sc,
  motif.name = valid_footprints,
  genome     = BSgenome.Btaurus.UCSC.bosTau9
)

# 4.2 Diagnostics and Visualizations Plotting
message("--> Generating diagnostic footprint tracks...")

p_footprint <- PlotFootprint(sc, features = "Runx1", group.by = "group") +
  theme_classic() +
  labs(title = "Runx1 Genomic Footprint Track", x = "Distance from Motif Center (bp)")

ggsave(
  filename = "Runx1_Footprint_Track.pdf", 
  plot     = p_footprint, 
  width    = 6, 
  height   = 5, 
  dpi      = 300
)

message("--> Advanced TF trajectory and footprint analysis successfully completed.")