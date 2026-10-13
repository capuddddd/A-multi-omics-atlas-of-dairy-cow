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
# Sequence names must also match the underlying indexed fragments; do not rename
# them here to hide an incompatible genome or assembly.
unknown_sequences <- setdiff(seqlevels(sc), seqlevels(BSgenome.Btaurus.UCSC.bosTau9))
if (length(unknown_sequences) > 0) {
  stop("Sequence names do not match bosTau9; reconcile genome, peaks and fragments: ",
       paste(unknown_sequences, collapse = ", "))
}

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
# 3. Select the assay used for footprinting
# ==========================================
DefaultAssay(sc) <- "peaks"

# ==========================================
# 4. Transcription Factor Footprinting & Plotting
# ==========================================
message("--> Initiating high-resolution TF genomic footprinting calculations...")

# Run footprinting analysis for prioritized lineages (Runx1 and others)
# Ensure the motif names exactly match the identifiers within names(pfm)
target_footprints <- c("Runx1", "SOX4", "Stat5a::Stat5b", "NR3C1")

# Footprint expects motif names stored in Signac, not PFM identifier keys.
available_motifs <- unname(unlist(GetMotifData(sc, assay = "peaks", slot = "motif.names")))
valid_footprints <- intersect(target_footprints, available_motifs)
missing_footprints <- setdiff(target_footprints, available_motifs)
if (length(missing_footprints) > 0) {
  warning("Requested motif names were not found: ", paste(missing_footprints, collapse = ", "))
}
if (length(valid_footprints) == 0) stop("None of the requested footprint motifs are available.")

sc <- Footprint(
  object     = sc,
  motif.name = valid_footprints,
  genome     = BSgenome.Btaurus.UCSC.bosTau9
)

# 4.2 Diagnostics and Visualizations Plotting
message("--> Generating diagnostic footprint tracks...")

p_footprint <- PlotFootprint(sc, features = valid_footprints, group.by = "group") +
  theme_classic() +
  labs(title = "Selected Motif Footprints", x = "Distance from Motif Center (bp)")

ggsave(
  filename = "selected_motif_footprints.pdf",
  plot     = p_footprint,
  width    = 6,
  height   = 5,
  dpi      = 300
)

message("--> Advanced TF trajectory and footprint analysis successfully completed.")
