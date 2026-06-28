#' @title Comprehensive scATAC-seq Downstream Analysis Pipeline
#' @description Implements QC filtering, Harmony batch correction, gene activity matrix estimation,
#' differential accessibility (DA) testing, motif enrichment, genomic foot printing, and vectorized SNP overlap.
#' @author Chen Ao
#' @date 2026-06-28

library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(GenomicRanges)
library(tidyverse)
library(harmony)
library(future)

# ==========================================
# 1. High-Performance Context Initialization
# ==========================================

# Use multisession instead of deprecated multiprocess
plan("multisession", workers = 4)
options(future.globals.maxSize = 50 * 1024^3) # 50 GB RAM allocation allocation

# ==========================================
# 2. Quality Control & Filtering Module
# ==========================================

#' Execute Standard Signac QC Metrics Derivation
#' @param seurat_obj Pre-merged scATAC Seurat object.
#' @param annotations GRanges reference object containing gene models (e.g., UCD1.3).
#' @return A Seurat object appended with QC score columns.
compute_atac_qc_metrics <- function(seurat_obj, annotations) {
  message("--> Injecting genome annotations and computing nucleosome vs TSS scoring...")
  
  Annotation(seurat_obj) <- annotations
  
  # Compute signal metrics profiles
  seurat_obj <- NucleosomeSignal(object = seurat_obj)
  seurat_obj <- TSSEnrichment(object = seurat_obj, fast = FALSE)
  
  # Derive additional summary proportions
  seurat_obj$pct_reads_in_peaks <- (seurat_obj$peak_region_fragments / seurat_obj$passed_filters) * 100
  seurat_obj$blacklist_ratio    <- seurat_obj$blacklist_region_fragments / seurat_obj$peak_region_fragments
  
  return(seurat_obj)
}

# ==========================================
# 3. Dimensional Reduction & Integration
# ==========================================

#' Run TF-IDF, SVD, Harmony Integration, and Neighbor Clustering
#' @param seurat_obj Subsetted scATAC Seurat object.
#' @param batch_var Latent metadata column name used for batch integration (e.g., "orig.ident").
#' @return A fully integrated, embedded, and clustered Seurat object.
run_atac_integration <- function(seurat_obj, batch_var = "orig.ident") {
  message("--> Executing TF-IDF normalization and SVD reduction...")
  DefaultAssay(seurat_obj) <- "peaks"
  
  seurat_obj <- RunTFIDF(seurat_obj)
  seurat_obj <- FindTopFeatures(seurat_obj, min.cutoff = "q0")
  seurat_obj <- RunSVD(seurat_obj, reduction.key = "LSI_", reduction.name = "lsi")
  
  message("--> Running Harmony batch effect correction...")
  seurat_obj <- RunHarmony(
    object          = seurat_obj,
    group.by.vars   = batch_var,
    reduction       = "lsi",
    assay.use       = "peaks",
    project.dim     = FALSE,
    reduction.save  = "harmony"
  )
  
  message("--> Computing UMAP embeddings and clustering structures...")
  seurat_obj <- RunUMAP(seurat_obj, dims = 2:30, reduction = "harmony")
  seurat_obj <- FindNeighbors(seurat_obj, reduction = "harmony", dims = 2:30)
  seurat_obj <- FindClusters(seurat_obj, verbose = FALSE, algorithm = 3, resolution = 0.2)
  
  return(seurat_obj)
}

# ==========================================
# 4. Multimodal Gene Activity Estimation
# ==========================================

#' Construct Quantified Gene Activity Matrix Proxy Layer
#' @param seurat_obj Processed Seurat object.
#' @return Seurat object containing an independent log-normalized "RNA" assay slot.
estimate_gene_activities <- function(seurat_obj) {
  message("--> Quantifying gene activities across peak profiles...")
  gene_activities <- GeneActivity(seurat_obj)
  
  seurat_obj[["RNA"]] <- CreateAssayObject(counts = gene_activities)
  seurat_obj <- NormalizeData(
    object               = seurat_obj,
    assay                = "RNA",
    normalization.method = "LogNormalize",
    scale.factor         = median(seurat_obj$nCount_RNA)
  )
  return(seurat_obj)
}

# ==========================================
# 5. Advanced Differential & Motif Analysis
# ==========================================

#' Run Differential Accessibility (DA) Testing and Downstream Motif Enrichment
#' @param seurat_obj Fully annotated Seurat workspace container.
#' @param cluster_column Column name storing target cell identity clusters.
#' @param ident_1 Group name 1 for pairwise comparison testing.
#' @param ident_2 Group name 2 for pairwise comparison testing.
#' @return A list containing DA peaks and significantly enriched motif matrices.
run_da_and_motif_enrichment <- function(seurat_obj, cluster_column, ident_1, ident_2) {
  message(paste("--> Running DA testing via LR for:", ident_1, "vs", ident_2))
  Idents(seurat_obj)       <- cluster_column
  DefaultAssay(seurat_obj) <- "peaks"
  
  da_peaks <- FindMarkers(
    object      = seurat_obj,
    ident.1     = ident_1,
    ident.2     = ident_2,
    test.use    = "LR",
    min.pct     = 0.05,
    latent.vars = "nCount_peaks"
  )
  
  # Select stringently increased target elements
  target_features <- rownames(da_peaks[da_peaks$p_val < 0.05 & da_peaks$avg_log2FC > 0.5, ])
  
  message("--> Executing motif enrichment analysis on top hyper-accessible peaks...")
  enriched_motifs <- FindMotifs(object = seurat_obj, features = target_features)
  
  return(list(da_peaks = da_peaks, enriched_motifs = enriched_motifs))
}

# ==========================================
# 6. High-Throughput Vectorized Overlap Modules
# ==========================================

#' Ultra-Fast Vectorized Chromosome Intersections for Genomic Variant Coordinates
#' @param snp_df Data frame containing raw variant records (columns: chr, position).
#' @param seurat_obj Target Seurat reference containing target peaks.
#' @return Filtered subset of data frame entries residing within defined chromatin peak boundaries.
intersect_snps_with_peaks <- function(snp_df, seurat_obj) {
  message("--> Performing high-performance vectorized overlaps via GenomicRanges...")
  
  # Parse standard peak structures
  peak_regions_df <- data.frame(region = rownames(seurat_obj)) %>%
    tidyr::separate(region, into = c("chr", "start", "end"), sep = "-", convert = TRUE)
  
  peaks_gr <- makeGRangesFromDataFrame(peak_regions_df)
  
  # Construct GRanges matrix structure for target SNPs directly to bypass slow row loops
  snps_gr <- GRanges(
    seqnames = snp_df$chr,
    ranges   = IRanges(start = snp_df$position, end = snp_df$position)
  )
  
  # Execute lightning-fast geometric vector intersection checks
  overlaps <- findOverlaps(query = snps_gr, subject = peaks_gr)
  
  filtered_snp_df <- snp_df[unique(queryHits(overlaps)), ] %>%
    mutate(snp_pos = paste(chr, position, sep = "-"))
  
  return(filtered_snp_df)
}