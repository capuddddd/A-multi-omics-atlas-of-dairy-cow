#' @title Single-Cell Quality Control and Downstream Analysis Pipeline
#' @description Prepares Seurat objects by removing ambient RNA contamination, merging 
#' metadata, filtering low-quality cells, and running downstream dimension reduction.
#' @author Chen Ao
#' @date 2026-06-28

library(dplyr)
library(Seurat)
library(patchwork)
library(DoubletFinder)
library(cowplot)
library(clustree)
library(presto)
library(RColorBrewer)
library(ggpubr)
library(ggplot2)
library(harmony)
library(egg)
library(decontX)

# ==========================================
# 1. Helper Functions
# ==========================================

#' Merge multi-lane or split metadata columns (e.g., from DoubletFinder) into a single vector
merge_metadata_columns <- function(meta_data, pattern) {
  target_cols <- grep(pattern, colnames(meta_data), value = TRUE)
  combined_vector <- NULL
  
  for (col in target_cols) {
    non_na_data <- meta_data[!is.na(meta_data[[col]]), col]
    if (is.null(combined_vector)) {
      combined_vector <- non_na_data
    } else {
      combined_vector <- c(combined_vector, non_na_data)
    }
  }
  return(combined_vector)
}

library(DoubletFinder)
library(Seurat)
library(dplyr)

# ==========================================
# 0. Pipeline Engine Function
# ==========================================

#' Run Complete DoubletFinder Workflow
#' @param rds_path String path to input Seurat object.
#' @param output_path String path to save processed Seurat object.
#' @param pcs Integer vector specifying PCs to use for embedding. Default 1:21.
#' @param base_doublet_rate Numeric percentage representing expected baseline doublet rate. Default 0.075 (7.5%).
run_doublet_finder_pipeline <- function(rds_path, output_path, pcs = 1:21, base_doublet_rate = 0.075) {
  
  if (!file.exists(rds_path)) {
    stop(paste("Input file does not exist:", rds_path))
  }
  
  message("--> Loading Seurat workspace data...")
  sc <- readRDS(rds_path)
  
  # 0.1 Parameter Sweep (pK Identification without ground-truth)
  message("--> Running parameter sweep distributions...")
  sweep_res_list <- paramSweep(sc, PCs = pcs, sct = FALSE)
  sweep_stats    <- summarizeSweep(sweep_res_list, GT = FALSE)
  bcmvn          <- find.pK(sweep_stats)
  
  # Dynamically extract pK value at maximum mean BC metric
  optimal_pK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  message(paste("--> Dynamic Optimal pK identified as:", optimal_pK))
  
  # 0.2 Homotypic Doublet Proportion Calculations
  message("--> Modeling homotypic doublet proportions...")
  homotypic_prop <- modelHomotypic(sc@meta.data$seurat_clusters)
  
  nExp_poi     <- round(base_doublet_rate * nrow(sc@meta.data))
  nExp_poi_adj <- round(nExp_poi * (1 - homotypic_prop))
  
  # 0.3 First Pass: Standard Detection Unadjusted (Generates baseline pANN scores)
  message("--> Executing DoubletFinder Pass 1 (Unadjusted)...")
  sc <- doubletFinder(sc, PCs = pcs, pN = 0.25, pK = optimal_pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  
  # Extract names dynamically to bypass hardcoded column index dependencies
  pANN_col_name <- grep(paste0("^pANN_0.25_", optimal_pK), colnames(sc@meta.data), value = TRUE)
  DF_raw_name   <- grep(paste0("^DF.classifications_0.25_", optimal_pK, "_", nExp_poi, "$"), colnames(sc@meta.data), value = TRUE)
  
  # 0.4 Second Pass: Refinement via Homotypic Adjustment Strategy
  message("--> Executing DoubletFinder Pass 2 (Homotypic Adjusted)...")
  sc <- doubletFinder(sc, PCs = pcs, pN = 0.25, pK = optimal_pK, nExp = nExp_poi_adj, reuse.pANN = pANN_col_name, sct = FALSE)
  
  DF_adj_name <- grep(paste0("^DF.classifications_0.25_", optimal_pK, "_", nExp_poi_adj, "$"), colnames(sc@meta.data), value = TRUE)
  
  # 0.5 Tiered Confidence Scoring Evaluation
  message("--> Consolidating tiered doublet classifications...")
  sc$DF_hi.lo <- sc@meta.data[[DF_adj_name]] # Start with adjusted predictions as base
  
  # If cell is classified as a Doublet under loose parameters but Singlet under stricter ones -> Low Confidence
  sc$DF_hi.lo[sc@meta.data[[DF_raw_name]] == "Doublet" & sc@meta.data[[DF_adj_name]] == "Singlet"] <- "Doublet-Low Confidence"
  sc$DF_hi.lo[sc$DF_hi.lo == "Doublet"] <- "Doublet-High Confidence"
  
  # Output descriptive execution table logs
  message("Final script doublet distribution classifications summary:")
  print(table(sc$DF_hi.lo))
  
  # 1.6 File Serialization Persistence Output
  message(paste("--> Saving processed object workspace ->", output_path))
  saveRDS(sc, output_path)
  message("DoubletFinder pipeline completed successfully.")
}

# ==========================================
# 2. Metadata Consolidation & Ambient RNA QC
# ==========================================

# NOTE: Execute this workflow block with your active workspace 'sc' Seurat object.
# consolidate_and_filter_contamination <- function(sc, output_dir = ".") {

# 2.1 Merge split DoubletFinder metrics into unified metadata columns
meta_data <- sc@meta.data
pANN_combined <- merge_metadata_columns(meta_data, "^pANN")
DF_combined <- merge_metadata_columns(meta_data, "^DF.classifications")

sc <- AddMetaData(sc, metadata = pANN_combined, col.name = "pANN_combined")
sc <- AddMetaData(sc, metadata = DF_combined, col.name = "DF_combined")

# Clean up split intermediate columns to keep metadata tidy
cols_to_remove <- grep("^pANN_|^DF.classifications", colnames(sc@meta.data), value = TRUE)
sc@meta.data <- sc@meta.data %>% select(-all_of(cols_to_remove))

# 2.2 Run decontX to remove ambient RNA background matrix contamination
counts <- GetAssayData(sc, assay = "RNA", layer = "counts")
decontX_results <- decontX(counts) 
sc$Contamination <- decontX_results$contamination

# Update raw counts with decontX decontaminated expression matrix
sc <- SetAssayData(sc, assay = "RNA", layer = "counts", new.data = decontX_results$decontXcounts)

# 2.3 Strict QC Filtering Thresholds
sc <- subset(sc, subset = Contamination < 0.3 & pANN_combined < 0.4 & DF_combined != "Doublet-High Confidence")

# 2.4 Diagnostic Plot Generation
p1 <- FeaturePlot(sc, features = "pANN_combined", raster = FALSE) +
  scale_color_viridis_c() +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  xlab("scVI_UMAP_1") +
  ylab("scVI_UMAP_2")

ggsave("pANN_feature_plot.pdf", egg::set_panel_size(p1, width = unit(5, "in"), height = unit(5, "in")), 
       width = 9, height = 9, units = 'in', dpi = 300)

# ==========================================
# 3. Main Seurat Analysis Pipeline Function
# ==========================================

#' Standardized downstream workflow for subpopulation refinement
run_seurat_analysis <- function(rds_path, output_dir = ".", sample_col = "sample", dpi = 300) {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  sc <- readRDS(rds_path)
  
  # 3.1 Mitochondrial and Ribosomal QC Assessment (Bos taurus Gene Set)
  mt_genes <- c("ATP6", "ATP8", "COX1", "COX2", "COX3", "CYTB", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6")
  mt_present <- intersect(mt_genes, rownames(sc))
  
  raw_counts <- GetAssayData(sc, assay = "RNA", layer = "counts")
  
  sc[["percent.mt"]] <- colSums(raw_counts[mt_present, , drop = FALSE]) / sc[[paste0("nCount_", "RNA")]] * 100
  
  rb_genes <- rownames(sc)[grep("^RP[SL]", rownames(sc))]
  sc[["percent.ribo"]] <- colSums(raw_counts[rb_genes, , drop = FALSE]) / sc[[paste0("nCount_", "RNA")]] * 100
  
  # Filter out structural mitochondrial and ribosomal expression from variable workspace
  sc <- subset(sc, features = setdiff(rownames(sc), c(mt_present, rb_genes)))
  
  # Export pre-filter distribution checks
  pdf(file.path(output_dir, "quality_vlnplot_clusters.pdf"), height = 10, width = 20)
  print(VlnPlot(sc, features = c("nFeature_RNA", "nCount_RNA", "percent.ribo"), pt.size = 0, group.by = "seurat_clusters", ncol = 3, raster = FALSE))
  dev.off()
  
  # Subset cells based on threshold settings
  sc <- subset(sc, subset = nFeature_RNA > 300 & nFeature_RNA < 20000 & percent.ribo < 20 & percent.mt < 20)
  
  # 3.2 Normalization and Scaling (with Cell Cycle Regression)
  sc <- NormalizeData(sc, normalization.method = "LogNormalize", scale.factor = 10000)
  sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 2000)
  
  # Verify presence of Cell Cycle Scores prior to regression execution
  regression_vars <- intersect(c("S.Score", "G2M.Score"), colnames(sc@meta.data))
  if (length(regression_vars) > 0) {
    sc <- ScaleData(sc, features = rownames(sc), vars.to.regress = regression_vars)
  } else {
    sc <- ScaleData(sc, features = rownames(sc))
  }
  
  # 3.3 Linear Dimension Reduction and Batch Integration via Harmony
  sc <- RunPCA(sc, features = VariableFeatures(object = sc))
  sc <- RunHarmony(sc, reduction = "pca", group.by.vars = sample_col, reduction.save = "harmony")
  
  # Determine optimal Harmony PCs based on cumulative variance standard deviations
  pct <- sc[["harmony"]]@stdev / sum(sc[["harmony"]]@stdev) * 100
  cumu <- cumsum(pct)
  co1 <- which(cumu > 90 & pct < 5)[1]
  co2 <- sort(which((pct[1:(length(pct) - 1)] - pct[2:length(pct)]) > 0.1), decreasing = TRUE)[1] + 1
  pcs <- min(co1, co2)
  bestpc <- 1:pcs
  
  # 3.4 Non-linear Dimension Reduction and Clustering Graph Optimization
  sc <- RunUMAP(sc, reduction = "harmony", dims = bestpc, reduction.name = "umap", min.dist = 0.4)
  sc <- FindNeighbors(sc, reduction = "harmony", dims = bestpc)
  
  # Clustree evaluation across variable resolutions
  sc <- FindClusters(sc, resolution = seq(0.1, 0.5, by = 0.1))
  pdf(file.path(output_dir, "Harmony_clustree.pdf"))
  print(clustree(sc))
  dev.off()
  
  # Finalize targeted structural clustering assignment
  sc <- FindClusters(sc, resolution = 0.1)
  
  # 3.5 Generate and Save Plots
  p_orig <- DimPlot(sc, reduction = "umap", group.by = "orig.ident", raster = FALSE)
  ggsave(file.path(output_dir, "umap_orig_ident.pdf"), egg::set_panel_size(p_orig, width = unit(5, "in"), height = unit(5, "in")), 
         width = 9, height = 9, units = 'in', dpi = dpi)
  
  p_clust <- DimPlot(sc, reduction = "umap", group.by = "seurat_clusters", label = TRUE, raster = FALSE)
  ggsave(file.path(output_dir, "umap_seurat_clusters.pdf"), egg::set_panel_size(p_clust, width = unit(5, "in"), height = unit(5, "in")),
         width = 9, height = 9, units = "in", dpi = dpi)
  
  # Save clean structured downstream reference data
  saveRDS(sc, file.path(output_dir, "merge_processed.rds"))
}

library(Seurat)
library(clusterProfiler)
library(org.Bt.eg.db)
library(tidyverse)

# ==========================================
# 4. Functional Enrichment Pipeline Engine
# ==========================================

#' Run Comparative Enrichment and Generate Structured Visualizations
#' @param seurat_obj Seurat workspace object.
#' @param group_col Metadata column string used to split sample cohorts. Default is 'group'.
#' @param output_dir Output string path for saving metrics and graphics.
run_comparative_enrichment <- function(seurat_obj, group_col = "group", output_dir = ".") {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # 4.1 Calculate Differential Marker Genes Across Defined Cluster Identities
  message("--> Calculating all cluster-specific markers...")
  Idents(seurat_obj) <- group_col
  obj_markers <- FindAllMarkers(seurat_obj, only.pos = TRUE)
  saveRDS(obj_markers, file.path(output_dir, "obj.markers.rds"))
  
  if (nrow(obj_markers) == 0) {
    stop("No significant markers identified. Halting downstream enrichment processing steps.")
  }
  
  # 4.2 Multi-Omic Identifier Key Transformations (Bos taurus annotation)
  message("--> Converting Bovine SYMBOL IDs to ENTREZID mappings...")
  symbol_to_entrez <- mapIds(
    org.Bt.eg.db, 
    keys    = obj_markers$gene, 
    keytype = "SYMBOL", 
    column  = "ENTREZID"
  )
  
  bitr_df <- bitr(
    obj_markers$gene, 
    fromType = "SYMBOL", 
    toType   = "ENTREZID", 
    OrgDb    = "org.Bt.eg.db"
  )
  
  # Structural validation merge
  merged_markers <- merge(obj_markers, bitr_df, by.x = "gene", by.y = "SYMBOL")
  
  # Split lists by identifier groups
  gc_sample <- split(merged_markers$ENTREZID, merged_markers$cluster)
  
  # 4.3 ClusterProfiler Comparative Pathway Analysis Execution
  message("--> Processing comparative KEGG functional arrays...")
  # NOTE: To toggle Gene Ontology biological processes, uncomment the block line sequence below:
  # xx <- compareCluster(gc_sample, fun = "enrichGO", OrgDb = "org.Bt.eg.db", ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
  
  comp_kegg <- compareCluster(
    gc_sample, 
    fun          = "enrichKEGG", 
    organism     = "bta", 
    pvalueCutoff = 1, 
    qvalueCutoff = 1
  )
  
  res_matrix <- as.data.frame(comp_kegg@compareClusterResult)
  if (nrow(res_matrix) == 0) {
    stop("Comparative output returned empty dataframes. Review input threshold selections.")
  }
  
  # 4.4 Translate Clustered Enrichment ENTREZID Lists Back Into Clean SYMBOL Annotations
  message("--> Remapping nested ENTREZID lists back to explicit SYMBOL markers...")
  res_matrix$geneID <- sapply(res_matrix$geneID, function(gene_string) {
    entrez_vector <- unlist(strsplit(as.character(gene_string), split = "/"))
    mapped_symbols <- names(symbol_to_entrez[symbol_to_entrez %in% entrez_vector])
    clean_symbols <- paste(unique(mapped_symbols), collapse = "/")
    return(clean_symbols)
  })
  
  write.csv(res_matrix, file.path(output_dir, "comparative_kegg_results.csv"), row.names = FALSE)
  
  # ==========================================
  # Graphics Generation Step
  # ==========================================
  
  message("--> Generating publication-ready ggplot outputs...")
  
  # Filter top entries per experimental condition (e.g. Dry vs Peak Lactation Periods)
  target_clusters <- intersect(c("Dry", "Lac"), unique(res_matrix$Cluster))
  
  plot_data <- res_matrix %>%
    filter(Cluster %in% target_clusters) %>%
    group_by(Cluster) %>%
    top_n(n = 5, wt = -pvalue) %>%
    ungroup() %>%
    arrange(Cluster, desc(pvalue))
  
  # Freeze row assignments to factor ordering sequences to avoid alphabetical re-sorting
  plot_data$Description <- factor(plot_data$Description, levels = unique(plot_data$Description))
  
  # Establish visualization theme profiles
  custom_theme <- theme_classic() +
    theme(
      axis.title          = element_text(size = 13),
      axis.text           = element_text(size = 11),
      axis.text.y         = element_blank(), 
      axis.ticks.length.y = unit(0, "cm"),
      plot.title          = element_text(size = 13, hjust = 0.5, face = "bold"),
      legend.position     = "none",
      plot.margin         = margin(t = 6, r = 10, l = 6, b = 6)
    )
  
  # Auto-calculate dynamic vector lengths for alternating text parameters
  color_assignments <- ifelse(plot_data$Cluster == "Dry", "#6bb9d2", "#d55640")
  
  p_bar <- ggplot(data = plot_data, aes(x = -log10(pvalue), y = Description, fill = Cluster)) +
    geom_bar(stat = "identity", width = 0.6, alpha = 0.8) +
    scale_fill_manual(values = c("Dry" = "#6bb9d2", "Lac" = "#d55640")) +
    scale_x_continuous(expand = c(0, 0)) + 
    labs(
      x     = "-Log10(pvalue)", 
      y     = "Enriched Dynamic Conditions (Lac / Dry)", 
      title = "KEGG Pathway Enrichment Matrix"
    ) +
    # Layer 1: Pathway Description text labels
    geom_text(
      aes(x = 0.05, label = Description), 
      hjust = 0, 
      vjust = -0.5, 
      size  = 3.5, 
      color = "black", 
      fontface = "bold"
    ) +
    # Layer 2: Mapped Enriched Gene ID strings
    geom_text(
      aes(x = 0.05, label = geneID), 
      hjust = 0, 
      vjust = 1.5, 
      size  = 2.8, 
      color = color_assignments
    ) +
    custom_theme

  ggsave(
    filename = file.path(output_dir, "markers_kegg_enrichment.pdf"), 
    plot     = p_bar, 
    width    = 8, 
    height   = 6, 
    dpi      = 300
  )
  
  message("Pipeline execution completed successfully.")
}