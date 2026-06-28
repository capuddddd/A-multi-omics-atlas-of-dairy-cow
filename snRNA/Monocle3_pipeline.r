#' @title Modular Single-Cell Trajectory & Pseudotime Analysis Pipeline (Monocle3)
#' @description Comprehensive workflow wrapping Monocle3, ClusterGVis, and custom GAM/Spline fitting
#' for scRNA-seq and scATAC-seq multimodal datasets.
#' @author Chen Ao
#' @date 2026-06-28

library(Seurat)
library(monocle3)
library(tidyverse)
library(patchwork)
library(SeuratWrappers)
library(ClusterGVis)
library(reshape2)
library(ComplexHeatmap)

# ==========================================
# 1. Monocle3 Core Pipeline Functions
# ==========================================

#' Initialize and Align CellDataSet from Seurat Object
#' @param seurat_obj Seurat object context.
#' @param assay Active target assay. Default 'RNA'.
#' @param alignment_field Covariate formula covariate string (e.g. 'sample').
#' @return Preprocessed and aligned cell_data_set object.
initialize_monocle_cds <- function(seurat_obj, assay = "RNA", alignment_field = "sample") {
  message("--> Building CellDataSet from matrix counts...")
  counts_mat <- GetAssayData(seurat_obj, assay = assay, layer = "counts")
  
  # Remove low-frequency filtering artifacts
  counts_mat <- counts_mat[rowSums(counts_mat > 0) >= 3, ]
  
  gene_anno <- data.frame(gene_short_name = rownames(counts_mat), row.names = rownames(counts_mat))
  
  cds <- new_cell_data_set(
    expression_data = counts_mat,
    cell_metadata   = seurat_obj@meta.data,
    gene_metadata   = gene_anno
  )
  
  # Structural alignment and batch correction
  cds <- preprocess_cds(cds, num_dim = 50, residual_model_formula_str = paste0("~", alignment_field))
  cds <- align_cds(cds, alignment_group = alignment_field, useNames = FALSE)
  cds <- reduce_dimension(cds, preprocess_method = "PCA", reduction_method = "UMAP")
  
  return(cds)
}

#' Transfer Native Integrated Embeddings into Monocle Slot Matrix
#' @param cds Monocle cell_data_set object.
#' @param seurat_obj Target Seurat reference container.
#' @param reduction Character string name of reduction layer (e.g. 'umap').
#' @return Modifed cell_data_set with direct coordinate spatial projections.
sync_embeddings_to_cds <- function(cds, seurat_obj, reduction = "umap") {
  message("--> Syncing native integrated coordinate matrices to CDS slot...")
  seurat_embeds <- Embeddings(seurat_obj, reduction = reduction)
  seurat_embeds <- seurat_embeds[rownames(cds@int_colData$reducedDims$UMAP), ]
  cds@int_colData$reducedDims$UMAP <- seurat_embeds
  return(cds)
}

# ==========================================
# 2. Advanced Mathematical Pseudotime Routing
# ==========================================

#' Extract Target Root Principal Nodes Dynamically Based on Max Cluster Identity
#' @param cds Preprocessed cell_data_set.
#' @param celltype_col Annotation field used to evaluate starting cluster.
#' @param target_identity Cluster subset name serving as root origin.
get_earliest_principal_node <- function(cds, celltype_col = "celltype", target_identity = NULL) {
  if (is.null(target_identity)) stop("Must specify target identity for root mapping calculation.")
  
  cell_indices <- which(colData(cds)[, celltype_col] == target_identity)
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  
  root_node <- igraph::V(principal_graph(cds)[["UMAP"]])$name[
    as.numeric(names(which.max(table(closest_vertex[cell_indices, ]))))
  ]
  return(root_node)
}

# ==========================================
# 3. Dynamic Visualizations & Kinetic Modeling
# ==========================================

#' Batch Plot and Fit Gene Kinetics Across Pseudotime Curves via GAM Modeling
#' @param seurat_obj Seurat matrix container.
#' @param gene_list Explicit vector array of features to plot.
#' @param output_dir Character string target output folder destination.
plot_gam_pseudotime_kinetics <- function(seurat_obj, gene_list, output_dir = "./gam_plots") {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  if (!"pseudotime" %in% colnames(seurat_obj@meta.data)) {
    stop("Pseudotime vector is missing from metadata slot container hierarchy.")
  }
  
  message("--> Starting batch calculation of Generalized Additive Models (GAM)...")
  for (gene in gene_list) {
    if (!gene %in% rownames(seurat_obj)) {
      warning(paste("Feature skipped. Missing from object:", gene))
      next
    }
    
    expr_vals <- FetchData(seurat_obj, vars = gene)[, 1]
    plot_df   <- data.frame(Pseudotime = seurat_obj$pseudotime, Expression = expr_vals)
    
    p_gam <- ggplot(plot_df, aes(x = Pseudotime, y = Expression)) +
      geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = TRUE, color = "#2B5C8F", fill = "#D4E2F0") +
      theme_classic() +
      labs(title = gene, x = NULL, y = NULL) +
      theme(
        plot.title   = element_text(hjust = 0.5, size = 16, face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.8),
        axis.text.x  = element_text(size = 11, color = "black"),
        axis.text.y  = element_text(size = 11, color = "black")
      )
    
    ggsave(file.path(output_dir, paste0(gene, "_pseudotime_gam.pdf")), plot = p_gam, width = 4, height = 3)
  }
}

#' Standardized Smooth Spline Matrix Generation for High-Throughput Heatmaps
#' @param seurat_obj Seurat data structure.
#' @param target_genes Vector list containing target trajectory DEGs.
#' @param output_path File path destination for pdf object generation.
plot_spline_trajectory_heatmap <- function(seurat_obj, target_genes, output_path = "spline_heatmap.pdf") {
  message("--> Building smooth spline normalization matrices...")
  
  # Intersect arrays to secure matching data indexes
  valid_genes <- intersect(target_genes, rownames(seurat_obj))
  sorted_cells <- colnames(seurat_obj)[order(seurat_obj$pseudotime)]
  
  raw_mat <- as.matrix(GetAssayData(seurat_obj, assay = "RNA", layer = "data")[valid_genes, sorted_cells])
  
  # Apply 3-degree of freedom smoothing spline fits row-by-row
  smoothed_mat <- t(apply(raw_mat, 1, function(row_vals) {
    smooth.spline(row_vals, df = 3)$y
  }))
  
  # High-performance Z-Score scaling centering normalization rows
  scaled_mat <- t(apply(smoothed_mat, 1, function(row_vals) {
    (row_vals - mean(row_vals)) / sd(row_vals)
  }))
  rownames(scaled_mat) <- valid_genes
  
  # Render final matrix output utilizing ComplexHeatmap engine blocks
  pdf(output_path, width = 8, height = 10)
  ht <- Heatmap(
    scaled_mat,
    name              = "Z-Score Expression",
    cluster_columns   = FALSE,
    show_column_names = FALSE,
    cluster_rows      = TRUE,
    show_row_dend     = FALSE,
    col               = colorRampPalette(c("#5bacdb", "white", "#ea7ea3"))(100),
    row_names_gp      = gpar(fontsize = 8)
  )
  draw(ht)
  dev.off()
}