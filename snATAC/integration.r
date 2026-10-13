#' @title Cross-Modal Label Transfer & Co-Embedding Pipeline (scRNA -> scATAC)
#' @description Performs single-cell multi-omic integration, label transfer from
#' scRNA-seq reference to scATAC-seq query, RNA data imputation, and integrated co-embedding.
#' @author Chen Ao
#' @date 2026-06-28

library(Seurat)
library(tidyverse)
library(cowplot)

# ==========================================
# 1. Core Label Transfer Module
# ==========================================

#' Run Cross-Modal Anchor Finding and Label Transfer
#' @param ref_obj Seurat reference object (scRNA-seq)
#' @param query_obj Seurat query object (scATAC-seq gene activity matrix)
#' @param ref_assay Assay name for reference object. Default is "RNA".
#' @param query_assay Assay name for query object. Default is "RNA".
#' @param weight_reduction Computational dimensional reduction for anchoring. Default is "lsi".
#' @return A list containing the updated query object and calculated transfer anchors.
run_label_transfer <- function(ref_obj, query_obj, ref_assay = "RNA", query_assay = "RNA", weight_reduction = "lsi") {

  message("--> Finding cross-modal transfer anchors via CCA...")
  transfer_anchors <- FindTransferAnchors(
    reference       = ref_obj,
    query           = query_obj,
    features        = VariableFeatures(object = ref_obj),
    reference.assay = ref_assay,
    query.assay     = query_assay,
    reduction       = "cca"
  )

  message("--> Projecting reference cell type labels onto query dataset...")
  celltype_predictions <- TransferData(
    anchorset        = transfer_anchors,
    refdata          = ref_obj$celltype,
    dims             = 1:30,
    weight.reduction = query_obj[[weight_reduction]]
  )

  query_obj <- AddMetaData(query_obj, metadata = celltype_predictions)

  return(list(query = query_obj, anchors = transfer_anchors))
}

# ==========================================
# 2. Cross-Modal Co-Embedding Module
# ==========================================

#' Construct Cross-Modal Joint Co-Embedding Layer
#' @param ref_obj Seurat reference object.
#' @param query_obj Processed Seurat query object.
#' @param anchors Transfer anchors object generated from run_label_transfer.
#' @param weight_reduction Reduction slot used for anchor weights. Default is "lsi".
#' @return An integrated co-embedded Seurat object.
run_coembedding <- function(ref_obj, query_obj, anchors, weight_reduction = "lsi") {

  message("--> Performing RNA data imputation using transfer anchors...")
  genes_use <- VariableFeatures(ref_obj)
  ref_data  <- GetAssayData(ref_obj, assay = "RNA", layer = "data")[genes_use, ]

  imputation <- TransferData(
    anchorset        = anchors,
    refdata          = ref_data,
    weight.reduction = query_obj[[weight_reduction]],
    dims             = 1:30
  )

  # Isolate imputed matrix inside an independent slot to prevent baseline RNA corruption
  query_obj[["ImputedRNA"]] <- CreateAssayObject(data = imputation)
  DefaultAssay(query_obj)   <- "ImputedRNA"

  # Subset features to ensure uniform dimensional structural merging
  ref_sub   <- subset(ref_obj, features = genes_use)
  query_sub <- subset(query_obj, features = genes_use)

  message("--> Merging cross-modal layers and executing dimension reduction...")
  coembed <- merge(x = ref_sub, y = query_sub)

  # Unified celltype mapping (Prioritize ground-truth; fall back to predictions)
  coembed$celltype <- ifelse(!is.na(coembed$celltype), coembed$celltype, coembed$predicted.id)

  # Define technical modalities based on cohort metadata patterns
  coembed$tech <- case_when(
    coembed$sample %in% c("MT007D","MT790D","MT986D","MT007L","MT790L","MT986L") ~ "scRNA",
    coembed$sample %in% c("MA007D","MA790D","MA986D","MA007L","MA790L","MA986L") ~ "scATAC",
    TRUE ~ "Unknown"
  )

  # Standard downstream reduction for co-embedding visualization
  coembed <- ScaleData(coembed, features = genes_use, do.scale = FALSE, verbose = FALSE)
  coembed <- RunPCA(coembed, features = genes_use, verbose = FALSE)
  coembed <- RunUMAP(coembed, dims = 1:30, verbose = FALSE)

  return(coembed)
}

# ==========================================
# 3. Pipeline Sandbox Execution Run
# ==========================================

# File path allocations
REF_PATH   <- Sys.getenv("ATLAS_RNA_REFERENCE")
QUERY_PATH <- Sys.getenv("ATLAS_ATAC_QUERY")

if (!file.exists(REF_PATH) || !file.exists(QUERY_PATH)) {
  stop("Set ATLAS_RNA_REFERENCE and ATLAS_ATAC_QUERY to existing RDS files.")
}

scrna <- readRDS(REF_PATH)
sc    <- readRDS(QUERY_PATH)

# Run integration flow
transfer_results <- run_label_transfer(scrna, sc, weight_reduction = "lsi")
sc               <- transfer_results$query
coembed          <- run_coembedding(scrna, sc, transfer_results$anchors, weight_reduction = "lsi")

# ==========================================
# 4. Data Visualization Exports
# ==========================================

message("--> Exporting high-resolution diagnostic plots...")

# 4.1 Query dataset multi-dimensional integration overview
p1 <- DimPlot(sc, reduction = "umap", group.by = "sample", label = TRUE) + theme(legend.position = "none")
p2 <- DimPlot(sc, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + theme(legend.position = "none")
p3 <- DimPlot(sc, reduction = "umap", group.by = "group", label = TRUE, repel = TRUE) + theme(legend.position = "none")
p4 <- DimPlot(sc, reduction = "umap", group.by = "predicted.id", label = TRUE)
combined_plot <- plot_grid(p1, p2, p3, p4, ncol = 2, nrow = 2)
ggsave("query_integration_check.pdf", combined_plot, width = 12, height = 10)

# 4.2 Cross-modal joint co-embedding UMAP distributions
pe1 <- DimPlot(coembed, reduction = "umap", group.by = "tech", raster = TRUE, colors.ggplot = TRUE)
pe2 <- DimPlot(coembed, reduction = "umap", group.by = "celltype", label = TRUE, raster = TRUE)
coembed_plot <- plot_grid(pe1, pe2, ncol = 2, nrow = 1)
ggsave("umap_atac_rna_coembed.pdf", coembed_plot, width = 15, height = 7, dpi = 300)

# 4.3 Confusion matrix alignment heatmap (Evaluated only if query carries native annotations)
if (!all(is.na(sc$celltype))) {
  predictions <- table(sc$celltype, sc$predicted.id)
  predictions <- as.data.frame(predictions / rowSums(predictions))

  p_heat <- ggplot(predictions, aes(Var1, Var2, fill = Freq)) +
    geom_tile(color = "white") +
    scale_fill_gradientn(name = "Fraction of cells", colors = c("#ffffc8", "#f7a583", "#7d0025")) +
    labs(x = "Cell type annotation (RNA)", y = "Predicted cell type label (ATAC)") +
    theme_cowplot() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      panel.grid.major = element_blank()
    )
  ggsave("label_transfer_confusion_matrix.pdf", p_heat, width = 10, height = 8, dpi = 300)
}

message("--> Complete multi-omic analytics workflow finished successfully.")
