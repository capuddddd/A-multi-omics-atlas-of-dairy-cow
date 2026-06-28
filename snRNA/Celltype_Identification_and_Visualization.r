#' @title Single-Cell Cell Type Identification and Visualization Pipeline
#' @description Unified script for cell type identity assignment, custom dot plots, 
#' marker gene validation heatmaps, and faceted expression analysis.
#' @author Chen Ao
#' @date 2026-06-28

library(Seurat)
library(dplyr)
library(ggplot2)
library(forcats)
library(ggh4x)
library(RColorBrewer)

# ==========================================
# 1. Configuration & Marker Gene Lists
# ==========================================

# Factor level orders for visualization
MAMMARY_ORDERS <- c("LumSec", "Myoepithelial", "LumHR", "Macrophages", "Fibroblasts", "Vascular", "Lymphatic", "Myofibroblasts", "T cells")
LIVER_ORDERS   <- c("Hep", "Vascular", "Kuffer", "T cells", "NK cells", "HSC", "B cells", "Plasma", "Chol", "DC")

# Mammary Gland canonical markers
MAMMARY_MARKERS <- list(
  LumSec   = c("ELF5", "NFIB", "SLC28A3", "CCL28", "COBL"),
  Myoepi   = c("DST", "MYLK", "MYH11", "KRT17"),
  LumHR    = c("ESR1", "ESRRG", "ERBB4", "TTC6"),
  Macro    = c("CD74", "BLNK", "ZEB2", "CD68", "C1QA", "CSF1R"),
  Fibro    = c("COL1A1", "DCN", "COL3A1"),
  Vascular = c("PECAM1", "VWF", "MECOM"),
  Lymphatic= c("PROX1", "MMRN1", "KLHL4"),
  Myofibro = c("EBF1", "ADGRL3"),
  T        = c("CTSW", "IKZF1", "CD3E"),
  Proli    = c("TOP2A", "MKI67"),
  Myo      = c("ACTA2", "TAGLN"),
  NKT      = c("CTSW", "NKG7", "CD3E")
)

# Liver canonical markers (Bos taurus / Mammalian)
LIVER_MARKERS <- list(
  Hep      = c("PAH", "CPS1", "CYP3A74"),
  Vascular = c("STAB1", "STAB2", "VWF"),
  Kuffer   = c("C1QA", "ZEB2", "CSF1R"),
  T        = c("CD3E", "SKAP1", "THEMIS"),
  NK       = c("KLRB1", "KLRK1", "CTSW"),
  HSC      = c("RELN", "DCN", "COL3A1"),
  B        = c("MS4A1", "BANK1", "EBF1"),
  Plasma   = c("JCHAIN", "ITGA4", "DOCK2"),
  Chol     = c("CFTR", "PKHD1", "APOA1"),
  DC       = c("CLEC9A", "CADM1", "ID2")
)

# ==========================================
# 2. Visualization Functions
# ==========================================

#' Generate Structured DotPlot with Custom Continuous Scales
plot_custom_dotplot <- function(seurat_obj, gene_list, output_path, width = 10, height = 5) {
  # Generate base dotplot data layout matrix
  p_base <- DotPlot(seurat_obj, features = gene_list) + RotatedAxis()
  plot_data <- p_base$data
  
  # Ensure strict preservation of gene and identity ordering arrays
  plot_data$features.plot <- fct_inorder(as.factor(plot_data$features.plot))
  plot_data$id            <- fct_inorder(as.factor(plot_data$id))
  
  p_custom <- ggplot(plot_data, aes(x = features.plot, y = id)) +
    geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
    geom_point(aes(size = pct.exp, color = avg.exp.scaled), shape = 21, color = "black", stroke = 1) +
    guides(size = guide_legend(override.aes = list(shape = 21, colour = "black", fill = NA))) +
    scale_color_gradientn(
      colors = c("#5bacdb", "#62b8e3", "#c5e1f0", "white", "#f0b4c9", "#ea87aa", "#ea7ea3"),
      values = c(0, 0.25, 0.5, 0.75, 1)
    ) +
    scale_size(range = c(0, 5)) +
    theme_bw() +
    theme(
      panel.background = element_blank(),
      panel.border     = element_rect(fill = NA, color = "black"),
      panel.grid.major = element_line(color = "grey80"),
      axis.title       = element_blank(),
      axis.text.x      = element_text(color = "black", size = 12, angle = 90, hjust = 1, vjust = 0.5),
      axis.text.y      = element_text(color = "black", size = 12)
    )
  
  ggsave(output_path, plot = p_custom, width = width, height = height)
}

#' Generate Grouped Condition DotPlot across Dynamic Cellular States
plot_condition_dotplot <- function(seurat_obj, markers, cluster_order, output_path, width = 10, height = 6) {
  # Create combined state classification field
  seurat_obj$celltype2 <- paste(seurat_obj$celltype, seurat_obj$group, sep = "_")
  
  # Extract scaling base coordinates
  p_base <- DotPlot(seurat_obj, features = markers, group.by = "celltype2")
  plot_data <- p_base$data
  
  # Structural metadata alignment validation join 
  anno <- distinct(data.frame(id = seurat_obj$celltype2, group = seurat_obj$group))
  plot_data <- left_join(plot_data, anno, by = "id")
  
  # Filter missing identities and re-level factors
  plot_data <- plot_data %>% filter(id %in% cluster_order)
  plot_data$id <- factor(plot_data$id, levels = cluster_order)
  
  p_cond <- ggplot(plot_data, aes(x = features.plot, y = id, size = pct.exp, fill = avg.exp.scaled)) +
    geom_point(shape = 21, colour = "black", stroke = 1) +
    guides(size = guide_legend(override.aes = list(shape = 21, colour = "black", fill = NA))) +
    scale_fill_gradientn(colours = c("#61b4cf", "white", "#AE402F")) +
    theme_bw() +
    theme(
      panel.background = element_blank(),
      panel.border     = element_rect(fill = NA, color = "black"),
      panel.grid.major = element_line(color = "grey80"),
      axis.title       = element_blank(),
      axis.text.y      = element_text(color = "black", size = 12),
      axis.text.x      = element_text(color = "black", size = 12, angle = 45, hjust = 1, vjust = 1)
    )
  
  ggsave(output_path, plot = p_cond, width = width, height = height)
}