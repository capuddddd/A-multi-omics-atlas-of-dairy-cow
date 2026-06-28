#' @title Multi-Condition CellChat Interaction Comparison and Visualization Pipeline
#' @description Comprehensive workflow for multi-sample comparative cell-cell communication 
#' utilizing CellChat, including custom network centrality heatmaps and customized ggplot tile plots.
#' @author Chen Ao
#' @date 2026-06-28

library(Seurat)
library(tidyverse)
library(CellChat)
library(NMF)
library(ggalluvial)
library(patchwork)
library(ComplexHeatmap)
library(reshape2)

# ==========================================
# 1. Custom CellChat Customization Functions
# ==========================================

#' Optimized Custom Heatmap Function for Signaling Role Analysis
#' @description Tweaked version of netAnalysis_signalingRole_heatmap supporting strict row/column clustering toggles.
my_netAnalysis_signalingRole_heatmap <- function(object, signaling = NULL, pattern = c("outgoing", "incoming", "all"), slot.name = "netP",
                                                 color.use = NULL, color.heatmap = "BuGn", title = NULL, width = 10, height = 8, 
                                                 font.size = 8, font.size.title = 10, cluster.rows = FALSE, cluster.cols = FALSE) {
  pattern <- match.arg(pattern)
  if (length(slot(object, slot.name)$centr) == 0) {
    stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores!")
  }
  
  centr <- slot(object, slot.name)$centr
  outgoing <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
  incoming <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
  dimnames(outgoing) <- list(levels(object@idents), names(centr))
  dimnames(incoming) <- dimnames(outgoing)
  
  for (i in 1:length(centr)) {
    outgoing[, i] <- centr[[i]]$outdeg
    incoming[, i] <- centr[[i]]$indeg
  }
  
  if (pattern == "outgoing") {
    mat <- t(outgoing)
    legend.name <- "Outgoing"
  } else if (pattern == "incoming") {
    mat <- t(incoming)
    legend.name <- "Incoming"
  } else if (pattern == "all") {
    mat <- t(outgoing + incoming)
    legend.name <- "Overall"
  }
  
  title <- if (is.null(title)) paste0(legend.name, " signaling patterns") else paste0(legend.name, " signaling patterns - ", title)

  if (!is.null(signaling)) {
    mat1 <- mat[rownames(mat) %in% signaling, , drop = FALSE]
    mat <- matrix(0, nrow = length(signaling), ncol = ncol(mat))
    idx <- match(rownames(mat1), signaling)
    mat[idx[!is.na(idx)], ] <- mat1
    dimnames(mat) <- list(signaling, colnames(mat1))
  }
  
  mat.ori <- mat
  mat <- sweep(mat, 1L, apply(mat, 1, max), '/', check.margin = FALSE)
  mat[mat == 0] <- NA

  if (is.null(color.use)) {
    color.use <- scPalette(length(colnames(mat)))
  }
  color.heatmap.use <- grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, name = color.heatmap)))(102)

  df <- data.frame(group = colnames(mat))
  rownames(df) <- colnames(mat)
  names(color.use) <- colnames(mat)
  
  col_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use), which = "column",
                                      show_legend = FALSE, show_annotation_name = FALSE,
                                      simple_anno_size = grid::unit(0.2, "cm"))
  ha2 <- HeatmapAnnotation(Strength = anno_barplot(colSums(mat.ori), border = FALSE, gp = gpar(fill = color.use, col = color.use)), show_annotation_name = FALSE)

  pSum <- rowSums(mat.ori)
  pSum.original <- pSum
  pSum <- -1 / log(pSum)
  pSum[is.na(pSum)] <- 0
  idx1 <- which(is.infinite(pSum) | pSum < 0)
  if (length(idx1) > 0) {
    values.assign <- seq(max(pSum) * 1.1, max(pSum) * 1.5, length.out = length(idx1))
    position <- sort(pSum.original[idx1], index.return = TRUE)$ix
    pSum[idx1] <- values.assign[match(1:length(idx1), position)]
  }

  ha1 <- rowAnnotation(Strength = anno_barplot(pSum, border = FALSE), show_annotation_name = FALSE)

  legend.break <- if (min(mat, na.rm = TRUE) == max(mat, na.rm = TRUE)) max(mat, na.rm = TRUE) else c(round(min(mat, na.rm = TRUE), digits = 1), round(max(mat, na.rm = TRUE), digits = 1))
  
  ht1 <- Heatmap(mat, col = color.heatmap.use, na_col = "white", name = "Relative strength",
                 bottom_annotation = col_annotation, top_annotation = ha2, right_annotation = ha1,
                 cluster_rows = cluster.rows, cluster_columns = cluster.cols,
                 row_names_side = "left", row_names_rot = 0, row_names_gp = gpar(fontsize = font.size), column_names_gp = gpar(fontsize = font.size),
                 width = unit(width, "cm"), height = unit(height, "cm"),
                 column_title = title, column_title_gp = gpar(fontsize = font.size.title), column_names_rot = 90,
                 heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "plain"), title_position = "leftcenter-rot",
                                             border = NA, at = legend.break, legend_height = unit(20, "mm"), labels_gp = gpar(fontsize = 8), grid_width = unit(2, "mm"))
  )
  return(ht1)
}

# ==========================================
# 2. Main Processing Function
# ==========================================

#' Standard Execution Pipeline for Multi-Condition Comparison
run_comparative_cellchat <- function(rds_path, output_dir, gene_translation_df = NULL) {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # 2.1 Load and Split Seurat Object by Metadata Group Conditons
  sc <- readRDS(rds_path)
  
  # Drop unused identity levels from subset workspace
  if ("subtype" %in% colnames(sc@meta.data)) {
    sc$subtype <- droplevels(sc$subtype)
  }
  
  Idents(sc) <- "group"
  sco.A <- subset(sc, idents = 'Dry')
  sco.B <- subset(sc, idents = 'Lac')
  
  # Apply Gene Translation Matrix mappings if provided
  if (!is.null(gene_translation_df)) {
    sco.A <- RenameGenesSeurat_v2(obj = sco.A, newnames = gene_translation_df$x1, de.assay = "RNA")
    sco.B <- RenameGenesSeurat_v2(obj = sco.B, newnames = gene_translation_df$x1, de.assay = "RNA")
  }
  
  # 2.2 Initialize Independent Base CellChat Arrays
  construct_individual_cellchat <- function(seurat_sub) {
    cc <- createCellChat(GetAssayData(seurat_sub, assay = "RNA", layer = "data"), meta = seurat_sub@meta.data, group.by = "celltype")
    cc@DB <- CellChatDB.human
    cc <- subsetData(cc)
    cc <- identifyOverExpressedGenes(cc)
    cc <- identifyOverExpressedInteractions(cc)
    cc <- smoothData(cc, adj = PPI.human)
    cc <- computeCommunProb(cc, raw.use = TRUE, population.size = TRUE)
    cc <- filterCommunication(cc, min.cells = 10)
    cc <- computeCommunProbPathway(cc)
    cc <- aggregateNet(cc)
    cc <- netAnalysis_computeCentrality(cc, slot.name = "netP")
    return(cc)
  }
  
  message("Executing CellChat pipeline for Condition: Dry")
  sca <- construct_individual_cellchat(sco.A)
  message("Executing CellChat pipeline for Condition: Lac")
  scb <- construct_individual_cellchat(sco.B)
  
  saveRDS(sca, file.path(output_dir, "sca.rds"))
  saveRDS(scb, file.path(output_dir, "scb.rds"))
  
  # 2.3 Cross-Dataset Merge Analytics
  cco.list <- list(Dry = sca, Lac = scb)
  cellchat <- mergeCellChat(cco.list, add.names = names(cco.list))
  saveRDS(cellchat, file.path(output_dir, "cellchat_merged.rds"))
  
  # --- 2.4 Structural Figure Exports ---
  
  # Interaction Number vs Strength Metrics Comparison
  gg1 <- compareInteractions(cellchat, show.legend = FALSE, group = c(1,2), measure = "count")
  gg2 <- compareInteractions(cellchat, show.legend = FALSE, group = c(1,2), measure = "weight")
  ggsave(file.path(output_dir, "Overview_number_strength.png"), plot = (gg1 + gg2), width = 6, height = 4)
  
  # Differential Interaction Nets (Circle Layouts)
  pdf(file.path(output_dir, "Diff_number_strength_net.pdf"))
  netVisual_diffInteraction(cellchat, weight.scale = TRUE, measure = "count")
  netVisual_diffInteraction(cellchat, weight.scale = TRUE, measure = "weight")
  dev.off()
  
  # Differential Heatmaps
  png(file.path(output_dir, "Diff_number_strength_heatmap.png"), width = 2000, height = 1000, res = 300)
  h1 <- netVisual_heatmap(cellchat, measure = "count")
  h2 <- netVisual_heatmap(cellchat, measure = "weight")
  print(h1 + h2)
  dev.off()
  
  # Quantified Interaction Weight Networks Across Full Object Matrices
  pdf(file.path(output_dir, "Counts_Compare_net.pdf"))
  weight.max <- getMaxWeight(cco.list, attribute = c("idents", "count"))
  par(mfrow = c(1, 2))
  for (i in 1:length(cco.list)) {
    netVisual_circle(cco.list[[i]]@net$count, weight.scale = TRUE, label.edge = FALSE,
                     edge.weight.max = weight.max[2], edge.width.max = 14,
                     title.name = paste0("Number of interactions - ", names(cco.list)[i]))
  }
  dev.off()
  
  # Target Specific Structural Network Filtering Visualizations
  pdf(file.path(output_dir, "Counts_Compare_select.pdf"))
  selected_cells <- intersect(c(
    "Liver sinusoidal endothelial cells", "Kuffer cells", "Hepatocytes", "Cytotoxic CD8+ T cells", 
    "Cytotoxic CD4+ T cells", "NK cells", "Stellate cells", "B cells", "γδT cells", 
    "Cholangiocytes", "Dendritic cells", "Blood vasular endothelial cells", "Cycling immune cells"
  ), levels(cellchat@idents))
  
  if(length(selected_cells) > 1) {
    count1 <- cco.list[[1]]@net$count[selected_cells, selected_cells, drop = FALSE]
    count2 <- cco.list[[2]]@net$count[selected_cells, selected_cells, drop = FALSE]
    weight.max.sel <- max(max(count1), max(count2))
    par(mfrow = c(1, 2))
    netVisual_circle(count1, weight.scale = TRUE, label.edge = TRUE, edge.weight.max = weight.max.sel, edge.width.max = 12, title.name = paste0("Interactions - ", names(cco.list)[1]))
    netVisual_circle(count2, weight.scale = TRUE, label.edge = TRUE, edge.weight.max = weight.max.sel, edge.width.max = 12, title.name = paste0("Interactions - ", names(cco.list)[2]))
  }
  dev.off()
  
  # Pathway Signaling Group Strengths Comparative Ranks
  col_palette <- c("#ffc0cb", "#61b4cf")
  gg_rank1 <- rankNet(cellchat, mode = "comparison", color.use = col_palette, stacked = TRUE, do.stat = TRUE)
  gg_rank2 <- rankNet(cellchat, mode = "comparison", color.use = col_palette, stacked = FALSE, do.stat = TRUE)
  ggsave(file.path(output_dir, "Compare_pathway_strength.png"), plot = (gg_rank1 + gg_rank2), width = 10, height = 6)
  
  # Outgoing vs Incoming Role Dynamics Mapped to Unified ComplexHeatmaps
  pathway.union <- union(cco.list[[1]]@netP$pathways, cco.list[[2]]@netP$pathways)
  
  pdf(file.path(output_dir, "Compare_signal_pattern_all.pdf"), width = 10, height = 8)
  ht_all1 <- netAnalysis_signalingRole_heatmap(cco.list[[1]], pattern = "all", signaling = pathway.union, title = names(cco.list)[1], width = 8, height = 16)
  ht_all2 <- netAnalysis_signalingRole_heatmap(cco.list[[2]], pattern = "all", signaling = pathway.union, title = names(cco.list)[2], width = 8, height = 16)
  draw(ht_all1 + ht_all2, ht_gap = unit(0.5, "cm"))
  dev.off()
  
  png(file.path(output_dir, "Compare_signal_pattern_outgoing.png"), width = 3000, height = 2500, res = 300)
  ht_out1 <- netAnalysis_signalingRole_heatmap(cco.list[[1]], pattern = "outgoing", signaling = pathway.union, title = names(cco.list)[1], width = 8, height = 14)
  ht_out2 <- netAnalysis_signalingRole_heatmap(cco.list[[2]], pattern = "outgoing", signaling = pathway.union, title = names(cco.list)[2], width = 8, height = 14)
  draw(ht_out1 + ht_out2, ht_gap = unit(0.5, "cm"))
  dev.off()
  
  png(file.path(output_dir, "Compare_signal_pattern_incoming.png"), width = 3000, height = 2500, res = 300)
  ht_in1 <- my_netAnalysis_signalingRole_heatmap(cco.list[[1]], pattern = "incoming", signaling = pathway.union, title = names(cco.list)[1], width = 8, height = 14)
  ht_in2 <- my_netAnalysis_signalingRole_heatmap(cco.list[[2]], pattern = "incoming", signaling = pathway.union, title = names(cco.list)[2], width = 8, height = 14)
  draw(ht_in1 + ht_in2, ht_gap = unit(0.5, "cm"))
  dev.off()
  
  # Targeted Functional Bubble Expression Profiles
  p_bub1 <- netVisual_bubble(cellchat, sources.use = 1:8, targets.use = 1:8, comparison = c(1, 2), angle.x = 45)
  ggsave(file.path(output_dir, "Compare_test_bubble1.pdf"), plot = p_bub1, width = 12, height = 10)
  
  # Condition-Specific Enriched Signaling Patterns (DEGs Mapping Approach)
  pos.dataset <- "Lac"
  features.name <- paste0(pos.dataset, ".merged")
  cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.05, thresh.p = 0.05, group.DE.combined = FALSE) 
  
  net_deg <- netMappingDEG(cellchat, features.name = features.name, variable.all = TRUE)
  net.up <- subsetCommunication(cellchat, net = net_deg, datasets = "Lac", ligand.logFC = 0.05, receptor.logFC = NULL)
  net.down <- subsetCommunication(cellchat, net = net_deg, datasets = "Dry", ligand.logFC = -0.05, receptor.logFC = NULL)
  
  if(nrow(net.up) > 0) {
    gg_up <- netVisual_bubble(cellchat, pairLR.use = net.up[, "interaction_name", drop = FALSE], sources.use = 1:9, targets.use = 1:9, comparison = c(1, 2), angle.x = 90, remove.isolate = TRUE, title.name = "Up-regulated signaling")
    ggsave(file.path(output_dir, "Up_regulated_bubble.pdf"), plot = gg_up, width = 15, height = 12)
  }
}

# ==========================================
# 3. Custom Ggplot Matrix Tile Plot Call
# ==========================================

#' Customized Tile Heatmap Plotting Routine Matrix Outputs
plot_custom_tile_heatmap <- function(cellchat_obj, output_path) {
  count_matrix <- cellchat_obj@net$count
  weight_matrix <- cellchat_obj@net$weight
  
  count_long <- melt(count_matrix)
  colnames(count_long) <- c("Source", "Target", "Count")
  weight_long <- melt(weight_matrix)
  colnames(weight_long) <- c("Source", "Target", "Weight")
  
  merged_data <- merge(count_long, weight_long, by = c("Source", "Target"))
  
  p_tile <- ggplot(data = merged_data, aes(x = Target, y = Source)) +
    geom_tile(aes(fill = Weight, width = Count / 200, height = Count / 200)) +
    scale_fill_gradient2(
      low = "white", mid = "#FF7F7F", high = "#8B0000",
      midpoint = median(merged_data$Weight), limits = c(3, 17)
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(title = "Interaction Heatmap Matrix", x = "Target Cell", y = "Source Cell", fill = "Weight") +
    geom_vline(xintercept = seq(1.5, length(unique(merged_data$Target)) - 0.5, by = 1), size = 0.3) +
    geom_hline(yintercept = seq(1.5, length(unique(merged_data$Source)) - 0.5, by = 1), size = 0.3)
  
  ggsave(output_path, plot = p_tile, width = 8, height = 7)
}