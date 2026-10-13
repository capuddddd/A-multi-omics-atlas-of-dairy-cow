#' @title Merge Single-Cell h5 Datasets
#' @description Load multiple 10X .h5 files, create Seurat objects, and merge them into a single dataset.
#' @author Chen Ao
#' @date 2026-06-28

library(Seurat)

# ==========================================
# 1. Core Functions
# ==========================================

#' Load h5 datasets and create a list of Seurat objects
load_h5_datasets <- function(data_dir, pattern = "\\.h5$") {
  if (!dir.exists(data_dir)) {
    stop(paste("Data directory does not exist:", data_dir))
  }

  fs <- list.files(data_dir, pattern = pattern, full.names = TRUE)
  if (length(fs) == 0) {
    stop(paste("No .h5 files found in:", data_dir))
  }

  sce_list <- lapply(fs, function(file_path) {
    file_name <- basename(file_path)
    message("Processing: ", file_name)

    # Read 10X h5 matrix
    counts <- Read10X_h5(file_path)

    # Keep the complete file stem; underscores may be part of a sample ID.
    sample_id <- sub("\\.h5$", "", file_name)

    if (is.list(counts)) {
      if (!"Gene Expression" %in% names(counts)) stop("Missing Gene Expression: ", file_path)
      counts <- counts[["Gene Expression"]]
    }

    # Create Seurat object
    sce <- CreateSeuratObject(counts = counts, project = sample_id)
    sce$sample <- sample_id

    return(sce)
  })

  return(sce_list)
}

#' Merge Seurat object list and synchronize identities
merge_seurat_lists <- function(sce_list, cell_ids, project_name) {
  if (length(sce_list) != length(cell_ids)) {
    stop("The number of Seurat objects does not match the length of cell_ids.")
  }

  if (length(sce_list) == 0 || anyNA(cell_ids) || any(!nzchar(cell_ids)) || anyDuplicated(cell_ids)) {
    stop("Provide a nonempty object list with unique, nonempty sample IDs.")
  }
  for (i in seq_along(sce_list)) sce_list[[i]]$sample <- cell_ids[[i]]
  if (length(sce_list) == 1) {
    merged_sc <- RenameCells(sce_list[[1]], add.cell.id = cell_ids[[1]])
  } else {
    merged_sc <- merge(x = sce_list[[1]], y = sce_list[-1],
                       add.cell.ids = cell_ids, project = project_name)
  }
  Idents(merged_sc) <- "sample"

  return(merged_sc)
}

# ==========================================
# 2. Execution Pipeline
# ==========================================

# Run with: Rscript snRNA/Merge_Single-Cell_h5_Datasets.r INPUT_DIR OUTPUT_RDS
main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 2) stop("Usage: Rscript Merge_Single-Cell_h5_Datasets.r INPUT_DIR OUTPUT_RDS")
  sce_list <- load_h5_datasets(args[[1]])
  cell_ids <- vapply(sce_list, function(sce) as.character(sce$sample[[1]]), character(1))
  merged_sc <- merge_seurat_lists(sce_list, cell_ids, "MT_integrated")
  dir.create(dirname(args[[2]]), recursive = TRUE, showWarnings = FALSE)
  saveRDS(merged_sc, args[[2]])
  message("Saved merged object: ", args[[2]])
}

if (sys.nframe() == 0) main()
