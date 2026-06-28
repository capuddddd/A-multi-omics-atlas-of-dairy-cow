#' @title Merge Single-Cell h5 Datasets
#' @description Load multiple 10X .h5 files, create Seurat objects, and merge them into a single dataset.
#' @author Chen Ao
#' @date 2026-06-28

library(Seurat)
library(stringr)

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
    
    # Extract sample ID (prefix before the first underscore)
    sample_id <- str_split(file_name, "_", simplify = TRUE)[1, 1]
    
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
  
  # Merge all Seurat objects in the list
  merged_sc <- merge(
    x = sce_list[[1]],
    y = sce_list[-1],
    add.cell.ids = cell_ids,
    project = project_name
  )
  
  # Align active idents with cell_ids and update metadata
  names(cell_ids) <- levels(merged_sc)
  merged_sc <- RenameIdents(merged_sc, cell_ids)
  merged_sc$sample <- Idents(merged_sc)
  
  return(merged_sc)
}

# ==========================================
# 2. Execution Pipeline
# ==========================================

# Define directories (Replace with your actual paths or configuration setup)
DATA_DIR   <- "./data/Mammary_gland/merge_raw_h5" 
OUTPUT_DIR <- "./results/mammary_gland"

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# Run processing workflow
message("--- Starting Single-Cell Data Merging ---")
sce_list <- load_h5_datasets(data_dir = DATA_DIR)

cell_ids <- c("MT007A", "MT007B", "MT790A", "MT790B", "MT986A", "MT986B")
merged_sc <- merge_seurat_lists(
  sce_list     = sce_list, 
  cell_ids     = cell_ids, 
  project_name = "MT_integrated"
)

# Save output object
output_file <- file.path(OUTPUT_DIR, "Mamm_merge.rds")
saveRDS(merged_sc, file = output_file)
message("Pipeline completed successfully. Output saved to: ", output_file)