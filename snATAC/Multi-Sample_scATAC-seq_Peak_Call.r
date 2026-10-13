#' @title Automated Multi-Sample scATAC-seq Peak Call Synthesis & Quant Map Pipeline
#' @description Merges distinct genomic coordinates into a single consensus peak dictionary,
#' quantifies Tn5 transposase insertion insertions across cohorts, and initializes integrated workspace slots.
#' @author Chen Ao
#' @date 2026-06-28

library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(GenomicRanges)
library(tidyverse)
library(future)

# ==========================================
# 1. High-Performance Runtime Profiles
# ==========================================

# Use multisession instead of deprecated multiprocess
plan("multisession", workers = 4)
options(future.globals.maxSize = 50 * 1024^3) # 50 GB RAM allocation allocation

# Directory mapping profiles
SET_DIR    <- Sys.getenv("ATLAS_ATAC_DATA_DIR")
if (!dir.exists(SET_DIR)) stop("Set ATLAS_ATAC_DATA_DIR to the input directory.")
DATA_PATHS <- c(
  "./Liver/230346A_GA007A_A/outs",
  "./Liver/230346A_GA007B_A/outs",
  "./Liver/230346A_GA790A_A/outs",
  "./Liver/230346A_GA790B_A/outs",
  "./Liver/230346A_GA986A_A/outs",
  "./Liver/230346A_GA986B_A/outs"
)
SAMPLE_IDS <- c("LA007A", "LA007B", "LA790A", "LA790B", "LA986A", "LA986B")

setwd(SET_DIR)

# ==========================================
# 2. Master Coordinate Synthesis Functions
# ==========================================

#' Extract and Merge Consensus Peaks Across Experimental Batches
#' @param paths_vector Character vector containing absolute directory target locations.
#' @param min_len Minimum acceptable fragment segment width filter.
#' @param max_len Maximum acceptable fragment segment width filter.
#' @param mode Coordinate intersection execution mode ("reduce" or "disjoin").
#' @return A filtered genomic GRanges data structure.
create_consensus_peaks <- function(paths_vector, min_len = 20, max_len = 10000, mode = "reduce") {
  message("--> Parsing dynamic target bed maps across cohorts...")

  gr_list <- lapply(paths_vector, function(dir_path) {
    bed_file <- list.files(path = dir_path, pattern = "_peaks.bed$", full.names = TRUE)
    if (length(bed_file) != 1) stop(paste("Missing targeted peak file within slot destination:", dir_path))

    bed_data <- read.table(bed_file, col.names = c("chr", "start", "end"))
    # BED starts are zero-based; GRanges starts are one-based.
    return(makeGRangesFromDataFrame(bed_data, starts.in.df.are.0based = TRUE))
  })

  # Structural collapsing
  unified_ranges <- Reduce(f = c, x = gr_list)

  message(paste("--> Synthesizing peak intervals utilizing mode set to:", mode))
  consensus_peaks <- switch(mode,
    "reduce"  = GenomicRanges::reduce(unified_ranges),
    "disjoin" = GenomicRanges::disjoin(unified_ranges),
    stop("Requested structural operation profile is invalid.")
  )

  # Filter length criteria constraints
  peak_widths <- width(consensus_peaks)
  consensus_peaks <- consensus_peaks[peak_widths >= min_len & peak_widths <= max_len]

  # Strict Bos taurus autosomal filtering restriction
  bovine_chroms <- c(as.character(1:29), "X")
  consensus_peaks <- consensus_peaks[seqnames(consensus_peaks) %in% bovine_chroms]

  return(consensus_peaks)
}

# ==========================================
# 3. Seurat Matrix Matrix Compilation Module
# ==========================================

#' Process Matrix Quantifications for an Isolated Run Sample
#' @param data_dir Source path target character string.
#' @param project_id Named cohort assignment label.
#' @param master_peaks Consensus reference genomic interval profile.
#' @param cutoff_metric Quantitative minimum sequencing depth filtering assignment threshold.
#' @return Fully instantiated Seurat container object.
quantify_sample_features <- function(data_dir, project_id, master_peaks, cutoff_metric = 500) {
  message(paste("--> Processing Workspace Entity Identification Profile:", project_id))

  frag_path <- list.files(path = data_dir, pattern = "fragments\\.tsv\\.gz$", full.names = TRUE)
  meta_path <- list.files(path = data_dir, pattern = "singlecell\\.csv$", full.names = TRUE)

  if (length(frag_path) != 1 || length(meta_path) != 1) {
    stop(paste("Prerequisite resource artifacts missing from path destination:", data_dir))
  }

  # Load and trim metrics frame metadata slots
  meta_df <- read.csv(meta_path, header = TRUE, row.names = 1, stringsAsFactors = FALSE)[-1, ]
  meta_df <- meta_df[meta_df$passed_filters > cutoff_metric, ]

  # Fragment indexing guard confirmation checklist
  if (!file.exists(paste0(frag_path, ".tbi"))) {
    warning(paste("Tabix index was absent for reference path target. Indexing now:", frag_path))
    Signac::IndexFragments(frag_path)
  }

  fragment_instance <- CreateFragmentObject(path = frag_path, cells = rownames(meta_df))

  # Quantify insertion count densities over consensus vectors
  counts_matrix <- FeatureMatrix(
    fragments = fragment_instance,
    features  = master_peaks,
    cells     = rownames(meta_df)
  )

  chromatin_assay <- CreateChromatinAssay(counts = counts_matrix, fragments = fragment_instance)

  seurat_obj <- CreateSeuratObject(
    counts    = chromatin_assay,
    assay     = "peaks",
    project   = project_id,
    meta.data = meta_df
  )

  # Ensure clean identifier headers
  seurat_obj <- RenameCells(seurat_obj, add.cell.id = project_id)
  return(seurat_obj)
}

# ==========================================
# 4. Pipeline Execution Script Flow
# ==========================================

# 4.1 Generate Consensus Genomic Coordinates Dictionary
consensus_peaks <- create_consensus_peaks(DATA_PATHS, mode = "reduce")

# 4.2 Map and Compile Features Across Individual Files
seurat_list <- purrr::map2(
  .x           = DATA_PATHS,
  .y           = SAMPLE_IDS,
  .f           = ~ quantify_sample_features(data_dir = .x, project_id = .y, master_peaks = consensus_peaks)
)

# 4.3 Consolidate Matrix Collections into a Single Unified Object Workspace
message("--> Consolidating individual runs into a single workspace...")
combined_sc_atac <- Reduce(f = function(x, y) merge(x, y), x = seurat_list)

# Validate current operational metrics summary report
print(combined_sc_atac)
message("--> Data structuring sequence finalized successfully.")
