# Spatiotemporal and Cis-Regulatory Dynamics in Bovine Tissues: A Multi-Omics Atlas

This repository contains the official production-ready core pipeline and analytical scripts for our manuscript submitted to *Nature Communications*. The workflow integrates single-cell RNA-seq (scRNA-seq), single-nucleus ATAC-seq (snATAC-seq), and spatial transcriptomics (10x Genomics Visium) to profile the cis-regulatory dynamics and cellular heterogeneity in *Bos taurus* mammary gland and liver tissues.

---

## 🗺️ Pipeline Architecture Overview

The multi-omics analysis platform is structured into four main sequential modules:

1. **Spatial Preprocessing & Topology Analysis** (`spatial_pipeline.py`)
   * Visium data QC (filtering mitochondrial/ribosomal variants via NCBI style mapping).
   * Dimensionality reduction, Leiden clustering, and differential expression matrix parsing.
   * Spatial neighborhood enrichment and network centrality scoring via `Squidpy`.
2. **Cell-Type Spatial Deconvolution** (`run_tangram.py`)
   * Deep learning-based projection of annotated scRNA-seq references onto spatial coordinates via `Tangram`.
   * Mapping cell-type specific probability densities across local histological niches.
3. **Advanced Transcription Factor Mining** (`tf_advanced_analysis.R`)
   * Calculation of per-cell transcription factor (TF) accessibility deviations via `chromvar`.
   * Identification of trajectory-associated dynamic TFs via Spearman correlations along pseudotime vectors.
   * High-resolution genomic footprinting track plotting (`Signac`/`Footprint`) using custom local JASPAR2024 PFM structures.
4. **Gene Regulatory Network (GRN) Inference** (`run_scenicplus.py`)
   * Enhancer-driven GRN modeling based on cross-modal associations (TF-to-Region & Region-to-Gene) via `SCENIC+`.
   * Assembly and evaluation of multi-omic dynamic regulons (eRegulons).
