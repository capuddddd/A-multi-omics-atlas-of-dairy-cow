"""
SCENIC+ Single-Cell Multi-Omics Gene Regulatory Network Inference Pipeline
Environment Requirement: scenicplus (Python >= 3.8)
Author: Chen Ao
Date: 2026-06-28
"""


def main():
    """Run the legacy template after validating the required upstream version."""
    import os
    import scanpy as sc
    import numpy as np
    import pandas as pd

    # Import the core SCENIC+ components.
    import scenicplus
    from scenicplus.scenicplusclass import create_SCENICPLUS_object
    from scenicplus.preprocessing.filtering import filter_genes, filter_regions
    from scenicplus.wrappers.run_scenicplus import run_scenicplus

    # ==========================================
    # 1. Initialize the environment and parameters
    # ==========================================
    WORKDIR = "./scenicplus_output"
    if not os.path.exists(WORKDIR):
        os.makedirs(WORKDIR)

    os.chdir(WORKDIR)

    # Set the number of parallel workers.
    N_CORES = 16

    # ==========================================
    # 2. Load the multi-omics input matrices
    # ==========================================
    print("--> Loading processed scRNA-seq and scATAC-seq matrices...")
    adata_rna = sc.read_h5ad("../rna_processed.h5ad")
    adata_atac = sc.read_h5ad("../atac_processed.h5ad")

    # Represent genomic coordinates in chr-start-end format.
    # Convert "chr1:1000-2000" to "chr1-1000-2000" when necessary.
    adata_atac.var_names = adata_atac.var_names.str.replace(":", "-")

    # ==========================================
    # 3. Configure candidate regulatory regions and associations
    # ==========================================
    # SCENIC+ requires species-appropriate motif resources; bovine/vertebrate examples follow.
    # Download the corresponding reference resources before running:
    # 1. cisTarget.regions.feather (region scoring matrix)
    # 2. motifs-v10nr_clust-nr.vertebrates-m0.001-o0.0.tbl (motif annotation table)

    db_annot_fn = "../motifs-v10nr_clust-nr.vertebrates-m0.001-o0.0.tbl"
    ctx_db_fn = "../bosTau9_screen_v10_nr.vertebrates.regions.feather" # Match the database to the species and assembly.

    # ==========================================
    # 4. Initialize the SCENIC+ object
    # ==========================================
    print("--> Initializing SCENIC+ master container object...")

    scplus_obj = create_SCENICPLUS_object(
        type_of_data={'RNA': 'AnnData', 'ATAC': 'AnnData'},
        data={'RNA': adata_rna, 'ATAC': adata_atac},
        multi_ome_mode=True,
        key_to_group_cells='cell_type' # Cell-type or cluster annotation column.
    )

    # Filter genes and regions with low cell coverage.
    filter_genes(scplus_obj, min_cells=5)
    filter_regions(scplus_obj, min_cells=5)

    print(scplus_obj)

    # ==========================================
    # 5. Run the SCENIC+ network inference wrapper
    # ==========================================
    print("--> Running SCENIC+ GRN computation wrapper (This will take hours)...")

    # The wrapper is intended to perform the following steps:
    # 1. Identify TF-binding regions through cisTarget motif enrichment.
    # 2. Model relationships between region accessibility and gene expression (region-to-gene).
    # 3. Model relationships between TF expression and region accessibility (TF-to-region).
    # 4. Assemble enhancer-driven regulons (eRegulons).

    try:
        run_scenicplus(
            scplus_obj=scplus_obj,
            variable_gene_keys=['highly_variable'],
            pfm_matrix_path="../jaspar2024.txt", # Local JASPAR2024 matrix file.
            species='bos_taurus',                # Target species.
            assembly='bosTau9',                  # Genome assembly.
            tf_file="../tf_names.txt",           # Species-specific list of known transcription factors.
            save_path="scplus_obj_final.pkl",
            biomart_host="http://www.ensembl.org",
            upstream=150000,                     # Maximum upstream enhancer search distance.
            downstream=150000,                   # Maximum downstream enhancer search distance.
            calculate_TF_eRegulon_specificity=True,
            calculate_Gex_eRegulon_specificity=True,
            n_cores=N_CORES
        )
        print("--> SCENIC+ Pipeline executed successfully!")
    except Exception as e:
        raise RuntimeError("SCENIC+ legacy template failed; see docs/KNOWN_LIMITATIONS.md") from e


if __name__ == "__main__":
    main()
