"""Visium QC, clustering, spatial summaries and marker plots.

Inputs are resolved from the invocation directory; outputs use --work_dir.
The fixed marker panel and QC thresholds must be reviewed for each dataset.
"""

import argparse
from pipeline_utils import input_path, output_directory, select_library
import re

def main(work_dir, adata_path, library_id=None):
    adata_path = input_path(adata_path, directory=True)
    output = output_directory(work_dir)
    import scanpy as sc
    import squidpy as sq
    from matplotlib import pyplot as plt

    # 1. Load Visium spatial transcriptomics data
    print("--> Loading Visium data...")
    adata = sq.read.visium(adata_path)
    adata.var_names_make_unique()
    library = select_library(adata, library_id)

    # 2. Quality control and mitochondrial/ribosomal gene filtering
    MT = ["ATP6", "ATP8", "COX1", "COX2", "COX3", "CYTB", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6"]
    adata.var["mt"] = adata.var_names.isin(set(MT))

    rb_genes = [gene for gene in adata.var_names if re.match(r'^RP[SL]', gene)]
    adata = adata[:, ~adata.var_names.isin(rb_genes)].copy()

    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
    adata = adata[:, ~adata.var_names.isin(MT)].copy()

    sc.pp.filter_cells(adata, min_counts=50)
    sc.pp.filter_genes(adata, min_cells=10)

    # 3. Normalize expression and reduce dimensionality
    # Preserve counts for mapping and log expression for differential testing.
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    adata.raw = adata.copy()
    sc.pp.scale(adata, max_value=10)

    sc.tl.pca(adata, svd_solver="arpack")
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=20)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=1.0)

    # 4. Analyze spatial heterogeneity and topology with Squidpy
    # 4.1 Visualize clusters in UMAP and spatial coordinates
    sc.pl.umap(adata, color=["leiden"], wspace=0.4, show=False)
    plt.savefig(output / "UMAP_leiden.pdf", bbox_inches='tight')
    plt.close()

    sq.pl.spatial_scatter(adata, color="leiden", wspace=0.4, library_id=library, show=False)
    plt.savefig(output / "Spatial_leiden_scatter.pdf", bbox_inches='tight')
    plt.close()

    # 4.2 Compute spatial centrality and neighborhood enrichment
    sq.gr.spatial_neighbors(adata, coord_type="generic", delaunay=True)
    sq.gr.centrality_scores(adata, cluster_key="leiden")
    sq.pl.centrality_scores(adata, cluster_key="leiden", show=False)
    plt.savefig(output / "centrality_scores.pdf", bbox_inches='tight')
    plt.close()

    sq.gr.spatial_neighbors(adata)
    sq.gr.nhood_enrichment(adata, cluster_key="leiden")
    sq.pl.nhood_enrichment(adata, cluster_key="leiden", show=False)
    plt.savefig(output / "neighborhood_enrichment.pdf", bbox_inches='tight')
    plt.close()

    # 5. Plot liver cell-type marker expression
    genelist = {
        "Hepatocytes": ["CYP1A2", "CYP3A74", "ALB"],
        "LSECs": ["CLEC4G", "STAB1", "STAB2", "CLEC3B"],
        "Kupffer cells": ["C1QA", "C1QB", "BOLA-DRA"],
        "Hepatic Stellate": ["DCN", "AEBP1", "IGFBP3"],
        "Cholangiocytes": ["APOA1", "ALDOB"]
    }

    # Plot marker expression from adata, grouped by Leiden clusters.
    sc.pl.dotplot(adata, genelist, groupby='leiden', dendrogram=True, cmap='bwr', show=False)
    plt.gcf().set_size_inches(10, 5)
    plt.savefig(output / 'marker_dotplot.pdf', format='pdf', bbox_inches='tight')
    plt.close()

    # 6. Differential gene expression analysis
    sc.tl.rank_genes_groups(adata, groupby='leiden', use_raw=True, method="wilcoxon")

    sc.pl.rank_genes_groups_matrixplot(adata, groupby='leiden', values_to_plot='logfoldchanges', cmap='bwr', show=False)
    plt.savefig(output / "deg_heatmap.pdf", bbox_inches='tight')
    plt.close()

    sc.pl.rank_genes_groups_dotplot(adata, groupby="leiden", standard_scale="var", n_genes=5, cmap='bwr', show=False)
    plt.gcf().set_size_inches(20, 12)
    plt.savefig(output / 'celltype_deg_dotplot.pdf', format='pdf', bbox_inches='tight')
    plt.close()

    # 7. Save processed data
    adata.write(output / "processed_data.h5ad")
    print("--> Pipeline finished successfully. Data saved to 'processed_data.h5ad'.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Spatial Transcriptomics Analysis Pipeline")
    parser.add_argument("--work_dir", type=str, required=True, help="Working directory path")
    parser.add_argument("--adata_path", type=str, required=True, help="Path to raw Visium data directory")
    parser.add_argument("--library_id", help="Visium library ID; inferred for one library")
    args = parser.parse_args()

    main(args.work_dir, args.adata_path, args.library_id)
