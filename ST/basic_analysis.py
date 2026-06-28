import argparse
import os
import re
import scanpy as sc
import squidpy as sq
from matplotlib import pyplot as plt

def main(work_dir, adata_path):
    os.chdir(work_dir)

    try:
        # 1. 载入空间转录组数据 (Visium)
        print("--> Loading Visium data...")
        adata = sq.read.visium(adata_path)
        adata.var_names_make_unique()

        # 2. 质控与过滤 (过滤线粒体与核糖体基因)
        MT = ["ATP6", "ATP8", "COX1", "COX2", "COX3", "CYTB", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6"]
        adata.var["mt"] = adata.var_names.isin(set(MT))
        
        rb_genes = [gene for gene in adata.var_names if re.match(r'^RP[SL]', gene)]
        adata = adata[:, ~adata.var_names.isin(rb_genes)]
        
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
        adata = adata[:, ~adata.var_names.isin(MT)]

        sc.pp.filter_cells(adata, min_counts=50)
        sc.pp.filter_genes(adata, min_cells=10)

        # 3. 数据标准化与降维
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        sc.pp.scale(adata, max_value=10)

        sc.tl.pca(adata, svd_solver="arpack")
        sc.pp.neighbors(adata, n_neighbors=10, n_pcs=20)
        sc.tl.umap(adata)
        sc.tl.leiden(adata, resolution=1.0)

        # 4. 空间异质性与高级拓扑分析 (Squidpy)
        # 4.1 基础聚类 UMAP 与空间斑块映射
        sc.pl.umap(adata, color=["leiden"], wspace=0.4, show=False)
        plt.savefig("UMAP_leiden.pdf", bbox_inches='tight')
        plt.close()

        sq.pl.spatial_scatter(adata, color="leiden", wspace=0.4, library_id="spatial", show=False)
        plt.savefig("Spatial_leiden_scatter.pdf", bbox_inches='tight')
        plt.close()

        # 4.2 空间中心度与邻域富集分析
        sq.gr.spatial_neighbors(adata, coord_type="generic", delaunay=True)
        sq.gr.centrality_scores(adata, cluster_key="leiden")
        sq.pl.centrality_scores(adata, cluster_key="leiden", show=False)
        plt.savefig("centrality_scores.pdf", bbox_inches='tight')
        plt.close()

        sq.gr.spatial_neighbors(adata)
        sq.gr.nhood_enrichment(adata, cluster_key="leiden")
        sq.pl.nhood_enrichment(adata, cluster_key="leiden", show=False)
        plt.savefig("neighborhood_enrichment.pdf", bbox_inches='tight')
        plt.close()

        # 5. 肝脏特异性 Marker 基因气泡图 (Dotplot)
        genelist = {
            "Hepatocytes": ["CYP1A2", "CYP3A74", "ALB"],
            "LSECs": ["CLEC4G", "STAB1", "STAB2", "CLEC3B"],
            "Kupffer cells": ["C1QA", "C1QB", "BOLA-DRA"],
            "Hepatic Stellate": ["DCN", "AEBP1", "IGFBP3"],
            "Cholangiocytes": ["APOA1", "ALDOB"]
        }
        
        # 使用当前正确的 adata 对象，依据 leiden 聚类群绘制
        sc.pl.dotplot(adata, genelist, groupby='leiden', dendrogram=True, cmap='bwr', show=False)
        plt.gcf().set_size_inches(10, 5)
        plt.savefig('marker_dotplot.pdf', format='pdf', bbox_inches='tight')
        plt.close()

        # 6. 差异表达基因 (DEG) 分析
        sc.tl.rank_genes_groups(adata, groupby='leiden', use_raw=False, method="wilcoxon")
        
        sc.pl.rank_genes_groups_matrixplot(adata, groupby='leiden', values_to_plot='logfoldchanges', cmap='bwr', show=False)
        plt.savefig("deg_heatmap.pdf", bbox_inches='tight')
        plt.close()

        sc.pl.rank_genes_groups_dotplot(adata, groupby="leiden", standard_scale="var", n_genes=5, cmap='bwr', show=False)
        plt.gcf().set_size_inches(20, 12)
        plt.savefig('celltype_deg_dotplot.pdf', format='pdf', bbox_inches='tight')
        plt.close()

        # 7. 持久化保存
        adata.write("processed_data.h5ad")
        print("--> Pipeline finished successfully. Data saved to 'processed_data.h5ad'.")

    except Exception as e:
        print(f"Error occurred: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Spatial Transcriptomics Analysis Pipeline")
    parser.add_argument("--work_dir", type=str, required=True, help="Working directory path")
    parser.add_argument("--adata_path", type=str, required=True, help="Path to raw Visium data directory")
    args = parser.parse_args()
    
    main(args.work_dir, args.adata_path)