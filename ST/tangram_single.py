"""
Tangram Spatial Cell Type Mapping & Deconvolution Pipeline
Environment Requirement: tangram-sc, scanpy, pytorch
Author: Chen Ao
Date: 2026-06-28
"""

import argparse
import os
import scanpy as sc
import tangram as tg
import pandas as pd
import numpy as np
import torch
from matplotlib import pyplot as plt

def main(work_dir, sc_path, st_path, annotation_key):
    # 设置工作目录
    os.chdir(work_dir)
    
    # 自动检测 GPU 加速
    device = "cuda:0" if torch.cuda.is_available() else "cpu"
    print(f"--> Using device: {device}")

    try:
        # 1. 载入单细胞与空间数据
        print("--> Loading reference scRNA-seq and target Spatial data...")
        adata_sc = sc.read_h5ad(sc_path)  # 单细胞源数据
        adata_st = sc.read_h5ad(st_path)  # 空间组学靶数据（上一步处理好的 processed_data.h5ad）

        # 确保数据中不包含非稀疏的大型非必需层，保持 counts 干净
        if 'counts' in adata_sc.layers:
            adata_sc.X = adata_sc.layers['counts'].copy()
        if 'counts' in adata_st.layers:
            adata_st.X = adata_st.layers['counts'].copy()

        # 2. 筛选共享的训练特征基因 (Marker Genes)
        print("--> Selecting training marker genes for alignment...")
        # 挑选 scRNA-seq 中的高变基因作为特征引导
        if 'highly_variable' not in adata_sc.var:
            sc.pp.highly_variable_genes(adata_sc, flavor='seurat', n_top_genes=2000, subset=False)
        
        # 提取两者的交集基因并进行切片对齐 (Tangram 核心预处理)
        tg.pp_adatas(adata_sc, adata_st, genes=None, lookup_genes=None)

        # 3. 运行 Tangram 空间映射模型
        print(f"--> Training Tangram mapping model based on cluster annotation: '{annotation_key}'...")
        # clusters_map_all 模式会将 scRNA 中的每个单细胞单独映射，随后按 annotation_key 聚合到空转 Spot
        cv_adata = tg.map_cells_to_space(
            adata_cells=adata_sc,
            adata_space=adata_st,
            mode="clusters_map_all",
            cluster_label=annotation_key,
            device=device
        )

        # 4. 提取反褶积预测丰度矩阵并回写至空转对象
        print("--> Extracting and projecting cell-type spatial density matrices...")
        # 从模型结果中投影每种细胞类型的空间概率密度
        tg.project_cell_types(cv_adata, adata_st, section_name="spatial")

        # 5. 可视化：保存主要细胞类型的空间丰度图
        print("--> Plotting cell-type spatial distribution tracks...")
        # 获取推断出的所有细胞类型名称
        cell_types = list(cv_adata.obs[annotation_key].unique())
        print(f"Detected cell types mapping: {cell_types}")

        # 使用 Scanpy 绘制这些推断细胞类型在空间 Spot 上的连续丰度权重
        # 预测出的细胞类型丰度会以当前名字自动回写进 adata_st.obs 中
        sc.pl.spatial(
            adata_st, 
            color=cell_types[:6],  # 默认展示前6个主要细胞类型，可根据需求调整
            cmap="viridis", 
            library_id="spatial", 
            show=False
        )
        plt.gcf().set_size_inches(15, 10)
        plt.savefig("tangram_celltype_spatial_density.pdf", bbox_inches='tight')
        plt.close()

        # 6. 数据持久化保存
        # 将带有 Tangram 反褶积丰度结果（存储在 obs 中）的空间对象重新写入新文件
        output_st_path = "spatial_with_tangram.h5ad"
        adata_st.write(output_st_path)
        print(f"--> Pipeline successfully completed. Results saved to '{output_st_path}'.")

    except Exception as e:
        print(f"Error occurred during Tangram alignment: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Tangram Spatial Deconvolution Pipeline")
    parser.add_argument("--work_dir", type=str, required=True, help="Working directory path")
    parser.add_argument("--sc_path", type=str, required=True, help="Path to reference scRNA-seq h5ad file")
    parser.add_argument("--st_path", type=str, required=True, help="Path to target Spatial h5ad file")
    parser.add_argument("--annotation_key", type=str, default="celltype", help="The obs column name in scRNA-seq representing cell types")
    
    args = parser.parse_args()
    main(args.work_dir, args.sc_path, args.st_path, args.annotation_key)