"""
SCENIC+ Single-Cell Multi-Omics Gene Regulatory Network Inference Pipeline
Environment Requirement: scenicplus (Python >= 3.8)
Author: Chen Ao
Date: 2026-06-28
"""

import os
import scanpy as sc
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

# 导入 SCENIC+ 核心核心组件
import scenicplus
from scenicplus.scenicplusclass import create_SCENICPLUS_object
from scenicplus.preprocessing.filtering import filter_genes, filter_regions
from scenicplus.wrappers.run_scenicplus import run_scenicplus

# ==========================================
# 1. 环境与参数初始化
# ==========================================
WORKDIR = "./scenicplus_output"
if not os.path.exists(WORKDIR):
    os.makedirs(WORKDIR)

os.chdir(WORKDIR)

# 设置并行线程数
N_CORES = 16 

# ==========================================
# 2. 读取多组学输入矩阵
# ==========================================
print("--> Loading processed scRNA-seq and scATAC-seq matrices...")
adata_rna = sc.read_h5ad("../rna_processed.h5ad")
adata_atac = sc.read_h5ad("../atac_processed.h5ad")

# 确保基因组坐标格式为符合标准的 chr-start-end
# 若原格式为 "chr1:1000-2000"，需将其转换为 "chr1-1000-2000"
adata_atac.var_names = adata_atac.var_names.str.replace(":", "-")

# ==========================================
# 3. 构建候选调控区域与空间关联
# ==========================================
# SCENIC+ 需要下载物种对应的基序数据 (Motif database)，此处以牛/脊椎动物基序数据库为例
# 请提前下载好对应的以下三个基础组件文件：
# 1. cisTarget.regions.feather (区域打分矩阵)
# 2. motifs-v10nr_clust-nr.vertebrates-m0.001-o0.0.tbl (Motif对齐注释表)

db_annot_fn = "../motifs-v10nr_clust-nr.vertebrates-m0.001-o0.0.tbl"
ctx_db_fn = "../bosTau9_screen_v10_nr.vertebrates.regions.feather" # 对应你的物种版本

# ==========================================
# 4. 初始化创建 SCENIC+ 核心对象
# ==========================================
print("--> Initializing SCENIC+ master container object...")

scplus_obj = create_SCENICPLUS_object(
    type_of_data={'RNA': 'AnnData', 'ATAC': 'AnnData'},
    data={'RNA': adata_rna, 'ATAC': adata_atac},
    multi_ome_mode=True,
    key_to_group_cells='cell_type' # 你的细胞类型或聚类群标志列名
)

# 过滤低质量基因和区域
filter_genes(scplus_obj, min_cells=5)
filter_regions(scplus_obj, min_cells=5)

print(scplus_obj)

# ==========================================
# 5. 核心：一步法运行 SCENIC+ 计算网络
# ==========================================
print("--> Running SCENIC+ GRN computation wrapper (This will take hours)...")

# 该函数将依次执行：
# 1. 寻找 TF 结合区域 (Cistarget Motif enrichment)
# 2. 估计 区域可及性 - 基因表达 之间的相关性 (Region-to-Gene models)
# 3. 估计 转录因子表达 - 区域可及性 之间的相关性 (TF-to-Region models)
# 4. 组装整合生成增强子调控子 (eRegulons)

try:
    run_scenicplus(
        scplus_obj=scplus_obj,
        variable_gene_keys=['highly_variable'],
        pfm_matrix_path="../jaspar2024.txt", # 使用你本地的 JASPAR2024 矩阵表
        species='bos_taurus',                # 设定对应目标物种
        assembly='bosTau9',                  # 基因组版本
        tf_file="../tf_names.txt",           # 物种特异性已知转录因子白名单列表
        save_path="scplus_obj_final.pkl",
        biomart_host="http://www.ensembl.org",
        upstream=150000,                     # 寻找增强子的上游最大搜索跨度
        downstream=150000,                   # 寻找增强子的下游最大搜索跨度
        calculate_TF_eRegulon_specificity=True,
        calculate_Gex_eRegulon_specificity=True,
        n_cores=N_CORES
    )
    print("--> SCENIC+ Pipeline executed successfully!")
except Exception as e:
    print(f"--> Pipeline failed with error: {str(e)}")