# Bovine multi-omics analysis scripts

奶牛组织单核 RNA、单核 ATAC 和 Visium 空间转录组分析代码。
当前仓库是研究脚本集合：部分文件提供函数，部分依赖已准备好的 Seurat 对象。
没有随仓库提供原始数据、完整参考文件或验证环境，因此不能据此宣称整套流程可直接复现论文结果。

## 代码导航

| 文件 | 用途 | 运行方式 |
| --- | --- | --- |
| `snRNA/Merge_Single-Cell_h5_Datasets.r` | 合并 10x H5 数据 | Rscript 命令行 |
| `snRNA/QC_and_Basic_Analysis.r` | DoubletFinder、去污染、降维、富集 | source 后调用函数 |
| `snRNA/Celltype_Identification_and_Visualization.r` | 细胞注释与作图 | source 后调用函数 |
| `snRNA/CellChat_Inter.r` | 细胞通讯比较 | source 后调用函数 |
| `snRNA/Monocle3_pipeline.r` | 轨迹与拟时序作图 | source 后调用函数 |
| `snATAC/Multi-Sample_scATAC-seq_Peak_Call.r` | 共识 peaks 和片段计数 | 配置输入后运行 |
| `snATAC/QC_basic_analysis.r` | ATAC QC、整合、基因活性 | source 后调用函数 |
| `snATAC/integration.r` | RNA 到 ATAC 标签转移 | 配置环境变量后运行 |
| `snATAC/Findmotif.r` | motif、chromVAR、footprint | 依赖会话中的 `sc` 和 bosTau9 匹配检查 |
| `snATAC/heatmap_RNA&ATAC.r` | TF 活性与表达热图 | 依赖会话中的 `sc` |
| `ST/basic_analysis.py` | Visium QC、聚类、空间统计 | Python 命令行 |
| `ST/tangram_single.py` | 单细胞到空间映射 | Python 命令行 |
| `snATAC/scenicplus.py` | SCENIC+ 旧接口示例 | 尚未验证，见限制说明 |

## 环境

空间分析依赖 `scanpy`、`squidpy`、`matplotlib`、`igraph`、`leidenalg`；
映射额外依赖 `tangram-sc` 和 `torch`。SCENIC+ 应使用单独环境，并按其对应版本的官方教程准备输入。
R 依赖列表见 [DEPENDENCIES.md](docs/DEPENDENCIES.md)。
尚无实际运行环境的版本锁定文件；请在成功运行后保存 `pip freeze` 和 R 的 `sessionInfo()`，
不要把未经测试的版本列表当作可复现环境。

## 运行示例

从仓库根目录运行；输入路径相对于启动命令时的目录解析，`--work_dir` 是输出目录。
输出使用固定文件名，重复运行请使用新的目录，以保留之前结果。

```bash
Rscript snRNA/Merge_Single-Cell_h5_Datasets.r data/h5 results/merged.rds

python3 ST/basic_analysis.py \
  --adata_path data/visium/sample1 --work_dir results/spatial/sample1

python3 ST/tangram_single.py \
  --sc_path data/reference.h5ad \
  --st_path results/spatial/sample1/processed_data.h5ad \
  --annotation_key celltype --work_dir results/tangram/sample1 --seed 0
```

一个 Visium library 时自动推断 ID；多个时使用 `--library_id` 明确选择。
空间预处理保存 `layers["counts"]` 和 log-normalized `raw`；差异表达使用 `raw`。
Tangram 优先读取 counts 层，保存映射矩阵、CSV 和 `obsm["tangram_ct_pred"]`。
这些预测是映射权重，不应直接解释为绝对细胞数。

ATAC 输入通过环境变量配置：

```bash
export ATLAS_ATAC_DATA_DIR="$PWD/data/atac"
export ATLAS_RNA_REFERENCE="$PWD/data/reference.rds"
export ATLAS_ATAC_QUERY="$PWD/data/query.rds"
```

Peak 脚本中的 `DATA_PATHS`、`SAMPLE_IDS` 仍需按实际样本配置。
基因组组装、染色体命名、注释与片段文件必须一致。

## 验证和贡献

```bash
python3 -m unittest discover -s tests -v
python3 -m compileall -q ST snATAC
```

单元测试验证文件路径与 library 选择，不代表科学分析已验证。
修改 QC、物种注释、归一化、差异分析或模型参数时，必须用代表性数据比较前后结果。
更多未解决的问题见 [KNOWN_LIMITATIONS.md](docs/KNOWN_LIMITATIONS.md)。
许可证和数据来源状态见 [DATA_AND_LICENSE.md](docs/DATA_AND_LICENSE.md)。
