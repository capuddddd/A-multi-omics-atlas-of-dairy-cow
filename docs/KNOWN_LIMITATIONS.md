# 已知限制与需要验证的内容

本次整理完成的是代码可读性、明确输入输出及可确认缺陷的修复。
没有输入数据和 R/生物信息学运行环境，尚未进行端到端计算或结果一致性验证。

- SCENIC+ 示例使用旧版模块及构造方式，RNA/ATAC AnnData 不足以证明符合对应版本输入要求；
  motif 数据库路径已定义但未传入流程。需按官方教程准备 cisTopic、motif enrichment 等数据并重写对应入口。
  本次只移除导入时自动执行和全局屏蔽警告，并让失败产生非零退出码。
- `Findmotif.r` 已通过 Signac 的 motif.names 选择 footprint 目标，替换原 `HasMotif`/标量 `||` 筛选。
  不再强行修改染色体名；输入 sequence names 与 bosTau9 不符会报错。
  基因组组装、fragment 坐标与目标 motif 的大小写仍需真实数据确认。
- `Multi-Sample_scATAC-seq_Peak_Call.r` 仍保留原样本示例；
  删除 singlecell.csv 第一行、染色体只选 1–29/X、片段索引方法均需对真实文件确认。
  BED 到 GRanges 的坐标转换已明确为零基到一基，该修复可能改变原结果。
- DoubletFinder 元数据合并现在保留细胞索引并拒绝冲突；若有多轮不同 pANN 或分类，
  应明确选择目标列，不能静默拼接。污染过滤函数现在显式调用，且保留合并后的列。
- 多样本 H5 合并使用完整文件名作为 sample ID，避免下划线截断；
  若命名包含测序 lane，应由实际实验样本表定义合并关系。
- Visium 的线粒体列表和肝脏 marker panel、最小 counts、Leiden resolution 都是原研究参数；
  需要按组织、基因 ID 和样本重新确认。固定 20 PCs 不适合很小的数据集。
  保存 counts/raw 并改用 log expression 进行差异分析会改变旧版差异表达结果，应复核。
- 多个 Visium library 的空间邻居计算仍需确认 library_key 与跨样本连接策略。
- Monocle3 使用内部 slot，拟时序缺失值、无限值、零方差与 spline 自变量等情况仍需验证。
- 部分 R 作图依赖会话中的对象或预设元数据，README 已区分函数脚本和会话片段。
- 尚无成功运行环境的依赖锁、原始数据来源/访问号、参考文件校验和以及论文引用信息。

官方参考：

- [Tangram mapping API](https://github.com/broadinstitute/Tangram/blob/master/tangram/mapping_utils.py)
- [Tangram annotations API](https://github.com/broadinstitute/Tangram/blob/master/tangram/utils.py)
- [SCENIC+](https://github.com/aertslab/scenicplus)

- [Signac motif implementation](https://github.com/stuart-lab/signac/blob/master/R/motifs.R)
- [Signac footprint implementation](https://github.com/stuart-lab/signac/blob/master/R/footprinting.R)
