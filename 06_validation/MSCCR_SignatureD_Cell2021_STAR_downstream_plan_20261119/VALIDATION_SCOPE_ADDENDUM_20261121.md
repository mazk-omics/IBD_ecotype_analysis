# Signature-D validation 补充范围

本补充文件冻结 2026-11-21 新增的验证范围。Cell 2021 STAR Methods 和官方 EcoTyper 实现是方法学主依据；上一轮 `paper_aligned_recovery_characterization_20261012` 只作为代码、输出结构和数值回归基线，旧模型结果不并入本轮推断。

## 1. 上一轮可复用内容

- Seurat v3.1.3 数值 smoke test、`FindMarkers(test.use="wilcox", min.pct=0.1, logfc.threshold=0.05)` 和默认全基因 Bonferroni `p_val_adj`；
- `<5` assigned cells 的 state 排除规则、逐数据集完整表达读取和不可评价分母表；
- n1–n7、基因外部归属与 discovery 归属一致性、七指标等权平均 rank、每 state top 10 markers；
- 以完整 dataset 为留出单位的 LOOCV，held-out dataset 不参与 marker 选择；
- DiseaseGroup 与 Region3 两版 n4、表型 1,000 次 assignment-label permutation、UMAP 和报告结构；
- 旧结果只用于程序回归测试。本轮 state、cell、dataset 或模型身份不同时，不要求数值复现旧结果。

## 2. 本轮单细胞主验证

四个正式 scRNA-seq 数据集分别完成 state recovery、100 次 recovery shuffle、全转录组 marker 检验、n1–n7 排名和 dataset-level LOOCV。所有步骤使用本轮冻结的 Signature-D state、细胞清单、父类型映射和 gene identity；不根据 recovery 或 marker 结果重新筛 state、marker、细胞或数据集。

### 2.1 Known phenotype enrichment

按 Cell 2021 原法，对预先冻结的 phenotype × state 统计分配细胞或样本数，打乱 state assignment labels 1,000 次，计算 Z，并以 `Z > 1.96` 报告原文轨道。Disease、Inflammation、Region3 和 nonIBD/IBD contrasts 在运行前冻结；患者级受限置换或混合模型作为独立 IBD extension，不能替代原文轨道。

### 2.2 Cell-state representation

参照 Cell 2021 Figure S2G–S2H 的 tumor-type representation，在本项目中分别输出 Disease、Inflammation、Region3、dataset 和 IBD/nonIBD 的 state representation。保留每个 parent type 的输入、正式分配和 Unassigned 分母，不以 pooled 比例掩盖缺失 strata。该部分保持原图统计结构，但 tumor type 到 IBD phenotype 的映射标记为 `IBD_EXTENSION`。

### 2.3 UMAP

主图按 Cell 2021 方法使用正式 state abundance `H*` 的 Euclidean distance 和 R `umap` 0.2.0.0 默认参数，逐 parent type 显示并冻结随机种子。基于全转录组表达的细胞 UMAP 仅作为补充可视化，不能替代 H* UMAP，也不能用于重新划分 state。

### 2.4 State annotation

先完成不使用 phenotype 的 marker 排名，再进行注释。每个 state 保存 discovery markers、四队列 external markers、方向一致性、已知功能 marker、文献证据和数据库版本。手工或文献 marker 不进入 LOOCV 准确率；相互矛盾的证据保留，不强制命名。注释结果与 marker validation 分表输出，避免循环论证。

## 3. nonIBD 验证

nonIBD scRNA-seq 和 bulk RNA-seq 均为正式验证层。运行前分别冻结队列、组织部位、样本/供者、健康定义、疾病排除、表达单位和可评价父类型。

- nonIBD scRNA-seq：投影到同一冻结 state，不重新 NMF；完成 recovery、state representation、known phenotype enrichment、marker 支持和可评价分母审计；
- nonIBD bulk：使用同一 Signature-D/EcoTyper 模型恢复 state/ecotype abundance，报告 state 覆盖、正式 assignment、组织部位和 cohort 差异；
- IBD 与 nonIBD 的比较按 cohort 先分析再汇总，不能把平台或组织部位差异当作疾病效应；
- nonIBD 阴性结果保留，不因恢复率低而修改模型或删除 state。

## 4. bulk RNA validation

bulk validation 同时包含 IBD 和 nonIBD 队列。每个队列独立完成输入尺度审计、模型 recovery、state/ecotype coverage、已知 phenotype enrichment 和方向一致性，再做预先指定的跨队列汇总。发现队列不回流为外部验证；用于构建 Signature-D 的样本不得计入独立 validation 分母。

## 5. 放行顺序

1. G5：四数据集 100-shuffle recovery 完成；
2. G6：四数据集 n1–n7 与完整 dataset-level LOOCV；
3. G6A：H* UMAP、marker/literature annotation 和不可评价分母审计；
4. G7：Cell 2021 exact known phenotype enrichment；
5. G7A：S2G–S2H 对应的 Disease/Region/IBD-nonIBD representation；
6. G11A：IBD 与 nonIBD bulk validation；
7. G11B：独立 nonIBD scRNA validation；
8. 所有适用 gate 通过后才可将 validation package 标记为 COMPLETE。
