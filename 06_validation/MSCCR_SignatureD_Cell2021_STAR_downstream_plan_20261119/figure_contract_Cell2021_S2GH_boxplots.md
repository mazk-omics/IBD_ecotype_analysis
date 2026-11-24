# Cell 2021 Fig. S2G/S2H boxplot contract

- 结论：用 discovery bulk 队列比较不同亲本细胞类型和不同 IBD 表型情境中的 cell-state representation specificity。
- 主数据：discovery bulk 的样本级 EcoTyper state assignment；每个样本在每个亲本细胞类型内已经是一个 dominant state。四套外部 scRNA 数据只用于 validation 扩展，不替代本主图。
- 检验：对每个 state × 表型水平构建 2 × 2 表，执行单侧 Fisher 精确检验；在每个表型水平内跨全部 cell states 做 BH 校正。
- 展示值：保留全部可评估组合，包括非富集和不显著组合；使用 -log10(BH q-value)，上限 10。不得先按 OddsRatio 或显著性筛选后再求平均。
- 主表型轴：Disease 与 Region3 分开计算和作图，禁止机械交叉成 Disease × Region3；Inflammation 作为次级图。
- G：每个点为一个 cell state；先在当前单一表型轴的所有水平中平均 enrichment，再按亲本细胞类型画箱线图。
- H：每个点为一个亲本细胞类型；先在该类型的 states 中平均 enrichment，再按当前 discovery 表型水平画箱线图。
- 排序：G 先按 lineage 连续分组，再在 lineage 内按 enrichment 中位数递减；H 按组内 enrichment 中位数递减。
- 版式：复刻 Cell 2021 Fig. S2G/S2H。G、H 横向并列；灰色箱线图；彩色散点；G 底部 lineage 色带及双图例；H 右侧 cell-type 图例。
- 输出：每个表型轴分别输出 PDF、SVG、600 dpi TIFF、300 dpi PNG 及 source-data TSV。
- 限制：本项目有 20 个亲本细胞类型，横轴比原文更密；保持原文视觉语法但增加画布宽度。
- 解释边界：由于 discovery 中 Disease 与 Region3 的构成不完全平衡，两个边际 Fisher 结果分别表示疾病或部位关联，不能解释为互相校正后的独立效应。
