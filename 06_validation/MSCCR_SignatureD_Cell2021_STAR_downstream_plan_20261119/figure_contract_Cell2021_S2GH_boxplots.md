# Cell 2021 Fig. S2G/S2H boxplot contract

- 结论：比较不同亲本细胞类型和不同数据集表型情境中的 cell-state representation specificity。
- 数据：每个样本在每个亲本细胞类型内的 dominant state；单侧 Fisher 精确检验；同一数据集和表型轴内跨细胞类型 BH 校正。
- 正向富集：仅保留 OddsRatio > 1 且可评估的组合；展示值为 -log10(BH q-value)，上限 10。
- G：每个点为一个 cell state；先在数据集×表型水平中平均正向 enrichment，再按亲本细胞类型画箱线图。
- H：每个点为一个亲本细胞类型；先在该类型的 states 中平均正向 enrichment，再按数据集×表型水平画箱线图。
- 排序：G 和 H 均按组内 enrichment 中位数递减。
- 版式：复刻 Cell 2021 Fig. S2G/S2H。G、H 横向并列；灰色箱线图；彩色散点；G 底部 lineage 色带及双图例；H 右侧 cell-type 图例。
- 输出：每个表型轴分别输出 PDF、SVG、600 dpi TIFF、300 dpi PNG 及 source-data TSV。
- 限制：本项目有 20 个亲本细胞类型，横轴比原文更密；保持原文视觉语法但增加画布宽度。
