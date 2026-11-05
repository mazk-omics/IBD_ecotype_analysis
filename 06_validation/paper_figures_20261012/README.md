# 上一轮 EcoTyper 下游图件与 QA 代码包

本目录归档冻结模型的研究图件代码。图件展示上一轮 BayesPrism 敏感性分析，不重新 NMF、不修改 state/ecotype，也不能用网络边推断通讯、空间邻近或因果关系。

## 默认入口

| 顺序 | 脚本 | 用途 |
|---:|---|---|
| 1 | `scripts/prepare_paper_figures.py` | 冻结输入、统计表与逐图 source data |
| 2 | `scripts/compute_networks_paper_igraph122.R` | R umap 0.2.0.0 与 igraph 1.2.2 数值坐标 |
| 3 | `scripts/reuse_existing_figures.py` | 逐字节复用已冻结官方图件 |
| 4 | `scripts/plot_missing_panels.py` | 绘制 Fig2B、Fig4B 与 Fig5A 缺失面板 |
| 5 | `scripts/plot_Fig4B_cell2021_style_20261017.R` | Fig4B 的 Cell 2021 风格补充版本 |
| 6 | `scripts/check_figures.py` | PDF 字体、碰撞、对齐和输入哈希 QA |

`export_official_colors.R`、`save_umap_config.R` 与 `install_paper_umap.R` 保存实际颜色和历史环境准备过程。`compute_umap_networks.R` 是早期组合版本；最终数值入口为 `compute_networks_paper_igraph122.R`。

`scripts/figure_qa/` 是当前 Python QA 工具；`scripts/nature_figure_qa/` 保留当时使用的完整副本及 R 对齐脚本。两套同名文件的哈希相同并非重复计算，而是历史运行目录布局。

## 图件契约

- 每张新图保存 PDF、SVG、PNG 和 source data。
- PDF/SVG 保留可编辑文字，字体下限 5 pt。
- Fig5A 的 385 个比较使用同一 BH family。
- 细胞比例色块的 0–1 是跨 ecotype 的显示缩放，不是实际百分比。
- UMAP 是冻结 H* 的样本级展示，不是单细胞 UMAP 或独立身份验证。
- Fig4B 是冻结 Jaccard 共现网络，0 边不画，但节点和非零边不得删减。

方法对照见 `reports/阅读入口与方法对照.md`，机器可读契约见 `reports/figure_contract.json`，最终 QA 汇总见 `reports/figure_QA_status.tsv` 与 `reports/final_plot_source_validation.json`。

`scripts_manifest_repository_20261105.tsv` 记录本次进入 Git 的全部图件脚本 SHA256。统一阅读入口：`docs/methods/ECOTYPER_DOWNSTREAM_ANALYSIS.md`。
