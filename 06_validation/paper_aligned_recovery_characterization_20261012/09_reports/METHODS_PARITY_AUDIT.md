# 方法对照审计

本轮针对冻结53状态/7生态型，只读原件，不重新discovery或修改recovery assignment。所有涉及上游表达的结果均为BayesPrism敏感性分析，不宣称CIBERSORTx数值等价。

两版七指标和两个留出折已计算，具体差异见 n4_disease_vs_region_comparison.tsv 与 marker_LOOCV_summary.tsv。

| Module | Validation | Parity class | Unresolved |
|---|---|---|---|
| Model identity | PASS | PAPER_ADAPTED_IBD | Reference donor overlap not globally provable from normalized cohort IDs |
| Recovery evidence | PARTIAL | PAPER_ADAPTED_IBD | Official P absent; Z evidence is not cell-lineage accuracy |
| External DEG | PASS | PAPER_ADAPTED_IBD | Native Bonferroni p_val_adj kept; not relabeled as BH |
| Seven metrics | PASS | PAPER_ADAPTED_IBD | Original rank/tie/BH-family code unavailable; two datasets restrict variability |
| Marker LOOCV | PASS | PAPER_ADAPTED_IBD | State assignment itself uses model markers; not independent identity truth |
| State phenotype | PASS | PAPER_ADAPTED_IBD | Cell dependence; sample-dominant unit is extra adaptation; HC surgical controls |
| Region representation | PASS | PAPER_ADAPTED_IBD | Repeated patients; undefined original exact family |
| Co-association | PASS | PAPER_ADAPTED_IBD | Original P tail not explicit; retain Z, no significance stars |
| Global detection | BLOCKED | BLOCKED | Original math missing in accessible BioC; joint criterion / membership null not confirmed |
| Discrete ecotype features | PASS | PAPER_ADAPTED_IBD | Discovery-associated; repeated biopsies not independent patient outcomes |
| Hallmark ecotype features | PASS | PAPER_ADAPTED_IBD | Paper MSigDB release not confirmed; not original CIBERSORTx upstream |
| State biological cards | PARTIAL | PAPER_ADAPTED_IBD | Mast/Glial cross-lineage concerns preserved; overlap is not formal GSEA |
| Ligand-receptor | BLOCKED | BLOCKED | Cannot yet claim original-equivalent LR; FigureG not fabricated |
| Existing bulk/nonIBD | PARTIAL | PAPER_ADAPTED_IBD | Full new metadata join / all-run manifest not completed this batch |


主要依据：[Luca et al., Cell 2021 STAR Methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC8526411/)。全局CE与LR阻断只影响各自模块，没有用替代检验冒充原文。服务器环境规格、函数默认源码、数值烟雾测试及输入校验保存在00_audit。
