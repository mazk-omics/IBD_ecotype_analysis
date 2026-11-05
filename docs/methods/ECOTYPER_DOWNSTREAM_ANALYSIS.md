# EcoTyper 下游验证、解释与绘图流程

本文整理上一轮 EcoTyper 敏感性分析的下游代码、方法边界和图件验收入口。上一轮模型包含冻结的 53 个 cell states 与 7 个 ecotypes；上游表达来自 BayesPrism-derived G-like 矩阵，因此本目录中的结论和图件不能直接写成当前 CIBERSORTx HiRes 正式路线的结果。

## 阅读顺序

1. 先读方法一致性审计：`06_validation/paper_aligned_recovery_characterization_20261012/09_reports/METHODS_PARITY_AUDIT.md`。
2. 再读完整方法矩阵：`06_validation/paper_aligned_recovery_characterization_20261012/00_audit/method_parity_matrix.tsv`。
3. 按代码入口执行或审阅：`06_validation/paper_aligned_recovery_characterization_20261012/README.md`。
4. 图件方法、源数据和 QA：`06_validation/paper_figures_20261012/README.md`。
5. 新增或修改脚本遵循：`docs/methods/ECOTYPER_DOWNSTREAM_CODE_STYLE.md` 与 `Profile/CODE_STYLE.md`。

## 分析层级

| 层级 | 主要问题 | 核心脚本 | 结果解释边界 |
|---|---|---|---|
| 输入与模型身份 | 冻结状态、生态型、细胞及完整 counts 是否对应 | `preflight_and_full_counts.py` | 只核验身份，不重新 discovery 或 recovery |
| Recovery 审计 | 状态恢复、生态型共关联、疾病/炎症/肠段是否可复核 | `audit_recovery_and_phenotypes.py` | 状态 Z 不是细胞谱系准确率；重复活检不是独立患者 |
| 外部 marker | 同一父类型内的状态差异表达 | `external_markers_seurat313_parallel4.R` | Seurat 3.1.3 原生 Bonferroni 保留，不改称 BH |
| 七指标与 LOOCV | n1-n7 排名及留数据集评价 | `rank_markers_n4_and_LOOCV_parallel4_v3.py` | 只有两个外部队列；状态分配本身使用模型 marker |
| Marker 验收 | 训练/留出隔离、rank、维度与冻结模型完整性 | `check_final_marker_outputs.py` | 不把不可评价折补成 0 |
| 表型与肠段 | 状态/生态型的疾病、炎症和 Region3 表征 | `audit_recovery_and_phenotypes.py` | 明确统计单位、多重检验族和依赖结构 |
| Hallmark | 50 个 Hallmark 与 dominant ecotype 的描述 | `hallmark_bulk_features.py` | MSigDB 7.4；上游仍是 BayesPrism 敏感性路线 |
| 证据卡与报告 | 合并恢复、marker、表型与来源 | `assemble_characterization_reports.py` | BLOCKED/PARTIAL 不得自动升级为 PASS |
| 研究图件 | UMAP、ecotype 网络、特征总览 | `paper_figures_20261012/scripts/` | 网络是共现，不是通讯、空间邻近或因果 |
| 图件 QA | 字体、碰撞、对齐、源文件和导出格式 | `paper_figures_20261012/scripts/check_figures.py` | 排版 PASS 不替代统计或生物学验证 |

## 方法学验收状态

- `PASS`：本项目适配后的计算和独立核验均已完成，但仍需保留 `PAPER_ADAPTED_IBD` 标签。
- `PARTIAL`：已有结果可描述，但原文次数、独立性或正式统计尚不完整。
- `BLOCKED`：原文关键定义不足，不能用相似方法替代后宣称复现。上一轮全局 ecotype detection 与 ligand-receptor 模块属于此类。
- 所有上游表达相关结果都必须写明“BayesPrism sensitivity”，不能与正在运行的 CIBERSORTx HiRes 正式分析混写。

完整状态见 `method_parity_matrix.tsv`。该表同时记录公式、参数、实际输入、软件版本、适配和未解决问题。

## 推荐执行顺序

以下命令在项目根目录执行；历史脚本保留了服务器实际路径，正式复跑前先核对脚本开头的环境和路径。

```bash
python 06_validation/paper_aligned_recovery_characterization_20261012/scripts/preflight_and_full_counts.py
python 06_validation/paper_aligned_recovery_characterization_20261012/scripts/audit_recovery_and_phenotypes.py
Rscript 06_validation/paper_aligned_recovery_characterization_20261012/scripts/external_markers_seurat313_parallel4.R
python 06_validation/paper_aligned_recovery_characterization_20261012/scripts/rank_markers_n4_and_LOOCV_parallel4_v3.py
python 06_validation/paper_aligned_recovery_characterization_20261012/scripts/check_final_marker_outputs.py
python 06_validation/paper_aligned_recovery_characterization_20261012/scripts/assemble_characterization_reports.py
python 06_validation/paper_aligned_recovery_characterization_20261012/scripts/plot_verified_characterization.py
```

图件流程：

```bash
python 06_validation/paper_figures_20261012/scripts/prepare_paper_figures.py
Rscript 06_validation/paper_figures_20261012/scripts/compute_networks_paper_igraph122.R
python 06_validation/paper_figures_20261012/scripts/plot_missing_panels.py
python 06_validation/paper_figures_20261012/scripts/check_figures.py
```

## 验证要求

- 冻结模型、恢复输入和关键结果使用 SHA256，整理前后不得变化。
- 样本、供者、父类型、状态和生态型 ID 在合并前显式检查集合与顺序。
- 缺失状态、未分配样本和不可评价折保留为缺失或明确类别，不静默改成 0。
- P 值、Q 值、置换次数、随机种子、统计单位和多重检验族写入脚本或结果表。
- 每张新图同时保存 PDF、SVG、PNG 与逐图 source data；PDF 通过文字、碰撞和面板对齐检查。
- 图件 QA 只证明文件与版式满足约定，不证明分析结论正确。

代码归档本身可运行：

```bash
python scripts/python/audit_ecotyper_downstream_archive_20261105.py
```

该检查核对两套 GitHub LF 版本脚本清单的 SHA256、Python 语法和必要方法文档，不读取患者级结果或重跑统计。上一轮实际执行脚本的原始字节哈希仍单独保存在 `00_audit/scripts_manifest_final.tsv`，不能用仓库换行规范化后的哈希覆盖。

## 当前正式路线如何复用

当前 Signature-D TPM → CIBERSORTx HiRes → EcoTyper 运行完成后，可以复用本流程的结构，但不能直接复用上一轮的状态编号、生态型成员、marker 排名或图件源表。应先冻结新模型并重新生成：

1. model/state/ecotype identity manifest；
2. 新的 recovery 与表型审计输入；
3. 新模型对应的 marker、LOOCV 和多重检验族；
4. 新的逐图 source data 与 QA；
5. 标明 CIBERSORTx HiRes 正式上游，不再沿用 BayesPrism sensitivity 标签。

## Git 边界

Git 跟踪代码、方法说明、小型审计表、环境说明和 QA 摘要。原始数据、患者级矩阵、大型中间文件、完整运行环境、可重新生成图件和未发表主要结果继续留在研究存储中。
