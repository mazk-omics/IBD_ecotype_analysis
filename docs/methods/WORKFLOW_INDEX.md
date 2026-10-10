# 研究流程与代码阅读入口

本索引用于 2026-10-10 代码归档后的逐步审阅。归档没有重新计算数据或统一历史脚本的参数；实际输入矩阵、临床注释、冻结模型和输出保留在研究存储中。运行前核对每个脚本开头的输入路径、环境、阶段与输出目录。

| 阶段 | 按执行顺序阅读的入口 | 检查重点 |
|---|---|---|
| IBDverse 原始结构与元数据 | `scripts/python/inspect_ibdverse_h5ad.py` → `audit_ibdverse_reference_inputs.py` | counts/归一化层、样本与供者连接、冲突隔离；服务器版与本地历史差异见来源清单 |
| 参考导出 | `scripts/python/three_reference_pilots_v1/`、`joint20_v2/build_joint20_donor_site_reference_v2.py` | 细胞类型映射、donor/site/sample 分组、训练与留出供者边界、CPM 尺度 |
| 留出验证与网页输入 | `scripts/python/generate_joint20_holdout_mixtures.py`、`prepare_job13_equalrna8_validation.py`、`prepare_job13_natural28_validation.py` | 细胞数比例与 RNA 比例的区别、基因顺序、独立留出边界 |
| MSCCR 小规模输入与 S-mode | `scripts/python/prepare_msccr_platform120_pilot.py` → `prepare_msccr_joint20v2_web.py` → `prepare_gj20v2_smode_source_subset.py`；`review_cibersortx_job17_msccr_smode.py` 与 `review_cibersortx_job20_msccr_repeat.py` | discovery 样本范围、样本一致性、源参考对应关系、重复运行差异 |
| HiRes 与接口审计 | `scripts/python/prepare_msccr120_hires_upload.py` → `prepare_msccr120_hires_gene_batches.py` → `review_cibersortx_job21_hires.py` / `job23` / `job24` → `check_cibersortx_hires_for_ecotyper.py` | 网页基因限制、块间重叠、类型表达是否可变；不能仅凭文件可合并判定表达可靠 |
| BayesPrism 敏感性 discovery | `scripts/R/ecotyper/01_prepare_core12_sensitivity_inputs.R`、`scripts/bash/ecotyper/01_launch_core12_BC_parallel.sh`、`05_EcoTyper/config/MSCCR_CORE12/` | Z/G-like 转换、presorted 类型与样本顺序、Region3 实际分组标准化；官方源码补丁见 `docs/software_versions/` |
| Bulk recovery | `scripts/recovery/gse72819_bp_20261009/`、`06_validation/01_bulk_progression/` 与 `06_validation/02_bulk_treatment/` 中各 `scripts/` | 数据来源、基因身份、冻结模型、官方 recovery 输入规范；不要把敏感性恢复等同独立验证 |
| 单细胞 recovery | `scripts/bash/run_ibdverse_recovery_pilot_20261009.sh`、`run_gse214695_recovery_pilot_20261009.sh`，配套 `prepare_*` 与 `summarize_*` | 部位/谱系覆盖、输入尺度、细胞注释、冻结模型及运行完成检查 |
| 初步状态解释 | `scripts/python/prepare_paper_annotation_inputs_20261010.py` → `scripts/R/interpretation_20261010/interpret_states_ecotypes_paper_20261010.R` → `check_paper_annotation_outputs_20261010.py` | 论文方法的实际适配、数据库版本、统计阈值、候选命名证据与限制 |
| nonIBD/健康敏感性恢复 | `scripts/python/prepare_nonibd_recovery_20261010.py` → `scripts/bash/run_nonibd_recovery_20261010.sh` → `summarize_nonibd_recovery_20261010.py` 与 `scripts/R/nonibd_coassociation_20261010.R` | 全量上游是否包含对照，外部队列与同源对照的证据区分 |
| 图件与交付检查 | `scripts/python/plot_*recovery_nature.py`、`assemble_*recovery_figures.py`、`present_*20261010.py`、`check_recovery_delivery_20261010.py` | 数据源、图件一致性与输出检查；最终研究图件仍保存在研究目录 |

## 环境与版本

- 解释环境的 Conda 导出、显式包列表和 R sessionInfo：`environments/conda/` 与 `environments/R/`。
- 官方 EcoTyper 基线、项目 HEAD 和实际工作树补丁：`docs/software_versions/ecotyper_server_20261010.json`。
- 安装/修复环境的原脚本保留在 `scripts/archived/20261010_research_folder_sources/06_validation/ecotyper_original_method_annotation_20261010/`。这是历史执行步骤，不表示需要重复执行全部修复。
- 当同名源码不同，服务器当前实际文件作为运行来源；本地历史副本放在 `scripts/archived/20261010_local_pre_server_sync/`，差异映射见 `docs/provenance/local_server_reconciliation_20261010.tsv`。

## 尚未由此次归档验证的内容

语法和归档检查不替代分析正确性与生物学验证。当前 CIBERSORTx 正式 HiRes 路线、独立 recovery 证据及后续临床关联仍需按研究计划推进；本轮 Git 收尾没有重新运行这些计算，也没有修改输入或统计参数。
