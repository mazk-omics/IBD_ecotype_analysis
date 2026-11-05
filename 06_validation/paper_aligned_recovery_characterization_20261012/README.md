# 上一轮 EcoTyper 下游表征代码包

本目录归档冻结 53 个 states / 7 个 ecotypes 的下游验证与表征。上游为 BayesPrism-derived G-like 敏感性路线；代码不适用于直接宣称当前 CIBERSORTx HiRes 正式路线的结果。

Git 只纳入脚本、方法一致性审计和小型验收记录。患者级矩阵、完整结果表、日志和大型中间文件继续保存在研究存储中。

## 默认入口

| 顺序 | 脚本 | 用途 |
|---:|---|---|
| 1 | `scripts/preflight_and_full_counts.py` | 冻结模型、细胞身份、原始 counts 与输入哈希预审 |
| 2 | `scripts/audit_recovery_and_phenotypes.py` | recovery、共关联、疾病、炎症和 Region3 核验 |
| 3 | `scripts/external_markers_seurat313_parallel4.R` | Seurat 3.1.3 全基因状态 DEG 与分层 DEG |
| 4 | `scripts/rank_markers_n4_and_LOOCV_parallel4_v3.py` | 最终两版 n4、七指标和两折 LOOCV |
| 5 | `scripts/check_final_marker_outputs.py` | 训练/留出隔离、维度、rank 与冻结模型验收 |
| 6 | `scripts/compare_n4_final.py` | Disease 与 Region3 两版 n4 的成对比较 |
| 7 | `scripts/hallmark_bulk_features.py` | 50 Hallmark 的 bulk ecotype 特征 |
| 8 | `scripts/assemble_characterization_reports.py` | 合并证据卡与方法边界 |
| 9 | `scripts/plot_verified_characterization.py` | 只绘制已验收的状态和生态型结果 |
| 10 | `scripts/plot_n4_comparison.py` | 两版 n4 的方法比较图 |

## 历史版本

- `external_markers_seurat313.R`、`rank_markers_n4_and_LOOCV.py`、`parallel4.py`、`parallel4_v2.py`：早期或串行版本，保留用于追溯。
- `pilot_B_rank_and_LOOCV.py`：只用于 B 细胞 smoke test。
- `finalize_marker_reports*.sh` 与 `run_marker_characterization.sh`：服务器实际调度记录；最终计算对应 `finalize_marker_reports_v4.sh`。
- `serve_analysis_transfer.py`：历史单连接传输辅助，不是分析入口。

`00_audit/scripts_manifest_final.tsv` 保存历史执行环境中的原始字节 SHA256。GitHub 归档仅机械统一了换行符为 LF，未改变算法、参数或命令；统一换行后的仓库哈希见 `scripts_manifest_repository_20261105.tsv`，并由 `scripts/python/audit_ecotyper_downstream_archive_20261105.py` 复核。

## 方法边界

- `PASS` 表示 IBD 适配后的实现和核验通过，不表示逐字复现原论文全部代码。
- recovery Z、共关联和 marker LOOCV 分别回答不同问题，不能合成“验证准确率”。
- 30 次 shuffle 的既有状态 Z 不写成原论文 100 次等价实现。
- 全局 ecotype detection 与 ligand-receptor 原版定义不足，保持 `BLOCKED`。
- 详细公式、参数、版本和 unresolved 项见 `00_audit/method_parity_matrix.tsv`。

统一阅读入口：`docs/methods/ECOTYPER_DOWNSTREAM_ANALYSIS.md`。
