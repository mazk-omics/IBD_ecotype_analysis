# MSCCR Control：当前 85-state 模型的内部 nonIBD 覆盖审计

## 分析身份

本分析只回答“当前 85-state/8-ecotype 模型能否在 MSCCR Control 中恢复”的内部覆盖问题。它不把 MSCCR Control 宣称为独立外部验证，也不重新验证 discovery CD/UC 样本。

- 表达尺度：线性 CPM。
- CPM 分母：每个样本全部原始 Ensembl count 行的文库总量；模型基因筛选发生在 CPM 计算之后。
- 未重建 TPM。
- 未重跑 discovery CD/UC。
- MSCCR Control-only 数据不能单独计算 Cell 2021 已知表型富集 Z；所有样本表型相同，置换后表型计数不变。

## 输入重建和预检查

原有 Control CPM 只包含上一轮 12-type 模型的 4,553 个基因，对当前 20-type 模型的最低父类型覆盖率约为 27.9%，因此未进入 recovery。

`scripts/47_prepare_MSCCR_Control_CPM_current85.R` 从同一原始 counts 和冻结基因身份表重建当前模型基因并集输入：

- 461 个 Control、Non-inflamed 样本；243 位供者。
- Ileum 121、Colon 115、Rectum 225 个样本。
- 当前模型基因并集 12,960 个，成功映射 12,750 个。
- 20 个父类型最低模型基因覆盖率 96.6%。
- 三个部位均通过 `scripts/41_preflight_GSE57945_bulk_recovery.R` 的 75% 停止线。

## 官方 EcoTyper recovery 验收

`scripts/48_run_MSCCR_Control_bulk_recovery_current85.R` 按 Ileum、Colon、Rectum 顺序调用未修改的 `EcoTyper_recovery_bulk.R`，每个部位使用 20 threads。`scripts/49_finalize_MSCCR_Control_bulk_recovery_current85.R` 只读复核结果，没有重复运行 recovery。

- 60/60 个 `state_abundances.txt`。
- 60/60 个 `state_assignment.txt`。
- 3/3 个 `ecotype_abundance.txt`。
- 3/3 个 `ecotype_assignment.txt`。
- 三个部位状态均为 `ACCEPTED`。

首次包装启动失败保留在 `output/36_MSCCR_Control_bulk_recovery_current85`：外层系统 R 4.5 将 `R_HOME/R_LIBS_SITE` 传给 conda R 4.1，后者误载入 R 4.5 编译的 `data.table`。`retry1` 显式清除传给 conda 子进程的 R 环境变量后通过。该问题不属于输入数据或 EcoTyper 算法失败。

## Cell state 覆盖

| 部位 | 样本 | 供者 | 观察到的正式 state | 至少 5 个样本 | 分配比例至少 10% |
|---|---:|---:|---:|---:|---:|
| Ileum | 121 | 121 | 85/85 | 73 | 42 |
| Colon | 115 | 110 | 82/85 | 74 | 48 |
| Rectum | 225 | 225 | 84/85 | 78 | 49 |

未观察到的部位-state 组合只有 4 个：

- Colon：`CD8_T_S02`、`CD8_T_S04`、`Fibroblast_S05`。
- Rectum：`CD8_T_S04`。

因此 82/85 个 state 在三个部位均出现，84/85 个至少在两个部位出现；`CD8_T_S04` 只在 Ileum Control 中出现。

## Ecotype 覆盖

8 个 ecotype 在三个部位均被观察到。完整分母来自官方内部 `initial_ecotype_assignment.txt`；正式 `ecotype_assignment.txt` 会按官方代码过滤 `Ecotype == Unassigned` 的样本。

| 部位 | 完整样本 | 正式分配 | Unassigned | `Q < 0.25` | `Q < 0.05` | 观察到的 ecotype |
|---|---:|---:|---:|---:|---:|---:|
| Ileum | 121 | 121 | 0 | 121 | 108 | 8/8 |
| Colon | 115 | 114 | 1 | 114 | 104 | 8/8 |
| Rectum | 225 | 222 | 3 | 223 | 193 | 8/8 |

这里 `Q < 0.25` 是 EcoTyper `ecotypes_recovery.R` 的正式分配阈值；`Q < 0.05` 仅作为更严格的描述性指标。Rectum 中有 1 个样本满足 `Q < 0.25`，但未满足 `AssignedToEcotypeStates`，因此仍为 Unassigned。

## 与外部验证的整合

`scripts/51_integrate_nonIBD_bulk_scRNA_MSCCR_evidence.R` 将 GSE57945 外部 bulk、四套 scRNA nonIBD 验证与 MSCCR 内部覆盖并列，不把 MSCCR 覆盖计入显著性验证：

- 正式 state：85。
- GSE57945 nonIBD `Z > 1.96`：12。
- 至少一套 scRNA nonIBD `Z > 1.96`：31。
- 外部 bulk 与至少一套 scRNA 同时显著：3。
- MSCCR 三部位均观察到：82。
- MSCCR 至少两部位观察到：84。

严格跨模态同时支持的 3 个 state 为：

| State | GSE57945 bulk Z | 显著 scRNA 数据集 | MSCCR 覆盖 |
|---|---:|---|---|
| `CD8_T_S06` | 2.476 | SCP1884 | 3/3 部位 |
| `Fibroblast_S02` | 2.116 | GSE214695、IBDverse E-MTAB-16999、SCP1884、SCP259 | 3/3 部位 |
| `Ileal_absorptive_S02` | 2.630 | SCP1884 | 3/3 部位 |

## 结论边界

MSCCR Control 支持绝大多数 state 和全部 ecotype 在非 IBD、非炎症组织中的可恢复性，说明这些程序不是 IBD 样本专属。但 MSCCR Control 与 discovery 来源相同，且本步骤没有同流程的 MSCCR IBD recovery 对照，因此不能替代 GSE57945 的外部 bulk 表型富集，也不能把“观察到”解释为“nonIBD 显著富集”。

主要结果文件位于服务器：

- `output/35_MSCCR_Control_CPM_current85_preflight/reconstructed_current_model_gene_union`
- `output/36_MSCCR_Control_bulk_recovery_current85_retry1`
- `output/37_MSCCR_Control_bulk_recovery_summary`
- `output/38_integrated_nonIBD_validation_evidence`
