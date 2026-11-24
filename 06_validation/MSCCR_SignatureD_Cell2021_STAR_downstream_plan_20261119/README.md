# Signature-D：Cell 2021 STAR Methods 对齐的下游分析方案

本目录是本轮 `Signature-D TPM -> CIBERSORTx HiRes -> EcoTyper` 正式模型的下游分析预注册方案。它替代旧包中“先做描述性关联、再把外部 recovery 放到可选阶段”的顺序，但当前只冻结方法，不运行服务器任务，也不改动已验收模型。

## 核心原则

1. `CELL2021_EXACT`：能够从 Cell 2021 STAR Methods 原样执行的步骤，公式、阈值、置换次数、统计单位和软件版本均固定。
2. `OFFICIAL_SOFTWARE`：EcoTyper 官方代码的实际实现。若与论文文字不同，必须并列复算并记录差异，不能静默选择。
3. `IBD_EXTENSION`：重复活检、Disease、Inflammation、Region3、跨疾病和跨部位保守性等 IBD 专属问题，独立于原文复现主链。
4. `SUPPLEMENTARY_SENSITIVITY`：BayesPrism、细胞类型删减、替代多重校正等只作为敏感性分析。
5. `NOT_APPLICABLE`：肿瘤生存、ICI、癌症突变特征等没有同构 IBD 终点时明确不执行，不用 OR 或横断面比较冒充原文 Cox/ICI 分析。

## 文件

- `CELL2021_STAR_METHODS_DOWNSTREAM_PLAN.md`：完整分析顺序、公式、参数、统计单位、输出与停止规则。
- `CELL2021_STAR_METHODS_CROSSWALK.tsv`：STAR Methods 模块与本项目实现的一一对应表。
- `VALIDATION_SCOPE_ADDENDUM_20261121.md`：冻结 marker LOOCV、known phenotype enrichment、S2G–S2H 对应 representation、H* UMAP、state annotation，以及 IBD/nonIBD 单细胞和 bulk validation 的新增正式范围。
- `scripts/16_prepare_state_annotation_literature_database_review.R`：在不修改冻结模型和首轮标签的前提下，为 85 个 state 建立逐条文献、数据库、组织背景和冲突证据审阅表。
- `inputs/state_annotation_source_register.tsv`：登记正式使用、仅用于解释和待接入的注释来源；版本、DOI/登录号、用途与限制必须显式保存。
- `scripts/17_cross_reference_state_annotation_databases.R`：将冻结 marker 与 CellMarker 2.0、PanglaoDB 和 HPA v26 交叉核对，分别保存肠道优先命中、广义命中和 HPA 表达排序，不自动改名。
- `inputs/state_annotation_database_celltype_map.tsv`：显式记录 20 个母细胞类型在三个数据库中的对应规则和解释限制。
- `STATE_ANNOTATION_EVIDENCE_METHOD_20261123.md`：冻结 state 后的正式注释证据层、命名边界、A-D 文献证据分级和文章表述规则。
- `scripts/18_integrate_high_confidence_literature_review.R`：将 26 个高注释置信度 state 的逐条原始文献证据、冲突和保留/改名决定写入 v1 审阅表；其余 state 保持待审。
- `inputs/state_annotation_literature_evidence_high_confidence.tsv`：首批逐 state 文献证据记录，包含来源 ID、参考 state、组织/疾病背景、证据等级和命名边界。
- `scripts/19_integrate_remaining_literature_review.R`：审阅其余 59 个 state，输出无待审项的完整 v2 表和按谱系证据等级审计。
- `inputs/state_annotation_literature_evidence_remaining.tsv`：其余中低置信 state 的正向、缺失与冲突证据；允许正式结论为未解析或跨谱系，而非强行命名。
- `output/10_state_annotation/literature_database_review/state_annotation_manuscript_label_table_v2.tsv`：面向文章的精简标签表，按 A-B、C、D 分为正文可解释、描述性程序和补充材料未解析/冲突三层。
- `scripts/20_align_annotation_interpretability_with_cell2021.R`：按 Cell 2021 的实际注释口径，将 A-B 级定义为文献锚定命名状态、C 级定义为 marker 驱动功能程序，并将 A-C 共同计入可解释状态。
- `output/10_state_annotation/literature_database_review/state_annotation_manuscript_label_table_v3_cell2021_aligned.tsv`：修正后的文章标签表；55/85 个 state 可解释，其中 25 个文献锚定、30 个功能程序。
- `scripts/21_cross_reference_published_state_signatures.R`：复用上一轮冻结的 published-signature 资源，对当前 85 个 state 的可用 marker 做同谱系重叠审计；只生成升级候选，不自动改变证据等级。
- `scripts/22_run_known_phenotype_enrichment_exact.R`：复现 G7 的逐细胞置换算法；当前输入只含 IBD，故正式身份为 `CELL2021_STATISTICAL_FRAMEWORK_IBD_WITHIN_DISEASE_EXTENSION`，待 nonIBD recovery 后才形成原文生物学对照轨道。
- `scripts/23_run_known_phenotype_enrichment_patient_aware.R`：G7 IBD 扩展；Disease 在供者层置换，Inflammation/Region3 在供者内置换样本标签，不替代原文逐细胞轨道。
- `scripts/24_run_cell_state_representation.R`：G7A；按样本 dominant state 和 Fisher 检验输出 Disease、Inflammation、Region3、dataset representation，并将当前 IBD/nonIBD 层记为待 G11B。
- `scripts/25_plot_marker_validation_and_z.R`：用 R 输出 marker LODO 热图、四队列 recovery Z、Stouffer Meta Z，以及 exact/患者感知 phenotype Z 图和逐图 source data。
- `scripts/26_plot_cell_state_representation.R`：用 R 输出 pooled 与逐数据集的 Disease、Inflammation、Region3、dataset dominant-state Fisher representation 热图；无 nonIBD 时不伪造对照图。
- `scripts/27_summarize_G7_G7A_validation.R`：汇总 recovery、marker LODO、exact/患者感知 phenotype enrichment 和 dominant-state representation，输出分母清楚的中文审阅报告。
- `scripts/52_run_discovery_state_representation_Cell2021.R`：严格使用 discovery bulk 样本级 state assignment，按 Cell 2021 的单侧 Fisher 和表型水平内跨全部 state 的 BH 校正计算 representation specificity；Disease 与 Region3 独立计算，禁止机械交叉。
- `scripts/53_plot_discovery_Cell2021_S2GH_boxplots.R`：从 discovery representation 结果分别输出 Disease 和 Region3 两套 Cell 2021 Fig. S2G/S2H 主箱线图，Inflammation 为次级图；四套 scRNA 图不替代此主图。
- `scripts/28_audit_and_freeze_nonIBD_scRNA_cells.R`：从四队列健康供者中按同一 Joint20 规则冻结 nonIBD eligible pools 和 recovery 细胞；不修改既有 IBD 输入或 85-state 模型。
- `METHOD_IDENTITY_CORRECTION_20261123.md`：记录 G7/G7A 当前结果是算法对齐的 IBD 内部扩展，而非正常/疾病完整复现，并冻结 nonIBD 补齐边界。
- `NONIBD_BULK_VALIDATION_PLAN_20261123.md`：冻结 MSCCR Control、GSE57945 及其他候选 nonIBD bulk 在当前 85-state/8-ecotype 模型中的输入尺度、独立性和 recovery 方案。
- `MSCCR_CONTROL_INTERNAL_NONIBD_COVERAGE_20261124.md`：记录当前模型在 461 个 MSCCR Control 样本中的 CPM 重建、官方 recovery、85-state/8-ecotype 覆盖、环境失败审计和外部证据整合。
- `RELEASE_GATES.tsv`：每一阶段开始和放行所需的可审计条件。

## 当前结论

旧包 `MSCCR_SignatureD_downstream_ready_20261114` 仍可作为描述性分析和绘图原型，但在本方案实施前，不应作为“完全符合 Cell 2021”的主分析 runner。BayesPrism 对照不再位于原文主链；外部 scRNA state recovery、marker LOOCV 和 scRNA ecotype detection 被提前为核心验证。

## 实施进度

- Phase 0 本地方法证据审计已实现：`scripts/00_audit_local_method_provenance.R`。
- Phase 0 正式输入冻结已实现：`scripts/00_freeze_inputs.R`，等待在一次获授权的持久服务器会话中读取完整正式数值输出。
- Phase 1 三层 state-QC 只读影响审计已实现：`scripts/01_audit_discovery_state_qc.R`；它只重算 marker、AFI 和 dropout 判定，不修改模型。
- 本地证据已确认正式配置跳过 EcoTyper Step 7；官方源码将 Step 7 定义为 advanced dropout QC。因此当前不修改模型，但 G2 三层 state-QC 放行门槛保持未通过。
- 官方 `1-Jaccard` ecotype partition 保持主结果；STAR 文字版 Euclidean comparison 降为非阻断敏感性，不作为重跑理由。
- 尚未启动 dropout 复算、state 删除、ecotype 重聚类或任何外部 recovery。

## 2026-11-24 更新

- 四套 nonIBD scRNA recovery、marker validation、GSE57945 外部 bulk recovery 和 Cell 2021 已知表型富集已完成；旧的“尚未启动外部 recovery”仅描述本方案初建时状态。
- GSE57945 与四套 scRNA 的 nonIBD 显著性一致性已逐数据集计算，不跨数据集合并样本。
- MSCCR Control 已按完整原始文库分母重建当前模型基因并集的线性 CPM；未重建 TPM、未重跑 discovery CD/UC。
- MSCCR Control 三部位官方 bulk recovery 已验收：60/60 state abundance、60/60 state assignment、3/3 ecotype abundance 和 3/3 ecotype assignment。
- 85 个正式 state 中 82 个在 MSCCR 三部位均观察到，84 个至少在两个部位观察到；8 个 ecotype 在三个部位均观察到。
- 外部 GSE57945 bulk 与至少一套 scRNA 同时 `Z > 1.96` 的 state 为 `CD8_T_S06`、`Fibroblast_S02` 和 `Ileal_absorptive_S02`。MSCCR 只作为内部覆盖证据，不计入该显著性判定。
- discovery bulk 的 Fig. S2G/S2H representation 已按原文统计单位补齐：Disease 和 Region3 分开计算与作图，全部可评估组合均进入箱线图，不按 OR 或显著性预筛选；不再保留 Disease × Region3 机械交叉图。
