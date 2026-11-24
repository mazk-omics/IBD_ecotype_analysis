# Cell 2021 STAR Methods 对齐的 Signature-D 下游分析方案

版本：2026-11-19。当前文件只冻结方法和决策边界，不读取运行中间文件，不连接服务器，不根据结果调参。

## 1. 目标与适用边界

主要目标不是寻找最好看的亚组差异，而是回答三个按顺序排列的问题：

1. 本轮正式模型的 cell states 是否按 Cell 2021 的质量控制和外部 scRNA-seq recovery 方法可复现？
2. 冻结的 ecotypes 是否按 Cell 2021 的单细胞共现检验可检测，并对替代聚类方法稳定？
3. 在前两项证据成立后，这些 states/ecotypes 是否在 IBD 的 Disease、Inflammation、Region3、样本和队列之间保守，并与临床特征相关？

“完全对齐”在本项目中指：可迁移的算法严格按 STAR Methods 执行；EcoTyper 官方代码与论文文字不一致时双轨审计；癌症特异终点明确标为不适用；所有 IBD 必需改动单独标为扩展。不能把 IBD 改写后的统计称为原文复现。

## 2. 分析轨道

| 轨道 | 定义 | 是否进入主证据链 |
|---|---|---|
| `CELL2021_EXACT` | STAR Methods 明文公式和参数的直接实现 | 是 |
| `OFFICIAL_SOFTWARE` | 冻结版本 EcoTyper 官方 R 代码的实际结果 | 是，且与原文文字逐项比对 |
| `IBD_EXTENSION` | 为重复活检和 IBD 表型增加的患者级模型与保守性检验 | 是，但不得写成 Cell 2021 原法 |
| `SUPPLEMENTARY_SENSITIVITY` | BayesPrism、删减高风险细胞类型、替代校正等 | 否 |
| `NOT_APPLICABLE` | 无同构 IBD 数据的肿瘤专属分析 | 否，不以近似替代 |

## 3. 冻结对象与再现性要求

正式运行前一次性冻结：

- 正式 EcoTyper 目录、`status.txt`、config、console、所有 state/ecotype 数值表及 SHA256；
- CIBERSORTx fractions、20 类 HiRes 表达矩阵、正式 TPM、样本注释和 PatientID；
- EcoTyper 源码逐文件 SHA256、容器镜像 digest、R 版本、包版本、BLAS、线程数和 RNGkind；
- 外部 scRNA-seq 四队列的原始计数、规范化矩阵、cell/donor/sample 注释、父细胞类型映射和基因 ID 映射；
- discovery、validation 和 sensitivity 数据集角色。当前外部单细胞队列已参与 signature/purification 评价，因此 state/ecotype recovery 应表述为“未参与 state/ecotype discovery 的正交终点验证”，不能表述为完全未接触的独立验证；
- 所有随机步骤的种子。建议总种子 `20261119`，各模块使用固定派生种子并保存实际抽样 ID。

禁止行为：根据 recovery、P 值、silhouette 或图形观感调整 state 数、ecotype 数、marker 数、细胞类型映射、样本排除或阈值。任何方法改变必须产生新版本目录，旧结果保持可复现。

## 4. Phase 0：正式输出与模型身份审计

### 4.1 输入完整性

- 样本集合、顺序和 ID 一致；重复 ID、缺失注释、非有限值、负 abundance、列和异常均单独列出。
- `fraction`、`state abundance`、`ecotype abundance` 和 HiRes expression 分开审计，不能互相替代。
- 核实 fraction 表的所有度量列，包括 `P-value`、`Correlation`、`RMSE` 和 `Absolute score`，不得误当细胞类型。
- 冻结官方 `AssignedEcotype`；另可计算 `DominantEcotype` 供描述，但并列保存且不得混用。
- 对 abundance 并列最大值不使用 `ties.method = "first"`。并列最大时离散标签设为 `NA/Tied`，连续 abundance 保留。

### 4.2 模型质量控制闭环

逐父类型复核 Cell 2021 三层 state QC：

1. marker 数：基因在至少 50% discovery 样本中非零，且该 state 相对其余 states 的 log2FC 最高；`n <= 10` 的 state 删除；
2. AFI：删除 `AFI >= 1`。原文对单个 epithelial state 的事后例外不能迁移到 IBD；本项目不设例外；
3. dropout：每个 marker 计算相同数值的最大样本比例和平均 log2 表达；分别在 states 间标准化，取两项平均 Z；删除 mean Z `> 1.96`；
4. 被删除 state 的原始 winning assignments 设为 unassigned，不能在剩余 states 中重新 argmax；
5. 保存 initial states、每层删除原因、最终 states、W、H、marker 和 assignment 的映射表。

若正式运行没有完成任一 QC 层，停止下游主分析。补做 QC 会改变模型身份，必须建立新模型版本，不能在原目录内覆盖。

## 5. Phase 1：state discovery 与 ecotype discovery 的原文复核

### 5.1 state discovery

- 每个父类型默认选择相对离散度最高的 1,000 个基因；分箱、relative dispersion 和基因集合与官方输出逐项核对。
- 每个基因在数据集/预设标准化层内 `log2` 后均值 0、方差 1，再做 posneg transform。
- NMF 使用 KL divergence；rank `2:20`；每个 rank 默认 50 次 restart；选择最低 RMSE 模型。
- H 每个样本列在父类型内归一到和为 1。
- rank 选择严格执行：在至少两个连续 rank 的 CC 高于 0.95 后，寻找首次低于 0.95 的 crossing，再取 crossing 两侧最接近 0.95 的 rank。IBD 不复制论文对 epithelial/neutrophil 的人工例外。
- 当前项目的 Disease/Region3 校正属于输入预处理适配，必须保存未校正、仅 Region3、Disease+Region3 三套设计的预定敏感性；正式主模型只使用已经冻结的设计，不根据下游结果更换。

### 5.2 ecotype discovery

1. 每个样本、每个父类型只把最丰富的保留 state 记为 1，其余为 0；QC 删除 state 的原 winning sample 保持 unassigned；
2. state 两两计算 Jaccard；超几何检验 `P > 0.01` 时令 Jaccard 为 0；
3. `CELL2021_EXACT` 按 STAR Methods 对过滤后的 Jaccard 行使用 Euclidean distance、average linkage；
4. `OFFICIAL_SOFTWARE` 按冻结官方 `ecotypes.R` 的 `as.dist(1 - Jaccard)`、average linkage重算；
5. 两轨均在合法 K 范围内最大化平均 silhouette，并删除成员数 `<= 2` 的簇；
6. 比较两轨的 K、ARI、成员 Jaccard 和逐 state 映射。正式软件 `1-Jaccard` 结果始终作为本项目主模型；Euclidean 结果是非阻断敏感性审计，不触发自动重跑。只有在差异明显时才写入局限性；若要改用 STAR 文字分区，必须另建新模型版本；
7. 按原文做 Louvain 和 k-medoids 稳健性：两者 K 搜索 `2:30`，最大化平均 silhouette，比较前删除少于 3 states 的簇；
8. 另做患者级稳定性：每位患者随机一份活检，至少 1,000 次固定模型重聚类，保存每次样本 ID、K、ARI、state membership stability。此项为 `IBD_EXTENSION`。

### 5.3 ecotype abundance 与 assignment

- 每个 ecotype 对成员 state abundance 取平均，再在每个样本内跨 ecotypes 归一到和为 1；
- 对每个样本/ecotype，用双侧 unequal-variance t test 比较该 ecotype 成员 states 与其他 ecotypes 的 states；
- 对候选最高 abundance ecotype 的 P 值按样本做 BH；只有 `q < 0.25` 且样本至少有一个 CE 成员 state 的有效 assignment 时才正式分配，否则 `Unassigned`；
- 下游连续分析优先使用 ecotype abundance；正式 assignment 用于复现原文类别分析；强制 argmax 只作描述性图，不能替代正式 assignment。

## 6. Phase 2：外部 scRNA-seq state recovery

外部队列固定为当前可评价的 `SCP1884`、`GSE214695`、`SCP259` 和 `IBDverse/E-MTAB-16999`。每个队列先冻结父类型映射、供者、样本、Disease、Inflammation、Region3 和可评价原因。

### 6.1 recovery

- 只使用 discovery 冻结的 W；对每个外部表达矩阵 M* 迭代更新 H'，不重新拟合 W；
- 输入按队列内每个基因 `log2` 并缩放到单位方差；分析单位是单细胞，且只在预先映射的父细胞类型内 recovery；
- 删除 discovery QC 未保留 states，得到 H''；每个 cell 分给父类型内最大权重 state；
- 每个父类型超过 2,500 cells 时，用冻结种子随机抽 2,500；保存实际 cell IDs；不得为了提高结果而反复抽样挑最好的一次；
- 每个队列分别保存输入覆盖率、W 基因覆盖率、总 cells、可分配 cells、每 state cells、每 donor cells 和失败原因。

### 6.2 recovery 显著性

- state marker 与 QC step 1 相同，但不使用“至少 50% 样本非零”的过滤；
- normalized counts 做 `log2`，每个基因在父类型 cells 中均值 0、方差 1；
- 对 state s 的每个 marker，计算 assigned-s cells 对其余 states cells 的平均 fold change，再对 marker 取平均得到 `FC_s`；
- 独立打乱每个 marker 在 cells 间的表达并重新求解 Equation 8，共 100 次；
- `Z_s = (FC_s - mean(FC_s,null)) / sd(FC_s,null)`；
- 各 scRNA 数据集 Z 用等权 Stouffer 合并。原文判定使用 `Z > 1.96`；原文未规定 recovery 层 BH，因此 exact 表保留名义阈值，另附 BH 结果但不反写原文判定；
- 分开报告“有输入”“有 assignment”“finite Z”“Z>1.96”，不能合并成一个成功率。

## 7. Phase 3：外部 marker 排名与整队列 LOOCV

### 7.1 七指标 marker 排名

在每个 dataset × parent cell type 中先排除 `<5` assigned cells 的 state 及其 cells。对每个 gene-state 计算：

- `n1`：基因平均 TPM/CPM > 0 的数据集数；
- `n2`：基因被分配给该 state 的数据集数；
- `n3 = n2/n1`；
- `n4`：至少一个数据集把基因分给该 state 的独立疾病层数；本项目以预先冻结的 IBD disease strata 替代原文 tumor types，并明确标记为 `IBD_EXTENSION`；
- `n5`：Seurat v3.1.3 `FindMarkers(min.pct=0.1, logfc.threshold=0.05)` 中 `Q<0.05` 的数据集数；
- `n6`：各数据集 nominal P 转带方向双侧 Z，Stouffer 合并，meta P 做 BH 后再转带方向 Z；
- `n7`：各可评价数据集 log2FC 的平均值。

基因先按最大 log2FC 分配给 state；平局用跨数据集平均 log2FC 最大者；若该归属与 EcoTyper state 不一致则排除。`n1:n7` 分别转 rank 后等权平均，得到 composite marker score。手工 marker 不进入验证准确率，只能在注释表单独标记。

### 7.2 LOOCV

- 留出单位是完整 scRNA-seq dataset，不是 donor；
- 每次只在其余数据集重算全部七指标和每 state top 10 markers；
- held-out 数据中每个基因 `log2` 并在父类型所有 cells 中缩放到单位方差；
- 对每个真实/recovered state，计算 top 10 markers 的平均表达，最高者为预测 state；
- 输出逐 dataset、逐 parent、逐 state 是否正确和总体 LOOCV rate；缺少 10 genes、少于 5 cells 或无 recovery 的 state 均保留在分母审计表中，不悄悄删除。

## 8. Phase 4：表型 recovery 与 discovery-external concordance

### 8.1 Cell 2021 exact 表型富集

对每个预先定义的 phenotype P 和 state s，统计被分配至 s 的 P 样本/cells 数 `N_s`；将 state assignment labels 打乱 1,000 次得到 null；计算：

`Z_s = (N_s - mean(N_s,null)) / sd(N_s,null)`。

原文使用 `Z > 1.96` 为富集。该轨道不把简单 Wilcoxon 或 Fisher 当成替代。discovery 与 external scRNA 使用同一预先固定 contrast，随后报告方向一致率和预定的 concordance 检验。

### 8.2 IBD 患者级扩展

原文置换默认样本/cell 可交换，不适合 MSCCR 重复活检作为主要推断。因此另外运行：

- 每位患者一份冻结抽样；
- PatientID 内/疾病与部位 strata 内受限置换；
- 连续 abundance 的混合模型或 GEE：`abundance ~ Disease + Inflammation + Region3 + covariates + (1|PatientID)`；
- state 和 ecotype 分开定义 BH family，Disease、Inflammation、Region3 不混成一个事后 family；
- 同时报告效应量、CI、患者数、活检数和可评价 strata，不只报告 P/Q。

Exact 与 IBD extension 并列表述；若二者不同，以患者级结果作为 IBD 生物学推断，原文轨道仅说明算法复现。

## 9. Phase 5：scRNA-seq ecotype detection

### 9.1 四矩阵 co-occurrence index

state assignment 必须在不知道 ecotype membership 的条件下完成。以生物学 sample 为列、ecotype 成员 states 为行：

- `F1`：分母仅为能分配到 EcoTyper states 的同父类型 cells；
- `F2`：分母为该 sample 中该父类型全部 cells；
- `C1/C2`：分别对 F1/F2 的 state 行做 Spearman；
- `C3/C4`：把 F1/F2 的零替换为 NA 后再做 Spearman；
- `C = mean(C1,C2,C3,C4)`，逐元素平均并保存每个元素实际可用矩阵数。

### 9.2 单个 ecotype 检验

令 C 对角线为 NA。对 ecotype i 内每个 state s，计算其与同 ecotype 其他 states 的均值 `mu1`，以及与其余 ecotypes states 的均值 `mu2`；`Delta_s = mu1 - mu2`；`AvgDif_ref = mean(Delta_s)`。独立置换 C 的每一行 1,000 次，并按原文把 NA 对角位置与新对角位置交换，得到 `AvgDif_null`：

`Z_i = (AvgDif_ref - mean(AvgDif_null)) / sd(AvgDif_null)`。

Z 转为 P；原文以单个 ecotype `P<0.05` 报告，不做 BH。补充表可给 BH，但不能替换原文 exact 结果。

### 9.3 全局联合检验

按 STAR Methods 对每个 ecotype/state 计算同 ecotype mean co-occurrence，逐行置换并重复 1,000,000 次；每次统计有多少 ecotypes 的置换 co-occurrence 同时达到或超过 reference。输出 joint empirical P、Monte Carlo 分辨率和计数。若公式、比较方向或 NA 交换不能从原文公式与历史官方实现完全复原，则该模块不得用近似检验冒充；必须先完成小矩阵人工例和原始结果回归测试，才能放行。

### 9.4 IBD 保守性扩展

- 分别在 CD、UC、nonIBD（若存在）、Inflamed/Non-inflamed 和 Region3 中重建 C，但 state/ecotype membership 始终冻结；
- leave-one-Disease、leave-one-Region3、leave-one-dataset 重算 detection，不重新选择模型；
- patient bootstrap 计算 AvgDif、成员相关和 detection P 的分布；
- “跨疾病/跨部位保守”只在预先规定的可评价 strata 中报告，不能以 pooled 显著掩盖某一层完全缺失。

## 10. Phase 6：原文 feature analysis 与 IBD 临床扩展

### 10.1 可原样迁移的 feature analysis

- Hallmark：bulk TPM `log2`；每个基因跨样本缩放到单位方差；基因集内取平均；
- 连续特征：按正式 ecotype assignment 分组，某 ecotype 对其余 ecotypes 做双侧 Wilcoxon；效应为组内均值减其余组均值；
- 离散特征：按原文比较两组间连续 ecotype abundance，双侧 Wilcoxon；
- 所有预先列出的 features 构成一个明确的 BH family；`q<0.05`；输出 family 大小、缺失和效应单位；
- CIBERSORTx 细胞比例仅按最高 ecotype 类分组后取均值并用于组成显示；显示缩放值不能用于推断。

### 10.2 IBD 临床扩展

- Disease、Inflammation、Region3 同时进入模型，PatientID 作为随机截距或 cluster；
- 年龄、性别、用药、采样时间和 cohort 只有在预先核实可比性后进入；
- treatment response 只用治疗前样本，逐 cohort recovery 后分析；有 time/event 才使用 Cox，无 time/event 不用二分类 OR 冒充 HR；
- progression/response 的每个 endpoint 单独冻结 eligible sample、事件数、缺失规则和 BH family；
- 结果表必须区分 association、prediction 和 validation；未做校准和外部预测时不得称预测模型。

### 10.3 原文肿瘤专属模块

pan-cancer survival meta、ICI 跨治疗 Liptak 排名、癌症突变 signature、肿瘤 stage 和正常-肿瘤场效应在无同构 IBD 数据时标记 `NOT_APPLICABLE`。这不是缺失结果，也不以现有横断面 IBD 分析替代。

## 11. Phase 7：ligand-receptor、空间与功能解释

若要声称复现 Cell 2021 ligand-receptor enrichment，必须严格使用：

- HiRes 纯化表达；每个 state 内 CE-vs-other CE 的 Welch t test；BH；`q<0.05` 且 `log2FC>0.1`；
- 原文 ligand-receptor 数据库的冻结版本；
- 每个 state pair 随机抽取相同数量基因 100 次形成 null；
- state-pair Z 通过 Stouffer 汇总到 ecotype pair。

数据库版本、基因 universe、方向和 undefined pair 未闭环前，该模块保持阻断。普通相关网络、CellChat 或 bulk 共现不得替代。空间共定位需要独立空间数据；没有空间数据时标为 `NOT_APPLICABLE`。

## 12. Phase 8：敏感性分析

以下均不改变主模型：

- BayesPrism 与 CIBERSORTx 的 fraction、marker 和 ecotype 对照；
- 完整 20 类、去 ILC、去 ILC+DC、去预先定义高风险稀有类型；先投影到冻结模型，不重新 NMF；
- STAR 文字距离与官方代码距离；
- nominal P 与补充 BH；
- 一活检/患者、全活检混合模型和受限置换；
- 不同可评价 scRNA 队列组合的 leave-one-dataset-out。

任何敏感性结果较好都不能替换主结果；若提示模型身份应改变，只能提出下一版模型，不覆盖本版。

## 13. 图表计划

| 图 | 内容 | 方法身份 |
|---|---|---|
| Figure 1 | rank/CC、三层 state QC、正式 state heatmaps、UMAP | Cell 2021 exact + official audit |
| Figure 2 | 外部 recovery Z、逐队列覆盖、七指标 markers、整队列 LOOCV、表型 concordance | Cell 2021 exact |
| Figure 3 | Jaccard、ecotype network、silhouette、Louvain/k-medoids、STAR-vs-code partition | Cell 2021 exact + official audit |
| Figure 4 | scRNA C 矩阵、单 ecotype AvgDif/P、1e6 global detection | Cell 2021 exact |
| Figure 5 | Disease/Inflammation/Region3、Hallmark、组成与临床特征 | Cell 2021 feature method + IBD extension |
| Figure 6 | leave-one-disease/region/dataset、patient bootstrap、风险细胞类型敏感性 | IBD extension |
| Supplement | BayesPrism 对照、替代 BH、全部审计与失败项 | sensitivity/audit |

所有图由 R 生成，PDF/SVG/TIFF 同源，source data 单独导出。图不能从截图或图形对象反推数值。

## 14. 实施脚本顺序

正式实现时使用线性 R 研究脚本，不建立通用 pipeline 框架：

1. `00_freeze_inputs.R`
2. `01_audit_discovery_state_qc.R`
3. `02_audit_ecotype_discovery.R`
4. `03_recover_scrna_states.R`
5. `04_test_state_recovery.R`
6. `05_rank_markers_and_lodo.R`
7. `06_test_phenotype_concordance.R`
8. `07_test_scrna_ecotype_detection.R`
9. `08_test_ibd_conservation.R`
10. `09_characterize_ecotypes.R`
11. `10_run_clinical_extensions.R`
12. `11_run_sensitivity_analyses.R`
13. `12_build_figures.R`
14. `13_validate_delivery.R`

每个脚本顶部列明输入和输出，参数靠近实际调用，保存中间对象、种子、抽样 ID、sessionInfo 和错误表。耗时 recovery 可分父类型运行，但聚合前必须校验模型 hash、输入 hash 和参数完全一致。

## 15. 放行与停止规则

- Phase 0–1 未通过，不运行或解释后续主分析；
- discovery QC 缺失会改变 state 身份时，停止并新建模型版本；
- 外部数据缺少父类型、供者或样本映射时，只把对应 state/dataset 标为不可评价，不以零补齐；
- recovery 结果差不触发自动重跑、改 marker、删基因、改 state 或挑选抽样；
- silhouette 低、ecotype detection 不显著或保守性不足均是结果，不是自动调参理由；
- 只有输入、模型、代码、环境、随机性、输出和统计回归测试全部闭环，才可把模块状态标记为 `COMPLETE`；
- 主结论需同时报告负结果、不可评价分母和模型局限，不以敏感性分析覆盖主结果。

## 16. 对旧下游包的处置

`MSCCR_SignatureD_downstream_ready_20261114` 中可复用的是输入审计框架、统一长表、Hallmark 计算、R 图形 QA 和 source-data 合同。需要重写或后移的是：

- biopsy-level Wilcoxon/Fisher 不能替代原文 phenotype permutation；
- 强制 `max.col(..., ties.method="first")` 必须移除；
- fraction 度量列必须补入 `Absolute score`；
- external scRNA recovery、marker LOOCV 和 global ecotype detection 从 optional 提升为核心；
- BayesPrism 对照移到补充敏感性；
- 原文不适用模块明确标记，不再用“近似足够”作为完成标准。
