# State annotation v4 审阅记录

本轮在冻结的 85 个 state 和 8 个 ecotype 上继续注释，不重跑 CIBERSORTx、NMF 或 ecotype discovery，也不删除成员。注释对象是可重复的转录状态或功能程序；只有证据充分时才解释为稳定细胞亚型。

## 1. v4 注释优先级主表

`scripts/54_build_state_annotation_priority_master.R` 合并以下证据：

- 当前文章标签和 A-D 文献证据等级；
- published-signature 同谱系重叠；
- 四套 scRNA recovery 及 Stouffer Meta Z；
- GSE57945、MSCCR nonIBD 和跨模态验证；
- Disease、Region3 marker LODO；
- marker 通路富集；
- discovery bulk 的 inflamed/noninflamed representation。

该脚本只生成审阅优先级，不自动改变标签或等级。85 个 state 的等级仍为 A=7、B=18、C=30、D=30。

## 2. 首批高优先级审阅

首批对象为 10 个同时满足“D 级、ecotype 成员、恢复证据较高、需要人工判断”的 state。逐条记录位于 `inputs/state_annotation_priority_review_v4.tsv`，并由 `scripts/55_integrate_priority_state_annotation_review.R` 写入完整证据表。

- 10 个 state 全部保持 D 级，没有证据升级；
- 8 个保留为未解析的母谱系归属程序，仅放补充材料；
- `NK_S03` 和 `Plasma_S05` 的 marker 更接近炎症/损伤上皮信号，记录为跨谱系敏感性信号；
- 所有 10 个 state 均保留在原 ecotype 定义中，`ModelDefinitionChanged=FALSE`。

这里的“跨谱系”是解释警示，不等同于已经证明的污染或错误细胞亚型。是否影响 E2 需要在后续 ecotype 敏感性分析中检验，不能通过事后删成员代替。

## 3. Inflamed/noninflamed 供者感知敏感性

正式 inflammation 分析使用 discovery bulk，完全沿用 Fig. S2G/S2H 的 dominant-state、单侧 Fisher 和表型水平内 BH 方法。`scripts/56_integrate_inflammation_patient_sensitivity.R` 仅复用外部 scRNA known-phenotype 结果，对逐细胞 exact 轨道和供者感知轨道进行补充敏感性汇总，不重新计算 EcoTyper。

- 供者感知轨道中共有 11 个 state 在至少一个可评估数据集显著；
- 其中 10 个与 discovery bulk 的 inflamed/noninflamed 方向一致；
- `NK_S03` 和 `Plasma_S05` 在逐细胞 exact 轨道出现 inflamed 证据，但未通过供者感知复现；
- 因此二者不据此提升注释等级，也不作为主文 inflammation 结论。

供者感知比较只在 SCP1884 和 SCP259 可评估。IBDverse 的逐细胞 exact 结果不能自动视为独立供者复现；本结果作为补充敏感性证据，不替代 discovery 样本级 Cell 2021 Fisher/BH representation。

## 4. 第二批解释边界审阅

第二批覆盖 11 个恢复较稳定但仍未解析的 D 级 state，以及 15 个属于 ecotype 的 C 级功能程序。决策表位于 `inputs/state_annotation_priority_review_v4_batch2.tsv`，由 `scripts/57_integrate_state_annotation_review_batch2.R` 验证并整合。

- 11 个 D 级 state 全部保持补充材料使用，不因 recovery Meta Z 较高而升级；
- 15 个 C 级 state 全部保持功能程序，不解释为稳定细胞亚型；
- `Enteroendocrine_S02`、`Monocyte_Macrophage_S02` 和 `Paneth_S04` 保留跨谱系敏感性标记；
- `B_S01`、`DC_S03`、`Endothelial_S03` 和 `Mast_S01` 可作为功能程序，但保留谱系解释警示；
- 标签、证据等级、state 和 ecotype 定义均未改变。

## 5. 当前文章使用边界

- A-B：可作为正文命名状态使用；
- C：只作为 marker 支持的功能程序使用，不宣称为稳定亚型；
- D：只放补充材料，保留未解析或跨谱系警示；
- state 和 ecotype 的原始定义保持不变；任何成员排除只能作为预先声明的敏感性轨道。

服务器验收状态为 `PASS_STATE_ANNOTATION_PRIORITY_REVIEW_V4_BATCH1`、`PASS_STATE_ANNOTATION_INFLAMMATION_SENSITIVITY` 和 `PASS_STATE_ANNOTATION_PRIORITY_REVIEW_V4_BATCH2`。
