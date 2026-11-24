# G7/G7A 方法身份更正（2026-11-23）

## 结论

`output/11_known_phenotype_exact` 已复现 Cell 2021 的逐细胞标签置换算法、1,000 次 shuffle 和 `Z > 1.96` 阈值，但输入只含 IBD 细胞，因此它不是原文正常组织与疾病组织对照的完整生物学复现。该结果正式改称 `CELL2021_STATISTICAL_FRAMEWORK_IBD_WITHIN_DISEASE_EXTENSION`。

`output/12_cell_state_representation` 中 Disease、Inflammation、Region3 和 dataset 的 dominant-state Fisher 结果同样属于 IBD 内部扩展；在 nonIBD recovery 完成前，不得把它写成 Figure S2G-S2H 对应的 IBD/nonIBD representation。

本次更正不改变任何原始结果、置换参数、state assignment、85-state 模型或 79-state ecotype membership，只更正结果的生物学身份与证据级别。

## 补齐方案

1. 从 SCP1884、GSE214695、SCP259 和 IBDverse 的健康供者中按既有 Joint20 映射冻结 nonIBD 细胞；
2. 每个 dataset × parent type 使用固定种子 1234，最多抽取 2,500 个细胞；
3. 使用同一冻结 85-state 模型做 recovery，不重新 NMF、不修改 rank、不按结果筛 state；
4. 分队列、分父类型比较 nonIBD 与 IBD，保留 Unassigned、低细胞数和不可评价分母；
5. 再按 Cell 2021 的 1,000 次逐细胞标签置换计算 phenotype enrichment，并按样本 dominant state 做 one-sided Fisher representation；
6. 四套 scRNA 数据曾参与 Signature-D 表达纯化验证，但未参与本轮 state/ecotype discovery，因此表述为“正交于 state discovery 的外部 recovery”，不表述为完全未接触的独立验证。

## 解释边界

- 既有 IBD-only phenotype 结果保留为补充分析，不作为正常/疾病主验证；
- nonIBD 阴性结果不得触发模型重训、state 删除或健康样本重抽；
- cohort、部位和细胞量差异必须分层展示，不能把 pooled 差异直接解释为疾病效应；
- 低支持父类型预先按输入分母标记，不能在查看 enrichment 结果后调整阈值。
