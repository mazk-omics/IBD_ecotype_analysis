# 当前 85-state / 8-ecotype 模型的 nonIBD bulk 验证方案

## 1. 已有数据资产

| 队列 | nonIBD 样本 | 当前角色 | 当前输入状态 |
|---|---:|---|---|
| MSCCR / GSE193677 | 461 个 Control、Non-inflamed 活检；243 位患者；Ileum 121、Colon 115、Rectum 225 | 内部投影与流程一致性检查；参与过上游反卷积，不能作为完全独立验证 | 旧模型存在 CPM 输入；当前主路线需按与 CD/UC discovery 相同的 Ensembl75 protein-coding gene-length 流程重建 TPM |
| GSE57945 | 42 个 nonIBD ileum bulk | 外部 nonIBD bulk 主验证 | 已有 canonical metadata 和线性 RPKM；需重新审计当前 85-state 模型基因覆盖后 recovery |
| GSE235236 | GEO 报告 8 个 HC sigmoid bulk | 外部结肠候选；与 GSE214695 同研究，不作为完全独立第二证据 | 尚未正式导入；需核对 RSEM 字段、健康定义、gene ID 和表达尺度 |
| GSE93624 | 含 RISK nonIBD ileum | 候选敏感性队列；可能与 GSE57945 供者/研究重叠 | 当前表达含负值且尺度身份未闭环，暂不进入 recovery |
| 用户新增 bulk | 待上传 | 新的外部验证层 | 上传后先冻结样本、疾病、部位、表达单位和供者重叠 |

旧目录 `nonIBD_healthy_recovery_20261010` 使用 12 父类型、53-state、7-ecotype 模型。旧结果不能直接并入当前模型；可复用其样本筛选、源文件审计、供者边界和表达预处理证据。

## 2. 正式输入原则

1. 当前模型、rank、state 和 ecotype membership 全部冻结；nonIBD 不重新 NMF，也不按结果删除 state；
2. MSCCR Control 使用与 2025 CD/UC 主路线相同的原始 counts、Ensembl75 protein-coding 过滤、gene length 和 HGNC 映射重建 TPM；不把旧 CPM 静默替换成当前主输入；
3. 外部队列保留其可审计的线性表达单位。RPKM、TPM、芯片标准化表达分别审计并独立 recovery，不做跨平台拼接矩阵；
4. 每个队列独立输出 gene coverage、样本分母、state abundance、state assignment、ecotype abundance、ecotype assignment 和 Unassigned；
5. 同一患者多部位样本保持患者身份，部位结果分层报告；样本数和患者数同时作为分母；
6. nonIBD 阴性结果保留，不触发模型重训或样本后验筛选。

## 3. 主分析与证据等级

- MSCCR：在 Ileum、Colon、Rectum 内分别比较 IBD 与 nonIBD，并补充患者去重/等权敏感性；定位为内部一致性证据；
- GSE57945：以当前模型 recovery 覆盖和 nonIBD state/ecotype representation 为外部 bulk 证据；与 RISK 相关队列先做供者/样本指纹排重；
- 其他外部队列：逐队列分析后再汇总，不将平台或部位差异解释为疾病效应；
- 跨队列汇总：state/ecotype 方向可作 meta-analysis，但只对彼此独立且结局定义一致的队列合并；共享供者或同研究队列作单独敏感性。

## 4. 与单细胞 nonIBD 的关系

nonIBD scRNA recovery 用于检验细胞层 state assignment、marker 支持和 phenotype enrichment；nonIBD bulk 用于检验组织层 state/ecotype abundance 与 representation。两者是互补证据，不能相互替代，也不能混用分母。
