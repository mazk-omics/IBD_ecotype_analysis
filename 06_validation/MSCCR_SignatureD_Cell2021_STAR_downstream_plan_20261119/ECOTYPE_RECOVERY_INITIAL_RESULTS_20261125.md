# EcoType recovery 与 discovery bulk 初步结果

## 结论

冻结模型保持不变：79 个 cell states 归入 E1-E8。四个预先指定的单细胞 atlas 合并后，按 Cell 2021 STAR Methods 计算一个 pooled F1/F2 与共现矩阵 C。单个 EcoType 检验中，E1、E2、E3、E4、E5、E7 显著，E6、E8 不显著；1,000,000 次联合置换中，没有一次使 8 个 EcoType 的成员内共现均值同时达到观察值，联合经验 P<1e-6。

Discovery bulk 共 2,025 个样本、919 名患者。EcoTyper official assignment 接受 1,937 个样本，88 个未分配；E1-E8 均有正式分配样本。Discovery bulk 用于模型内部重建和生物学表征，不称为独立验证。

## pooled scRNA EcoType 检验

四个 atlas 共 299 个生物学样本、134,810 个输入细胞、64,612 个正式分配细胞。先在样本和父细胞类型内计算状态比例：F1 分母为可分配到所有正式 states 的细胞，F2 分母为该父类型全部输入细胞。F1、F2 及各自 zero-to-NA 版本产生四个 Spearman 矩阵，严格平均为 C。

| EcoType | Z | Cell 2021 单侧 P | 1,000 次经验 P | 结论 |
|---|---:|---:|---:|---|
| E1 | 5.602 | 1.06e-8 | 0.000999 | 显著 |
| E2 | 9.265 | 9.72e-21 | 0.000999 | 显著 |
| E3 | 6.530 | 3.30e-11 | 0.000999 | 显著 |
| E4 | 6.350 | 1.08e-10 | 0.000999 | 显著 |
| E5 | 2.387 | 0.00850 | 0.00999 | 显著但较弱 |
| E6 | 0.895 | 0.185 | 0.176 | 不显著 |
| E7 | 4.650 | 1.66e-6 | 0.000999 | 显著 |
| E8 | 0.782 | 0.217 | 0.218 | 不显著 |

原文单个 EcoType 展示口径为置换零分布得到 Z 后转单侧 P，P<0.05 标星且不以 BH 替代。项目同时保存经验 P、BH8 和此前逐数据集 32 项敏感性结果。联合检验使用成员内共现均值，不等同于单个检验的成员内减成员外 AvgDif，因此“整体 P<1e-6”与“单个 6/8 显著”可以同时成立。

正式 Figure 4D 风格输出只画一张 pooled scRNA 热图。四个 atlas 的独立矩阵保留为异质性审计，不铺成四张主图。

## discovery bulk 重建

| EcoType | Official assigned samples | 占全部样本比例 |
|---|---:|---:|
| E1 | 366 | 18.1% |
| E2 | 269 | 13.3% |
| E3 | 178 | 8.8% |
| E4 | 347 | 17.1% |
| E5 | 71 | 3.5% |
| E6 | 137 | 6.8% |
| E7 | 221 | 10.9% |
| E8 | 348 | 17.2% |

连续 abundance 使用全部 2,025 个样本。Disease、Inflammation、Region3 分开分析：

- 炎症轴最强。E2、E3、E7 在 inflamed 中较高；E1、E4、E5、E6、E8 在 non-inflamed 中较高。452 名同时具有两类活检的患者内配对分析中，8 个方向均通过 BH8。
- 疾病轴较弱。患者中位数敏感性分析支持 E4 在 CD 较高、E3 在 UC 较高；E6 的 UC 趋势未通过 BH8。
- 部位轴次于炎症。连续 abundance 的一对其余检验提示 E2 偏 ileum、E3 偏 rectum、E5/E8 偏 colon；E1、E2、E3、E5、E8 的 Region3 全局检验通过 BH8。

最高丰度类别的 Fisher 富集与连续 abundance 回答的问题不同，不能要求每一项完全一致。表型结论以连续 abundance、患者敏感性和效应方向共同判断，不以单个极小 P 值直接命名机制。

## 方法边界与后续

- 本轮没有重新计算、删改或重聚类任何 cell state / EcoType。
- pooled scRNA 使用所有预先指定 atlas、样本和可定义成员，没有按 P 值挑选数据。
- discovery bulk 不是独立 validation；外部证据来自 pooled scRNA。
- 仍需继续 EcoType 成员功能注释、配体-受体分析、临床关联和空间转录组 recovery。外部 bulk recovery 可作为后续独立验证补充，但不混入本轮 discovery 表征。
