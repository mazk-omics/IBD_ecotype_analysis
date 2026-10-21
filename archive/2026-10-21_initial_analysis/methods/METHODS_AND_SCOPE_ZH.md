# IBD EcoTyper 初步分析归档：方法与范围

本目录是公开仓库中的可追溯代码归档，不是结果发布包。归档代码保持原字节不变；患者级、样本级、原始表达、模型对象、结果表和科研图件均未纳入。

## 1. 冻结 BayesPrism CORE12 / REGION_IBD 敏感性路线

- 分析对象为冻结的 2025 份 CD+UC 活检，涉及 919 位患者。
- 冻结模型包含 53 个 cell states 与 7 个 ecotypes。
- 上游表达来自 BayesPrism G-like 路线，只能解释为敏感性分析；不得称为 CIBERSORTx 正式主线结果。
- 官方 EcoTyper discovery、state/ecotype recovery、单细胞共现、表型和临床审计代码按原运行目录镜像到 code/。
- 旧 BayesPrism 注释和对应性审计与 Joint20 路线分别存放，不能合并为同一发现模型。

## 2. 官方 EcoTyper recovery 与审计

- 使用冻结 discovery 模型，不重新拟合状态或生态型，不改变既有 assignment。
- recovery、marker/LOOCV、表型富集、区域代表性、ecotype 共现和临床关联均保留各自统计口径。
- 置换、Fisher、BH 家族、重复患者处理和外部队列限制以归档脚本为准。
- GlobalDetection 和 LigandReceptor 模块未达到原文等价复现条件，未用替代分析冒充完成。

## 3. Joint20 v2 Signature 与伪混合验证

- Joint20 v2 参考包含 20 个父类型、9,294 个训练细胞、32 位供者和 212 个 donor/site/sample/type 汇总参考 profile。
- CIBERSORTx 网页端使用 RNA-Seq reference samples + phenotype classes 模式，不按单细胞模式上传。
- 外部伪混合验证使用 GSE214695 和 IBDverse；每数据集 100 个混合物、每个 1,000 个细胞、有放回抽样。
- BayesPrism 使用原始 counts；同批 CIBERSORTx 输入使用 CPM。两种估计量的尺度和基因空间不同，不能称为纯算法单变量比较。
- 已完成 Signature 构建以及伪混合/Job25/Job26 配对复核；正式 HiRes 和正式 Discovery 不在本归档范围。

## 4. Joint20 matched BayesPrism 分支

- 与主线冻结 CD+UC 样本及 Joint20 父类型来源匹配。
- BayesPrism 去卷积和 QC 已完成；正式 EcoTyper NMF 未启动，记录为 ecotyper_nmf_started=false。
- 该分支不能替代或重命名为正式 CIBERSORTx Discovery，也不能复用旧 53 states / 7 ecotypes 作为新结果。

## 5. 临床与外部结局审计

- 临床分析脚本保留疾病、炎症、部位及重复活检处理的预设口径。
- 外部结局分析保持队列、终点、有效分母和协变量敏感性分析分离。
- 本公开归档不包含患者级明细、样本 ID、模型输出、效应量、P 值或图件。

## 6. 复现前提

- 需要合法取得各公开数据集和本项目受限输入，并按原脚本约定准备路径、基因标识及 metadata。
- 多个脚本保留服务器绝对路径，这是运行现场的原始记录；复现时应通过独立配置映射路径，不要修改本归档快照后再声称字节一致。
- CIBERSORTx 网页任务需要相应账号、网页输入和原任务输出；本仓库不包含认证状态或大体积上传/下载文件。
- 服务器 EcoTyper 软件树在审计时不是干净官方 checkout；必须结合记录的 revision 和 dirty 状态解释。
