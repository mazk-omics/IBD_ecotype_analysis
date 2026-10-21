# 2026-10-21 IBD EcoTyper 初步分析安全归档

这是一次面向公开 GitHub 仓库的代码与方法归档。它用于记录已经完成的分析流程、运行边界、软件版本和源文件指纹，不构成未发表结果的公开发布。

## 目录

- code/：从本机已完成分析目录复制的字节一致代码快照，保留原相对路径。
- configs/：公开安全的分析范围和软件版本摘要。
- methods/：方法、路线边界、代码版本说明和复现前提。
- provenance/：代码来源映射、执行状态、本机/服务器审计摘要、未入库结果的 SHA-256 指纹。
- limitations/：禁止性解释、未完成模块和复现限制。

## 明确未入库

未提交任何患者级或样本级表、原始/标准化表达矩阵、模型对象、CIBERSORTx 上传/下载文件、科研结果表、图件、文献全文、缓存、环境或密钥。预提交审计中识别出的 29 个汇总结果文件和 35 张主汇报 PNG 均继续保留在本地，不在本分支公开。

## 使用方式

先阅读 configs/analysis_scope.yml、methods/METHODS_AND_SCOPE_ZH.md 与 limitations/KNOWN_LIMITATIONS_ZH.md。复现时依据 provenance/code_source_map.tsv 核对代码 SHA-256，并自行合法准备未随仓库分发的输入数据。
