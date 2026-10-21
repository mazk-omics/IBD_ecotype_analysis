# 本机与服务器来源审计摘要

- 本机仓库：D:\IBD_ecotype\IBD_ecotype_analysis；归档基线为 origin/main 的 7f96eb8edca6d73114b69500d19fddea019977d3。
- 服务器项目根：/home/mazekai/IBD_EcoTyper；服务器仓库审计时落后本机一笔提交且存在大量未跟踪分析输出。
- 审计时没有发现正在运行的 BayesPrism、EcoTyper 或 CIBERSORTx 计算进程。
- 服务器大型目录包括约 3.4 GiB 的 CIBERSORTx、9.1 GiB 的 EcoTyper 和约 16 GiB 的 Joint20 matched sensitivity 结果；这些内容没有复制到 Git。
- Joint20 v2 构建摘要记录 20 types、212 profiles、9,294 cells、32 donors；网页输入模式是 RNA-Seq reference samples + phenotype classes。
- Joint20 matched BayesPrism 完成去卷积和 QC，但 completion status 明确 ecotyper_nmf_started=false。
- paper-aligned recovery 的 GlobalDetection 与 LigandReceptor 为 BLOCKED；早期非零日志已有最终成功重跑覆盖。
- 本地和服务器若干关键小表已做 SHA-256 对照；结果内容不进入公开归档，只保留本机源文件指纹。
