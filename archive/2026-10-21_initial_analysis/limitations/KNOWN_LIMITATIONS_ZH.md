# 已知限制与禁止性解释

1. 本归档不发布科学结果；不能从归档存在推断某个临床或功能关联已经显著。
2. 冻结 53-state / 7-ecotype 路线是 BayesPrism 上游敏感性分析，不是 CIBERSORTx 正式主线。
3. Joint20 matched BayesPrism 仅完成去卷积与 QC，EcoTyper NMF 未启动。
4. CIBERSORTx HiRes 与正式 Discovery 尚未完成，相关实验脚本、缓存、环境和局部 Job24 产物未归档。
5. paper-aligned recovery 的 GlobalDetection 与 LigandReceptor 仍为 BLOCKED；全部请求模块并未完成。
6. recovery Z 不是表型一致性、细胞身份真值或独立临床验证。
7. ecotype 共现不代表物理接触、配体受体通讯或因果机制。
8. 重复活检不能按独立患者解释；患者级或聚类处理必须按各脚本的预设口径执行。
9. RISK 等外部队列存在量纲/来源限制；GSE73661 的最终可评价范围受 W0 样本限制。
10. 早期非零退出日志已被后续成功重跑取代；不能将早期尝试作为最终失败，也不能删除其审计意义。
11. 服务器 EcoTyper revision 记录了实验性 Region3 scaling 修改，不应声明为干净官方版本。
12. 本仓库不含原始数据、大矩阵、模型对象或网页任务文件，因此不能仅凭本提交端到端重跑。
