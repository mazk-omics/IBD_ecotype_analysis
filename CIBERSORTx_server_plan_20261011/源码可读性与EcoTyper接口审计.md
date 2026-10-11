# CIBERSORTx 源码可读性与 EcoTyper 接口审计

日期：2026-10-11。范围：本机已下载的四个官方公开应用层，以及官方论文和 EcoTyper 调用代码。未连接服务器，未执行容器或编译程序，未修改研究数据或认证机制。

这是第一轮结构检查记录。后续已完成三个镜像的全部层扫描及核心数值函数反编译，最新范围、结果和限制请阅读[完整镜像调查与数值源码恢复](完整镜像调查与数值源码恢复.md)。下文保留第一轮的判断，不作为后续工作尚未完成的声明。

## 当前判断

**部分官方 R 实现可直接阅读；Fractions 和 HiRes 的主程序是编译后的 C++ ELF 文件，不能通过解压得到完整原始源码。** 两个主程序仍保留符号表，后续可以进一步定位算法函数、检查调用关系和尝试恢复伪代码。但本轮只完成静态结构审计，尚未反编译，也没有得到可以替代官方 HiRes 的完整实现。

“解压镜像”“反编译程序”和“恢复作者原始源码”是三个不同层次。此次已经能识别应用层中的可读代码和编译程序；保留符号表有利于继续调查，却不能恢复原来的注释、全部变量名、构建配置，也不保证数值等价。镜像公开下载与程序运行授权同样是两件事，当前调查没有解决账号审批的问题。

## 检查过程与证据

1. 复用 2026-10-10 保存的四个 COPY 应用层，重新核对各层 SHA256；压缩文件合计 1,751,604 字节。检查对象是小型应用层，不是整个容器的所有依赖层。
2. 在内存中读取文件，识别 ELF 头、节表、动态依赖、函数符号数量和编译单元名称。没有输出全部二进制字符串。
3. 对 R 文本记录包名、相关函数名和预设数值调用的出现次数，没有复制原始程序或字面量到报告。
4. 对照官方论文区分群体表达推断与样本级表达推断，再检查官方 EcoTyper 的 HiRes 调用接口。

复现脚本：[audit_cibersortx_source_structure_20261011.py](../scripts/python/audit_cibersortx_source_structure_20261011.py)。运行环境：本机 `.tools/cibersortx_review` Python 环境，仅使用标准库。

结构证据：[application_structure_audit.json](application_structure_audit.json)、[readable_R_algorithm_structure.json](readable_R_algorithm_structure.json)、[ELF_section_inventory.tsv](ELF_section_inventory.tsv)。每个应用文件均记录大小和 SHA256。

## 实际找到的模块

| 模块 | 实际形式 | 本次确认的内容 | 对项目的帮助 |
|---|---|---|---|
| CIBERSORTxHiRes | Linux x86-64 ELF，3,829,200 字节 | 有 `.symtab` 和 `.dynsym`，无 `.debug_info`；编译单元包含 `CIBERSORTxHiRes.cpp`、`GETable.cpp`、`LinearRegression.cpp` | 能定位部分核心结构，但尚不能恢复样本级表达算法 |
| CIBERSORTxFractions | Linux x86-64 ELF，3,374,088 字节 | 同样有符号表、无调试信息；编译单元包含 `SVRWrapper.cpp` 和 `svm.cpp` | 可进一步核对比例求解与 R 辅助程序的关系 |
| CIBERSORTxGEP.R | 可读 R 文本，1,087 行 | 有 GEP 相关函数定义；出现 `nnls`、`lm`、`normalize.quantiles` 等调用标记 | 可用于理解群体表达推断；不能仅凭这些内容认定已得到 HiRes 核心 |
| run_SVR.R | 可读 R 文本 | 出现 `svm` 调用，加载 `e1071` 等包 | 提供比例求解辅助实现的线索 |
| run_adjust_w_pseudo.R、run_combatba.R | 可读 R 文本 | 出现 `combatba` 调用 | 提供批次校正辅助实现的线索 |
| run_optimization.R | 可读 R 文本 | 定义 `performOptimization` | 需要进一步核对调用者、目标函数和约束 |
| 分位数归一化、热图相关 R 文件 | 可读 R 文本 | 找到归一化和绘图模块 | 有助于接口核对和结果呈现，不等于表达推断核心 |

HiRes 与 Fractions 分别记录到 12,311 和 10,619 个函数符号，其中包含第三方库和模板生成的函数，不能将这些数量理解为同等数量的独立算法步骤。HiRes 中存在 `Rscript`、`run_optimization` 和 `CIBERSORTxGEP` 标记，提示原生主程序可能调用外部 R 模块；这属于静态线索，尚未用控制流或运行实验确认。

R 文件中的数值调用次数是文本匹配统计，可能包括注释或字符串，并非运行时调用次数。没有匹配到某个名称也不能证明相应功能不存在。

## 为什么可读的 GEP 不能直接顶替 HiRes

原论文区分两种表达输出：group mode 推断每个细胞类型的群体代表表达；high-resolution mode 推断每个细胞类型在各样本中的表达。后者还涉及逐基因排序、窗口推断、平滑和样本顺序恢复。把一个 NNLS 群体表达结果复制到所有样本，并不实现样本级推断。[CIBERSORTx 原论文及 Methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC6610714/)

因此，目前的可读 R 部分可以减少算法理解的工作量，但不能据此宣布已经拿到完整 HiRes，也不能将它直接作为正式 EcoTyper discovery 的输入。

## EcoTyper 接口已经可以核对哪些内容

官方 `pipeline/csx_hires_worker.R` 在 Docker 和 Singularity 路线中都调用 HiRes，并传入账号与 token。输入包括 `--mixture`、`--sigmatrix`、`--classes`、`--cibresults`，另外设置线程数和关闭热图。代码将已有类别及比例结果复制到 worker 输出目录；它本身是调用包装，不包含 HiRes 的数学实现。[官方 EcoTyper worker](https://github.com/digitalcytometry/ecotyper/blob/master/pipeline/csx_hires_worker.R)

本项目需要核对的目标输出仍是各细胞类型的“基因 × 样本”表达矩阵，样本名称、基因名称、细胞类型名称、表达尺度及输入比例必须能与下游流程对应。当前结构审计未验证这些输出，也未证明旧 Docker 与现在网页结果数值一致。

检查所对应的 HiRes 镜像摘要为：

`cibersortx/hires@sha256:e8da6850311d163e33a343d29a0d2ffc8b18c1ec4604995b8285c7bc2017c83e`

来源：[官方 HiRes 镜像页面](https://hub.docker.com/r/cibersortx/hires)。旧镜像与当前网页运行环境不同，需要分别记录版本，不能只按工具名称视为同一实现。

## 下一步调查的落点

有价值的后续技术工作是对 HiRes 的**数值算法部分**建立函数与数据流映射，核对它怎样调用 GEP/优化辅助模块、怎样输出各细胞类型的样本表达；随后才谈是否能恢复或独立实现缺失部分。认证逻辑不属于这一算法调查范围。

若需要自行实现，必须用现有网页 HiRes 的小规模结果核对输入输出、表达变化和误差，明确标为独立实现，而不能预先声称与官方一致。此次没有启动这一重构，也没有运行任何替代分析。

四个应用层中没有找到 LICENSE/COPYING 文件名，但本次没有检查整个容器及全部分发材料，因此不能据此判断完整许可内容。原始供应商代码和二进制没有加入项目 Git；报告与审计脚本仅保留结构证据。
