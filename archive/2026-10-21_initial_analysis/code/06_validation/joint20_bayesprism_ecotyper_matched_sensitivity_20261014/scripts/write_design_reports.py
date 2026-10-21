# 汇总已通过的输入审计和首次后台状态，生成本次交付文件。
# 只写独立敏感性目录；长任务状态以日志文件为准，不假定完成。
from pathlib import Path
import csv
import json
import sys
from datetime import datetime

# Step1. 读取本轮实际审计、PID与冻结证据
out = Path(sys.argv[1])
audit = json.loads((out / '00_audit/input_audit.json').read_text())
local = json.loads((out / '00_audit/local_mainline_evidence.json').read_text())
snapshot_time = datetime.now().astimezone().isoformat()
stage_path = out / 'logs/current_stage.txt'
stage = stage_path.read_text().strip() if stage_path.exists() else 'NOT_STARTED'
pid_path = out / 'logs/formal_R.pid'
pid = pid_path.read_text().strip() if pid_path.exists() else 'NOT_STARTED'
launcher_path = out / 'logs/launcher.pid'
launcher = launcher_path.read_text().strip() if launcher_path.exists() else 'NOT_STARTED'
started_path = out / 'logs/started_at.txt'
started = started_path.read_text().strip() if started_path.exists() else 'NOT_STARTED'
command_path = out / 'logs/formal_command.txt'
command = command_path.read_text().strip() if command_path.exists() else 'NOT_STARTED'
session_path = out / 'logs/tmux_session.txt'
session = session_path.read_text().strip() if session_path.exists() else 'NOT_STARTED'
exit_path = out / 'logs/exit_status.txt'
exit_status = exit_path.read_text().strip() if exit_path.exists() else 'RUNNING; no exit code yet'
retained_path = out / '02_bayesprism/retained_genes.txt'
retained = len(retained_path.read_text().splitlines()) if retained_path.exists() else 'pending new.prism'

# Step2. 冻结配对设计对照，未决项保留WAIT_MAINLINE
rows = [
 ['discovery_samples', '2025 frozen Job16/metadata/CPM IDs', str(audit['samples']), 'EXACT_SET_AND_ORDER', audit['ordered_sample_sha256']],
 ['CD', str(local['disease_counts']['CD']), str(audit['disease_counts']['CD']), 'EXACT', 'frozen metadata'],
 ['UC', str(local['disease_counts']['UC']), str(audit['disease_counts']['UC']), 'EXACT', 'frozen metadata'],
 ['patients', str(audit['patients']), str(audit['patients']), 'EXACT', 'PatientID; biopsies are not independent patients'],
 ['metadata', 'ID/PatientID/Disease/Inflammation/Region3', 'same values and sample order', 'EXACT', 'local-server annotation byte identical; legacy metadata checked'],
 ['NonIBD_controls', 'excluded', 'all 461 excluded', 'EXACT', 'separate from four frozen technical IBD exclusions'],
 ['technical_exclusions', '3 CD and 1 UC', str(audit['frozen_technical_IBD_exclusions']), 'EXACT', 'original technical_exclusion_samples.tsv, no new exclusions'],
 ['inflammation', 'Inflamed and Non-inflamed', 'both retained', 'EXACT', json.dumps(audit['inflammation_counts'])],
 ['training_source', '9294 GSE282122 cells, 32 donors, 212 CPM profiles', 'same 9294 original cells and 32 donors, raw counts', 'SOURCE_MATCHED_REPRESENTATION_DIFFERS', 'not algorithm-only identical numeric input'],
 ['reference_parent_types', '20 Joint20 v2 types', '20 same labels and frozen order', 'EXACT', '00_audit/frozen_types.txt'],
 ['input_gene_space', '4474 signature genes; 17040 bulk CPM background', f"{audit['bulk_mapped_genes']} mapped raw-count genes; {retained} prism genes", 'MODEL_SPECIFIC', 'same identity bridge; BP cleanup/outlier filter differs'],
 ['gene_identity', 'audited one-to-one Ensembl75/HGNC bridge', 'same accepted mapping', 'MATCHED_IDENTITY_POLICY', 'no new aliases or duplicate aggregation'],
 ['bulk_scale', 'nonlog CPM/full original library denominator', 'original numeric counts, fractional values preserved', 'MODEL_SPECIFIC', 'no global rounding or depth replacement'],
 ['purified_expression', 'complete 2025 HiRes not available', 'initial posterior Z to conditional G-like', 'WAIT_HIRES', 'not equal expression units or model'],
 ['NMF_parent_types', 'NOT_FROZEN', 'WAIT_MAINLINE', 'HOLD', 'do not substitute Core12 or force all 20'],
 ['Region3_scaling', 'NOT_FROZEN; historical templates disagree', 'WAIT_MAINLINE', 'HOLD', 'no legacy Region3 patch automatically reused'],
 ['NMF_gene_filter_rank_restarts_seed', 'NOT_FROZEN', 'WAIT_MAINLINE', 'HOLD', 'matched template includes historical values only as provenance'],
 ['state_QC_ecotype_clustering', 'NOT_FROZEN', 'WAIT_MAINLINE', 'HOLD', 'official presorted interface inspected; no NMF launched'],
 ['legacy_53_states_7_ecotypes', 'separate frozen branch B', 'never reused as branch C results', 'FROZEN', 'legacy output paths read only'],
]
with (out / 'MATCHED_DESIGN_MANIFEST.tsv').open('w', newline='') as handle:
    writer = csv.writer(handle, delimiter='\t')
    writer.writerow(['design_item', 'CIBERSORTx_mainline', 'Joint20_BayesPrism', 'match_status', 'evidence_or_limitation'])
    writer.writerows(rows)

# Step3. 记录审计限制和当前后台运行证据
text = f'''# Joint20 BayesPrism输入及配对设计审计

审计状态：{audit['status']}。检查时间：{snapshot_time}。

实际主线冻结发现队列为CD {audit['disease_counts']['CD']}份、UC {audit['disease_counts']['UC']}份，共{audit['samples']}份活检、{audit['patients']}位患者。metadata字段完整，SampleID无重复，集合及顺序与完整Job16 fractions、本机冻结manifest和服务器CPM表头完全一致。顺序SHA256：`{audit['ordered_sample_sha256']}`。注释文件本机/服务器字节一致。

原始2490份中排除461份Non-IBD Control，另沿用主线已冻结的4份技术排除（3 CD、1 UC），故原始CD 1157/UC 872与实际发现CD 1154/UC 871的差异得到解释。四个ID与原技术清单精确匹配，没有新增IBD过滤。保留Inflamed和Non-inflamed：{json.dumps(audit['inflammation_counts'])}。区域：{json.dumps(audit['region_counts'])}。排除明细及原技术清单见00_audit。

复用benchmark已审计reference：9294训练细胞、32供者、20类，与原selected_cell_manifest逐字段一致；现存reference counts与原始selected_training_counts在原gene symbol索引上逐元素一致。原始计数非负整数、barcode唯一、训练split未变。25703个非零symbol，无新增重复合并；同源clean_reference.rds原样复制，保持Rb/Mrp/other_Rb/chrM/MALAT1/chrX/chrY清理和exp.cells=1。训练来源包含的非IBD细胞未被删除。CIBERSORTx使用212个汇总CPM profiles构建4474×20 signature；BayesPrism使用原始9294细胞counts，因此仅训练来源/父类型/下游设计匹配。

Bulk来源固定为GSE193677_MSCCR_Biopsy_counts.txt.gz，原始文件SHA256与上次benchmark一致。以相同一对一Ensembl75/HGNC桥提取{audit['bulk_mapped_genes']}个基因，非负有限，样本库深度>0。使用完整56632基因库分母回算主线CPM，最大绝对差{audit['mainline_CPM_roundtrip_max_abs_error']:.12g}，符合原文本三位小数导出精度；现存MSCCR120 counts子集与本轮对应值逐项完全相同。最终new.prism基因数：{retained}，outlier过滤按完整队列计算，不能强行复用120例的最终基因名单。

已知局限保留：{audit['fractional_samples']}份样本包含小数，映射输入共{audit['mapped_fractional_entries']}个非整数项；未将输入舍入、改CPM或替换深度。pre-prism若内部逐项整数化，最大样本丢失比例上界为{audit['max_preprism_truncation_loss_bound_fraction']:.8g}。BayesPrism内部rmultinom的整数size限制仍存在，保留小数不等于精确抽样小数counts。实际prism保留小数数量另见prism_sample_audit.tsv。

服务器BayesPrism版本实测2.2.3，使用现有bayesprism环境，外层R seed=20261013、n.cores=8、update.gibbs=TRUE、key=NULL、outlier.cut=0.01、outlier.fraction=0.1、pseudo.min=1e-8；Gibbs及优化保留安装默认：chain.length=1000、burn.in=500、thinning=2、Gibbs内部seed=123、alpha=1，优化MAP/sigma=2。外层R seed与内部采样seed分别记录，均与benchmark一致；线程由run.prism显式n.cores=8决定。实际默认和函数源码快照保存在00_audit/02_bayesprism，推断完成后保存对象control_param。未新试marker/Ig过滤或组织分层。只读120例对象的输出契约检查与R脚本解析见logs/runtime_contract.stdout.log。

启动前资源：40 CPU、内存约512 GiB/可用511 GiB、磁盘剩余约173 GiB，无当前用户R推断/tmux长任务。按2025×16970×20双精度Z估计单数组约5.12 GiB，8 worker复制及reference/后验聚合按80 GiB级内存预算，完整对象+候选表达及中间结果按30 GiB磁盘预算；这是容量保守估算，不是实测峰值。保留8核、BLAS/OMP单线程，启动前即时资源快照及完成后time -v分别保存，未占满40核。若运行时新增共享负载，需复核。

未发现完全匹配9294细胞/20类/2025样本的既有全量结果；已有Joint20 MSCCR120及模拟benchmark不能替代，旧版多数据集2490样本结果也不匹配。本轮建立新的独立输出分支。候选结果路径清单见logs/existing_result_inventory.txt。

主线状态：完整2025 Job16 NS fractions已存在；校正及HiRes开发只覆盖120例，Job21/23/24局部HiRes及独立研究实现不能替代完整官方HiRes。完整配对HiRes、正式NMF父类型、gene filtering、rank/restarts/seed/QC/ecotype和Region3决定未冻结。BayesPrism可先运行；最终EcoTyper NMF保持HOLD。

旧版53 states/7 ecotypes、Joint20 benchmark和CIBERSORTx主线只读。所有新增分析输出、脚本与日志均在本独立目录；未更改软件源码、训练选择、正式参考或既有结果。

启动前按用户要求读取GitHub仓库main中的scripts/R/bayesprism_full/02_run_full_bayesprism.R与scripts/R/bayesprism_epithelial/16_build_core12_ecotyper_glike.R。参考其全对象检查点和Z归一化接口，但未沿用旧323144细胞参考、2490样本、两阶段/Core12选择；本轮推断参数以最近Joint20 benchmark已审计脚本为准。
'''
(out / 'JOINT20_BP_INPUT_AUDIT.md').write_text(text, encoding='utf-8')
launch = f'''# 全量MSCCR IBD BayesPrism后台启动记录

状态快照：{snapshot_time}。当前阶段：{stage}。总体退出状态：{exit_status}。

- tmux session：`{session}`；命令：`tmux new-session -d -s {session} bash {out}/scripts/run_background.sh`。
- 启动时间：{started}；launcher PID：{launcher}；正式R PID：{pid}。
- 正式命令：`{command}`。实际tmux命令及pane PID另见logs/tmux_launch_command.txt和tmux_first_status.txt。
- 证据：logs/first_runtime_check.txt（进程及子进程）、formal.stdout.log（需要出现Run Gibbs/Initial Gibbs等真实运算信息）、formal.stderr.log；不能仅以tmux存在判断已经运算。
- 使用8个BayesPrism worker，OMP/BLAS/MKL/NUMEXPR各1线程；不自动重试失败任务。
- 命令、脚本及输入SHA256见scripts、00_audit/input_SHA256.tsv、logs/script_SHA256.txt。启动前资源快照见logs/resource_prelaunch.txt。

完整结果将写02_bayesprism/BayesPrism_result.rds，保留Z及各后验。成功后后台自动执行prepare_glike_qc.R，只导出候选G-like与数值QC；marker/crosstalk审阅及正式配对HiRes等待后续完成。不会启动EcoTyper NMF。

formal_exit_status.txt和postprocessing_exit_status.txt分别记录两阶段退出码，exit_status.txt记录整体退出码，time -v保留资源用量。任一阶段失败即停止且不重启；完成后写completed_at.txt并正常退出tmux任务。没有退出码表示本快照时仍在运行，不应写成已完成。本轮仅在确认真实进入BayesPrism计算后结束前台SSH，tmux负责继续运行。

接管查看：`tmux attach -t {session}`。运行完成或失败后tmux窗口可能正常消失，最终状态以logs退出码和完成标记为准。
'''
(out / 'BAYESPRISM_FULL_IBD_LAUNCH_LOG.md').write_text(launch, encoding='utf-8')
print('REPORTS_WRITTEN', stage, 'R_PID', pid)
