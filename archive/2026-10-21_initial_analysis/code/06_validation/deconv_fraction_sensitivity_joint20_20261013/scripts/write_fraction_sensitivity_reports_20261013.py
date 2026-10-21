# 汇总实际数值与失败边界；外部已知真值验证为主要证据，开发验证和真实bulk分别报告。
# 本机科学Python运行；读取本轮审计与推断表，生成中文交付报告、执行状态和文件清单。
from pathlib import Path
import hashlib
import json
import shutil
import numpy as np
import pandas as pd

# Step1. 读入实际完成的结果，不以计划代替输出
project = Path(__file__).resolve().parents[2]
out = Path((project / '06_validation/fraction_sensitivity_current_path.txt').read_text(encoding='utf-8'))
reports = out / '07_reports'
types = (out / '00_audit/frozen_types.txt').read_text().splitlines()
reference = json.loads((out / '01_reference/reference_build_audit.json').read_text())
summary = json.loads((reports / 'numeric_summary.json').read_text())
metrics = pd.read_csv(out / '04_MSCCR120_concordance/type_concordance_metrics.tsv', sep='\t')
main = metrics.loc[metrics.group.eq('All120')].copy()
external_summaries = []
independence_rows = []
blocked = []

# Step2. 每个外部数据集独立报告，避免把不同覆盖类型的MAE分母混在一起
for dataset in ['GSE214695', 'IBDverse']:
    folder = out / '08_external_pseudobulk' / dataset
    result = json.loads((folder / '04_accuracy/accuracy_summary.json').read_text())
    independence = json.loads((folder / '00_audit/source_independence_audit.json').read_text())
    generation = json.loads((folder / '02_simulation/generation_summary.json').read_text())
    settings = json.loads((folder / '00_audit/endpoint_prespecification.json').read_text())
    table = pd.read_csv(folder / '04_accuracy/accuracy_metrics_all_scales.tsv', sep='\t')
    primary = table.loc[table.truth_scale.eq('Retained_gene_RNA_PRIMARY') & table.method.eq('BayesPrism_final') & table.cell_type.isin(result['evaluable_type_names'])].copy()
    primary = primary.sort_values('MAE_pp', ascending=False)
    result_row = result['primary_result']
    cell_row = result['cell_count_result']
    lines = [f'# {dataset} 外部 pseudo-bulk benchmark', '', f"BayesPrism外部准确性计算 **COMPLETE**；官方CIBERSORTx同混合物配对 **BLOCKED**。100份混合物，每份1000个有放回细胞抽样，完整20类推断，独立作者可评价类型{result['evaluable_types']}类。", '', f"主要终点：final theta相对实际保留基因RNA-contribution truth，可评价类型MAE **{result_row['MAE_pp']:.3f}个百分点**，RMSE **{result_row['RMSE_pp']:.3f}个百分点**，pooled Pearson **{result_row['Pearson_r']:.3f}**、Spearman **{result_row['Spearman_rho']:.3f}**。", '', f"cell-count secondary：MAE **{cell_row['MAE_pp']:.3f}个百分点**，RMSE **{cell_row['RMSE_pp']:.3f}个百分点**，pooled Pearson **{cell_row['Pearson_r']:.3f}**。两种真值针对不同量，不以较好者替换预定主要终点。", '', '## 来源与独立性', '', independence['source'], '', f"全候选作者映射池{independence['eligible_cells']:,}个细胞、{independence['eligible_donors']}位供者；正式100个混合物实际抽中{generation['unique_cells']:,}个不同细胞、{generation['source_donors']}位供者。完整原始基因计数可取得；实际读取抽中的原始细胞行，非既往recovery小试验CPM。", '', 'Joint20参考只用GSE282122。两个外部来源文件/研究accession不参与当前训练，细胞标识及样本ID核对记录见来源审计。跨研究患者匿名编号没有可用身份交叉表，不能仅据编号不同宣称生物学donor-disjoint已完全证明。数据此前已用于参考小试验或recovery查看，本轮不是未触碰的前瞻性外部验收。', '', '作者label/annotation作父类型真值，没有使用EcoTyper恢复状态。歧义或不属于Joint20的细胞在生成前明确排除；这是闭集、可映射细胞组成的accuracy实验，不等价于含所有未知类型的真实组织。', '', '## 模拟方法和分母', '', '各可评价类型高斯抽样，μ为外部供者组成的等供者均值，SD=max(0.25,3×供者SD)。负值置零，归一化后多项式分配1000次抽样，类型内有放回抽原作者细胞。真值使用实际抽到的次数而非仅用目标权重。线性累加真正counts，重复gene symbols相加；BP输入未经CPM转换或人为深度舍入。', '', f"正式种子{settings['seed_formal']}，smoke种子{settings['seed_smoke']}。BayesPrism使用同一冻结参考、相同默认Gibbs参数、update.gibbs=TRUE；实际保留{result['retained_genes']:,}个基因。", '', '每个类型按实际保留基因的分子贡献计算RNA truth，同时保存cell-count truth和完整基因RNA truth。所有输出保留Full20，不删除类型、不重归一化可评价类型。跨混合物复用细胞和供者，100个模拟不是100位独立生物学供者；pooled correlation也不替代逐类型评价。', '', f"不可提供明确非零真值的类型：{', '.join(result['unsupported_types'])}。这些类型没有进入设计混合物，预测质量作为设计中的unsupported-type false-positive control；不宣称原组织中该类型不存在。其预测质量均值{result['unsupported_fraction_mean'] * 100:.3f}%，最大{result['unsupported_fraction_max'] * 100:.3f}%。", '', '## 逐类型主要结果', '', '| 类型 | Pearson | Spearman | MAE pp | RMSE pp | signed bias pp |', '|---|---:|---:|---:|---:|---:|']
    for row in primary.itertuples():
        lines.append(f'| {row.cell_type} | {row.Pearson_r:.3f} | {row.Spearman_rho:.3f} | {row.MAE_pp:.3f} | {row.RMSE_pp:.3f} | {row.mean_signed_error_pp:.3f} |')
    lines += ['', '主要误差类型：' + '；'.join(f'{r.cell_type} signed bias {r.mean_signed_error_pp:+.3f} pp、MAE {r.MAE_pp:.3f} pp' for r in primary.head(3).itertuples()) + '。总体相关不能代替这些系统偏差。', '', f"逐类型Pearson中位数{primary.Pearson_r.median():.3f}、Spearman中位数{primary.Spearman_rho.median():.3f}。脚本外层R set.seed(20261013)与Gibbs内部默认seed=123分别记录；内部chain.length=1000、burn.in=500、thinning=2、alpha=1，见installed_run_control.txt。未将外层seed误称为内部Gibbs seed。"]
    lines += ['', '## 官方CIBERSORTx配对边界', '', '用户已选择取回同批输入后自行运行官方网页，并确认先保留无batch correction主验证。目前尚未取得官方外部预测，配对比较标记BLOCKED（等待用户官方输出）；服务器没有可用Docker/token，当前可访问网页未登录。没有把研究性SVR或Job14/15旧预测冒充本批对照。精简上传目录已保存同一counts的CPM、100个SampleID、冻结Job13 signature、relative/QN关闭/no batch correction/100 permutations说明，总计48.72 MB。', '', '## 证据文件', '', '- 00_audit：作者映射、Full20覆盖、供者组成、来源独立性及预定终点。', '- 02_simulation：100000条抽样manifest、raw pseudo-bulk及两套生成真值。精确RNA真值保存后删除本轮可再生成的大细胞counts缓存；完整源文件不变，抽样manifest与pseudo-bulk counts保留。', '- 03_bayesprism：initial/final theta、CV、实际过滤基因、参数与版本；没有另存大Z数组的完整结果对象。', '- 04_accuracy：三种真值尺度×两种theta的全部数值、逐样本误差和unsupported质量。', '', f"来源：[GSE214695 GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE214695)。" if dataset == 'GSE214695' else '来源：[IBDverse原始研究](https://www.nature.com/articles/s41588-026-02634-7)、[E-MTAB-16999](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-16999)。']
    text = '\n'.join(lines) + '\n'
    (folder / f'{dataset}_EXTERNAL_BENCHMARK_REPORT.md').write_text(text, encoding='utf-8')
    (reports / f'{dataset}_EXTERNAL_BENCHMARK_REPORT.md').write_text(text, encoding='utf-8')
    external_summaries.append(result)
    independence_rows.append({'dataset': dataset, 'candidate_cells': independence['eligible_cells'], 'candidate_donors': independence['eligible_donors'], 'source_dataset_disjoint': True, 'biological_donor_linkage': 'UNKNOWN', 'evaluable_types': result['evaluable_types'], 'unsupported': '|'.join(result['unsupported_types']), 'official_pair': 'BLOCKED'})
    blocked.append({'module': dataset + '_official_CIBERSORTx_paired', 'status': 'BLOCKED', 'reason': 'Awaiting user-run official outputs; current agent runtime unavailable; same-mixture inputs retrieved and verified', 'prepared_inputs': str(out / '09_CIBERSORTx_upload_small' / (dataset + '_pseudo_bulk100_CPM.txt'))})

# Step3. 真实bulk只能报告一致性；配对患者bootstrap区间附于完整表
lines = ['# MSCCR120 Real-bulk fraction concordance', '', f"完整Joint20、同120位患者、每肠段40例，final theta为主结果。类型间Pearson中位数 **{summary['Pearson_median_across_types']:.3f}**、Spearman中位数 **{summary['Spearman_median_across_types']:.3f}**；所有样本×类型平均绝对差 **{summary['real_overall_mean_absolute_difference_pp']:.3f}个百分点**；每样本total variation中位数 **{summary['median_sample_total_variation'] * 100:.2f}%**。", '', '没有独立真实组成，以上指标不能称accuracy。Full20全部直接比较，没有重归一化；共同类型轴也是完整20类，排除类型质量为0。BP衡量实际保留基因reads贡献；CIBERSORTx为优化signature及S-mode后的relative estimates。相同训练细胞及类型体系不消除算法、基因分母、供者汇总和平台校正差别。', '', '主要结果附2000次、按肠段分层的患者配对bootstrap百分位95%区间。每位患者一例；各细胞类型的相关/误差不被当作相互独立证据。BA limits为bias±1.96×样本差异SD的描述性界限，不是准确性区间或预定临床容许范围。', '', 'CD/UC、Inflamed/Non-inflamed、Colon/Rectum/Ileum均已分层，全部类型和样本保留。高丰度定义为两方法均值的平均≥5%，低丰度<5%；接近零阈值为0.01个百分点。initial theta和Job17重复运行敏感性单独保存。', '', '| 类型 | Pearson (95% CI) | Spearman | bias pp | median abs pp | RMSE pp | BA LoA pp |', '|---|---|---:|---:|---:|---:|---|']
for row in main.itertuples():
    lines.append(f'| {row.cell_type} | {row.Pearson_r:.3f} [{row.Pearson_r_CI_low:.3f}, {row.Pearson_r_CI_high:.3f}] | {row.Spearman_rho:.3f} | {row.mean_difference_pp:.3f} | {row.median_absolute_difference_pp:.3f} | {row.RMSE_pp:.3f} | [{row.BA_lower_LoA_pp:.3f}, {row.BA_upper_LoA_pp:.3f}] |')
weak = main.sort_values('Spearman_rho').head(5)
large = main.sort_values('RMSE_pp', ascending=False).head(5)
lines += ['', '## 重点差异', '', '最低Spearman的五类：' + ', '.join(weak.cell_type) + '。', 'RMSE最高的五类：' + ', '.join(large.cell_type) + '。低相关与高绝对误差是不同问题，需同时检查丰度、零值和同类谱系互相分配。', '', '异常样本依完整Full20 total variation排序，见sample_discrepancies.tsv。逐类型数值含两方法零值/接近零频率、全部统计量及主要CI，见type_concordance_metrics.tsv。', '', '图件：20类总体及三个肠段热图、全部20类散点、Bland–Altman、均值偏差与配对绝对差异。图件用完整数据，附PDF/SVG/PNG和渲染QA。', '', 'Job20参数由用户提供的运行摘录确认：relative、S-mode、SourceSubset1060、QN关闭、100 permutations。完整网页build和上传输入哈希不可核验；Job17/20重复输出已存在实质差异，不能自动归因于BayesPrism。']
(reports / 'MSCCR120_CONCORDANCE_REPORT.md').write_text('\n'.join(lines) + '\n', encoding='utf-8')

# Step4. 参考来源与辅助开发验证，不能替代外部accuracy
capacity = pd.read_csv(out / '01_reference/reference_type_capacity.tsv', sep='\t')
reference_lines = ['# Joint20 reference provenance', '', f"参考状态 **EXACT_TRAINING_CELLS_AND_TYPES**：全部{reference['reference_cells']}训练细胞与原H5AD的counts、barcode、供者、样本和作者标签逐一核验；只使用GSE282122，32位供者，完整20类。", '', 'Job13由212个donor/site/sample/type pooled CPM profiles构建4474基因signature；本轮BP直接用对应9294个原始单细胞counts。SourceSubset1060为每个原参考profile抽5个细胞的训练子集，只用于Job20的S-mode，既不是全部训练细胞也不是212列原构建输入。此差别保留为方法限制。', '', f"参考源基因{reference['source_genes']:,}，训练细胞中非零符号基因{reference['reference_genes']:,}。MSCCR原始counts经已验证一对一Ensembl/HGNC桥保留{reference['bulk_mapped_genes']:,}个有训练表达的基因，不限制为signature4474基因。原counts中{reference['raw_fractional_entries_preserved']:,}个非整数条目原样保留，没有CPM反算或舍入。", '', '原共同训练细胞全部来自冻结training split，未纳入Job14/15或最终保留供者CD12/UC17。旧323144细胞两阶段BayesPrism仅用于清理/运行参数来源；本轮未重用旧fractions，未改标签、合并类型或两阶段解上皮。', '', '细胞、供者、区域、sample、profile、训练划分与S-mode用途见reference_provenance.tsv和signature_reference_provenance.tsv；每类数量及供者容量见reference_type_capacity.tsv。外部研究与训练source文件不同，但跨研究匿名供者缺乏患者级交叉表，生物学donor overlap不能仅靠ID命名空间完全排除。']
(reports / 'REFERENCE_PROVENANCE.md').write_text('\n'.join(reference_lines) + '\n', encoding='utf-8')
auxiliary = pd.read_csv(out / '00_audit/existing_pseudobulk_audit.tsv', sep='\t')
aux_lines = ['# Job14/15 辅助开发集 benchmark', '', '这是开发集辅助证据，不能代替两个主要外部数据集。既有CIBERSORTx结果和报告复用，没有重新运行Job14/15。', '', f"复核已有官方MAE：Job14等RNA8例 **{auxiliary.loc[auxiliary.job.eq(14), 'historical_official_MAE_pp'].iloc[0]:.3f}个百分点**；Job15自然28例 **{auxiliary.loc[auxiliary.job.eq(15), 'historical_official_MAE_pp'].iloc[0]:.3f}个百分点**。", '', '自然28例有原始cell membership、counts池、种子与cell-count/UMI真值，已精确重建并校验历史CPM；与本轮训练无细胞/供者重叠，来自同一GSE282122，已经用于signature开发评价，属于development benchmark。', '', '等RNA8例为逐细胞CPM平均，其原始counts简单求和得到的是自然28例中的配对自然样本，并非相同mixture。没有人为指定深度或四舍五入CPM来冒充raw counts，严格count-based配对标记BLOCKED。', '', '替代设计：后续可从真实UMI counts对每细胞做固定、可追踪的分子下采样，使每个入选细胞贡献相同真实分子数，再累加整数counts。它是新实验，必须记录新的抽样manifest、seed和真值，并用同一新mixture取得新的官方CIBERSORTx预测；不能复用Job14旧预测或称为历史mixture精确重建。', '', '辅助自然28例主要终点仍是历史全导出基因UMI真值；BP实际保留基因RNA truth与cell-count truth分别作尺度敏感性。这与外部主要终点的预定定义分别记录，不能混用。']
for row in summary['natural28_primary_accuracy']:
    aux_lines.append(f"- {row['method']}：MAE {row['MAE_pp']:.3f} pp，RMSE {row['RMSE_pp']:.3f} pp，Pearson {row['Pearson_r']:.3f}。")
(reports / 'PSEUDOBULK_BENCHMARK_REPORT.md').write_text('\n'.join(aux_lines) + '\n', encoding='utf-8')
blocked.append({'module': 'Job14_strict_count_based_pair', 'status': 'BLOCKED', 'reason': 'Historical mean per-cell CPM cannot be equated to raw-count summation; new controlled-downsampling design needs new official predictions'})
blocked.append({'module': 'cross_study_biological_donor_identity', 'status': 'PARTIAL', 'reason': 'Distinct source cohorts and accessions; no individual-level cross-study donor identity crosswalk'})
blocked.append({'module': 'Job20_runtime_provenance', 'status': 'PARTIAL', 'reason': 'Official fractions and user runtime excerpt available; exact web build and originally uploaded mixture/source hashes unavailable'})
pd.DataFrame(blocked).to_csv(reports / 'BLOCKED_AND_PARTIAL_ITEMS.tsv', sep='\t', index=False)
pd.DataFrame(independence_rows).to_csv(reports / 'EXTERNAL_REFERENCE_INDEPENDENCE_AND_COVERAGE.tsv', sep='\t', index=False)

# Step5. 总报告按新的优先级组织实际计算，冻结模型只给风险建议
lines = ['# IBD EcoTyper fraction sensitivity final report', '', '主要证据为两个外部作者注释数据集的已知真值accuracy。MSCCR120提供独立的真实bulk方法concordance；Job14/15仅作辅助开发证据。完整推断都使用同一9294-cell Joint20参考和冻结20标签，53 states、7 ecotypes、discovery未修改。', '', '| 实验 | 已评价混合物/样本 | 主要MAE | 评价性质 |', '|---|---:|---:|---|']
for item in external_summaries:
    lines.append(f"| {item['dataset']} | 100×1000细胞 | {item['primary_result']['MAE_pp']:.3f} pp / {item['evaluable_types']}可评价类 | BayesPrism相对保留基因RNA真值的accuracy |")
lines += [f"| MSCCR120 | 120位患者 / Full20 | {summary['real_overall_mean_absolute_difference_pp']:.3f} pp平均绝对方法差 | Real-bulk concordance，无独立真值 |", '', '外部不同数据集的类型覆盖和生物学背景不同，表中MAE分母分别为18与16类，不能据这两个数直接判断某肠段“更准确”。Full20 supplemental MAE另存，不能通过加入零真值列降低误差后替代主要终点。', '', '官方CIBERSORTx同外部mixture配对受runtime限制BLOCKED，已经提供同批CPM和参数清单。现有Job20只与相同MSCCR120配对；Job14/15只与相同开发混合物配对。未生成研究性SVR、HiRes或G-like结果。', '', '主要风险谱系建议由外部逐类型MAE、signed bias、未支持类型分配质量以及真实bulk差异共同判读。MSCCR最低相关类型：' + ', '.join(weak.cell_type) + '；最大RMSE类型：' + ', '.join(large.cell_type) + '。后续应针对这些父谱系加强独立marker、亚型覆盖及状态归属审计，不能凭fractions相近宣称G-like、NMF状态或ecotype跨方法等价。', '', '不确定性：真实bulk主分析附患者分层bootstrap95%CI；模拟的100个混合物复用来源细胞和供者，不作为100个独立生物学重复。作者类型注释是可追踪的独立参考标签，仍受作者标注、粗类型映射及歧义排除限制；跨研究患者级身份交叉表缺失，明确报告donor独立性核验边界。', '', '具体类型数值、truth定义、基因分母、initial/final sensitivity、来源哈希、完整命令和退出码见各独立报告及execution_status.json。大counts/RDS保存在单一服务器结果路径；本地取回了数值、图件、审计与日志，未覆盖历史官方结果。']
(reports / 'FINAL_SENSITIVITY_REPORT.md').write_text('\n'.join(lines) + '\n', encoding='utf-8')
with (reports / 'FINAL_SENSITIVITY_REPORT.md').open('a', encoding='utf-8') as handle:
    handle.write('\n## 实际验证暴露的系统偏差\n\nGSE214695的Tuft平均高估4.54 pp；IBDverse的Tuft平均高估17.03 pp、Stem/TA平均低估13.22 pp。真实MSCCR120也出现Tuft高于Job20、Stem/TA低于Job20的同方向差异。外部真值支持将这两类列为重点限制，但真实bulk没有组成真值，不能据此推断Job20必然正确或两者共享同一误差机制。冻结的53 states与7 ecotypes没有因这些结果调整。\n')
parity = (out / '00_audit/Cell2021_methods_parity.md').read_text(encoding='utf-8')
parity += '\n## 主验证修订\n\nGSE214695与IBDverse是主要外部pseudo-bulk。每数据集100份、每份1000细胞、有放回、Gaussian高扰动负值截零并归一化，沿用Cell2021核心流程。μ/σ改为外部作者类型的供者组成，而非癌种bulk推断；类型覆盖为18/16类，BP raw counts与官方CIBERSORTx CPM输入有尺度差异，主要endpoint为实际保留基因RNA贡献，另保存Cell2021式cell-count truth。外部官方CIBERSORTx配对BLOCKED。Job14/15已降级为开发辅助证据。\n'
(reports / 'METHODS_PARITY_AUDIT.md').write_text(parity, encoding='utf-8')

# Step6. 实际脚本和文件manifest，状态值严格采用约定枚举
for filename in ['analyze_fraction_sensitivity_joint20_20261013.py', 'plot_fraction_sensitivity_joint20_20261013.py', 'audit_external_pseudobulk_joint20_20261013.py', 'generate_external_pseudobulk_joint20_20261013.py', 'evaluate_external_pseudobulk_joint20_20261013.py', 'write_fraction_sensitivity_reports_20261013.py']:
    shutil.copyfile(project / 'scripts/python' / filename, out / 'scripts' / filename)
status = {'reference': 'COMPLETE', 'MSCCR120_BayesPrism': 'COMPLETE', 'MSCCR120_concordance': 'COMPLETE', 'GSE214695_external_BayesPrism_accuracy': 'COMPLETE', 'IBDverse_external_BayesPrism_accuracy': 'COMPLETE', 'external_official_CIBERSORTx_paired': 'BLOCKED', 'Job14_strict_count_pair': 'BLOCKED', 'Job15_auxiliary': 'COMPLETE', 'reports': 'COMPLETE', 'cross_study_donor_identity': 'PARTIAL', 'SSH_connection_count': 1, 'SSH_session_id': 69993, 'reconnected': False, 'new_CIBERSORTx_run': False, 'new_signature_training': False, 'HiRes_Glike_run': False, 'frozen_states': 53, 'frozen_ecotypes': 7, 'discovery_modified': False, 'blocked_items': blocked, 'server_root': '/home/mazekai/IBD_EcoTyper/06_validation/deconv_fraction_sensitivity_joint20_20261013', 'primary_experiments': ['GSE214695 external pseudo-bulk', 'IBDverse external pseudo-bulk'], 'notes': 'Official paired external fraction inference blocked; large source counts and full BayesPrism RDS remain on server with path manifest. Local delivery contains numerical results, scripts, provenance, logs and figures.'}
(out / 'execution_status.json').write_text(json.dumps(status, indent=2), encoding='utf-8')
manifest = []
for path in sorted(out.rglob('*')):
    if path.is_file() and path.name != 'local_output_manifest.tsv':
        with path.open('rb') as handle:
            digest = hashlib.file_digest(handle, 'sha256').hexdigest()
        manifest.append({'relative_path': str(path.relative_to(out)), 'bytes': path.stat().st_size, 'sha256': digest})
pd.DataFrame(manifest).to_csv(out / 'local_output_manifest.tsv', sep='\t', index=False)
print('REPORTS_COMPLETE', reports)
