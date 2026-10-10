# 审计 Job7 等 RNA 对照，按样本和类型名称匹配已知比例，并与 Job6 的同一批细胞配对比较。
# 在本机项目 Python 环境运行；读 Downloads、冻结留出真值和 Job6 输出，写入 Job7_review。
# 原始输出保留副本；等 RNA 对照不等于自然组织 bulk，8 个混合物也不是 8 位独立供者。
import hashlib
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import pypdfium2 as pdfium
from pypdf import PdfReader

# Step1. 冻结本次结果与 PDF，防止下载目录后续替换文件
project = Path('D:/IBD_ecotype/IBD_ecotype_analysis')
base = project / 'CIBERSORTx_signature_audit/main_candidate_validation_20261006'
truth_dir = base / 'server_outputs/validation_joint20_20261006'
out = base / 'Job7_review'
originals = out / 'original_outputs'
originals.mkdir(parents=True, exist_ok=True)
file_audit = {}
for name in ['CIBERSORTx_Job7_Results.txt', 'CIBERSORTx_Job7_Results.pdf', 'CIBERSORTx_Job7_stackedChart.pdf']:
    source = Path('C:/Users/马泽凯/Downloads') / name
    target = originals / name
    source_hash = hashlib.sha256(source.read_bytes()).hexdigest()
    if target.exists():
        assert hashlib.sha256(target.read_bytes()).hexdigest() == source_hash, f'冻结副本不同：{name}'
    else:
        shutil.copy2(source, target)
    assert hashlib.sha256(target.read_bytes()).hexdigest() == source_hash
    file_audit[name] = {'bytes': source.stat().st_size, 'sha256': source_hash}

# Step2. 对照预设的 8 个样本，按名称重排预测与真值，不依赖网页列顺序
job7 = pd.read_csv(originals / 'CIBERSORTx_Job7_Results.txt', sep='\t', index_col=0)
job6 = pd.read_csv(base / 'Job6_review/original_outputs/CIBERSORTx_Job6_Results.txt', sep='\t', index_col=0)
truth = pd.read_csv(truth_dir / 'ground_truth_equalRNA_fractions.tsv', sep='\t', index_col=0)
cell_truth = pd.read_csv(truth_dir / 'ground_truth_cell_fractions.tsv', sep='\t', index_col=0)
RNA_truth = pd.read_csv(truth_dir / 'ground_truth_export_UMI_fractions.tsv', sep='\t', index_col=0)
pairs = pd.read_csv(truth_dir / 'paired_equalRNA_samples.tsv', sep='\t').set_index('sample').loc[truth.index]
design = pd.read_csv(truth_dir / 'mixture_design.tsv', sep='\t').set_index('sample').loc[pairs.paired_main_sample].copy()
design.index = truth.index
membership = pd.read_csv(truth_dir / 'mixture_cell_membership.tsv', sep='\t')
membership = membership.loc[membership['sample'].isin(pairs.paired_main_sample)]
assert job7.index.is_unique and job7.columns.is_unique and truth.index.is_unique and truth.columns.is_unique
assert pairs.index.is_unique and pairs.paired_main_sample.is_unique
assert set(job7.index) == set(truth.index) and set(job7.columns) == set(truth.columns) | {'P-value', 'Correlation', 'RMSE'}
pred7 = job7.loc[truth.index, truth.columns]
pred6 = job6.loc[pairs.paired_main_sample, truth.columns].copy()
pred6.index = truth.index
truth6_cell = cell_truth.loc[pairs.paired_main_sample, truth.columns].copy()
truth6_cell.index = truth.index
truth6_RNA = RNA_truth.loc[pairs.paired_main_sample, truth.columns].copy()
truth6_RNA.index = truth.index
assert pred7.shape == (8, 20) and np.isfinite(pred7.to_numpy()).all()
assert pred7.ge(0).all().all() and pred7.le(1).all().all() and np.allclose(pred7.sum(axis=1), 1)
assert np.allclose(truth.sum(axis=1), 1) and np.allclose(truth, truth6_cell)
assert not membership.duplicated(['sample', 'cell_id']).any()
counts = pd.crosstab(membership['sample'], membership.pilot_type).reindex(index=pairs.paired_main_sample, columns=truth.columns, fill_value=0)
recounted_truth = counts.div(counts.sum(axis=1), axis=0)
assert counts.sum(axis=1).eq(200).all() and np.allclose(recounted_truth.to_numpy(), truth.to_numpy())
pred7.to_csv(out / 'predicted_fractions_aligned.tsv', sep='\t')
print('Matched paired mixtures:', pred7.shape, 'finite, nonnegative, normalized; cell truth independently recounted', flush=True)

# Step3. 检查本地等 RNA 输入与 signature，确认真值单位和基因覆盖
equal_input = pd.read_csv(base / 'web_upload/GJ20_Holdout8_EqualRNA_Mixture_CPM.txt', sep='\t', index_col=0)
signature = pd.read_csv(base / 'web_upload/GJ20_Signature20.txt', sep='\t', index_col=0)
assert equal_input.index.is_unique and equal_input.columns.is_unique and set(equal_input.columns) == set(truth.index)
assert equal_input.shape == (23633, 8) and np.isfinite(equal_input.to_numpy()).all() and equal_input.ge(0).all().all()
assert np.allclose(equal_input.sum(axis=0), 1e6, rtol=1e-6)
assert set(signature.columns) == set(truth.columns) and signature.index.is_unique
overlap = signature.index.intersection(equal_input.index)
assert len(overlap) == len(signature.index) == 5816
input_hash = hashlib.sha256((base / 'web_upload/GJ20_Holdout8_EqualRNA_Mixture_CPM.txt').read_bytes()).hexdigest()
signature_hash = hashlib.sha256((base / 'web_upload/GJ20_Signature20.txt').read_bytes()).hexdigest()
export_summary = json.loads((truth_dir / 'summary.json').read_text(encoding='utf-8'))
assert input_hash == export_summary['files']['GJ20_Holdout8_EqualRNA_Mixture_CPM.txt']['sha256']
print('Local input verified:', equal_input.shape, 'signature overlap:', len(overlap), flush=True)

# Step4. 在同一组 8 个样本比较 MAE、总变差距离和逐类型偏差
# Job6 的自然 RNA 真值不同；同时列出细胞目标和 RNA 目标，不能把两个目标混成同一准确率。
error7 = pred7 - truth
error6_cell = pred6 - truth6_cell
error6_RNA = pred6 - truth6_RNA
sample_metrics = design[['region', 'scenario', 'n_cells', 'n_donors']].copy()
sample_metrics['paired_main_sample'] = pairs.paired_main_sample
sample_metrics['Job6_cell_MAE_pp'] = error6_cell.abs().mean(axis=1) * 100
sample_metrics['Job6_RNA_MAE_pp'] = error6_RNA.abs().mean(axis=1) * 100
sample_metrics['Job7_equalRNA_MAE_pp'] = error7.abs().mean(axis=1) * 100
sample_metrics['Job7_minus_Job6_cell_MAE_pp'] = sample_metrics.Job7_equalRNA_MAE_pp - sample_metrics.Job6_cell_MAE_pp
sample_metrics['Job6_cell_TV'] = error6_cell.abs().sum(axis=1) / 2
sample_metrics['Job6_RNA_TV'] = error6_RNA.abs().sum(axis=1) / 2
sample_metrics['Job7_equalRNA_TV'] = error7.abs().sum(axis=1) / 2
sample_metrics['Job7_fraction_pearson'] = pred7.corrwith(truth, axis=1)
sample_metrics[['Job7_web_P_value', 'Job7_web_expression_correlation', 'Job7_web_expression_RMSE']] = job7.loc[truth.index, ['P-value', 'Correlation', 'RMSE']].to_numpy()
sample_metrics.to_csv(out / 'paired_sample_accuracy_metrics.tsv', sep='\t')
metric_columns = ['Job6_cell_MAE_pp', 'Job6_RNA_MAE_pp', 'Job7_equalRNA_MAE_pp', 'Job7_minus_Job6_cell_MAE_pp', 'Job6_cell_TV', 'Job6_RNA_TV', 'Job7_equalRNA_TV']
region_metrics = sample_metrics.groupby('region')[metric_columns].mean()
region_metrics['n_pairs'] = sample_metrics.groupby('region').size()
region_metrics.to_csv(out / 'paired_region_accuracy_metrics.tsv', sep='\t')
type_metrics = pd.DataFrame({'truth_equalRNA_mean': truth.mean(), 'Job6_cell_truth_mean': truth6_cell.mean(), 'Job6_RNA_truth_mean': truth6_RNA.mean(),
    'Job6_predicted_mean': pred6.mean(), 'Job7_predicted_mean': pred7.mean(), 'Job6_cell_MAE_pp': error6_cell.abs().mean() * 100,
    'Job6_RNA_MAE_pp': error6_RNA.abs().mean() * 100, 'Job7_equalRNA_MAE_pp': error7.abs().mean() * 100,
    'Job7_equalRNA_bias_pp': error7.mean() * 100, 'Job7_zero_when_present': (pred7.eq(0) & truth.gt(0)).sum(), 'true_present_samples': truth.gt(0).sum()})
type_metrics.index.name = 'cell_type'
type_metrics.to_csv(out / 'paired_cell_type_accuracy_metrics.tsv', sep='\t')
print(region_metrics.to_string(), flush=True)

# Step5. 保存 160 个逐项对照和原本缺失类型的预测，避免均值掩盖上皮偏差
long_table = truth.rename_axis(index='sample', columns='cell_type').stack().rename('truth_equalRNA').reset_index()
for name, table in [('Job6_RNA_truth', truth6_RNA), ('Job6_predicted', pred6), ('Job7_predicted', pred7)]:
    values = table.rename_axis(index='sample', columns='cell_type').stack().rename(name).reset_index()
    long_table = long_table.merge(values, on=['sample', 'cell_type'], validate='one_to_one')
long_table = long_table.merge(sample_metrics[['region', 'scenario', 'paired_main_sample']].reset_index(), on='sample', validate='many_to_one')
long_table['Job7_bias_pp'] = (long_table.Job7_predicted - long_table.truth_equalRNA) * 100
long_table.to_csv(out / 'all_160_paired_fraction_comparisons.tsv', sep='\t', index=False)
long_table.loc[long_table.truth_equalRNA.eq(0)].to_csv(out / 'absent_type_predictions.tsv', sep='\t', index=False)
region_type = long_table.groupby(['region', 'cell_type']).agg(truth_equalRNA=('truth_equalRNA', 'mean'), Job6_RNA_truth=('Job6_RNA_truth', 'mean'), Job6_predicted=('Job6_predicted', 'mean'), Job7_predicted=('Job7_predicted', 'mean'), Job7_bias_pp=('Job7_bias_pp', 'mean'))
region_type.to_csv(out / 'paired_region_cell_type_mean_fractions.tsv', sep='\t')
epithelial = long_table.loc[long_table.cell_type.isin(['Paneth', 'Goblet', 'Stem_TA_progenitor', 'Tuft', 'Ileal_absorptive', 'Colonic_absorptive', 'Enteroendocrine'])]
epithelial.to_csv(out / 'paired_epithelial_fraction_comparisons.tsv', sep='\t', index=False)
print(region_type.loc['Ileum'].loc[['Paneth', 'Goblet', 'Stem_TA_progenitor', 'Tuft', 'Ileal_absorptive', 'Enteroendocrine']].to_string(), flush=True)

# Step6. 完整提取并渲染两份 PDF，后续以 TXT 精确值审计，PDF 只核对显示是否对应
render_dir = out / 'pdf_review'
render_dir.mkdir(exist_ok=True)
pdf_texts = {}
for name in ['CIBERSORTx_Job7_Results.pdf', 'CIBERSORTx_Job7_stackedChart.pdf']:
    pdf_path = originals / name
    reader = PdfReader(pdf_path)
    pdf_texts[name] = {'pages': len(reader.pages), 'texts': [page.extract_text() for page in reader.pages]}
    document = pdfium.PdfDocument(str(pdf_path))
    for page_number in range(len(document)):
        page = document[page_number]
        bitmap = page.render(scale=2.4)
        bitmap.to_pil().save(render_dir / f'{pdf_path.stem}_page{page_number + 1}.png')
        bitmap.close()
        page.close()
    document.close()
table_text = pdf_texts['CIBERSORTx_Job7_Results.pdf']['texts'][0]
assert all(sample in table_text for sample in truth.index)
display_rows = [line.split() for line in table_text.splitlines() if line.split() and line.split()[0] in truth.index]
display_table = pd.DataFrame([list(map(float, row[1:])) for row in display_rows], index=[row[0] for row in display_rows], columns=list(job7.columns[:-1]))
display_table = display_table.loc[truth.index]
assert np.allclose(display_table, job7.loc[truth.index, display_table.columns], atol=0.00051)
(render_dir / 'pdf_text_extraction.json').write_text(json.dumps(pdf_texts, ensure_ascii=False, indent=2), encoding='utf-8')

# Step7. 写入审计摘要；用户未提供 Job7 日志，不能独立确认网页参数和求解器收敛
ileum = design.region.eq('Ileum')
summary = {'status': 'numeric_audit_complete_pdf_visual_review_pending', 'job': 'Job7', 'samples': 8, 'cell_types': 20, 'all_160_pairs_matched': True,
    'all_fractions_finite_nonnegative_and_sum_to_one': True, 'cell_truth_recounted_from_membership': True, 'equalRNA_truth_equals_paired_cell_truth': True,
    'Job6_paired_cell_MAE_pp': float(error6_cell.abs().to_numpy().mean() * 100), 'Job6_paired_RNA_MAE_pp': float(error6_RNA.abs().to_numpy().mean() * 100),
    'Job7_equalRNA_MAE_pp': float(error7.abs().to_numpy().mean() * 100), 'region_metrics': region_metrics.to_dict(orient='index'),
    'ileum_mean_Paneth_truth_percent': float(truth.loc[ileum, 'Paneth'].mean() * 100),
    'ileum_mean_Job6_Paneth_percent': float(pred6.loc[ileum, 'Paneth'].mean() * 100), 'ileum_mean_Job7_Paneth_percent': float(pred7.loc[ileum, 'Paneth'].mean() * 100),
    'ileum_Goblet_zero_when_present': int((pred7.loc[ileum, 'Goblet'].eq(0) & truth.loc[ileum, 'Goblet'].gt(0)).sum()),
    'ileum_Stem_TA_zero_when_present': int((pred7.loc[ileum, 'Stem_TA_progenitor'].eq(0) & truth.loc[ileum, 'Stem_TA_progenitor'].gt(0)).sum()),
    'web_expression_correlation_median': float(job7.Correlation.median()), 'web_P_values_reported_zero': int(job7['P-value'].eq(0).sum()),
    'run_settings_independently_confirmed': False, 'solver_convergence_confirmed': False, 'Job7_runtime_log_provided': False,
    'local_equalRNA_input_sha256': input_hash, 'local_signature_sha256': signature_hash, 'original_outputs': file_audit,
    'simulation_samples_are_independent_donors': False, 'natural_bulk_accuracy_validated': False,
    'decision': 'Equal-RNA control improves paired cell-fraction error but retains major ileal epithelial bias; current Joint20 remains a baseline, not an accepted validated whole-cohort primary reference.',
    'decision_is_project_judgment_not_predefined_statistical_pass_threshold': True, 'MSCCR_deconvolution_run': False, 'new_server_work_performed': False}
visual_path = render_dir / 'pdf_visual_review.json'
if visual_path.exists():
    visual_review = json.loads(visual_path.read_text(encoding='utf-8'))
    assert visual_review['source_sha256'] == {name: info['sha256'] for name, info in file_audit.items() if name.endswith('.pdf')}
    summary['status'] = 'review_completed'
    summary['PDF_visual_review'] = visual_review
summary['PDF_table_matches_TXT_at_display_precision'] = True
(out / 'review_summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')

# Step8. 保存中文审阅结论，让读者直接看到配对关系、单位、主要偏差和适用边界
report = [
    '# Job7 等 RNA 配对对照审阅（2026-10-06）', '',
    '结论：等 RNA 后整体细胞比例误差下降，但回肠 Paneth 仍明显高估，Goblet 与 Stem/TA 仍严重低估。当前 Joint20 保留为对照基线，暂不接受为已经验证的 MSCCR 全队列主参考。这里的接受判断是项目判断，未预设统计学通过阈值。', '',
    '## 1. 本次核对了什么', '',
    '- Job7 有 8 个样本、20 个类型，160 项预测均按名称匹配；比例有限、非负且每个样本合计为 1。',
    '- 8 个样本对应预先固定的 Job6 配对名单：回肠与结直肠各 balanced 1/2、epithelial_high 1、immune_high 1。',
    '- 每个混合物仍用原来的 200 个细胞。Job7 输入先将每个细胞转成 CPM，再平均，使每个细胞贡献相同 RNA 总量；已知 RNA 贡献因此等于实际细胞数比例。真值直接从 cell membership 重新计数，与保存真值一致，未使用 signature 或预测比例生成。',
    '- 配对样本共涉及 1,423 个不同细胞，混合物之间存在复用；8 个混合物不能视为 8 位独立患者，不做独立患者显著性检验。训练供者隔离和最终留出供者未使用的结论沿用此前完成的输入审计。',
    '- 本地输入为 23,633 基因 × 8 个样本，CPM 列和为 1,000,000；signature 的 5,816 个基因全部覆盖，输入 hash 与服务器导出记录一致。',
    '- 两份 PDF 的完整页面已渲染供视觉检查；表格 PDF 的 8 行、20 类比例与 TXT 在三位小数显示精度内一致。网页 PDF 未显示 RMSE，数值审计以 TXT 为准。',
    '- 此次未提供 Job7 runtime log 或参数截图，不能独立确认实际使用的 signature 文件、batch correction、QN、置换次数及求解器收敛。下述机制比较以按既定指引保持其他参数为前提。', '',
    '## 2. 同一批样本的误差比较', '',
    'MAE 为每个样本 20 类比例绝对误差的平均，单位是百分点。Job6 有细胞比例和自然 RNA 贡献两个不同目标；Job7 的等 RNA 目标与细胞比例相同。不能用 Job6 全部 28 个样本的指标与 Job7 的 8 个样本直接比较。', '',
    '| 配对范围 | 样本数 | Job6 对细胞比例 MAE | Job6 对自然 RNA 比例 MAE | Job7 对等 RNA 比例 MAE |',
    '|---|---:|---:|---:|---:|'
]
for region, row in region_metrics.iterrows():
    label = {'Ileum': '回肠', 'ColonRectum': '结直肠'}[region]
    report.append(f'| {label} | {int(row.n_pairs)} | {row.Job6_cell_MAE_pp:.2f} | {row.Job6_RNA_MAE_pp:.2f} | {row.Job7_equalRNA_MAE_pp:.2f} |')
report.append(f'| 合计 | 8 | {summary["Job6_paired_cell_MAE_pp"]:.2f} | {summary["Job6_paired_RNA_MAE_pp"]:.2f} | {summary["Job7_equalRNA_MAE_pp"]:.2f} |')
report += ['', '8 个样本的细胞目标 MAE 都下降，但不代表每个类型都改善。总变差距离 TV 能反映整组成分差异：回肠 Job6 细胞目标 TV 为 0.565，Job7 为 0.324；结直肠分别为 0.335 和 0.144。平均到 20 列后的 MAE 较小，仍可能掩盖集中在少数上皮类型上的大偏差。', '',
    '## 3. 回肠偏差仍然保留', '',
    '| 回肠配对场景 | 等 RNA / 细胞真值 Paneth | Job6 自然 RNA 真值 Paneth | Job6 预测 | Job7 预测 | Job7 偏差（百分点） |',
    '|---|---:|---:|---:|---:|---:|']
for row in long_table.loc[long_table.region.eq('Ileum') & long_table.cell_type.eq('Paneth')].itertuples():
    report.append(f'| {row.paired_main_sample} | {row.truth_equalRNA * 100:.2f}% | {row.Job6_RNA_truth * 100:.2f}% | {row.Job6_predicted * 100:.2f}% | {row.Job7_predicted * 100:.2f}% | +{row.Job7_bias_pp:.2f} |')
report += ['', '注意：表中 4 个样本均来自同一组配对细胞；Job6 与 Job7 的 RNA 真值不同。Job6 预测值减小本身不等于校准通过，必须和各自真值比较。', '',
    '| 类型 | 回肠等 RNA 真值均值 | Job7 预测均值 | 平均偏差（百分点） |', '|---|---:|---:|---:|']
for name in ['Paneth', 'Goblet', 'Stem_TA_progenitor', 'Tuft', 'Ileal_absorptive', 'Enteroendocrine']:
    row = region_type.loc[('Ileum', name)]
    report.append(f'| {name} | {row.truth_equalRNA * 100:.2f}% | {row.Job7_predicted * 100:.2f}% | {row.Job7_bias_pp:+.2f} |')
report += ['',
    'Goblet 和 Stem/TA 在 4 个回肠样本中均实际存在，但前三个样本都被预测为零（各 3/4）。Tuft 从 Job6 四个配对样本均为零改善为四个都非零，但仍有低估。等 RNA 处理没有同时恢复所有上皮类型。比例受总和为 1 的约束，不能据此把每一类低估都认定为直接错分给 Paneth。', '',
    '## 4. 结直肠改善较明显，仍有类型偏差', '',
    '4 个结直肠样本的等 RNA MAE 为 1.44 个百分点。Goblet 真值均值 6.88%、预测 8.48%；Stem/TA 真值 7.63%、预测 9.91%；Colonic_absorptive 真值 7.13%、预测 4.50%；Tuft 真值 8.13%、预测 6.11%。这些描述来自 4 个混合物，不能视为已经通过独立结直肠队列验证。', '',
    '基线设计有意不放入结直肠 Paneth 和 Ileal_absorptive；前者预测最大仅 0.0133%，后者最大 0.0229%；回肠本来缺失的 Colonic_absorptive 均预测为零。这支持本组场景中的缺失类型识别，但不能推出任意背景下的假阳性率。', '',
    '## 5. 结果支持什么，不支持什么', '',
    '- 支持：自然 RNA 产量差异会影响细胞比例误差，等 RNA 控制使多个类型和整体误差改善。',
    '- 支持：RNA 产量不同不能完全解释回肠问题；在等 RNA 输入下仍有明显的 Paneth 高估与 Goblet、Stem/TA 低估。',
    '- 待检验：参考表达幅度、回肠/非回肠 Paneth 混合、部位构成、marker 共线性及算法求解/归一化适配各自贡献多少。此前训练 Paneth 的 8 个回肠、42 个非回肠组成和 marker 差异是候选线索，这次对照不能单独证明因果。',
    '- 不支持：据此认定 EcoTyper 或 CIBERSORTx 全面不适合肠道；这次仍是对某一参考和模拟输入的审计，尚未测试 downstream cell-state expression 的准确性。',
    '- 不支持：把等 RNA 的改善当作自然 bulk 细胞比例已经准确。真实组织并不会让每个细胞贡献相同 RNA，跨平台和真实组织偏差仍需验证。', '',
    f'Job7 网页表达拟合 Correlation 中位数为 {summary["web_expression_correlation_median"]:.3f}，8 个 P-value 均显示为 0。前者是表达拟合指标，后者是有限置换输出，均不等于细胞比例准确率。', '',
    '## 6. 本日收尾与后续工作顺序', '',
    '今天完成 Job7 数值审计、原始文件冻结、PDF 检查和 Job6 配对比较；未访问服务器、未提交新网页任务、未运行 MSCCR 或 recovery。今天在这里结束。', '',
    '后续第一优先级建议使用现有 Ileum19 signature，在相同的回肠留出细胞和样本上做对照，同时保留自然 counts 与等 RNA 两种输入。它可以判断回肠参考是否改善主要偏差，但两份训练参考的供者和细胞也不同，因此不是只改变部位的纯因果实验。随后再决定 Joint20 v2 的部位/供者抽样、Paneth 类型层级和表达标度；保留最终留出的 CD12、UC17 供者用于冻结方案后的检验。', '',
    'MSCCR CD+UC 2,025 样本的网页输入已经完成格式与基因覆盖核对，可继续保留；当前 Joint20 的 MSCCR 运行若开展，应标注为探索性基线，不能以输入合规替代参考准确性验证。', '',
    '详细数据：paired_sample_accuracy_metrics.tsv、paired_region_accuracy_metrics.tsv、paired_cell_type_accuracy_metrics.tsv、all_160_paired_fraction_comparisons.tsv、paired_epithelial_fraction_comparisons.tsv、absent_type_predictions.tsv；原始 TXT/PDF 保存在 original_outputs。'
]
(out / 'Job7等RNA配对对照审阅.md').write_text('\n'.join(report) + '\n', encoding='utf-8')
print(json.dumps(summary, ensure_ascii=True, indent=2), flush=True)
