# 审阅新版 Joint20 的 Job14，用冻结的 8 个等 RNA 混合物与旧 Job7 按样本、类型配对比较。
# 在本机 cibersortx_review 环境运行；读结果、日志与独立计数真值，写本地审阅目录，不访问服务器或改训练/留出集。
# 误差单位为百分点；这里是开发模拟验证，不把复用细胞的 8 个样本当成独立患者。
from pathlib import Path
import hashlib
import json
import shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pypdf import PdfReader
import pypdfium2 as pdfium

# Step1. 保存原件，核对日志和冻结的输入轴
root = Path('D:/IBD_ecotype/IBD_ecotype_analysis')
prepared = root / 'CIBERSORTx_ready_20261007/07_Joint20v2_EqualRNA8_validation/核对记录_勿上传'
out = root / 'CIBERSORTx_signature_audit/Joint20v2_Job14_review_20261007'
original = out / 'original_outputs'
qa = out / 'PDF_review'
original.mkdir(parents=True, exist_ok=True)
qa.mkdir(exist_ok=True)
sources = [Path('C:/Users/马泽凯/Downloads') / name for name in ['CIBERSORTx_Job14_Results.txt', 'CIBERSORTx_Job14_Results.pdf', 'CIBERSORTx_Job14_stackedChart.pdf']]
sources.append(Path('C:/Users/马泽凯/.codex/attachments/f9ad3969-b1cc-40fe-87f0-052bc8d901d6/已粘贴的文本.txt'))
integrity = []
for source in sources:
    name = source.name if source.name != '已粘贴的文本.txt' else 'CIBERSORTx_Job14_runtime_log.txt'
    destination = original / name
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    if destination.exists():
        assert hashlib.sha256(destination.read_bytes()).hexdigest() == digest, '不覆盖不同的网页原件'
    else:
        shutil.copyfile(source, destination)
    integrity.append({'source': str(source), 'archived_file': name, 'bytes': source.stat().st_size, 'sha256': digest})
pd.DataFrame(integrity).to_csv(out / 'original_integrity.tsv', sep='\t', index=False)
log = (original / 'CIBERSORTx_Job14_runtime_log.txt').read_text(encoding='utf-8-sig')
assert 'GJ20v2_Signature20.txt' in log and 'GJ20_Holdout8_EqualRNA_Mixture_CPM.txt' in log
assert 'perm: 100' in log and 'QN: FALSE' in log and 'All done.' in log
assert 'rmbatchSmode: TRUE' not in log and 'rmbatchBmode: TRUE' not in log
truth = pd.read_csv(prepared / 'ground_truth_equalRNA8_fractions.tsv', sep='\t', index_col=0)
before = pd.read_csv(prepared / 'Job7_same_samples_predictions.tsv', sep='\t', index_col=0)
raw = pd.read_csv(original / 'CIBERSORTx_Job14_Results.txt', sep='\t', index_col=0)
assert set(raw.index) == set(before.index) == set(truth.index)
assert set(raw.columns) == set(truth.columns) | {'P-value', 'Correlation', 'RMSE'}
before = before.loc[truth.index, truth.columns]
after = raw.loc[truth.index, truth.columns]
assert before.shape == after.shape == truth.shape == (8, 20)
assert np.isfinite(after.to_numpy()).all() and after.ge(0).all().all() and after.le(1).all().all()
assert np.allclose(after.sum(axis=1), 1)
design = pd.read_csv(prepared / 'equalRNA8_design.tsv', sep='\t').set_index('sample').loc[truth.index]
after.to_csv(out / 'Job14_predictions_aligned.tsv', sep='\t', index_label='sample')
print('Aligned 8 frozen mixtures x 20 types; Job14 completed with prescribed input names/QN/permutations', flush=True)

# Step2. 配对计算整体与分区域误差，关注每类偏差而不只看总均值
old_error = before - truth
new_error = after - truth
sample_metrics = design[['region', 'scenario', 'n_cells', 'n_donors']].copy()
sample_metrics['Job7_MAE_pp'] = old_error.abs().mean(axis=1) * 100
sample_metrics['Job14_MAE_pp'] = new_error.abs().mean(axis=1) * 100
sample_metrics['MAE_change_pp'] = sample_metrics.Job14_MAE_pp - sample_metrics.Job7_MAE_pp
sample_metrics['Job7_TV'] = old_error.abs().sum(axis=1) / 2
sample_metrics['Job14_TV'] = new_error.abs().sum(axis=1) / 2
sample_metrics.to_csv(out / 'paired_sample_metrics.tsv', sep='\t', index_label='sample')
group_rows = []
type_rows = []
for group in ['All', 'Ileum', 'ColonRectum']:
    samples = truth.index if group == 'All' else design.index[design.region.eq(group)]
    old = before.loc[samples]
    new = after.loc[samples]
    known = truth.loc[samples]
    old_abs = (old - known).abs()
    new_abs = (new - known).abs()
    group_rows.append({'group': group, 'samples': len(samples), 'Job7_MAE_pp': old_abs.to_numpy().mean() * 100, 'Job14_MAE_pp': new_abs.to_numpy().mean() * 100, 'improved_samples': int(new_abs.mean(axis=1).lt(old_abs.mean(axis=1)).sum())})
    for cell_type in truth.columns:
        present = known[cell_type].gt(0)
        type_rows.append({'group': group, 'cell_type': cell_type, 'truth_mean_percent': known[cell_type].mean() * 100, 'Job7_mean_percent': old[cell_type].mean() * 100, 'Job14_mean_percent': new[cell_type].mean() * 100, 'Job7_MAE_pp': old_abs[cell_type].mean() * 100, 'Job14_MAE_pp': new_abs[cell_type].mean() * 100, 'Job7_bias_pp': (old[cell_type] - known[cell_type]).mean() * 100, 'Job14_bias_pp': (new[cell_type] - known[cell_type]).mean() * 100, 'present_samples': int(present.sum()), 'Job7_false_zeros': int((old[cell_type].eq(0) & present).sum()), 'Job14_false_zeros': int((new[cell_type].eq(0) & present).sum()), 'Job7_present_below1pct': int((old[cell_type].lt(0.01) & present).sum()), 'Job14_present_below1pct': int((new[cell_type].lt(0.01) & present).sum())})
groups = pd.DataFrame(group_rows)
types = pd.DataFrame(type_rows)
groups.to_csv(out / 'paired_region_metrics.tsv', sep='\t', index=False)
types.to_csv(out / 'paired_type_metrics.tsv', sep='\t', index=False)
print(groups.to_string(index=False), flush=True)

# 保存全部 160 个对照，避免均值抵消高估/低估；1% 仅作低估描述，不是检测合格阈值。
comparisons = []
for sample in truth.index:
    for cell_type in truth.columns:
        comparisons.append({'sample': sample, 'region': design.at[sample, 'region'], 'scenario': design.at[sample, 'scenario'], 'cell_type': cell_type, 'truth_percent': truth.at[sample, cell_type] * 100, 'Job7_percent': before.at[sample, cell_type] * 100, 'Job14_percent': after.at[sample, cell_type] * 100, 'Job7_error_pp': old_error.at[sample, cell_type] * 100, 'Job14_error_pp': new_error.at[sample, cell_type] * 100})
comparisons = pd.DataFrame(comparisons)
comparisons.to_csv(out / 'all_160_paired_comparisons.tsv', sep='\t', index=False)
comparisons.loc[comparisons.truth_percent.eq(0)].to_csv(out / 'absent_type_predictions.tsv', sep='\t', index=False)
focus = ['Paneth', 'Goblet', 'Stem_TA_progenitor', 'Enteroendocrine', 'Tuft', 'Ileal_absorptive', 'Colonic_absorptive']
comparisons.loc[comparisons.cell_type.isin(focus)].to_csv(out / 'key_epithelial_sample_comparisons.tsv', sep='\t', index=False)

# Step3. 绘制真值、旧版和新版的区域均值；共用坐标便于查看改善与代价
fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), sharex=True)
for axis, region in zip(axes, ['Ileum', 'ColonRectum']):
    table = types.loc[types.group.eq(region)].set_index('cell_type').loc[focus]
    position = np.arange(len(focus))
    axis.barh(position - 0.24, table.truth_mean_percent, height=0.22, label='Known equal-RNA fraction', color='#333333')
    axis.barh(position, table.Job7_mean_percent, height=0.22, label='Job7: original Joint20', color='#AFB9C7')
    axis.barh(position + 0.24, table.Job14_mean_percent, height=0.22, label='Job14: donor/site Joint20 v2', color='#008A86')
    axis.set_yticks(position, [name.replace('_', ' ') for name in focus], fontsize=9)
    axis.invert_yaxis()
    axis.set_title(f'{region}: 4 matched development mixtures')
    axis.set_xlabel('Mean fraction (%)')
    axis.set_xlim(0, 35)
    axis.grid(axis='x', alpha=0.2)
    axis.set_axisbelow(True)
axes[0].legend(fontsize=8, loc='lower right')
fig.suptitle('Joint20 v2 improves ileal Paneth/Goblet; Stem/TA remains underestimated')
fig.tight_layout()
fig.savefig(out / 'Job14_vs_Job7_key_types.png', dpi=180)
plt.close(fig)

# Step4. 全页检查原 PDF，确认表格与 TXT 的显示舍入一致
pdf_texts = {}
for name in ['CIBERSORTx_Job14_Results.pdf', 'CIBERSORTx_Job14_stackedChart.pdf']:
    reader = PdfReader(original / name)
    text = '\n'.join(page.extract_text() or '' for page in reader.pages)
    pdf_texts[name] = text
    (qa / f'{Path(name).stem}_text.txt').write_text(text, encoding='utf-8')
    document = pdfium.PdfDocument(str(original / name))
    for page_number in range(len(document)):
        document[page_number].render(scale=1.6).to_pil().save(qa / f'{Path(name).stem}_page{page_number + 1}.png')
table_text = pdf_texts['CIBERSORTx_Job14_Results.pdf']
header = next(line for line in table_text.splitlines() if line.startswith('Input Sample ')).split()[2:]
pdf_rows = [line.split() for line in table_text.splitlines() if line.split() and line.split()[0] in after.index]
pdf_table = pd.DataFrame([[float(x) for x in row[1:]] for row in pdf_rows], index=[row[0] for row in pdf_rows], columns=header)
assert len(pdf_table) == 8 and set(truth.columns).issubset(pdf_table.columns)
pdf_error = (pdf_table.loc[after.index, after.columns] - after).abs().to_numpy().max()
assert pdf_error <= 0.00050001

# Step5. 保存实际验证结论与限制，不把表达拟合的 Correlation/P-value 当比例准确率
summary = {'job': 'Job14', 'baseline': 'Job7', 'status': 'review_completed', 'samples': 8, 'types': 20, 'truth': 'Frozen cell-count fraction equals equal-RNA contribution by mixture construction', 'groups': groups.to_dict(orient='records'), 'Job14_reported_zero_P_values': int(raw['P-value'].eq(0).sum()), 'Job14_expression_correlation_median': float(raw.Correlation.median()), 'PDF_TXT_max_fraction_rounding_error': float(pdf_error), 'log_complete': 'All done.' in log, 'explicit_error_or_warning_in_log': any(term in log.lower() for term in ['error:', 'warning:']), 'QN_disabled_verified': 'QN: FALSE' in log, 'permutations100_verified': 'perm: 100' in log, 'batch_correction_enable_option_seen': False, 'Job7_actual_parameters_verified': False, 'final_holdout_used': False, 'natural_mixture_validation_done_for_v2': False, 'MSCCR_run_by_this_review': False, 'HiRes_run_by_this_review': False, 'decision': 'Meaningful ileal improvement; retain v2 for natural-mixture validation, not yet a validated final reference; prioritize Stem/TA and EEC/Tuft residual errors'}
(out / 'review_summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
print('Saved:', out, flush=True)
