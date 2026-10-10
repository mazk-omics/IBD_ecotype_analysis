# 审阅 Job15 新版 Joint20 在原 28 个自然混合物中的比例，与 Job6 同样本配对，并联系 Job14 等 RNA 对照。
# 在本机 cibersortx_review 环境运行；读网页原件与冻结真值，写审阅表和图，不访问服务器、不改变参考或留出。
# 自然 counts 贡献与细胞数是真值的不同口径；误差分开保存，不将模拟样本视为独立患者。
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

# Step1. 保存网页原件，核对输入与结果轴
root = Path('D:/IBD_ecotype/IBD_ecotype_analysis')
prepared = root / 'CIBERSORTx_ready_20261007/08_Joint20v2_Natural28_validation/核对记录_勿上传'
out = root / 'CIBERSORTx_signature_audit/Joint20v2_Job15_review_20261007'
original = out / 'original_outputs'
qa = out / 'PDF_review'
original.mkdir(parents=True, exist_ok=True)
qa.mkdir(exist_ok=True)
sources = [Path('C:/Users/马泽凯/Downloads') / name for name in ['CIBERSORTx_Job15_Results.txt', 'CIBERSORTx_Job15_Results.pdf', 'CIBERSORTx_Job15_stackedChart.pdf']]
sources.append(Path('C:/Users/马泽凯/.codex/attachments/a3f2270d-a5d8-4a4b-b980-3f4fcce86628/已粘贴的文本.txt'))
integrity = []
for source in sources:
    name = source.name if source.name != '已粘贴的文本.txt' else 'CIBERSORTx_Job15_runtime_log.txt'
    destination = original / name
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    if destination.exists():
        assert hashlib.sha256(destination.read_bytes()).hexdigest() == digest, '不覆盖不同内容的网页原件'
    else:
        shutil.copyfile(source, destination)
    integrity.append({'source': str(source), 'archived_file': name, 'bytes': source.stat().st_size, 'sha256': digest})
pd.DataFrame(integrity).to_csv(out / 'original_integrity.tsv', sep='\t', index=False)
log = (original / 'CIBERSORTx_Job15_runtime_log.txt').read_text(encoding='utf-8-sig')
assert 'GJ20v2_Signature20.txt' in log and 'GJ20_Holdout28_Mixture_CPM.txt' in log
assert 'perm: 100' in log and 'QN: FALSE' in log and 'All done.' in log
assert 'rmbatchSmode: TRUE' not in log and 'rmbatchBmode: TRUE' not in log
truth_rna = pd.read_csv(prepared / 'ground_truth_export_UMI_fractions.tsv', sep='\t', index_col=0)
truth_cells = pd.read_csv(prepared / 'ground_truth_cell_fractions.tsv', sep='\t', index_col=0)
before = pd.read_csv(prepared / 'Job6_same_samples_predictions.tsv', sep='\t', index_col=0)
raw = pd.read_csv(original / 'CIBERSORTx_Job15_Results.txt', sep='\t', index_col=0)
assert raw.index.is_unique and raw.columns.is_unique
assert set(raw.index) == set(before.index) == set(truth_rna.index) == set(truth_cells.index)
assert set(raw.columns) == set(truth_rna.columns) | {'P-value', 'Correlation', 'RMSE'}
before = before.loc[truth_rna.index, truth_rna.columns]
after = raw.loc[truth_rna.index, truth_rna.columns]
truth_cells = truth_cells.loc[truth_rna.index, truth_rna.columns]
assert after.shape == before.shape == truth_rna.shape == truth_cells.shape == (28, 20)
assert np.isfinite(after.to_numpy()).all() and after.ge(0).all().all() and after.le(1).all().all()
assert np.allclose(after.sum(axis=1), 1)
design = pd.read_csv(prepared / 'natural28_design.tsv', sep='\t', index_col=0).loc[truth_rna.index]
after.to_csv(out / 'Job15_predictions_aligned.tsv', sep='\t', index_label='sample')
print('Aligned 28 frozen mixtures × 20 types; completed log, QN FALSE, permutations 100', flush=True)

# Step2. 两种真值分别计算整体、区域、类型和场景误差；算法必需的分组循环保留
group_rows = []
type_rows = []
scenario_rows = []
sample_rows = []
for truth_name, truth in [('export_UMI', truth_rna), ('cell_count', truth_cells)]:
    old_error = before - truth
    new_error = after - truth
    sample_metrics = design[['region', 'scenario', 'n_cells', 'n_donors']].copy()
    sample_metrics['truth_basis'] = truth_name
    sample_metrics['Job6_MAE_pp'] = old_error.abs().mean(axis=1) * 100
    sample_metrics['Job15_MAE_pp'] = new_error.abs().mean(axis=1) * 100
    sample_metrics['MAE_change_pp'] = sample_metrics.Job15_MAE_pp - sample_metrics.Job6_MAE_pp
    sample_metrics['Job6_TV'] = old_error.abs().sum(axis=1) / 2
    sample_metrics['Job15_TV'] = new_error.abs().sum(axis=1) / 2
    sample_rows.append(sample_metrics.reset_index(names='sample'))
    for group in ['All', 'Ileum', 'ColonRectum', 'Joint']:
        samples = truth.index if group == 'All' else design.index[design.region.eq(group)]
        old = before.loc[samples]
        new = after.loc[samples]
        known = truth.loc[samples]
        old_abs = (old - known).abs()
        new_abs = (new - known).abs()
        group_rows.append({'truth_basis': truth_name, 'group': group, 'samples': len(samples), 'Job6_MAE_pp': old_abs.to_numpy().mean() * 100, 'Job15_MAE_pp': new_abs.to_numpy().mean() * 100, 'improved_samples': int(new_abs.mean(axis=1).lt(old_abs.mean(axis=1)).sum())})
        for cell_type in truth.columns:
            present = known[cell_type].gt(0)
            type_rows.append({'truth_basis': truth_name, 'group': group, 'cell_type': cell_type, 'truth_mean_percent': known[cell_type].mean() * 100, 'Job6_mean_percent': old[cell_type].mean() * 100, 'Job15_mean_percent': new[cell_type].mean() * 100, 'Job6_MAE_pp': old_abs[cell_type].mean() * 100, 'Job15_MAE_pp': new_abs[cell_type].mean() * 100, 'Job6_bias_pp': (old[cell_type] - known[cell_type]).mean() * 100, 'Job15_bias_pp': (new[cell_type] - known[cell_type]).mean() * 100, 'present_samples': int(present.sum()), 'Job6_false_zeros': int((old[cell_type].eq(0) & present).sum()), 'Job15_false_zeros': int((new[cell_type].eq(0) & present).sum())})
    for (region, scenario), block in sample_metrics.groupby(['region', 'scenario']):
        scenario_rows.append({'truth_basis': truth_name, 'region': region, 'scenario': scenario, 'samples': len(block), 'Job6_MAE_pp': block.Job6_MAE_pp.mean(), 'Job15_MAE_pp': block.Job15_MAE_pp.mean(), 'improved_samples': int(block.MAE_change_pp.lt(0).sum())})
groups = pd.DataFrame(group_rows)
types = pd.DataFrame(type_rows)
pd.concat(sample_rows, ignore_index=True).to_csv(out / 'paired_sample_metrics.tsv', sep='\t', index=False)
groups.to_csv(out / 'paired_region_metrics.tsv', sep='\t', index=False)
types.to_csv(out / 'paired_type_metrics.tsv', sep='\t', index=False)
pd.DataFrame(scenario_rows).to_csv(out / 'paired_scenario_metrics.tsv', sep='\t', index=False)
print(groups.to_string(index=False), flush=True)

# Step3. 全部逐比例对照和四个明确缺失类型，避免均值掩盖个别失败
comparisons = []
for sample in truth_rna.index:
    for cell_type in truth_rna.columns:
        comparisons.append({'sample': sample, 'region': design.at[sample, 'region'], 'scenario': design.at[sample, 'scenario'], 'cell_type': cell_type, 'RNA_truth_percent': truth_rna.at[sample, cell_type] * 100, 'cell_truth_percent': truth_cells.at[sample, cell_type] * 100, 'Job6_percent': before.at[sample, cell_type] * 100, 'Job15_percent': after.at[sample, cell_type] * 100, 'Job6_RNA_error_pp': (before.at[sample, cell_type] - truth_rna.at[sample, cell_type]) * 100, 'Job15_RNA_error_pp': (after.at[sample, cell_type] - truth_rna.at[sample, cell_type]) * 100})
comparisons = pd.DataFrame(comparisons)
comparisons.to_csv(out / 'all_560_paired_comparisons.tsv', sep='\t', index=False)
comparisons.loc[comparisons.RNA_truth_percent.eq(0)].to_csv(out / 'absent_type_predictions.tsv', sep='\t', index=False)
absent = design.loc[design.region.eq('Joint'), ['region', 'scenario', 'absent_type']].reset_index(names='sample')
absent = absent.merge(comparisons, left_on=['sample', 'absent_type'], right_on=['sample', 'cell_type'], validate='one_to_one')
assert absent.RNA_truth_percent.eq(0).all() and len(absent) == 4
absent.to_csv(out / 'four_absent_type_controls.tsv', sep='\t', index=False)
focus = ['Paneth', 'Goblet', 'Stem_TA_progenitor', 'Enteroendocrine', 'Tuft', 'Ileal_absorptive', 'Colonic_absorptive']
comparisons.loc[comparisons.cell_type.isin(focus)].to_csv(out / 'key_epithelial_sample_comparisons.tsv', sep='\t', index=False)
print(types.loc[types.truth_basis.eq('export_UMI') & types.group.isin(['Ileum', 'ColonRectum']) & types.cell_type.isin(focus)].to_string(index=False), flush=True)
print(absent[['sample', 'absent_type', 'Job6_percent', 'Job15_percent']].to_string(index=False), flush=True)

# Step4. 联系 Job14 同细胞等 RNA 对照；两种输入使用各自真值，不能直接混称比例
equal_dir = root / 'CIBERSORTx_signature_audit/Joint20v2_Job14_review_20261007'
equal_truth = pd.read_csv(root / 'CIBERSORTx_ready_20261007/07_Joint20v2_EqualRNA8_validation/核对记录_勿上传/ground_truth_equalRNA8_fractions.tsv', sep='\t', index_col=0)
equal_after = pd.read_csv(equal_dir / 'Job14_predictions_aligned.tsv', sep='\t', index_col=0).loc[equal_truth.index, equal_truth.columns]
base_names = equal_truth.index.str.removesuffix('_EqualRNA')
assert base_names.isin(after.index).all()
paired_input_rows = []
for equal_sample, natural_sample in zip(equal_truth.index, base_names):
    for cell_type in equal_truth.columns:
        paired_input_rows.append({'natural_sample': natural_sample, 'equalRNA_sample': equal_sample, 'region': design.at[natural_sample, 'region'], 'cell_type': cell_type, 'equalRNA_truth_percent': equal_truth.at[equal_sample, cell_type] * 100, 'natural_RNA_truth_percent': truth_rna.at[natural_sample, cell_type] * 100, 'Job14_equalRNA_percent': equal_after.at[equal_sample, cell_type] * 100, 'Job15_natural_percent': after.at[natural_sample, cell_type] * 100})
pd.DataFrame(paired_input_rows).to_csv(out / 'Job14_Job15_same_cell_8pair_comparisons.tsv', sep='\t', index=False)

# Step5. 科学比较图显示各区域 counts 贡献真值和两版预测
fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))
position = np.arange(len(focus))
for axis, region in zip(axes, ['Ileum', 'ColonRectum']):
    table = types.loc[types.truth_basis.eq('export_UMI') & types.group.eq(region)].set_index('cell_type').loc[focus]
    axis.barh(position - 0.24, table.truth_mean_percent, height=0.22, label='Measured counts contribution', color='#25364B')
    axis.barh(position, table.Job6_mean_percent, height=0.22, label='Job6: original Joint20', color='#AFB9C7')
    axis.barh(position + 0.24, table.Job15_mean_percent, height=0.22, label='Job15: donor/site Joint20 v2', color='#008A86')
    axis.set_yticks(position, [name.replace('_', ' ') for name in focus], fontsize=9)
    axis.invert_yaxis()
    axis.set_title(f'{region}: 12 matched natural mixtures')
    axis.set_xlabel('Mean fraction (%)')
    axis.grid(axis='x', alpha=0.2)
    axis.set_axisbelow(True)
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, fontsize=8, loc='upper center', bbox_to_anchor=(0.5, 0.935), ncol=3)
fig.suptitle('Joint20 v2: natural-mixture validation against frozen counts contributions')
fig.tight_layout(rect=(0, 0, 1, 0.87))
fig.savefig(out / 'Job15_vs_Job6_key_types.png', dpi=180)
plt.close(fig)

# Step6. 两份 PDF 全页渲染，与 TXT 检查显示舍入；循环只用于实际页数
pdf_texts = {}
page_counts = {}
for name in ['CIBERSORTx_Job15_Results.pdf', 'CIBERSORTx_Job15_stackedChart.pdf']:
    reader = PdfReader(original / name)
    text = '\n'.join(page.extract_text() or '' for page in reader.pages)
    pdf_texts[name] = text
    page_counts[name] = len(reader.pages)
    (qa / f'{Path(name).stem}_text.txt').write_text(text, encoding='utf-8')
    document = pdfium.PdfDocument(str(original / name))
    for page_number in range(len(document)):
        document[page_number].render(scale=1.6).to_pil().save(qa / f'{Path(name).stem}_page{page_number + 1}.png')
table_text = pdf_texts['CIBERSORTx_Job15_Results.pdf']
header = next(line for line in table_text.splitlines() if line.startswith('Input Sample ')).split()[2:]
pdf_rows = [line.split() for line in table_text.splitlines() if line.split() and line.split()[0] in after.index]
pdf_table = pd.DataFrame([[float(x) for x in row[1:]] for row in pdf_rows], index=[row[0] for row in pdf_rows], columns=header)
assert len(pdf_table) == 28 and set(truth_rna.columns).issubset(pdf_table.columns)
pdf_error = (pdf_table.loc[after.index, after.columns] - after).abs().to_numpy().max()
assert pdf_error <= 0.00050001
summary = {'job': 'Job15', 'baseline': 'Job6', 'status': 'numeric_review_completed_visual_review_pending', 'samples': 28, 'types': 20, 'groups': groups.to_dict(orient='records'), 'Job15_reported_zero_P_values': int(raw['P-value'].eq(0).sum()), 'Job15_expression_correlation_median': float(raw.Correlation.median()), 'Job15_expression_correlation_min': float(raw.Correlation.min()), 'PDF_pages': page_counts, 'PDF_TXT_max_fraction_rounding_error': float(pdf_error), 'log_complete': True, 'explicit_error_or_warning_in_log': any(term in log.lower() for term in ['error:', 'warning:']), 'QN_disabled_verified': True, 'permutations100_verified': True, 'batch_correction_enable_option_seen': False, 'Job6_baseline_convergence_limitation': 'Original downloaded log stops at P8%; earlier page showed max iteration warning; full28 output exists, solver convergence not independently confirmed', 'final_holdout_used': False, 'MSCCR_or_HiRes_run_by_this_review': False}
(out / 'review_summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
print('PDF page counts:', page_counts, '; maximum TXT rounding difference:', pdf_error, flush=True)
print('Saved:', out, flush=True)
