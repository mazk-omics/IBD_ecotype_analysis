# 审阅Job9回肠自然混合物，重点与Job8同细胞等RNA对照比较，Job6作辅助基线。
# 在本机.tools/cibersortx_review环境运行；读Downloads及冻结真值，写Job9_review。
# 主真值为原始counts贡献比例；不改输入、重抽供者或把12个混合物当独立患者。
from pathlib import Path
import hashlib
import json
import shutil
import importlib.metadata
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pypdf import PdfReader
import pypdfium2 as pdfium

# Step1. 保存原结果和日志，核对样本、19类轴及本地提交文件完整性。
project = Path('D:/IBD_ecotype/IBD_ecotype_analysis')
ready = project / 'CIBERSORTx_ready_20261007'
metadata = ready / '核对记录_勿上传'
out = project / 'CIBERSORTx_signature_audit/ileum_candidate_validation_20261007/Job9_review'
originals = out / 'original_outputs'
originals.mkdir(parents=True, exist_ok=True)
sources = [Path('C:/Users/马泽凯/Downloads') / name for name in ['CIBERSORTx_Job9_Results.txt', 'CIBERSORTx_Job9_Results.pdf', 'CIBERSORTx_Job9_stackedChart.pdf']]
sources.append(Path('C:/Users/马泽凯/.codex/attachments/666d322f-de0e-4fa9-8919-b4f09f0d32fe/已粘贴的文本.txt'))
file_audit = []
for source in sources:
    name = 'CIBERSORTx_Job9_runtime_log.txt' if source.name == '已粘贴的文本.txt' else source.name
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    target = originals / name
    if target.exists():
        assert hashlib.sha256(target.read_bytes()).hexdigest() == digest, '已有冻结副本不同，不覆盖。'
    else:
        shutil.copy2(source, target)
    file_audit.append({'name': name, 'source': str(source), 'sha256': digest, 'bytes': source.stat().st_size})
result = pd.read_csv(originals / 'CIBERSORTx_Job9_Results.txt', sep='\t', index_col=0)
truth = pd.read_csv(metadata / 'GI19_Natural12_ground_truth_RNA.tsv', sep='\t', index_col=0)
cell_truth = pd.read_csv(metadata / 'GI19_Natural12_ground_truth_cell.tsv', sep='\t', index_col=0)
joint20 = pd.read_csv(metadata / 'GI19_Natural12_Joint20_baseline_predictions20.tsv', sep='\t', index_col=0).loc[truth.index]
design = pd.read_csv(metadata / 'GI19_natural12_sample_design.tsv', sep='\t').set_index('sample').loc[truth.index]
membership = pd.read_csv(metadata / 'GI19_natural12_existing_cell_membership.tsv', sep='\t')
signature = pd.read_csv(ready / '01_Ileum19_validation/GI19_Signature19.txt', sep='\t', index_col=0)
mixture = pd.read_csv(ready / '01_Ileum19_validation/GI19_HoldoutNatural12_Mixture_CPM.txt', sep='\t', index_col=0)
assert result.index.is_unique and result.columns.is_unique and set(result.index) == set(truth.index)
assert set(result.columns) == set(truth.columns) | {'P-value', 'Correlation', 'RMSE'}
assert signature.columns.equals(truth.columns) and mixture.columns.equals(truth.index)
prediction = result.loc[truth.index, truth.columns]
baseline = joint20.loc[truth.index, truth.columns]
assert prediction.shape == (12, 19) and np.isfinite(prediction.to_numpy()).all() and prediction.ge(0).all().all()
assert np.allclose(prediction.sum(axis=1), 1) and np.allclose(truth.sum(axis=1), 1) and np.allclose(cell_truth.sum(axis=1), 1)
integrity = pd.read_csv(metadata / 'upload_file_integrity.tsv', sep='\t').set_index('relative_path')
for name in ['GI19_Signature19.txt', 'GI19_HoldoutNatural12_Mixture_CPM.txt']:
    path = ready / '01_Ileum19_validation' / name
    assert hashlib.sha256(path.read_bytes()).hexdigest() == integrity.loc[str(path.relative_to(ready)), 'sha256']
log = (originals / 'CIBERSORTx_Job9_runtime_log.txt').read_text(encoding='utf-8')
assert 'GI19_HoldoutNatural12_Mixture_CPM.txt' in log and 'GI19_Signature19.txt' in log
assert '[Options] QN: FALSE' in log and '[Options] perm: 100' in log and 'All done.' in log

# Step2. 从逐细胞组成重算细胞比例和UMI贡献，保留两种单位的区别。
assert not membership.duplicated(['sample', 'cell_id']).any()
counts = pd.crosstab(membership['sample'], membership.pilot_type).reindex(index=truth.index, columns=truth.columns, fill_value=0)
RNA_counts = membership.groupby(['sample', 'pilot_type']).export_library_counts.sum().unstack(fill_value=0).reindex(index=truth.index, columns=truth.columns, fill_value=0)
assert np.allclose(counts.div(counts.sum(axis=1), axis=0), cell_truth)
assert np.allclose(RNA_counts.div(RNA_counts.sum(axis=1), axis=0), truth)
assert counts.sum(axis=1).eq(200).all() and not set(membership.Patient) & {'CD12', 'UC17'}
assert mixture.index.is_unique and signature.index.is_unique and np.allclose(mixture.sum(axis=0), 1e6)
overlap = signature.index.intersection(mixture.index)
prediction.to_csv(out / 'predicted_fractions_aligned.tsv', sep='\t')
truth.to_csv(out / 'ground_truth_RNA19.tsv', sep='\t')
cell_truth.to_csv(out / 'ground_truth_cell19.tsv', sep='\t')
print('Matched:', prediction.shape, 'signature coverage:', len(overlap), '/', len(signature), 'donors:', sorted(membership.Patient.unique()), flush=True)

# Step3. 固定19类主轴比较，不重归一化；完整20类轴另保留额外类型的预测质量。
error = prediction - truth
baseline_error = baseline - truth
truth20 = truth.reindex(columns=joint20.columns, fill_value=0)
prediction20 = prediction.reindex(columns=joint20.columns, fill_value=0)
sample_metrics = design.copy()
sample_metrics['Job9_RNA_MAE19_pp'] = error.abs().mean(axis=1) * 100
sample_metrics['Job6_RNA_MAE19_pp'] = baseline_error.abs().mean(axis=1) * 100
sample_metrics['Job9_RNA_TV20'] = (prediction20 - truth20).abs().sum(axis=1) / 2
sample_metrics['Job6_RNA_TV20'] = (joint20 - truth20).abs().sum(axis=1) / 2
sample_metrics['Job9_cell_MAE19_pp'] = (prediction - cell_truth).abs().mean(axis=1) * 100
sample_metrics['Job6_cell_MAE19_pp'] = (baseline - cell_truth).abs().mean(axis=1) * 100
sample_metrics['Job9_fraction_Pearson_RNA'] = prediction.corrwith(truth, axis=1)
sample_metrics['Job9_web_expression_Correlation'] = result.loc[truth.index, 'Correlation']
sample_metrics['Job9_web_RMSE'] = result.loc[truth.index, 'RMSE']
sample_metrics.to_csv(out / 'paired_sample_accuracy_metrics.tsv', sep='\t')
type_metrics = pd.DataFrame({'RNA_truth_mean_percent': truth.mean() * 100, 'cell_truth_mean_percent': cell_truth.mean() * 100, 'Job9_mean_percent': prediction.mean() * 100, 'Job6_mean_percent': baseline.mean() * 100,
    'Job9_RNA_MAE_pp': error.abs().mean() * 100, 'Job6_RNA_MAE_pp': baseline_error.abs().mean() * 100, 'Job9_bias_pp': error.mean() * 100,
    'Job9_zero_when_present': (prediction.eq(0) & truth.gt(0)).sum(), 'true_present_samples': truth.gt(0).sum()})
type_metrics.index.name = 'cell_type'
type_metrics.to_csv(out / 'paired_cell_type_accuracy_metrics.tsv', sep='\t')
scenario_metrics = sample_metrics.groupby('scenario')[['Job9_RNA_MAE19_pp', 'Job6_RNA_MAE19_pp', 'Job9_RNA_TV20', 'Job6_RNA_TV20']].mean()
scenario_metrics.to_csv(out / 'scenario_accuracy_metrics.tsv', sep='\t')
long = truth.rename_axis(index='sample', columns='cell_type').stack().rename('RNA_truth').reset_index()
for name, frame in [('cell_truth', cell_truth), ('Job9', prediction), ('Job6', baseline)]:
    values = frame.rename_axis(index='sample', columns='cell_type').stack().rename(name).reset_index()
    long = long.merge(values, on=['sample', 'cell_type'], validate='one_to_one')
long['Job9_bias_pp'] = (long.Job9 - long.RNA_truth) * 100
long.to_csv(out / 'all_228_paired_fraction_comparisons.tsv', sep='\t', index=False)
print(scenario_metrics.round(3).to_string(), flush=True)
print(type_metrics.round(3).to_string(), flush=True)

# Step4. 比较相同细胞的自然与等RNA构造，各自使用对应真值；不能跨构造直接比较预测比例。
equal_dir = out.parent / 'Job8_review'
equal_truth = pd.read_csv(metadata / 'GI19_EqualRNA4_ground_truth_RNA.tsv', sep='\t', index_col=0)
equal_pred = pd.read_csv(equal_dir / 'predicted_fractions_aligned.tsv', sep='\t', index_col=0).loc[equal_truth.index, truth.columns]
pairs = pd.read_csv(metadata / 'GI19_equalRNA4_pairing.tsv', sep='\t').set_index('sample').loc[equal_truth.index]
equal_joint20 = pd.read_csv(metadata / 'GI19_EqualRNA4_Joint20_baseline_predictions20.tsv', sep='\t', index_col=0).loc[equal_truth.index, truth.columns]
equal_error = equal_pred - equal_truth
assert np.allclose(cell_truth.loc[pairs.paired_main_sample].to_numpy(), equal_truth.to_numpy())
paired_sample_metrics = pairs.copy()
paired_sample_metrics['Job9_natural_MAE19_pp'] = sample_metrics.loc[pairs.paired_main_sample, 'Job9_RNA_MAE19_pp'].to_numpy()
paired_sample_metrics['Job8_equalRNA_MAE19_pp'] = equal_error.abs().mean(axis=1) * 100
paired_sample_metrics['Job9_minus_Job8_MAE_pp'] = paired_sample_metrics.Job9_natural_MAE19_pp - paired_sample_metrics.Job8_equalRNA_MAE19_pp
paired_sample_metrics.to_csv(out / 'same_cells_natural_vs_equalRNA_sample_metrics.tsv', sep='\t')
paired_rows = []
for equal_sample, row in pairs.iterrows():
    natural_sample = row.paired_main_sample
    for cell_type in truth.columns:
        paired_rows.append({'natural_sample': natural_sample, 'equalRNA_sample': equal_sample, 'cell_type': cell_type,
            'natural_RNA_truth': truth.loc[natural_sample, cell_type], 'Job9_natural_prediction': prediction.loc[natural_sample, cell_type], 'Job9_natural_bias_pp': error.loc[natural_sample, cell_type] * 100,
            'equalRNA_truth': equal_truth.loc[equal_sample, cell_type], 'Job8_equalRNA_prediction': equal_pred.loc[equal_sample, cell_type], 'Job8_equalRNA_bias_pp': equal_error.loc[equal_sample, cell_type] * 100})
paired_constructs = pd.DataFrame(paired_rows)
paired_constructs.to_csv(out / 'same_cells_natural_vs_equalRNA_76_pairs.tsv', sep='\t', index=False)
construct_type_summary = paired_constructs.groupby('cell_type')[['natural_RNA_truth', 'Job9_natural_prediction', 'Job9_natural_bias_pp', 'equalRNA_truth', 'Job8_equalRNA_prediction', 'Job8_equalRNA_bias_pp']].mean()
construct_type_summary['Job9_natural_MAE_pp'] = paired_constructs.Job9_natural_bias_pp.abs().groupby(paired_constructs.cell_type).mean()
construct_type_summary['Job8_equalRNA_MAE_pp'] = paired_constructs.Job8_equalRNA_bias_pp.abs().groupby(paired_constructs.cell_type).mean()
construct_type_summary.to_csv(out / 'paired_construct_cell_type_means.tsv', sep='\t')
natural_four_MAE = error.loc[pairs.paired_main_sample].abs().to_numpy().mean() * 100
equal_four_MAE = equal_error.abs().to_numpy().mean() * 100

# Step5. 核查粗谱系和实际亚型覆盖，不能将回肠Stem/TA误差套用SCP259的抽样解释。
epithelial_types = ['Ileal_absorptive', 'Paneth', 'Goblet', 'Stem_TA_progenitor', 'Enteroendocrine', 'Tuft']
epithelium = pd.DataFrame({'RNA_truth': truth[epithelial_types].sum(axis=1), 'Job9': prediction[epithelial_types].sum(axis=1), 'Job6': baseline[epithelial_types].sum(axis=1), 'scenario': design.scenario})
epithelium['Job9_abs_error_pp'] = (epithelium.Job9 - epithelium.RNA_truth).abs() * 100
epithelium['Job6_abs_error_pp'] = (epithelium.Job6 - epithelium.RNA_truth).abs() * 100
epithelium.to_csv(out / 'epithelial_total_fractions.tsv', sep='\t')
epithelium.groupby('scenario').mean(numeric_only=True).to_csv(out / 'epithelial_scenario_means.tsv', sep='\t')
candidate = project / 'CIBERSORTx_signature_audit/three_candidate_comparison_20261006/gse_ileum'
markers = pd.read_csv(candidate / 'marker_panel_audit.tsv', sep='\t')
markers = markers.loc[markers.pilot_type.isin(epithelial_types)].copy()
markers['in_test_mixture'] = markers.gene.isin(mixture.index)
markers.to_csv(out / 'epithelial_marker_coverage_and_specificity.tsv', sep='\t', index=False)
train_dir = project / 'IBDverse_reference/three_reference_pilots_20261006/run_20261006_183211/gse_ileum'
training = pd.read_csv(train_dir / 'selected_cell_manifest.tsv', sep='\t')
assert not set(training.Patient) & set(membership.Patient)
training.groupby(['pilot_type', 'author_label']).size().rename('cells').to_csv(out / 'training_author_subtype_counts.tsv', sep='\t')
pool = pd.read_csv(project / 'CIBERSORTx_signature_audit/main_candidate_validation_20261006/server_outputs/validation_joint20_20261006/selected_pool_manifest.tsv', sep='\t', usecols=['cell_id', 'final_analysis'])
annotated = membership.merge(pool, on='cell_id', validate='many_to_one')
assert annotated.final_analysis.notna().all()
annotated.groupby(['sample', 'pilot_type', 'final_analysis']).size().rename('cells').to_csv(out / 'test_author_subtype_counts.tsv', sep='\t')
pilot_summary = json.loads((train_dir / 'summary.json').read_text(encoding='utf-8'))

# Step6. 提取完整PDF并渲染所有页，核对显示值；视觉结论在看图后单独记录。
render_dir = out / 'pdf_review'
render_dir.mkdir(exist_ok=True)
pdf_text = {}
for name in ['CIBERSORTx_Job9_Results.pdf', 'CIBERSORTx_Job9_stackedChart.pdf']:
    reader = PdfReader(originals / name)
    pdf_text[name] = {'pages': len(reader.pages), 'text': [page.extract_text() for page in reader.pages]}
    document = pdfium.PdfDocument(str(originals / name))
    for page_number in range(len(document)):
        page = document[page_number]
        bitmap = page.render(scale=2.4)
        bitmap.to_pil().save(render_dir / f'{Path(name).stem}_page{page_number + 1}.png')
        bitmap.close()
        page.close()
    document.close()
(render_dir / 'pdf_text_extraction.json').write_text(json.dumps(pdf_text, ensure_ascii=False, indent=2), encoding='utf-8')
display_lines = [line.split() for line in pdf_text['CIBERSORTx_Job9_Results.pdf']['text'][0].splitlines() if line.split() and line.split()[0] in truth.index]
display = pd.DataFrame([list(map(float, row[1:])) for row in display_lines], index=[row[0] for row in display_lines], columns=result.columns[:-1]).loc[truth.index]
assert np.allclose(display, result.loc[truth.index, display.columns], atol=0.00051)

# Step7. 统一类型和样本轴画图，避免用较高网页拟合相关性掩盖组成误差。
labels = ['B1', 'B2', 'B3', 'B4', 'E1', 'E2', 'E3', 'E4', 'I1', 'I2', 'I3', 'I4']
fig, axes = plt.subplots(1, 3, figsize=(17, 7), sharey=True)
for ax, frame, title in zip(axes, [truth, baseline, prediction], ['Known RNA contribution', 'Job6 Joint20 (same 19 types)', 'Job9 Ileum19']):
    im = ax.imshow(frame.T.to_numpy() * 100, aspect='auto', vmin=0, vmax=80, cmap='viridis')
    ax.set_title(title)
    ax.set_xticks(range(12), labels, rotation=45)
    ax.set_yticks(range(19), truth.columns)
    ax.set_xlabel('Balanced / epithelial-high / immune-high')
fig.subplots_adjust(left=0.15, right=0.87, bottom=0.14, wspace=0.08)
color_axis = fig.add_axes([0.90, 0.22, 0.015, 0.62])
fig.colorbar(im, cax=color_axis, label='RNA fraction (%)')
fig.savefig(out / 'truth_vs_Job6_vs_Job9.png', dpi=180, bbox_inches='tight')
plt.close(fig)
fig, ax = plt.subplots(figsize=(9, 7))
im = ax.imshow(error.T.to_numpy() * 100, aspect='auto', vmin=-50, vmax=50, cmap='RdBu_r')
ax.set_xticks(range(12), labels)
ax.set_yticks(range(19), truth.columns)
ax.set_title('Job9 minus known RNA contribution')
fig.colorbar(im, ax=ax, label='Bias (percentage points)')
fig.tight_layout()
fig.savefig(out / 'Job9_fraction_bias.png', dpi=180)
plt.close(fig)

# 相同四批细胞分列两种真值，直接展示改变RNA权重之后哪些误差持续存在。
natural_four_truth = truth.loc[pairs.paired_main_sample]
natural_four_prediction = prediction.loc[pairs.paired_main_sample]
fig, axes = plt.subplots(1, 4, figsize=(16, 7), sharey=True)
for ax, frame, title in zip(axes, [equal_truth, equal_pred, natural_four_truth, natural_four_prediction], ['Equal RNA truth', 'Job8 equal RNA', 'Natural RNA truth', 'Job9 natural RNA']):
    im = ax.imshow(frame.T.to_numpy() * 100, aspect='auto', vmin=0, vmax=85, cmap='viridis')
    ax.set_title(title)
    ax.set_xticks(range(4), ['B1', 'B2', 'E1', 'I1'])
    ax.set_yticks(range(19), truth.columns)
    ax.set_xlabel('Same four cell mixtures')
fig.subplots_adjust(left=0.16, right=0.88, bottom=0.12, wspace=0.12)
color_axis = fig.add_axes([0.91, 0.22, 0.015, 0.62])
fig.colorbar(im, cax=color_axis, label='Contribution (%)')
fig.savefig(out / 'same_cells_Job8_equalRNA_vs_Job9_natural.png', dpi=180, bbox_inches='tight')
plt.close(fig)

# Step8. 保存摘要和环境版本；不将无校正的模拟测试等同于真实MSCCR性能。
summary = {'status': 'numeric_review_complete_visual_pending', 'date': '2026-10-07', 'job': 'Job9', 'samples': 12, 'types': 19, 'reference': 'GSE282122 Ileum19', 'test': 'GSE282122 fixed development ileum holdout natural mixtures',
    'truth_recomputed_from_membership': True, 'training_test_donors_disjoint': True, 'development_donors': sorted(membership.Patient.unique().tolist()), 'simulation_samples_are_independent_patients': False,
    'signature_genes': len(signature), 'signature_genes_in_test': len(overlap), 'signature_overlap_percent': len(overlap) / len(signature) * 100,
    'Job9_RNA_MAE19_pp': float(error.abs().to_numpy().mean() * 100), 'Job6_RNA_MAE19_pp': float(baseline_error.abs().to_numpy().mean() * 100),
    'Job9_mean_TV20': float(sample_metrics.Job9_RNA_TV20.mean()), 'Job6_mean_TV20': float(sample_metrics.Job6_RNA_TV20.mean()),
    'samples_with_lower_MAE19': int(sample_metrics.Job9_RNA_MAE19_pp.lt(sample_metrics.Job6_RNA_MAE19_pp).sum()), 'types_with_lower_MAE': int(type_metrics.Job9_RNA_MAE_pp.lt(type_metrics.Job6_RNA_MAE_pp).sum()),
    'same_four_cells_natural_Job9_MAE19_pp': float(natural_four_MAE), 'same_four_cells_equalRNA_Job8_MAE19_pp': float(equal_four_MAE),
    'epithelial_total_MAE_Job9_pp': float(epithelium.Job9_abs_error_pp.mean()), 'epithelial_total_MAE_Job6_pp': float(epithelium.Job6_abs_error_pp.mean()),
    'web_Correlation_min': float(result.Correlation.min()), 'web_Correlation_max': float(result.Correlation.max()), 'web_Correlation_median': float(result.Correlation.median()), 'all_P_values_reported_zero': bool(result['P-value'].eq(0).all()),
    'completed_log_present': True, 'QN': False, 'permutations': 100, 'warnings_in_supplied_console': 'WARNING:' in log, 'separate_warning_panel_provided': False, 'batch_correction_explicitly_shown': False,
    'actual_uploaded_file_hash_independently_verified': False, 'PDF_table_matches_TXT_at_display_precision': True, 'reference_training_cells': pilot_summary['web_reference_cells'], 'reference_actual_training_donors': pilot_summary['actual_web_donors'],
    'original_files': file_audit, 'decision': 'Ileum19 improves on Joint20 in all 12 natural holdouts, but Paneth inflation and Stem/TA/EEC underestimation remain; retain as regional developmental candidate, not a validated production reference.'}
(out / 'review_summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
(out / 'audit_runtime_versions.json').write_text(json.dumps({name: importlib.metadata.version(name) for name in ['numpy', 'pandas', 'matplotlib', 'pypdf', 'pypdfium2']}, indent=2), encoding='utf-8')
print('SUMMARY:', json.dumps(summary, ensure_ascii=False, indent=2), flush=True)
print(construct_type_summary.loc[['Paneth', 'Stem_TA_progenitor', 'Enteroendocrine', 'Goblet', 'Tuft']].round(3).to_string(), flush=True)

# Step9. 保存以Job8/Job9配对比较为主的中文报告，辅以全部12个自然样本审计。
paired_table = ['| 同一批细胞 | Job8等RNA MAE百分点 | Job9自然混合物 MAE百分点 |', '|---|---:|---:|']
for name, row in paired_sample_metrics.iterrows():
    paired_table.append(f'| {row.paired_main_sample} | {row.Job8_equalRNA_MAE19_pp:.2f} | {row.Job9_natural_MAE19_pp:.2f} |')
focus_table = ['| 细胞类型 | 等RNA真值均值% | Job8预测均值% | 自然RNA真值均值% | Job9预测均值% | 等RNA MAE百分点 | 自然RNA MAE百分点 |', '|---|---:|---:|---:|---:|---:|---:|']
for cell_type, row in construct_type_summary.loc[['Paneth', 'Stem_TA_progenitor', 'Enteroendocrine', 'Goblet', 'Tuft']].iterrows():
    focus_table.append(f'| {cell_type} | {row.equalRNA_truth * 100:.2f} | {row.Job8_equalRNA_prediction * 100:.2f} | {row.natural_RNA_truth * 100:.2f} | {row.Job9_natural_prediction * 100:.2f} | {row.Job8_equalRNA_MAE_pp:.2f} | {row.Job9_natural_MAE_pp:.2f} |')
type_table = ['| 类型 | 自然RNA真值均值% | Job9均值% | Job9 MAE百分点 | Job6 MAE百分点 | Job9假零样本数 |', '|---|---:|---:|---:|---:|---:|']
for cell_type, row in type_metrics.iterrows():
    type_table.append(f'| {cell_type} | {row.RNA_truth_mean_percent:.2f} | {row.Job9_mean_percent:.2f} | {row.Job9_RNA_MAE_pp:.2f} | {row.Job6_RNA_MAE_pp:.2f} | {int(row.Job9_zero_when_present)}/12 |')
stem_training = training.loc[training.pilot_type.eq('Stem_TA_progenitor'), 'author_label'].value_counts()
report = f'''# Job9回肠自然混合物审阅：重点与Job8等RNA配对比较

日期：2026-10-07。结论：同细胞的等RNA测试比自然混合物更准确，但Paneth高估、Stem/TA与EEC低估在两种构造下持续存在。回肠专用参考相对联合参考有改进，保留为需要优化的候选，不据本轮宣布已适合正式MSCCR或EcoTyper输入。

## 1. 先明确比较对象和单位

Job9使用GI19_Signature19.txt预测GI19_HoldoutNatural12_Mixture_CPM.txt。自然混合物为各细胞原始counts相加，再转CPM；主真值为export_library_counts占比，属于RNA/UMI贡献近似，不能与细胞数比例直接混用。Job8先给每个细胞相同总表达权重，再混合，因此其已知RNA贡献与细胞数量比例相等。

两次采用同一回肠19类参考，按冻结配对记录只取Ileum_balanced_1、Ileum_balanced_2、Ileum_epithelial_high_1、Ileum_immune_high_1四批相同细胞比较。逐类型细胞数比例亦与等RNA真值吻合。不能拿Job9全部12个样本的均值直接对比Job8四个样本。其余8个自然样本作为扩展诊断单列。

本轮12个自然样本每个200细胞，使用CD10、CD9、H2、UC4四位开发留出供者。训练参考950细胞、15位实际入选供者，各19类50细胞；训练与测试供者无重叠。混合物复用细胞及供者，不是12名独立患者；最终留出供者CD12、UC17未使用。

## 2. 四对样本：等RNA与自然混合物的主要比较

同一四批细胞中，Job8等RNA MAE为 **{equal_four_MAE:.2f}** 个百分点；Job9自然混合物为 **{natural_four_MAE:.2f}**，增加{natural_four_MAE - equal_four_MAE:.2f}个百分点。各自与对应真值比较，而不是直接将两份预测比例相减。

{chr(10).join(paired_table)}

下面所有均值均只来自这四对样本，不混入其他自然样本：

{chr(10).join(focus_table)}

Paneth的偏差从等RNA平均高估7.86个百分点，变成自然构造高估19.81个百分点；Goblet在等RNA中均值接近真值，自然构造则低估。Tuft在等RNA中四个样本均检出，平均接近真值，自然构造有两个配对样本输出0，均值低估。

Stem/TA不是在自然构造中所有表现都更差：等RNA有3/4输出0，自然构造是2/4为0；两边仍总体低估，不能只用假零次数认定某一构造已解决问题。EEC在两边也均低估，尤其上皮富集样本。

这些比较提示RNA贡献不均与表达动态范围会影响结果，并可能放大部分偏差；但不能全部解释错误，因为等RNA下Paneth和Stem/TA仍明显失真，且自然构造已经用RNA贡献而非细胞数量作为真值。构造改变同时改变了目标比例和动态范围，四对样本不能精确分离每个因素的因果作用；Job8尚缺完整runtime log，也不能宣称所有实际设置已被独立核实相同。

## 3. 最直观的一对：上皮富集样本

Ileum_epithelial_high_1的Paneth：等RNA真值16%，Job8预测34.18%；自然RNA真值46.98%，Job9预测81.22%。自然真值本来就比细胞数比例高，但预测仍比其高出34.25个百分点，不能把81.22%当成仅由Paneth细胞RNA含量高造成的合理结果。

同一批细胞的Stem/TA：等RNA真值6%、自然RNA真值6.89%，两次都预测0。EEC：等RNA真值15%、预测5.23%；自然RNA真值18.81%、预测3.76%。同一细胞背景下，抬高Paneth的同时其他上皮类型明显偏低；这是组成分配的现象，不能据此追溯具体细胞被误识别成Paneth。

四对真值/预测比较图：`same_cells_Job8_equalRNA_vs_Job9_natural.png`。B1、B2、E1、I1对应上述四批细胞，两种构造分别显示自己的真值。

## 4. 12个自然混合物的扩展审计

Job9共同19类MAE **{summary['Job9_RNA_MAE19_pp']:.2f}** 个百分点，同样本Job6联合参考 **{summary['Job6_RNA_MAE19_pp']:.2f}**；12/12混合物均改善，19类中12类MAE较低。完整20类轴包含缺失的Colonic_absorptive真值0并保留Joint20额外预测质量，平均总变差分别{summary['Job9_mean_TV20']:.3f}、{summary['Job6_mean_TV20']:.3f}；没有重归一化共同类型。总变差不是单细胞错分率。

全部12样本平均Paneth真实RNA贡献22.72%、预测44.79%；Goblet真值10.56%、预测8.38%；Stem/TA真值10.23%、预测2.33%；EEC真值10.59%、预测3.31%。Stem/TA有7/12、Tuft有5/12为假零。注意这些是扩展12样本均值，不是前述4对样本的均值。

{chr(10).join(type_table)}

六类上皮加总的逐样本MAE为{summary['epithelial_total_MAE_Job9_pp']:.2f}个百分点，联合参考为{summary['epithelial_total_MAE_Job6_pp']:.2f}。上皮总量与内部细分都比联合参考有改善，但绝对偏差仍在。

## 5. 输入和参考覆盖：不能复用SCP259的原因判断

signature为5,706基因×19类，测试基因交集5,706，覆盖100%；本地输入SHA与冻结提交记录一致。结果12×19无缺失、无负数、每样本合计1；两套真值均从逐细胞清单重算一致。Job9日志明确QN关闭、100次置换且All done，但未显式给出batch correction模式，独立Warnings面板未提供；实际网站保存文件的哈希未单独核验。

网页拟合相关性中位数{summary['web_Correlation_median']:.4f}，范围{summary['web_Correlation_min']:.4f}–{summary['web_Correlation_max']:.4f}，P-value均报告0。它们与组成准确性不是同一指标；高相关不能推翻有真值的比例偏差。两份PDF各1页，结果PDF的三位小数显示与TXT一致，原PDF未修改；图中的Key是颜色图例，不是额外样本。

回肠训练Stem/TA的实际组成是{stem_training.to_dict()}，并非SCP259仅1个Stem的情形。对应LGR5、OLFM4、MKI67、TOP2A均已进入回肠signature并存在于测试文件；没有依据说“只是忘记加入Stem marker”。MKI67在参考中的最高类型为Goblet，也提醒增殖基因不必严格限定于Stem/TA。

Paneth的DEFA5、DEFA6均在参考和测试中存在，参考表达较强；这为检查高表达基因与参考谱尺度提供入口，当前并未通过删除基因或改变权重来证明其造成误差。应先诊断训练与留出供者各类型的表达谱分布、参考系数尺度和signature基因对拟合的贡献，不能把Paneth高估直接当成回肠真实组织的生物学结论。模拟设计的Paneth比例也不代表MSCCR真实患者的组成。

## 6. 对下一步的具体意义

1. 保留回肠专用参考，等RNA及自然构造相对联合参考均有改善，说明按部位限制参考值得继续研究；本轮不能证明正式MSCCR性能。
2. 优先处理Paneth对组成的高估和Stem/TA、EEC持续低估，先在这四对固定开发混合物上诊断，避免每次重抽细胞使比较失去解释力。
3. Tuft在等RNA中可恢复、自然构造中易低估，应检查低RNA贡献场景的检出，而不是立即删除该类。Stem/TA已有LGR5pos stem与TA覆盖，应重点检查供者间谱差异、合并后的稳定性与参考代表性，而非照搬SCP259成熟/未成熟配额。
4. 任何参考重建、基因权重或校正改动，都应一次只变一个因素，同时复查等RNA和自然构造；避免以等RNA变好就宣布自然bulk可用。最终留出供者留待锁定方案之后，当前不重抽、不运行新任务。
5. 当前只有fraction验证，没有审阅cell-type expression/G矩阵；不能据此推断presorted EcoTyper输入兼容性已解决。

## 7. 本地复查文件

重点：`same_cells_natural_vs_equalRNA_sample_metrics.tsv`、`same_cells_natural_vs_equalRNA_76_pairs.tsv`、`paired_construct_cell_type_means.tsv`和同细胞四栏对照图。

扩展审计：`paired_sample_accuracy_metrics.tsv`、`paired_cell_type_accuracy_metrics.tsv`、`all_228_paired_fraction_comparisons.tsv`、`epithelial_total_fractions.tsv`、marker和作者亚型覆盖表。

原结果与日志冻结于`original_outputs/`；摘要`review_summary.json`；PDF视觉审阅记录独立保存于`pdf_review/pdf_visual_review.json`；环境版本记录于`audit_runtime_versions.json`。脚本位于`scripts/python/review_cibersortx_job9_ileum_natural.py`，遵循项目线性研究脚本规范。

本轮仅本机审阅已有结果，没有访问服务器、修改原矩阵、代提交网页或开展recovery临床分析。
'''
(out / 'Job9回肠自然与等RNA配对审阅.md').write_text(report, encoding='utf-8')
print(paired_sample_metrics.round(3).to_string(), flush=True)
print('Report:', out / 'Job9回肠自然与等RNA配对审阅.md', flush=True)
