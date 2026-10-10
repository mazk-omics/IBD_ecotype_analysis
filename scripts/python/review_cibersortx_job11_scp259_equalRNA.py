# 审阅Job11结直肠等RNA结果，与Job10同细胞自然混合物配对，Job7联合参考作辅助。
# 在本机.tools/cibersortx_review环境运行；读Downloads/冻结真值，写Job11_review。
# 按四对固定细胞分别匹配真值；不重抽、不修改参考、不把模拟样本当独立患者。
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

# Step1. 保存原输出和日志，核对实际任务、四样本18类及本地输入哈希。
project = Path('D:/IBD_ecotype/IBD_ecotype_analysis')
ready = project / 'CIBERSORTx_ready_20261007'
metadata = ready / '核对记录_勿上传'
out = project / 'CIBERSORTx_signature_audit/colon_candidate_validation_20261007/Job11_review'
originals = out / 'original_outputs'
originals.mkdir(parents=True, exist_ok=True)
sources = [Path('C:/Users/马泽凯/Downloads') / name for name in ['CIBERSORTx_Job11_Results.txt', 'CIBERSORTx_Job11_Results.pdf', 'CIBERSORTx_Job11_stackedChart.pdf']]
sources.append(Path('C:/Users/马泽凯/.codex/attachments/541fc4c3-bb89-4762-a682-b42c595c5104/已粘贴的文本.txt'))
file_audit = []
for source in sources:
    name = 'CIBERSORTx_Job11_runtime_log.txt' if source.name == '已粘贴的文本.txt' else source.name
    target = originals / name
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    if target.exists():
        assert hashlib.sha256(target.read_bytes()).hexdigest() == digest, '已有冻结副本不同，不覆盖。'
    else:
        shutil.copy2(source, target)
    file_audit.append({'name': name, 'source': str(source), 'bytes': source.stat().st_size, 'sha256': digest})
result = pd.read_csv(originals / 'CIBERSORTx_Job11_Results.txt', sep='\t', index_col=0)
truth = pd.read_csv(metadata / 'SC18_EqualRNA4_ground_truth_RNA.tsv', sep='\t', index_col=0)
cell_truth = pd.read_csv(metadata / 'SC18_EqualRNA4_ground_truth_cell.tsv', sep='\t', index_col=0)
pairs = pd.read_csv(metadata / 'SC18_equalRNA4_pairing.tsv', sep='\t').set_index('sample').loc[truth.index]
signature = pd.read_csv(ready / '02_SCP259_ColonRectum_validation/SC18_Signature18.txt', sep='\t', index_col=0)
mixture = pd.read_csv(ready / '02_SCP259_ColonRectum_validation/SC18_HoldoutEqualRNA4_Mixture_CPM.txt', sep='\t', index_col=0)
assert result.index.is_unique and result.columns.is_unique and set(result.index) == set(truth.index)
assert set(result.columns) == set(truth.columns) | {'P-value', 'Correlation', 'RMSE'}
assert signature.columns.equals(truth.columns) and mixture.columns.equals(truth.index) and np.allclose(truth, cell_truth)
prediction = result.loc[truth.index, truth.columns]
assert prediction.shape == (4, 18) and np.isfinite(prediction.to_numpy()).all() and prediction.ge(0).all().all()
assert np.allclose(prediction.sum(axis=1), 1) and np.allclose(truth.sum(axis=1), 1)
assert signature.index.is_unique and mixture.index.is_unique and np.allclose(mixture.sum(axis=0), 1e6)
integrity = pd.read_csv(metadata / 'upload_file_integrity.tsv', sep='\t').set_index('relative_path')
for name in ['SC18_Signature18.txt', 'SC18_HoldoutEqualRNA4_Mixture_CPM.txt']:
    path = ready / '02_SCP259_ColonRectum_validation' / name
    assert hashlib.sha256(path.read_bytes()).hexdigest() == integrity.loc[str(path.relative_to(ready)), 'sha256']
log = (originals / 'CIBERSORTx_Job11_runtime_log.txt').read_text(encoding='utf-8')
assert 'SC18_HoldoutEqualRNA4_Mixture_CPM.txt' in log and 'SC18_Signature18.txt' in log
assert '[Options] QN: FALSE' in log and '[Options] perm: 100' in log and 'All done.' in log
overlap = signature.index.intersection(mixture.index)

# Step2. 同一四批细胞的自然RNA真值和等RNA真值分别重算，避免与Job10全部12样本混比。
membership = pd.read_csv(metadata / 'SC18_natural12_existing_cell_membership.tsv', sep='\t')
membership = membership.loc[membership['sample'].isin(pairs.paired_main_sample)].copy()
assert not membership.duplicated(['sample', 'cell_id']).any() and not set(membership.Patient) & {'CD12', 'UC17'}
counts = pd.crosstab(membership['sample'], membership.pilot_type).reindex(index=pairs.paired_main_sample, columns=truth.columns, fill_value=0)
recomputed_cell = counts.div(counts.sum(axis=1), axis=0)
assert counts.sum(axis=1).eq(200).all() and np.allclose(recomputed_cell.to_numpy(), truth.to_numpy())
RNA_counts = membership.groupby(['sample', 'pilot_type']).export_library_counts.sum().unstack(fill_value=0).reindex(index=pairs.paired_main_sample, columns=truth.columns, fill_value=0)
natural_truth = pd.read_csv(metadata / 'SC18_Natural12_ground_truth_RNA.tsv', sep='\t', index_col=0).loc[pairs.paired_main_sample, truth.columns]
assert np.allclose(RNA_counts.div(RNA_counts.sum(axis=1), axis=0), natural_truth)
natural_prediction = pd.read_csv(out.parent / 'Job10_review/predicted_fractions_aligned.tsv', sep='\t', index_col=0).loc[pairs.paired_main_sample, truth.columns]
# 仅改配对样本索引以进行逐项计算，不改变任何组成数值。
natural_truth.index = truth.index
natural_prediction.index = truth.index
prediction.to_csv(out / 'predicted_fractions_aligned.tsv', sep='\t')
truth.to_csv(out / 'ground_truth_equalRNA18.tsv', sep='\t')
natural_truth.to_csv(out / 'paired_natural_ground_truth_RNA18.tsv', sep='\t')
natural_prediction.to_csv(out / 'paired_Job10_natural_predictions18.tsv', sep='\t')
print('Matched four cell sets and eighteen types; signature coverage:', len(overlap), '/', len(signature), 'donors:', sorted(membership.Patient.unique()), flush=True)

# Step3. 分别计算构造内误差；同RNA单位不意味着两种构造的目标比例相同。
error = prediction - truth
natural_error = natural_prediction - natural_truth
joint20 = pd.read_csv(metadata / 'SC18_EqualRNA4_Joint20_baseline_predictions20.tsv', sep='\t', index_col=0).loc[truth.index]
baseline = joint20.loc[truth.index, truth.columns]
baseline_error = baseline - truth
truth20 = truth.reindex(columns=joint20.columns, fill_value=0)
prediction20 = prediction.reindex(columns=joint20.columns, fill_value=0)
sample_metrics = pairs.copy()
sample_metrics['Job11_equalRNA_MAE18_pp'] = error.abs().mean(axis=1) * 100
sample_metrics['Job10_natural_MAE18_pp'] = natural_error.abs().mean(axis=1) * 100
sample_metrics['Job11_minus_Job10_MAE_pp'] = sample_metrics.Job11_equalRNA_MAE18_pp - sample_metrics.Job10_natural_MAE18_pp
sample_metrics['Job7_equalRNA_MAE18_pp'] = baseline_error.abs().mean(axis=1) * 100
sample_metrics['Job11_TV'] = error.abs().sum(axis=1) / 2
sample_metrics['Job10_TV'] = natural_error.abs().sum(axis=1) / 2
sample_metrics['Job7_TV20'] = (joint20 - truth20).abs().sum(axis=1) / 2
sample_metrics['Job11_fraction_Pearson'] = prediction.corrwith(truth, axis=1)
sample_metrics['Job11_web_expression_Correlation'] = result.loc[truth.index, 'Correlation']
sample_metrics['Job11_web_RMSE'] = result.loc[truth.index, 'RMSE']
sample_metrics.to_csv(out / 'paired_sample_accuracy_metrics.tsv', sep='\t')
type_metrics = pd.DataFrame({'equalRNA_truth_mean_percent': truth.mean() * 100, 'Job11_equalRNA_mean_percent': prediction.mean() * 100,
    'natural_RNA_truth_mean_percent': natural_truth.mean() * 100, 'Job10_natural_mean_percent': natural_prediction.mean() * 100,
    'Job11_equalRNA_MAE_pp': error.abs().mean() * 100, 'Job10_natural_MAE_pp': natural_error.abs().mean() * 100, 'Job7_equalRNA_MAE_pp': baseline_error.abs().mean() * 100,
    'Job11_equalRNA_bias_pp': error.mean() * 100, 'Job10_natural_bias_pp': natural_error.mean() * 100,
    'Job11_zero_when_present': (prediction.eq(0) & truth.gt(0)).sum(), 'Job10_zero_when_present': (natural_prediction.eq(0) & natural_truth.gt(0)).sum()})
type_metrics.index.name = 'cell_type'
type_metrics.to_csv(out / 'paired_cell_type_accuracy_metrics.tsv', sep='\t')
long = truth.rename_axis(index='sample', columns='cell_type').stack().rename('equalRNA_truth').reset_index()
for name, frame in [('Job11', prediction), ('natural_RNA_truth', natural_truth), ('Job10', natural_prediction), ('Job7', baseline)]:
    values = frame.rename_axis(index='sample', columns='cell_type').stack().rename(name).reset_index()
    long = long.merge(values, on=['sample', 'cell_type'], validate='one_to_one')
long['Job11_equalRNA_bias_pp'] = (long.Job11 - long.equalRNA_truth) * 100
long['Job10_natural_bias_pp'] = (long.Job10 - long.natural_RNA_truth) * 100
long.to_csv(out / 'all_72_paired_fraction_comparisons.tsv', sep='\t', index=False)
print(sample_metrics.round(3).to_string(), flush=True)
print(type_metrics.round(3).to_string(), flush=True)

# Step4. 粗上皮总量与内部类型准确性分开；marker和训练亚型证据沿用已审阅的冻结表。
epithelial_types = ['Goblet', 'Stem_TA_progenitor', 'Enteroendocrine', 'Colonic_absorptive', 'Tuft']
epithelium = pd.DataFrame({'equalRNA_truth': truth[epithelial_types].sum(axis=1), 'Job11': prediction[epithelial_types].sum(axis=1), 'natural_RNA_truth': natural_truth[epithelial_types].sum(axis=1), 'Job10': natural_prediction[epithelial_types].sum(axis=1)})
epithelium['Job11_abs_error_pp'] = (epithelium.Job11 - epithelium.equalRNA_truth).abs() * 100
epithelium['Job10_abs_error_pp'] = (epithelium.Job10 - epithelium.natural_RNA_truth).abs() * 100
epithelium.to_csv(out / 'paired_epithelial_total_fractions.tsv', sep='\t')
for name in ['relevant_marker_coverage_and_specificity.tsv', 'SCP259_training_author_subtype_counts.tsv', 'SCP259_training_type_disease_counts.tsv']:
    shutil.copy2(out.parent / 'Job10_review' / name, out / name)
markers = pd.read_csv(out / 'relevant_marker_coverage_and_specificity.tsv', sep='\t')
markers['in_Job11_equalRNA_mixture'] = markers.gene.isin(mixture.index)
markers.to_csv(out / 'relevant_marker_coverage_and_specificity.tsv', sep='\t', index=False)

# Step5. 渲染两份原PDF所有页面，核对数值显示精度；原PDF保持不变。
render_dir = out / 'pdf_review'
render_dir.mkdir(exist_ok=True)
pdf_text = {}
for name in ['CIBERSORTx_Job11_Results.pdf', 'CIBERSORTx_Job11_stackedChart.pdf']:
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
display_lines = [line.split() for line in pdf_text['CIBERSORTx_Job11_Results.pdf']['text'][0].splitlines() if line.split() and line.split()[0] in truth.index]
display = pd.DataFrame([list(map(float, row[1:])) for row in display_lines], index=[row[0] for row in display_lines], columns=result.columns[:-1]).loc[truth.index]
assert np.allclose(display, result.loc[truth.index, display.columns], atol=0.00051)

# Step6. 四批细胞统一轴展示两种真值与预测；偏差图使用百分点。
fig, axes = plt.subplots(1, 4, figsize=(16, 7), sharey=True)
for ax, frame, title in zip(axes, [truth, prediction, natural_truth, natural_prediction], ['Equal RNA truth', 'Job11 SCP259 equal RNA', 'Natural RNA truth', 'Job10 SCP259 natural RNA']):
    im = ax.imshow(frame.T.to_numpy() * 100, aspect='auto', vmin=0, vmax=40, cmap='viridis')
    ax.set_title(title, fontsize=11)
    ax.set_xticks(range(4), ['B1', 'B2', 'E1', 'I1'])
    ax.set_yticks(range(18), truth.columns)
    ax.set_xlabel('Same four cell mixtures')
fig.subplots_adjust(left=0.16, right=0.88, bottom=0.12, wspace=0.12)
color_axis = fig.add_axes([0.91, 0.22, 0.015, 0.62])
fig.colorbar(im, cax=color_axis, label='Contribution (%)')
fig.savefig(out / 'same_cells_Job11_equalRNA_vs_Job10_natural.png', dpi=180, bbox_inches='tight')
plt.close(fig)
fig, axes = plt.subplots(1, 2, figsize=(10, 7), sharey=True)
for ax, frame, title in zip(axes, [error, natural_error], ['Job11 equal RNA bias', 'Job10 natural RNA bias']):
    im = ax.imshow(frame.T.to_numpy() * 100, aspect='auto', vmin=-30, vmax=30, cmap='RdBu_r')
    ax.set_title(title)
    ax.set_xticks(range(4), ['B1', 'B2', 'E1', 'I1'])
    ax.set_yticks(range(18), truth.columns)
fig.subplots_adjust(left=0.24, right=0.85, wspace=0.14)
color_axis = fig.add_axes([0.88, 0.22, 0.02, 0.62])
fig.colorbar(im, cax=color_axis, label='Bias (percentage points)')
fig.savefig(out / 'Job11_vs_Job10_paired_bias.png', dpi=180, bbox_inches='tight')
plt.close(fig)

# Step7. 保存摘要，设置缺失与跨数据集局限不依靠文件名补齐。
summary = {'status': 'numeric_review_complete_visual_pending', 'date': '2026-10-07', 'job': 'Job11', 'paired_job': 'Job10', 'samples': 4, 'types': 18,
    'reference_source': 'SCP259', 'test_source': 'GSE282122 fixed development ColonRectum cell mixtures', 'primary_comparison': 'Same four cell mixtures; equal-RNA and natural RNA predictions evaluated against construct-specific truth',
    'truth_recomputed_from_membership': True, 'development_donors': sorted(membership.Patient.unique().tolist()), 'simulation_samples_are_independent_patients': False,
    'Job11_equalRNA_MAE18_pp': float(error.abs().to_numpy().mean() * 100), 'Job10_paired_natural_MAE18_pp': float(natural_error.abs().to_numpy().mean() * 100), 'Job7_equalRNA_MAE18_pp': float(baseline_error.abs().to_numpy().mean() * 100),
    'Job11_mean_TV': float(sample_metrics.Job11_TV.mean()), 'Job10_paired_mean_TV': float(sample_metrics.Job10_TV.mean()), 'Job7_equalRNA_mean_TV20': float(sample_metrics.Job7_TV20.mean()),
    'samples_equalRNA_MAE_lower_than_natural': int(sample_metrics.Job11_equalRNA_MAE18_pp.lt(sample_metrics.Job10_natural_MAE18_pp).sum()),
    'samples_equalRNA_MAE_higher_than_joint': int(sample_metrics.Job11_equalRNA_MAE18_pp.gt(sample_metrics.Job7_equalRNA_MAE18_pp).sum()),
    'Job11_epithelial_total_MAE_pp': float(epithelium.Job11_abs_error_pp.mean()), 'Job10_paired_epithelial_total_MAE_pp': float(epithelium.Job10_abs_error_pp.mean()),
    'signature_genes': len(signature), 'signature_overlap_genes': len(overlap), 'signature_overlap_percent': len(overlap) / len(signature) * 100,
    'web_Correlation_min': float(result.Correlation.min()), 'web_Correlation_max': float(result.Correlation.max()), 'web_Correlation_median': float(result.Correlation.median()), 'all_P_values_reported_zero': bool(result['P-value'].eq(0).all()),
    'completed_log_present': True, 'QN': False, 'permutations': 100, 'warnings_in_supplied_console': 'WARNING:' in log, 'separate_warning_panel_provided': False, 'batch_correction_explicitly_shown': False,
    'actual_uploaded_file_hash_independently_verified': False, 'PDF_table_matches_TXT_at_display_precision': True, 'original_files': file_audit,
    'decision': 'Equal-RNA input does not resolve major SCP259 cross-dataset biases; retain as diagnostic reference, not the production MSCCR signature.'}
(out / 'review_summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
(out / 'audit_runtime_versions.json').write_text(json.dumps({name: importlib.metadata.version(name) for name in ['numpy', 'pandas', 'matplotlib', 'pypdf', 'pypdfium2']}, indent=2), encoding='utf-8')
print('SUMMARY:', json.dumps(summary, ensure_ascii=False, indent=2), flush=True)

# Step8. 保存以Job11/Job10同细胞配对为主的报告，辅助基线不能替代跨队列验证。
sample_table = ['| 同一批细胞 | Job11等RNA MAE百分点 | Job10自然RNA MAE百分点 | 联合参考等RNA MAE百分点 |', '|---|---:|---:|---:|']
for name, row in sample_metrics.iterrows():
    sample_table.append(f'| {row.paired_main_sample} | {row.Job11_equalRNA_MAE18_pp:.2f} | {row.Job10_natural_MAE18_pp:.2f} | {row.Job7_equalRNA_MAE18_pp:.2f} |')
type_table = ['| 类型 | 等RNA真值均值% | Job11预测均值% | 自然RNA真值均值% | Job10预测均值% | 等RNA MAE百分点 | 自然RNA MAE百分点 |', '|---|---:|---:|---:|---:|---:|---:|']
for cell_type, row in type_metrics.iterrows():
    type_table.append(f'| {cell_type} | {row.equalRNA_truth_mean_percent:.2f} | {row.Job11_equalRNA_mean_percent:.2f} | {row.natural_RNA_truth_mean_percent:.2f} | {row.Job10_natural_mean_percent:.2f} | {row.Job11_equalRNA_MAE_pp:.2f} | {row.Job10_natural_MAE_pp:.2f} |')
report = f'''# Job11结直肠等RNA与Job10自然混合物配对审阅

日期：2026-10-07。结论：等RNA构造没有消除这版SCP259参考的主要跨数据集偏差。Goblet、吸收细胞、CD8_T及单核/巨噬细胞持续低估，Tuft、ILC与Stem/TA持续高估。保留为诊断对照，不支持作为正式MSCCR参考；不是否定完整SCP259数据或CIBERSORTx算法。

## 1. 核对比较对象

Job11是SC18_HoldoutEqualRNA4_Mixture_CPM.txt、SC18_Signature18.txt；Job10是同参考的自然counts求和混合物。参考来自SCP259，测试来自GSE282122结直肠开发留出细胞，属于跨数据集测试，不是SCP259自身供者留出。

两次配对ColonRectum_balanced_1、balanced_2、epithelial_high_1、immune_high_1四批相同细胞，各200个。Job10其余8个样本不进入本报告主要比较，因此本报告自然MAE4.49不是之前全部12样本的4.54。

等RNA先赋每个细胞相同总表达权重，真值等于细胞数量比例；自然混合物求和原始counts再CPM，主真值是UMI/RNA贡献占比。两种构造各自使用自己的真值评价，不能把预测比例直接相减当成误差改变。两套真值已按逐细胞清单重算一致，供者仍是CD10、CD9、H2、UC11、UC15、UC4；不同混合物复用供者及细胞，不是4个独立患者。最终留出供者CD12、UC17未使用。

## 2. 总体变化较小，而且并非每个样本都改善

共同18类平均绝对误差：Job11等RNA **{summary['Job11_equalRNA_MAE18_pp']:.2f}** 个百分点，Job10对应自然混合物 **{summary['Job10_paired_natural_MAE18_pp']:.2f}**。仅下降{summary['Job10_paired_natural_MAE18_pp'] - summary['Job11_equalRNA_MAE18_pp']:.2f}个百分点，3个样本较低、balanced_2反而较高。平均总变差分别{summary['Job11_mean_TV']:.3f}、{summary['Job10_paired_mean_TV']:.3f}，是组成L1距离的一半，不代表单细胞错分率。

{chr(10).join(sample_table)}

辅助对照：联合参考在同四个等RNA样本、共同18类上的MAE为{summary['Job7_equalRNA_MAE18_pp']:.2f}个百分点，4/4均优于SCP259；共同类型不重新归一化。完整20类轴保留额外Paneth/Ileal_absorptive预测质量，联合参考平均总变差{summary['Job7_equalRNA_mean_TV20']:.3f}。这些额外类型在本轮真值为0。联合参考与测试共享GSE282122研究背景，虽供者留出但有域内优势，不能据此保证真实MSCCR同样排名。

## 3. 类型偏差在两种构造中持续存在

以下均值和MAE全部基于相同四批细胞，不混入Job10另外8个自然样本：

{chr(10).join(type_table)}

最明显的持续偏差：

- Goblet：等RNA真值6.88%、预测0.08%，自然RNA真值13.23%、预测1.54%。等RNA3/4、自然2/4样本为假零，等RNA没有恢复可靠检出。
- 结肠吸收细胞：等RNA7.12%→1.17%，自然10.84%→2.64%，两边都低估。
- Tuft：等RNA8.12%→19.90%，自然3.38%→19.33%，两边都高估。自然真值本来因RNA贡献较低而减少，但预测均值仍维持约20%；需要检查参考拟合中的竞争分配，不能将其解释成测试细胞真的出现大量Tuft。
- ILC：等RNA4.50%→15.46%，自然2.41%→11.42%，等RNA平均高估10.96个百分点，比自然的9.01更大。
- CD8_T：等RNA5.12%→0.29%，自然2.13%→0；等RNA假零从4/4降为2/4，但检出非零的两个值也只有约0.64%和0.50%，仍明显低估。
- 单核/巨噬细胞：等RNA6.12%→0.79%，自然4.70%→0.43%，两种构造各有2/4假零。
- Stem/TA：等RNA7.62%→12.32%，自然16.30%→25.00%，持续高估。该方向与回肠参考的Stem/TA低估不同，不能将全部肠道数据归为同一个错误模式。

## 4. 粗上皮趋势不能替代细分准确性

Goblet、Stem/TA、EEC、吸收细胞、Tuft加总：等RNA真值均值{epithelium.equalRNA_truth.mean() * 100:.2f}%、Job11均值{epithelium.Job11.mean() * 100:.2f}%；自然RNA真值{epithelium.natural_RNA_truth.mean() * 100:.2f}%、Job10均值{epithelium.Job10.mean() * 100:.2f}%。逐样本粗上皮MAE分别{summary['Job11_epithelial_total_MAE_pp']:.2f}、{summary['Job10_paired_epithelial_total_MAE_pp']:.2f}个百分点。

内部的Goblet/吸收细胞不足与Tuft/Stem过量可能在总量上相互抵消。堆叠图可见上皮富集的总体趋势，但不支持18类精细分配可靠；不存在单细胞预测标签，不能据比例反推某一类细胞被逐一识别成另一类。

## 5. 输入检查与结果可信范围

两次日志均列相同signature、100次置换、QN关闭，Job11有All done。输出4×18，无缺失、无负值、各行合计1；两份PDF各1页，表格显示值与TXT在三位小数精度一致。原始PDF及TXT均保留SHA核对副本，没有修改。

signature基因4,612，测试覆盖4,381，交集{summary['signature_overlap_percent']:.2f}%；不是整体基因交集过低。marker明细另保存，Goblet诊断基因及CD8A/CD8B在参考和测试中存在，不能由这次失败直接认定它们都未进入矩阵。

网页表达拟合Correlation中位数{summary['web_Correlation_median']:.3f}，范围{summary['web_Correlation_min']:.3f}–{summary['web_Correlation_max']:.3f}，P-value全部报告0。它们不等同于细胞比例准确性；balanced_2网页Correlation约0.925，但真实比例向量与预测向量的相关性约−0.155，说明两个Correlation概念不能混用。

两次日志均未显式显示batch correction设置，独立Warnings面板未提供；原提交计划为关闭，实际仍以Job Parameters为准。原文件名、运行参数和本地SHA吻合不等于独立核验网站文件哈希。因此这里描述的是当前两次输出，不声称已证明完全收敛或某一种生物学机制。

## 6. 对参考选择和后续诊断的意义

等RNA去掉每个细胞总RNA贡献的差异后，严重偏差仍在，所以“仅由RNA量不均造成”不足以解释SCP259结果。另一方面，两种构造改变了目标比例和动态范围，不能用4对混合物估计每个因素的独立因果效应。

现有900细胞SCP259简化参考每类50个，Goblet以43个Immature Goblet为主，Stem/TA仅1个Stem，其余主要为TA和前体。跨队列标签边界、内部亚型代表性、平台及疾病背景差异仍是待检验的方向。ILC中部分髓系共享marker的高表达也值得回查，但当前结果不证明细胞污染或注释错误。

优先顺序建议：先降低这版SC18作为主候选的优先级；若继续挽救，先做SCP259自身供者留出与跨队列对照，区分参考内部可辨识性与迁移问题。之后再固定这四对开发混合物，单独测试有匹配reference的校正或重新定义亚型覆盖；不要同时改标签、抽样和校正后难以解释改善。

回肠GI19和结直肠SC18都已完成等RNA与自然两种构造的审阅：前者相对联合参考有改善但Paneth/Stem/EEC问题持续，后者两种构造均偏差明显。不能跨部位直接比较MAE来判断哪种组织更容易，也不能把当前SCP259小试结果推广到完整图谱。正式MSCCR及下游G矩阵适配仍需另验；本轮没有进行recovery临床关联。

## 7. 文件入口

- `paired_sample_accuracy_metrics.tsv`：四对样本各自真值下的误差和辅助联合参考。
- `paired_cell_type_accuracy_metrics.tsv`：18类均值、MAE、方向及假零数。
- `all_72_paired_fraction_comparisons.tsv`：4×18逐项对照。
- `same_cells_Job11_equalRNA_vs_Job10_natural.png`、`Job11_vs_Job10_paired_bias.png`：统一轴组成与偏差图，B1/B2/E1/I1对应四对样本。
- `paired_epithelial_total_fractions.tsv`：粗上皮总量。
- marker、作者亚型及疾病覆盖表：沿用Job10审阅证据并补Job11基因交集。
- `original_outputs/`、`review_summary.json`、`pdf_review/pdf_visual_review.json`：原始副本、摘要及视觉审阅记录。

脚本：`scripts/python/review_cibersortx_job11_scp259_equalRNA.py`；环境：`.tools/cibersortx_review`，版本记录于`audit_runtime_versions.json`。本轮仅本机审阅，无服务器访问、重抽或新任务提交。
'''
(out / 'Job11结直肠等RNA与自然配对审阅.md').write_text(report, encoding='utf-8')
print('Report:', out / 'Job11结直肠等RNA与自然配对审阅.md', flush=True)
