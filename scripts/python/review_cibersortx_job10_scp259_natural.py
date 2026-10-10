# 审阅 Job10：SCP259 参考预测 GSE282122 结直肠固定留出混合物，与 Job6 同样本配对比较。
# 在本机 .tools/cibersortx_review 环境运行；读 Downloads/冻结真值，写 colon_candidate_validation_20261007/Job10_review。
# 主指标为 RNA 贡献比例误差；不改参考、真值、原 PDF，也不把模拟样本当独立患者。
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
import pypdfium2 as pdfium
from pypdf import PdfReader

# Step1. 保存原始文件和日志副本，核对输入记录、样本名称与18类轴。
project = Path('D:/IBD_ecotype/IBD_ecotype_analysis')
ready = project / 'CIBERSORTx_ready_20261007'
metadata = ready / '核对记录_勿上传'
out = project / 'CIBERSORTx_signature_audit/colon_candidate_validation_20261007/Job10_review'
originals = out / 'original_outputs'
originals.mkdir(parents=True, exist_ok=True)
sources = [Path('C:/Users/马泽凯/Downloads') / name for name in ['CIBERSORTx_Job10_Results.txt', 'CIBERSORTx_Job10_Results.pdf', 'CIBERSORTx_Job10_stackedChart.pdf']]
sources.append(Path('C:/Users/马泽凯/.codex/attachments/d840ffc5-775e-4848-b6e5-aab51cc3a107/已粘贴的文本.txt'))
file_audit = []
for source in sources:
    name = 'CIBERSORTx_Job10_runtime_log.txt' if source.name == '已粘贴的文本.txt' else source.name
    target = originals / name
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    if target.exists():
        assert hashlib.sha256(target.read_bytes()).hexdigest() == digest, '已有冻结副本不同，不覆盖。'
    else:
        shutil.copy2(source, target)
    file_audit.append({'name': name, 'source': str(source), 'bytes': source.stat().st_size, 'sha256': digest})
job10 = pd.read_csv(originals / 'CIBERSORTx_Job10_Results.txt', sep='\t', index_col=0)
truth = pd.read_csv(metadata / 'SC18_Natural12_ground_truth_RNA.tsv', sep='\t', index_col=0)
cell_truth = pd.read_csv(metadata / 'SC18_Natural12_ground_truth_cell.tsv', sep='\t', index_col=0)
job6 = pd.read_csv(metadata / 'SC18_Natural12_Joint20_baseline_predictions20.tsv', sep='\t', index_col=0)
design = pd.read_csv(metadata / 'SC18_natural12_sample_design.tsv', sep='\t').set_index('sample').loc[truth.index]
membership = pd.read_csv(metadata / 'SC18_natural12_existing_cell_membership.tsv', sep='\t')
signature = pd.read_csv(ready / '02_SCP259_ColonRectum_validation/SC18_Signature18.txt', sep='\t', index_col=0)
mixture = pd.read_csv(ready / '02_SCP259_ColonRectum_validation/SC18_HoldoutNatural12_Mixture_CPM.txt', sep='\t', index_col=0)
assert job10.index.is_unique and set(job10.index) == set(truth.index)
assert set(job10.columns) == set(truth.columns) | {'P-value', 'Correlation', 'RMSE'}
assert signature.columns.equals(truth.columns) and mixture.columns.equals(truth.index)
prediction = job10.loc[truth.index, truth.columns]
assert prediction.shape == (12, 18) and np.isfinite(prediction.to_numpy()).all() and prediction.ge(0).all().all()
assert np.allclose(prediction.sum(axis=1), 1) and np.allclose(truth.sum(axis=1), 1) and np.allclose(cell_truth.sum(axis=1), 1)
integrity = pd.read_csv(metadata / 'upload_file_integrity.tsv', sep='\t').set_index('relative_path')
for name in ['SC18_Signature18.txt', 'SC18_HoldoutNatural12_Mixture_CPM.txt']:
    path = ready / '02_SCP259_ColonRectum_validation' / name
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    assert digest == integrity.loc[str(path.relative_to(ready)), 'sha256']
log = (originals / 'CIBERSORTx_Job10_runtime_log.txt').read_text(encoding='utf-8')
assert 'SC18_HoldoutNatural12_Mixture_CPM.txt' in log and 'SC18_Signature18.txt' in log
assert '[Options] perm: 100' in log and '[Options] QN: FALSE' in log and 'All done.' in log

# Step2. 用冻结的逐细胞组成重新算两套真值，避免仅凭文件名判断已知比例。
assert not membership.duplicated(['sample', 'cell_id']).any()
counts = pd.crosstab(membership['sample'], membership.pilot_type).reindex(index=truth.index, columns=truth.columns, fill_value=0)
recomputed_cell = counts.div(counts.sum(axis=1), axis=0)
RNA_counts = membership.groupby(['sample', 'pilot_type']).export_library_counts.sum().unstack(fill_value=0).reindex(index=truth.index, columns=truth.columns, fill_value=0)
recomputed_RNA = RNA_counts.div(RNA_counts.sum(axis=1), axis=0)
assert np.allclose(recomputed_cell, cell_truth) and np.allclose(recomputed_RNA, truth)
assert counts.sum(axis=1).eq(200).all() and not set(membership.Patient) & {'CD12', 'UC17'}
assert mixture.index.is_unique and signature.index.is_unique and np.allclose(mixture.sum(axis=0), 1e6)
overlap = signature.index.intersection(mixture.index)
truth.to_csv(out / 'ground_truth_RNA18.tsv', sep='\t')
cell_truth.to_csv(out / 'ground_truth_cell18.tsv', sep='\t')
prediction.to_csv(out / 'predicted_fractions_aligned.tsv', sep='\t')
print('Matched:', prediction.shape, 'signature genes present:', len(overlap), '/', len(signature), 'donors:', membership.Patient.nunique(), flush=True)

# Step3. 同轴比较两份结果；18类主指标不重归一化，20类另保留缺失类型的误预测质量。
error = prediction - truth
baseline = job6.loc[truth.index, truth.columns]
baseline_error = baseline - truth
cell_error = prediction - cell_truth
types20 = job6.columns
truth20 = truth.reindex(columns=types20, fill_value=0)
prediction20 = prediction.reindex(columns=types20, fill_value=0)
error20 = prediction20 - truth20
baseline_error20 = job6.loc[truth.index] - truth20
sample_metrics = design.copy()
sample_metrics['Job10_RNA_MAE18_pp'] = error.abs().mean(axis=1) * 100
sample_metrics['Job6_RNA_MAE18_pp'] = baseline_error.abs().mean(axis=1) * 100
sample_metrics['Job10_RNA_MAE20_pp'] = error20.abs().mean(axis=1) * 100
sample_metrics['Job6_RNA_MAE20_pp'] = baseline_error20.abs().mean(axis=1) * 100
sample_metrics['Job10_RNA_TV'] = error20.abs().sum(axis=1) / 2
sample_metrics['Job6_RNA_TV'] = baseline_error20.abs().sum(axis=1) / 2
sample_metrics['Job10_cell_MAE18_pp'] = cell_error.abs().mean(axis=1) * 100
sample_metrics['Job6_cell_MAE18_pp'] = (baseline - cell_truth).abs().mean(axis=1) * 100
sample_metrics['Job10_fraction_Pearson_RNA'] = prediction.corrwith(truth, axis=1)
sample_metrics['Job10_web_expression_correlation'] = job10.loc[truth.index, 'Correlation']
sample_metrics['Job10_web_expression_RMSE'] = job10.loc[truth.index, 'RMSE']
sample_metrics.to_csv(out / 'paired_sample_accuracy_metrics.tsv', sep='\t')
type_metrics = pd.DataFrame({'RNA_truth_mean_percent': truth.mean() * 100, 'cell_truth_mean_percent': cell_truth.mean() * 100, 'Job10_mean_percent': prediction.mean() * 100, 'Job6_mean_percent': baseline.mean() * 100,
    'Job10_RNA_MAE_pp': error.abs().mean() * 100, 'Job6_RNA_MAE_pp': baseline_error.abs().mean() * 100, 'Job10_bias_pp': error.mean() * 100,
    'true_present_samples': truth.gt(0).sum(), 'Job10_zero_when_present': (prediction.eq(0) & truth.gt(0)).sum(), 'RNA_truth_min_percent': truth.min() * 100, 'RNA_truth_max_percent': truth.max() * 100})
type_metrics.index.name = 'cell_type'
type_metrics.to_csv(out / 'paired_cell_type_accuracy_metrics.tsv', sep='\t')
scenario_metrics = sample_metrics.groupby('scenario')[['Job10_RNA_MAE18_pp', 'Job6_RNA_MAE18_pp', 'Job10_RNA_TV', 'Job6_RNA_TV', 'Job10_cell_MAE18_pp', 'Job6_cell_MAE18_pp']].mean()
scenario_metrics.to_csv(out / 'scenario_accuracy_metrics.tsv', sep='\t')
long = truth.rename_axis(index='sample', columns='cell_type').stack().rename('RNA_truth').reset_index()
for name, frame in [('cell_truth', cell_truth), ('Job10', prediction), ('Job6', baseline)]:
    values = frame.rename_axis(index='sample', columns='cell_type').stack().rename(name).reset_index()
    long = long.merge(values, on=['sample', 'cell_type'], validate='one_to_one')
long['Job10_bias_pp'] = (long.Job10 - long.RNA_truth) * 100
long.to_csv(out / 'all_216_paired_fraction_comparisons.tsv', sep='\t', index=False)
print(scenario_metrics.to_string(), flush=True)
print(type_metrics.round(3).to_string(), flush=True)

# Step4. 检查粗谱系总量；总量准确不意味着谱系内部类型准确。
lineage_types = {'Epithelial': ['Goblet', 'Stem_TA_progenitor', 'Enteroendocrine', 'Colonic_absorptive', 'Tuft'],
    'Immune': ['B', 'Plasma', 'CD4_T', 'CD8_T', 'NK', 'ILC', 'Monocyte_Macrophage', 'DC', 'Mast'], 'Stromal_with_glia': ['Glial', 'Endothelial', 'Pericyte', 'Fibroblast']}
lineage_rows = []
for lineage, types in lineage_types.items():
    for sample in truth.index:
        lineage_rows.append({'sample': sample, 'scenario': design.loc[sample, 'scenario'], 'lineage': lineage, 'RNA_truth': truth.loc[sample, types].sum(), 'Job10': prediction.loc[sample, types].sum(), 'Job6': baseline.loc[sample, types].sum()})
lineages = pd.DataFrame(lineage_rows)
lineages['Job10_abs_error_pp'] = (lineages.Job10 - lineages.RNA_truth).abs() * 100
lineages.to_csv(out / 'coarse_lineage_fractions.tsv', sep='\t', index=False)
lineages.groupby(['scenario', 'lineage'])[['RNA_truth', 'Job10', 'Job6', 'Job10_abs_error_pp']].mean().to_csv(out / 'coarse_lineage_scenario_means.tsv', sep='\t')

# Step5. 保存与本轮偏差有关的既有 marker 证据和训练细胞覆盖；不据此宣称已证明因果。
candidate = project / 'CIBERSORTx_signature_audit/three_candidate_comparison_20261006/scp_rectum'
markers = pd.read_csv(candidate / 'marker_panel_audit.tsv', sep='\t')
markers = markers.loc[markers.pilot_type.isin(['Goblet', 'Tuft', 'Stem_TA_progenitor', 'Colonic_absorptive', 'CD8_T', 'ILC', 'Monocyte_Macrophage'])].copy()
markers['in_test_mixture'] = markers.gene.isin(mixture.index)
markers.to_csv(out / 'relevant_marker_coverage_and_specificity.tsv', sep='\t', index=False)
train_dir = project / 'IBDverse_reference/three_reference_pilots_20261006/run_20261006_183211/scp259'
training = pd.read_csv(train_dir / 'selected_cell_manifest.tsv', sep='\t')
training.groupby(['pilot_type', 'author_label']).size().rename('cells').to_csv(out / 'SCP259_training_author_subtype_counts.tsv', sep='\t')
training.groupby(['pilot_type', 'Disease']).size().rename('cells').to_csv(out / 'SCP259_training_type_disease_counts.tsv', sep='\t')
pilot_audit = json.loads((train_dir / 'summary.json').read_text(encoding='utf-8'))

# Step6. 提取和渲染完整 PDF；显示值与 TXT 比较，图像另做人工视觉审阅。
render_dir = out / 'pdf_review'
render_dir.mkdir(exist_ok=True)
pdf_texts = {}
for name in ['CIBERSORTx_Job10_Results.pdf', 'CIBERSORTx_Job10_stackedChart.pdf']:
    reader = PdfReader(originals / name)
    pdf_texts[name] = {'pages': len(reader.pages), 'text': [page.extract_text() for page in reader.pages]}
    document = pdfium.PdfDocument(str(originals / name))
    for i in range(len(document)):
        page = document[i]
        bitmap = page.render(scale=2.4)
        bitmap.to_pil().save(render_dir / f'{Path(name).stem}_page{i + 1}.png')
        bitmap.close()
        page.close()
    document.close()
(render_dir / 'pdf_text_extraction.json').write_text(json.dumps(pdf_texts, ensure_ascii=False, indent=2), encoding='utf-8')
text = pdf_texts['CIBERSORTx_Job10_Results.pdf']['text'][0]
display_rows = [line.split() for line in text.splitlines() if line.split() and line.split()[0] in truth.index]
display = pd.DataFrame([list(map(float, row[1:])) for row in display_rows], index=[row[0] for row in display_rows], columns=list(job10.columns[:-1])).loc[truth.index]
assert np.allclose(display, job10.loc[truth.index, display.columns], atol=0.00051)

# Step7. 画真值与预测对照，以及误差图；保持百分比和百分点单位分开。
plot_types = truth.columns.tolist()
fig, axes = plt.subplots(1, 3, figsize=(17, 7), sharey=True)
for ax, frame, title in zip(axes, [truth, baseline, prediction], ['Known RNA contribution', 'Job6 Joint20 (same 18 types)', 'Job10 SCP259 18 types']):
    im = ax.imshow(frame.T.to_numpy() * 100, aspect='auto', vmin=0, vmax=40, cmap='viridis')
    ax.set_title(title)
    ax.set_xticks(range(12), ['B1','B2','B3','B4','E1','E2','E3','E4','I1','I2','I3','I4'], rotation=45)
    ax.set_yticks(range(18), plot_types)
    ax.set_xlabel('Balanced / epithelial-high / immune-high')
fig.subplots_adjust(left=0.15, right=0.87, wspace=0.08, bottom=0.14)
color_axis = fig.add_axes([0.90, 0.22, 0.015, 0.62])
fig.colorbar(im, cax=color_axis, label='RNA fraction (%)')
fig.savefig(out / 'truth_vs_Job6_vs_Job10.png', dpi=180, bbox_inches='tight')
plt.close(fig)
fig, ax = plt.subplots(figsize=(9, 7))
im = ax.imshow(error.T.to_numpy() * 100, aspect='auto', vmin=-30, vmax=30, cmap='RdBu_r')
ax.set_xticks(range(12), ['B1','B2','B3','B4','E1','E2','E3','E4','I1','I2','I3','I4'])
ax.set_yticks(range(18), plot_types)
ax.set_title('Job10 minus known RNA contribution')
fig.colorbar(im, ax=ax, label='Bias (percentage points)')
fig.tight_layout()
fig.savefig(out / 'Job10_fraction_bias.png', dpi=180)
plt.close(fig)

# Step8. 保存摘要，实际网页版上传内容与未显示的设置不能靠本地文件替代证明。
summary = {'status': 'numeric_review_complete_visual_pending', 'job': 'Job10', 'samples': 12, 'types': 18, 'test_source': 'GSE282122 fixed development holdout mixtures', 'reference_source': 'SCP259',
    'natural_counts_then_CPM': True, 'truth_recomputed_from_membership': True, 'simulation_samples_are_independent_patients': False, 'development_donors': sorted(membership.Patient.unique().tolist()),
    'signature_overlap_genes': len(overlap), 'signature_genes': len(signature), 'signature_overlap_percent': len(overlap) / len(signature) * 100,
    'Job10_RNA_MAE18_pp': float(error.abs().to_numpy().mean() * 100), 'Job6_RNA_MAE18_pp': float(baseline_error.abs().to_numpy().mean() * 100),
    'Job10_RNA_MAE20_pp': float(error20.abs().to_numpy().mean() * 100), 'Job6_RNA_MAE20_pp': float(baseline_error20.abs().to_numpy().mean() * 100),
    'Job10_mean_TV': float(sample_metrics.Job10_RNA_TV.mean()), 'Job6_mean_TV': float(sample_metrics.Job6_RNA_TV.mean()),
    'samples_with_higher_MAE18': int(sample_metrics.Job10_RNA_MAE18_pp.gt(sample_metrics.Job6_RNA_MAE18_pp).sum()), 'types_with_lower_MAE': int(type_metrics.Job10_RNA_MAE_pp.lt(type_metrics.Job6_RNA_MAE_pp).sum()),
    'web_expression_correlation_min': float(job10.Correlation.min()), 'web_expression_correlation_max': float(job10.Correlation.max()), 'web_expression_correlation_median': float(job10.Correlation.median()),
    'all_P_values_reported_zero': bool(job10['P-value'].eq(0).all()), 'completed_log_present': True, 'permutations': 100, 'QN': False,
    'solver_warning_in_supplied_console_text': 'WARNING:' in log, 'separate_warning_panel_provided': False, 'batch_correction_setting_explicitly_shown': False, 'actual_uploaded_file_hash_independently_verified': False,
    'PDF_table_matches_TXT_at_display_precision': True, 'reference_training_cells': pilot_audit['web_reference_cells'], 'reference_training_donors': pilot_audit['actual_web_donors'],
    'original_files': file_audit, 'decision': 'Do not use this SCP259 pilot as the production MSCCR reference on this cross-dataset result; retain as developmental diagnostic comparison.'}
(out / 'review_summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
versions = {name: importlib.metadata.version(name) for name in ['numpy', 'pandas', 'matplotlib', 'pypdf', 'pypdfium2']}
(out / 'audit_runtime_versions.json').write_text(json.dumps(versions, indent=2), encoding='utf-8')
print('SUMMARY:', json.dumps(summary, ensure_ascii=False, indent=2), flush=True)

# Step9. 将本轮结论和限制写成中文报告；百分比为组成，百分点为预测误差。
type_table = ['| 细胞类型 | RNA真值均值% | Job10均值% | Job6均值% | Job10 MAE百分点 | Job6 MAE百分点 | 真值阳性但Job10为零 |', '|---|---:|---:|---:|---:|---:|---:|']
for cell_type, row in type_metrics.iterrows():
    type_table.append(f'| {cell_type} | {row.RNA_truth_mean_percent:.2f} | {row.Job10_mean_percent:.2f} | {row.Job6_mean_percent:.2f} | {row.Job10_RNA_MAE_pp:.2f} | {row.Job6_RNA_MAE_pp:.2f} | {int(row.Job10_zero_when_present)}/12 |')
scenario_table = ['| 混合情景 | Job10 MAE百分点 | Job6 MAE百分点 | Job10总变差 | Job6总变差 |', '|---|---:|---:|---:|---:|']
for scenario, row in scenario_metrics.iterrows():
    scenario_table.append(f'| {scenario} | {row.Job10_RNA_MAE18_pp:.2f} | {row.Job6_RNA_MAE18_pp:.2f} | {row.Job10_RNA_TV:.3f} | {row.Job6_RNA_TV:.3f} |')
epithelium = lineages.loc[lineages.lineage.eq('Epithelial')]
lineage_table = ['| 混合情景 | 上皮RNA真值均值% | Job10上皮均值% | 上皮总量MAE百分点 |', '|---|---:|---:|---:|']
for scenario, frame in epithelium.groupby('scenario'):
    lineage_table.append(f'| {scenario} | {frame.RNA_truth.mean() * 100:.2f} | {frame.Job10.mean() * 100:.2f} | {frame.Job10_abs_error_pp.mean():.2f} |')
report = f'''# Job10：SCP259结直肠参考的自然混合物审阅

审阅日期：2026-10-07。结论：计算完成，表格结构及真值匹配通过；当前小规模SCP259参考在这一轮跨数据集组成预测中偏差明显，审计不支持直接提升为正式MSCCR反卷积参考。此处没有预设的正式通过阈值，也不能据此否定SCP259原始数据集或CIBERSORTx算法。PDF视觉审阅状态单独记录于 `pdf_review/pdf_visual_review.json`。

## 1. 实际验证了什么

参考是SCP259的18类signature；测试是GSE282122的12个结直肠固定开发留出混合物，并非SCP259自身留出细胞。每个混合物200个细胞，balanced、epithelial_high、immune_high各4个；总体使用供者CD10、CD9、H2、UC11、UC15、UC4。混合物之间复用了细胞及供者，不能当12位独立患者进行显著性推断。严格最终留出供者CD12、UC17未使用。

自然混合物先相加逐细胞原始counts，再归一化为CPM。因此主真值为各类型的export_library_counts占总UMI比例；这是可计算的RNA贡献近似，不等同于真实组织的绝对RNA质量。细胞数量比例另行保留，只作补充。两套真值均已从冻结的逐细胞清单重算，与既有文件一致，没有由signature或预测fraction反推真值。

这轮原始日志指定SC18_HoldoutNatural12_Mixture_CPM.txt、SC18_Signature18.txt、100次置换、QN关闭，且有All done。完整输出12×18，无缺失、无负数，各样本比例合计1。日志未显式列出batch correction模式，原提交计划为关闭，实际设置仍以网页Job Parameters为准；独立Warnings面板未提供。本地文件SHA256与冻结提交记录相符，但这不等同于核验网站保存文件的哈希。

两份原PDF各1页，Results显示值与TXT在三位小数精度一致。stackedChart中的Key柱是颜色图例，不是第13个样本。原图部分长样本标签被裁切，不影响数值审计；TXT和本报告为数值判断依据。

## 2. 与Job6联合参考的同样本对照

18个共同类型不重新归一化：Job10平均绝对误差为 **{summary['Job10_RNA_MAE18_pp']:.2f}** 个百分点，Job6为 **{summary['Job6_RNA_MAE18_pp']:.2f}**，约{summary['Job10_RNA_MAE18_pp'] / summary['Job6_RNA_MAE18_pp']:.2f}倍。12/12混合物Job10误差更大；仅EEC、Endothelial、Fibroblast三类的MAE较Job6低。

为保留Job6额外Paneth/Ileal_absorptive的误预测质量，另在完整20类轴补SCP259缺失两类为0：MAE分别{summary['Job10_RNA_MAE20_pp']:.2f}、{summary['Job6_RNA_MAE20_pp']:.2f}个百分点，平均总变差分别{summary['Job10_mean_TV']:.3f}、{summary['Job6_mean_TV']:.3f}。这两类在本轮测试中的真实贡献均为0；不是故意丢弃未建模的阳性细胞。总变差为两份组成的L1距离的一半，不代表同等比例的单细胞被错分。

{chr(10).join(scenario_table)}

Job6与测试共享GSE282122研究背景，虽供者留出，但相对于跨队列SCP259有域内优势。因此本轮支持“哪份参考更贴合当前GSE留出数据”，不能直接证明在真实MSCCR上仍保持相同排名。

## 3. 主要偏差及各类型完整表

Tuft均值从真值3.07%变成19.78%，ILC从2.50%变成11.04%，Stem/TA从16.47%变成25.40%。Goblet从14.43%降为1.32%，结肠吸收细胞从9.75%降为1.92%，单核/巨噬细胞从4.01%降为0.16%。CD8_T在全部12个真值阳性的混合物中输出0；单核/巨噬细胞9/12、Goblet 4/12为0。

这些是组成分配的系统性偏差，不能由此断言具体的Goblet细胞被逐一识别成Tuft；输出只有混合物比例，没有单细胞预测标签。类型均值容易掩盖样本差异，逐样本与216对真值/预测值已另存。

{chr(10).join(type_table)}

网页表达拟合Correlation范围{summary['web_expression_correlation_min']:.3f}–{summary['web_expression_correlation_max']:.3f}，中位数{summary['web_expression_correlation_median']:.3f}；全部P-value显示0。这些指标不能替代组成误差，尤其有限100次置换的P值0不是数学意义的概率为0，也不是每类比例准确的证据。

## 4. 粗谱系总量与内部类型分配

将Goblet、Stem/TA、EEC、Colonic_absorptive、Tuft加总，上皮真实均值{epithelium.RNA_truth.mean() * 100:.2f}%，Job10均值{epithelium.Job10.mean() * 100:.2f}%，逐样本总量MAE **{epithelium.Job10_abs_error_pp.mean():.2f}** 个百分点。分情景的趋势仍可辨认：

{chr(10).join(lineage_table)}

这说明粗上皮/免疫丰度趋势尚有信息，但上皮内部细分明显失真，不能因堆叠图随上皮富集而变化就接受18类精细比例。间质加总含Glial只是本审计的粗分组约定；胶质并非经典间质细胞，逐类型数据仍保留，不据粗加总赋予新生物学定义。

## 5. 输入检查和抽样提供的改进线索

signature共4,612个基因，测试覆盖4,381个，重叠 **{summary['signature_overlap_percent']:.2f}%**；不是整体基因交集过低。Goblet的MUC2、CLCA1、SPINK4、FCGBP、ZG16、TFF3在signature及测试中均存在，且参考中的最高类型都为Goblet；Tuft的POU2F3、TRPM5、AVIL、SH2D6也存在并指向Tuft。因此“完全没有对应marker”不足以解释此次失真，存在基因不等同于跨队列表达谱可迁移。

实际网页参考900个细胞、24位供者，每类50个。Goblet中43个Immature Goblet、7个Goblet；吸收细胞中34个Immature Enterocytes、11个Enterocytes、5个Best4+ Enterocytes；Stem/TA中18个TA1、17个TA2、11个Enterocyte Progenitors、3个Secretory TA，仅1个Stem。大类型总数平衡并没有保证内部亚型具有足够代表性。SCP259来源为UC/健康，没有CD；样本中的Epi/LP表示组织分离组分，不等同于精确结肠/直肠位置。

Stem/TA的LGR5、MKI67、TOP2A未进入最终signature，在既有输出推断参考中对应均值为0；OLFM4虽保留，但不能单独证明整个Stem/TA合并类准确。这里描述的是抽样和推断输出，不能把它当原始完整SCP259这些基因均不存在的证据。既有诊断中Stem/TA与吸收细胞原尺度列相关性较高，ILC/髓系部分共享基因亦缺乏严格类型特异性；应检查合并层级和高权重共享基因，而不是盲目增加每类marker数量。

SCP259的既有signature条件数约53.72，优于联合参考约73.02，但本轮预测更差：条件数、marker存在和基因覆盖只是必要诊断，不能代替留出准确性。

## 6. 后续判断与执行顺序

1. 先审阅已准备的SC18_HoldoutEqualRNA4结果，将每细胞RNA量差异与参考迁移问题区分。自然混合物已按RNA贡献评价，不宜把当前全部误差归因于细胞比例/RNA比例混用。
2. 如继续挽救SCP259，优先做其自身供者留出与跨队列对照，保持样本、真值及设置固定；再单独比较有明确匹配reference的S-mode设置。作者的方法说明将S-mode用于单细胞参考的跨平台校正，并不保证所有类型均可修复；参见[作者方法说明](https://pmc.ncbi.nlm.nih.gov/articles/PMC7695353/)及[原论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC6610714/)。单独改变一个因素，避免把校正、标签重分层、抽样同时改变后无法解释改善。
3. 如果重建参考，优先对Goblet成熟度、Stem与各TA/吸收前体组成设清楚的小规模覆盖约束，再检查CD8/NK/ILC及髓系的共同marker和归并边界；不要求每个罕见亚型严格等量，也不靠强行添加marker来保证可辨识性。平台、疾病和解剖部位差异都是待验证解释，当前数据未区分它们的作用。
4. 当前SC18保留为诊断对照，不直接用于正式MSCCR或下游EcoTyper输入。GSE联合参考本轮相对更好，但既往上皮亚型问题仍在，不能仅因击败SC18就宣布其完成验证。此处也未审阅cell-type expression/G矩阵，不能据fraction结果推定其适合presorted流程。

## 7. 文件与复查入口

- `paired_sample_accuracy_metrics.tsv`：12个混合物的配对误差、两种真值评价、网页拟合指标。
- `paired_cell_type_accuracy_metrics.tsv`：18类均值、MAE、偏差和假零数量。
- `all_216_paired_fraction_comparisons.tsv`：12×18全量逐项比较。
- `coarse_lineage_fractions.tsv`：粗谱系逐样本贡献，另附分情景平均表。
- `relevant_marker_coverage_and_specificity.tsv`、`SCP259_training_author_subtype_counts.tsv`：marker和抽样证据。
- `truth_vs_Job6_vs_Job10.png`、`Job10_fraction_bias.png`：统一轴的真值/预测和偏差图。B1–B4、E1–E4、I1–I4依次对应三个情景。
- `original_outputs/`：三个原始网页结果及运行日志的SHA核对副本；`review_summary.json`：机器可读摘要。
- 脚本：`scripts/python/review_cibersortx_job10_scp259_natural.py`；运行环境：`.tools/cibersortx_review`，包版本记录于`audit_runtime_versions.json`。

本轮全部在本机完成，未连接服务器、重抽供者、修改参考/真值或代提交网页。
'''
(out / 'Job10结直肠自然混合物审阅.md').write_text(report, encoding='utf-8')
print('Report:', out / 'Job10结直肠自然混合物审阅.md', flush=True)
