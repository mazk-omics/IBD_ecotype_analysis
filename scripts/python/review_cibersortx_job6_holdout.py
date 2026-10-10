# 审阅 Job6 留出测试：把预测比例按样本 ID 和类型名称匹配到已知细胞、RNA 贡献比例。
# 在本机项目环境运行；读取冻结的 Job6 输出和留出真值，写入 Job6_review，不改原 signature。
# 不将网页的表达拟合 Correlation/RMSE 当成比例准确率；同一细胞可跨模拟样本复用，不做独立患者推断。
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

# Step1. 按名称匹配结果，核对行列、比例范围和归一化
project = Path('D:/IBD_ecotype/IBD_ecotype_analysis')
base = project / 'CIBERSORTx_signature_audit/main_candidate_validation_20261006'
out = base / 'Job6_review'
truth_dir = base / 'server_outputs/validation_joint20_20261006'
results = pd.read_csv(out / 'original_outputs/CIBERSORTx_Job6_Results.txt', sep='\t', index_col=0)
cell_truth = pd.read_csv(truth_dir / 'ground_truth_cell_fractions.tsv', sep='\t', index_col=0)
RNA_truth = pd.read_csv(truth_dir / 'ground_truth_export_UMI_fractions.tsv', sep='\t', index_col=0)
original_RNA_truth = pd.read_csv(truth_dir / 'ground_truth_original_UMI_fractions.tsv', sep='\t', index_col=0)
design = pd.read_csv(truth_dir / 'mixture_design.tsv', sep='\t').set_index('sample').loc[cell_truth.index]
assert results.index.is_unique and set(results.index) == set(cell_truth.index)
assert set(results.columns) == set(cell_truth.columns) | {'P-value', 'Correlation', 'RMSE'}
predicted = results.loc[cell_truth.index, cell_truth.columns]
RNA_truth = RNA_truth.loc[cell_truth.index, cell_truth.columns]
assert predicted.shape == (28, 20) and np.isfinite(predicted.to_numpy()).all()
assert predicted.ge(0).all().all() and predicted.le(1).all().all() and np.allclose(predicted.sum(axis=1), 1)
assert np.allclose(cell_truth.sum(axis=1), 1) and np.allclose(RNA_truth.sum(axis=1), 1)
assert np.allclose(RNA_truth, original_RNA_truth.loc[cell_truth.index, cell_truth.columns])
predicted.to_csv(out / 'predicted_fractions_aligned.tsv', sep='\t')
print('匹配成功:', predicted.shape, '结果表含完整 28 样本；各列按类型名称对应', flush=True)

# Step2. 逐样本和逐类型计算绝对偏差，单位为百分点
cell_error = predicted - cell_truth
RNA_error = predicted - RNA_truth
per_sample = design[['region', 'scenario', 'absent_type', 'n_cells', 'n_donors']].copy()
per_sample['cell_MAE_pp'] = cell_error.abs().mean(axis=1) * 100
per_sample['RNA_MAE_pp'] = RNA_error.abs().mean(axis=1) * 100
per_sample['cell_RMSE_pp'] = np.sqrt(cell_error.pow(2).mean(axis=1)) * 100
per_sample['RNA_RMSE_pp'] = np.sqrt(RNA_error.pow(2).mean(axis=1)) * 100
# TV 是两组成分布的总变差距离，避免平均到 20 列后掩盖集中在 Paneth 的大偏差。
per_sample['cell_TV'] = cell_error.abs().sum(axis=1) / 2
per_sample['RNA_TV'] = RNA_error.abs().sum(axis=1) / 2
per_sample['fraction_pearson_cell'] = predicted.corrwith(cell_truth, axis=1)
per_sample['fraction_pearson_RNA'] = predicted.corrwith(RNA_truth, axis=1)
per_sample[['web_P_value', 'web_expression_correlation', 'web_expression_RMSE']] = results.loc[cell_truth.index, ['P-value', 'Correlation', 'RMSE']].to_numpy()
per_sample.to_csv(out / 'sample_accuracy_metrics.tsv', sep='\t')
per_type = pd.DataFrame({'cell_truth_mean': cell_truth.mean(), 'RNA_truth_mean': RNA_truth.mean(), 'predicted_mean': predicted.mean(),
    'cell_MAE_pp': cell_error.abs().mean() * 100, 'RNA_MAE_pp': RNA_error.abs().mean() * 100,
    'cell_bias_pp': cell_error.mean() * 100, 'RNA_bias_pp': RNA_error.mean() * 100,
    'fraction_pearson_cell': predicted.corrwith(cell_truth), 'fraction_pearson_RNA': predicted.corrwith(RNA_truth),
    'fraction_spearman_cell': predicted.rank().corrwith(cell_truth.rank()), 'fraction_spearman_RNA': predicted.rank().corrwith(RNA_truth.rank()),
    'zero_predictions': predicted.eq(0).sum(), 'true_present_samples': cell_truth.gt(0).sum(),
    'predicted_zero_when_present': (predicted.eq(0) & cell_truth.gt(0)).sum()})
per_type.index.name = 'cell_type'
per_type.to_csv(out / 'cell_type_accuracy_metrics.tsv', sep='\t')
region_metrics = per_sample.groupby('region')[['cell_MAE_pp', 'RNA_MAE_pp', 'cell_TV', 'RNA_TV', 'fraction_pearson_cell', 'fraction_pearson_RNA']].mean()
region_metrics['n_samples'] = design.groupby('region').size()
region_metrics.to_csv(out / 'region_accuracy_metrics.tsv', sep='\t')
print(region_metrics.to_string(), flush=True)

# Step3. 保存每个样本每类的完整对照，缺失与稀有类型不混成一个总体准确率
long_table = cell_truth.rename_axis(index='sample', columns='cell_type').stack().rename('cell_truth').reset_index()
rna_long = RNA_truth.rename_axis(index='sample', columns='cell_type').stack().rename('RNA_truth').reset_index()
prediction_long = predicted.rename_axis(index='sample', columns='cell_type').stack().rename('predicted').reset_index()
long_table = long_table.merge(rna_long, on=['sample', 'cell_type'], validate='one_to_one')
long_table = long_table.merge(prediction_long, on=['sample', 'cell_type'], validate='one_to_one')
long_table = long_table.merge(design[['region', 'scenario']].reset_index(), on='sample', validate='many_to_one')
long_table.to_csv(out / 'all_560_fraction_comparisons.tsv', sep='\t', index=False)
absent_records = []
for sample, row in design.loc[design.absent_type.notna()].iterrows():
    name = row.absent_type
    assert cell_truth.loc[sample, name] == 0
    absent_records.append({'sample': sample, 'cell_type': name, 'cell_truth': 0, 'RNA_truth': 0, 'predicted': predicted.loc[sample, name], 'predicted_percent': predicted.loc[sample, name] * 100})
absent = pd.DataFrame(absent_records)
absent.to_csv(out / 'designed_absent_type_predictions.tsv', sep='\t', index=False)
rare = long_table.loc[long_table.cell_truth.gt(0) & long_table.cell_truth.le(0.01)].copy()  # 本项目描述性分组：实际细胞占比不超过 1%
rare.to_csv(out / 'rare_type_comparisons_up_to_1percent.tsv', sep='\t', index=False)
region_type = long_table.groupby(['region', 'cell_type']).agg(cell_truth=('cell_truth', 'mean'), RNA_truth=('RNA_truth', 'mean'), predicted=('predicted', 'mean'))
region_type.to_csv(out / 'region_cell_type_mean_fractions.tsv', sep='\t')

# Step4. 在原训练细胞中核对 Paneth 亚型表达，区别实测证据与因果推测
training_dir = project / 'IBDverse_reference/three_reference_pilots_20261006/run_20261006_183211/gse_joint'
training_cells = pd.read_csv(training_dir / 'selected_cell_manifest.tsv', sep='\t')
training_expression = pd.read_csv(training_dir / 'GSE282122_Joint20_scRNA_reference.txt', sep='\t', index_col=0)
assert np.array_equal(training_expression.columns, training_cells.web_column)
markers = ['DEFA5', 'DEFA6', 'REG3A', 'REG1A', 'LYZ', 'MUC2', 'TFF3', 'LGR5', 'OLFM4']
marker_counts = training_expression.loc[markers].T
marker_counts.index = training_cells.index
marker_CPM = marker_counts.div(training_cells.total_counts, axis=0) * 1e6
marker_CPM['pilot_type'] = training_cells.pilot_type
marker_CPM['author_label'] = training_cells.author_label
paneth_marker_CPM = marker_CPM.loc[training_cells.pilot_type.eq('Paneth')].groupby('author_label')[markers].mean()
paneth_marker_CPM.to_csv(out / 'training_Paneth_subtype_mean_marker_CPM.tsv', sep='\t')
training_cells.loc[training_cells.pilot_type.eq('Paneth')].groupby('author_label').agg(n_cells=('cell_id', 'size'), mean_total_counts=('total_counts', 'mean'), median_total_counts=('total_counts', 'median')).to_csv(out / 'training_Paneth_subtype_library_sizes.tsv', sep='\t')
marker_CPM.groupby('pilot_type')[markers].mean().to_csv(out / 'training_type_mean_marker_CPM.tsv', sep='\t')
del training_expression
print('Paneth 亚型 marker 差异（原始 counts / 每细胞完整库总量）：', flush=True)
print(paneth_marker_CPM[['DEFA5', 'DEFA6', 'REG3A']].to_string(), flush=True)

# Step5. 保存校验与结论，保留原始输出用于后续追溯
files = {}
for file in (out / 'original_outputs').iterdir():
    files[file.name] = {'bytes': file.stat().st_size, 'sha256': hashlib.sha256(file.read_bytes()).hexdigest()}
summary = {'status': 'review_completed', 'job': 'Job6', 'samples': 28, 'cell_types': 20,
    'cell_MAE_pp': float(cell_error.abs().to_numpy().mean() * 100), 'RNA_MAE_pp': float(RNA_error.abs().to_numpy().mean() * 100),
    'web_expression_correlation_median': float(results.Correlation.median()), 'web_P_values_reported_zero': int(results['P-value'].eq(0).sum()),
    'web_metrics_are_not_fraction_accuracy': True, 'all_fractions_finite_nonnegative_and_sum_to_one': True,
    'ileum_Paneth_mean_cell_percent': float(cell_truth.loc[design.region.eq('Ileum'), 'Paneth'].mean() * 100),
    'ileum_Paneth_mean_RNA_percent': float(RNA_truth.loc[design.region.eq('Ileum'), 'Paneth'].mean() * 100),
    'ileum_Paneth_mean_predicted_percent': float(predicted.loc[design.region.eq('Ileum'), 'Paneth'].mean() * 100),
    'decision': 'current uncorrected Joint20 signature not accepted as validated whole-cohort primary reference; investigate ileal Paneth scale and subtype composition first',
    'decision_is_project_judgment_not_predefined_statistical_pass_threshold': True,
    'Paneth_training_composition': {'Ileal paneth': 8, 'Non ileal paneth': 42},
    'equal_RNA_control_results_available': False, 'solver_convergence_confirmed': False,
    'previous_UI_max_iterations_warning_absent_from_downloaded_log': True,
    'MSCCR_deconvolution_run': False, 'clinical_analysis_performed': False, 'original_outputs': files}
(out / 'review_summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)
