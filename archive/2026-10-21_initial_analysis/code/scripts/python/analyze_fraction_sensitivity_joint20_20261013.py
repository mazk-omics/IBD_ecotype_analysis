# 配对真实120例比较一致性，并在精确恢复的自然28例上评价相对真值误差。
# 本机科学Python环境运行；只读取本轮BayesPrism和既有官方结果，输出数值表与患者bootstrap区间。
from pathlib import Path
import json
import sys
import numpy as np
import pandas as pd
from scipy import stats

# Step1. 按SampleID和20类冻结标签对齐，不重新归一化
project = Path(__file__).resolve().parents[2]
out = Path((project / '06_validation/fraction_sensitivity_current_path.txt').read_text(encoding='utf-8'))
types = (out / '00_audit/frozen_types.txt').read_text().splitlines()
annotation = pd.read_csv(out / '00_audit/MSCCR120_annotation.tsv', sep='\t').set_index('ID')
job20_path = project / 'CIBERSORTx_signature_audit/MSCCR_Joint20v2_Job20_review_20261009/original_outputs/CIBERSORTx_Job20_Adjusted.txt'
cx = pd.read_csv(job20_path, sep='\t', index_col=0).loc[annotation.index, types]
bp = pd.read_csv(out / '02_bayesprism_joint20_MSCCR120/theta_final.tsv', sep='\t', index_col=0).loc[annotation.index, types]
initial = pd.read_csv(out / '02_bayesprism_joint20_MSCCR120/theta_initial.tsv', sep='\t', index_col=0).loc[annotation.index, types]
assert bp.shape == (120, 20) and annotation.PatientID.nunique() == 120
assert np.isfinite(bp).all().all() and bp.ge(0).all().all() and np.allclose(bp.sum(axis=1), 1)
near_zero = 1e-4  # 0.01个百分点；两方法使用相同预定阈值。

# Step2. 主比较及分层，比例差异统一为百分点
comparisons = [('All120', bp, cx), ('Initial_vs_Job20', initial, cx)]
job17 = pd.read_csv(project / 'CIBERSORTx_signature_audit/MSCCR_Joint20v2_Job17_review_20261008/original_outputs/CIBERSORTx_Job17_Adjusted.txt', sep='\t', index_col=0).loc[annotation.index, types]
comparisons.append(('Final_vs_Job17', bp, job17))
for field in ['Disease', 'Inflammation', 'Region3']:
    for group, indices in annotation.groupby(field, observed=True).groups.items():
        comparisons.append((field + ':' + group, bp.loc[indices], cx.loc[indices]))
records = []
for group, a, b in comparisons:
    for cell_type in types:
        x, y = a[cell_type].to_numpy(), b[cell_type].to_numpy()
        delta = (x - y) * 100
        variable = np.std(x) > 0 and np.std(y) > 0
        r = stats.pearsonr(x, y).statistic if variable else np.nan
        rho = stats.spearmanr(x, y).statistic if variable else np.nan
        bias, sd = delta.mean(), delta.std(ddof=1)
        records.append({'group': group, 'cell_type': cell_type, 'n': len(x), 'Pearson_r': r, 'Spearman_rho': rho, 'BP_mean_percent': x.mean() * 100, 'CX_mean_percent': y.mean() * 100, 'mean_difference_pp': bias, 'median_absolute_difference_pp': np.median(np.abs(delta)), 'RMSE_pp': np.sqrt(np.mean(delta**2)), 'BA_bias_pp': bias, 'BA_lower_LoA_pp': bias - 1.96 * sd, 'BA_upper_LoA_pp': bias + 1.96 * sd, 'BP_zero_fraction': np.mean(x == 0), 'CX_zero_fraction': np.mean(y == 0), 'BP_near_zero_fraction': np.mean(x <= near_zero), 'CX_near_zero_fraction': np.mean(y <= near_zero)})
metrics = pd.DataFrame(records)

# Step3. 120个独立患者配对bootstrap；保持双方同一抽样索引
n_boot = 2000
rng = np.random.default_rng(20261013)
draw_blocks = []
for region in ['Colon', 'Rectum', 'Ileum']:
    region_indices = np.flatnonzero(annotation.Region3.eq(region))
    draw_blocks.append(rng.choice(region_indices, size=(n_boot, len(region_indices)), replace=True))
draws = np.concatenate(draw_blocks, axis=1)
x_boot = bp.to_numpy()[draws]
y_boot = cx.to_numpy()[draws]
dx = x_boot - x_boot.mean(axis=1, keepdims=True)
dy = y_boot - y_boot.mean(axis=1, keepdims=True)
pearson_boot = (dx * dy).sum(axis=1) / np.sqrt((dx**2).sum(axis=1) * (dy**2).sum(axis=1))
x_rank = stats.rankdata(x_boot, axis=1)
y_rank = stats.rankdata(y_boot, axis=1)
rx = x_rank - x_rank.mean(axis=1, keepdims=True)
ry = y_rank - y_rank.mean(axis=1, keepdims=True)
spearman_boot = (rx * ry).sum(axis=1) / np.sqrt((rx**2).sum(axis=1) * (ry**2).sum(axis=1))
delta_boot = (x_boot - y_boot) * 100
bootstrap_values = {'Pearson_r': pearson_boot, 'Spearman_rho': spearman_boot, 'mean_difference_pp': delta_boot.mean(axis=1), 'median_absolute_difference_pp': np.median(np.abs(delta_boot), axis=1), 'RMSE_pp': np.sqrt((delta_boot**2).mean(axis=1))}
for name, values in bootstrap_values.items():
    low, high = np.nanquantile(values, [0.025, 0.975], axis=0)
    metrics.loc[metrics.group.eq('All120'), name + '_CI_low'] = low
    metrics.loc[metrics.group.eq('All120'), name + '_CI_high'] = high
metrics.to_csv(out / '04_MSCCR120_concordance/type_concordance_metrics.tsv', sep='\t', index=False)
long = bp.stack().rename('BayesPrism_final').to_frame().join(cx.stack().rename('CIBERSORTx_Job20'))
long.index.names = ['sample', 'cell_type']
long = long.reset_index().merge(annotation.reset_index(), left_on='sample', right_on='ID', validate='many_to_one')
long['signed_difference_pp'] = (long.BayesPrism_final - long.CIBERSORTx_Job20) * 100
long.to_csv(out / '04_MSCCR120_concordance/paired_fractions_long.tsv', sep='\t', index=False)
delta = (bp - cx).abs()
sample_metrics = annotation.copy()
sample_metrics['total_variation'] = delta.sum(axis=1) / 2
sample_metrics['maximum_type_difference_pp'] = delta.max(axis=1) * 100
sample_metrics['largest_discrepancy_type'] = delta.idxmax(axis=1)
sample_metrics.sort_values('total_variation', ascending=False).to_csv(out / '04_MSCCR120_concordance/sample_discrepancies.tsv', sep='\t', index_label='sample')
main = metrics.loc[metrics.group.eq('All120')].copy()
main['abundance_group'] = np.where((main.BP_mean_percent + main.CX_mean_percent) / 2 >= 5, 'High_mean_ge5_percent', 'Low_mean_lt5_percent')
main.to_csv(out / '04_MSCCR120_concordance/abundance_and_lineage_summary.tsv', sep='\t', index=False)
real_summary = {'real_samples': 120, 'patients': 120, 'types': 20, 'renormalized': False, 'primary_theta': 'final', 'Pearson_median_across_types': float(main.Pearson_r.median()), 'Spearman_median_across_types': float(main.Spearman_rho.median()), 'median_sample_total_variation': float(sample_metrics.total_variation.median()), 'real_overall_mean_absolute_difference_pp': float(delta.to_numpy().mean() * 100), 'bootstrap_draws': n_boot, 'bootstrap_unit': 'patient, paired, region-stratified', 'near_zero_threshold': near_zero, 'natural28_primary_accuracy': []}
(out / '07_reports/numeric_summary.json').write_text(json.dumps(real_summary, indent=2), encoding='utf-8')
if not (out / '05_existing_pseudobulk_benchmark/theta_final.tsv').exists():
    print('MSCCR120_CONCORDANCE_COMPLETE; auxiliary natural28 still running')
    print(json.dumps(real_summary, indent=2))
    sys.exit(0)

# Step4. 自然28例按历史全基因UMI真值比较，同时报告细胞数与BP保留基因真值敏感性
truth_dir = out / '05_existing_pseudobulk_benchmark'
truth = pd.read_csv(truth_dir / 'ground_truth_export_UMI_fractions.tsv', sep='\t', index_col=0)[types]
cell_truth = pd.read_csv(truth_dir / 'ground_truth_cell_fractions.tsv', sep='\t', index_col=0).loc[truth.index, types]
bp_sim = pd.read_csv(truth_dir / 'theta_final.tsv', sep='\t', index_col=0).loc[truth.index, types]
bp_sim_initial = pd.read_csv(truth_dir / 'theta_initial.tsv', sep='\t', index_col=0).loc[truth.index, types]
cx_sim = pd.read_csv(project / 'CIBERSORTx_signature_audit/Joint20v2_Job15_review_20261007/original_outputs/CIBERSORTx_Job15_Results.txt', sep='\t', index_col=0).loc[truth.index, types]
retained_truth = pd.read_csv(truth_dir / 'ground_truth_BP_retained_gene_UMI_fractions.tsv', sep='\t', index_col=0).loc[truth.index, types]
accuracy_records = []
for truth_scale, target in [('Historical_export_UMI_primary', truth), ('Cell_count_secondary', cell_truth), ('BP_retained_gene_UMI_secondary', retained_truth)]:
    for method, prediction in [('BayesPrism_final', bp_sim), ('BayesPrism_initial', bp_sim_initial), ('CIBERSORTx_Job15', cx_sim)]:
        for cell_type in types + ['ALL_POOLED']:
            x = prediction.to_numpy().ravel() if cell_type == 'ALL_POOLED' else prediction[cell_type].to_numpy()
            y = target.to_numpy().ravel() if cell_type == 'ALL_POOLED' else target[cell_type].to_numpy()
            error = (x - y) * 100
            variable = np.std(x) > 0 and np.std(y) > 0
            accuracy_records.append({'truth_scale': truth_scale, 'method': method, 'cell_type': cell_type, 'n_mixtures': 28, 'Pearson_r': stats.pearsonr(x, y).statistic if variable else np.nan, 'Spearman_rho': stats.spearmanr(x, y).statistic if variable else np.nan, 'MAE_pp': np.abs(error).mean(), 'RMSE_pp': np.sqrt((error**2).mean()), 'mean_signed_error_pp': error.mean()})
accuracy = pd.DataFrame(accuracy_records)
accuracy.to_csv(truth_dir / 'natural28_accuracy_metrics.tsv', sep='\t', index=False)
bp_error = (bp_sim - truth).abs() * 100
cx_error = (cx_sim - truth).abs() * 100
paired_error = bp_error - cx_error
paired_error.to_csv(truth_dir / 'paired_absolute_error_difference_pp.tsv', sep='\t', index_label='sample')
design = pd.read_csv(truth_dir / 'mixture_design.tsv', sep='\t').set_index('sample')
scenario_records = []
for (region, scenario), ids in design.groupby(['region', 'scenario'], observed=True).groups.items():
    scenario_records.append({'region': region, 'scenario': scenario, 'n': len(ids), 'BP_MAE_pp': bp_error.loc[ids].to_numpy().mean(), 'CX_MAE_pp': cx_error.loc[ids].to_numpy().mean(), 'paired_difference_pp': paired_error.loc[ids].to_numpy().mean()})
pd.DataFrame(scenario_records).to_csv(truth_dir / 'scenario_accuracy_metrics.tsv', sep='\t', index=False)
# 模拟跨样本复用细胞，普通28样本bootstrap不构成供者级区间；只给配对误差分布和按场景描述。
pair_summary = pd.DataFrame({'cell_type': types, 'BP_MAE_pp': bp_error.mean(), 'CX_MAE_pp': cx_error.mean(), 'mean_paired_absolute_error_difference_pp': paired_error.mean(), 'median_paired_absolute_error_difference_pp': paired_error.median(), 'BP_better_mixtures': (paired_error < 0).sum()}).reset_index(drop=True)
pair_summary.to_csv(truth_dir / 'paired_type_accuracy_difference.tsv', sep='\t', index=False)
summary = {'real_samples': 120, 'patients': 120, 'types': 20, 'renormalized': False, 'primary_theta': 'final', 'Pearson_median_across_types': float(main.Pearson_r.median()), 'Spearman_median_across_types': float(main.Spearman_rho.median()), 'median_sample_total_variation': float(sample_metrics.total_variation.median()), 'real_overall_mean_absolute_difference_pp': float(delta.to_numpy().mean() * 100), 'bootstrap_draws': n_boot, 'bootstrap_unit': 'patient, paired, region-stratified', 'near_zero_threshold': near_zero, 'simulation_samples': 28, 'equalRNA8_raw_count_pair': 'BLOCKED', 'simulation_independent_external_validation': False, 'natural28_primary_accuracy': accuracy.loc[accuracy.truth_scale.eq('Historical_export_UMI_primary') & accuracy.cell_type.eq('ALL_POOLED')].to_dict('records')}
(out / '07_reports/numeric_summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
print(json.dumps(summary, indent=2))
