# 独立验证收缩阶段的分界、弱类型NNLS与OLS；分别比较官方可变、恒定和NA行。
# 本机科学Python环境运行，读取冻结面板和新R结果，保存完整差异及阶段种子敏感性。
from pathlib import Path
import hashlib
import itertools
import json
import numpy as np
import pandas as pd
from scipy.optimize import nnls

# Step1. 小无并列样例穷举8!个排列，核对AS89的小样本尾概率
project = Path(__file__).resolve().parents[2]
audit = project / 'CIBERSORTx_server_plan_20261011'
inputs = audit / 'core_benchmark_inputs'
checks_dir = audit / 'refinement_validation'
baseline = audit / 'refinement_original_run2'
toys = pd.read_csv(checks_dir / 'spearman_toys_for_enumeration.tsv', sep='\t')
permutations = np.asarray(list(itertools.permutations(range(1, 9))), dtype=float)
null_s = np.sum((permutations - np.arange(1, 9))**2, axis=1)
toy_checks = []
for row in toys.itertuples(index=False):
    y = np.asarray([int(value) for value in row.y_order.split(',')])
    observed_s = int(np.sum((y - np.arange(1, 9))**2))
    one_tail = np.mean(null_s <= observed_s) if observed_s <= (8**3 - 8) / 6 else np.mean(null_s >= observed_s)
    expected_p = min(2 * one_tail, 1)
    error = abs(expected_p - row.two_sided_p)
    toy_checks.append({'example': row.example, 'enumerated_permutations': len(permutations), 'expected_p': expected_p, 'error': error, 'passed': error < 1e-12})
assert all(row['passed'] for row in toy_checks)
pd.DataFrame(toy_checks).to_csv(checks_dir / 'AS89_vs_exhaustive_permutations.tsv', sep='\t', index=False)

# Step2. 在真实可检测类型屏蔽后，SciPy独立搜索新分界并重估弱类型残差
refit_inputs = pd.read_csv(checks_dir / 'refit_inputs.tsv', sep='\t')
weak_inputs = pd.read_csv(checks_dir / 'weak_inputs.tsv', sep='\t')
weak_coefficients = pd.read_csv(checks_dir / 'weak_coefficients.tsv', sep='\t')
split_checks = []
weak_checks = []
for gene, frame in refit_inputs.groupby('GeneSymbol', sort=False):
    f = frame.iloc[:, 4:].to_numpy()
    y = frame['mixture'].to_numpy()
    cuts = np.arange(30, 91)
    residual_squares = []
    for cut in cuts:
        beta1, residual1 = nnls(f[:cut], y[:cut])
        beta2, residual2 = nnls(f[cut:], y[cut:])
        residual_squares.append(residual1**2 + residual2**2)
    expected_cut = int(cuts[np.argmin(residual_squares)])
    recorded_cut = int(frame['best_cut'].iloc[0])
    split_checks.append({'GeneSymbol': gene, 'scipy_cut': expected_cut, 'r_cut': recorded_cut, 'passed': expected_cut == recorded_cut})
    weak_frame = weak_inputs.loc[weak_inputs['GeneSymbol'] == gene]
    expected_beta, residual = nnls(weak_frame.iloc[:, 3:].to_numpy(), weak_frame['residual'].to_numpy())
    recorded = weak_coefficients.loc[weak_coefficients['GeneSymbol'] == gene]
    assert recorded['CellType'].tolist() == weak_frame.columns[3:].tolist()
    relative_error = float(np.max(np.abs(expected_beta - recorded['coefficient'].to_numpy())) / max(1, np.max(np.abs(expected_beta))))
    weak_checks.append({'GeneSymbol': gene, 'relative_coefficient_error': relative_error, 'passed': relative_error < 1e-10})
assert all(row['passed'] for row in split_checks) and all(row['passed'] for row in weak_checks)
pd.DataFrame(split_checks).to_csv(checks_dir / 'refined_split_vs_scipy.tsv', sep='\t', index=False)
pd.DataFrame(weak_checks).to_csv(checks_dir / 'weak_residual_nnls_vs_scipy.tsv', sep='\t', index=False)

# Step3. 同一x与锚点y用独立最小二乘验证OLS，不只比较相关性
calibration_inputs = pd.read_csv(checks_dir / 'calibration_inputs.tsv', sep='\t')
ols_checks = []
for (gene, celltype), frame in calibration_inputs.groupby(['GeneSymbol', 'CellType'], sort=False):
    x = frame['smooth_x'].to_numpy()
    target = frame['target_y'].to_numpy()
    if np.std(x) == 0:
        beta = np.asarray([np.mean(target), 0.0])
    else:
        beta = np.linalg.lstsq(np.column_stack([np.ones(len(x)), x]), target, rcond=None)[0]
    recorded = np.asarray([frame['r_intercept'].iloc[0], frame['r_slope'].iloc[0]])
    coefficient_error = float(np.max(np.abs(beta - recorded)) / max(1, np.max(np.abs(beta))))
    prediction_error = float(np.max(np.abs((beta - recorded)[0] + (beta - recorded)[1] * x)) / max(1, np.max(np.abs(target))))
    ols_checks.append({'GeneSymbol': gene, 'CellType': celltype, 'relative_coefficient_error': coefficient_error, 'relative_prediction_error': prediction_error, 'passed': max(coefficient_error, prediction_error) < 1e-10})
assert all(row['passed'] for row in ols_checks)
pd.DataFrame(ols_checks).to_csv(checks_dir / 'OLS_vs_numpy_lstsq.tsv', sep='\t', index=False)

# Step4. 同一官方26个可变组合比较前后；NA与恒定行分别记录，不能人为套官方掩码
group = pd.read_csv(inputs / 'official_group_24genes.tsv', sep='\t', index_col=0)
pair_records = []
for celltype in group.columns:
    official = pd.read_csv(inputs / f'official_hires_{celltype}.tsv', sep='\t', index_col=0)
    core = pd.read_csv(audit / 'core_benchmark_original_run1' / f'CORE_ONLY_{celltype}.tsv', sep='\t', index_col=0)
    refined = pd.read_csv(baseline / f'REFINEMENT_ONLY_{celltype}.tsv', sep='\t', index_col=0)
    assert official.index.tolist() == core.index.tolist() == refined.index.tolist()
    assert official.columns.tolist() == refined.columns.tolist() == core.columns.tolist()
    for gene in official.index:
        y = official.loc[gene].to_numpy()
        x_core = core.loc[gene].to_numpy()
        x_refined = refined.loc[gene].to_numpy()
        finite = np.isfinite(y)
        if not finite.any():
            official_status = 'NA_all_samples'
        elif not finite.all():
            official_status = 'partial_NA'
        elif np.ptp(y) > 1e-8:
            official_status = 'variable'
        else:
            official_status = 'constant'
        predicted_variable = bool(np.ptp(x_refined) > 1e-8)
        core_mae = float(np.mean(np.abs(np.log2(x_core[finite] + 1) - np.log2(y[finite] + 1)))) if finite.any() else None
        refined_mae = float(np.mean(np.abs(np.log2(x_refined[finite] + 1) - np.log2(y[finite] + 1)))) if finite.any() else None
        core_cor = None
        refined_cor = None
        if official_status == 'variable':
            if np.ptp(x_core) > 1e-8:
                core_cor = float(np.corrcoef(np.log2(x_core + 1), np.log2(y + 1))[0, 1])
            if predicted_variable:
                refined_cor = float(np.corrcoef(np.log2(x_refined + 1), np.log2(y + 1))[0, 1])
        pair_records.append({'GeneSymbol': gene, 'CellType': celltype, 'official_status': official_status, 'predicted_variable': predicted_variable, 'core_log2_mae': core_mae, 'refined_log2_mae': refined_mae, 'core_log2_correlation': core_cor, 'refined_log2_correlation': refined_cor, 'max_absolute_error': float(np.max(np.abs(x_refined[finite] - y[finite]))) if finite.any() else None})
pairs = pd.DataFrame(pair_records)
pairs.to_csv(checks_dir / 'all_480_pairs_vs_official.tsv', sep='\t', index=False)
variable = pairs.loc[pairs['official_status'] == 'variable'].copy()
variable.sort_values('refined_log2_mae', ascending=False).to_csv(checks_dir / 'official_variable_pairs_prioritized.tsv', sep='\t', index=False)
pd.crosstab(pairs['official_status'], pairs['predicted_variable']).to_csv(checks_dir / 'official_status_vs_predicted_variation.tsv', sep='\t')

# Step5. 保留原始统计，仅变收缩阶段种子；检测残余误差是否只是重抽样造成
seed_records = []
for run in ['refinement_original_seed2', 'refinement_original_seed3']:
    for celltype in group.columns:
        reference = pd.read_csv(baseline / f'REFINEMENT_ONLY_{celltype}.tsv', sep='\t', index_col=0)
        alternate = pd.read_csv(audit / run / f'REFINEMENT_ONLY_{celltype}.tsv', sep='\t', index_col=0)
        assert reference.index.tolist() == alternate.index.tolist() and reference.columns.tolist() == alternate.columns.tolist()
        relevant = variable.loc[variable['CellType'] == celltype, 'GeneSymbol']
        for gene in relevant:
            official = pd.read_csv(inputs / f'official_hires_{celltype}.tsv', sep='\t', index_col=0).loc[gene].to_numpy()
            prediction = alternate.loc[gene].to_numpy()
            correlation = float(np.corrcoef(np.log2(prediction + 1), np.log2(official + 1))[0, 1]) if np.ptp(prediction) > 1e-8 else None
            seed_records.append({'run': run, 'GeneSymbol': gene, 'CellType': celltype, 'changed_type_variation': bool((np.ptp(reference.loc[gene]) > 1e-8) != (np.ptp(prediction) > 1e-8)), 'mean_log2_difference_from_baseline': float(np.mean(np.abs(np.log2(reference.loc[gene] + 1) - np.log2(prediction + 1)))), 'log2_mae_vs_official': float(np.mean(np.abs(np.log2(prediction + 1) - np.log2(official + 1)))), 'log2_correlation_vs_official': correlation})
seed_comparison = pd.DataFrame(seed_records)
seed_comparison.to_csv(checks_dir / 'refinement_seed_sensitivity.tsv', sep='\t', index=False)

# Step6. 哈希核对原件；保存运行范围、真实缺口和全部不一致组合
provenance = json.loads((inputs / 'input_provenance.json').read_text(encoding='utf-8'))
integrity = [{'path': item['path'], 'unchanged': hashlib.sha256(Path(item['path']).read_bytes()).hexdigest() == item['sha256']} for item in provenance['sources']]
assert all(row['unchanged'] for row in integrity)
pd.DataFrame(integrity).to_csv(checks_dir / 'source_integrity.tsv', sep='\t', index=False)
cor_checks = pd.read_csv(checks_dir / 'residual_spearman_vs_R_public.tsv', sep='\t')
worst_pairs = variable.nlargest(5, 'refined_log2_mae')[['GeneSymbol', 'CellType', 'refined_log2_mae', 'refined_log2_correlation']]
worst_pairs = worst_pairs.astype(object).where(pd.notna(worst_pairs), None)
summary = {'status': 'panel_refinement_validated_not_complete_official_hires', 'genes': 24, 'samples': 120, 'types': 20, 'nnls_fits_per_refinement_run': 15288, 'AS89_toy_cases_passed': len(toy_checks), 'enumerated_permutations_per_case': 40320, 'residual_spearman_R_checks_passed': int(cor_checks['passed'].sum()), 'refined_split_scipy_genes_passed': len(split_checks), 'weak_refit_scipy_genes_passed': len(weak_checks), 'OLS_pairs_passed': len(ols_checks), 'OLS_max_relative_prediction_error': max(row['relative_prediction_error'] for row in ols_checks), 'official_variable_pairs': len(variable), 'official_constant_pairs': int((pairs['official_status'] == 'constant').sum()), 'official_all_NA_pairs': int((pairs['official_status'] == 'NA_all_samples').sum()), 'refined_variable_official_variable_pairs': int(variable['predicted_variable'].sum()), 'median_variable_pair_core_log2_mae': float(variable['core_log2_mae'].median()), 'median_variable_pair_refined_log2_mae': float(variable['refined_log2_mae'].median()), 'median_variable_pair_core_log2_correlation': float(variable['core_log2_correlation'].median()), 'median_variable_pair_refined_log2_correlation_among_varying_predictions': float(variable['refined_log2_correlation'].median()), 'pairs_correlation_at_least_0_95': int((variable['refined_log2_correlation'] >= 0.95).sum()), 'official_variable_but_predicted_constant': variable.loc[~variable['predicted_variable'], ['GeneSymbol', 'CellType']].to_dict('records'), 'worst_5_by_log2_mae': worst_pairs.to_dict('records'), 'refinement_seed_changes_of_variation_status': int(seed_comparison['changed_type_variation'].sum()), 'official_source_files_unchanged': True, 'unresolved': ['Current website NA/reliability mask and optional range filter', 'Auxiliary Group-mode constant refill', 'Old numerical-version interpretation versus current website version', 'Full official random stream and upstream statistical randomness', 'Independent signature, fractions and S-mode'], 'clinical_use': 'Research prototype only; no formal discovery or recovery input'}
(checks_dir / 'refinement_validation_summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')
print(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False))
