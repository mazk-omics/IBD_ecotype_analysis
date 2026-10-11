# 对1249基因运行做独立数学核对：Spearman、重拟合分界、弱类型NNLS和OLS。
# 本机科学Python环境运行；读取R导出的真实节点，由SciPy/NumPy独立求解，结果与表达分开保存。
from pathlib import Path
import itertools
import json
import numpy as np
import pandas as pd
from scipy.optimize import nnls

# Step1. 小无并列样例穷举8!个排列，核对AS89的小样本尾概率
project = Path(__file__).resolve().parents[2]
audit = project / 'CIBERSORTx_server_plan_20261011'
checks_dir = audit / 'full_target_math_checks'
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
    # 移动全零设计行不会改变理论RSS；仅这种结构性并列允许不同索引，保留严格匹配标记。
    rss = np.asarray(residual_squares)
    objective_gap = float(rss[recorded_cut - 30] - rss.min())
    numerical_tolerance = 10 * np.finfo(float).eps * max(1, np.sum(y**2))
    moved_rows = f[min(expected_cut, recorded_cut):max(expected_cut, recorded_cut)]
    structural_tie = expected_cut != recorded_cut and bool((moved_rows == 0).all()) and abs(objective_gap) <= numerical_tolerance
    split_checks.append({'GeneSymbol': gene, 'scipy_cut': expected_cut, 'r_cut': recorded_cut, 'strict_cut_equal': expected_cut == recorded_cut, 'structural_zero_design_tie': structural_tie, 'rss_gap': objective_gap, 'rss_roundoff_tolerance': float(numerical_tolerance), 'passed': expected_cut == recorded_cut or structural_tie})
    weak_frame = weak_inputs.loc[weak_inputs['GeneSymbol'] == gene]
    expected_beta, residual = nnls(weak_frame.iloc[:, 3:].to_numpy(), weak_frame['residual'].to_numpy())
    recorded = weak_coefficients.loc[weak_coefficients['GeneSymbol'] == gene]
    assert recorded['CellType'].tolist() == weak_frame.columns[3:].tolist()
    relative_error = float(np.max(np.abs(expected_beta - recorded['coefficient'].to_numpy())) / max(1, np.max(np.abs(expected_beta))))
    weak_checks.append({'GeneSymbol': gene, 'relative_coefficient_error': relative_error, 'passed': relative_error < 1e-10})
pd.DataFrame(split_checks).to_csv(checks_dir / 'refined_split_vs_scipy.tsv', sep='\t', index=False)
pd.DataFrame(weak_checks).to_csv(checks_dir / 'weak_residual_nnls_vs_scipy.tsv', sep='\t', index=False)
print('分界严格不同：', sum(not row['strict_cut_equal'] for row in split_checks), '；结构性并列：', sum(row['structural_zero_design_tie'] for row in split_checks))
print('弱类型系数不同：', [row for row in weak_checks if not row['passed']])
assert all(row['passed'] for row in split_checks) and all(row['passed'] for row in weak_checks)

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

# Step4. 汇总全目标集合的检查；通过数学核对不等于与当前官方软件等价
cor_checks = pd.read_csv(checks_dir / 'residual_spearman_vs_R_public.tsv', sep='\t')
shrink_checks = pd.read_csv(checks_dir / 'shrink_invariants.tsv', sep='\t')
assert cor_checks['passed'].all() and shrink_checks['passed'].all()
summary = {'AS89_toy_cases': len(toy_checks), 'exhaustive_permutations_per_case': len(permutations), 'residual_Spearman_checks': len(cor_checks), 'shrink_gene_checks': len(shrink_checks), 'scipy_split_checks': len(split_checks), 'strict_cut_equal': int(sum(item['strict_cut_equal'] for item in split_checks)), 'structural_zero_design_ties': int(sum(item['structural_zero_design_tie'] for item in split_checks)), 'scipy_weak_NNLS_checks': len(weak_checks), 'numpy_OLS_checks': len(ols_checks), 'OLS_max_relative_prediction_error': max(item['relative_prediction_error'] for item in ols_checks), 'all_objective_and_formula_checks_passed': True, 'limits': '12 cut indices differ across solvers for structural zero-design ties. Independent formula checks do not establish official version, RNG or upstream parity.'}
(checks_dir / 'full_target_math_summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')
print(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False))
