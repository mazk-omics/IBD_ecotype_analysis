# 用SciPy/NumPy检查Group真实抽样的NNLS、log中位数、CV、SD和解析t检验。
# 科学Python环境运行；24个均匀冻结基因共720次拟合，验证公式而非官方随机流。
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.optimize import nnls
from scipy.stats import t, false_discovery_control

# Step1. 按原始样本顺序重建每次抽样设计，不从官方输出读取系数
project = Path(__file__).resolve().parents[2]
audit = project / 'CIBERSORTx_server_plan_20261011'
group_dir = audit / 'group_independent_run1'
checks = audit / 'group_math_checks'
mixture = pd.read_csv(project / 'CIBERSORTx_ready_20261008/04_MSCCR120_HiRes_gene_batches/MSCCR_Platform120_CPM.txt', sep='\t', index_col=0)
fractions = pd.read_csv(audit / 'core_benchmark_inputs/fractions_fixed.tsv', sep='\t', index_col=0)
assert mixture.columns.tolist() == fractions.index.tolist()
indices = pd.read_csv(checks / 'bootstrap_indices.tsv', sep='\t')
coefficients = pd.read_csv(checks / 'bootstrap_coefficients.tsv', sep='\t')
F = fractions.to_numpy()
records = []
for gene, rows in indices.groupby('GeneSymbol', sort=False):
    expected = coefficients.loc[coefficients['GeneSymbol'] == gene, fractions.columns].to_numpy()
    predicted = []
    for row in rows.iloc[:, 2:].to_numpy(dtype=int):
        chosen = row - 1
        assert len(set(chosen)) == 90 and np.all((chosen >= 0) & (chosen < 120))
        predicted.append(nnls(F[chosen], mixture.loc[gene].to_numpy()[chosen])[0])
    predicted = np.asarray(predicted)
    relative = abs(predicted - expected) / np.maximum(1, abs(expected))
    records.append({'GeneSymbol': gene, 'NNLS_calls': len(predicted), 'maximum_relative_coefficient_difference': float(relative.max())})
    assert relative.max() < 1e-8
pd.DataFrame(records).to_csv(checks / 'scipy_NNLS_check.tsv', sep='\t', index=False)

# Step2. 用冻结R抽样系数复算统计，解析SE与bootstrap SD分开验收
tables = {name: pd.read_csv(group_dir / f'{name}.tsv', sep='\t', index_col=0) for name in ['bootstrap_beta_before_q_guard', 'GEPs_CVs', 'GEPs_StdErrs', 'analytic_SE', 'GEPs_Pvals', 'GEPs_Qvals', 'GEPs', 'GEPs_Filtered']}
covariance_basis = np.linalg.inv(F.T @ F)
statistic_records = []
for gene, rows in coefficients.groupby('GeneSymbol', sort=False):
    a = rows[fractions.columns].to_numpy()
    beta = np.exp2(np.median(np.log2(a + 1), axis=0)) - 1
    cv = np.sqrt(np.expm1(np.std(np.log1p(a), axis=0, ddof=1)**2))
    sd = np.std(a, axis=0, ddof=1)
    variance = np.sum((mixture.loc[gene].to_numpy() - F @ beta)**2) / 100
    se = np.sqrt(variance * np.diag(covariance_basis))
    with np.errstate(invalid='ignore', divide='ignore'):
        p = 2 * t.sf(beta / se, df=100)
    for name, value in [('bootstrap_beta_before_q_guard', beta), ('GEPs_CVs', cv), ('GEPs_StdErrs', sd), ('analytic_SE', se), ('GEPs_Pvals', p)]:
        expected = tables[name].loc[gene].to_numpy()
        assert np.allclose(value, expected, rtol=1e-10, atol=1e-10, equal_nan=True)
        statistic_records.append({'GeneSymbol': gene, 'statistic': name, 'max_absolute_error': float(np.nanmax(abs(value - expected))) if np.isfinite(value).any() else 0})
pd.DataFrame(statistic_records).to_csv(checks / 'numpy_statistics_check.tsv', sep='\t', index=False)

# Step3. 全基因BH和过滤就地核对，结果不能只在检查点范围进行多重校正
for celltype in fractions.columns:
    p = tables['GEPs_Pvals'][celltype].to_numpy()
    tested = np.isfinite(p)
    q = np.ones(len(p))
    q[tested] = false_discovery_control(p[tested], method='bh')
    assert np.allclose(q, tables['GEPs_Qvals'][celltype], rtol=1e-10, atol=1e-10)
thresholds = pd.read_csv(group_dir / 'group_cv_thresholds.tsv', sep='\t').set_index('CellType')
for celltype in fractions.columns:
    expected_na = tables['GEPs_CVs'][celltype].to_numpy() >= thresholds.loc[celltype, 'cv_threshold']
    assert np.array_equal(expected_na, tables['GEPs_Filtered'][celltype].isna())
strong = tables['GEPs_Qvals'].min(axis=1).to_numpy() < 1e-5
expected = tables['bootstrap_beta_before_q_guard'].to_numpy().copy()
expected[strong] = np.where(tables['GEPs_Qvals'].to_numpy()[strong] > .25, 0, expected[strong])
assert np.allclose(expected, tables['GEPs'].to_numpy(), rtol=1e-12, atol=1e-12)
summary = {'all_formula_checks_passed': True, 'NNLS_checkpoints': sum(r['NNLS_calls'] for r in records), 'checked_genes': len(records), 'maximum_relative_NNLS_difference': max(r['maximum_relative_coefficient_difference'] for r in records), 'statistic_checks': len(statistic_records), 'full_BH_and_filter_pairs': int(tables['GEPs'].size), 'limit': 'Formula consistency only; independent RNG/full residual differ from cached legacy and web outputs.'}
(checks / 'group_math_summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
print(json.dumps(summary, indent=2))
