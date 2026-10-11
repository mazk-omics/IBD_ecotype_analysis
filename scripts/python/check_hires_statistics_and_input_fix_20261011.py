# 独立核对原始输入下的R求解、完整基因BH和HiRes统计；比较输入修正前后的同一冻结基准。
# 本机 cibersortx_review 环境运行，读取新派生结果并保存检查表，不修改官方输入或输出。
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from scipy.optimize import nnls
from scipy.stats import chi2, false_discovery_control, norm, t

# Step1. 同一24基因、120样本、20类型；原始输入和调整后导出分开记录
project = Path(__file__).resolve().parents[2]
audit = project / 'CIBERSORTx_server_plan_20261011'
inputs = audit / 'core_benchmark_inputs'
core = audit / 'core_benchmark_original_run1'
entry_dir = audit / 'entry_statistics_original_run1'
worker_dir = audit / 'worker_statistics_original_run1'
outdir = audit / 'statistics_validation'
outdir.mkdir(exist_ok=False)
mixture = pd.read_csv(audit / 'input_space_audit/mixture_original_24genes.tsv', sep='\t', index_col=0)
fractions = pd.read_csv(inputs / 'fractions_fixed.tsv', sep='\t', index_col=0)
r_nnls = pd.read_csv(core / 'core_global_nnls.tsv', sep='\t', index_col=0)
worker = pd.read_csv(worker_dir / 'worker_type_statistics.tsv', sep='\t')
genes = pd.read_csv(worker_dir / 'worker_gene_statistics.tsv', sep='\t')
assert mixture.columns.tolist() == fractions.index.tolist() and mixture.index.tolist() == r_nnls.index.tolist()
assert genes['GeneSymbol'].tolist() == mixture.index.tolist() and len(worker) == 24 * 20

# Step2. 用SciPy独立NNLS核对系数和两段分界搜索
solver_rows = []
split_rows = []
fraction_array = fractions.to_numpy()
for gene in mixture.index:
    y = mixture.loc[gene].to_numpy()
    beta, residual = nnls(fraction_array, y)
    relative_error = float(np.max(np.abs(beta - r_nnls.loc[gene].to_numpy())) / max(1, np.max(np.abs(beta))))
    solver_rows.append({'GeneSymbol': gene, 'relative_coefficient_error': relative_error, 'passed': relative_error < 1e-8})
    sample_order = np.lexsort((np.arange(len(y)), y))
    sorted_y = y[sample_order]
    sorted_fractions = fraction_array[sample_order]
    cuts = np.arange(30, 91)
    rss = []
    for cut in cuts:
        beta_left, residual_left = nnls(sorted_fractions[:cut], sorted_y[:cut])
        beta_right, residual_right = nnls(sorted_fractions[cut:], sorted_y[cut:])
        rss.append(residual_left**2 + residual_right**2)
    scipy_cut = int(cuts[np.argmin(rss)])
    r_cut = int(genes.loc[genes['GeneSymbol'] == gene, 'best_cut'].iloc[0])
    split_rows.append({'GeneSymbol': gene, 'r_cut': r_cut, 'scipy_cut': scipy_cut, 'passed': r_cut == scipy_cut})
assert all(row['passed'] for row in solver_rows) and all(row['passed'] for row in split_rows)
pd.DataFrame(solver_rows).to_csv(outdir / 'r_vs_scipy_original_nnls.tsv', sep='\t', index=False)
pd.DataFrame(split_rows).to_csv(outdir / 'r_vs_scipy_split.tsv', sep='\t', index=False)

# Step3. 独立计算全基因入口p、BH；NA不当成零，不在24基因中重新BH
input_provenance = json.loads((audit / 'input_space_audit/input_space_provenance.json').read_text(encoding='utf-8'))
full_mixture = pd.read_csv(input_provenance['original_source'], sep='\t', index_col=0)
entry_beta = pd.read_csv(entry_dir / 'entry_bootstrap_medians.tsv', sep='\t', index_col=0)
entry_p = pd.read_csv(entry_dir / 'entry_raw_pvalues.tsv', sep='\t', index_col=0)
entry_q = pd.read_csv(entry_dir / 'entry_bh_qvalues.tsv', sep='\t', index_col=0)
worker_q = pd.read_csv(entry_dir / 'entry_worker_qvalues.tsv', sep='\t', index_col=0)
assert entry_beta.index.tolist() == full_mixture.index.tolist() == entry_p.index.tolist() == entry_q.index.tolist()
assert full_mixture.columns.tolist() == fractions.index.tolist() and entry_beta.columns.tolist() == fractions.columns.tolist()
residuals = full_mixture.to_numpy() - entry_beta.to_numpy() @ fraction_array.T
residual_variance = np.sum(residuals**2, axis=1) / (len(fractions) - fractions.shape[1])
covariance_diagonal = np.diag(np.linalg.inv(fraction_array.T @ fraction_array))
se = np.sqrt(residual_variance[:, None] * covariance_diagonal[None, :])
with np.errstate(divide='ignore', invalid='ignore'):
    expected_p = 2 * t.sf(entry_beta.to_numpy() / se, df=len(fractions) - fractions.shape[1])
assert np.array_equal(np.isnan(expected_p), np.isnan(entry_p.to_numpy()))
p_error = float(np.nanmax(np.abs(expected_p - entry_p.to_numpy())))
assert p_error < 1e-10
bh_rows = []
for celltype in entry_p.columns:
    valid = entry_p[celltype].notna()
    expected_q = false_discovery_control(entry_p.loc[valid, celltype].to_numpy(), method='bh')
    error = float(np.max(np.abs(expected_q - entry_q.loc[valid, celltype].to_numpy())))
    bh_rows.append({'CellType': celltype, 'tested_genes': int(valid.sum()), 'max_q_error': error, 'passed': error < 1e-12})
assert all(row['passed'] for row in bh_rows)
strong_gene = (entry_q < 1e-5).any(axis=1)
expected_worker_q = entry_q.copy()
expected_worker_q.loc[strong_gene] = entry_q.loc[strong_gene].mask(entry_q.loc[strong_gene] > 0.01, 1)
assert np.allclose(expected_worker_q, worker_q, equal_nan=True, atol=1e-12, rtol=1e-12)
pd.DataFrame(bh_rows).to_csv(outdir / 'full_gene_bh_vs_scipy.tsv', sep='\t', index=False)

# Step4. Fisher和Z分布使用独立库核对，不将数值规则误写成独立统计证据
fisher_q = worker['entry_q'].to_numpy().copy()
fisher_q[fisher_q == 0] = 1e-100
with np.errstate(divide='ignore', invalid='ignore'):
    expected_global_p = chi2.sf(-2 * (np.log(worker['probability_nonpositive']) + np.log(fisher_q)), df=4)
expected_global_p[np.isnan(fisher_q) | (fisher_q == 1)] = 1
global_error = float(np.max(np.abs(expected_global_p - worker['global_p'])))
expected_z_p = 2 * norm.sf(np.abs(worker['difference_z']))
expected_z_p[expected_z_p == 0] = 1e-100
z_error = float(np.max(np.abs(expected_z_p - worker['z_p'])))
tested = worker['empirical_p'].notna()
expected_differential_p = chi2.sf(-2 * (np.log(worker.loc[tested, 'empirical_p']) + np.log(worker.loc[tested, 'z_p'])), df=4)
differential_error = float(np.max(np.abs(expected_differential_p - worker.loc[tested, 'differential_p'])))
assert max(global_error, z_error, differential_error) < 1e-12

# Step5. 同一误差指标比较修正前后；记录仍未实现的完整流程范围
old_group = pd.read_csv(audit / 'core_benchmark_run1/core_vs_official_group.tsv', sep='\t')
new_group = pd.read_csv(core / 'core_vs_official_group.tsv', sep='\t')
comparison = old_group.merge(new_group, on='CellType', suffixes=('_adjusted_export', '_original_mixture'))
comparison.to_csv(outdir / 'input_fix_vs_official_group.tsv', sep='\t', index=False)
old_mae = float(old_group['bootstrap_median_log2_mae'].median())
new_mae = float(new_group['bootstrap_median_log2_mae'].median())
provenance = json.loads((inputs / 'input_provenance.json').read_text(encoding='utf-8'))
integrity = [{'path': item['path'], 'unchanged': hashlib.sha256(Path(item['path']).read_bytes()).hexdigest() == item['sha256']} for item in provenance['sources']]
integrity.append({'path': input_provenance['original_source'], 'unchanged': hashlib.sha256(Path(input_provenance['original_source']).read_bytes()).hexdigest() == input_provenance['original_sha256']})
assert all(item['unchanged'] for item in integrity)
pd.DataFrame(integrity).to_csv(outdir / 'source_integrity.tsv', sep='\t', index=False)
summary = {'status': 'input_correction_and_statistics_verified_not_full_hires', 'panel_genes': len(mixture), 'samples': len(fractions), 'types': fractions.shape[1], 'full_entry_genes': len(entry_p), 'full_entry_pvalue_error': p_error, 'full_gene_bh_max_error': max(row['max_q_error'] for row in bh_rows), 'nnls_genes_passed': len(solver_rows), 'nnls_max_relative_error': max(row['relative_coefficient_error'] for row in solver_rows), 'split_genes_passed': len(split_rows), 'worker_detectable_panel_genes': int(genes['gene_detectable'].sum()), 'worker_differential_pairs_below_0_05': int((worker['differential_p'] < 0.05).sum()), 'fisher_global_max_error': global_error, 'fisher_differential_max_error': differential_error, 'z_p_max_error': z_error, 'old_median_type_group_log2_mae': old_mae, 'corrected_median_type_group_log2_mae': new_mae, 'mae_reduction_percent': (1 - new_mae / old_mae) * 100, 'all_official_source_files_unchanged': True, 'unresolved': ['Current website input routing not proven by numerical comparison alone', 'Original-versus-adjusted gene export filtering', 'Independent RNG versus official threaded RNG', 'Residual Spearman fallback, shrinkage, targeted refit, OLS, rare-type overrides and final NA/floor rules', 'Signature building, fractions and S-mode independent implementations', 'No complete HiRes or official-equivalent EcoTyper input produced']}
(outdir / 'validation_summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps(summary, ensure_ascii=False, indent=2))
