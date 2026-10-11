# 用 SciPy 独立求解同一24基因 NNLS，核对 R 数值原型的求解与真实数据最优性条件。
# 本机 cibersortx_review 环境运行；保存比较、源文件完整性和阶段耗时，不宣称官方 HiRes 等价。
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from scipy.optimize import nnls

# Step1. 读取固定基准和 R 输出，核对接口顺序
project = Path(__file__).resolve().parents[2]
audit = project / 'CIBERSORTx_server_plan_20261011'
inputs = audit / 'core_benchmark_inputs'
outputs = audit / 'core_benchmark_run1'
mixture = pd.read_csv(inputs / 'mixture_adjusted_24genes.tsv', sep='\t', index_col=0)
fractions = pd.read_csv(inputs / 'fractions_fixed.tsv', sep='\t', index_col=0)
reference = pd.read_csv(outputs / 'core_global_nnls.tsv', sep='\t', index_col=0)
assert reference.index.tolist() == mixture.index.tolist() and reference.columns.tolist() == fractions.columns.tolist()
assert mixture.columns.tolist() == fractions.index.tolist()

# Step2. 不同语言求解器比较系数与残差；满秩才比较系数唯一性
rank = int(np.linalg.matrix_rank(fractions.to_numpy()))
rows = []
for gene in mixture.index:
    predicted, residual = nnls(fractions.to_numpy(), mixture.loc[gene].to_numpy())
    r_predicted = reference.loc[gene].to_numpy()
    scale = max(1.0, float(np.max(np.abs(predicted))))
    coefficient_error = float(np.max(np.abs(predicted - r_predicted)))
    r_residual = float(np.linalg.norm(mixture.loc[gene].to_numpy() - fractions.to_numpy() @ r_predicted))
    rows.append({'GeneSymbol': gene, 'max_absolute_coefficient_error': coefficient_error, 'relative_coefficient_error': coefficient_error / scale, 'residual_norm_error': abs(residual - r_residual), 'passed': coefficient_error / scale < 1e-8 and abs(residual - r_residual) / max(1, residual) < 1e-8})
comparison = pd.DataFrame(rows)
assert rank == fractions.shape[1] and comparison['passed'].all()
comparison.to_csv(outputs / 'r_vs_scipy_nnls.tsv', sep='\t', index=False)

# Step3. 确认未修改官方文件，记录真正运行过的范围与缺口
provenance = json.loads((inputs / 'input_provenance.json').read_text(encoding='utf-8'))
integrity = [{'path': item['path'], 'unchanged': hashlib.sha256(Path(item['path']).read_bytes()).hexdigest() == item['sha256']} for item in provenance['sources']]
assert all(item['unchanged'] for item in integrity)
pd.DataFrame(integrity).to_csv(outputs / 'source_integrity_after_run.tsv', sep='\t', index=False)
timing = pd.read_csv(outputs / 'timings.tsv', sep='\t')
kkt = pd.read_csv(outputs / 'real_gene_kkt_checks.tsv', sep='\t')
hires = pd.read_csv(outputs / 'core_vs_official_hires.tsv', sep='\t')
group = pd.read_csv(outputs / 'core_vs_official_group.tsv', sep='\t')
summary = {'status': 'core_benchmark_passed_not_full_hires', 'genes': len(mixture), 'samples': fractions.shape[0], 'types': fractions.shape[1], 'fraction_rank': rank, 'fraction_condition_number': float(np.linalg.cond(fractions.to_numpy())), 'nnls_fits': int(timing['nnls_fits'].sum()), 'measured_core_seconds': timing[['global_nnls_seconds', 'bootstrap_seconds', 'window_seconds']].sum().to_dict(), 'cross_solver_genes_passed': int(comparison['passed'].sum()), 'max_relative_cross_solver_error': float(comparison['relative_coefficient_error'].max()), 'max_scaled_dual_violation': float(kkt['dual_violation_scaled'].max()), 'max_scaled_complementarity': float(kkt['complementarity_scaled'].max()), 'official_complete_gene_type_pairs': int(hires['n_official_complete'].sum()), 'official_varying_gene_type_pairs': int(hires['n_official_varying'].sum()), 'median_type_bootstrap_vs_official_group_log2_mae': float(group['bootstrap_median_log2_mae'].median()), 'all_source_files_unchanged': True, 'missing_stages': ['independent fractions and S-mode', 'full gene-space preprocessing and BH', 'detectability and differential filters', 'two-segment refinement and rare-type constant overrides', 'OLS calibration and final output rules', 'official random streams and current website version parity'], 'clinical_interpretation': 'No discovery or recovery inference from this uncalibrated core', 'full_runtime_estimate': 'Not provided: most HiRes statistical steps are not implemented'}
(outputs / 'benchmark_summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps(summary, ensure_ascii=False, indent=2))

# Step4. 核对按目标基因拆出的结果仍一致，支持以后固定全样本的基因分工
partition = audit / 'core_benchmark_partition_1_2'
partition_checks = []
for path in [outputs / 'core_global_nnls.tsv', outputs / 'core_global_bootstrap_median.tsv', *sorted(outputs.glob('CORE_ONLY_*.tsv'))]:
    full = pd.read_csv(path, sep='\t', index_col=0)
    partial = pd.read_csv(partition / path.name, sep='\t', index_col=0)
    assert partial.index.tolist() == full.index[:2].tolist() and partial.columns.tolist() == full.columns.tolist()
    error = float(np.max(np.abs(partial.to_numpy() - full.iloc[:2].to_numpy())))
    partition_checks.append({'file': path.name, 'max_absolute_difference': error, 'passed': error == 0})
assert all(item['passed'] for item in partition_checks)
pd.DataFrame(partition_checks).to_csv(outputs / 'target_gene_partition_checks.tsv', sep='\t', index=False)
print('固定种子和全样本下，目标基因分工一致性检查：', len(partition_checks), '项通过')

# Step5. 记录新静态证据和源版本，不把分析器临时类型解释说成原始源码
typed = audit / 'decompiled_numeric/CIBERSORTxHiRes_randomization_typed'
evidence_paths = [typed / 'typed_001fc958.c', typed / 'typed_001fc71a.c', typed / 'typed_00129f15.c']
evidence = {'method': 'Static numerical-only ABI interpretation; not source recovery or executable parity', 'binary_sha256': hashlib.sha256((audit / 'vendor_reference/hires/src/CIBERSORTxHiRes').read_bytes()).hexdigest(), 'randomization': 'Separate full-length permutation without replacement for each column; no shared whole-row permutation', 'spearman_branch': 'prho when n<1291 and no ties; t approximation otherwise; exact vendor RNG and prho not validated', 'exports': [{'path': str(path.relative_to(audit)), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()} for path in evidence_paths], 'resource_scenario': {'genes': 17040, 'samples': 2025, 'types': 20, 'double_cube_bytes': 17040 * 2025 * 20 * 8, 'default_window': 80, 'sliding_windows': 2025 - 80 + 1, 'only_global_and_sliding_nnls_fits': 17040 * (1 + 100 + (2025 - 80 + 1) * 5)}}
(audit / 'numeric_core_and_randomization_evidence.json').write_text(json.dumps(evidence, ensure_ascii=False, indent=2), encoding='utf-8')
