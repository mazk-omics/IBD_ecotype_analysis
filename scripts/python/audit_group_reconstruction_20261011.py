# 审计独立Group估计、官方CV条件下的过滤，以及独立Group接入HiRes后的变化。
# 科学Python环境运行；只读取官方原件和研究结果，输出诊断，不能当作全流程等价验证。
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control

# Step1. 完整17040基因用于Group审计；原始名称与官方点号别名严格一对一
project = Path(__file__).resolve().parents[2]
audit = project / 'CIBERSORTx_server_plan_20261011'
official_dir = project / 'CIBERSORTx_signature_audit/MSCCR_Joint20v2_Job24_HiRes_review_20261009/original_outputs'
group_dir = audit / 'group_independent_run1'
outdir = audit / 'group_independent_validation'
outdir.mkdir(exist_ok=True)
independent = pd.read_csv(group_dir / 'GEPs.tsv', sep='\t', index_col=0)
alias = independent.index.str.replace('-', '.', regex=False)
assert alias.is_unique
comparison = []
sources = []
official_tables = {}
for name in ['GEPs', 'GEPs_Filtered', 'GEPs_CVs', 'GEPs_Pvals', 'GEPs_Qvals', 'GEPs_StdErrs']:
    path = official_dir / f'CIBERSORTxGEP_Job24_{name}.txt'
    official = pd.read_csv(path, sep='\t', index_col=0)
    assert official.index.is_unique and set(official.index) == set(alias)
    official = official.loc[alias, independent.columns]
    official.index = independent.index
    official_tables[name] = official
    prediction = pd.read_csv(group_dir / f'{name}.tsv', sep='\t', index_col=0)
    assert prediction.index.equals(independent.index) and prediction.columns.equals(independent.columns)
    sources.append({'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()})
    for celltype in independent.columns:
        x = prediction[celltype].to_numpy()
        y = official[celltype].to_numpy()
        common = np.isfinite(x) & np.isfinite(y)
        error = abs(np.log2(x[common] + 1) - np.log2(y[common] + 1))
        log_x, log_y = np.log2(x[common] + 1), np.log2(y[common] + 1)
        correlation = float(np.corrcoef(log_x, log_y)[0, 1]) if np.std(log_x) > 0 and np.std(log_y) > 0 else np.nan
        comparison.append({'table': name, 'CellType': celltype, 'common_pairs': int(common.sum()), 'NA_mismatch': int(sum(np.isnan(x) != np.isnan(y))), 'median_log2_absolute_error': float(np.median(error)), 'mean_log2_absolute_error': float(np.mean(error)), 'log2_correlation': correlation})
pd.DataFrame(comparison).to_csv(outdir / 'group_table_comparisons.tsv', sep='\t', index=False)

# Step2. 只给官方CV，独立重算阈值；不读取官方NA来选择阈值
cv = official_tables['GEPs_CVs']
filtered = official_tables['GEPs_Filtered']
filter_records = []
for celltype in cv:
    values = cv[celltype].to_numpy()
    sorted_cv = np.sort(values)
    ycut = max(1, np.quantile(values[values < 999999], 0.75))
    endpoint = np.argmin(abs(sorted_cv - ycut)) + 1
    rank = np.arange(1, endpoint + 1) / endpoint
    threshold = sorted_cv[:endpoint][np.argmin((rank - ycut)**2 + sorted_cv[:endpoint]**2)]
    predicted_na = values >= threshold
    observed_na = filtered[celltype].isna().to_numpy()
    filter_records.append({'CellType': celltype, 'ycut': float(ycut), 'endpoint_rank': int(endpoint), 'threshold': float(threshold), 'official_NA': int(observed_na.sum()), 'conditional_filter_NA_mismatch': int(sum(predicted_na != observed_na))})
pd.DataFrame(filter_records).to_csv(outdir / 'official_CV_conditional_filter_check.tsv', sep='\t', index=False)
assert sum(item['conditional_filter_NA_mismatch'] for item in filter_records) == 0

# Step3. 官方P到Q的条件核对，缺失P单独记录；不将这些Q替换研究估计
pvalues = official_tables['GEPs_Pvals']
qvalues = official_tables['GEPs_Qvals']
bh_records = []
for celltype in pvalues:
    p = pvalues[celltype].to_numpy()
    tested = np.isfinite(p)
    q = np.ones(len(p))
    q[tested] = false_discovery_control(p[tested], method='bh')
    observed = qvalues[celltype].to_numpy()
    bh_records.append({'CellType': celltype, 'tested_genes': int(tested.sum()), 'max_absolute_difference': float(np.max(abs(q - observed))), 'differences_over_1e_10': int(sum(abs(q - observed) > 1e-10))})
pd.DataFrame(bh_records).to_csv(outdir / 'official_P_conditional_BH_check.tsv', sep='\t', index=False)

# 系数本身与强marker置零规则分开看；短向量版本仅作为来源诊断，绝不进入主流程
beta = pd.read_csv(group_dir / 'bootstrap_beta_before_q_guard.tsv', sep='\t', index_col=0)
independent_p = pd.read_csv(group_dir / 'GEPs_Pvals.tsv', sep='\t', index_col=0)
independent_q = pd.read_csv(group_dir / 'GEPs_Qvals.tsv', sep='\t', index_col=0)
legacy_p = pd.read_csv(audit / 'group_legacy_residual_diagnostic/recycled_residual_Pvals.tsv', sep='\t', index_col=0)
legacy_q = pd.read_csv(audit / 'group_legacy_residual_diagnostic/recycled_residual_Qvals.tsv', sep='\t', index_col=0)
guard_records = []
for celltype in independent.columns:
    before, observed = beta[celltype].to_numpy(), official_tables['GEPs'][celltype].to_numpy()
    tested = np.isfinite(pvalues[celltype])
    guard_records.append({'CellType': celltype, 'before_guard_mean_log2_error': float(np.mean(abs(np.log2(before + 1) - np.log2(observed + 1)))), 'before_guard_log2_correlation': float(np.corrcoef(np.log2(before + 1), np.log2(observed + 1))[0, 1]), 'full_residual_P_mean_absolute_error': float(np.mean(abs(independent_p.loc[tested, celltype] - pvalues.loc[tested, celltype]))), 'recycled_residual_P_mean_absolute_error': float(np.mean(abs(legacy_p.loc[tested, celltype] - pvalues.loc[tested, celltype])))})
pd.DataFrame(guard_records).to_csv(outdir / 'coefficient_and_residual_version_diagnostic.tsv', sep='\t', index=False)
strong_counts = {'official_Group': int((qvalues.min(axis=1) < 1e-5).sum()), 'independent_full_residual': int((independent_q.min(axis=1) < 1e-5).sum()), 'diagnostic_recycled_residual': int((legacy_q.min(axis=1) < 1e-5).sum())}

# Step4. 新Group作为唯一辅助表，比较完整1249基因；上游fraction仍是官方冻结输入
mapping = pd.read_csv(audit / 'full_target_benchmark_inputs/gene_name_mapping.tsv', sep='\t')
pair_tables = []
for celltype in independent.columns:
    path = official_dir / f'CIBERSORTxHiRes_Job24_{celltype}_Window60.txt'
    official = pd.read_csv(path, sep='\t', index_col=0)
    assert official.index.tolist() == mapping['official_output_gene'].tolist()
    official.index = mapping['GeneSymbol'].tolist()
    prediction = pd.read_csv(audit / 'auxiliary_independent_group_run1' / f'AUX_DEPENDENT_{celltype}.tsv', sep='\t', index_col=0)
    assert official.index.equals(prediction.index) and official.columns.equals(prediction.columns)
    x = prediction.to_numpy()
    y = official.to_numpy()
    official_na = np.isnan(y).all(axis=1)
    predicted_na = np.isnan(x).all(axis=1)
    official_variable = ~official_na & (np.ptp(np.nan_to_num(y), axis=1) > 1e-8)
    predicted_variable = ~predicted_na & (np.ptp(np.nan_to_num(x), axis=1) > 1e-8)
    covered = np.isfinite(x).all(axis=1) & np.isfinite(y).all(axis=1)
    error = np.full(len(x), np.nan)
    error[covered] = np.mean(abs(np.log2(x[covered] + 1) - np.log2(y[covered] + 1)), axis=1)
    pair_tables.append(pd.DataFrame({'GeneSymbol': official.index, 'CellType': celltype, 'official_NA': official_na, 'predicted_NA': predicted_na, 'official_variable': official_variable, 'predicted_variable': predicted_variable, 'all_samples_covered': covered, 'log2_MAE': error}))
pairs = pd.concat(pair_tables, ignore_index=True)
pairs.to_csv(outdir / 'hires_independent_group_pairs.tsv', sep='\t', index=False)
pairs.loc[pairs['official_NA'] != pairs['predicted_NA']].to_csv(outdir / 'hires_independent_group_NA_disagreements.tsv', sep='\t', index=False)
variable = pairs.loc[pairs['official_variable']]
summary = {'genes_Group': len(independent), 'types': len(independent.columns), 'conditional_filter_pairs': int(cv.size), 'conditional_filter_mismatch': sum(item['conditional_filter_NA_mismatch'] for item in filter_records), 'conditional_BH_max_error': max(item['max_absolute_difference'] for item in bh_records), 'strong_marker_gene_counts': strong_counts, 'before_guard_median_type_log2_correlation': float(np.median([r['before_guard_log2_correlation'] for r in guard_records])), 'Group_independent_filter_NA_mismatch': int(sum(item['NA_mismatch'] for item in comparison if item['table'] == 'GEPs_Filtered')), 'hires_pairs': len(pairs), 'official_NA': int(pairs['official_NA'].sum()), 'predicted_NA': int(pairs['predicted_NA'].sum()), 'hires_NA_missed': int((pairs['official_NA'] & ~pairs['predicted_NA']).sum()), 'hires_false_NA': int((~pairs['official_NA'] & pairs['predicted_NA']).sum()), 'official_variable_pairs': len(variable), 'covered_variable_pairs': int(variable['all_samples_covered'].sum()), 'variable_predictions_on_official_variable': int(variable['predicted_variable'].sum()), 'median_log2_MAE_on_covered_variable_pairs': float(variable['log2_MAE'].median()), 'sources': sources, 'limits': ['CV conditional formula agreement is not independent CV estimation agreement.', 'Independent Group uses full residual and paper zero guard; legacy cached code has dimension and guard differences.', 'Fixed official fractions remain upstream; no independent signature, fractions or S-mode.', 'HiRes worker and refinement reused unchanged; differences in NA do not demonstrate better biological accuracy.']}
(outdir / 'group_validation_summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')
print(json.dumps({k: v for k, v in summary.items() if k != 'sources'}, ensure_ascii=False, indent=2, allow_nan=False))
