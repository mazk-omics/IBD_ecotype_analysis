# 审计1249基因的HiRes研究输出：全上下文回填、名称映射、旧24基因保持性与数值一致性。
# 本机科学Python环境运行；读取新研究输出和官方原件，写全部组合及限制，不按表现筛选。
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

# Step1. 固定20类型和全部目标基因；辅助回填的外部依赖必须与独立数值部分分开报告
project = Path(__file__).resolve().parents[2]
audit = project / 'CIBERSORTx_server_plan_20261011'
inputs = audit / 'full_target_benchmark_inputs'
outdir = audit / 'full_target_validation'
outdir.mkdir(exist_ok=True)  # 只重写本脚本诊断表，不改研究表达与原件
mapping = pd.read_csv(inputs / 'gene_name_mapping.tsv', sep='\t')
fractions = pd.read_csv(inputs / 'fractions_fixed.tsv', sep='\t', index_col=0)
official_dir = project / 'CIBERSORTx_signature_audit/MSCCR_Joint20v2_Job24_HiRes_review_20261009/original_outputs'
frozen = pd.read_csv(audit / 'core_benchmark_inputs/gene_selection.tsv', sep='\t')['GeneSymbol'].tolist()
pair_tables = []
repeat_records = []
source_records = []

# Step2. 每类型用完整1249基因核对；统计公式向量化，NA缺失与恒定行独立标记
for celltype in fractions.columns:
    path = official_dir / f'CIBERSORTxHiRes_Job24_{celltype}_Window60.txt'
    official = pd.read_csv(path, sep='\t', index_col=0)
    refined = pd.read_csv(audit / 'refinement_full_target_run1' / f'REFINEMENT_ONLY_{celltype}.tsv', sep='\t', index_col=0)
    refilled = pd.read_csv(audit / 'auxiliary_full_target_run2' / f'AUX_DEPENDENT_{celltype}.tsv', sep='\t', index_col=0)
    assert official.index.tolist() == mapping['official_output_gene'].tolist()
    official.index = mapping['GeneSymbol'].tolist()
    assert refined.index.tolist() == refilled.index.tolist() == official.index.tolist()
    assert refined.columns.tolist() == refilled.columns.tolist() == official.columns.tolist() == fractions.index.tolist()
    source_records.append({'path': str(path.resolve()), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()})
    y = official.to_numpy()
    valid_y = np.isfinite(y)
    n_y = valid_y.sum(axis=1)
    y_range = np.max(np.where(valid_y, y, -np.inf), axis=1) - np.min(np.where(valid_y, y, np.inf), axis=1)
    status = np.where(n_y == 0, 'NA_all_samples', np.where(n_y < y.shape[1], 'partial_NA', np.where(y_range > 1e-8, 'variable', 'constant')))
    for stage, prediction in [('refinement', refined), ('auxiliary_dependent_refill', refilled)]:
        x = prediction.to_numpy()
        valid_x = np.isfinite(x)
        n_x = valid_x.sum(axis=1)
        common = valid_x & valid_y
        n_common = common.sum(axis=1)
        x_range = np.max(np.where(valid_x, x, -np.inf), axis=1) - np.min(np.where(valid_x, x, np.inf), axis=1)
        error = np.sum(np.where(common, np.abs(np.log2(np.where(valid_x, x, 0) + 1) - np.log2(np.where(valid_y, y, 0) + 1)), 0), axis=1)
        mae = np.divide(error, n_common, out=np.full(len(y), np.nan), where=n_common > 0)
        log_x = np.log2(np.where(common, x, 0) + 1)
        log_y = np.log2(np.where(common, y, 0) + 1)
        mean_x = np.divide(np.sum(log_x, axis=1), n_common, out=np.zeros(len(y)), where=n_common > 0)
        mean_y = np.divide(np.sum(log_y, axis=1), n_common, out=np.zeros(len(y)), where=n_common > 0)
        centered_x = np.where(common, log_x - mean_x[:, None], 0)
        centered_y = np.where(common, log_y - mean_y[:, None], 0)
        denominator = np.sqrt(np.sum(centered_x**2, axis=1) * np.sum(centered_y**2, axis=1))
        correlation = np.divide(np.sum(centered_x * centered_y, axis=1), denominator, out=np.full(len(y), np.nan), where=(denominator > 0) & (status == 'variable') & (x_range > 1e-8))
        pair_tables.append(pd.DataFrame({'GeneSymbol': official.index, 'CellType': celltype, 'stage': stage, 'official_status': status, 'predicted_all_NA': n_x == 0, 'predicted_variable': (n_x == y.shape[1]) & (x_range > 1e-8), 'common_samples': n_common, 'log2_mae': mae, 'log2_correlation': correlation}))
    previous = pd.read_csv(audit / 'refinement_canonical24_run1' / f'REFINEMENT_ONLY_{celltype}.tsv', sep='\t', index_col=0)
    difference = refined.loc[frozen].to_numpy() - previous.loc[frozen].to_numpy()
    historical = pd.read_csv(audit / 'refinement_original_run2' / f'REFINEMENT_ONLY_{celltype}.tsv', sep='\t', index_col=0)
    legacy_difference = refined.loc[frozen].to_numpy() - historical.loc[frozen].to_numpy()
    repeat_records.append({'CellType': celltype, 'max_absolute_difference_from_canonical_24': float(np.max(np.abs(difference))), 'exactly_equal': bool((difference == 0).all()), 'max_absolute_difference_from_legacy_seed_24': float(np.max(np.abs(legacy_difference)))})
assert all(item['exactly_equal'] for item in repeat_records)
pd.DataFrame(repeat_records).to_csv(outdir / 'frozen_24_vs_expanded_context.tsv', sep='\t', index=False)
pairs = pd.concat(pair_tables, ignore_index=True)
pairs.to_csv(outdir / 'all_gene_type_pairs_both_stages.tsv', sep='\t', index=False)

# Step3. 报告覆盖和漏检，不把预测NA或恒定的失败组合从分母中抹掉
summaries = []
for stage, frame in pairs.groupby('stage', sort=False):
    variable = frame.loc[frame['official_status'] == 'variable']
    covered = variable.loc[variable['common_samples'] == 120]
    official_na = frame['official_status'] == 'NA_all_samples'
    false_na = (~official_na) & frame['predicted_all_NA']
    na_missed = official_na & ~frame['predicted_all_NA']
    summaries.append({'stage': stage, 'pairs': len(frame), 'official_variable_pairs': len(variable), 'variable_pairs_all_samples_covered': len(covered), 'variable_predictions_on_official_variable': int(variable['predicted_variable'].sum()), 'correlation_at_least_0_95_on_official_variable': int((variable['log2_correlation'] >= 0.95).sum()), 'median_variable_log2_mae_on_covered_pairs': float(covered['log2_mae'].median()), 'official_NA_pairs': int(official_na.sum()), 'predicted_NA_pairs': int(frame['predicted_all_NA'].sum()), 'NA_missed': int(na_missed.sum()), 'false_NA': int(false_na.sum()), 'official_constant_predicted_variable': int(((frame['official_status'] == 'constant') & frame['predicted_variable']).sum())})
    frame.loc[false_na | na_missed].to_csv(outdir / f'{stage}_NA_disagreements.tsv', sep='\t', index=False)
    variable.sort_values('log2_mae', ascending=False).to_csv(outdir / f'{stage}_variable_pairs_by_error.tsv', sep='\t', index=False)
    pd.crosstab(frame['official_status'], [frame['predicted_all_NA'], frame['predicted_variable']]).to_csv(outdir / f'{stage}_status_confusion.tsv', sep='\t')

# 漏掉变化的组合按实际分支拆开，辅助NA不同也附上前面统计，便于逐条修复
worker = pd.read_csv(audit / 'worker_full_target_run2/worker_type_statistics.tsv', sep='\t')
refinement_statistics = pd.read_csv(audit / 'refinement_full_target_run1/refinement_type_statistics.tsv', sep='\t')
trace = pairs.loc[pairs['stage'] == 'refinement'].merge(worker, on=['GeneSymbol', 'CellType'], validate='one_to_one')
trace = trace.merge(refinement_statistics, on=['GeneSymbol', 'CellType'], validate='one_to_one')
missed = trace.loc[(trace['official_status'] == 'variable') & ~trace['predicted_variable']].copy()
missed['local_branch'] = np.where(missed['global_p'] >= 0.05, 'excluded_by_global_detection', np.where(missed['weak_type_refitted'], 'weak_variation_constant_refit', 'other'))
assert not (missed['local_branch'] == 'other').any()
missed.to_csv(outdir / 'missed_variation_branch_trace.tsv', sep='\t', index=False)
mask = pd.read_csv(outdir / 'auxiliary_dependent_refill_NA_disagreements.tsv', sep='\t')
mask.merge(worker, on=['GeneSymbol', 'CellType'], validate='one_to_one').merge(refinement_statistics, on=['GeneSymbol', 'CellType'], validate='one_to_one').to_csv(outdir / 'NA_disagreement_branch_trace.tsv', sep='\t', index=False)

# Step4. 输入来源哈希复核；原始研究文件不变，官方aux依赖单独登记
provenance = json.loads((inputs / 'input_provenance.json').read_text(encoding='utf-8'))
for item in provenance['sources']:
    assert hashlib.sha256(Path(item['path']).read_bytes()).hexdigest() == item['sha256']
aux_path = official_dir / 'CIBERSORTxGEP_Job24_GEPs_Filtered.txt'
source_records.append({'path': str(aux_path.resolve()), 'sha256': hashlib.sha256(aux_path.read_bytes()).hexdigest()})
worker_summary = pd.read_csv(audit / 'worker_full_target_run2/worker_statistics_summary.tsv', sep='\t').to_dict('records')
refinement_summary = pd.read_csv(audit / 'refinement_full_target_run1/refinement_summary.tsv', sep='\t').to_dict('records')
summary = {'status': 'Full target numerical benchmark; auxiliary refill still depends on official Group input', 'genes': 1249, 'samples': 120, 'types': 20, 'gene_name_aliases': provenance['gene_name_aliases'], 'frozen_24_exactly_preserved_types': len(repeat_records), 'canonical_seed_correction': 'Historical 24-gene worker seed indices were inherited from Adjusted export; canonical indices now use the original mixture. Same expression and thresholds.', 'missed_variation_branches': missed['local_branch'].value_counts().to_dict(), 'stages': summaries, 'worker_run': worker_summary, 'refinement_run': refinement_summary, 'sources': source_records, 'limits': ['Not independent full workflow: fixed official fractions and supplied official Group reliability output.', 'No vendor binary executed and no official authentication logic modified.', 'Raw numerical stage still has detection/shrinkage/randomness/version discrepancies.', 'Do not use these benchmark outputs as formal EcoTyper discovery/recovery.']}
(outdir / 'full_target_validation_summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')
print(json.dumps({k: v for k, v in summary.items() if k != 'sources'}, ensure_ascii=False, indent=2, allow_nan=False))
