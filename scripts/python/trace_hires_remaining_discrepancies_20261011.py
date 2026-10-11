# 定位冻结面板剩余的不一致：检测、弱变化回填、随机性与官方NA之间的关系。
# 本机科学Python环境运行；只读取研究检查点摘要与官方输出，不改表达、不套官方掩码。
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

# Step1. 合并全480个组合的状态与两个数值阶段，避免只挑相关性高的行
project = Path(__file__).resolve().parents[2]
audit = project / 'CIBERSORTx_server_plan_20261011'
outdir = audit / 'refinement_validation'
official_dir = project / 'CIBERSORTx_signature_audit/MSCCR_Joint20v2_Job24_HiRes_review_20261009/original_outputs'
pairs = pd.read_csv(outdir / 'all_480_pairs_vs_official.tsv', sep='\t')
worker = pd.read_csv(audit / 'worker_statistics_original_run1/worker_type_statistics.tsv', sep='\t')
refinement = pd.read_csv(audit / 'refinement_original_run2/refinement_type_statistics.tsv', sep='\t')
trace = pairs.merge(worker, on=['GeneSymbol', 'CellType'], validate='one_to_one')
trace = trace.merge(refinement, on=['GeneSymbol', 'CellType'], validate='one_to_one')
trace['active_global_p_below_0_05'] = trace['global_p'] < 0.05
trace['constant_after_refinement'] = np.abs(trace['segment1_final'] - trace['segment2_final']) <= 1e-5
assert len(trace) == 480
trace.to_csv(outdir / 'all_pairs_stage_trace.tsv', sep='\t', index=False)

# Step2. 官方Group过滤不能直接当HiRes掩码；全表保存反例，供后续外层规则定位
source_paths = [official_dir / f'CIBERSORTxGEP_Job24_{suffix}.txt' for suffix in ['GEPs_Filtered', 'GEPs_Qvals', 'GEPs_Pvals']]
filtered = pd.read_csv(source_paths[0], sep='\t', index_col=0)
group_q = pd.read_csv(source_paths[1], sep='\t', index_col=0)
group_p = pd.read_csv(source_paths[2], sep='\t', index_col=0)
trace['official_group_filtered_NA'] = [pd.isna(filtered.loc[g, c]) for g, c in zip(trace['GeneSymbol'], trace['CellType'])]
trace['official_group_q'] = [group_q.loc[g, c] for g, c in zip(trace['GeneSymbol'], trace['CellType'])]
trace['official_group_p'] = [group_p.loc[g, c] for g, c in zip(trace['GeneSymbol'], trace['CellType'])]
pd.crosstab(trace['official_status'], trace['official_group_filtered_NA']).to_csv(outdir / 'official_hires_status_vs_group_filter.tsv', sep='\t')
trace.loc[trace['official_group_filtered_NA'] & (trace['official_status'] != 'NA_all_samples')].to_csv(outdir / 'group_NA_but_hires_finite.tsv', sep='\t', index=False)
assert int(((trace['official_status'] == 'NA_all_samples') & ~trace['official_group_filtered_NA']).sum()) == 0

# Step3. 优先列出漏检和错误变化，再看低相关性；不按官方表现改变阈值或种子
missed = trace.loc[(trace['official_status'] == 'variable') & ~trace['predicted_variable']].copy()
extra = trace.loc[(trace['official_status'] == 'constant') & trace['predicted_variable']].copy()
missed['local_branch'] = np.where(~missed['active_global_p_below_0_05'], 'excluded_by_global_detection', np.where(missed['weak_type_refitted'], 'weak_variation_residual_refit_constant', 'other_constant_branch'))
missed.to_csv(outdir / 'official_variation_missed_branch_trace.tsv', sep='\t', index=False)
extra.to_csv(outdir / 'official_constant_but_predicted_variable.tsv', sep='\t', index=False)
trace.loc[trace['official_status'] == 'variable'].sort_values('refined_log2_correlation', na_position='first').to_csv(outdir / 'variable_pairs_sorted_by_correlation.tsv', sep='\t', index=False)
seeds = pd.read_csv(outdir / 'refinement_seed_sensitivity.tsv', sep='\t')
c1qb = seeds.loc[(seeds['GeneSymbol'] == 'C1QB') & (seeds['CellType'] == 'Plasma')]
baseline_c1qb = trace.loc[(trace['GeneSymbol'] == 'C1QB') & (trace['CellType'] == 'Plasma')].iloc[0]

# Step4. 同种子两次结果应一致；来源哈希与阶段种子限制一并保存
group = pd.read_csv(audit / 'core_benchmark_inputs/official_group_24genes.tsv', sep='\t', index_col=0)
repeat_checks = []
for celltype in group.columns:
    x = pd.read_csv(audit / 'refinement_original_run1' / f'REFINEMENT_ONLY_{celltype}.tsv', sep='\t', index_col=0)
    y = pd.read_csv(audit / 'refinement_original_run2' / f'REFINEMENT_ONLY_{celltype}.tsv', sep='\t', index_col=0)
    assert x.index.tolist() == y.index.tolist() and x.columns.tolist() == y.columns.tolist()
    repeat_checks.append({'CellType': celltype, 'exactly_equal': bool(np.array_equal(x.to_numpy(), y.to_numpy())), 'max_absolute_difference': float(np.max(np.abs(x.to_numpy() - y.to_numpy())))})
assert all(item['exactly_equal'] for item in repeat_checks)
pd.DataFrame(repeat_checks).to_csv(outdir / 'same_seed_repeatability.tsv', sep='\t', index=False)

# Step5. 旧外层回填有整类型全1保护；24基因面板不能代表完整1249基因的这个判断
# 仅检验该规则是否解释NA，完整官方结果只作诊断上下文，不用于修改研究预测。
refill_records = []
context_sources = []
for celltype in group.columns:
    prediction = pd.read_csv(audit / 'refinement_original_run2' / f'REFINEMENT_ONLY_{celltype}.tsv', sep='\t', index_col=0)
    full_path = official_dir / f'CIBERSORTxHiRes_Job24_{celltype}_Window60.txt'
    official_full = pd.read_csv(full_path, sep='\t', index_col=0)
    context_sources.append({'path': str(full_path.resolve()), 'sha256': hashlib.sha256(full_path.read_bytes()).hexdigest()})
    panel_all_one = bool((prediction.to_numpy() == 1).all())
    official_full_all_one = bool((official_full.to_numpy() == 1).all())
    for gene in prediction.index:
        candidate_refill = bool(prediction.loc[gene].max() == 1)
        auxiliary_na = bool(pd.isna(filtered.loc[gene, celltype]))
        observed_na = bool(trace.loc[(trace['GeneSymbol'] == gene) & (trace['CellType'] == celltype), 'official_status'].iloc[0] == 'NA_all_samples')
        refill_records.append({'GeneSymbol': gene, 'CellType': celltype, 'panel_type_all_one': panel_all_one, 'official_full_type_all_one': official_full_all_one, 'candidate_refill_row_max_one': candidate_refill, 'group_filtered_NA': auxiliary_na, 'predicted_NA_panel_context': candidate_refill and auxiliary_na and not panel_all_one, 'diagnostic_NA_official_full_context': candidate_refill and auxiliary_na and not official_full_all_one, 'observed_official_NA': observed_na})
refill = pd.DataFrame(refill_records)
refill.to_csv(outdir / 'outer_refill_context_diagnostic.tsv', sep='\t', index=False)
context_summary = {'panel_guard_NA_count': int(refill['predicted_NA_panel_context'].sum()), 'official_full_context_guard_NA_count': int(refill['diagnostic_NA_official_full_context'].sum()), 'official_full_context_guard_mismatches': refill.loc[refill['diagnostic_NA_official_full_context'] != refill['observed_official_NA'], ['GeneSymbol', 'CellType']].to_dict('records'), 'source_files': context_sources, 'limits': 'Official full-context all-one guard is an oracle diagnostic, not an independent prediction. No output mask or expression changed.'}
(outdir / 'outer_refill_context_summary.json').write_text(json.dumps(context_summary, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')
source_records = [{'path': str(path.resolve()), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()} for path in source_paths]
summary = {'panel_gene_type_pairs': len(trace), 'missed_variable_pairs': missed[['GeneSymbol', 'CellType', 'local_branch']].to_dict('records'), 'extra_variable_pairs': extra[['GeneSymbol', 'CellType']].to_dict('records'), 'official_all_NA_pairs': int((trace['official_status'] == 'NA_all_samples').sum()), 'group_filtered_NA_pairs': int(trace['official_group_filtered_NA'].sum()), 'group_filtered_NA_but_hires_finite': int((trace['official_group_filtered_NA'] & (trace['official_status'] != 'NA_all_samples')).sum()), 'C1QB_Plasma_baseline_log2_correlation': float(baseline_c1qb['refined_log2_correlation']), 'C1QB_Plasma_alternate_stage_seed_correlations': c1qb['log2_correlation_vs_official'].tolist(), 'same_seed_exact_repeat_types': len(repeat_checks), 'source_files': source_records, 'interpretation': 'All HiRes NA pairs in this panel are Group-filtered NA, but 20 other Group-filtered NA pairs are finite in HiRes; direct mask reuse fails.', 'limits': 'Only refinement seeds varied; initial worker and full-gene entry statistics remain fixed. No independent full workflow or current-webpage equivalence.'}
(outdir / 'remaining_discrepancies_summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')
print(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False))
