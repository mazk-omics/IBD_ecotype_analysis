# 按BayesPrism实际保留基因计算外部RNA真值，保存细胞数尺度和完整20类准确性结果。
# 服务器ibd_refaudit运行；只读冻结抽样manifest和counts，不运行研究性SVR或EcoTyper。
from pathlib import Path
import json
import sys
import numpy as np
import pandas as pd
from scipy import sparse, stats

# Step1. 通过细胞与基因清单精确恢复实际分母
parent = Path(sys.argv[1])
dataset = sys.argv[2]
folder = parent / '08_external_pseudobulk' / dataset
input = folder / '02_simulation'
output = folder / '04_accuracy'
types = (parent / '00_audit/frozen_types.txt').read_text().splitlines()
settings = json.loads((folder / '00_audit/endpoint_prespecification.json').read_text())
coverage = pd.read_csv(folder / '00_audit/Full20_type_coverage.tsv', sep='\t').set_index('cell_type')
evaluable = coverage.index[coverage.evaluable].tolist()
members = pd.read_csv(input / 'cell_sampling_manifest.tsv.gz', sep='\t')
pool = pd.read_csv(input / 'unique_source_cell_manifest.tsv.gz', sep='\t')
source_genes = pd.Index((input / 'source_genes.tsv').read_text().splitlines())
counts = sparse.load_npz(input / 'unique_source_cell_counts_csr.npz')
retained_genes = pd.Index((folder / '03_bayesprism/retained_genes.txt').read_text().splitlines())
positions = source_genes.get_indexer(retained_genes)
assert (positions >= 0).all()
cell_rows = pd.Index(pool.CellID).get_indexer(members.CellID)
assert (cell_rows >= 0).all() and len(members) == 100000
cell_reads = counts[:, positions].sum(axis=1).A1
members['actual_retained_gene_reads'] = cell_reads[cell_rows]
RNA_truth = members.groupby(['sample', 'CellType'], observed=True).actual_retained_gene_reads.sum().unstack(fill_value=0).reindex(columns=types, fill_value=0)
RNA_truth = RNA_truth.div(RNA_truth.sum(axis=1), axis=0)
RNA_truth.to_csv(output / 'truth_BP_actual_retained_gene_RNA_fraction.tsv', sep='\t', index_label='sample')
cell_truth = pd.read_csv(input / 'truth_cell_count_fraction.tsv', sep='\t', index_col=0).loc[RNA_truth.index, types]
full_truth = pd.read_csv(input / 'truth_full_gene_RNA_fraction.tsv', sep='\t', index_col=0).loc[RNA_truth.index, types]
initial = pd.read_csv(folder / '03_bayesprism/theta_initial.tsv', sep='\t', index_col=0).loc[RNA_truth.index, types]
final = pd.read_csv(folder / '03_bayesprism/theta_final.tsv', sep='\t', index_col=0).loc[RNA_truth.index, types]
assert final.shape == (100, 20) and np.isfinite(final).all().all() and final.ge(0).all().all() and np.allclose(final.sum(axis=1), 1)
assert np.allclose(RNA_truth.sum(axis=1), 1)
denominators = members.groupby('sample', observed=True).actual_retained_gene_reads.sum().to_frame('retained_gene_counts')
denominators['full_gene_counts'] = pd.read_csv(input / 'pseudo_bulk_raw_counts.tsv.gz', sep='\t', index_col=0).sum().reindex(denominators.index)
denominators['retained_RNA_fraction_of_full_library'] = denominators.retained_gene_counts / denominators.full_gene_counts
denominators.to_csv(output / 'actual_gene_denominator_audit.tsv', sep='\t', index_label='sample')

# Step2. 真值尺度分开评价；可评价类型不重归一化，完整20轴保留缺失类型质量
records = []
long_records = []
for truth_scale, truth in [('Retained_gene_RNA_PRIMARY', RNA_truth), ('Cell_count_SECONDARY', cell_truth), ('Full_gene_RNA_SECONDARY', full_truth)]:
    for method, prediction in [('BayesPrism_final', final), ('BayesPrism_initial', initial)]:
        for cell_type in types + ['POOLED_EVALUABLE', 'POOLED_FULL20']:
            if cell_type.startswith('POOLED'):
                axis = evaluable if cell_type == 'POOLED_EVALUABLE' else types
                x, y = prediction[axis].to_numpy().ravel(), truth[axis].to_numpy().ravel()
            else:
                x, y = prediction[cell_type].to_numpy(), truth[cell_type].to_numpy()
            delta = (x - y) * 100
            variable = x.std() > 0 and y.std() > 0
            records.append({'dataset': dataset, 'truth_scale': truth_scale, 'method': method, 'cell_type': cell_type, 'evaluable': cell_type in evaluable or cell_type == 'POOLED_EVALUABLE', 'n_mixtures': 100, 'Pearson_r': stats.pearsonr(x, y).statistic if variable else np.nan, 'Spearman_rho': stats.spearmanr(x, y).statistic if variable else np.nan, 'MAE_pp': np.abs(delta).mean(), 'RMSE_pp': np.sqrt(np.mean(delta**2)), 'mean_signed_error_pp': delta.mean()})
        predicted_long = prediction.stack().rename('predicted_fraction')
        truth_long = truth.stack().rename('truth_fraction')
        predicted_long.index.names = ['sample', 'cell_type']
        truth_long.index.names = ['sample', 'cell_type']
        assert predicted_long.index.is_unique and truth_long.index.is_unique
        long = predicted_long.to_frame().join(truth_long)
        long = long.reset_index()
        long['method'] = method
        long['truth_scale'] = truth_scale
        long_records.append(long)
metrics = pd.DataFrame(records)
metrics.to_csv(output / 'accuracy_metrics_all_scales.tsv', sep='\t', index=False)
pd.concat(long_records, ignore_index=True).to_csv(output / 'predictions_and_truth_long.tsv', sep='\t', index=False)
unsupported = [t for t in types if t not in evaluable]
unsupported_mass = final[unsupported].sum(axis=1).rename('unsupported_type_fraction')
unsupported_mass.to_csv(output / 'unsupported_type_predicted_mass.tsv', sep='\t', index_label='sample')
(final - RNA_truth).mul(100).to_csv(output / 'primary_signed_errors_pp.tsv', sep='\t', index_label='sample')
(final - cell_truth).mul(100).to_csv(output / 'cell_count_signed_errors_pp.tsv', sep='\t', index_label='sample')
summary = {'dataset': dataset, 'status': 'COMPLETE', 'mixtures': 100, 'draws_per_mixture': 1000, 'reference_types': 20, 'evaluable_types': len(evaluable), 'evaluable_type_names': evaluable, 'unsupported_types': unsupported, 'renormalized': False, 'primary_endpoint': settings['primary_accuracy_endpoint'], 'primary_result': metrics.loc[metrics.truth_scale.eq('Retained_gene_RNA_PRIMARY') & metrics.method.eq('BayesPrism_final') & metrics.cell_type.eq('POOLED_EVALUABLE')].iloc[0].to_dict(), 'cell_count_result': metrics.loc[metrics.truth_scale.eq('Cell_count_SECONDARY') & metrics.method.eq('BayesPrism_final') & metrics.cell_type.eq('POOLED_EVALUABLE')].iloc[0].to_dict(), 'retained_genes': len(retained_genes), 'unsupported_fraction_mean': float(unsupported_mass.mean()), 'unsupported_fraction_max': float(unsupported_mass.max()), 'official_CIBERSORTx_paired': 'BLOCKED', 'official_reason': 'No available official runtime and token; inputs prepared for same mixtures; research SVR not used', 'simulation_samples_not_independent_biological_donors': True, 'cross_study_donor_linkage': 'UNKNOWN; accession and source context disjoint from GSE282122'}
(output / 'accuracy_summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
# 精确RNA真值和误差已落盘；删除本轮可再生成的大稀疏细胞缓存，保留原源文件、抽样manifest和pseudo-bulk counts。
cache_path = input / 'unique_source_cell_counts_csr.npz'
cache_bytes = cache_path.stat().st_size
cache_path.unlink()
(output / 'compact_cache_audit.json').write_text(json.dumps({'removed_generated_cell_counts_cache': str(cache_path), 'bytes_saved': cache_bytes, 'original_raw_sources_modified': False, 'sampling_manifest_retained': True, 'pseudo_bulk_raw_counts_retained': True, 'RNA_truth_retained': True, 'regeneration_basis': 'unchanged source counts + exact sampling manifest'}, indent=2), encoding='utf-8')
(folder / 'execution_status.json').write_text(json.dumps({'audit': 'COMPLETE', 'smoke': 'COMPLETE', 'generation': 'COMPLETE', 'BayesPrism': 'COMPLETE', 'accuracy': 'COMPLETE', 'official_CIBERSORTx_paired': 'BLOCKED', 'biological_cross_study_donor_independence': 'PARTIAL'}, indent=2), encoding='utf-8')
print(json.dumps(summary, indent=2), flush=True)
