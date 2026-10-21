# 只读复算Joint20外部验证和真实bulk一致性，补充最小的组成误差与分层效应审计。
# 本机科学Python运行；结果只写本次closeout目录，不更改旧truth、预测、signature或模型。
from pathlib import Path
import hashlib
import json
import shutil
import numpy as np
import pandas as pd
from scipy import stats

# Step1. 冻结输入清单与基本参考轴，保护已完成的其他任务
out = Path(__file__).resolve().parents[1]
root = out.parent
project = root.parents[1]
types_path = root / '00_audit/frozen_types.txt'
types = types_path.read_text().splitlines()
signature_path = project / 'CIBERSORTx_ready_20261007/09_MSCCR_Joint20v2_Fractions/GJ20v2_Signature20.txt'
signature = pd.read_csv(signature_path, sep='\t', index_col=0)
assert signature.shape == (4474, 20) and signature.columns.tolist() == types and signature.index.is_unique
assert np.isfinite(signature).all().all() and signature.ge(0).all().all()
inputs = [types_path, signature_path]
audit = [{'check': 'Frozen signature and type axis', 'status': 'PASS', 'evidence': '4474 unique genes; 20 unique types; original type order; finite nonnegative expression'}]
capacity_path = root / '01_reference/reference_type_capacity.tsv'
reference_path = root / '01_reference/reference_cells.tsv'
capacity = pd.read_csv(capacity_path, sep='\t').set_index('cell_type').loc[types]
reference = pd.read_csv(reference_path, sep='\t')
assert len(reference) == 9294 and reference.Patient.nunique() == 32 and reference.cell_id.is_unique
assert reference.groupby('coarse_identity').size().reindex(types).eq(capacity.cells).all()
inputs += [capacity_path, reference_path]
metrics = []
zero_records = []
error_pairs = []
details = []
cohorts = {}

# Step2. 逐一核对同批输入、实际抽样、作者映射和三个真值分母
for name, job in [('GSE214695', 25), ('IBDverse', 26)]:
    folder = root / '08_external_pseudobulk' / name
    review = folder / f'06_CIBERSORTx_Job{job}_review_20261014'
    coverage_path = folder / '00_audit/Full20_type_coverage.tsv'
    mapping_path = folder / '00_audit/author_type_mapping.tsv'
    member_path = folder / '02_simulation/cell_sampling_manifest.tsv.gz'
    cell_path = folder / '02_simulation/truth_cell_count_fraction.tsv'
    full_path = folder / '02_simulation/truth_full_gene_RNA_fraction.tsv'
    RNA_path = folder / '04_accuracy/truth_BP_actual_retained_gene_RNA_fraction.tsv'
    denominator_path = folder / '04_accuracy/actual_gene_denominator_audit.tsv'
    final_path = folder / '03_bayesprism/theta_final.tsv'
    initial_path = folder / '03_bayesprism/theta_initial.tsv'
    official_path = review / 'original_outputs' / f'CIBERSORTx_Job{job}_Results.txt'
    cpm_path = root / '09_CIBERSORTx_upload_small' / (name + '_pseudo_bulk100_CPM.txt')
    gene_path = folder / '03_bayesprism/retained_genes.txt'
    gene_audit_path = folder / '03_bayesprism/gene_filter_audit.tsv'
    source_gene_path = folder / '02_simulation/source_genes.tsv'
    old_metrics_path = review / 'paired_accuracy_metrics_all_scales.tsv'
    inputs += [coverage_path, mapping_path, member_path, cell_path, full_path, RNA_path, denominator_path, final_path, initial_path, official_path, cpm_path, gene_path, gene_audit_path, source_gene_path, old_metrics_path]
    coverage = pd.read_csv(coverage_path, sep='\t').set_index('cell_type').loc[types]
    active = coverage.index[coverage.evaluable].tolist()
    assert len(active) == (18 if name == 'GSE214695' else 16)
    members = pd.read_csv(member_path, sep='\t')
    mapping = pd.read_csv(mapping_path, sep='\t').set_index('AuthorLabel')
    assert mapping.index.is_unique and len(members) == 100000 and members.groupby('sample').size().eq(1000).all()
    assert not members.duplicated(['sample', 'draw_number']).any()
    assert members.AuthorLabel.map(mapping.CellType).eq(members.CellType).all()
    assert set(members.CellType) == set(active)
    cell = pd.read_csv(cell_path, sep='\t', index_col=0)
    assert cell.shape == (100, 20) and cell.index.is_unique and cell.columns.tolist() == types
    counts = members.groupby(['sample', 'CellType'], observed=True).size().unstack(fill_value=0).reindex(index=cell.index, columns=types, fill_value=0)
    assert np.allclose(counts.to_numpy() / 1000, cell, atol=1e-12, rtol=0)
    full = pd.read_csv(full_path, sep='\t', index_col=0).loc[cell.index, types]
    RNA = pd.read_csv(RNA_path, sep='\t', index_col=0).loc[cell.index, types]
    denominator = pd.read_csv(denominator_path, sep='\t', index_col=0).loc[cell.index]
    assert (denominator.retained_gene_counts > 0).all() and (denominator.full_gene_counts >= denominator.retained_gene_counts).all()
    for truth in [cell, full, RNA]:
        assert np.isfinite(truth).all().all() and truth.ge(0).all().all() and np.allclose(truth.sum(axis=1), 1, atol=1e-10, rtol=0)
        assert truth[coverage.index[~coverage.evaluable]].eq(0).all().all()
    full_read_array = full.mul(denominator.full_gene_counts, axis=0).to_numpy()
    retained_read_array = RNA.mul(denominator.retained_gene_counts, axis=0).to_numpy()
    full_integer_residual = np.abs(full_read_array - np.rint(full_read_array)).max()
    retained_integer_residual = np.abs(retained_read_array - np.rint(retained_read_array)).max()
    assert full_integer_residual < 1e-6 and retained_integer_residual < 1e-6
    official = pd.read_csv(official_path, sep='\t', index_col=0)
    assert official.shape == (100, 23) and official.index.is_unique and set(official.index) == set(cell.index)
    assert [x for x in official.columns if x in types] == types
    predictions = {'CIBERSORTx': official.loc[cell.index, types], 'BayesPrism_final': pd.read_csv(final_path, sep='\t', index_col=0).loc[cell.index, types], 'BayesPrism_initial': pd.read_csv(initial_path, sep='\t', index_col=0).loc[cell.index, types]}
    for prediction in predictions.values():
        assert prediction.shape == (100, 20) and np.isfinite(prediction).all().all() and prediction.ge(0).all().all() and np.allclose(prediction.sum(axis=1), 1, atol=1e-10, rtol=0)
    cpm = pd.read_csv(cpm_path, sep='\t', index_col=0)
    assert cpm.index.is_unique and cpm.columns.tolist() == cell.index.tolist() and np.isfinite(cpm).all().all() and cpm.ge(0).all().all()
    assert np.allclose(cpm.sum(axis=0), 1e6, atol=0.1, rtol=0)
    retained = pd.Index(gene_path.read_text().splitlines())
    source_genes = pd.Index(source_gene_path.read_text().splitlines())
    gene_audit = pd.read_csv(gene_audit_path, sep='\t').set_index('gene')
    assert retained.is_unique and source_genes.is_unique and set(retained) <= set(source_genes)
    assert set(gene_audit.index[gene_audit.retained_in_prism]) == set(retained)
    audit += [{'check': name + ' pairing, mapping and Full20 retention', 'status': 'PASS', 'evidence': f'100 matching samples, {len(active)} eligible types, no axis shift or renormalization'}, {'check': name + ' cell-count truth from actual manifest', 'status': 'PASS', 'evidence': '100000 draws reconstructed; every sample exactly 1000; exact saved cell truth'}, {'check': name + ' RNA truth algebra, denominators and gene membership', 'status': 'PASS', 'evidence': f'All row sums 1; full/retained implied type reads integer within {max(full_integer_residual, retained_integer_residual):.3g}; actual retained genes {len(retained)}'}, {'check': name + ' independent source-count reaggregation', 'status': 'PARTIAL', 'evidence': 'Source counts remain remote; local cell-level cache removed after original truth export. Reviewed generator/evaluator, source manifest and preserved denominator tables; no new remote reaggregation performed.'}, {'check': name + ' official runtime provenance', 'status': 'PARTIAL', 'evidence': 'Parameters explicitly user-confirmed; original outputs and local inputs hashed; exact web build, full log and upload receipt hash missing'}]
    read_retention = RNA.mul(denominator.retained_gene_counts, axis=0).sum() / full.mul(denominator.full_gene_counts, axis=0).sum().replace(0, np.nan)
    old = pd.read_csv(old_metrics_path, sep='\t').set_index(['truth_scale', 'method', 'cell_type'])
    cohort_metrics = []

    # 每个类型和尺度独立计算；相关恒零保留NA，不将失败当成r=0或删除
    for scale, truth in [('Retained_gene_RNA_PRIMARY', RNA), ('Cell_count_SECONDARY', cell), ('Full_gene_RNA_SECONDARY', full)]:
        for method, prediction in predictions.items():
            error = (prediction - truth) * 100
            abs_energy = error.abs().sum().sum()
            square_energy = np.square(error.to_numpy()).sum()
            for cell_type in types + ['POOLED_EVALUABLE', 'POOLED_FULL20']:
                axis = active if cell_type == 'POOLED_EVALUABLE' else types
                x = prediction[axis].to_numpy().ravel() if cell_type.startswith('POOLED') else prediction[cell_type].to_numpy()
                y = truth[axis].to_numpy().ravel() if cell_type.startswith('POOLED') else truth[cell_type].to_numpy()
                delta = (x - y) * 100
                variable = np.std(x) > 0 and np.std(y) > 0
                row = {'dataset': name, 'method': method, 'truth_scale': scale, 'cell_type': cell_type, 'evaluable': cell_type in active or cell_type == 'POOLED_EVALUABLE', 'n': len(x), 'Pearson_r': stats.pearsonr(x, y).statistic if variable else np.nan, 'Spearman_rho': stats.spearmanr(x, y).statistic if variable else np.nan, 'MAE_pp': np.abs(delta).mean(), 'RMSE_pp': np.sqrt(np.square(delta).mean()), 'bias_pp': delta.mean(), 'prediction_mean_percent': x.mean() * 100, 'truth_mean_percent': y.mean() * 100, 'slope_prediction_on_truth': np.cov(y, x, ddof=0)[0, 1] / np.var(y) if np.var(y) > 0 else np.nan, 'absolute_error_share_Full20': np.abs(delta).sum() / abs_energy if not cell_type.startswith('POOLED') else np.nan, 'squared_error_share_Full20': np.square(delta).sum() / square_energy if not cell_type.startswith('POOLED') else np.nan}
                original_method = f'CIBERSORTx_Job{job}' if method == 'CIBERSORTx' else method
                prior = old.loc[(scale, original_method, cell_type)]
                assert np.allclose([row[k] for k in ['MAE_pp', 'RMSE_pp', 'bias_pp', 'Pearson_r', 'Spearman_rho']], [prior[k] for k in ['MAE_pp', 'RMSE_pp', 'mean_signed_error_pp', 'Pearson_r', 'Spearman_rho']], atol=1e-10, rtol=1e-10, equal_nan=True)
                metrics.append(row)
                cohort_metrics.append(row)
            corr = error.corr()
            for cell_type in ['Plasma', 'ILC', 'Stem_TA_progenitor', 'Tuft']:
                peers = corr[cell_type].drop(cell_type).dropna().sort_values().head(3)
                for peer, correlation in peers.items():
                    error_pairs.append({'dataset': name, 'method': method, 'truth_scale': scale, 'cell_type': cell_type, 'other_type': peer, 'error_correlation': correlation, 'interpretation': 'Descriptive closure-dependent error covariance; not proof of read reassignment'})
    for method, prediction in predictions.items():
        for cell_type in types:
            y, x = cell[cell_type].to_numpy(), prediction[cell_type].to_numpy()
            n_positive, n_absent = int((y > 0).sum()), int((y == 0).sum())
            zero_records.append({'dataset': name, 'method': method, 'cell_type': cell_type, 'n_true_positive': n_positive, 'n_true_absent': n_absent, 'false_negative_exact_zero': int(((y > 0) & (x == 0)).sum()), 'false_positive_exact_positive': int(((y == 0) & (x > 0)).sum()), 'false_positive_above_0p01pp': int(((y == 0) & (x > 1e-4)).sum()), 'mean_prediction_when_absent_percent': float(x[y == 0].mean() * 100) if n_absent else np.nan, 'max_prediction_when_absent_percent': float(x[y == 0].max() * 100) if n_absent else np.nan, 'retained_fraction_of_type_RNA': read_retention[cell_type]})
    saved = review / 'zero_prediction_audit.tsv'
    old_zero = pd.read_csv(saved, sep='\t').set_index('cell_type')
    inputs.append(saved)
    for row in zero_records[-60:]:
        if row['method'] == 'CIBERSORTx':
            assert row['false_negative_exact_zero'] == old_zero.loc[row['cell_type'], 'CX_zero_when_cell_truth_positive']
    audit.append({'check': name + ' all reported metrics and exact-zero counts', 'status': 'PASS', 'evidence': '198 method/scale/type or pooled metric rows reproduced within 1e-10; exact-zero false negatives reproduced'})
    cohorts[name] = {'coverage': coverage, 'members': members, 'cell': cell, 'full': full, 'RNA': RNA, 'predictions': predictions, 'denominator': denominator, 'retained': retained, 'gene_audit': gene_audit, 'cpm': cpm, 'mapping': mapping}
    print(name, 'PASS: exact cell truth, axes, denominator algebra and 198 prior metric rows')

# Step3. 异常谱系的marker、共线性和基因过滤，仅作已有材料诊断
metric_table = pd.DataFrame(metrics)
zero_table = pd.DataFrame(zero_records)
metric_table.to_csv(out / 'EXTERNAL_RECOMPUTED_METRICS.tsv', sep='\t', index=False)
zero_table.to_csv(out / 'DETECTION_AND_GENE_DENOMINATORS.tsv', sep='\t', index=False)
pd.DataFrame(error_pairs).to_csv(out / 'COMPOSITION_ERROR_COVARIANCE.tsv', sep='\t', index=False)
pair_rows = []
log_corr = np.log1p(signature).corr()
raw_corr = signature.corr()
unit = signature.to_numpy() / np.sqrt(np.square(signature.to_numpy()).sum(axis=0))
condition = float(np.linalg.cond(unit))
for i, cell_type in enumerate(types):
    for j in range(i + 1, len(types)):
        peer = types[j]
        pair_rows.append({'type_1': cell_type, 'type_2': peer, 'Pearson_raw_signature': raw_corr.loc[cell_type, peer], 'Pearson_log1p_signature': log_corr.loc[cell_type, peer], 'cosine_signature': float(unit[:, i] @ unit[:, j]), 'unit_L2_column_condition_number': condition})
pd.DataFrame(pair_rows).to_csv(out / 'SIGNATURE_COLLINEARITY.tsv', sep='\t', index=False)
audit.append({'check': 'Exact signature alias or rank loss', 'status': 'PASS', 'evidence': f'Full column rank {np.linalg.matrix_rank(unit)}/20; unit-L2 condition number {condition:.6f}; pairwise similarities descriptive, not unique mechanism evidence'})
marker_rows = []
for cell_type, markers in {'Plasma': ['IGKC', 'IGHA1', 'IGLC2', 'IGLC3', 'JCHAIN', 'MZB1'], 'ILC': ['IL7R', 'KIT', 'NCR2', 'RORC'], 'Stem_TA_progenitor': ['OLFM4', 'LGR5', 'MKI67'], 'Tuft': ['POU2F3', 'TRPM5', 'SH2D6']}.items():
    for gene in markers:
        for name, data in cohorts.items():
            a = data['gene_audit']
            marker_rows.append({'dataset': name, 'cell_type': cell_type, 'gene': gene, 'in_signature': gene in signature.index, 'signature_expression': signature.loc[gene, cell_type] if gene in signature.index else np.nan, 'in_clean_BP_reference': bool(a.loc[gene, 'in_frozen_clean_reference']) if gene in a.index else False, 'in_BP_prism': gene in data['retained'], 'filter_stage': 'new.prism mixture outlier/align step' if gene in a.index and a.loc[gene, 'in_frozen_clean_reference'] and gene not in data['retained'] else 'retained or not in input'})
pd.DataFrame(marker_rows).to_csv(out / 'KEY_MARKER_AUDIT.tsv', sep='\t', index=False)

# Step4. MSCCR120逐样本和类型复算，保留Job17/20重复差异
annotation_path = root / '00_audit/MSCCR120_annotation.tsv'
bulk_final_path = root / '02_bayesprism_joint20_MSCCR120/theta_final.tsv'
bulk_initial_path = root / '02_bayesprism_joint20_MSCCR120/theta_initial.tsv'
job20_path = project / 'CIBERSORTx_signature_audit/MSCCR_Joint20v2_Job20_review_20261009/original_outputs/CIBERSORTx_Job20_Adjusted.txt'
job17_path = project / 'CIBERSORTx_signature_audit/MSCCR_Joint20v2_Job17_review_20261008/original_outputs/CIBERSORTx_Job17_Adjusted.txt'
prior_bulk_path = root / '04_MSCCR120_concordance/type_concordance_metrics.tsv'
sample_qc_path = root / '02_bayesprism_joint20_MSCCR120/sample_QC.tsv'
inputs += [annotation_path, bulk_final_path, bulk_initial_path, job20_path, job17_path, prior_bulk_path, sample_qc_path]
annotation = pd.read_csv(annotation_path, sep='\t').set_index('ID')
assert len(annotation) == 120 and annotation.PatientID.nunique() == 120 and annotation.groupby('Region3').size().eq(40).all()
bulk = {'BayesPrism_final': pd.read_csv(bulk_final_path, sep='\t', index_col=0).loc[annotation.index, types], 'BayesPrism_initial': pd.read_csv(bulk_initial_path, sep='\t', index_col=0).loc[annotation.index, types], 'CIBERSORTx_Job20': pd.read_csv(job20_path, sep='\t', index_col=0).loc[annotation.index, types], 'CIBERSORTx_Job17': pd.read_csv(job17_path, sep='\t', index_col=0).loc[annotation.index, types]}
for frame in bulk.values():
    assert frame.shape == (120, 20) and np.isfinite(frame).all().all() and frame.ge(0).all().all() and np.allclose(frame.sum(axis=1), 1, atol=1e-10)
bulk_rows = []
for label, left, right in [('BP_final_vs_Job20', 'BayesPrism_final', 'CIBERSORTx_Job20'), ('BP_initial_vs_Job20', 'BayesPrism_initial', 'CIBERSORTx_Job20'), ('Job20_vs_Job17', 'CIBERSORTx_Job20', 'CIBERSORTx_Job17')]:
    for cell_type in types:
        x, y = bulk[left][cell_type], bulk[right][cell_type]
        delta = (x - y) * 100
        bulk_rows.append({'comparison': label, 'cell_type': cell_type, 'Pearson_r': stats.pearsonr(x, y).statistic, 'Spearman_rho': stats.spearmanr(x, y).statistic, 'bias_pp': delta.mean(), 'MAE_pp': delta.abs().mean(), 'median_absolute_difference_pp': delta.abs().median(), 'RMSE_pp': np.sqrt(np.square(delta).mean()), 'BA_low_pp': delta.mean() - 1.96 * delta.std(ddof=1), 'BA_high_pp': delta.mean() + 1.96 * delta.std(ddof=1)})
bulk_metrics = pd.DataFrame(bulk_rows)
prior = pd.read_csv(prior_bulk_path, sep='\t').query("group == 'All120'").set_index('cell_type').loc[types]
current = bulk_metrics.loc[bulk_metrics.comparison.eq('BP_final_vs_Job20')].set_index('cell_type').loc[types]
assert np.allclose(current[['Pearson_r', 'Spearman_rho', 'bias_pp', 'RMSE_pp']].to_numpy(), prior[['Pearson_r', 'Spearman_rho', 'mean_difference_pp', 'RMSE_pp']].to_numpy(), atol=1e-10)
delta = (bulk['BayesPrism_final'] - bulk['CIBERSORTx_Job20']).abs()
tv = delta.sum(axis=1) / 2
repeat_tv = (bulk['CIBERSORTx_Job20'] - bulk['CIBERSORTx_Job17']).abs().sum(axis=1) / 2
bulk_summary = {'Pearson_median_across_types': current.Pearson_r.median(), 'Spearman_median_across_types': current.Spearman_rho.median(), 'mean_absolute_difference_pp': delta.to_numpy().mean() * 100, 'median_sample_total_variation': tv.median(), 'Job20_vs_Job17_mean_absolute_difference_pp': float((bulk['CIBERSORTx_Job20'] - bulk['CIBERSORTx_Job17']).abs().to_numpy().mean() * 100), 'Job20_vs_Job17_median_total_variation': repeat_tv.median(), 'n_independent_patients': 120}
bulk_metrics.to_csv(out / 'MSCCR120_RECOMPUTED_CONCORDANCE.tsv', sep='\t', index=False)

# Step5. 最小分层方向/效应比较，不新建临床因果模型或对p值筛选
effect_rows = []
for scope in ['All120', 'Colon', 'Rectum', 'Ileum']:
    ids = annotation.index if scope == 'All120' else annotation.index[annotation.Region3.eq(scope)]
    for field, group1, group0 in [('Disease', 'UC', 'CD'), ('Inflammation', 'Inflamed', 'Non-inflamed')]:
        group1_ids = annotation.loc[ids].index[annotation.loc[ids, field].eq(group1)]
        group0_ids = annotation.loc[ids].index[annotation.loc[ids, field].eq(group0)]
        assert len(group1_ids) > 0 and len(group0_ids) > 0
        for cell_type in types:
            values = {m: (frame.loc[group1_ids, cell_type].mean() - frame.loc[group0_ids, cell_type].mean()) * 100 for m, frame in bulk.items()}
            effect_rows.append({'scope': scope, 'contrast': group1 + '_minus_' + group0, 'cell_type': cell_type, 'n_group1': len(group1_ids), 'n_group0': len(group0_ids), **values, 'BP_final_Job20_sign_agree': np.sign(values['BayesPrism_final']) == np.sign(values['CIBERSORTx_Job20']), 'Job20_Job17_sign_agree': np.sign(values['CIBERSORTx_Job20']) == np.sign(values['CIBERSORTx_Job17']), 'meaning': 'Descriptive mean contrast in pp; not adjusted clinical or causal effect; near-zero signs can differ trivially'})
for group1, group0 in [('Ileum', 'Colon'), ('Rectum', 'Colon')]:
    ids1, ids0 = annotation.index[annotation.Region3.eq(group1)], annotation.index[annotation.Region3.eq(group0)]
    for cell_type in types:
        values = {m: (frame.loc[ids1, cell_type].mean() - frame.loc[ids0, cell_type].mean()) * 100 for m, frame in bulk.items()}
        effect_rows.append({'scope': 'Region3_marginal', 'contrast': group1 + '_minus_' + group0, 'cell_type': cell_type, 'n_group1': len(ids1), 'n_group0': len(ids0), **values, 'BP_final_Job20_sign_agree': np.sign(values['BayesPrism_final']) == np.sign(values['CIBERSORTx_Job20']), 'Job20_Job17_sign_agree': np.sign(values['CIBERSORTx_Job20']) == np.sign(values['CIBERSORTx_Job17']), 'meaning': 'Marginal anatomical contrast; disease/inflammation distributions differ; not anatomical matching experiment'})
effects = pd.DataFrame(effect_rows)
effects.to_csv(out / 'MSCCR120_DESCRIPTIVE_EFFECTS.tsv', sep='\t', index=False)
pd.crosstab([annotation.Region3, annotation.Disease], annotation.Inflammation).to_csv(out / 'MSCCR120_METADATA_COVERAGE.tsv', sep='\t')
ref_regions = reference.groupby(['coarse_identity', 'Site'], observed=True).agg(cells=('cell_id', 'size'), donors=('Patient', 'nunique'), profiles=('profile_id', 'nunique')).reset_index()
ref_regions.to_csv(out / 'REFERENCE_ANATOMICAL_COVERAGE.tsv', sep='\t', index=False)
sample_qc = pd.read_csv(sample_qc_path, sep='\t')
fraction_bound = sample_qc.fractional_raw_entries / sample_qc.mapped_raw_counts
bulk_summary.update({'raw_library_min': float(sample_qc.full_raw_library_counts.min()), 'raw_library_max': float(sample_qc.full_raw_library_counts.max()), 'samples_with_fractional_counts': int(sample_qc.fractional_raw_entries.gt(0).sum()), 'fractional_entries': int(sample_qc.fractional_raw_entries.sum()), 'max_pre_prism_read_loss_upper_bound_fraction_if_truncated': float(fraction_bound.max()), 'pooled_pre_prism_read_loss_upper_bound_fraction_if_truncated': float(sample_qc.fractional_raw_entries.sum() / sample_qc.mapped_raw_counts.sum())})
audit += [{'check': 'MSCCR120 core concordance and Job17 repeat', 'status': 'PASS', 'evidence': '120 distinct patients, 40 per Region3; all reported r/rho/bias/RMSE reproduced; repeat differences retained'}, {'check': 'MSCCR120 fractional count origin and sampler precision', 'status': 'PARTIAL', 'evidence': 'Cached GEO metadata distinguishes counts from adjcounts and describes featureCounts raw counts; selected path/hash and library-scale totals support count input. Specific fractional-assignment mechanism is not documented; BP Gibbs uses rmultinom integer sizes, so fractional input preservation is not proof of exact fractional sampling.'}, {'check': 'Job20 S-mode source representation', 'status': 'PARTIAL', 'evidence': 'Same training source, 212 profiles each represented by 5 of 1060 training cells, all 20 types; source subset is not literal 212-profile signature construction input; no new correction efficacy experiment'}, {'check': 'Scope and clinical effect interpretation', 'status': 'PASS', 'evidence': 'Only descriptive contrasts; no new signature, BP fit, S-mode, NMF, discovery or remote SSH'}]

# Step6. 输入哈希与轻量输出，独立来源复核不假装已经重新读取远程counts
for folder in [root / 'scripts', root / '08_external_pseudobulk/GSE214695/03_bayesprism', root / '08_external_pseudobulk/IBDverse/03_bayesprism']:
    inputs += [p for p in folder.iterdir() if p.is_file() and p.suffix in ['.py', '.R', '.sh', '.txt']]
inputs += [root / '00_audit/server_input_manifest.tsv', root / '00_audit/local_input_manifest.tsv', root / '00_audit/signature_reference_provenance.tsv', root / '01_reference/reference_build_audit.json', project / 'CIBERSORTx_ready_20261008/01_MSCCR_platform120_pilot/核对记录_勿上传/source_subset_cell_manifest.tsv', project / 'CIBERSORTx_ready_20261008/01_MSCCR_platform120_pilot/核对记录_勿上传/source_subset_summary.json', project / 'CIBERSORTx_ready_20261007/06_Joint20v2_DonorSite_reference/核对记录_勿上传/profile_manifest.tsv', project / 'CIBERSORTx_signature_audit/MSCCR_Joint20v2_Job20_review_20261009/original_outputs/Job20_runtime_excerpt_user_provided.txt', root.parent / 'paper_aligned_recovery_characterization_20261012/00_audit/PMC8526411_bioc.xml', root.parent / 'paper_aligned_recovery_characterization_20261012/00_audit/prior_frozen_states.tsv', root.parent / 'paper_aligned_recovery_characterization_20261012/00_audit/prior_frozen_ecotypes.tsv']
input_rows = []
for path in sorted(set(inputs)):
    assert path.exists(), str(path)
    input_rows.append({'path': str(path), 'bytes': path.stat().st_size, 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()})
pd.DataFrame(input_rows).to_csv(out / 'INPUT_MANIFEST.tsv', sep='\t', index=False)
pd.DataFrame(audit).to_csv(out / 'AUDIT_CHECKS.tsv', sep='\t', index=False)
(out / 'numeric_audit_summary.json').write_text(json.dumps({'external_metric_rows_reproduced': len(metrics), 'external_type_rows': 120, 'signature_unit_L2_condition_number': condition, 'MSCCR120': bulk_summary, 'all_checks': audit, 'status': 'NUMERICAL_AUDIT_PASS_WITH_PROVENANCE_LIMITATIONS', 'new_remote_connections': 0, 'new_deconvolution_runs': 0}, indent=2), encoding='utf-8')
print(json.dumps(bulk_summary, indent=2))
print('All core saved numerical results reproduced. No prior input or output modified.')
