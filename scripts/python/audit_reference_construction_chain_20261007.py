# 审计三个参考从抽样 counts 到网页 inferred_refsample、signature 的表达变化，不修改输入。
# 在本地 .tools/cibersortx_review 环境、D:/IBD_ecotype 下运行；输出构建链条、marker 和训练内诊断。
# NNLS 只用于定位表达匹配问题，不替代 CIBERSORTx，也不算独立留出验证。
from pathlib import Path
import hashlib
import json
import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.optimize import nnls
from audit_cibersortx_joint_signature import PANELS

# Step1. 固定已有参考和审计输出，保留原始文件身份
root = Path(__file__).resolve().parents[2]
pilots = root / 'IBDverse_reference/three_reference_pilots_20261006/run_20261006_183211'
audits = root / 'CIBERSORTx_signature_audit/three_candidate_comparison_20261006'
out = root / 'CIBERSORTx_signature_audit/construction_chain_20261007'
out.mkdir(parents=True, exist_ok=True)
candidate_names = [('gse_joint', 'gse_joint'), ('gse_ileum', 'gse_ileum'), ('scp259', 'scp_rectum')]
chain_records = []
marker_records = []
identity_records = []
integrity_records = []
candidate_records = []

# 三个候选按同一口径逐个处理；按类型循环是分组计算所需。
for pilot_name, audit_name in candidate_names:
    source = pilots / pilot_name
    prior = audits / audit_name
    destination = out / pilot_name
    destination.mkdir(exist_ok=True)
    meta = json.loads((source / 'summary.json').read_text(encoding='utf-8'))
    manifest = pd.read_csv(source / 'selected_cell_manifest.tsv', sep='\t')
    genes = pd.read_csv(source / 'gene_export_manifest.tsv', sep='\t', keep_default_na=False)
    signature_path = next(prior.glob('*_signature_gene_symbols_restored.txt'))
    reference_path = next(prior.glob('*_inferred_refsample_gene_symbols_restored.txt'))
    signature = pd.read_csv(signature_path, sep='\t', index_col=0)
    inferred = pd.read_csv(reference_path, sep='\t', index_col=0)
    source_paths = [source / 'selected_counts.h5', source / 'selected_cell_manifest.tsv', source / 'gene_export_manifest.tsv', source / meta['web_input_filename'], signature_path, reference_path]
    before = {}
    for path in source_paths:
        with path.open('rb') as handle:
            before[str(path)] = hashlib.file_digest(handle, 'sha256').hexdigest()

    # Step2. 从只读 H5 重建上传 counts，并逐个数值核对网页输入
    with h5py.File(source / 'selected_counts.h5', 'r') as handle:
        counts = csr_matrix((handle['counts/data'][:], handle['counts/indices'][:], handle['counts/indptr'][:]), shape=tuple(handle['counts'].attrs['shape']))
        assert handle['cell_id'].asstr()[:].tolist() == manifest.cell_id.tolist()
        assert handle['web_column'].asstr()[:].tolist() == manifest.web_column.tolist()
        assert handle['gene_symbol'].asstr()[:].tolist() == genes.gene_symbol.tolist()
    web_mask = genes.web_included.to_numpy(bool)
    web_genes = genes.loc[web_mask, 'gene_symbol'].tolist()
    exported = counts[:, web_mask].toarray().astype(float)
    web = pd.read_csv(source / meta['web_input_filename'], sep='\t', index_col=0)
    assert web.index.tolist() == web_genes == inferred.index.tolist()
    assert web.columns.tolist() == manifest.web_column.tolist()
    assert np.array_equal(exported.T, web.to_numpy())
    assert np.isfinite(exported).all() and (exported >= 0).all()
    del web
    original_libraries = np.asarray(counts.sum(axis=1)).ravel()
    exported_libraries = exported.sum(axis=1)
    assert np.allclose(original_libraries, manifest.extracted_total_counts.to_numpy(), rtol=0, atol=.01)
    equal_cell_cpm = exported / exported_libraries[:, None] * 1e6
    type_names = signature.columns.tolist()
    assert set(type_names) == set(manifest.pilot_type)
    reference_labels = inferred.columns.str.replace(r'\.\d+$', '', regex=True)
    raw_mean = pd.DataFrame(index=web_genes, columns=type_names, dtype=float)
    raw_pooled = raw_mean.copy()
    inferred_mean = raw_mean.copy()
    detection = raw_mean.copy()
    mito_mask = genes.gene_symbol.str.startswith('MT-').to_numpy()
    raw_mito_counts = np.asarray(counts[:, mito_mask].sum(axis=1)).ravel()
    inferred_all_zero = (inferred.to_numpy() == 0).all(axis=1)

    # Step3. 检查类型覆盖、供者 RNA 权重、表达过滤和重新归一化
    for cell_type in type_names:
        use = manifest.pilot_type.eq(cell_type).to_numpy()
        selected = manifest.loc[use].copy()
        raw_mean[cell_type] = equal_cell_cpm[use].mean(axis=0)
        raw_pooled[cell_type] = exported[use].sum(axis=0) / exported_libraries[use].sum() * 1e6
        detection[cell_type] = (exported[use] > 0).mean(axis=0)
        inferred_mean[cell_type] = inferred.loc[:, reference_labels == cell_type].mean(axis=1)
        lost = inferred_mean[cell_type].eq(0).to_numpy() & (raw_mean[cell_type].to_numpy() > 0)
        retained = inferred_mean[cell_type].to_numpy() > 0
        selected['export_library_counts'] = exported_libraries[use]
        donor = selected.groupby('individual_id', observed=True).agg(n_cells=('cell_id', 'size'), counts=('export_library_counts', 'sum'))
        donor['count_share'] = donor['counts'] / donor['counts'].sum()
        donor.to_csv(destination / f'{cell_type}_donor_counts_and_RNA_weight.tsv', sep='\t')
        mito_percent = raw_mito_counts[use] / original_libraries[use] * 100
        signature_locations = raw_mean.index.get_indexer(signature.index)
        source_on_signature = raw_mean[cell_type].to_numpy()[signature_locations]
        inferred_on_signature = inferred_mean[cell_type].to_numpy()[signature_locations]
        active = (source_on_signature > 0) & (inferred_on_signature > 0)
        ratios = inferred_on_signature[active] / source_on_signature[active]
        record = {'candidate': pilot_name, 'cell_type': cell_type, 'n_cells': int(use.sum()), 'n_donors': int(len(donor)),
            'n_author_labels': int(selected.author_label.nunique()), 'site_counts': json.dumps(selected.Site.value_counts().to_dict(), ensure_ascii=False),
            'author_label_counts': json.dumps(selected.author_label.value_counts().to_dict(), ensure_ascii=False),
            'max_donor_cell_share': float(donor.n_cells.max() / use.sum()), 'max_donor_export_count_share': float(donor.count_share.max()),
            'median_library_counts': float(np.median(original_libraries[use])), 'median_detected_genes': float(selected.extracted_n_genes.median()),
            'median_mito_percent': float(np.median(mito_percent)), 'n_mito_above20pct_descriptive_only': int((mito_percent > 20).sum()),
            'export_retained_original_count_fraction': float(exported_libraries[use].sum() / original_libraries[use].sum()),
            'source_detected_genes': int((raw_mean[cell_type] > 0).sum()), 'inferred_nonzero_genes': int(retained.sum()),
            'source_expressed_but_inferred_zero_genes': int(lost.sum()),
            'lost_equal_cell_CPM_mass_fraction': float(raw_mean[cell_type].to_numpy()[lost].sum() / 1e6),
            'lost_pooled_count_CPM_mass_fraction': float(raw_pooled[cell_type].to_numpy()[lost].sum() / 1e6),
            'active_inferred_vs_source_median_ratio_on_signature': float(np.median(ratios)) if len(ratios) else None,
            'source_vs_inferred_pearson_on_signature': float(np.corrcoef(source_on_signature, inferred_on_signature)[0, 1]),
            'source_vs_inferred_log1p_pearson_on_signature': float(np.corrcoef(np.log1p(source_on_signature), np.log1p(inferred_on_signature))[0, 1]),
            'signature_column_sum_subset': float(signature[cell_type].sum())}
        chain_records.append(record)
        # 20% 仅是描述分布的标记，未据此删除细胞，也不是统一的肠道 QC 阈值。
        for gene in PANELS[cell_type]:
            in_source = gene in raw_mean.index
            row = {'candidate': pilot_name, 'cell_type': cell_type, 'gene': gene, 'in_source_export': in_source, 'in_signature': gene in signature.index}
            if in_source:
                row.update(source_mean_raw_count=float(exported[use, web_genes.index(gene)].mean()), source_detection_fraction=float(detection.loc[gene, cell_type]),
                    source_equal_cell_CPM=float(raw_mean.loc[gene, cell_type]), source_pooled_count_CPM=float(raw_pooled.loc[gene, cell_type]),
                    inferred_CPM=float(inferred_mean.loc[gene, cell_type]), inferred_highest_type=inferred_mean.loc[gene].idxmax())
            if gene in signature.index:
                row.update(signature_expression=float(signature.at[gene, cell_type]), signature_highest_type=signature.loc[gene].idxmax())
            marker_records.append(row)
    raw_mean.to_csv(destination / 'source_equal_cell_mean_CPM.tsv', sep='\t')
    raw_pooled.to_csv(destination / 'source_pooled_counts_CPM.tsv', sep='\t')
    inferred_mean.to_csv(destination / 'web_inferred_type_mean_CPM.tsv', sep='\t')
    gene_loss = pd.DataFrame({'gene': web_genes, 'all_web_inferred_profiles_zero': inferred_all_zero, 'max_source_detection_fraction': detection.max(axis=1).to_numpy(), 'max_source_mean_CPM': raw_mean.max(axis=1).to_numpy()})
    gene_loss.to_csv(destination / 'inferred_gene_loss_diagnostics.tsv', sep='\t', index=False)
    detection_bins = pd.cut(gene_loss.max_source_detection_fraction, [0, .1, .25, .5, .75, 1], include_lowest=True)
    gene_loss.groupby(detection_bins, observed=True).agg(n_genes=('gene', 'size'), inferred_all_zero=('all_web_inferred_profiles_zero', 'sum')).to_csv(destination / 'loss_by_source_detection_bin.tsv', sep='\t')

    # Step4. 用 NNLS 做训练内表达匹配诊断；与网页算法、独立留出结果分别保存
    matrix = signature.to_numpy(float)
    matrix_scale = np.max(matrix)
    for profile_name, profiles in [('ideal_signature_self_check', signature), ('training_equal_cell_mean', raw_mean.loc[signature.index]), ('training_pooled_counts', raw_pooled.loc[signature.index])]:
        fitted = pd.DataFrame(index=type_names, columns=type_names, dtype=float)
        for cell_type in type_names:
            coefficient, residual = nnls(matrix / matrix_scale, profiles[cell_type].to_numpy() / matrix_scale)
            assert coefficient.sum() > 0
            fractions = coefficient / coefficient.sum()
            fitted[cell_type] = fractions
            identity_records.append({'candidate': pilot_name, 'diagnostic': profile_name, 'true_training_label': cell_type,
                'target_fraction': float(fractions[type_names.index(cell_type)]), 'highest_predicted_type': type_names[int(np.argmax(fractions))],
                'highest_fraction': float(fractions.max()), 'scaled_residual_norm': float(residual), 'independent_validation': False})
        fitted.index.name = 'predicted_type'
        fitted.to_csv(destination / f'NNLS_{profile_name}.tsv', sep='\t')
    candidate_records.append({'candidate': pilot_name, 'n_cells': len(manifest), 'n_types': len(type_names), 'web_genes': len(web_genes),
        'all_inferred_zero_genes': int(inferred_all_zero.sum()), 'uploaded_counts_exactly_match_H5': True,
        'raw_count_max': float(exported.max()), 'inferred_column_sum_min': float(inferred.sum(axis=0).min()), 'inferred_column_sum_max': float(inferred.sum(axis=0).max())})
    for path in source_paths:
        with path.open('rb') as handle:
            after = hashlib.file_digest(handle, 'sha256').hexdigest()
        assert before[str(path)] == after
        integrity_records.append({'path': str(path), 'sha256_before': before[str(path)], 'sha256_after': after, 'unchanged': True})
    print(pilot_name, 'counts exact; cells/types/genes:', len(manifest), len(type_names), len(web_genes), flush=True)

# Step5. 保存跨参考可读结果与明确的诊断边界
chain = pd.DataFrame(chain_records)
markers = pd.DataFrame(marker_records)
identity = pd.DataFrame(identity_records)
chain.to_csv(out / 'type_construction_chain.tsv', sep='\t', index=False)
markers.to_csv(out / 'marker_construction_chain.tsv', sep='\t', index=False)
identity.to_csv(out / 'NNLS_training_diagnostics.tsv', sep='\t', index=False)
pd.DataFrame(candidate_records).to_csv(out / 'candidate_structure_checks.tsv', sep='\t', index=False)
pd.DataFrame(integrity_records).to_csv(out / 'input_integrity.tsv', sep='\t', index=False)
summary = {'status': 'construction_chain_audit_completed', 'originals_unchanged': True, 'candidates': candidate_records,
    'training_NNLS_is_not_holdout_validation': True, 'no_new_web_jobs_or_remote_runs': True,
    'full_heldout_single_cell_counts_available_locally': False, 'final_holdout_used': False}
(out / 'audit_summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
print(chain.loc[chain.cell_type.isin(['Paneth', 'Stem_TA_progenitor', 'Goblet', 'Tuft', 'ILC']), ['candidate', 'cell_type', 'median_mito_percent', 'lost_equal_cell_CPM_mass_fraction', 'active_inferred_vs_source_median_ratio_on_signature']].to_string(index=False))
print('Output:', out)
