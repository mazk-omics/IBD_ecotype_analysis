# 审阅 Job13 供者/部位伪 bulk signature，比较旧 Job4 的 marker、类型区分信息与输入汇总。
# 在本机 cibersortx_review 环境运行；读取下载文件和冻结 v2 参考，写本地审阅目录，不登录服务器、不改参考或真值。
# 这里的“5 倍特异基因”是矩阵数值诊断，不能冒充网页 DEG/q-value，也不能替代留出比例验证。
from pathlib import Path
import hashlib
import json
import shutil
import numpy as np
import pandas as pd
from pypdf import PdfReader
import pypdfium2 as pdfium

# Step1. 保存网页原件，并核对同一组 20 类型
root = Path('D:/IBD_ecotype/IBD_ecotype_analysis')
out = root / 'CIBERSORTx_signature_audit/Joint20v2_Job13_review_20261007'
original = out / 'original_outputs'
original.mkdir(parents=True, exist_ok=True)
stem = 'CIBERSORTx_Job13_GJ20v2_PhenotypeClasses.CIBERSORTx_Job13_GJ20v2_DonorSite_CPM.bm.K999'
integrity = []
for extension in ['txt', 'pdf']:
    source = Path('C:/Users/马泽凯/Downloads') / f'{stem}.{extension}'
    destination = original / source.name
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    if destination.exists():
        assert hashlib.sha256(destination.read_bytes()).hexdigest() == digest, '不覆盖内容不同的网页原件'
    else:
        shutil.copyfile(source, destination)
    integrity.append({'source': str(source), 'archived_file': source.name, 'bytes': source.stat().st_size, 'sha256': digest})
pd.DataFrame(integrity).to_csv(out / 'original_integrity.tsv', sep='\t', index=False)
new = pd.read_csv(original / f'{stem}.txt', sep='\t', index_col=0)
old = pd.read_csv(root / 'CIBERSORTx_signature_audit/three_candidate_comparison_20261006/gse_joint/GSE282122_Joint20_signature_gene_symbols_restored.txt', sep='\t', index_col=0)
reference_dir = root / 'CIBERSORTx_ready_20261007/06_Joint20v2_DonorSite_reference'
reference = pd.read_csv(reference_dir / 'GJ20v2_DonorSite_CPM.txt', sep='\t', index_col=0)
profiles = pd.read_csv(reference_dir / '核对记录_勿上传/profile_manifest.tsv', sep='\t')
assert new.index.is_unique and new.columns.is_unique and set(new.columns) == set(old.columns) == set(profiles.coarse_identity)
assert np.isfinite(new.to_numpy()).all() and new.ge(0).all().all()
assert new.index.isin(reference.index).all()
old = old.loc[:, new.columns]
print('Old/new signature shape:', old.shape, new.shape, flush=True)

# Step2. 按名字核对经典 marker，避免用全图红色面积判断某个基因是否存在
markers = ['DEFA5', 'DEFA6', 'REG3A', 'REG1A', 'REG1B', 'LYZ', 'PRSS2', 'PLA2G2A', 'MMP7', 'REG3G', 'ITLN2', 'SOX9', 'MUC2', 'TFF3', 'SPINK4', 'LGR5', 'OLFM4', 'MKI67', 'TOP2A', 'CHGA', 'CHGB', 'NEUROG3', 'SOX10', 'S100B', 'POU2F3', 'TRPM5']
marker_rows = []
for gene in markers:
    row = {'gene': gene, 'in_v2_input': gene in reference.index, 'in_Job4_signature': gene in old.index, 'in_Job13_signature': gene in new.index}
    for version, matrix in [('Job4', old), ('Job13', new)]:
        row[f'{version}_Paneth'] = matrix.at[gene, 'Paneth'] if gene in matrix.index else np.nan
        row[f'{version}_highest_type'] = matrix.loc[gene].idxmax() if gene in matrix.index else 'NOT_SELECTED'
    marker_rows.append(row)
marker_comparison = pd.DataFrame(marker_rows)
marker_comparison.to_csv(out / 'marker_comparison.tsv', sep='\t', index=False, na_rep='NA')
print(marker_comparison.head(11).to_string(index=False), flush=True)

# Step3. 衡量矩阵中可区分各类型的基因数，保存完整判定依据
# 最高列、2 倍和 5 倍比值只是相同数值规则；不是统计显著性或每类型网页实际分配的 DEG 数。
diagnostics = []
for version, matrix in [('Job4', old), ('Job13', new)]:
    values = matrix.to_numpy()
    for cell_type in matrix.columns:
        max_other = matrix.drop(columns=cell_type).max(axis=1)
        ratio = matrix[cell_type] / max_other
        ranked = pd.DataFrame({'expression': matrix[cell_type], 'max_other_expression': max_other, 'ratio_max_other': ratio})
        ranked.sort_values('ratio_max_other', ascending=False).to_csv(out / f'{version}_{cell_type}_gene_specificity.tsv', sep='\t', index_label='gene')
        diagnostics.append({'version': version, 'cell_type': cell_type, 'highest_expression_genes': int(matrix.idxmax(axis=1).eq(cell_type).sum()), 'ratio_ge2_genes': int(ratio.ge(2).sum()), 'ratio_ge5_genes': int(ratio.ge(5).sum()), 'top3_fraction_of_selected_column_sum': float(matrix[cell_type].nlargest(3).sum() / matrix[cell_type].sum()), 'raw_condition_number': float(np.linalg.cond(values)), 'unit_L2_column_condition_number': float(np.linalg.cond(values / np.linalg.norm(values, axis=0)))})
diagnostics = pd.DataFrame(diagnostics)
diagnostics.to_csv(out / 'type_diagnostics_old_new.tsv', sep='\t', index=False)
gene_sets = pd.DataFrame({'gene': old.index.union(new.index)})
gene_sets['in_Job4'] = gene_sets.gene.isin(old.index)
gene_sets['in_Job13'] = gene_sets.gene.isin(new.index)
gene_sets.to_csv(out / 'signature_gene_membership.tsv', sep='\t', index=False)

# Step4. 反查参考与输出数值，核对是否有类型错位或非 log 转换异常
# 检查“输入每类型逐基因 median+1”的数值吻合；这是本次文件证据，不宣称读过网页源码。
aggregation_checks = []
for cell_type in new.columns:
    columns = profiles.loc[profiles.coarse_identity.eq(cell_type), 'profile_id']
    type_reference = reference.loc[new.index, columns]
    predicted_median = type_reference.median(axis=1) + 1
    predicted_mean = type_reference.mean(axis=1) + 1
    aggregation_checks.append({'cell_type': cell_type, 'input_profiles': len(columns), 'median_plus1_max_abs_error': float((predicted_median - new[cell_type]).abs().max()), 'mean_plus1_max_abs_error': float((predicted_mean - new[cell_type]).abs().max()), 'matches_median_plus1': bool(np.allclose(predicted_median, new[cell_type], rtol=1e-6, atol=1e-6))})
aggregation_checks = pd.DataFrame(aggregation_checks)
aggregation_checks.to_csv(out / 'reference_to_signature_aggregation_check.tsv', sep='\t', index=False)
assert aggregation_checks.matches_median_plus1.all()
paneth_profiles = profiles.loc[profiles.coarse_identity.eq('Paneth'), ['profile_id', 'Patient', 'Disease', 'region', 'Site', 'sample_id', 'selected_cells']].copy()
paneth_markers = reference.loc[['DEFA5', 'DEFA6', 'LYZ', 'REG3A', 'MUC2', 'OLFM4'], paneth_profiles.profile_id].T
paneth_markers.index.name = 'profile_id'
paneth_profiles = paneth_profiles.merge(paneth_markers.reset_index(), on='profile_id', validate='one_to_one')
paneth_profiles.to_csv(out / 'Paneth_six_reference_profiles.tsv', sep='\t', index=False)
paneth_profiles.groupby('region')[['DEFA5', 'DEFA6', 'LYZ', 'REG3A']].agg(['mean', 'median', 'min', 'max']).to_csv(out / 'Paneth_region_marker_summary.tsv', sep='\t')

# Step5. 保存类型相关性，说明强 marker 保留和多类型区分信息下降可以同时发生
pair_rows = []
common_genes = old.index.intersection(new.index)
for a, b in [('Paneth', 'Goblet'), ('Paneth', 'Stem_TA_progenitor'), ('Goblet', 'Stem_TA_progenitor'), ('Stem_TA_progenitor', 'Colonic_absorptive'), ('Stem_TA_progenitor', 'Ileal_absorptive')]:
    pair_rows.append({'type_a': a, 'type_b': b, 'Job4_log2p1_pearson': float(np.log2(old[[a, b]] + 1).corr().iloc[0, 1]), 'Job13_log2p1_pearson': float(np.log2(new[[a, b]] + 1).corr().iloc[0, 1]), 'common_genes': len(common_genes), 'Job4_common_log2p1_pearson': float(np.log2(old.loc[common_genes, [a, b]] + 1).corr().iloc[0, 1]), 'Job13_common_log2p1_pearson': float(np.log2(new.loc[common_genes, [a, b]] + 1).corr().iloc[0, 1])})
pd.DataFrame(pair_rows).to_csv(out / 'epithelial_pair_similarity.tsv', sep='\t', index=False)

# Step6. 阅读完整单页 PDF；只渲染审阅，不修改原图
pdf = PdfReader(original / f'{stem}.pdf')
(out / 'Job13_pdf_text.txt').write_text('\n'.join(page.extract_text() or '' for page in pdf.pages), encoding='utf-8')
document = pdfium.PdfDocument(str(original / f'{stem}.pdf'))
for page_number in range(len(document)):
    document[page_number].render(scale=1.4).to_pil().save(out / f'Job13_page_{page_number + 1}.png')

# Step7. 归档为新候选，并保存可以复核的结果摘要
candidate_dir = root / 'signature_matrix/gse_joint_v2'
candidate_dir.mkdir(exist_ok=True)
for extension in ['txt', 'pdf']:
    destination = candidate_dir / f'{stem}.{extension}'
    if destination.exists():
        assert destination.read_bytes() == (original / f'{stem}.{extension}').read_bytes()
    else:
        shutil.copyfile(original / f'{stem}.{extension}', destination)
summary = {'Job4_shape': list(old.shape), 'Job13_shape': list(new.shape), 'gene_overlap': len(old.index.intersection(new.index)), 'types_unchanged': True, 'all_markers_DEFA5_DEFA6_REG3A_REG1A_LYZ_retained': all(g in new.index for g in ['DEFA5', 'DEFA6', 'REG3A', 'REG1A', 'LYZ']), 'all_type_values_match_input_median_plus1': bool(aggregation_checks.matches_median_plus1.all()), 'old_raw_condition_number': float(diagnostics.loc[diagnostics.version.eq('Job4'), 'raw_condition_number'].iloc[0]), 'new_raw_condition_number': float(diagnostics.loc[diagnostics.version.eq('Job13'), 'raw_condition_number'].iloc[0]), 'pdf_pages': len(pdf.pages), 'web_runtime_log_available': False, 'web_DEG_qvalues_available': False, 'deconvolution_validation_done': False, 'server_access_this_review': False}
(out / 'audit_summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
print('Review complete:', out, flush=True)
