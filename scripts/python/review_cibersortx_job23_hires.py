# 审阅Job23第二块HiRes并对照Job21：核对命名修正、120例fraction和共同基因的表达一致性。
# 在本机cibersortx_review环境运行；保存下载原件、逐类型和逐基因审计表，不合并正式G矩阵、不启动EcoTyper。
from pathlib import Path
from zipfile import ZipFile
import hashlib
import json
import shutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Step1. 归档并安全解压，原件已有时只核对，不覆盖不同内容
project = Path('D:/IBD_ecotype/IBD_ecotype_analysis')
source_zip = Path('C:/Users/马泽凯/Downloads/CIBERSORTx_Job23_output.zip')
inputs = project / 'CIBERSORTx_ready_20261008/04_MSCCR120_HiRes_gene_batches'
prior = project / 'CIBERSORTx_signature_audit/MSCCR_Joint20v2_Job21_HiRes_review_20261009/original_outputs'
outdir = project / 'CIBERSORTx_signature_audit/MSCCR_Joint20v2_Job23_HiRes_review_20261009'
originals = outdir / 'original_outputs'
originals.mkdir(parents=True, exist_ok=True)
zip_hash = hashlib.sha256(source_zip.read_bytes()).hexdigest()
archived_zip = outdir / source_zip.name
if archived_zip.exists():
    assert hashlib.sha256(archived_zip.read_bytes()).hexdigest() == zip_hash
else:
    shutil.copy2(source_zip, archived_zip)
integrity = []
with ZipFile(archived_zip) as archive:
    for item in archive.infolist():
        path = (originals / item.filename).resolve()
        assert path.is_relative_to(originals.resolve()), f'ZIP路径越界：{item.filename}'
        if item.is_dir():
            continue
        content = archive.read(item.filename)
        if path.exists():
            assert path.read_bytes() == content
        else:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(content)
        integrity.append({'file': item.filename, 'bytes': len(content), 'sha256': hashlib.sha256(content).hexdigest()})
pd.DataFrame(integrity).to_csv(outdir / 'original_integrity.tsv', sep='\t', index=False)
runtime = (originals / 'CIBERSORTx_Job23_runtime_log.txt').read_text(encoding='utf-8')
error = (originals / 'CIBERSORTx_Job23_error_log.txt').read_text(encoding='utf-8')
runtime21 = (prior / 'CIBERSORTx_Job21_runtime_log.txt').read_text(encoding='utf-8')
assert 'Imputing high-resolution cell type GEPs...done.' in runtime
assert 'MSCCR120_HiRes_B02_WebNames_Genes1249.txt' in runtime and 'rmbatchSmode: TRUE' in runtime and 'QN: FALSE' in runtime
options23 = [line for line in runtime.splitlines() if line.startswith('>[Options]')]
options21 = [line for line in runtime21.splitlines() if line.startswith('>[Options]')]
(outdir / 'runtime_options_comparison.txt').write_text('Job21\n' + '\n'.join(options21) + '\n\nJob23\n' + '\n'.join(options23) + '\n', encoding='utf-8')
print('Job23完成；警告：', error.strip())

# Step2. 核对冻结120例和20类，测量两块重新估计的fraction差异
annotation = pd.read_csv(inputs / '核对记录_勿上传/MSCCR120_annotation.tsv', sep='\t')
signature = pd.read_csv(inputs / 'GJ20v2_Signature20.txt', sep='\t', index_col=0)
types = signature.columns.tolist()
sample_ids = annotation['ID'].tolist()
fractions23 = pd.read_csv(originals / 'CIBERSORTxGEP_Job23_Fractions-Adjusted.txt', sep='\t', index_col=0)
fractions21 = pd.read_csv(prior / 'CIBERSORTxGEP_Job21_Fractions-Adjusted.txt', sep='\t', index_col=0)
for fractions in [fractions23, fractions21]:
    assert fractions.index.tolist() == sample_ids and fractions.columns[:20].tolist() == types
    assert np.isfinite(fractions[types].to_numpy()).all() and fractions[types].ge(0).all().all()
    assert np.allclose(fractions[types].sum(axis=1), 1, atol=1e-8)
fraction_change = fractions23[types] - fractions21[types]
reallocation = fraction_change.abs().sum(axis=1) * 50  # 两个总和为1的分布之间需要重新分配的比例，单位百分比。
fraction_type_summary = pd.DataFrame({'Job21_mean_percent': fractions21[types].mean() * 100, 'Job23_mean_percent': fractions23[types].mean() * 100, 'change_percentage_points': fraction_change.mean() * 100, 'median_absolute_change_percentage_points': fraction_change.abs().median() * 100})
fraction_type_summary.to_csv(outdir / 'Job23_vs_Job21_fraction_type_summary.tsv', sep='\t', index_label='CellType')
fraction_change.mul(100).to_csv(outdir / 'Job23_vs_Job21_fraction_change_percentage_points.tsv', sep='\t', index_label='sample')
pd.DataFrame({'ID': sample_ids, 'fraction_reallocation_percent': reallocation.to_numpy()}).merge(annotation, on='ID', validate='one_to_one').to_csv(outdir / 'sample_fraction_reallocation.tsv', sep='\t', index=False)
fractions23.reset_index().rename(columns={'Mixture': 'ID'}).merge(annotation, on='ID', validate='one_to_one').to_csv(outdir / 'Job23_fraction_with_annotation.tsv', sep='\t', index=False)
print('Fraction重新分配百分比，中位数/最大值：', reallocation.median(), reallocation.max())

# Step3. 比较上游Adjusted signature及Group表达，不能把它们当作样本级G
adjusted23 = pd.read_csv(originals / 'CIBERSORTx_Job23_sigmatrix_Adjusted.txt', sep='\t', index_col=0)
adjusted21 = pd.read_csv(prior / 'CIBERSORTx_Job21_sigmatrix_Adjusted.txt', sep='\t', index_col=0)
assert adjusted23.index.tolist() == adjusted21.index.tolist() and adjusted23.columns.tolist() == adjusted21.columns.tolist() == types
adjusted_difference = (adjusted23 - adjusted21).abs()
adjusted_rows = []
for cell_type in types:
    x = adjusted21[cell_type].to_numpy()
    y = adjusted23[cell_type].to_numpy()
    adjusted_rows.append({'CellType': cell_type, 'n_values_changed_over_1e_8': int((np.abs(x - y) > 1e-8).sum()), 'max_absolute_change': float(np.abs(x - y).max()), 'log2_pearson': float(np.corrcoef(np.log2(x + 1), np.log2(y + 1))[0, 1])})
pd.DataFrame(adjusted_rows).to_csv(outdir / 'adjusted_signature_comparison.tsv', sep='\t', index=False)
group23 = pd.read_csv(originals / 'CIBERSORTxGEP_Job23_GEPs_Filtered.txt', sep='\t', index_col=0)
group21 = pd.read_csv(prior / 'CIBERSORTxGEP_Job21_GEPs_Filtered.txt', sep='\t', index_col=0)
assert group23.shape == (17040, 20) and group23.index.tolist() == group21.index.tolist()
upstream_rows = []
for suffix in ['GEPs.txt', 'GEPs_Filtered.txt', 'Mixture-Adjusted.txt', 'Weights.txt']:
    upstream21 = pd.read_csv(prior / f'CIBERSORTxGEP_Job21_{suffix}', sep='\t', index_col=0)
    upstream23 = pd.read_csv(originals / f'CIBERSORTxGEP_Job23_{suffix}', sep='\t', index_col=0)
    assert upstream21.index.tolist() == upstream23.index.tolist() and upstream21.columns.tolist() == upstream23.columns.tolist()
    upstream_rows.append({'output': suffix, 'rows': upstream23.shape[0], 'columns': upstream23.shape[1], 'numeric_identical': bool(np.array_equal(upstream21.to_numpy(), upstream23.to_numpy(), equal_nan=True))})
pd.DataFrame(upstream_rows).to_csv(outdir / 'upstream_numeric_equality.tsv', sep='\t', index=False)

# Step4. 核对修改后的第二块所有请求基因是否实际加载
request = pd.read_csv(inputs / '先提交这两批/MSCCR120_HiRes_B02_WebNames_Genes1249.txt')['GeneSymbol'].tolist()
anchors = pd.read_csv(inputs / '核对记录_勿上传/共同100基因_选择依据.tsv', sep='\t')
assert len(request) == len(set(request)) == 1249 and len(anchors) == 100
type_rows = []
gene_rows = []
comparison_rows = []
common_type_rows = []
diagnostic_rows = []
output_gene_order = None

# Step5. 逐类型审阅表达，并按同一基因的120例样本比较两块
# 1/NA的含义沿用ZIP README；只对两块都完整且有变异的基因计算相关性，常数不能算作重复成功。
for cell_type in types:
    expression23 = pd.read_csv(originals / f'CIBERSORTxHiRes_Job23_{cell_type}_Window60.txt', sep='\t', index_col=0)
    expression21 = pd.read_csv(prior / f'CIBERSORTxHiRes_Job21_{cell_type}_Window60.txt', sep='\t', index_col=0)
    assert expression23.columns.tolist() == expression21.columns.tolist() == sample_ids and expression23.index.is_unique
    if output_gene_order is None:
        output_gene_order = expression23.index.tolist()
    assert expression23.index.tolist() == output_gene_order
    values = expression23.to_numpy(dtype=float)
    assert not np.isinf(values).any() and np.nanmin(values) >= 0
    all_one23 = expression23.eq(1).all(axis=1)
    all_na23 = expression23.isna().all(axis=1)
    complete23 = expression23.notna().all(axis=1)
    span23 = expression23.max(axis=1) - expression23.min(axis=1)
    varying23 = complete23 & span23.gt(1e-8)
    all_one21 = expression21.eq(1).all(axis=1)
    all_na21 = expression21.isna().all(axis=1)
    complete21 = expression21.notna().all(axis=1)
    varying21 = complete21 & (expression21.max(axis=1) - expression21.min(axis=1)).gt(1e-8)
    own_anchors = anchors.loc[anchors['selection_type'].eq(cell_type), 'GeneSymbol'].tolist()
    type_rows.append({'CellType': cell_type, 'Job21_output_genes': len(expression21), 'Job23_output_genes': len(expression23), 'Job21_complete_varying_genes': int(varying21.sum()), 'Job23_complete_varying_genes': int(varying23.sum()), 'Job23_all_one_genes': int(all_one23.sum()), 'Job23_all_NA_genes': int(all_na23.sum()), 'Job23_partial_NA_genes': int((~complete23 & ~all_na23).sum()), 'Job23_stable_non_one_genes': int((complete23 & ~all_one23 & ~varying23).sum()), 'Job23_own_five_anchors_varying': int(varying23.reindex(own_anchors).sum()), 'Job21_group_positive_all17040': int(group21[cell_type].gt(0).sum()), 'Job23_group_positive_all17040': int(group23[cell_type].gt(0).sum())})
    for gene in expression23.index:
        gene_rows.append({'CellType': cell_type, 'GeneSymbol': gene, 'all_one': bool(all_one23[gene]), 'all_NA': bool(all_na23[gene]), 'complete': bool(complete23[gene]), 'varying': bool(varying23[gene]), 'mean': float(expression23.loc[gene].mean()), 'span': float(span23[gene])})
    common_genes = [gene for gene in expression21.index if gene in expression23.index]
    assert set(common_genes) == set(anchors['GeneSymbol']), '两块实际共同基因应恰好为100个共同诊断基因。'
    these_rows = []
    for gene in common_genes:
        x = expression21.loc[gene].to_numpy(dtype=float)
        y = expression23.loc[gene].to_numpy(dtype=float)
        both_varying = bool(varying21[gene] and varying23[gene])
        pearson_log = np.nan
        pearson_linear = np.nan
        median_symmetric_difference = np.nan
        if both_varying:
            pearson_log = float(np.corrcoef(np.log2(x + 1), np.log2(y + 1))[0, 1])
            pearson_linear = float(np.corrcoef(x, y)[0, 1])
            median_symmetric_difference = float(np.median(2 * np.abs(x - y) / np.maximum(np.abs(x) + np.abs(y), 1e-12)) * 100)
        state21 = 'all_NA' if all_na21[gene] else 'all_one' if all_one21[gene] else 'varying' if varying21[gene] else 'other'
        state23 = 'all_NA' if all_na23[gene] else 'all_one' if all_one23[gene] else 'varying' if varying23[gene] else 'other'
        row = {'CellType': cell_type, 'GeneSymbol': gene, 'selection_type': anchors.set_index('GeneSymbol').loc[gene, 'selection_type'], 'state_Job21': state21, 'state_Job23': state23, 'same_state': state21 == state23, 'both_varying': both_varying, 'pearson_log2': pearson_log, 'pearson_linear': pearson_linear, 'median_symmetric_difference_percent': median_symmetric_difference, 'Job21_mean': float(np.nanmean(x)) if not np.isnan(x).all() else np.nan, 'Job23_mean': float(np.nanmean(y)) if not np.isnan(y).all() else np.nan}
        for job, data in [('Job21', x), ('Job23', y)]:
            finite = data[np.isfinite(data)]
            row[f'{job}_CV'] = float(finite.std(ddof=1) / finite.mean()) if len(finite) > 1 and finite.mean() > 0 else np.nan
            row[f'{job}_log2_span'] = float(np.log2(finite.max() + 1) - np.log2(finite.min() + 1)) if len(finite) else np.nan
        these_rows.append(row)
        comparison_rows.append(row)
    shared = pd.DataFrame(these_rows)
    valid = shared.loc[shared['both_varying']]
    own_shared = shared.loc[shared['selection_type'].eq(cell_type)]
    common_type_rows.append({'CellType': cell_type, 'common_genes': len(common_genes), 'Job21_varying_common_genes': int(varying21.reindex(common_genes).sum()), 'Job23_varying_common_genes': int(varying23.reindex(common_genes).sum()), 'both_varying_common_genes': len(valid), 'common_state_agreement_percent': float(shared['same_state'].mean() * 100), 'median_pearson_log2_both_varying': float(valid['pearson_log2'].median()), 'min_pearson_log2_both_varying': float(valid['pearson_log2'].min()), 'median_symmetric_difference_percent_both_varying': float(valid['median_symmetric_difference_percent'].median()), 'both_varying_genes_pearson_below_0_8': int(valid['pearson_log2'].lt(0.8).sum()), 'own_five_same_state': int(own_shared['same_state'].sum()), 'own_five_both_varying': int(own_shared['both_varying'].sum())})
    for gene in ['CD8A', 'CD8B', 'CD3D', 'VWF', 'CLDN5', 'SOX10', 'PLP1', 'DEFA5', 'DEFA6', 'MKI67', 'OLFM4', 'LGR5', 'MUC2', 'CHGA', 'POU2F3']:
        if gene in expression23.index:
            diagnostic_rows.append({'CellType': cell_type, 'GeneSymbol': gene, 'state_Job23': 'all_NA' if all_na23[gene] else 'all_one' if all_one23[gene] else 'varying' if varying23[gene] else 'other', 'mean': float(expression23.loc[gene].mean()), 'span': float(span23[gene])})
type_summary = pd.DataFrame(type_rows)
common_summary = pd.DataFrame(common_type_rows)
comparison = pd.DataFrame(comparison_rows)
missing_genes = [gene for gene in request if gene not in output_gene_order]
unexpected_genes = [gene for gene in output_gene_order if gene not in request]
assert len(output_gene_order) == 1249 and not missing_genes and not unexpected_genes
type_summary.to_csv(outdir / 'HiRes_type_usability_comparison.tsv', sep='\t', index=False)
pd.DataFrame(gene_rows).to_csv(outdir / 'HiRes_gene_status.tsv', sep='\t', index=False)
comparison.to_csv(outdir / 'shared100_gene_expression_comparison.tsv', sep='\t', index=False)
comparison.loc[comparison['CellType'].eq(comparison['selection_type'])].to_csv(outdir / 'own_type_100_anchors_comparison.tsv', sep='\t', index=False)
common_summary.to_csv(outdir / 'shared100_type_consistency_summary.tsv', sep='\t', index=False)
pd.DataFrame(diagnostic_rows).to_csv(outdir / 'diagnostic_marker_status.tsv', sep='\t', index=False)
print(type_summary[['CellType', 'Job21_complete_varying_genes', 'Job23_complete_varying_genes', 'Job23_own_five_anchors_varying']].to_string(index=False))
print(common_summary[['CellType', 'both_varying_common_genes', 'median_pearson_log2_both_varying', 'common_state_agreement_percent']].to_string(index=False))

# Step6. 保存图和审阅边界，当前两块不导出为正式合并表达
fig, axes = plt.subplots(1, 2, figsize=(13, 7))
positions = np.arange(len(types))
axes[0].barh(positions - 0.18, type_summary['Job21_complete_varying_genes'], height=0.36, label='Job21 / block 1', color='#94b7c4')
axes[0].barh(positions + 0.18, type_summary['Job23_complete_varying_genes'], height=0.36, label='Job23 / block 2', color='#267997')
axes[0].set_yticks(positions, types)
axes[0].invert_yaxis()
axes[0].set_xlabel('Complete nonconstant genes (different block genes)')
axes[0].legend(loc='lower right', fontsize=8)
axes[1].barh(positions, common_summary['median_pearson_log2_both_varying'], color='#267997')
axes[1].set_yticks(positions, types)
axes[1].invert_yaxis()
axes[1].set_xlim(-1, 1)
axes[1].axvline(0.8, color='#b96537', linestyle='--', linewidth=1)
axes[1].set_xlabel('Median log-expression Pearson r (varying in both blocks)')
fig.suptitle('Job23 vs Job21: output availability and shared-gene consistency')
fig.tight_layout()
fig.savefig(outdir / 'Job23_vs_Job21_HiRes_comparison.png', dpi=170)
plt.close(fig)
valid = comparison.loc[comparison['both_varying']]
summary = {'status': 'numeric_review_complete', 'source_zip_sha256': zip_hash, 'archive_members': len(integrity), 'runtime_completed': True, 'warning': error.strip(), 'requested_genes': len(request), 'hires_loaded_genes': len(output_gene_order), 'missing_genes': missing_genes, 'unexpected_genes': unexpected_genes, 'corrected_batch02_request_matching_verified': True, 'samples': len(sample_ids), 'types': len(types), 'window': 60, 'fraction_P_value_unique': fractions23['P-value'].unique().tolist(), 'fraction_correlation_median': float(fractions23['Correlation'].median()), 'fraction_RMSE_median': float(fractions23['RMSE'].median()), 'median_fraction_reallocation_vs_Job21_percent': float(reallocation.median()), 'max_fraction_reallocation_vs_Job21_percent': float(reallocation.max()), 'adjusted_signature_max_absolute_difference': float(adjusted_difference.max().max()), 'adjusted_signature_values_changed_over_1e_8': int(adjusted_difference.gt(1e-8).sum().sum()), 'shared_gene_type_pairs': len(comparison), 'shared_gene_type_pairs_both_varying': len(valid), 'shared_both_varying_median_pearson_log2': float(valid['pearson_log2'].median()), 'shared_state_agreement_percent': float(comparison['same_state'].mean() * 100), 'shared_both_varying_median_symmetric_difference_percent': float(valid['median_symmetric_difference_percent'].median()), 'all_one_types': type_summary.loc[type_summary['Job23_all_one_genes'].eq(len(output_gene_order)), 'CellType'].tolist(), 'full_transcriptome_hires_expression_available': False, 'formal_merged_matrix_exported': False, 'ecotyper_started': False}
(outdir / 'review_summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps(summary, ensure_ascii=False, indent=2))

# Step7. 为可能的第一块重跑准备命名清单，不改变表达矩阵或原B01，也不代替用户提交任务
request1 = pd.read_csv(inputs / '先提交这两批/MSCCR120_HiRes_B01_Genes1249.txt')['GeneSymbol'].tolist()
request1_web = [gene.replace('-', '.') for gene in request1]
assert len(request1_web) == len(set(request1_web)) == 1249 and set(request1_web).issubset(group23.index)
request1_text = 'GeneSymbol\n' + '\n'.join(request1_web) + '\n'
corrected1 = inputs / '先提交这两批/MSCCR120_HiRes_B01_WebNames_Genes1249.txt'
if corrected1.exists():
    assert corrected1.read_text(encoding='utf-8') == request1_text
else:
    corrected1.write_text(request1_text, encoding='utf-8')
pd.DataFrame({'original_gene': request1, 'web_request_gene': request1_web}).to_csv(outdir / 'Batch01_fixed_gene_name_mapping.tsv', sep='\t', index=False)
print('已准备修正命名的B01；尚未运行：', corrected1)
