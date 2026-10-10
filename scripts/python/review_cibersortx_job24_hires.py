# 审阅Job24第一块重跑并对照Job23第二块：核对命名修正、120例fraction和共同基因的表达一致性。
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
source_zip = Path('C:/Users/马泽凯/Downloads/CIBERSORTx_Job24_output.zip')
inputs = project / 'CIBERSORTx_ready_20261008/04_MSCCR120_HiRes_gene_batches'
prior = project / 'CIBERSORTx_signature_audit/MSCCR_Joint20v2_Job23_HiRes_review_20261009/original_outputs'
outdir = project / 'CIBERSORTx_signature_audit/MSCCR_Joint20v2_Job24_HiRes_review_20261009'
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
runtime = (originals / 'CIBERSORTx_Job24_runtime_log.txt').read_text(encoding='utf-8')
error = (originals / 'CIBERSORTx_Job24_error_log.txt').read_text(encoding='utf-8')
runtime23 = (prior / 'CIBERSORTx_Job23_runtime_log.txt').read_text(encoding='utf-8')
assert 'Imputing high-resolution cell type GEPs...done.' in runtime
assert 'MSCCR120_HiRes_B01_WebNames_Genes1249.txt' in runtime and 'rmbatchSmode: TRUE' in runtime and 'QN: FALSE' in runtime
options24 = [line for line in runtime.splitlines() if line.startswith('>[Options]')]
options23 = [line for line in runtime23.splitlines() if line.startswith('>[Options]')]
(outdir / 'runtime_options_comparison.txt').write_text('Job23\n' + '\n'.join(options23) + '\n\nJob24\n' + '\n'.join(options24) + '\n', encoding='utf-8')
print('Job24完成；警告：', error.strip())

# Step2. 核对冻结120例和20类，测量两块重新估计的fraction差异
annotation = pd.read_csv(inputs / '核对记录_勿上传/MSCCR120_annotation.tsv', sep='\t')
signature = pd.read_csv(inputs / 'GJ20v2_Signature20.txt', sep='\t', index_col=0)
types = signature.columns.tolist()
sample_ids = annotation['ID'].tolist()
fractions24 = pd.read_csv(originals / 'CIBERSORTxGEP_Job24_Fractions-Adjusted.txt', sep='\t', index_col=0)
fractions23 = pd.read_csv(prior / 'CIBERSORTxGEP_Job23_Fractions-Adjusted.txt', sep='\t', index_col=0)
for fractions in [fractions24, fractions23]:
    assert fractions.index.tolist() == sample_ids and fractions.columns[:20].tolist() == types
    assert np.isfinite(fractions[types].to_numpy()).all() and fractions[types].ge(0).all().all()
    assert np.allclose(fractions[types].sum(axis=1), 1, atol=1e-8)
fraction_change = fractions24[types] - fractions23[types]
reallocation = fraction_change.abs().sum(axis=1) * 50  # 两个总和为1的分布之间需要重新分配的比例，单位百分比。
fraction_type_summary = pd.DataFrame({'Job23_mean_percent': fractions23[types].mean() * 100, 'Job24_mean_percent': fractions24[types].mean() * 100, 'change_percentage_points': fraction_change.mean() * 100, 'median_absolute_change_percentage_points': fraction_change.abs().median() * 100})
fraction_type_summary.to_csv(outdir / 'Job24_vs_Job23_fraction_type_summary.tsv', sep='\t', index_label='CellType')
fraction_change.mul(100).to_csv(outdir / 'Job24_vs_Job23_fraction_change_percentage_points.tsv', sep='\t', index_label='sample')
pd.DataFrame({'ID': sample_ids, 'fraction_reallocation_percent': reallocation.to_numpy()}).merge(annotation, on='ID', validate='one_to_one').to_csv(outdir / 'sample_fraction_reallocation.tsv', sep='\t', index=False)
fractions24.reset_index().rename(columns={'Mixture': 'ID'}).merge(annotation, on='ID', validate='one_to_one').to_csv(outdir / 'Job24_fraction_with_annotation.tsv', sep='\t', index=False)
print('Fraction重新分配百分比，中位数/最大值：', reallocation.median(), reallocation.max())

# Step3. 比较上游Adjusted signature及Group表达，不能把它们当作样本级G
adjusted24 = pd.read_csv(originals / 'CIBERSORTx_Job24_sigmatrix_Adjusted.txt', sep='\t', index_col=0)
adjusted23 = pd.read_csv(prior / 'CIBERSORTx_Job23_sigmatrix_Adjusted.txt', sep='\t', index_col=0)
assert adjusted24.index.tolist() == adjusted23.index.tolist() and adjusted24.columns.tolist() == adjusted23.columns.tolist() == types
adjusted_difference = (adjusted24 - adjusted23).abs()
adjusted_rows = []
for cell_type in types:
    x = adjusted23[cell_type].to_numpy()
    y = adjusted24[cell_type].to_numpy()
    adjusted_rows.append({'CellType': cell_type, 'n_values_changed_over_1e_8': int((np.abs(x - y) > 1e-8).sum()), 'max_absolute_change': float(np.abs(x - y).max()), 'log2_pearson': float(np.corrcoef(np.log2(x + 1), np.log2(y + 1))[0, 1])})
pd.DataFrame(adjusted_rows).to_csv(outdir / 'adjusted_signature_comparison.tsv', sep='\t', index=False)
group24 = pd.read_csv(originals / 'CIBERSORTxGEP_Job24_GEPs_Filtered.txt', sep='\t', index_col=0)
group23 = pd.read_csv(prior / 'CIBERSORTxGEP_Job23_GEPs_Filtered.txt', sep='\t', index_col=0)
assert group24.shape == (17040, 20) and group24.index.tolist() == group23.index.tolist()
upstream_rows = []
for suffix in ['GEPs.txt', 'GEPs_Filtered.txt', 'Mixture-Adjusted.txt', 'Weights.txt']:
    upstream23 = pd.read_csv(prior / f'CIBERSORTxGEP_Job23_{suffix}', sep='\t', index_col=0)
    upstream24 = pd.read_csv(originals / f'CIBERSORTxGEP_Job24_{suffix}', sep='\t', index_col=0)
    assert upstream23.index.tolist() == upstream24.index.tolist() and upstream23.columns.tolist() == upstream24.columns.tolist()
    upstream_rows.append({'output': suffix, 'rows': upstream24.shape[0], 'columns': upstream24.shape[1], 'numeric_identical': bool(np.array_equal(upstream23.to_numpy(), upstream24.to_numpy(), equal_nan=True))})
pd.DataFrame(upstream_rows).to_csv(outdir / 'upstream_numeric_equality.tsv', sep='\t', index=False)

# Step4. 核对修改后的命名修正后的第一块所有请求基因是否实际加载
request = pd.read_csv(inputs / '先提交这两批/MSCCR120_HiRes_B01_WebNames_Genes1249.txt')['GeneSymbol'].tolist()
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
    expression24 = pd.read_csv(originals / f'CIBERSORTxHiRes_Job24_{cell_type}_Window60.txt', sep='\t', index_col=0)
    expression23 = pd.read_csv(prior / f'CIBERSORTxHiRes_Job23_{cell_type}_Window60.txt', sep='\t', index_col=0)
    assert expression24.columns.tolist() == expression23.columns.tolist() == sample_ids and expression24.index.is_unique
    if output_gene_order is None:
        output_gene_order = expression24.index.tolist()
    assert expression24.index.tolist() == output_gene_order
    values = expression24.to_numpy(dtype=float)
    assert not np.isinf(values).any() and np.nanmin(values) >= 0
    all_one24 = expression24.eq(1).all(axis=1)
    all_na24 = expression24.isna().all(axis=1)
    complete24 = expression24.notna().all(axis=1)
    span24 = expression24.max(axis=1) - expression24.min(axis=1)
    varying24 = complete24 & span24.gt(1e-8)
    all_one23 = expression23.eq(1).all(axis=1)
    all_na23 = expression23.isna().all(axis=1)
    complete23 = expression23.notna().all(axis=1)
    varying23 = complete23 & (expression23.max(axis=1) - expression23.min(axis=1)).gt(1e-8)
    own_anchors = anchors.loc[anchors['selection_type'].eq(cell_type), 'GeneSymbol'].tolist()
    type_rows.append({'CellType': cell_type, 'Job23_output_genes': len(expression23), 'Job24_output_genes': len(expression24), 'Job23_complete_varying_genes': int(varying23.sum()), 'Job24_complete_varying_genes': int(varying24.sum()), 'Job24_all_one_genes': int(all_one24.sum()), 'Job24_all_NA_genes': int(all_na24.sum()), 'Job24_partial_NA_genes': int((~complete24 & ~all_na24).sum()), 'Job24_stable_non_one_genes': int((complete24 & ~all_one24 & ~varying24).sum()), 'Job24_own_five_anchors_varying': int(varying24.reindex(own_anchors).sum()), 'Job23_group_positive_all17040': int(group23[cell_type].gt(0).sum()), 'Job24_group_positive_all17040': int(group24[cell_type].gt(0).sum())})
    for gene in expression24.index:
        gene_rows.append({'CellType': cell_type, 'GeneSymbol': gene, 'all_one': bool(all_one24[gene]), 'all_NA': bool(all_na24[gene]), 'complete': bool(complete24[gene]), 'varying': bool(varying24[gene]), 'mean': float(expression24.loc[gene].mean()), 'span': float(span24[gene])})
    common_genes = [gene for gene in expression23.index if gene in expression24.index]
    assert set(common_genes) == set(anchors['GeneSymbol']), '两块实际共同基因应恰好为100个共同诊断基因。'
    these_rows = []
    for gene in common_genes:
        x = expression23.loc[gene].to_numpy(dtype=float)
        y = expression24.loc[gene].to_numpy(dtype=float)
        both_varying = bool(varying23[gene] and varying24[gene])
        pearson_log = np.nan
        pearson_linear = np.nan
        median_symmetric_difference = np.nan
        if both_varying:
            pearson_log = float(np.corrcoef(np.log2(x + 1), np.log2(y + 1))[0, 1])
            pearson_linear = float(np.corrcoef(x, y)[0, 1])
            median_symmetric_difference = float(np.median(2 * np.abs(x - y) / np.maximum(np.abs(x) + np.abs(y), 1e-12)) * 100)
        state23 = 'all_NA' if all_na23[gene] else 'all_one' if all_one23[gene] else 'varying' if varying23[gene] else 'other'
        state24 = 'all_NA' if all_na24[gene] else 'all_one' if all_one24[gene] else 'varying' if varying24[gene] else 'other'
        row = {'CellType': cell_type, 'GeneSymbol': gene, 'selection_type': anchors.set_index('GeneSymbol').loc[gene, 'selection_type'], 'state_Job23': state23, 'state_Job24': state24, 'same_state': state23 == state24, 'both_varying': both_varying, 'pearson_log2': pearson_log, 'pearson_linear': pearson_linear, 'median_symmetric_difference_percent': median_symmetric_difference, 'Job23_mean': float(np.nanmean(x)) if not np.isnan(x).all() else np.nan, 'Job24_mean': float(np.nanmean(y)) if not np.isnan(y).all() else np.nan}
        for job, data in [('Job23', x), ('Job24', y)]:
            finite = data[np.isfinite(data)]
            row[f'{job}_CV'] = float(finite.std(ddof=1) / finite.mean()) if len(finite) > 1 and finite.mean() > 0 else np.nan
            row[f'{job}_log2_span'] = float(np.log2(finite.max() + 1) - np.log2(finite.min() + 1)) if len(finite) else np.nan
        these_rows.append(row)
        comparison_rows.append(row)
    shared = pd.DataFrame(these_rows)
    valid = shared.loc[shared['both_varying']]
    own_shared = shared.loc[shared['selection_type'].eq(cell_type)]
    common_type_rows.append({'CellType': cell_type, 'common_genes': len(common_genes), 'Job23_varying_common_genes': int(varying23.reindex(common_genes).sum()), 'Job24_varying_common_genes': int(varying24.reindex(common_genes).sum()), 'both_varying_common_genes': len(valid), 'common_state_agreement_percent': float(shared['same_state'].mean() * 100), 'median_pearson_log2_both_varying': float(valid['pearson_log2'].median()), 'min_pearson_log2_both_varying': float(valid['pearson_log2'].min()), 'median_symmetric_difference_percent_both_varying': float(valid['median_symmetric_difference_percent'].median()), 'both_varying_genes_pearson_below_0_8': int(valid['pearson_log2'].lt(0.8).sum()), 'own_five_same_state': int(own_shared['same_state'].sum()), 'own_five_both_varying': int(own_shared['both_varying'].sum())})
    for gene in ['CD8A', 'CD8B', 'CD3D', 'VWF', 'CLDN5', 'SOX10', 'PLP1', 'DEFA5', 'DEFA6', 'MKI67', 'OLFM4', 'LGR5', 'MUC2', 'CHGA', 'POU2F3']:
        if gene in expression24.index:
            diagnostic_rows.append({'CellType': cell_type, 'GeneSymbol': gene, 'state_Job24': 'all_NA' if all_na24[gene] else 'all_one' if all_one24[gene] else 'varying' if varying24[gene] else 'other', 'mean': float(expression24.loc[gene].mean()), 'span': float(span24[gene])})
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
print(type_summary[['CellType', 'Job23_complete_varying_genes', 'Job24_complete_varying_genes', 'Job24_own_five_anchors_varying']].to_string(index=False))
print(common_summary[['CellType', 'both_varying_common_genes', 'median_pearson_log2_both_varying', 'common_state_agreement_percent']].to_string(index=False))

# Step6. 保存图和审阅边界，当前两块不导出为正式合并表达
fig, axes = plt.subplots(1, 2, figsize=(13, 7))
positions = np.arange(len(types))
axes[0].barh(positions - 0.18, type_summary['Job23_complete_varying_genes'], height=0.36, label='Job23 / block 2', color='#94b7c4')
axes[0].barh(positions + 0.18, type_summary['Job24_complete_varying_genes'], height=0.36, label='Job24 / corrected block 1', color='#267997')
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
fig.suptitle('Job24 vs Job23: output availability and shared-gene consistency')
fig.tight_layout()
fig.savefig(outdir / 'Job24_vs_Job23_HiRes_comparison.png', dpi=170)
plt.close(fig)
valid = comparison.loc[comparison['both_varying']]
summary = {'status': 'numeric_review_complete', 'source_zip_sha256': zip_hash, 'archive_members': len(integrity), 'runtime_completed': True, 'warning': error.strip(), 'requested_genes': len(request), 'hires_loaded_genes': len(output_gene_order), 'missing_genes': missing_genes, 'unexpected_genes': unexpected_genes, 'corrected_batch01_request_matching_verified': True, 'samples': len(sample_ids), 'types': len(types), 'window': 60, 'fraction_P_value_unique': fractions24['P-value'].unique().tolist(), 'fraction_correlation_median': float(fractions24['Correlation'].median()), 'fraction_RMSE_median': float(fractions24['RMSE'].median()), 'median_fraction_reallocation_vs_Job23_percent': float(reallocation.median()), 'max_fraction_reallocation_vs_Job23_percent': float(reallocation.max()), 'adjusted_signature_max_absolute_difference': float(adjusted_difference.max().max()), 'adjusted_signature_values_changed_over_1e_8': int(adjusted_difference.gt(1e-8).sum().sum()), 'shared_gene_type_pairs': len(comparison), 'shared_gene_type_pairs_both_varying': len(valid), 'shared_both_varying_median_pearson_log2': float(valid['pearson_log2'].median()), 'shared_state_agreement_percent': float(comparison['same_state'].mean() * 100), 'shared_both_varying_median_symmetric_difference_percent': float(valid['median_symmetric_difference_percent'].median()), 'all_one_types': type_summary.loc[type_summary['Job24_all_one_genes'].eq(len(output_gene_order)), 'CellType'].tolist(), 'full_transcriptome_hires_expression_available': False, 'formal_merged_matrix_exported': False, 'ecotyper_started': False}
(outdir / 'review_summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps(summary, ensure_ascii=False, indent=2))

# Step7. 对照旧第一块：只比较实际共同的1,155基因，新增94基因单独统计
# Job21与Job24的实际加载集合不同，因此这不是严格相同输入的重复实验。
old_block = project / 'CIBERSORTx_signature_audit/MSCCR_Joint20v2_Job21_HiRes_review_20261009/original_outputs'
repeat_rows = []
repeat_gene_rows = []
for cell_type in types:
    old = pd.read_csv(old_block / f'CIBERSORTxHiRes_Job21_{cell_type}_Window60.txt', sep='\t', index_col=0)
    new = pd.read_csv(originals / f'CIBERSORTxHiRes_Job24_{cell_type}_Window60.txt', sep='\t', index_col=0)
    assert old.columns.tolist() == new.columns.tolist() == sample_ids
    common = [gene for gene in old.index if gene in new.index]
    extra = [gene for gene in new.index if gene not in old.index]
    assert len(common) == 1155 and len(extra) == 94
    extra_varying = new.loc[extra].notna().all(axis=1) & (new.loc[extra].max(axis=1) - new.loc[extra].min(axis=1)).gt(1e-8)
    old_varying = old.notna().all(axis=1) & (old.max(axis=1) - old.min(axis=1)).gt(1e-8)
    new_varying = new.notna().all(axis=1) & (new.max(axis=1) - new.min(axis=1)).gt(1e-8)
    for gene in common:
        x = old.loc[gene].to_numpy(dtype=float)
        y = new.loc[gene].to_numpy(dtype=float)
        old_state = 'all_NA' if np.isnan(x).all() else 'all_one' if (x == 1).all() else 'varying' if old_varying[gene] else 'other'
        new_state = 'all_NA' if np.isnan(y).all() else 'all_one' if (y == 1).all() else 'varying' if new_varying[gene] else 'other'
        comparable = bool(old_varying[gene] and new_varying[gene])
        r = float(np.corrcoef(np.log2(x + 1), np.log2(y + 1))[0, 1]) if comparable else np.nan
        relative = float(np.median(2 * np.abs(x - y) / np.maximum(np.abs(x) + np.abs(y), 1e-12)) * 100) if comparable else np.nan
        repeat_gene_rows.append({'CellType': cell_type, 'GeneSymbol': gene, 'state_Job21': old_state, 'state_Job24': new_state, 'same_state': old_state == new_state, 'both_varying': comparable, 'pearson_log2': r, 'median_symmetric_difference_percent': relative, 'exact_numeric_match': bool(np.array_equal(x, y, equal_nan=True)), 'Job21_mean': float(np.nanmean(x)) if not np.isnan(x).all() else np.nan, 'Job24_mean': float(np.nanmean(y)) if not np.isnan(y).all() else np.nan})
    repeat_rows.append({'CellType': cell_type, 'Job21_varying': int(old_varying.sum()), 'Job24_varying': int(new_varying.sum()), 'Job24_new94_varying': int(extra_varying.sum()), 'Job24_common1155_varying': int(new_varying.loc[common].sum())})
repeat_comparison = pd.DataFrame(repeat_gene_rows)
repeat_type_summary = pd.DataFrame(repeat_rows)
for cell_type in types:
    selected = repeat_comparison.loc[repeat_comparison['CellType'].eq(cell_type)]
    valid_repeat = selected.loc[selected['both_varying']]
    repeat_type_summary.loc[repeat_type_summary.CellType.eq(cell_type), 'Both_varying'] = len(valid_repeat)
    repeat_type_summary.loc[repeat_type_summary.CellType.eq(cell_type), 'Median_log2_r'] = valid_repeat['pearson_log2'].median()
    repeat_type_summary.loc[repeat_type_summary.CellType.eq(cell_type), 'Below_0_8_count'] = int(valid_repeat['pearson_log2'].lt(0.8).sum())
repeat_type_summary.to_csv(outdir / 'Job24_vs_Job21_same_block_type_summary.tsv', sep='\t', index=False)
repeat_comparison.to_csv(outdir / 'Job24_vs_Job21_common1155_gene_comparison.tsv', sep='\t', index=False)
pd.DataFrame({'GeneSymbol': extra}).to_csv(outdir / 'Job24_restored94_genes.tsv', sep='\t', index=False)
print('第一块重跑对照：')
print(repeat_type_summary.to_string(index=False))
valid_repeat = repeat_comparison.loc[repeat_comparison['both_varying']]
summary['Job24_vs_Job23_both_varying_below_0_8_count'] = int(valid['pearson_log2'].lt(0.8).sum())
summary['Job24_vs_Job23_own_type_both_varying_count'] = int((comparison['both_varying'] & comparison['CellType'].eq(comparison['selection_type'])).sum())
summary['Job21_shared_genes'] = len(common)
summary['restored_genes_vs_Job21'] = len(extra)
summary['Job21_both_varying_pairs'] = len(valid_repeat)
summary['Job21_both_varying_median_log2_r'] = float(valid_repeat['pearson_log2'].median())
summary['Job21_both_varying_below_0_8_count'] = int(valid_repeat['pearson_log2'].lt(0.8).sum())
(outdir / 'review_summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')

# Step8. 将关键例外的样本信号画出来，区分接近恒定与排序改变
fig, axes = plt.subplots(2, 3, figsize=(12, 7))
examples = [('Plasma', 'JCHAIN'), ('Goblet', 'ITLN1'), ('Paneth', 'DEFA5')]
for column, (cell_type, gene) in enumerate(examples):
    current = pd.read_csv(originals / f'CIBERSORTxHiRes_Job24_{cell_type}_Window60.txt', sep='\t', index_col=0).loc[gene]
    for row, (job, folder) in enumerate([('Job21', old_block), ('Job23', prior)]):
        earlier = pd.read_csv(folder / f'CIBERSORTxHiRes_{job}_{cell_type}_Window60.txt', sep='\t', index_col=0).loc[gene]
        x = np.log2(earlier.to_numpy() + 1)
        y = np.log2(current.to_numpy() + 1)
        r = np.corrcoef(x, y)[0, 1]
        axes[row, column].scatter(x, y, s=14, alpha=0.65, color='#267997')
        axes[row, column].set_title(f'{cell_type}: {gene}\nr = {r:.3f}')
        axes[row, column].set_xlabel(f'{job}: log2(expression + 1)')
        axes[row, column].set_ylabel('Job24: log2(expression + 1)')
fig.suptitle('Same 120 samples: key HiRes expression comparisons')
fig.tight_layout()
fig.savefig(outdir / 'Job24_key_gene_repeat_comparison.png', dpi=170)
plt.close(fig)

