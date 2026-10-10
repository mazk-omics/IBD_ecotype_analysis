# 审阅Job21第一批HiRes：保存原件，核对样本、请求基因与每类可推断/可变表达。
# 本机cibersortx_review环境运行；读取下载ZIP与冻结输入，写审阅表及修正命名的第二批清单，不启动EcoTyper或网页任务。
from pathlib import Path
from zipfile import ZipFile
import hashlib
import json
import shutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Step1. 归档压缩包并安全解压，保留所有原始输出
project = Path('D:/IBD_ecotype/IBD_ecotype_analysis')
source_zip = Path('C:/Users/马泽凯/Downloads/CIBERSORTx_Job21_output.zip')
inputs = project / 'CIBERSORTx_ready_20261008/04_MSCCR120_HiRes_gene_batches'
outdir = project / 'CIBERSORTx_signature_audit/MSCCR_Joint20v2_Job21_HiRes_review_20261009'
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
        output_path = (originals / item.filename).resolve()
        assert output_path.is_relative_to(originals.resolve()), f'ZIP路径越界：{item.filename}'
        if item.is_dir():
            continue
        content = archive.read(item.filename)
        if output_path.exists():
            assert output_path.read_bytes() == content
        else:
            output_path.parent.mkdir(parents=True, exist_ok=True)
            output_path.write_bytes(content)
        integrity.append({'file': item.filename, 'bytes': len(content), 'sha256': hashlib.sha256(content).hexdigest()})
pd.DataFrame(integrity).to_csv(outdir / 'original_integrity.tsv', sep='\t', index=False)
runtime = (originals / 'CIBERSORTx_Job21_runtime_log.txt').read_text(encoding='utf-8')
error = (originals / 'CIBERSORTx_Job21_error_log.txt').read_text(encoding='utf-8')
readme = (originals / 'README_HiRes.txt').read_text(encoding='utf-8')
assert 'Imputing high-resolution cell type GEPs...done.' in runtime
assert 'MSCCR120_HiRes_B01_Genes1249.txt' in runtime and 'rmbatchSmode: TRUE' in runtime and 'QN: FALSE' in runtime
print('Job21日志确认第一批HiRes完成；错误/警告日志：', error.strip())

# Step2. 核对同一120例、20类与本次新估计的fraction
annotation = pd.read_csv(inputs / '核对记录_勿上传/MSCCR120_annotation.tsv', sep='\t')
mixture = pd.read_csv(inputs / 'MSCCR_Platform120_CPM.txt', sep='\t', index_col=0)
signature = pd.read_csv(inputs / 'GJ20v2_Signature20.txt', sep='\t', index_col=0)
types = signature.columns.tolist()
fractions = pd.read_csv(originals / 'CIBERSORTxGEP_Job21_Fractions-Adjusted.txt', sep='\t', index_col=0)
assert fractions.index.tolist() == annotation['ID'].tolist() == mixture.columns.tolist()
assert fractions.columns[:20].tolist() == types
assert np.isfinite(fractions[types].to_numpy()).all() and fractions[types].ge(0).all().all()
assert np.allclose(fractions[types].sum(axis=1), 1, atol=1e-8)
fraction_records = fractions.reset_index().rename(columns={'Mixture': 'ID'}).merge(annotation, on='ID', how='left', validate='one_to_one')
fraction_records.to_csv(outdir / 'Job21_fraction_with_annotation.tsv', sep='\t', index=False)
job20 = pd.read_csv(project / 'CIBERSORTx_signature_audit/MSCCR_Joint20v2_Job20_review_20261009/original_outputs/CIBERSORTx_Job20_Adjusted.txt', sep='\t', index_col=0)
assert fractions.index.tolist() == job20.index.tolist()
fraction_change = fractions[types] - job20[types]
fraction_type_comparison = pd.DataFrame({'Job20_mean_percent': job20[types].mean() * 100, 'Job21_mean_percent': fractions[types].mean() * 100, 'change_percentage_points': fraction_change.mean() * 100})
fraction_type_comparison.to_csv(outdir / 'Job21_vs_Job20_fraction_type_summary.tsv', sep='\t', index_label='CellType')
fraction_change.mul(100).to_csv(outdir / 'Job21_vs_Job20_fraction_change_percentage_points.tsv', sep='\t', index_label='sample')
reallocation = fraction_change.abs().sum(axis=1) * 50
print('Fraction相关性/RMSE中位数：', fractions[['Correlation', 'RMSE']].median().to_dict())
print('Job21 P-value字段取值：', fractions['P-value'].unique().tolist())

# Step3. 检查样本级表达，区分全1、全NA、稳定非1及可变基因
# ZIP内README将1与NA解释为表达证据或统计功效不足；不补零、不把全1行当有效表达。
request = pd.read_csv(inputs / '先提交这两批/MSCCR120_HiRes_B01_Genes1249.txt')['GeneSymbol'].tolist()
anchors = pd.read_csv(inputs / '核对记录_勿上传/共同100基因_选择依据.tsv', sep='\t')
group_filtered = pd.read_csv(originals / 'CIBERSORTxGEP_Job21_GEPs_Filtered.txt', sep='\t', index_col=0)
group_unfiltered = pd.read_csv(originals / 'CIBERSORTxGEP_Job21_GEPs.txt', sep='\t', index_col=0)
assert group_filtered.shape == (17040, 20) and group_filtered.columns.tolist() == types
type_rows = []
gene_rows = []
anchor_rows = []
output_gene_order = None
for cell_type in types:
    expression_path = originals / f'CIBERSORTxHiRes_Job21_{cell_type}_Window60.txt'
    expression = pd.read_csv(expression_path, sep='\t', index_col=0)
    assert expression.columns.tolist() == annotation['ID'].tolist() and expression.index.is_unique
    if output_gene_order is None:
        output_gene_order = expression.index.tolist()
    assert expression.index.tolist() == output_gene_order
    values = expression.to_numpy(dtype=float)
    assert not np.isinf(values).any() and np.nanmin(values) >= 0
    all_one = expression.eq(1).all(axis=1)
    all_na = expression.isna().all(axis=1)
    complete = expression.notna().all(axis=1)
    row_mean = expression.mean(axis=1)
    row_span = expression.max(axis=1) - expression.min(axis=1)
    varying = complete & row_span.gt(1e-8)  # 仅排除数值常数，不代表生物学显著或适合NMF。
    stable_non_one = complete & ~all_one & ~varying
    type_rows.append({'CellType': cell_type, 'output_genes': len(expression), 'all_one_genes': int(all_one.sum()), 'all_NA_genes': int(all_na.sum()), 'partial_NA_genes': int((~complete & ~all_na).sum()), 'stable_non_one_genes': int(stable_non_one.sum()), 'complete_varying_genes': int(varying.sum()), 'own_five_anchors_varying': int(varying.reindex(anchors.loc[anchors['selection_type'].eq(cell_type), 'GeneSymbol']).sum()), 'group_filtered_positive_genes_all17040': int(group_filtered[cell_type].gt(0).sum()), 'group_filtered_positive_genes_in_this_request': int(group_filtered.loc[expression.index, cell_type].gt(0).sum())})
    for gene in expression.index:
        gene_rows.append({'CellType': cell_type, 'GeneSymbol': gene, 'all_one': bool(all_one[gene]), 'all_NA': bool(all_na[gene]), 'complete': bool(complete[gene]), 'varying': bool(varying[gene]), 'mean': row_mean[gene], 'span': row_span[gene]})
    for row in anchors.itertuples(index=False):
        gene = row.GeneSymbol
        if gene in expression.index:
            anchor_rows.append({'CellType': cell_type, 'GeneSymbol': gene, 'selection_type': row.selection_type, 'all_one': bool(all_one[gene]), 'all_NA': bool(all_na[gene]), 'varying': bool(varying[gene]), 'mean': row_mean[gene], 'span': row_span[gene]})
type_summary = pd.DataFrame(type_rows)
type_summary.to_csv(outdir / 'HiRes_type_usability_summary.tsv', sep='\t', index=False)
pd.DataFrame(gene_rows).to_csv(outdir / 'HiRes_gene_status.tsv', sep='\t', index=False)
pd.DataFrame(anchor_rows).to_csv(outdir / 'shared100_gene_status_by_type.tsv', sep='\t', index=False)
print(type_summary[['CellType', 'all_one_genes', 'all_NA_genes', 'complete_varying_genes']].to_string(index=False))

# Step4. 追踪请求中缺失的基因与网页命名变化，保存可核对的对应关系
missing_genes = [gene for gene in request if gene not in output_gene_order]
missing_table = pd.DataFrame({'original_request_gene': missing_genes})
missing_table['web_group_gene'] = missing_table['original_request_gene'].str.replace('-', '.', regex=False)
missing_table['in_original_mixture'] = missing_table['original_request_gene'].isin(mixture.index)
missing_table['in_group_unfiltered'] = missing_table['web_group_gene'].isin(group_unfiltered.index)
assert len(missing_genes) == 94 and all('-' in gene for gene in missing_genes)
assert missing_table['in_original_mixture'].all() and missing_table['in_group_unfiltered'].all()
assert len(output_gene_order) == 1155 and set(anchors['GeneSymbol']).issubset(output_gene_order)
missing_table.to_csv(outdir / 'missing94_request_gene_name_mapping.tsv', sep='\t', index=False)

# 第二批只改请求清单的网页命名，保留原文件及三个冻结矩阵。
# 替换规则由本次94个缺失基因的Group输出支持；实际能否全部进入HiRes仍需下一次输出确认。
request2_path = inputs / '先提交这两批/MSCCR120_HiRes_B02_Genes1249.txt'
request2 = pd.read_csv(request2_path)['GeneSymbol'].tolist()
webnames_all = mixture.index.str.replace('-', '.', regex=False)
assert webnames_all.is_unique and set(webnames_all) == set(group_unfiltered.index)
request2_web = [gene.replace('-', '.') for gene in request2]
assert len(request2_web) == 1249 and len(set(request2_web)) == 1249
assert set(request2_web).issubset(group_unfiltered.index)
request2_fixed_path = inputs / '先提交这两批/MSCCR120_HiRes_B02_WebNames_Genes1249.txt'
request2_text = 'GeneSymbol\n' + '\n'.join(request2_web) + '\n'
if request2_fixed_path.exists():
    assert request2_fixed_path.read_text(encoding='utf-8') == request2_text
else:
    request2_fixed_path.write_text(request2_text, encoding='utf-8')
pd.DataFrame({'original_gene': request2, 'web_request_gene': request2_web}).to_csv(outdir / 'Batch02_fixed_gene_name_mapping.tsv', sep='\t', index=False)

# Step5. 保存审阅状态和简洁图件，不据第一批合并或启动discovery
fig, ax = plt.subplots(figsize=(9, 6))
ax.barh(type_summary['CellType'], type_summary['complete_varying_genes'], color='#247b9a')
ax.invert_yaxis()
ax.set_xlabel('Complete nonconstant genes in Job21 (1155 requested genes loaded)')
ax.set_title('Job21: nonconstant sample-level expression by cell type')
fig.tight_layout()
fig.savefig(outdir / 'HiRes_variable_gene_counts.png', dpi=160)
plt.close(fig)
summary = {'status': 'numeric_review_complete', 'source_zip_sha256': zip_hash, 'archive_members': len(integrity), 'runtime_completed': True, 'warning': error.strip(), 'requested_genes': len(request), 'hires_loaded_genes': len(output_gene_order), 'missing_genes': len(missing_genes), 'missing_gene_hyphen_to_dot_mapping_supported': True, 'corrected_request_web_parsing_verified': False, 'samples': 120, 'types': 20, 'window': 60, 'fraction_P_value_unique': fractions['P-value'].unique().tolist(), 'fraction_correlation_median': float(fractions['Correlation'].median()), 'fraction_RMSE_median': float(fractions['RMSE'].median()), 'median_fraction_reallocation_vs_Job20_percent': float(reallocation.median()), 'all_one_types': type_summary.loc[type_summary['all_one_genes'].eq(1155), 'CellType'].tolist(), 'full_transcriptome_hires_expression_available': False, 'cross_batch_consistency_verified': False, 'ecotyper_started': False, 'corrected_batch02_request': str(request2_fixed_path)}
(outdir / 'review_summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
print('审阅目录：', outdir)
print('命名修正后的第二批清单：', request2_fixed_path)
