# 审阅 MSCCR Job16：核对计算完整性、HTML/TXT一致性、各部位分布和全零类型，区分PDF导出错误与求解失败。
# 在本机 cibersortx_review 环境运行；读取用户回传4个原件、固定discovery注释和现有CPM输入，保存描述性审阅。
# 真实bulk没有比例金标准，不计算“准确率”或临床关联；不改变矩阵、删除样本、改校正模式或登录服务器。
from pathlib import Path
import hashlib
import json
import shutil
import re
import html
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Step1. 保存全部原件，与已提交的样本和类型核对
root = Path('D:/IBD_ecotype/IBD_ecotype_analysis')
ready = root / 'CIBERSORTx_ready_20261007/09_MSCCR_Joint20v2_Fractions'
out = root / 'CIBERSORTx_signature_audit/MSCCR_Joint20v2_Job16_review_20261008'
original = out / 'original_outputs'
original.mkdir(parents=True, exist_ok=True)
integrity = []
for name in ['CIBERSORTx_Job16_error_log.txt', 'CIBERSORTx_Job16_runtime_log.txt', 'CIBERSORTx_Job16_Results.txt', 'CIBERSORTx_Job16_Results.html']:
    source = Path('C:/Users/马泽凯/Downloads') / name
    destination = original / name
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    if destination.exists():
        assert hashlib.sha256(destination.read_bytes()).hexdigest() == digest, '不覆盖不同的网页原件'
    else:
        shutil.copyfile(source, destination)
    integrity.append({'source': str(source), 'archived_file': name, 'bytes': source.stat().st_size, 'sha256': digest})
pd.DataFrame(integrity).to_csv(out / 'original_integrity.tsv', sep='\t', index=False)
log = (original / 'CIBERSORTx_Job16_runtime_log.txt').read_text(encoding='utf-8-sig')
error = (original / 'CIBERSORTx_Job16_error_log.txt').read_text(encoding='utf-8-sig')
assert 'GJ20v2_Signature20.txt' in log and 'MSCCR_CD_UC_2025_GJ20_CPM.txt' in log
assert 'perm: 100' in log and 'QN: FALSE' in log and 'All done.' in log
assert 'rmbatchSmode: TRUE' not in log and 'rmbatchBmode: TRUE' not in log
annotation = pd.read_csv(ready / '核对记录_勿上传/MSCCR_CD_UC_2025_annotation.tsv', sep='\t', dtype=str)
signature = pd.read_csv(ready / 'GJ20v2_Signature20.txt', sep='\t', index_col=0)
raw = pd.read_csv(original / 'CIBERSORTx_Job16_Results.txt', sep='\t', index_col=0)
assert raw.index.is_unique and raw.columns.is_unique and set(raw.index) == set(annotation.ID)
assert raw.index.tolist() == annotation.ID.tolist()
assert set(raw.columns) == set(signature.columns) | {'P-value', 'Correlation', 'RMSE'}
fractions = raw.loc[annotation.ID, signature.columns]
assert fractions.shape == (2025, 20) and np.isfinite(raw.to_numpy()).all()
assert fractions.ge(0).all().all() and fractions.le(1).all().all() and np.allclose(fractions.sum(axis=1), 1)
joined = annotation.set_index('ID').join(raw, validate='one_to_one')
joined.to_csv(out / 'Job16_results_with_frozen_annotation.tsv', sep='\t', index_label='sample')
fractions.to_csv(out / 'Job16_fractions_aligned.tsv', sep='\t', index_label='sample')
all_zero = fractions.columns[fractions.eq(0).all()].tolist()
print('Complete 2025 × 20, identical frozen sample order; all-zero types:', all_zero, flush=True)

# Step2. 只解析网页表格的数据，不执行HTML中的脚本或外部资源
html_source = (original / 'CIBERSORTx_Job16_Results.html').read_text(encoding='utf-8-sig')
table = re.search(r'<table\b[^>]*id="table_id"[^>]*>(.*?)</table>', html_source, flags=re.S).group(1)
html_rows = []
for row in re.findall(r'<tr\b[^>]*>(.*?)</tr>', table, flags=re.S):
    fields = re.findall(r'<t[dh]\b[^>]*>(.*?)</t[dh]>', row, flags=re.S)
    html_rows.append([html.unescape(re.sub(r'<[^>]+>', '', field)).strip() for field in fields])
display = pd.DataFrame(html_rows[1:], columns=html_rows[0]).set_index('Input Sample').astype(float)
assert display.shape == (2025, 22) and display.index.tolist() == raw.index.tolist()
html_difference = (display.loc[raw.index, signature.columns] - fractions).abs().to_numpy().max()
assert html_difference <= 0.00050001
assert 'MpdfException' in error and 'pcre.backtrack_limit' in error and 'inc.htmlpdf.php' in error
print('HTML matches all 40500 fractions within display rounding; error is in HTML→PDF export', flush=True)

# Step3. 保存区域和总体分布，不按P值删样本，不开展疾病/临床关联
region_means = joined.groupby('Region3')[signature.columns].mean().T * 100
region_medians = joined.groupby('Region3')[signature.columns].median().T * 100
region_zeros = fractions.eq(0).join(joined.Region3).groupby('Region3')[signature.columns].mean().T * 100
region_means.to_csv(out / 'regional_mean_percent.tsv', sep='\t', index_label='cell_type')
region_medians.to_csv(out / 'regional_median_percent.tsv', sep='\t', index_label='cell_type')
region_zeros.to_csv(out / 'regional_zero_percent.tsv', sep='\t', index_label='cell_type')
type_summary = pd.DataFrame({'mean_percent': fractions.mean()*100, 'median_percent': fractions.median()*100, 'p95_percent': fractions.quantile(.95)*100, 'max_percent': fractions.max()*100, 'zero_samples': fractions.eq(0).sum(), 'zero_fraction': fractions.eq(0).mean(), 'all_zero': fractions.eq(0).all()})
type_summary.to_csv(out / 'overall_type_distribution.tsv', sep='\t', index_label='cell_type')
fit = raw[['P-value', 'Correlation', 'RMSE']].describe(percentiles=[.01, .05, .25, .5, .75, .95, .99])
fit.to_csv(out / 'expression_fit_distribution.tsv', sep='\t', index_label='statistic')
joined.groupby('Region3')[['Correlation', 'RMSE']].describe().to_csv(out / 'regional_expression_fit.tsv', sep='\t')
raw.sort_values('Correlation').head(30).join(annotation.set_index('ID')).to_csv(out / 'lowest30_expression_correlations.tsv', sep='\t', index_label='sample')
print('Correlation median:', raw.Correlation.median(), '; P-value counts:', raw['P-value'].value_counts().to_dict(), flush=True)
print(region_means.round(3).to_string(), flush=True)

# Step4. 有限检查代表基因及两个全零类型的signature基因覆盖，确认是否只是bulk缺基因
markers = ['IGKC', 'IGHG1', 'IGHA1', 'JCHAIN', 'MZB1', 'SDC1', 'DERL3', 'XBP1', 'IL7R', 'KIT', 'KLRB1', 'RORC', 'DEFA5', 'DEFA6', 'MUC2', 'OLFM4', 'LGR5']
marker_blocks = []
bulk_gene_nonzero = {}
for block in pd.read_csv(ready / 'MSCCR_CD_UC_2025_GJ20_CPM.txt', sep='\t', index_col=0, chunksize=500):
    hits = block.index.intersection(markers)
    if len(hits):
        marker_blocks.append(block.loc[hits])
    bulk_gene_nonzero.update(dict(zip(block.index, block.gt(0).any(axis=1))))
marker_expression = pd.concat(marker_blocks)
marker_rows = []
for gene in markers:
    if gene not in marker_expression.index:
        marker_rows.append({'gene': gene, 'in_bulk': False})
        continue
    expression = marker_expression.loc[gene]
    marker_rows.append({'gene': gene, 'in_bulk': True, 'in_signature': gene in signature.index, 'median_CPM': expression.median(), 'mean_CPM': expression.mean(), 'positive_samples': int(expression.gt(0).sum()), 'max_CPM': expression.max()})
marker_table = pd.DataFrame(marker_rows)
marker_table.to_csv(out / 'marker_bulk_expression_check.tsv', sep='\t', index=False)
marker_expression.T.join(annotation.set_index('ID')).to_csv(out / 'marker_CPM_with_annotation.tsv', sep='\t', index_label='sample')
zero_support = []
for cell_type in all_zero:
    other_max = signature.drop(columns=cell_type).max(axis=1)
    advantage = signature[cell_type].div(other_max.replace(0, np.nan))
    candidates = advantage.sort_values(ascending=False).head(20).index
    evidence = pd.DataFrame({'reference_expression': signature.loc[candidates, cell_type], 'reference_fold_over_other_max': advantage.loc[candidates], 'in_bulk': candidates.isin(bulk_gene_nonzero), 'nonzero_in_bulk': [bulk_gene_nonzero.get(gene, False) for gene in candidates]}, index=candidates)
    evidence.to_csv(out / f'{cell_type}_reference_top20_expression_advantage.tsv', sep='\t', index_label='gene')
    zero_support.append({'cell_type': cell_type, 'top20_expression_advantage_genes': len(evidence), 'in_bulk': int(evidence.in_bulk.sum()), 'nonzero_in_bulk': int(evidence.nonzero_in_bulk.sum())})
pd.DataFrame(zero_support).to_csv(out / 'all_zero_type_reference_coverage.tsv', sep='\t', index=False)
print(marker_table.to_string(index=False), flush=True)
print('All-zero type gene coverage:', zero_support, flush=True)

# 原件里的IGKC信号极低，只记录参考尺度与既有基因身份依据，不猜测映射或删marker。
plasma_markers = marker_table.loc[marker_table.gene.isin(['IGKC', 'IGHA1', 'IGHG1', 'JCHAIN', 'MZB1', 'DERL3', 'XBP1'])].copy()
plasma_markers['Plasma_signature_expression'] = plasma_markers.gene.map(signature.Plasma)
plasma_markers.to_csv(out / 'Plasma_marker_reference_bulk_scale.tsv', sep='\t', index=False)
identity_root = root / 'CIBERSORTx_signature_audit/MSCCR_web_input_adaptation_20261006/server_outputs/audit'
identity = pd.read_csv(identity_root / 'msccr_ibd_joint_web_20261006/exported_gene_identity.tsv', sep='\t')
identity.loc[identity.astype(str).eq('IGKC').any(axis=1)].to_csv(out / 'IGKC_existing_export_identity.tsv', sep='\t', index=False)

# Step5. 一张描述性图并列显示比例与零值，不将零值当组织中不存在
regions = ['Ileum', 'Colon', 'Rectum']
fig, axes = plt.subplots(1, 2, figsize=(10, 8))
for axis, frame, title, maximum, cmap in [(axes[0], region_means, 'Mean inferred fraction (%)', 75, 'viridis'), (axes[1], region_zeros, 'Samples inferred as zero (%)', 100, 'magma')]:
    image = axis.imshow(frame.loc[signature.columns, regions].to_numpy(), vmin=0, vmax=maximum, cmap=cmap, aspect='auto')
    axis.set_xticks(range(3), regions)
    axis.set_yticks(range(20), [name.replace('_', ' ') for name in signature.columns], fontsize=8)
    axis.set_title(title, fontsize=11)
    for row_number, cell_type in enumerate(signature.columns):
        for column_number, region in enumerate(regions):
            value = frame.at[cell_type, region]
            axis.text(column_number, row_number, f'{value:.1f}', ha='center', va='center', fontsize=8, color='white' if value < maximum*.55 else 'black')
    fig.colorbar(image, ax=axis, fraction=.045, pad=.04)
fig.suptitle('MSCCR Job16: complete output, but Plasma/ILC fractions are zero in every sample', fontsize=11)
fig.tight_layout()
fig.savefig(out / 'Job16_regional_distribution.png', dpi=180)
plt.close(fig)

# Step6. 保存完成状态与表达推断限制；没有虚构比例准确率或实际运行时长
summary = {'job': 'Job16', 'status': 'numeric_review_completed_visual_review_pending', 'sample_count': 2025, 'type_count': 20, 'patient_count': int(annotation.PatientID.nunique()), 'sample_order_matches_frozen_discovery': True, 'fractions_finite_nonnegative_rows_sum1': True, 'log_complete': True, 'max_iteration_warning_in_runtime_log': 'reaching max number' in log.lower(), 'error_stage': 'Website HTML-to-PDF export: mPDF PCRE backtrack limit, not fractions solver', 'HTML_fraction_count_checked': int(fractions.size), 'HTML_TXT_max_rounding_difference': float(html_difference), 'QN_disabled': True, 'permutations': 100, 'batch_correction_enable_option_seen': False, 'all_zero_types': all_zero, 'median_expression_correlation': float(raw.Correlation.median()), 'median_expression_RMSE': float(raw.RMSE.median()), 'P_value_counts': {str(k): int(v) for k, v in raw['P-value'].value_counts().items()}, 'true_bulk_proportion_ground_truth_available': False, 'actual_runtime_measured': False, 'HiRes_started': False, 'EcoTyper_started': False, 'server_accessed': False, 'decision': 'Keep complete NS baseline; do not treat all20 as validated or feed unchanged all-zero Plasma/ILC fractions to formal HiRes. Prioritize one matched cross-platform correction path, not another broad marker/signature search.'}
(out / 'review_summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
print('Saved:', out, flush=True)
