# 审阅Job17真实MSCCR的120例S-mode诊断，与Job16同样本未校正结果逐例比较；不做临床关联。
# 本机cibersortx_review环境运行，读取用户3份结果、冻结120例注释/CPM和原v2 signature，保存对照表、图和来源记录。
# 没有真实比例金标准；表达拟合、marker相关和非零比例均不能单独作为准确率，来源子集配套限制仍保留。
from pathlib import Path
import hashlib
import json
import shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pypdf import PdfReader
import pypdfium2 as pdfium

# Step1. 保存原件并核对冻结样本、类型和数值
root = Path('D:/IBD_ecotype/IBD_ecotype_analysis')
ready = root / 'CIBERSORTx_ready_20261008/01_MSCCR_platform120_pilot'
out = root / 'CIBERSORTx_signature_audit/MSCCR_Joint20v2_Job17_review_20261008'
original = out / 'original_outputs'
qa = out / 'QA'
original.mkdir(parents=True, exist_ok=True)
qa.mkdir(exist_ok=True)
checks = []
for name in ['CIBERSORTx_Job17_Adjusted.txt', 'CIBERSORTx_Job17_Results.pdf', 'CIBERSORTx_Job17_stackedChart.pdf']:
    source = Path('C:/Users/马泽凯/Downloads') / name
    destination = original / name
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    if destination.exists():
        assert hashlib.sha256(destination.read_bytes()).hexdigest() == digest, '不覆盖不同原件'
    else:
        shutil.copyfile(source, destination)
    checks.append({'source': str(source), 'archived_file': name, 'bytes': source.stat().st_size, 'sha256': digest})
pd.DataFrame(checks).to_csv(out / 'original_integrity.tsv', sep='\t', index=False)
annotation = pd.read_csv(ready / '核对记录_勿上传/selected_annotation.tsv', sep='\t', dtype=str).set_index('ID')
signature = pd.read_csv(ready / 'GJ20v2_Signature20.txt', sep='\t', index_col=0)
baseline = pd.read_csv(ready / '核对记录_勿上传/Job16_NS_same120_results.tsv', sep='\t', index_col=0)
adjusted = pd.read_csv(original / 'CIBERSORTx_Job17_Adjusted.txt', sep='\t', index_col=0)
assert adjusted.index.is_unique and adjusted.columns.is_unique
assert adjusted.index.tolist() == annotation.index.tolist() == baseline.index.tolist()
assert adjusted.columns.tolist() == baseline.columns.tolist() == signature.columns.tolist() + ['P-value', 'Correlation', 'RMSE']
types = signature.columns.tolist()
before = baseline[types]
after = adjusted[types]
assert after.shape == before.shape == (120, 20) and np.isfinite(adjusted.to_numpy()).all()
assert after.ge(0).all().all() and after.le(1).all().all() and np.allclose(after.sum(axis=1), 1)
assert annotation.PatientID.nunique() == 120 and annotation.Region3.value_counts().eq(40).all()
annotation.join(adjusted, validate='one_to_one').to_csv(out / 'Job17_results_with_frozen_annotation.tsv', sep='\t', index_label='sample')
print('Aligned 120 distinct donors × 20 types; fraction row sums:', after.sum(axis=1).min(), after.sum(axis=1).max(), flush=True)

# Step2. 核对PDF显示值与TXT，渲染原图供人工检查；PDF以三位小数显示而TXT保留更多精度
display_text = PdfReader(original / 'CIBERSORTx_Job17_Results.pdf').pages[0].extract_text()
display_rows = [line.split() for line in display_text.splitlines() if line.startswith('MSCCR_reGRID_')]
assert len(display_rows) == 120 and all(len(row) == 23 for row in display_rows)
display = pd.DataFrame(display_rows, columns=['sample'] + types + ['P-value', 'Correlation']).set_index('sample').astype(float)
assert display.index.tolist() == adjusted.index.tolist()
pdf_difference = (display - adjusted[types + ['P-value', 'Correlation']]).abs().to_numpy().max()
assert pdf_difference <= 0.00050001
for name in ['CIBERSORTx_Job17_Results.pdf', 'CIBERSORTx_Job17_stackedChart.pdf']:
    document = pdfium.PdfDocument(original / name)
    page = document[0]
    image = page.render(scale=min(2, 2000 / page.get_width())).to_pil()
    image.save(qa / (Path(name).stem + '.png'))
    document.close()
(qa / 'CIBERSORTx_Job17_Results_text.txt').write_text(display_text, encoding='utf-8')
print('PDF/TXT maximum display-rounding difference:', pdf_difference, flush=True)

# Step3. 逐样本保存拟合变化与比例重分配；改变量不是相对真值的误差
paired = annotation.copy()
paired['NS_expression_correlation'] = baseline.Correlation
paired['S_expression_correlation'] = adjusted.Correlation
paired['correlation_change'] = adjusted.Correlation - baseline.Correlation
paired['NS_RMSE'] = baseline.RMSE
paired['S_RMSE'] = adjusted.RMSE
paired['RMSE_change'] = adjusted.RMSE - baseline.RMSE
paired['fraction_reallocation_L1_half'] = (after - before).abs().sum(axis=1) / 2
paired.to_csv(out / 'paired_sample_metrics.tsv', sep='\t', index_label='sample')
changes = (after - before) * 100
changes.to_csv(out / 'paired_fraction_change_percentage_points.tsv', sep='\t', index_label='sample')
type_summary = pd.DataFrame({'NS_mean_percent': before.mean()*100, 'S_mean_percent': after.mean()*100, 'S_median_percent': after.median()*100, 'S_p95_percent': after.quantile(.95)*100, 'S_max_percent': after.max()*100, 'NS_zero_samples': before.eq(0).sum(), 'S_zero_samples': after.eq(0).sum(), 'S_below_0_1_percent_samples': after.lt(.001).sum(), 'mean_change_percentage_points': changes.mean()})
type_summary.to_csv(out / 'paired_type_distribution.tsv', sep='\t', index_label='cell_type')
regions = ['Ileum', 'Colon', 'Rectum']
region_rows = []
for region in regions:  # 三个部位各40例，保留配对比较，不汇总临床结局
    mask = annotation.Region3.eq(region)
    region_rows.append({'region': region, 'samples': int(mask.sum()), 'NS_correlation_median': baseline.loc[mask, 'Correlation'].median(), 'S_correlation_median': adjusted.loc[mask, 'Correlation'].median(), 'NS_RMSE_median': baseline.loc[mask, 'RMSE'].median(), 'S_RMSE_median': adjusted.loc[mask, 'RMSE'].median(), 'correlation_improved_samples': int(paired.loc[mask, 'correlation_change'].gt(0).sum()), 'RMSE_improved_samples': int(paired.loc[mask, 'RMSE_change'].lt(0).sum()), 'S_Plasma_zero_samples': int(after.loc[mask, 'Plasma'].eq(0).sum()), 'S_ILC_zero_samples': int(after.loc[mask, 'ILC'].eq(0).sum())})
region_fit = pd.DataFrame(region_rows)
region_fit.to_csv(out / 'paired_regional_fit.tsv', sep='\t', index=False)
regional_before = before.join(annotation.Region3).groupby('Region3')[types].mean().T[regions] * 100
regional_after = after.join(annotation.Region3).groupby('Region3')[types].mean().T[regions] * 100
regional_zero = after.eq(0).join(annotation.Region3).groupby('Region3')[types].mean().T[regions] * 100
regional_before.to_csv(out / 'regional_NS_mean_percent.tsv', sep='\t', index_label='cell_type')
regional_after.to_csv(out / 'regional_S_mean_percent.tsv', sep='\t', index_label='cell_type')
regional_zero.to_csv(out / 'regional_S_zero_percent.tsv', sep='\t', index_label='cell_type')
print(region_fit.round(4).to_string(index=False), flush=True)
print(type_summary.round(3).to_string(), flush=True)

# Step4. 同样本代表marker的描述性关联，只看支持与冲突，不把相关系数当准确率
bulk = pd.read_csv(ready / 'MSCCR_Platform120_CPM.txt', sep='\t', index_col=0)
assert bulk.columns.tolist() == after.index.tolist()
marker_pairs = {'Plasma': ['JCHAIN', 'MZB1', 'DERL3', 'XBP1', 'IGHA1', 'IGKC'], 'ILC': ['IL7R', 'KLRB1', 'RORC', 'KIT'], 'Paneth': ['DEFA5', 'DEFA6', 'LYZ', 'REG3A'], 'Glial': ['SOX10', 'S100B', 'PLP1'], 'Tuft': ['POU2F3', 'TRPM5']}
marker_rows = []
for cell_type, markers in marker_pairs.items():
    for gene in markers:
        if gene not in bulk.index:
            marker_rows.append({'cell_type': cell_type, 'gene': gene, 'region': 'All120', 'in_bulk': False})
            continue
        expression = bulk.loc[gene]
        for region in ['All120'] + regions:
            mask = pd.Series(True, index=after.index) if region == 'All120' else annotation.Region3.eq(region)
            gene_values = expression.loc[mask]
            old_values = before.loc[mask, cell_type]
            new_values = after.loc[mask, cell_type]
            old_r = old_values.corr(gene_values, method='spearman') if old_values.nunique() > 1 and gene_values.nunique() > 1 else np.nan
            new_r = new_values.corr(gene_values, method='spearman') if new_values.nunique() > 1 and gene_values.nunique() > 1 else np.nan
            marker_rows.append({'cell_type': cell_type, 'gene': gene, 'region': region, 'in_bulk': True, 'gene_median_CPM': gene_values.median(), 'gene_nonzero_samples': int(gene_values.gt(0).sum()), 'NS_spearman': old_r, 'S_spearman': new_r})
marker_support = pd.DataFrame(marker_rows)
marker_support.to_csv(out / 'marker_fraction_descriptive_correlations.tsv', sep='\t', index=False, na_rep='NA')
print('Representative marker associations (All120; not accuracy):', flush=True)
print(marker_support.loc[marker_support.region.eq('All120')].round(3).to_string(index=False), flush=True)

# Step5. 绘制同尺度部位均值和拟合配对图，避免原堆叠图缺少部位注释
plt.rcParams.update({'font.size': 10})
fig, axes = plt.subplots(2, 2, figsize=(12, 13), gridspec_kw={'height_ratios': [2.4, 1]}, constrained_layout=True)
color_max = np.ceil(max(regional_before.to_numpy().max(), regional_after.to_numpy().max()) / 10) * 10
for ax, values, title in zip(axes[0], [regional_before, regional_after], ['Job16: no correction, same 120 samples', 'Job17: source-subset S-mode diagnostic']):
    heat = ax.imshow(values, vmin=0, vmax=color_max, cmap='YlOrRd', aspect='auto')
    ax.set_xticks(range(3), regions)
    ax.set_yticks(range(20), types)
    ax.set_title(title)
    for row in range(20):
        for col in range(3):
            value = values.iloc[row, col]
            ax.text(col, row, f'{value:.1f}', ha='center', va='center', color='white' if value > color_max*.55 else 'black', fontsize=9)
fig.colorbar(heat, ax=list(axes[0]), shrink=.6, label='Mean predicted relative fraction (%)')
colors = {'Ileum': '#c97a22', 'Colon': '#3465a4', 'Rectum': '#429066'}
for region in regions:
    mask = annotation.Region3.eq(region)
    axes[1, 0].scatter(baseline.loc[mask, 'Correlation'], adjusted.loc[mask, 'Correlation'], s=22, alpha=.75, label=region, color=colors[region])
    axes[1, 1].scatter(baseline.loc[mask, 'RMSE'], adjusted.loc[mask, 'RMSE'], s=22, alpha=.75, color=colors[region])
for ax in axes[1]:
    ax.plot([0, 1.1], [0, 1.1], '--', color='gray', lw=1)
    ax.set_xlim(0, 1.05)
    ax.set_ylim(0, 1.05)
    ax.set_xlabel('Job16, no correction')
    ax.set_ylabel('Job17, S-mode diagnostic')
axes[1, 0].set_title('Expression correlation: above line = higher fit')
axes[1, 0].legend(fontsize=9)
axes[1, 1].set_title('RMSE: below line = lower residual')
fig.savefig(out / 'Job17_vs_Job16_same120.png', dpi=180)
plt.close(fig)

# Step6. 保存诊断结论边界；设置由用户确认，未收到runtime日志或校正后的signature表达文件
summary = {'status': 'descriptive_review_complete', 'samples': 120, 'types': 20, 'distinct_patients': 120, 'exact_frozen_order': True, 'pdf_values_match_txt_with_display_rounding': True, 'pdf_maximum_absolute_rounding_difference': float(pdf_difference), 'settings_confirmed_by_user': {'signature': 'GJ20v2_Signature20.txt', 'mixture': 'MSCCR_Platform120_CPM.txt', 'source': 'GJ20v2_SourceSubset1060_scRNA_CPM.txt', 'S_mode': True, 'QN': False, 'relative': True, 'permutations': 100}, 'runtime_log_provided': False, 'adjusted_signature_expression_provided': False, 'NS_correlation_median': float(baseline.Correlation.median()), 'S_correlation_median': float(adjusted.Correlation.median()), 'NS_RMSE_median': float(baseline.RMSE.median()), 'S_RMSE_median': float(adjusted.RMSE.median()), 'samples_with_higher_correlation': int(paired.correlation_change.gt(0).sum()), 'samples_with_lower_RMSE': int(paired.RMSE_change.lt(0).sum()), 'NS_all_zero_types': before.columns[before.eq(0).all()].tolist(), 'S_all_zero_types': after.columns[after.eq(0).all()].tolist(), 'S_Plasma_positive_samples': int(after.Plasma.gt(0).sum()), 'S_ILC_positive_samples': int(after.ILC.gt(0).sum()), 'median_fraction_reallocation_L1_half': float(paired.fraction_reallocation_L1_half.median()), 'true_fraction_gold_standard_available': False, 'exact_author_protocol_reference_pairing_verified': False, 'formal_full_cohort_correction_accepted': False, 'HiRes_started': False, 'EcoTyper_started': False, 'clinical_association_performed': False, 'interpretation': 'Strong improvement in expression fit and resolution of whole-column zeros. Preserve as source-subset diagnostic; per-type accuracy and exact standard reference pairing remain unverified.'}
summary['full_cohort_candidate_S_mode_test_supported'] = True
summary['unconditional_acceptance_of_all_20_types'] = False
(out / 'review_summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
print('REVIEW_SAVED:', out, flush=True)
