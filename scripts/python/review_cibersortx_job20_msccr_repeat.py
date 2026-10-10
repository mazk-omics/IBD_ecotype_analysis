# 审阅Job20的120例比例重跑，与Job17同样本结果比较；不进行临床关联或HiRes计算。
# 本机cibersortx_review环境运行，读取Downloads原件和冻结注释，保存独立审计表与PDF预览。
# 配对差异是重复运行的差异，不是相对真实细胞比例的误差；参数确认采用用户提供的日志摘录。
from pathlib import Path
import hashlib
import json
import shutil
import numpy as np
import pandas as pd
from pypdf import PdfReader
import pypdfium2 as pdfium

# Step1. 留存原件，核对样本与类型
root = Path('D:/IBD_ecotype/IBD_ecotype_analysis')
out = root / 'CIBERSORTx_signature_audit/MSCCR_Joint20v2_Job20_review_20261009'
original = out / 'original_outputs'
qa = out / 'QA'
original.mkdir(parents=True, exist_ok=True)
qa.mkdir(exist_ok=True)
source_dir = Path('C:/Users/马泽凯/Downloads')
integrity = []
for name in ['CIBERSORTx_Job20_Adjusted.txt', 'CIBERSORTx_Job20_Results.pdf', 'CIBERSORTx_Job20_stackedChart.pdf']:
    source = source_dir / name
    destination = original / name
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    if destination.exists():
        assert hashlib.sha256(destination.read_bytes()).hexdigest() == digest, '不覆盖不同原件'
    else:
        shutil.copy2(source, destination)
    integrity.append({'source': str(source), 'archived_file': name, 'bytes': source.stat().st_size, 'sha256': digest})
pd.DataFrame(integrity).to_csv(out / 'original_integrity.tsv', sep='\t', index=False)
ready = root / 'CIBERSORTx_ready_20261008/01_MSCCR_platform120_pilot'
annotation = pd.read_csv(ready / '核对记录_勿上传/selected_annotation.tsv', sep='\t', dtype=str).set_index('ID')
signature_types = pd.read_csv(ready / 'GJ20v2_Signature20.txt', sep='\t', nrows=0).columns[1:].tolist()
job17 = pd.read_csv(root / 'CIBERSORTx_signature_audit/MSCCR_Joint20v2_Job17_review_20261008/original_outputs/CIBERSORTx_Job17_Adjusted.txt', sep='\t', index_col=0)
job20 = pd.read_csv(original / 'CIBERSORTx_Job20_Adjusted.txt', sep='\t', index_col=0)
assert job20.index.is_unique and job20.columns.is_unique
assert job20.index.tolist() == job17.index.tolist() == annotation.index.tolist()
assert job20.columns.tolist() == job17.columns.tolist() == signature_types + ['P-value', 'Correlation', 'RMSE']
before = job17[signature_types]
after = job20[signature_types]
assert after.shape == (120, 20) and annotation.PatientID.nunique() == 120
assert np.isfinite(job20.to_numpy()).all() and after.ge(0).all().all() and after.le(1).all().all()
assert np.allclose(after.sum(axis=1), 1, atol=1e-8)
assert job20['P-value'].between(0, 1).all() and job20.Correlation.between(-1, 1).all() and job20.RMSE.ge(0).all()
annotation.join(job20, validate='one_to_one').to_csv(out / 'Job20_results_with_frozen_annotation.tsv', sep='\t', index_label='sample')
print('Job20：', job20.shape, '有效比例：', after.shape, '最大行和误差：', abs(after.sum(axis=1)-1).max())

# Step2. 核对PDF显示值，并渲染两份图件
reader = PdfReader(original / 'CIBERSORTx_Job20_Results.pdf')
display_text = '\n'.join(page.extract_text() or '' for page in reader.pages)
display_rows = [line.split() for line in display_text.splitlines() if line.startswith('MSCCR_reGRID_')]
assert len(display_rows) == 120 and all(len(row) == 23 for row in display_rows)
display = pd.DataFrame(display_rows, columns=['sample'] + signature_types + ['P-value', 'Correlation']).set_index('sample').astype(float)
assert display.index.tolist() == job20.index.tolist()
pdf_difference = (display - job20[signature_types + ['P-value', 'Correlation']]).abs().to_numpy().max()
assert pdf_difference <= 0.00050001
(qa / 'Job20_Results_text.txt').write_text(display_text, encoding='utf-8')
for name in ['CIBERSORTx_Job20_Results.pdf', 'CIBERSORTx_Job20_stackedChart.pdf']:
    document = pdfium.PdfDocument(original / name)
    assert len(document) == 1, '多页图件需要逐页审阅'
    page = document[0]
    page.render(scale=min(2, 2000 / page.get_width())).to_pil().save(qa / (Path(name).stem + '.png'))
    document.close()
print('PDF与TXT最大显示差异：', pdf_difference)

# Step3. 配对比较重复运行的拟合和比例变化
changes = (after - before) * 100
changes.to_csv(out / 'Job20_minus_Job17_percentage_points.tsv', sep='\t', index_label='sample')
paired = annotation.copy()
paired['Job17_correlation'] = job17.Correlation
paired['Job20_correlation'] = job20.Correlation
paired['correlation_change'] = job20.Correlation - job17.Correlation
paired['Job17_RMSE'] = job17.RMSE
paired['Job20_RMSE'] = job20.RMSE
paired['RMSE_change'] = job20.RMSE - job17.RMSE
paired['fraction_reallocation_L1_half'] = (after-before).abs().sum(axis=1)/2
paired['largest_type_change_percentage_points'] = changes.abs().max(axis=1)
paired['largest_changed_type'] = changes.abs().idxmax(axis=1)
paired.to_csv(out / 'paired_sample_metrics.tsv', sep='\t', index_label='sample')
type_summary = pd.DataFrame({'Job17_mean_percent': before.mean()*100, 'Job20_mean_percent': after.mean()*100, 'mean_change_percentage_points': changes.mean(), 'median_abs_change_percentage_points': changes.abs().median(), 'p95_abs_change_percentage_points': changes.abs().quantile(.95), 'max_abs_change_percentage_points': changes.abs().max(), 'Job17_zero_samples': before.eq(0).sum(), 'Job20_zero_samples': after.eq(0).sum(), 'Job20_median_percent': after.median()*100, 'Job20_max_percent': after.max()*100})
type_summary['paired_spearman'] = [before[cell_type].corr(after[cell_type], method='spearman') for cell_type in signature_types]
type_summary.to_csv(out / 'paired_type_distribution.tsv', sep='\t', index_label='cell_type')
regions = ['Ileum', 'Colon', 'Rectum']
region_rows = []
for region in regions:
    mask = annotation.Region3.eq(region)
    region_rows.append({'region': region, 'samples': int(mask.sum()), 'Job17_correlation_median': job17.loc[mask, 'Correlation'].median(), 'Job20_correlation_median': job20.loc[mask, 'Correlation'].median(), 'Job17_RMSE_median': job17.loc[mask, 'RMSE'].median(), 'Job20_RMSE_median': job20.loc[mask, 'RMSE'].median(), 'median_reallocation_fraction': paired.loc[mask, 'fraction_reallocation_L1_half'].median(), 'Job20_Plasma_zero_samples': int(after.loc[mask, 'Plasma'].eq(0).sum()), 'Job20_ILC_zero_samples': int(after.loc[mask, 'ILC'].eq(0).sum())})
region_summary = pd.DataFrame(region_rows)
region_summary.to_csv(out / 'paired_regional_fit.tsv', sep='\t', index=False)
regional20 = after.join(annotation.Region3).groupby('Region3')[signature_types].mean().T[regions]*100
regional17 = before.join(annotation.Region3).groupby('Region3')[signature_types].mean().T[regions]*100
regional20.to_csv(out / 'regional_Job20_mean_percent.tsv', sep='\t', index_label='cell_type')
(regional20-regional17).to_csv(out / 'regional_mean_change_percentage_points.tsv', sep='\t', index_label='cell_type')

# Step4. 核对用户日志摘录中的设置，保留随机机制与实际网页输入哈希的未知边界
log_path = original / 'Job20_runtime_excerpt_user_provided.txt'
settings_verified = False
if log_path.exists():
    log_text = log_path.read_text(encoding='utf-8')
    expected_settings = ['[Options] perm: 100', '[Options] rmbatchSmode: TRUE', '[Options] QN: FALSE', 'Run mode (relative or absolute): relative', 'GJ20v2_Signature20.txt', 'MSCCR_Platform120_CPM.txt', 'GJ20v2_SourceSubset1060_scRNA_CPM.txt', 'All done.']
    assert all(value in log_text for value in expected_settings)
    settings_verified = True
summary = {'status': 'numeric_review_complete_visual_review_pending', 'output_kind': '120-sample relative cell fractions, not HiRes gene expression', 'samples': 120, 'types': 20, 'distinct_patients': 120, 'same_frozen_sample_order_as_Job17': True, 'fractions_finite_nonnegative_and_sum_to_one': True, 'fraction_row_sum_max_abs_error': float(abs(after.sum(axis=1)-1).max()), 'pdf_max_abs_display_rounding_difference': float(pdf_difference), 'Job17_correlation_median': float(job17.Correlation.median()), 'Job20_correlation_median': float(job20.Correlation.median()), 'Job17_RMSE_median': float(job17.RMSE.median()), 'Job20_RMSE_median': float(job20.RMSE.median()), 'Job20_P_below_0_05_samples': int(job20['P-value'].lt(.05).sum()), 'Job20_all_zero_types': after.columns[after.eq(0).all()].tolist(), 'Job20_Plasma_positive_samples': int(after.Plasma.gt(0).sum()), 'Job20_ILC_positive_samples': int(after.ILC.gt(0).sum()), 'median_reallocation_percent': float(paired.fraction_reallocation_L1_half.median()*100), 'p95_reallocation_percent': float(paired.fraction_reallocation_L1_half.quantile(.95)*100), 'maximum_absolute_type_change_percentage_points': float(changes.abs().to_numpy().max()), 'minimum_type_paired_spearman': float(type_summary.paired_spearman.min()), 'runtime_log_provided': settings_verified, 'runtime_log_scope': 'user-provided text excerpt', 'settings_verified_from_log': settings_verified, 'web_input_file_hashes_verified': False, 'cause_of_repeat_differences_confirmed': False, 'HiRes_expression_provided': False, 'true_fraction_gold_standard_available': False, 'clinical_association_performed': False}
(out / 'review_summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps(summary, ensure_ascii=False, indent=2))
print(type_summary.round(4).to_string())
print(region_summary.round(4).to_string(index=False))
print(regional20.round(3).to_string())
print('审阅目录：', out)
