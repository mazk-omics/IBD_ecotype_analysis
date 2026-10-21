# 为对话交接说明填入已完成的数值表，不重新拟合模型或改变评价终点。
# 本机科学Python运行；读取已保存的配对结果与参考容量，输出完整可复制Markdown。
from pathlib import Path
import pandas as pd

# Step1. 读入本轮新增的交接说明与已经核验的结果
project = Path(__file__).resolve().parents[2]
root = project / '06_validation/deconv_fraction_sensitivity_joint20_20261013'
path = root / '07_reports/CHATGPT_DISCUSSION_BRIEF_20261014.md'
text = path.read_text(encoding='utf-8')
summary = pd.read_csv(root / '07_reports/EXTERNAL_OFFICIAL_PAIRED_SUMMARY_20261014.tsv', sep='\t')
labels = {'Retained_gene_RNA_PRIMARY': 'BP-retained RNA', 'Cell_count_SECONDARY': 'Cell-count', 'Full_gene_RNA_SECONDARY': 'Full-gene RNA'}

# Step2. 汇总表明确保留initial/final和三种truth，避免只选有利结果
lines = ['| Dataset | Truth | Method | MAE pp | RMSE pp | Pooled r | Pooled rho |', '|---|---|---|---:|---:|---:|---:|']
for row in summary.itertuples():
    lines.append(f'| {row.dataset} | {labels[row.truth_scale]} | {row.method} | {row.MAE_pp:.3f} | {row.RMSE_pp:.3f} | {row.Pearson_r:.3f} | {row.Spearman_rho:.3f} |')
text = text.replace('{{EXTERNAL_METRICS}}', '\n'.join(lines))
capacity = pd.read_csv(root / '01_reference/reference_type_capacity.tsv', sep='\t')
lines = ['| Cell type | Cells | Donors | Profiles |', '|---|---:|---:|---:|']
for row in capacity.itertuples():
    lines.append(f'| {row.cell_type} | {row.cells} | {row.donors} | {row.profiles} |')
text = text.replace('{{CAPACITY_TABLE}}', '\n'.join(lines))

# Step3. 逐类型数量真值附表及相关中位数，恒零失败单独标出
medians = ['| Dataset | Method | Evaluable types | Valid r | Median r | Median rho |', '|---|---|---:|---:|---:|---:|']
for name, job, placeholder in [('GSE214695', 25, '{{GSE_TYPE_TABLE}}'), ('IBDverse', 26, '{{IBD_TYPE_TABLE}}')]:
    folder = root / '08_external_pseudobulk' / name / f'06_CIBERSORTx_Job{job}_review_20261014'
    metrics = pd.read_csv(folder / 'paired_accuracy_metrics_all_scales.tsv', sep='\t')
    selected = metrics.loc[metrics.truth_scale.eq('Cell_count_SECONDARY') & metrics.evaluable & ~metrics.cell_type.str.startswith('POOLED')]
    for method in ['BayesPrism_final', f'CIBERSORTx_Job{job}']:
        block = selected.loc[selected.method.eq(method)]
        medians.append(f'| {name} | {method} | {len(block)} | {block.Pearson_r.notna().sum()} | {block.Pearson_r.median():.3f} | {block.Spearman_rho.median():.3f} |')
    cx = selected.loc[selected.method.eq(f'CIBERSORTx_Job{job}')].set_index('cell_type')
    bp = selected.loc[selected.method.eq('BayesPrism_final')].set_index('cell_type')
    lines = ['| Cell type | CX r | CX rho | CX MAE pp | BP MAE pp | CX bias pp |', '|---|---:|---:|---:|---:|---:|']
    for cell_type, row in cx.iterrows():
        r = f'{row.Pearson_r:.3f}' if pd.notna(row.Pearson_r) else 'NA'
        rho = f'{row.Spearman_rho:.3f}' if pd.notna(row.Spearman_rho) else 'NA'
        lines.append(f'| {cell_type} | {r} | {rho} | {row.MAE_pp:.3f} | {bp.loc[cell_type, "MAE_pp"]:.3f} | {row.mean_signed_error_pp:+.3f} |')
    text = text.replace(placeholder, '\n'.join(lines))
text = text.replace('{{CORRELATION_MEDIANS}}', '\n'.join(medians))
assert '{{' not in text and '}}' not in text
assert len(summary) == 18 and len(capacity) == 20
path.write_text(text, encoding='utf-8')
print('Complete discussion brief:', path)
print('Characters:', len(text), 'Bytes:', path.stat().st_size)
