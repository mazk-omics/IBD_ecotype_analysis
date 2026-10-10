# 审计 Job8 回肠 19 类等 RNA 测试，和 Job7 在同一批细胞、相同 19 类轴上配对比较。
# 在本机项目 Python 环境运行；读 Downloads、已冻结真值与 Job7，写入回肠验证目录。
# 只分析已有输出，不改 signature、留出细胞或原始 PDF；未提供日志的参数和收敛状态保持未确认。
import hashlib
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import pypdfium2 as pdfium
from pypdf import PdfReader

# Step1. 保存原始文件副本，核对后续重跑不会替换已有结果
project = Path('D:/IBD_ecotype/IBD_ecotype_analysis')
ready = project / 'CIBERSORTx_ready_20261007'
metadata = ready / '核对记录_勿上传'
old = project / 'CIBERSORTx_signature_audit/main_candidate_validation_20261006/Job7_review'
out = project / 'CIBERSORTx_signature_audit/ileum_candidate_validation_20261007/Job8_review'
originals = out / 'original_outputs'
originals.mkdir(parents=True, exist_ok=True)
file_audit = {}
for name in ['CIBERSORTx_Job8_Results.txt', 'CIBERSORTx_Job8_Results.pdf', 'CIBERSORTx_Job8_stackedChart.pdf']:
    source = Path('C:/Users/马泽凯/Downloads') / name
    target = originals / name
    source_hash = hashlib.sha256(source.read_bytes()).hexdigest()
    if target.exists():
        assert hashlib.sha256(target.read_bytes()).hexdigest() == source_hash, f'冻结副本不一致：{name}'
    else:
        shutil.copy2(source, target)
    assert hashlib.sha256(target.read_bytes()).hexdigest() == source_hash
    file_audit[name] = {'bytes': source.stat().st_size, 'sha256': source_hash}

# Step2. 按名称匹配 4 个样本与 19 类，核实输出类型与预期回肠候选一致
job8 = pd.read_csv(originals / 'CIBERSORTx_Job8_Results.txt', sep='\t', index_col=0)
job7 = pd.read_csv(old / 'original_outputs/CIBERSORTx_Job7_Results.txt', sep='\t', index_col=0)
truth = pd.read_csv(metadata / 'GI19_EqualRNA4_ground_truth_cell.tsv', sep='\t', index_col=0)
RNA_truth = pd.read_csv(metadata / 'GI19_EqualRNA4_ground_truth_RNA.tsv', sep='\t', index_col=0)
pairs = pd.read_csv(metadata / 'GI19_EqualRNA4_pairing.tsv', sep='\t').set_index('sample').loc[truth.index]
signature = pd.read_csv(ready / '01_Ileum19_validation/GI19_Signature19.txt', sep='\t', index_col=0)
integrity = pd.read_csv(metadata / 'upload_file_integrity.tsv', sep='\t').set_index('relative_path')
for name in ['GI19_Signature19.txt', 'GI19_HoldoutEqualRNA4_Mixture_CPM.txt']:
    path = ready / '01_Ileum19_validation' / name
    with path.open('rb') as handle:
        digest = hashlib.file_digest(handle, 'sha256').hexdigest()
    assert digest == integrity.loc[str(path.relative_to(ready)), 'sha256']
assert job8.index.is_unique and job8.columns.is_unique and set(job8.index) == set(truth.index)
assert set(job8.columns) == set(truth.columns) | {'P-value', 'Correlation', 'RMSE'}
assert signature.columns.equals(truth.columns) and np.allclose(truth, RNA_truth) and np.allclose(truth.sum(axis=1), 1)
pred8 = job8.loc[truth.index, truth.columns]
pred7 = job7.loc[truth.index, truth.columns]
omitted = job7.columns.difference(list(truth.columns) + ['P-value', 'Correlation', 'RMSE'])
assert omitted.tolist() == ['Colonic_absorptive'] and job7.loc[truth.index, omitted].eq(0).all().all()
assert pred8.shape == (4, 19) and np.isfinite(pred8.to_numpy()).all() and pred8.ge(0).all().all() and pred8.le(1).all().all()
assert np.allclose(pred8.sum(axis=1), 1) and np.allclose(pred7.sum(axis=1), 1)
pred8.to_csv(out / 'predicted_fractions_aligned.tsv', sep='\t')
truth.to_csv(out / 'ground_truth_equalRNA19.tsv', sep='\t')
print('Matched four paired samples and nineteen common types; omitted Joint20 fraction is exactly zero', flush=True)

# Step3. 配对比较比例误差；不使用网页表达拟合 Correlation 代替比例准确率
error8 = pred8 - truth
error7 = pred7 - truth
sample_metrics = pairs.copy()
sample_metrics['Job7_MAE_common19_pp'] = error7.abs().mean(axis=1) * 100
sample_metrics['Job8_MAE_common19_pp'] = error8.abs().mean(axis=1) * 100
sample_metrics['Job8_minus_Job7_MAE_pp'] = sample_metrics.Job8_MAE_common19_pp - sample_metrics.Job7_MAE_common19_pp
sample_metrics['Job7_TV'] = error7.abs().sum(axis=1) / 2
sample_metrics['Job8_TV'] = error8.abs().sum(axis=1) / 2
sample_metrics['Job7_fraction_pearson'] = pred7.corrwith(truth, axis=1)
sample_metrics['Job8_fraction_pearson'] = pred8.corrwith(truth, axis=1)
sample_metrics[['Job8_web_P_value', 'Job8_web_expression_correlation', 'Job8_web_expression_RMSE']] = job8.loc[truth.index, ['P-value', 'Correlation', 'RMSE']].to_numpy()
sample_metrics.to_csv(out / 'paired_sample_accuracy_metrics.tsv', sep='\t')
type_metrics = pd.DataFrame({'truth_mean_percent': truth.mean() * 100, 'Job7_mean_percent': pred7.mean() * 100, 'Job8_mean_percent': pred8.mean() * 100,
    'Job7_MAE_pp': error7.abs().mean() * 100, 'Job8_MAE_pp': error8.abs().mean() * 100, 'Job7_bias_pp': error7.mean() * 100, 'Job8_bias_pp': error8.mean() * 100,
    'Job7_zero_when_present': (pred7.eq(0) & truth.gt(0)).sum(), 'Job8_zero_when_present': (pred8.eq(0) & truth.gt(0)).sum(), 'true_present_samples': truth.gt(0).sum()})
type_metrics['Job8_minus_Job7_MAE_pp'] = type_metrics.Job8_MAE_pp - type_metrics.Job7_MAE_pp
type_metrics.index.name = 'cell_type'
type_metrics.to_csv(out / 'paired_cell_type_accuracy_metrics.tsv', sep='\t')
long_table = truth.rename_axis(index='sample', columns='cell_type').stack().rename('truth').reset_index()
for name, table in [('Job7_predicted', pred7), ('Job8_predicted', pred8)]:
    values = table.rename_axis(index='sample', columns='cell_type').stack().rename(name).reset_index()
    long_table = long_table.merge(values, on=['sample', 'cell_type'], validate='one_to_one')
long_table['Job8_bias_pp'] = (long_table.Job8_predicted - long_table.truth) * 100
long_table.to_csv(out / 'all_76_paired_fraction_comparisons.tsv', sep='\t', index=False)
print(sample_metrics[['Job7_MAE_common19_pp', 'Job8_MAE_common19_pp', 'Job7_TV', 'Job8_TV']].to_string(), flush=True)
print(type_metrics.to_string(), flush=True)

# 生物学解释先核对本地证据：标志基因已入选，不把 Stem/TA 为零简单归因于漏掉 marker。
candidate_dir = project / 'CIBERSORTx_signature_audit/three_candidate_comparison_20261006/gse_ileum'
marker_audit = pd.read_csv(candidate_dir / 'marker_panel_audit.tsv', sep='\t')
selected_markers = marker_audit.loc[marker_audit.gene.isin(['DEFA5', 'DEFA6', 'LYZ', 'MUC2', 'LGR5', 'OLFM4', 'MKI67', 'TOP2A'])]
selected_markers.to_csv(out / 'ileum_relevant_marker_specificity.tsv', sep='\t', index=False)
source_dir = project / 'CIBERSORTx_signature_audit/main_candidate_validation_20261006/server_outputs/validation_joint20_20261006'
pool = pd.read_csv(source_dir / 'selected_pool_manifest.tsv', sep='\t', usecols=['cell_id', 'final_analysis'])
membership = pd.read_csv(metadata / 'GI19_natural12_existing_cell_membership.tsv', sep='\t')
stem_members = membership.loc[membership['sample'].isin(pairs.paired_main_sample) & membership.pilot_type.eq('Stem_TA_progenitor')]
stem_members = stem_members.merge(pool, on='cell_id', validate='many_to_one')
assert stem_members.final_analysis.notna().all()
stem_counts = pd.crosstab(stem_members['sample'], stem_members.final_analysis).reindex(pairs.paired_main_sample, fill_value=0)
training = pd.read_csv(project / 'IBDverse_reference/three_reference_pilots_20261006/run_20261006_183211/gse_ileum/selected_cell_manifest.tsv', sep='\t')
training_stem = training.loc[training.pilot_type.eq('Stem_TA_progenitor'), 'author_label'].value_counts()
stem_counts.loc['Ileum19_training'] = training_stem.reindex(stem_counts.columns, fill_value=0)
stem_counts.to_csv(out / 'StemTA_author_subtype_counts.tsv', sep='\t')

# Step4. 提取并渲染两份完整 PDF，核对表格显示精度，不重新导出或修改源 PDF
render_dir = out / 'pdf_review'
render_dir.mkdir(exist_ok=True)
pdf_texts = {}
for name in ['CIBERSORTx_Job8_Results.pdf', 'CIBERSORTx_Job8_stackedChart.pdf']:
    pdf_path = originals / name
    reader = PdfReader(pdf_path)
    pdf_texts[name] = {'pages': len(reader.pages), 'texts': [page.extract_text() for page in reader.pages]}
    document = pdfium.PdfDocument(str(pdf_path))
    for page_number in range(len(document)):
        page = document[page_number]
        bitmap = page.render(scale=2.4)
        bitmap.to_pil().save(render_dir / f'{pdf_path.stem}_page{page_number + 1}.png')
        bitmap.close()
        page.close()
    document.close()
table_text = pdf_texts['CIBERSORTx_Job8_Results.pdf']['texts'][0]
display_rows = [line.split() for line in table_text.splitlines() if line.split() and line.split()[0] in truth.index]
display = pd.DataFrame([list(map(float, row[1:])) for row in display_rows], index=[row[0] for row in display_rows], columns=list(job8.columns[:-1])).loc[truth.index]
assert np.allclose(display, job8.loc[truth.index, display.columns], atol=0.00051)
(render_dir / 'pdf_text_extraction.json').write_text(json.dumps(pdf_texts, ensure_ascii=False, indent=2), encoding='utf-8')

# Step5. 保存审计摘要；实际网页 signature 文件与参数仍需要 runtime log 确认
MAE7 = float(error7.abs().to_numpy().mean() * 100)
MAE8 = float(error8.abs().to_numpy().mean() * 100)
summary = {'status': 'numeric_audit_complete_pdf_visual_review_pending', 'date': '2026-10-07', 'job': 'Job8', 'samples': 4, 'common_cell_types': 19,
    'observed_sample_and_type_labels_match_GI19_EqualRNA4': True, 'all_76_pairs_aligned': True, 'all_fractions_valid_and_sum_to_one': True,
    'Job7_MAE_common19_pp': MAE7, 'Job8_MAE_common19_pp': MAE8, 'descriptive_relative_MAE_reduction_percent': (MAE7 - MAE8) / MAE7 * 100,
    'all_four_sample_MAEs_improved': bool(sample_metrics.Job8_minus_Job7_MAE_pp.lt(0).all()), 'types_with_lower_MAE': int(type_metrics.Job8_minus_Job7_MAE_pp.lt(0).sum()),
    'Job7_mean_TV': float(sample_metrics.Job7_TV.mean()), 'Job8_mean_TV': float(sample_metrics.Job8_TV.mean()),
    'Paneth_mean_truth_percent': float(truth.Paneth.mean() * 100), 'Paneth_Job7_mean_percent': float(pred7.Paneth.mean() * 100), 'Paneth_Job8_mean_percent': float(pred8.Paneth.mean() * 100),
    'Stem_TA_Job8_zero_when_present': int(type_metrics.loc['Stem_TA_progenitor', 'Job8_zero_when_present']),
    'web_expression_correlation_median': float(job8.Correlation.median()), 'web_P_values_reported_zero': int(job8['P-value'].eq(0).sum()),
    'run_settings_and_actual_uploaded_signature_independently_confirmed': False, 'solver_convergence_confirmed': False,
    'natural_bulk_accuracy_validated': False, 'original_outputs': file_audit,
    'decision': 'Ileum19 improves these four equal-RNA holdouts compared with Joint20; retain as regional candidate, with remaining Stem/TA, Paneth, DC and EEC errors; await natural-count validation.',
    'simulation_samples_are_independent_patients': False, 'new_sampling_or_server_work_performed': False}
visual_path = render_dir / 'pdf_visual_review.json'
if visual_path.exists():
    visual_review = json.loads(visual_path.read_text(encoding='utf-8'))
    assert visual_review['source_sha256'] == {name: info['sha256'] for name, info in file_audit.items() if name.endswith('.pdf')}
    summary['status'] = 'review_completed'
    summary['PDF_visual_review'] = visual_review
summary['PDF_table_matches_TXT_at_display_precision'] = True
(out / 'review_summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')

# Step6. 写出中文报告，改善与残留误差同时展示，不按平均值给候选宣布通过
report = ['# Job8 回肠等 RNA 对照审阅（2026-10-07）', '',
    '结论：在相同的四个等 RNA 留出样本上，Ileum19 比 Joint20 改善了整体比例误差，尤其是 Paneth、Goblet、Tuft 和 NK。Stem/TA 仍在三个样本为零，Paneth 在上皮富集样本仍明显高估，DC 和 Enteroendocrine 的误差增大。保留为优先回肠候选，不能认定为自然 bulk 或所有细胞类型均已通过验证。', '',
    '## 配对范围与完整性', '',
    'Job8 包含 4 个回肠等 RNA 样本、19 类，样本名单与输出类型和既定 GI19_EqualRNA4 一致；76 项比例全部按名称对齐，有限、非负、和为 1。保存真值来自已审计的细胞名单，等 RNA 贡献等于细胞数比例，未使用 signature 或预测 fraction 构造真值。输入与 signature 的本地 SHA256 与准备清单一致。', '',
    'Job7 原本有 20 类；本次仅使用相同 4 个样本、共有 19 类，删去的 Colonic_absorptive 在这些样本中预测和真值均为零。没有对 Joint20 预测重新归一化。两边 MAE 均以 19 类为分母，因此 Job7 数值与上次按 20 类报告的 3.24 个百分点不同。', '',
    '用户未提供 Job8 日志和参数截图，实际上传的 signature 文件、batch correction、QN、置换次数及求解器收敛未独立核验。参考效果比较以按原指引保持其他设置为前提。Ileum19 与 Joint20 的训练细胞、供者、特征集合和类型数均不同，不能把改善纯因果归因于只改变部位。', '',
    '## 同一批样本、相同类型轴的误差', '',
    '| 样本 | Job7 MAE（百分点） | Job8 MAE（百分点） | Job7 TV | Job8 TV |', '|---|---:|---:|---:|---:|']
for sample, row in sample_metrics.iterrows():
    report.append(f'| {sample} | {row.Job7_MAE_common19_pp:.2f} | {row.Job8_MAE_common19_pp:.2f} | {row.Job7_TV:.3f} | {row.Job8_TV:.3f} |')
report += [f'| 平均 | {MAE7:.2f} | {MAE8:.2f} | {sample_metrics.Job7_TV.mean():.3f} | {sample_metrics.Job8_TV.mean():.3f} |', '',
    f'四个样本都改善，共有 {summary["types_with_lower_MAE"]}/19 类的平均绝对误差下降；整体 MAE 描述性下降 {summary["descriptive_relative_MAE_reduction_percent"]:.1f}%。这些混合物共享供者和部分细胞，不是四位独立患者，不做配对患者显著性推断。', '',
    '## 代表类型的均值与误差', '', '| 类型 | 真值均值 | Job7 均值 | Job8 均值 | Job7 MAE（百分点） | Job8 MAE（百分点） |', '|---|---:|---:|---:|---:|---:|']
for name in ['Paneth', 'Goblet', 'Stem_TA_progenitor', 'Tuft', 'Ileal_absorptive', 'DC', 'Enteroendocrine', 'NK', 'CD4_T', 'Monocyte_Macrophage']:
    row = type_metrics.loc[name]
    report.append(f'| {name} | {row.truth_mean_percent:.2f}% | {row.Job7_mean_percent:.2f}% | {row.Job8_mean_percent:.2f}% | {row.Job7_MAE_pp:.2f} | {row.Job8_MAE_pp:.2f} |')
report += ['', '均值接近不能替代逐样本准确性，因此同时提供逐类型 MAE 和全部 76 项对照。Paneth 平均偏差虽缩小，仍有实际真值的约两倍；DC、EEC 与 Stem/TA 仍需重点检查。', '',
    '## Paneth 的逐样本变化', '', '| 场景 | 真值 | Job7 | Job8 |', '|---|---:|---:|---:|']
for row in long_table.loc[long_table.cell_type.eq('Paneth')].itertuples():
    report.append(f'| {row.sample} | {row.truth * 100:.2f}% | {row.Job7_predicted * 100:.2f}% | {row.Job8_predicted * 100:.2f}% |')
report += ['',
    '回肠免疫富集样本的 Paneth 已接近真值；上皮富集样本仍为 34.18%，真值为 16%。Goblet 从前三个样本为零改善为四个均非零；Tuft 四个样本也更接近真值。Stem/TA 则仍是 balanced 1/2、epithelial_high 1 三个样本为零，仅 immune_high 1 预测 2.34%（真值 4%）。EEC 在上皮富集样本仅预测 5.23%（真值 15%）。这些是不同类型的残留校准问题，不能简单认为全被 Paneth 吸收。', '',
    f'网页 Correlation 中位数 {summary["web_expression_correlation_median"]:.3f}，四个 P-value 均显示为 0；表达拟合指标不能代替上述比例准确率或证明求解器收敛。', '',
    '## 下一步', '',
    '继续提交已经准备的 GI19_HoldoutNatural12_Mixture_CPM.txt，使用相同 GI19_Signature19.txt，relative、QN disabled、batch correction disabled、100 permutations，并保存完整日志和参数。等 RNA 改善还不能证明自然 RNA 贡献或真实 MSCCR 细胞比例已准确；自然 counts 输入是下一项必要对照。', '',
    'SCP259 的结直肠两项测试可按原计划随后进行。当前不需要因这四个结果立即改 signature、删基因、加细胞类型或运行新的抽样；先完成已固定测试，保留未使用的 CD12、UC17 最终检验供者。', '',
    '原始 TXT/PDF 在 original_outputs；全部结果见 all_76_paired_fraction_comparisons.tsv、paired_sample_accuracy_metrics.tsv、paired_cell_type_accuracy_metrics.tsv。新写研究脚本为 scripts/python/review_cibersortx_job8_ileum_equalRNA.py。']
report += ['', '## 回肠生物学的解释与证据边界', '',
    '人类肠道图谱显示部位间存在显著的细胞组成和转录程序差异，同一上皮谱系沿分化路径的表达也随部位变化。这使混合回肠和非回肠来源的参考可能发生表达代表性或幅度失配。[Organization of the human intestine at single-cell resolution](https://www.nature.com/articles/s41586-023-05915-x)。', '',
    '人类健康肠上皮研究还发现，常用 stem marker 并不特异标记所有人类 stem 细胞，LYZ 也不是人类 Paneth 独有标志。这支持检查多个 marker 的联合特异性，不能套用单一标志或鼠肠标准。[A Proximal-to-Distal Survey of Healthy Adult Human Small Intestine and Colon Epithelium](https://pmc.ncbi.nlm.nih.gov/articles/PMC9043569/)。这些文献提供机制背景，不直接证明本次偏差的原因。', '',
    '本地实测：Ileum19 中 LGR5、OLFM4、MKI67、TOP2A 都在 signature 内。既有审计的 LGR5 目标/其他最大值比（均加 1）约 1.24，OLFM4 约 2.05；MKI67 的最高值在 Goblet 列。这说明某些常用 marker 在当前参考中的类型区分并不很强，但不能据此把整个特征矩阵的误差归因于单个基因。', '',
    'Stem/TA 合并类的训练组成是 23 个 Ileal LGR5pos stem、27 个 Ileal TA；四个留出样本的 stem/TA 数分别为 5/6、5/5、8/4、3/5。最后一个样本有非零预测，而更富 stem 的上皮样本仍为零，简单的 stem/TA 数量比例变化并不足以完整解释结果。这里使用原作者 final_analysis 标签；minor 列中的 Ileal_enterocyte 对应较粗分组，不能误当成把成熟 enterocyte 加入了 Stem/TA。', '',
    '当前更合理的待检验机制是：回肠部位与谱系/功能状态差异影响了参考表达代表性，加上标度、共线性、注释及算法适配造成比例偏差。等 RNA 已控制每细胞总 RNA 贡献，但没有消除基因组成和状态差异。由于真值来自同一批已标注细胞，偏差仍是需要解释和改进的预测误差；不能用回肠特殊来跳过验证。先完成自然 counts 对照，再决定亚型平衡或参考重建。']
(out / 'Job8回肠等RNA对照审阅.md').write_text('\n'.join(report) + '\n', encoding='utf-8')
print(json.dumps(summary, ensure_ascii=True, indent=2), flush=True)
