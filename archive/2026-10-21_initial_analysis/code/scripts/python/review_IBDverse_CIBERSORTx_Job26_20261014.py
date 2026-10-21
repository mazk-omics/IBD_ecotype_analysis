# 将用户提供的官方Job26与同批100份外部模拟真值、BayesPrism配对，分别报告三种真值尺度与逐类型偏差。
# 本机科学Python运行；只新增审阅目录，历史输入和旧报告保留，所有误差统一为百分点。
from pathlib import Path
import hashlib
import json
import shutil
import numpy as np
import pandas as pd
from scipy import stats

# Step1. 保存官方原件与参数来源，严格核对100个样本和20类
project = Path(__file__).resolve().parents[2]
root = project / '06_validation/deconv_fraction_sensitivity_joint20_20261013'
dataset = root / '08_external_pseudobulk/IBDverse'
out = dataset / '06_CIBERSORTx_Job26_review_20261014'
out.mkdir(exist_ok=False)
(out / 'original_outputs').mkdir()
types = (root / '00_audit/frozen_types.txt').read_text().splitlines()
audit = json.loads((dataset / '04_accuracy/accuracy_summary.json').read_text())
evaluable = audit['evaluable_type_names']
unsupported = audit['unsupported_types']
sources = [Path.home() / 'Downloads' / name for name in ['CIBERSORTx_Job26_Results.txt', 'CIBERSORTx_Job26_Results.pdf', 'CIBERSORTx_Job26_stackedChart.pdf']]
inventory = []
for path in sources:
    shutil.copyfile(path, out / 'original_outputs' / path.name)
    inventory.append({'path': str(path), 'bytes': path.stat().st_size, 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()})
pd.DataFrame(inventory).to_csv(out / 'official_input_inventory.tsv', sep='\t', index=False)
official = pd.read_csv(sources[0], sep='\t', index_col=0)
cell = pd.read_csv(dataset / '02_simulation/truth_cell_count_fraction.tsv', sep='\t', index_col=0)[types]
RNA = pd.read_csv(dataset / '04_accuracy/truth_BP_actual_retained_gene_RNA_fraction.tsv', sep='\t', index_col=0).loc[cell.index, types]
full_RNA = pd.read_csv(dataset / '02_simulation/truth_full_gene_RNA_fraction.tsv', sep='\t', index_col=0).loc[cell.index, types]
initial = pd.read_csv(dataset / '03_bayesprism/theta_initial.tsv', sep='\t', index_col=0).loc[cell.index, types]
final = pd.read_csv(dataset / '03_bayesprism/theta_final.tsv', sep='\t', index_col=0).loc[cell.index, types]
assert official.shape == (100, 23) and official.index.is_unique and set(official.index) == set(cell.index)
assert set(official.columns) == set(types + ['P-value', 'Correlation', 'RMSE'])
cx = official.loc[cell.index, types]
assert np.isfinite(cx).all().all() and cx.ge(0).all().all() and np.allclose(cx.sum(axis=1), 1, atol=1e-8)
settings = {'job': 26, 'official_output_received': True, 'samples': 100, 'types': 20, 'evaluable_types': 16, 'sample_ID_exact_set_match': True, 'renormalized': False, 'signature': 'GJ20v2_Signature20.txt', 'relative': True, 'QN': False, 'batch_correction': False, 'permutations': 100, 'parameter_evidence': 'User explicitly confirmed all listed settings in this chat on 2026-10-14', 'official_runtime_log_available': False, 'exact_official_web_build': 'UNKNOWN', 'uploaded_server_file_hash_attestation': 'No official upload hash; local prepared inputs and signature hashes retained', 'SSH_used': False}
(out / 'official_QC_and_settings.json').write_text(json.dumps(settings, indent=2), encoding='utf-8')
print('Job26 QC:', cx.shape, 'row-sum range', cx.sum(axis=1).min(), cx.sum(axis=1).max())

# Step2. 三种真值和initial/final分别评价，可评价16类不重新归一化
metrics = []
paired_samples = []
paired_types = []
long_tables = []
for scale, truth in [('Retained_gene_RNA_PRIMARY', RNA), ('Cell_count_SECONDARY', cell), ('Full_gene_RNA_SECONDARY', full_RNA)]:
    for method, prediction in [('BayesPrism_final', final), ('BayesPrism_initial', initial), ('CIBERSORTx_Job26', cx)]:
        for name in types + ['POOLED_EVALUABLE', 'POOLED_FULL20']:
            axis = evaluable if name == 'POOLED_EVALUABLE' else types
            x = prediction[axis].to_numpy().ravel() if name.startswith('POOLED') else prediction[name].to_numpy()
            y = truth[axis].to_numpy().ravel() if name.startswith('POOLED') else truth[name].to_numpy()
            error = (x - y) * 100
            variable = x.std() > 0 and y.std() > 0
            metrics.append({'truth_scale': scale, 'method': method, 'cell_type': name, 'evaluable': name in evaluable or name == 'POOLED_EVALUABLE', 'n_mixtures': 100, 'Pearson_r': stats.pearsonr(x, y).statistic if variable else np.nan, 'Spearman_rho': stats.spearmanr(x, y).statistic if variable else np.nan, 'MAE_pp': np.abs(error).mean(), 'RMSE_pp': np.sqrt(np.mean(error**2)), 'mean_signed_error_pp': error.mean(), 'predicted_mean_percent': x.mean() * 100, 'truth_mean_percent': y.mean() * 100})
        long = prediction.stack().rename('predicted_fraction')
        target = truth.stack().rename('truth_fraction')
        long.index.names = ['sample', 'cell_type']
        target.index.names = ['sample', 'cell_type']
        assert long.index.is_unique and target.index.is_unique
        joined = long.to_frame().join(target).reset_index()
        joined['truth_scale'] = scale
        joined['method'] = method
        long_tables.append(joined)
    bp_error = (final - truth).abs() * 100
    cx_error = (cx - truth).abs() * 100
    sample = pd.DataFrame({'BP_MAE_pp': bp_error[evaluable].mean(axis=1), 'CX_MAE_pp': cx_error[evaluable].mean(axis=1)})
    sample['CX_minus_BP_MAE_pp'] = sample.CX_MAE_pp - sample.BP_MAE_pp
    sample['truth_scale'] = scale
    paired_samples.append(sample.reset_index(names='sample'))
    paired = pd.DataFrame({'cell_type': types, 'BP_MAE_pp': bp_error.mean().to_numpy(), 'CX_MAE_pp': cx_error.mean().to_numpy()})
    paired['CX_minus_BP_MAE_pp'] = paired.CX_MAE_pp - paired.BP_MAE_pp
    paired['evaluable'] = paired.cell_type.isin(evaluable)
    paired['truth_scale'] = scale
    paired_types.append(paired)
metrics = pd.DataFrame(metrics)
samples = pd.concat(paired_samples, ignore_index=True)
type_pairs = pd.concat(paired_types, ignore_index=True)
metrics.to_csv(out / 'paired_accuracy_metrics_all_scales.tsv', sep='\t', index=False)
samples.to_csv(out / 'paired_sample_MAE.tsv', sep='\t', index=False)
type_pairs.to_csv(out / 'paired_type_MAE.tsv', sep='\t', index=False)
pd.concat(long_tables, ignore_index=True).to_csv(out / 'predictions_and_truth_long.tsv', sep='\t', index=False)

# Step3. 表达拟合与比例准确性分开；审计零预测、未支持类型和实际基因分母
diagnostics = official.loc[cell.index, ['P-value', 'Correlation', 'RMSE']]
diagnostics.to_csv(out / 'official_expression_fit_diagnostics.tsv', sep='\t', index_label='sample')
zero = pd.DataFrame({'cell_type': types, 'CX_zero_predictions': cx.eq(0).sum().to_numpy(), 'BP_zero_predictions': final.eq(0).sum().to_numpy(), 'cell_truth_positive_mixtures': cell.gt(0).sum().to_numpy(), 'CX_zero_when_cell_truth_positive': (cx.eq(0) & cell.gt(0)).sum().to_numpy(), 'BP_zero_when_cell_truth_positive': (final.eq(0) & cell.gt(0)).sum().to_numpy()})
zero.to_csv(out / 'zero_prediction_audit.tsv', sep='\t', index=False)
absent = pd.DataFrame({'CIBERSORTx_Job26': cx[unsupported].sum(axis=1), 'BayesPrism_final': final[unsupported].sum(axis=1)})
absent.to_csv(out / 'unsupported_type_predicted_mass.tsv', sep='\t', index_label='sample')
denominator = pd.read_csv(dataset / '04_accuracy/actual_gene_denominator_audit.tsv', sep='\t', index_col=0).loc[cell.index]
retained_counts = RNA.mul(denominator.retained_gene_counts, axis=0).sum()
full_counts = full_RNA.mul(denominator.full_gene_counts, axis=0).sum()
gene_retention = pd.DataFrame({'full_gene_counts_over_all_simulations': full_counts, 'BP_retained_gene_counts_over_all_simulations': retained_counts, 'BP_retained_fraction_of_type_RNA': retained_counts / full_counts})
gene_retention.to_csv(out / 'gene_denominator_by_type.tsv', sep='\t', index_label='cell_type')
pooled = metrics.loc[metrics.cell_type.eq('POOLED_EVALUABLE')].copy()
pooled.to_csv(out / 'summary_16_evaluable_types.tsv', sep='\t', index=False)
summary = {'status': 'COMPLETE', 'official_job': 26, 'mixtures': 100, 'evaluable_types': 16, 'inference_types': 20, 'parameters_user_confirmed': True, 'results': pooled.to_dict('records'), 'unsupported_mass_mean': absent.mean().to_dict(), 'official_expression_fit_correlation_median': float(diagnostics.Correlation.median()), 'official_P_value_zero_count': int(diagnostics['P-value'].eq(0).sum()), 'simulation_samples_not_independent_biological_donors': True, 'primary_endpoint_unchanged': True, 'comparative_ranking_depends_on_truth_gene_denominator': False, 'pairing_evidence': 'User supplied official IBDverse results; exact prepared SampleIDs and Full20 matched; upload hash not independently available'}
(out / 'paired_summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')

# Step4. 报告保留原定主要终点，并明确总体误差与逐类型误差的区别
lines = ['# IBDverse 官方 CIBERSORTx Job26 配对验证', '', '2026-10-14：官方TXT与两份官方PDF已归档；100个SampleID及完整20类完全匹配，非负且每行和为1，没有重新归一化或剔除样本。用户在本对话确认原GJ20v2 signature、Relative、QN关闭、batch correction关闭、100 permutations；完整运行日志、网页build和服务器上传hash未取得。', '', '主要终点仍按2026-10-13预定：相对BayesPrism实际保留基因的RNA-contribution truth。cell-count与完整基因RNA truth并列作为预定次要尺度。主要指标仅用16个作者可评价类型，保留完整20类推断；不以缺失类型零列降低主要误差。', '', '| 真值 | 方法 | MAE pp | RMSE pp | pooled Pearson | pooled Spearman |', '|---|---|---:|---:|---:|---:|']
for row in pooled.loc[~pooled.method.eq('BayesPrism_initial')].itertuples():
    lines.append(f'| {row.truth_scale} | {row.method} | {row.MAE_pp:.3f} | {row.RMSE_pp:.3f} | {row.Pearson_r:.3f} | {row.Spearman_rho:.3f} |')
lines += ['', '## 判读', '', 'CIBERSORTx Job26在三种真值尺度的总体MAE均低于BayesPrism final。主要RNA尺度的优势集中在Tuft与Stem/TA：BP对Tuft高估、对Stem/TA低估。逐类型结果仍需完整保留，不把总体更低误差解释为每类均更好。BP保留基因truth是方法特异性尺度；完整基因RNA与cell-count分别保存，预定主要终点不因结果更改。', '', f'Plasma全部模拟RNA中仅{gene_retention.loc["Plasma", "BP_retained_fraction_of_type_RNA"] * 100:.2f}%落入BP保留基因集合；各类型保留程度不同，真值变换不是统一常数缩放。该分母差异是实际可核验的现象，不能自动将Plasma相对保留基因truth的高估都解释为算法失效。具体见gene_denominator_by_type.tsv。', '', '## 逐类型 cell-count 与主要RNA结果', '', '| 类型 | BP RNA MAE | CX RNA MAE | BP cell MAE | CX cell MAE | CX RNA signed bias |', '|---|---:|---:|---:|---:|---:|']
for name in evaluable:
    r = metrics.loc[metrics.cell_type.eq(name)].set_index(['truth_scale', 'method'])
    lines.append(f'| {name} | {r.loc[("Retained_gene_RNA_PRIMARY", "BayesPrism_final"), "MAE_pp"]:.3f} | {r.loc[("Retained_gene_RNA_PRIMARY", "CIBERSORTx_Job26"), "MAE_pp"]:.3f} | {r.loc[("Cell_count_SECONDARY", "BayesPrism_final"), "MAE_pp"]:.3f} | {r.loc[("Cell_count_SECONDARY", "CIBERSORTx_Job26"), "MAE_pp"]:.3f} | {r.loc[("Retained_gene_RNA_PRIMARY", "CIBERSORTx_Job26"), "mean_signed_error_pp"]:+.3f} |')
stem = zero.set_index('cell_type').loc['Stem_TA_progenitor']
ilc = zero.set_index('cell_type').loc['ILC']
primary_pairs = type_pairs.loc[type_pairs.truth_scale.eq('Retained_gene_RNA_PRIMARY') & type_pairs.evaluable]
lines += ['', f'主要RNA尺度上，CX有{int(primary_pairs.CX_minus_BP_MAE_pp.lt(0).sum())}/16类MAE更低，BP有{int(primary_pairs.CX_minus_BP_MAE_pp.gt(0).sum())}/16类更低；总体优势被少数大误差谱系影响。', '', f'ILC：CX在100份中全部预测为0；其中{int(ilc.cell_truth_positive_mixtures)}份实际含ILC。该系统漏检独立于Tuft/Stem改善，需要作为限制报告。']
lines += ['', f'Stem/TA：官方预测{int(stem.CX_zero_predictions)}/100份为0；在{int(stem.cell_truth_positive_mixtures)}份实际含该类细胞的混合物中，{int(stem.CX_zero_when_cell_truth_positive)}份被预测为0。不能因全局表达重建相关高而忽略该漏检。', '', '## 配对分布与限制', '']
for scale, block in samples.groupby('truth_scale', sort=False):
    lines.append(f'- {scale}：BP每混合物MAE较低{int(block.CX_minus_BP_MAE_pp.gt(0).sum())}/100，CX较低{int(block.CX_minus_BP_MAE_pp.lt(0).sum())}/100。仅描述模拟分布，不作供者独立性显著性检验。')
lines += ['', f'未设计的NK、Colonic_absorptive、Paneth和Glial预测质量：CX均值{absent.CIBERSORTx_Job26.mean() * 100:.3f}%、最大{absent.CIBERSORTx_Job26.max() * 100:.3f}%；BP均值{absent.BayesPrism_final.mean() * 100:.3f}%。这些类型不进入主要16类终点，完整20类指标另存。', '', f'官方输出表达重建Correlation中位数{diagnostics.Correlation.median():.3f}；100份显示P-value=0是有限置换输出，不代表数学概率为0，更不证明各细胞比例准确。官方RMSE为表达拟合量，未与本报告比例RMSE混用。全部100份保留。', '', '100个模拟复用来源细胞和供者，不是100个独立生物学重复；不据普通100样本bootstrap制造生物学置信区间。跨研究供者身份核验限制及作者类型映射限制沿用原审计。GSE214695 Job25配对已完成；两个数据集都已取得官方无校正结果。此次仅本地审阅，无SSH、无新signature、无HiRes/G-like、无冻结states/ecotypes修改。', '', '原件与SHA256：official_input_inventory.tsv；全套逐类型Pearson/Spearman/MAE/RMSE/signed bias：paired_accuracy_metrics_all_scales.tsv；逐样本与逐类型配对误差：paired_sample_MAE.tsv、paired_type_MAE.tsv。']
(out / 'IBDverse_JOB26_PAIRED_REPORT.md').write_text('\n'.join(lines) + '\n', encoding='utf-8')
shutil.copyfile(Path(__file__), out / Path(__file__).name)
print(pooled[['truth_scale', 'method', 'MAE_pp', 'RMSE_pp', 'Pearson_r', 'Spearman_rho']].to_string(index=False))
print('REPORT', out / 'IBDverse_JOB26_PAIRED_REPORT.md')
