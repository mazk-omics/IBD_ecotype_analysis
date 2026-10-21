# 汇总已完成的BayesPrism、表达QC和EcoTyper等待项，保留昨日启动记录不覆盖。
# 服务器ibd_refaudit环境运行，输出本次复核表及中文报告；不启动发现分析。
from pathlib import Path
import hashlib
import json
from datetime import datetime
import numpy as np
import pandas as pd

# Step1. 核对退出码、冻结样本和完整theta
root = Path('/home/mazekai/IBD_EcoTyper/06_validation/joint20_bayesprism_ecotyper_matched_sensitivity_20261014')
out = root / '05_completion_review_20261015'
exit_codes = {name: int((root / 'logs' / name).read_text()) for name in ['formal_exit_status.txt', 'postprocessing_exit_status.txt', 'exit_status.txt']}
assert all(value == 0 for value in exit_codes.values())
annotation = pd.read_csv(root / '00_audit/MSCCR_CD_UC_discovery_manifest.tsv', sep='\t', dtype=str)
frozen = json.loads((root / '00_audit/input_audit.json').read_text())
ids = annotation.ID.tolist()
assert len(ids) == len(set(ids)) == frozen['samples']
assert hashlib.sha256(('\n'.join(ids) + '\n').encode()).hexdigest() == frozen['ordered_sample_sha256']
theta = pd.read_csv(root / '02_bayesprism/theta_final.tsv', sep='\t', index_col=0)
assert theta.index.tolist() == ids and theta.shape == (len(ids), 20)
assert np.isfinite(theta.to_numpy()).all() and (theta.to_numpy() >= 0).all()
row_sum_error = float(np.abs(theta.sum(axis=1) - 1).max())
assert row_sum_error < 1e-6
assert 'EXPRESSION_READINESS_QC_COMPLETE' in (out / 'expression_review_status.txt').read_text()

# Step2. 汇总预定低mass标记及质量守恒，不据此排除样本或父类型
mass = pd.read_csv(root / '03_expression_qc/posterior_mass_by_sample_type.tsv', sep='\t')
mass = mass.merge(annotation[['ID', 'PatientID', 'Disease', 'Inflammation', 'Region3']], on='ID', validate='many_to_one')
assert len(mass) == len(ids) * 20
grouped = mass.groupby(['cell_type', 'Disease', 'Inflammation', 'Region3'], observed=True).agg(samples=('ID', 'size'), patients=('PatientID', 'nunique'), median_mass=('posterior_mass', 'median'), n_mass_lt100=('mass_lt100', 'sum'), n_mass_lt1000=('mass_lt1000', 'sum')).reset_index()
grouped['pct_mass_lt100'] = grouped.n_mass_lt100 / grouped.samples * 100
grouped['pct_mass_lt1000'] = grouped.n_mass_lt1000 / grouped.samples * 100
grouped.to_csv(out / 'low_mass_by_disease_inflammation_region.tsv', sep='\t', index=False)
by_sample = mass.groupby(['ID', 'PatientID', 'Disease', 'Inflammation', 'Region3'], observed=True).agg(types_mass_lt100=('mass_lt100', 'sum'), types_mass_lt1000=('mass_lt1000', 'sum')).reset_index()
by_sample.to_csv(out / 'sample_information_flags.tsv', sep='\t', index=False)
conservation = pd.read_csv(root / '03_expression_qc/Z_mass_conservation.tsv', sep='\t')
assert conservation.within_known_integer_sampler_bound.all() and conservation.max_integer_gene_abs_residual.max() <= 1e-6
type_qc = pd.read_csv(out / 'parent_type_information_QC.tsv', sep='\t')
assert len(type_qc) == 20 and type_qc.sample_order_exact.all() and type_qc.finite_nonnegative.all()
panel_review = pd.read_csv(out / 'cross_parent_marker_panel_review.tsv', sep='\t')
own = panel_review.loc[panel_review.expected_parent.eq(panel_review.observed_parent)].copy()
tops = panel_review.loc[panel_review.observed_rank.eq(1), ['expected_parent', 'observed_parent']].groupby('expected_parent').observed_parent.agg(','.join)
own['highest_scoring_parent_all_samples'] = own.expected_parent.map(tops)
own.to_csv(out / 'own_parent_marker_rank_summary.tsv', sep='\t', index=False)
missing_markers = pd.read_csv(out / 'marker_gene_coverage.tsv', sep='\t')
missing_markers = missing_markers.loc[~missing_markers.in_retained_genes, 'gene'].tolist()

# Step3. 保存EcoTyper候选接口清单；20个候选文件不代表选择20类做NMF
manifest = type_qc[['cell_type', 'genes', 'samples', 'sample_order_exact', 'finite_nonnegative', 'zero_mass_samples', 'pct_mass_lt100', 'pct_mass_lt1000']].copy()
manifest['candidate_RDS'] = [str(root / '03_expression_qc/candidate_glike_rds' / f'{ct}.rds') for ct in manifest.cell_type]
manifest['matrix_orientation'] = 'gene_by_sample'
manifest['expression_scale'] = 'initial_Z / sum_gene_initial_Z * 1e6; unlogged conditional CPM-like'
manifest['official_presorted_text_export'] = 'PENDING_FROZEN_PARENT_LIST'
manifest['discovery_release'] = 'HOLD_MAINLINE_CONFIG_AND_BIOLOGICAL_REVIEW'
manifest.to_csv(out / 'ECOTYPER_CANDIDATE_INPUT_MANIFEST.tsv', sep='\t', index=False)
status = {'checked_at': datetime.now().astimezone().isoformat(), 'formal_completed_at': (root / 'logs/current_stage.txt').read_text().strip(), 'all_stages_completed_at': (root / 'logs/completed_at.txt').read_text().strip(), 'exit_codes': exit_codes, 'samples': len(ids), 'disease_counts': annotation.Disease.value_counts().to_dict(), 'patients': int(annotation.PatientID.nunique()), 'types': 20, 'genes_per_type': 16195, 'theta_max_row_sum_error': row_sum_error, 'zero_mass_sample_type_pairs': int(mass.zero_mass.sum()), 'mass_conservation_all_pass': True, 'max_integer_gene_Z_residual': float(conservation.max_integer_gene_abs_residual.max()), 'max_fractional_loss_bound_counts': float(conservation.fractional_loss_bound.max()), 'max_abs_Z_fraction_minus_initial_theta': float(pd.read_csv(root / '03_expression_qc/Z_theta_consistency.tsv', sep='\t').iloc[:, 1].max()), 'candidate_matrices_numerically_valid': True, 'ecotyper_nmf_started': False, 'configuration': 'HOLD_PENDING_MAINLINE_FREEZE; user confirmed no mainline config path on 2026-10-15', 'full_HiRes': 'No matched full-2025 input found in reviewed project paths', 'sample_or_type_exclusions_added': 0}
(out / 'completion_status.json').write_text(json.dumps(status, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
lines = ['| 父类型 | mass中位数 | mass<100 | mass<1000 | 最终fraction中位数 |', '|---|---:|---:|---:|---:|']
for row in type_qc.itertuples():
    lines.append(f'| {row.cell_type} | {row.mass_median:.2f} | {row.n_mass_lt100}/{row.samples} | {row.n_mass_lt1000}/{row.samples} | {row.median_final_fraction:.6g} |')
marker_lines = ['| 面板对应类型 | 本类型marker排序 | 参考中本类型排序 | 全样本G-like最高类型 |', '|---|---:|---:|---|']
for row in own.itertuples():
    marker_lines.append(f'| {row.expected_parent} | {row.observed_rank} | {row.reference_rank} | {row.highest_scoring_parent_all_samples} |')
report = f'''# BayesPrism完成复核与EcoTyper衔接状态（2026-10-15）

检查时间：{status['checked_at']}。正式推断与候选表达数值QC均已完成，formal/postprocessing/overall退出码全部为0。正式阶段结束于2026-10-15 11:34:52，全流程结束于11:37:20。正式阶段耗时14:39:45，完整流程约14小时42分钟；已正常退出，当前未遗留本任务R或tmux进程。原stdout中的预计时间仅为采样阶段估计，现以实际完成标记为准。

## 已完成及数值核验

- 原CD 1154、UC 871，共2025份活检、919位患者；SampleID顺序SHA256与昨日冻结manifest完全一致，没有新增样本/父类型排除。
- 完整BayesPrism对象已保存（约11 GB）；初始/最终theta及CV、reference update、初始posterior Z均已保存。实际控制参数：8核、1000迭代、500 burn-in、thinning=2、内部seed=123、alpha=1、MAP/sigma=2；外层R seed=20261013。
- 20个候选G-like矩阵均独立回读，轴为16195 genes × 2025 samples，样本及基因顺序正确，符号唯一，有限非负，每样本和为1e6；zero posterior mass为0/40500。
- 最终theta行和最大误差{row_sum_error:.3g}。全部质量守恒检查符合预定整数采样容差；整数基因最大残差{status['max_integer_gene_Z_residual']:.3g}。Z质量比例与初始theta最大绝对差{status['max_abs_Z_fraction_minus_initial_theta']:.3g}，两者定义相关但非强制逐值相等。
- 输入小数未舍入；内部rmultinom整数化限制仍保留。质量守恒通过不证明细胞类型生物学特异性、低mass可靠性或与HiRes等价。
- 候选表达为初始Z按sample/type归一化的G-like，未做log或Region3标准化，未使用最终theta重新构造表达。EcoTyper源码会自行进行log2或log2(x+1)处理，因此导出前不能重复log。

## 低信息量父类型

昨天预定的mass<100和<1000标记仅用于审计，未据此自动删样本、删类型或确定NMF父类型。分组表已按Disease/Inflammation/Region3及PatientID汇总，不能只按全队列平均值选择父类型。

{chr(10).join(lines)}

ILC全部2025样本mass<100，不能因G-like被归一化到1e6就视为有充足信息；Stem/TA、Paneth、NK、DC等也需重点审阅。组织特异类型在不匹配部位低mass不自动代表算法失败。没有将旧版Core12列表套用到本次结果。

## Marker与跨父类型程序分布

使用昨日计划中的marker面板。ILC的IL7R/KIT、DC的CTSS以及T/NK等存在共享表达，属于背景审阅，不是专有身份判据。面板分数定义为各marker平均G-like的log2(x+1)均值；同时保存mass>=1000分层、原始Z contribution和同基因空间的参考counts汇总CPM。最高类型或排序偏移只是审阅线索，不能单凭该表宣称污染或错配。

{chr(10).join(marker_lines)}

完整基因与面板表见marker_distribution_by_type.tsv、cross_parent_marker_panel_review.tsv。20个候选表达的数值合格不等于全部20类均适合NMF；本轮未自动选择或排除父类型，也未依据尚未运行的NMF结果筛选。

本次18/20个面板在对应父类型最高。CD4的泛T面板在CD8较高、DC的共享CTSS在Monocyte/Macrophage较高；这两种排序在训练参考中也存在，不能据此直接判定新出现的错配。预定marker中缺失{', '.join(missing_markers)}，已报告而未补造。该初筛与模型参考相关，不能代替独立recovery或证明无跨谱系串扰。

另需独立复核Tuft：全队列最终RNA fraction中位数约{type_qc.loc[type_qc.cell_type.eq('Tuft'), 'median_final_fraction'].iloc[0] * 100:.2f}%。这不是组织中Tuft细胞数比例的估计真值；marker方向正确也不能验证此丰度。当前只列为生物学解释和外部验证事项，不为改善结果调整参考、重新分配fraction或删类型。低mass类型仍可出现大量可变基因，不能据此抵消信息量不足。

## EcoTyper当前门槛与准备

用户本轮确认没有已冻结的正式主线配置。服务器仍仅有旧MSCCR_CORE12配置；完整2025例匹配HiRes未在已审阅路径发现。原MATCHED_ECOTYPER_DISCOVERY_CONFIG.yml保持HOLD。必须先共同冻结：NMF父类型、Region3决定、基因选择和表达处理、rank/restarts/seed、状态QC、ecotype策略及实际源码。

当前EcoTyper HEAD仍为4a1c71498bca1802a9b374a3259bc7d3f0836fc8，工作树仍包含实验性Region3修改，入口和feature extraction SHA256与昨日一致。没有默认沿用该修改，也没有改动软件源码。已读源码确认presorted入口按每类型TXT读取gene × sample矩阵；Step2同时运行预处理和NMF，不能在设置未定时仅为检查接口就启动Step2。

候选接口清单已保存ECOTYPER_CANDIDATE_INPUT_MANIFEST.tsv。冻结父类型后可运行scripts/export_frozen_ecotyper_inputs.R，把选定候选RDS导出为每类型TXT及同序annotation；本轮不制造全20类正式TXT目录，避免被误作已冻结discovery输入。该导出脚本不会运行EcoTyper。尚缺的父类型冻结清单和设计信息必须先提供，当前无正式NMF任务。

最终状态：**BayesPrism完成；数值QC通过；低mass/生物学特异性需复核；EcoTyper正式NMF等待主线设计冻结。** 旧版53 states/7 ecotypes、Joint20 benchmark及CIBERSORTx主线均未修改。
'''
(out / 'COMPLETION_AND_ECOTYPER_READINESS.md').write_text(report, encoding='utf-8')
print(json.dumps(status, ensure_ascii=False))
print(type_qc[['cell_type', 'mass_median', 'n_mass_lt100', 'n_mass_lt1000']].to_string(index=False))
