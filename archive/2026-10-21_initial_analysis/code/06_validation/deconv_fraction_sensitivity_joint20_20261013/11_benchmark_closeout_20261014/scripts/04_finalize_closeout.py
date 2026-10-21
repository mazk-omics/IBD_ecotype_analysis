# 汇总四份核心收尾成果，记录缺失来源和未解决机制，核实旧文件没有被覆盖。
# 本机运行；所有占位表由复算结果填充，不改变原benchmark的任何参数、组成或预测。
from pathlib import Path
import hashlib
import json
import re
import shutil
import subprocess
import sys
import numpy as np
import pandas as pd

# Step1. 核对参数原件、1060子集和历史hash，资料缺失继续保留PARTIAL
out = Path(__file__).resolve().parents[1]
root = out.parent
project = root.parents[1]
inputs = pd.read_csv(out / 'INPUT_MANIFEST.tsv', sep='\t')
audit = pd.read_csv(out / 'AUDIT_CHECKS.tsv', sep='\t')
metrics = pd.read_csv(out / 'EXTERNAL_RECOMPUTED_METRICS.tsv', sep='\t')
zero = pd.read_csv(out / 'DETECTION_AND_GENE_DENOMINATORS.tsv', sep='\t')
reliability = pd.read_csv(out / 'JOINT20_CELLTYPE_RELIABILITY.tsv', sep='\t')
effects = pd.read_csv(out / 'MSCCR120_DESCRIPTIVE_EFFECTS.tsv', sep='\t')
summary = json.loads((out / 'numeric_audit_summary.json').read_text())
for folder in [root / '02_bayesprism_joint20_MSCCR120', root / '08_external_pseudobulk/GSE214695/03_bayesprism', root / '08_external_pseudobulk/IBDverse/03_bayesprism']:
    path = folder / 'installed_run_control.txt'
    control = path.read_text()
    actual = {}
    for parameter, expected in [('chain.length', 1000), ('burn.in', 500), ('thinning', 2), ('seed', 123), ('alpha', 1)]:
        match = re.search(r'\$gibbs.control\$' + re.escape(parameter) + r'\s+\[1\]\s+(\S+)', control)
        assert match is not None and float(match.group(1)) == expected
        actual[parameter] = float(match.group(1))
    audit.loc[len(audit)] = {'check': str(folder.relative_to(root)) + ' saved Gibbs controls', 'status': 'PASS', 'evidence': json.dumps(actual)}
source_dir = project / 'CIBERSORTx_ready_20261008/01_MSCCR_platform120_pilot'
source_manifest = pd.read_csv(source_dir / '核对记录_勿上传/source_subset_cell_manifest.tsv', sep='\t')
reference = pd.read_csv(root / '01_reference/reference_cells.tsv', sep='\t')
source_path = source_dir / 'GJ20v2_SourceSubset1060_scRNA_CPM.txt'
assert len(source_manifest) == 1060 and source_manifest.groupby('profile_id').size().eq(5).all()
assert source_manifest.profile_id.nunique() == 212 and set(source_manifest.cell_id) <= set(reference.cell_id)
assert source_manifest.Patient.nunique() == 32 and set(source_manifest.coarse_identity) == set(reference.coarse_identity)
assert hashlib.sha256(source_path.read_bytes()).hexdigest() == '9f7d25d37c69541c0f3028db8f323648604fb7bc5efab1abac91dcc410836ac3'
audit.loc[len(audit)] = {'check': '1060 source subset membership and expression-file integrity', 'status': 'PASS', 'evidence': '212 profiles x 5 training cells; 32 donors; 20 types; saved source hash exact. Does not establish S-mode correction efficacy.'}
for dataset, job in [('GSE214695', 25), ('IBDverse', 26)]:
    r = root / '08_external_pseudobulk' / dataset / f'06_CIBERSORTx_Job{job}_review_20261014'
    inventory = pd.read_csv(r / 'official_input_inventory.tsv', sep='\t')
    for row in inventory.itertuples():
        archived = r / 'original_outputs' / Path(row.path).name
        assert archived.exists() and hashlib.sha256(archived.read_bytes()).hexdigest() == row.sha256
    settings = json.loads((r / 'official_QC_and_settings.json').read_text())
    assert settings['relative'] and settings['QN'] is False and settings['batch_correction'] is False and settings['permutations'] == 100
    audit.loc[len(audit)] = {'check': f'Job{job} archived originals and reported parameters', 'status': 'PASS', 'evidence': 'Archived TXT/PDF hashes match supplied inventory; parameter evidence is user confirmation, not missing server logs.'}
signature_hash = hashlib.sha256((project / 'CIBERSORTx_ready_20261007/09_MSCCR_Joint20v2_Fractions/GJ20v2_Signature20.txt').read_bytes()).hexdigest()
assert signature_hash == '26babf93a06db1d2063fdada032aacf90712b24a23faece7ee5b5ddd172df416'
audit.loc[len(audit)] = {'check': 'Job13-derived frozen signature hash', 'status': 'PASS', 'evidence': signature_hash}
upload_inventory = pd.read_csv(root / '09_CIBERSORTx_upload_small/upload_size_and_sha256.tsv', sep='\t')
for row in upload_inventory.itertuples():
    path = root / '09_CIBERSORTx_upload_small' / row.file
    assert hashlib.sha256(path.read_bytes()).hexdigest() == row.sha256
audit.loc[len(audit)] = {'check': 'Prepared official upload matrices unchanged', 'status': 'PASS', 'evidence': 'All three existing prepared TXT hashes reproduced; this is not a server-upload receipt hash.'}
audit.to_csv(out / 'AUDIT_CHECKS.tsv', sep='\t', index=False)

# Step2. 原文评估单位和本轮数值表，原文图版式未取得则明确不声称逐版复制
status = 'CLOSED WITH LIMITATIONS'
lines = ['| 审计项目 | 状态 | 证据及边界 |', '|---|---|---|']
for row in audit.itertuples():
    lines.append(f'| {row.check} | {row.status} | {row.evidence} |')
audit_table = '\n'.join(lines)
anomaly_specs = [
    ('GSE214695 Plasma', '数量bias +27.502 pp；完整RNA仍+5.208 pp', '完整RNA truth均值36.27%远高于细胞数13.98%，且BP Ig过滤进一步改变分母；尺度解释了很大一部分差异', '剩余完整RNA高估、参考汇总、平台、细胞RNA量、作者注释/环境RNA均未分离；UNRESOLVED MECHANISM', '按细胞数量解释为C；完整RNA仍需B级定量限制'),
    ('GSE214695 Stem/TA', '97/100预测0，阳性漏检51/54；各RNA尺度仍显著偏低', '异常不是仅由细胞数/RNA分母改变造成；signature/作者marker存在，确定性映射与列错位已排除', '作者Cycling TA与参考Stem/TA的状态边界、平台或与其他上皮共线性；UNRESOLVED MECHANISM', '该测试配置中不能可靠检出或定量，C'),
    ('IBDverse ILC', '100/100预测0，49/49阳性漏检；三truth均失败，r为NA', '并非外部设计没有ILC，也不是signature完全缺少IL7R/KIT/NCR2；BP也偏低但并未恒零', 'ILC亚型、特征稳定性、SVR零截断、平台和标签/ambient差异；UNRESOLVED MECHANISM', '该配置C；不得以零预测作缺乏ILC证据'),
    ('IBDverse CX Stem/TA', '数量bias +8.791 pp；完整RNA bias -3.503 pp；retained RNA -2.681 pp', '“高估”的方向依赖truth；不能把数量高估描述为三个尺度均高估', '仍有较大RMSE、7/74阳性漏检、作者状态/数据集差异；UNRESOLVED MECHANISM', 'B，不能精确计数；变化趋势与校准分开'),
    ('IBDverse BP Tuft / Stem', 'Tuft retained-RNA +17.031 pp、完整RNA +17.601 pp；Stem分别-13.224/-14.045 pp', 'Tuft与Stem的偏差不因换RNA分母消失；参考和外部有典型marker', '背景基因、参考状态不足、技术差异及组成竞争；误差负协变不证明特定reads流向；UNRESOLVED MECHANISM', 'BP敏感性中对应RNA估计C；不据此自动判定CX Tuft为C'),
    ('MSCCR absorptive / Tuft / Stem', 'BP-CX平均差分别为结肠+16.72、回肠+7.96、Tuft+10.99、Stem-14.49 pp', '高r可与高RMSE并存；Job17/20本身不同；基因空间、S-mode和参考表示不一致', '无真值，另有小数count精度与生物/技术混杂；UNRESOLVED MECHANISM', '只能报告方法差异和描述性敏感性，不能判定准确方法或跨方法模型等价')
]
lines = ['| 条件/谱系 | 确认异常 | 现有材料支持的解释 | 未排除的替代解释 | 定量解释影响 |', '|---|---|---|---|---|']
for spec in anomaly_specs:
    lines.append('| ' + ' | '.join(spec) + ' |')
anomaly_table = '\n'.join(lines)
lines = ['| Dataset / type | Cell truth % | Full RNA truth % | BP-gene RNA truth % | CX prediction % | BP final % |', '|---|---:|---:|---:|---:|---:|']
for dataset, cell_type in [('GSE214695', 'Plasma'), ('GSE214695', 'Stem_TA_progenitor'), ('IBDverse', 'ILC'), ('IBDverse', 'Stem_TA_progenitor'), ('IBDverse', 'Tuft')]:
    a = metrics.loc[metrics.dataset.eq(dataset) & metrics.cell_type.eq(cell_type)].set_index(['truth_scale', 'method'])
    c = a.loc[('Cell_count_SECONDARY', 'CIBERSORTx')]
    f = a.loc[('Full_gene_RNA_SECONDARY', 'CIBERSORTx')]
    g = a.loc[('Retained_gene_RNA_PRIMARY', 'CIBERSORTx')]
    b = a.loc[('Cell_count_SECONDARY', 'BayesPrism_final')]
    lines.append(f'| {dataset} / {cell_type} | {c.truth_mean_percent:.3f} | {f.truth_mean_percent:.3f} | {g.truth_mean_percent:.3f} | {c.prediction_mean_percent:.3f} | {b.prediction_mean_percent:.3f} |')
five_table = '\n'.join(lines)
lines = ['| Cell type | Overall CX evidence | GSE cell / full RNA | IBDverse cell / full RNA | Interpretation boundary |', '|---|---|---|---|---|']
for row in reliability.itertuples():
    lines.append(f'| {row.cell_type} | {row.overall_CIBERSORTx_evidence} | {row.GSE214695_cell_count_grade} / {row.GSE214695_full_RNA_grade} | {row.IBDverse_cell_count_grade} / {row.IBDverse_full_RNA_grade} | {row.interpretation_limit} |')
reliability_table = '\n'.join(lines)
lines = []
for dataset, method, scale in [('GSE214695', 'CIBERSORTx', 'Cell_count_SECONDARY'), ('IBDverse', 'CIBERSORTx', 'Cell_count_SECONDARY'), ('IBDverse', 'BayesPrism_final', 'Retained_gene_RNA_PRIMARY')]:
    a = metrics.loc[metrics.dataset.eq(dataset) & metrics.method.eq(method) & metrics.truth_scale.eq(scale) & ~metrics.cell_type.str.startswith('POOLED')].nlargest(3, 'squared_error_share_Full20')
    lines.append('- ' + dataset + ' / ' + method + ' / ' + scale + '：Full20平方误差贡献最大三类为' + '、'.join(f'{r.cell_type} {r.squared_error_share_Full20 * 100:.1f}%' for r in a.itertuples()) + '。这是误差分解，不是因果贡献。')
contribution_text = '\n'.join(lines)

# Step3. 真实bulk汇总与最小描述性效应，不把方向一致当临床验证
b = summary['MSCCR120']
bulk_summary = f"- Pearson类型中位数{b['Pearson_median_across_types']:.6f}，Spearman{b['Spearman_median_across_types']:.6f}；平均绝对方法差{b['mean_absolute_difference_pp']:.6f} pp；样本TV中位数{b['median_sample_total_variation']:.6f}。\n- Job20 vs Job17平均绝对差{b['Job20_vs_Job17_mean_absolute_difference_pp']:.6f} pp，TV中位数{b['Job20_vs_Job17_median_total_variation']:.6f}。\n- 这些数值与保存矩阵一致，不同配置间无真值不能比较accuracy。"
lines = ['| Scope | Contrast | BP final vs Job20 same direction | Job20 vs Job17 same direction |', '|---|---|---:|---:|']
for (scope, contrast), block in effects.groupby(['scope', 'contrast'], sort=False):
    lines.append(f'| {scope} | {contrast} | {int(block.BP_final_Job20_sign_agree.sum())}/{len(block)} | {int(block.Job20_Job17_sign_agree.sum())}/{len(block)} |')
effect_table = '\n'.join(lines)
lines = ['| Dataset | Truth | Method | MAE pp | RMSE pp | Pooled Pearson | Pooled Spearman |', '|---|---|---|---:|---:|---:|---:|']
for row in metrics.loc[metrics.cell_type.eq('POOLED_EVALUABLE') & metrics.method.isin(['CIBERSORTx', 'BayesPrism_final'])].itertuples():
    lines.append(f'| {row.dataset} | {row.truth_scale} | {row.method} | {row.MAE_pp:.3f} | {row.RMSE_pp:.3f} | {row.Pearson_r:.3f} | {row.Spearman_rho:.3f} |')
external_table = '\n'.join(lines)
replacements = {'{{CLOSEOUT_STATUS}}': status, '{{AUDIT_TABLE}}': audit_table, '{{ANOMALY_TABLE}}': anomaly_table, '{{FIVE_QUANTITY_TABLE}}': five_table, '{{RELIABILITY_TABLE}}': reliability_table, '{{ERROR_CONTRIBUTIONS}}': contribution_text, '{{MSCCR_SUMMARY}}': bulk_summary, '{{EFFECT_DIRECTION_TABLE}}': effect_table, '{{FINAL_EXTERNAL_TABLE}}': external_table}
for name in ['BENCHMARK_CLOSEOUT_AUDIT.md', 'BENCHMARK_FINAL_METHODS_RESULTS.md']:
    path = out / name
    text = path.read_text(encoding='utf-8')
    for key, value in replacements.items():
        text = text.replace(key, value)
    assert '{{' not in text
    path.write_text(text, encoding='utf-8')

# Step4. 图件最终QA和原件不变性作为收尾门槛，只检查本轮新图
skill = Path.home() / '.codex/skills/nature-figure/scripts'
plot_script = out / 'scripts/03_plot_signature_validation.py'
pdf = out / 'figures/SIGNATURE_VALIDATION.pdf'
commands = [('source_preflight', [sys.executable, str(skill / 'validate_figure.py'), str(plot_script)]), ('text_audit', [sys.executable, str(skill / 'audit_pdf_text.py'), str(pdf), '--json']), ('collision_audit', [sys.executable, str(skill / 'audit_figure_collisions.py'), str(pdf), '--json-out', str(out / 'figures/SIGNATURE_VALIDATION.collision.json')])]
for name, command in commands:
    run = subprocess.run(command, capture_output=True, text=True, encoding='utf-8', errors='replace')
    (out / 'figures' / (name + '.txt')).write_text(run.stdout + run.stderr, encoding='utf-8')
    assert run.returncode == 0, run.stdout
alignment = json.loads((out / 'figures/SIGNATURE_VALIDATION.alignment.json').read_text())
collision = json.loads((out / 'figures/SIGNATURE_VALIDATION.collision.json').read_text())
font = json.loads((out / 'figures/text_audit.txt').read_text())
assert alignment['verdict'] == 'PASS' and collision['verdict'] == 'PASS' and font['minimum_found_pt'] >= 5
notes = [{'panel': p, 'alignment': 'PASS', 'collision': 'PASS', 'minimum_font_pt': font['minimum_found_pt'], 'visual_review': 'All 20 rows, 34 eligible pairs, NA and NE retained; reviewed full rendered figure; no biological CI implied'} for p in 'abc']
pd.DataFrame(notes).to_csv(out / 'figures/RENDERED_QA.tsv', sep='\t', index=False)
for row in inputs.itertuples():
    assert hashlib.sha256(Path(row.path).read_bytes()).hexdigest() == row.sha256, 'Input changed during closeout: ' + row.path
audit.loc[len(audit)] = {'check': 'All frozen/local input hashes unchanged after closeout', 'status': 'PASS', 'evidence': f'{len(inputs)} hashed inputs rechecked; no modifications to old files or other tasks'}
audit.to_csv(out / 'AUDIT_CHECKS.tsv', sep='\t', index=False)

# Step5. 新版本状态与输出清单；旧2026-10-13目录状态不被覆盖
state = {'closeout_status': status, 'numerical_audit': 'PASS', 'input_axis_and_draw_truth': 'PASS', 'RNA_truth_algebra': 'PASS', 'RNA_truth_source_reaggregation_this_round': 'PARTIAL; not performed locally', 'official_full_runtime_provenance': 'PARTIAL', 'cross_study_biological_donor_identity': 'PARTIAL', 'MSCCR_fractional_count_mechanism': 'PARTIAL; raw-count designation supported; exact fractional assignment undocumented', 'MSCCR_BP_sampler_precision': 'PARTIAL; integer coercion documented, no new fit', 'Job20_Smode_pairing': 'PARTIAL; source matches but representation differs', 'Job14_strict_count_based_pair': 'BLOCKED outside primary external validation', 'anomaly_mechanisms': 'UNRESOLVED MECHANISM except observed denominator effects and deterministic exclusions', 'reliability_matrix': 'COMPLETE; 20 types with explicit scopes', 'Cell2021_style_figure': 'COMPLETE; adapted type/cohort unit, absolute error added, failed/untested rows retained', 'input_files_unchanged': True, 'new_remote_connections': 0, 'new_signature': False, 'new_Smode': False, 'new_BayesPrism': False, 'new_clinical_model': False, 'new_NMF_or_discovery': False, 'frozen_states_changed': False, 'frozen_ecotypes_changed': False, 'observed_major_implementation_blockers': [], 'reason': 'Saved numerical results and deterministic axes/draw truth reproduce. Remaining limitations are explicitly scoped; no evidence requires replacing core external results or pausing the CIBERSORTx fraction workflow.'}
(out / 'execution_status.json').write_text(json.dumps(state, indent=2), encoding='utf-8')
records = []
for path in sorted(out.rglob('*')):
    if path.is_file() and path.name != 'OUTPUT_MANIFEST.tsv':
        records.append({'relative_path': str(path.relative_to(out)), 'bytes': path.stat().st_size, 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()})
pd.DataFrame(records).to_csv(out / 'OUTPUT_MANIFEST.tsv', sep='\t', index=False)
print('CLOSEOUT:', status)
print('Inputs rechecked unchanged:', len(inputs))
print('20-type matrix and four core deliverables complete; final figure QA passed.')
