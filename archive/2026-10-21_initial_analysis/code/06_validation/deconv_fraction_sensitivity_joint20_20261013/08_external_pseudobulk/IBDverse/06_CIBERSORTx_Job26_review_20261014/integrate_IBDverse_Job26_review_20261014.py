# 把已完成的Job26配对审阅登记到当前状态；旧状态先留档，历史计算数值不改。
# 本机运行；新增2026-10-14补充报告，两套外部官方配对均完成，不访问服务器。
from pathlib import Path
import json
import shutil
import hashlib
import pandas as pd

# Step1. 留存更新前状态，避免把历史BLOCKED记录抹掉
project = Path(__file__).resolve().parents[2]
root = project / '06_validation/deconv_fraction_sensitivity_joint20_20261013'
dataset = root / '08_external_pseudobulk/IBDverse'
review = dataset / '06_CIBERSORTx_Job26_review_20261014'
archive = review / 'prior_status_before_Job26'
archive.mkdir(exist_ok=False)
paths = [root / 'execution_status.json', root / '07_reports/BLOCKED_AND_PARTIAL_ITEMS.tsv', root / '07_reports/EXTERNAL_REFERENCE_INDEPENDENCE_AND_COVERAGE.tsv', dataset / 'execution_status.json', dataset / '04_accuracy/accuracy_summary.json']
for path in paths:
    shutil.copyfile(path, archive / (path.parent.name + '__' + path.name))

# Step2. 登记IBDverse完成状态；运行参数由用户确认，完整日志仍未提供
blocked = pd.read_csv(root / '07_reports/BLOCKED_AND_PARTIAL_ITEMS.tsv', sep='\t')
blocked = blocked.loc[~blocked.module.eq('IBDverse_official_CIBERSORTx_paired')].copy()
blocked.loc[len(blocked) + 1] = {'module': 'Job26_full_runtime_provenance', 'status': 'PARTIAL', 'reason': 'User confirmed signature, relative, QN off, no batch correction, 100 permutations; full official runtime log, web build and upload hashes unavailable', 'prepared_inputs': ''}
blocked.to_csv(root / '07_reports/BLOCKED_AND_PARTIAL_ITEMS.tsv', sep='\t', index=False)
status = json.loads((root / 'execution_status.json').read_text())
status.update({'IBDverse_official_CIBERSORTx_paired': 'COMPLETE', 'IBDverse_official_job': 26, 'external_official_CIBERSORTx_paired': 'COMPLETE', 'official_external_jobs_received': [25, 26], 'new_CIBERSORTx_run': True, 'new_official_run_operator': 'user; agent did not submit a new job', 'Job26_parameter_confirmation': 'user confirmed in chat 2026-10-14', 'Job26_Smode_sensitivity': 'NOT_RUN', 'blocked_items': blocked.fillna('').to_dict('records'), 'latest_review_date': '2026-10-14', 'server_reports_synchronized': False, 'notes': 'GSE214695 Job25 and IBDverse Job26 official paired analyses complete locally. Both no batch correction, parameters user-confirmed. Server retains 2026-10-13 snapshot; no new SSH connection.'})
(root / 'execution_status.json').write_text(json.dumps(status, indent=2), encoding='utf-8')
for path in [dataset / 'execution_status.json', dataset / '04_accuracy/accuracy_summary.json']:
    state = json.loads(path.read_text())
    state['official_CIBERSORTx_paired'] = 'COMPLETE'
    state['official_job'] = 26
    state['official_reason'] = 'User supplied and confirmed Job26; exact 100 SampleIDs and 20 classes matched; local paired evaluation complete'
    path.write_text(json.dumps(state, indent=2), encoding='utf-8')
coverage = pd.read_csv(root / '07_reports/EXTERNAL_REFERENCE_INDEPENDENCE_AND_COVERAGE.tsv', sep='\t')
coverage.loc[coverage.dataset.eq('IBDverse'), 'official_pair'] = 'COMPLETE: Job26'
coverage.to_csv(root / '07_reports/EXTERNAL_REFERENCE_INDEPENDENCE_AND_COVERAGE.tsv', sep='\t', index=False)

# Step3. 历史报告追加有日期的更新，保留原定终点和实际排序反转
report = review / 'IBDverse_JOB26_PAIRED_REPORT.md'
shutil.copyfile(report, root / '07_reports' / report.name)
addendum = '\n## 2026-10-14 Job26更新\n\nIBDverse官方配对已完成，先前BLOCKED为2026-10-13历史状态。16类可评价轴上，BP/CX的MAE：保留基因RNA真值3.294/2.247 pp；cell-count真值4.730/4.358 pp；完整基因RNA真值3.685/1.920 pp。CX总体误差在三种尺度均低，但主要RNA尺度仍有11/16类是BP误差更低；CX优势主要来自Tuft和Stem/TA。详细报告为IBDverse_JOB26_PAIRED_REPORT.md。GSE214695 Job25及IBDverse Job26官方配对均已完成。原Job26无校正结果保留；S-mode尚未运行，若追加只能先标记探索性敏感性分析，因为Job13的212个汇总profiles与现有1060个原训练单细胞子集不是逐列相同的构建输入。\n'
for path in [root / '07_reports/FINAL_SENSITIVITY_REPORT.md', root / '07_reports/IBDverse_EXTERNAL_BENCHMARK_REPORT.md', dataset / 'IBDverse_EXTERNAL_BENCHMARK_REPORT.md']:
    with path.open('a', encoding='utf-8') as handle:
        handle.write(addendum)
shutil.copyfile(Path(__file__), review / Path(__file__).name)
# Step4. 两套数据独立汇总，不能把18类和16类分母直接合并
joint = []
for name, job in [('GSE214695', 25), ('IBDverse', 26)]:
    folder = root / '08_external_pseudobulk' / name / f'06_CIBERSORTx_Job{job}_review_20261014'
    table = pd.read_csv(folder / 'paired_accuracy_metrics_all_scales.tsv', sep='\t')
    table = table.loc[table.cell_type.eq('POOLED_EVALUABLE')].copy()
    table.insert(0, 'dataset', name)
    joint.append(table)
joint = pd.concat(joint, ignore_index=True)
joint.to_csv(root / '07_reports/EXTERNAL_OFFICIAL_PAIRED_SUMMARY_20261014.tsv', sep='\t', index=False)
lines = ['# 两套外部官方配对验证：2026-10-14更新', '', 'GSE214695 Job25与IBDverse Job26均已完成配对。用户分别确认原GJ20v2 signature、Relative、QN关闭、batch correction关闭、100 permutations。各100份×1000次细胞抽样；可评价类型分别18、16类。Full20推断保留，不重新归一化，不剔除预测差的样本。', '', '| 数据集 | 真值尺度 | 方法 | MAE pp | RMSE pp | pooled Pearson | pooled Spearman |', '|---|---|---|---:|---:|---:|---:|']
for row in joint.loc[~joint.method.eq('BayesPrism_initial')].itertuples():
    lines.append(f'| {row.dataset} | {row.truth_scale} | {row.method} | {row.MAE_pp:.3f} | {row.RMSE_pp:.3f} | {row.Pearson_r:.3f} | {row.Spearman_rho:.3f} |')
lines += ['', 'GSE214695：BP在保留基因RNA与cell-count尺度总体误差较低，CX在完整基因RNA尺度较低。IBDverse：CX三种尺度总体误差均较低；主要RNA尺度BP仍有11/16类MAE较低，CX改善主要集中于Tuft与Stem/TA，而ILC全部预测为0（49份有真值）。不据总体MAE把任一方法宣称为所有谱系均更好。', '', '保留基因RNA是预定主要终点，对BP估计量有方法特异性；两个次要尺度一并报告，不根据结果切换主要终点。18/16类覆盖与背景不同，两个数据集的误差不直接合并。100个混合物复用细胞和供者，不能充当100个独立供者。跨研究患者身份核验及完整官方运行日志限制仍为PARTIAL。', '', '详细报告：GSE214695_JOB25_PAIRED_REPORT.md、IBDverse_JOB26_PAIRED_REPORT.md。此汇总为最新本地结果，2026-10-13压缩包和服务器报告属于历史快照；本次无新SSH、无新signature或EcoTyper discovery。']
(root / '07_reports/EXTERNAL_OFFICIAL_PAIRED_SUMMARY_20261014.md').write_text('\n'.join(lines) + '\n', encoding='utf-8')
records = []
for path in sorted(root.rglob('*')):
    if path.is_file() and path.name != 'local_output_manifest.tsv':
        records.append({'relative_path': str(path.relative_to(root)), 'bytes': path.stat().st_size, 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()})
pd.DataFrame(records).to_csv(root / 'local_output_manifest.tsv', sep='\t', index=False)
print('Job25 and Job26 paired COMPLETE; old status archived; no SSH')
