# 把已完成的Job25配对审阅登记到当前状态；旧状态先留档，历史计算数值不改。
# 本机运行；新增2026-10-14补充报告，IBDverse继续等待官方结果，不访问服务器。
from pathlib import Path
import json
import shutil
import hashlib
import pandas as pd

# Step1. 留存更新前状态，避免把历史BLOCKED记录抹掉
project = Path(__file__).resolve().parents[2]
root = project / '06_validation/deconv_fraction_sensitivity_joint20_20261013'
dataset = root / '08_external_pseudobulk/GSE214695'
review = dataset / '06_CIBERSORTx_Job25_review_20261014'
archive = review / 'prior_status_20261013'
archive.mkdir(exist_ok=False)
paths = [root / 'execution_status.json', root / '07_reports/BLOCKED_AND_PARTIAL_ITEMS.tsv', root / '07_reports/EXTERNAL_REFERENCE_INDEPENDENCE_AND_COVERAGE.tsv', dataset / 'execution_status.json', dataset / '04_accuracy/accuracy_summary.json']
for path in paths:
    shutil.copyfile(path, archive / (path.parent.name + '__' + path.name))

# Step2. 完成状态仅针对GSE214695；运行参数由用户确认，完整日志仍未提供
blocked = pd.read_csv(root / '07_reports/BLOCKED_AND_PARTIAL_ITEMS.tsv', sep='\t')
blocked = blocked.loc[~blocked.module.eq('GSE214695_official_CIBERSORTx_paired')].copy()
blocked.loc[len(blocked) + 1] = {'module': 'Job25_full_runtime_provenance', 'status': 'PARTIAL', 'reason': 'User confirmed signature, relative, QN off, no batch correction, 100 permutations; full official runtime log, web build and upload hashes unavailable', 'prepared_inputs': ''}
blocked.to_csv(root / '07_reports/BLOCKED_AND_PARTIAL_ITEMS.tsv', sep='\t', index=False)
status = json.loads((root / 'execution_status.json').read_text())
status.update({'GSE214695_official_CIBERSORTx_paired': 'COMPLETE', 'GSE214695_official_job': 25, 'IBDverse_official_CIBERSORTx_paired': 'BLOCKED', 'external_official_CIBERSORTx_paired': 'PARTIAL', 'official_external_jobs_received': [25], 'new_CIBERSORTx_run': True, 'new_official_run_operator': 'user; agent did not submit a new job', 'Job25_parameter_confirmation': 'user confirmed in chat 2026-10-14', 'Job25_Smode_sensitivity': 'NOT_RUN', 'blocked_items': blocked.fillna('').to_dict('records'), 'latest_review_date': '2026-10-14', 'server_reports_synchronized': False, 'notes': 'Job25 paired comparison completed locally; IBDverse official output pending. Server retains 2026-10-13 snapshot; no new SSH connection.'})
(root / 'execution_status.json').write_text(json.dumps(status, indent=2), encoding='utf-8')
for path in [dataset / 'execution_status.json', dataset / '04_accuracy/accuracy_summary.json']:
    state = json.loads(path.read_text())
    state['official_CIBERSORTx_paired'] = 'COMPLETE'
    state['official_job'] = 25
    state['official_reason'] = 'User supplied and confirmed Job25; exact 100 SampleIDs and 20 classes matched; local paired evaluation complete'
    path.write_text(json.dumps(state, indent=2), encoding='utf-8')
coverage = pd.read_csv(root / '07_reports/EXTERNAL_REFERENCE_INDEPENDENCE_AND_COVERAGE.tsv', sep='\t')
coverage.loc[coverage.dataset.eq('GSE214695'), 'official_pair'] = 'COMPLETE: Job25'
coverage.to_csv(root / '07_reports/EXTERNAL_REFERENCE_INDEPENDENCE_AND_COVERAGE.tsv', sep='\t', index=False)

# Step3. 历史报告追加有日期的更新，保留原定终点和实际排序反转
report = review / 'GSE214695_JOB25_PAIRED_REPORT.md'
shutil.copyfile(report, root / '07_reports' / report.name)
addendum = '\n## 2026-10-14 Job25更新\n\nGSE214695官方配对已完成，先前BLOCKED为2026-10-13历史状态。18类可评价轴上，BP/CX的MAE：保留基因RNA真值1.340/3.273 pp；cell-count真值2.434/4.140 pp；完整基因RNA真值2.577/1.499 pp。比较排序随基因分母改变，不能宣称一种方法全面更好。详细报告为GSE214695_JOB25_PAIRED_REPORT.md。IBDverse仍等待官方输出。原Job25无校正结果保留；S-mode尚未运行，若追加只能先标记探索性敏感性分析，因为Job13的212个汇总profiles与现有1060个原训练单细胞子集不是逐列相同的构建输入。\n'
for path in [root / '07_reports/FINAL_SENSITIVITY_REPORT.md', root / '07_reports/GSE214695_EXTERNAL_BENCHMARK_REPORT.md', dataset / 'GSE214695_EXTERNAL_BENCHMARK_REPORT.md']:
    with path.open('a', encoding='utf-8') as handle:
        handle.write(addendum)
shutil.copyfile(Path(__file__), review / Path(__file__).name)
records = []
for path in sorted(root.rglob('*')):
    if path.is_file() and path.name != 'local_output_manifest.tsv':
        records.append({'relative_path': str(path.relative_to(root)), 'bytes': path.stat().st_size, 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()})
pd.DataFrame(records).to_csv(root / 'local_output_manifest.tsv', sep='\t', index=False)
print('Job25 COMPLETE; IBDverse BLOCKED; old status archived; no SSH')
