# 将已核验的小型统计表嵌入审计报告，解码最后的只读运行快照。
# 只写本次本地报告目录，不接触正式运行目录或服务器。
from pathlib import Path
import base64
import json
import pandas as pd

# Step1. 解码已有SSH输出中的最终快照
out = Path('D:/IBD_ecotype/IBD_ecotype_analysis/analysis_reports/joint20_v2_sampling_bp_reference_audit_20261014')
snapshot = json.loads(base64.b64decode((out / 'final_runtime_snapshot.b64').read_text()))
(out / 'final_runtime_snapshot.json').write_text(json.dumps(snapshot, indent=2) + '\n', encoding='utf-8')
previous = json.loads((out / 'remote_readonly_audit.json').read_text(encoding='utf-8'))
assert snapshot['main_PID_exists'] and not snapshot['formal_exit_exists']
assert all(snapshot['workers'][pid]['utime_ticks'] > previous['proc_ticks'][pid]['utime_ticks'] for pid in snapshot['workers'])

# Step2. 格式化疾病、炎症及肠段的实际覆盖
overall_lines = ['|维度/类别|细胞|供者|样本|profiles|', '|---|---:|---:|---:|---:|']
for field in ['Disease', 'Inflammation', 'Site']:
    table = pd.read_csv(out / f'overall_{field}.tsv', sep='\t')
    for row in table.itertuples():
        overall_lines.append(f'|{field}: {row.audit_value}|{row.cells}|{row.donors}|{row.samples}|{row.profiles}|')
summary = pd.read_csv(out / 'type_summary.tsv', sep='\t')
type_lines = ['|父类型|细胞|供者|样本/profiles|I/CR|CD/UC/H|Infl/NInfl|最大供者%|', '|---|---:|---:|---:|---:|---:|---:|---:|']
for row in summary.itertuples():
    type_lines.append(f'|{row.coarse_identity}|{row.cells}|{row.donors}|{row.samples}/{row.profiles}|{row.region_Ileum}/{row.region_ColonRectum}|{row.Disease_CD}/{row.Disease_UC}/{row.Disease_Healthy}|{row.Inflammation_Inflamed}/{row.Inflammation_Non_Inflamed}|{row.max_donor_fraction * 100:.1f}|')

# Step3. 完成报告并校验关键总量与输入完整性
path = out / 'JOINT20_V2_SAMPLING_AND_BP_REFERENCE_AUDIT.md'
report = path.read_text(encoding='utf-8')
report = report.replace('{{OVERALL_TABLE}}', '\n'.join(overall_lines)).replace('{{TYPE_TABLE}}', '\n'.join(type_lines))
report = report.replace('审计日期：2026-10-14，Asia/Shanghai。', f'审计日期：2026-10-14，Asia/Shanghai；最终进程快照：{snapshot["timestamp"]}。')
assert '{{' not in report and summary.cells.sum() == 9294 and summary.profiles.sum() == 212
assert pd.read_csv(out / 'live_input_hash_verification.tsv', sep='\t').match.all()
assert pd.read_csv(out / 'sampleID_match_2025.tsv', sep='\t').exact_ID_match.all()
path.write_text(report, encoding='utf-8')
print('REPORT_READY', path)
print('FINAL_RUNTIME', snapshot['timestamp'], 'main_PID 99062', '8 workers CPU ticks increased')
print('INPUT_HASHES_PASS', len(pd.read_csv(out / 'live_input_hash_verification.tsv', sep='\t')))
