# 整理 2026-10-07 通过一个持久 SSH 会话读取的 recovery 历史记录，不重新运行服务器分析。
# 在本机项目目录运行；读取已回传的 JSON 证据，保存可阅读的服务器小文件副本和输入校验汇总。
# 原服务器脚本仅作为历史证据保留，不按新代码规范改写，也不下载大型表达矩阵。
from pathlib import Path
import csv
import io
import json

# Step1. 读入已回传证据，保留此次检查范围。
project = Path(__file__).resolve().parents[3]
outdir = project / 'IBD_ecotype_analysis/06_validation/recovery_inventory_20261007'
evidence = json.loads((outdir / 'source_evidence.json').read_text(encoding='utf-8'))
software = json.loads((outdir / 'software_details.json').read_text(encoding='utf-8'))
print('历史小文件数量:', len(evidence['files']))
print('已复核表达输入:', list(evidence['matrix_checks']))

# Step2. 保存原有脚本、日志与小型审计表，路径保持与服务器项目一致。
snapshot = outdir / 'server_snapshot'
for relative, content in evidence['files'].items():
    destination = snapshot / relative
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(content, encoding='utf-8')

for absolute, content in software['discovery_models'].items():
    relative = Path(absolute).as_posix().split('/IBD_EcoTyper/', 1)[1]
    destination = snapshot / relative
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(content, encoding='utf-8')

# Step3. 比较服务器现有输入与旧校验和；这里只复核已读取证据，不再次计算表达数据。
md5_records = {}
for relative, content in evidence['files'].items():
    if '06_ecotyper_recovery_input' not in relative or 'md5' not in relative:
        continue
    rows = list(csv.DictReader(io.StringIO(content), delimiter='\t'))
    for row in rows:
        md5_records[row['File']] = row['MD5']

checks = evidence['matrix_checks']
for cohort, check in checks.items():
    assert check['genes'] == check['unique_genes']
    assert check['samples'] == check['unique_samples'] == check['annotation_rows']
    assert check['exact_order_match_first_column']
    assert check['row_widths'] == [check['samples'] + 1]
    assert check['nonfinite_values'] == check['invalid_numeric_entries'] == 0
    assert check['expression_md5'] == md5_records[check['expression_path']]
    assert check['annotation_md5'] == md5_records[check['annotation_path']]
    check['matches_original_expression_md5'] = True
    check['matches_original_annotation_md5'] = True
    print(cohort, 'genes=', check['genes'], 'samples=', check['samples'], '顺序及旧 MD5 均一致')

# Step4. 保存本次复核结论；输入准备完成不等于已运行状态或生态型 recovery。
summary = {'checked_date': '2026-10-07', 'server_root': evidence['root'], 'matrix_checks': checks, 'project_recovery_directories': evidence['recovery_directories'], 'internal_recovery_directories': software['internal_recovery_directories'], 'tracked_example_recovery_files': len(evidence['tracked_software_recovery_files']), 'ibd_recovery_results_found_in_checked_paths': False, 'signature_gene_compatibility_with_final_discovery_checked': False, 'remote_analysis_rerun': False, 'ssh_sessions_created': 1, 'ssh_session_closed': True}
(outdir / 'input_recheck_summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
print('历史证据副本:', snapshot)
print('复核汇总:', outdir / 'input_recheck_summary.json')
