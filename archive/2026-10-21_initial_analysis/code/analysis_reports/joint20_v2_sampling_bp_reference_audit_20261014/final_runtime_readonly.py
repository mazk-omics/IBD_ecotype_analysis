# 经现有SSH会话只读确认细胞metadata与进程累计CPU；仅打印快照。
from pathlib import Path
from datetime import datetime
import json
import base64
import hashlib
import pandas as pd

# Step1. 对齐现有参考身份，category先转string再处理缺失
root = Path('/home/mazekai/IBD_EcoTyper')
run = root / '06_validation/joint20_bayesprism_ecotyper_matched_sensitivity_20261014'
cells = pd.read_csv(run / '01_reference/reference_cells.tsv', sep='\t')
fields = ['Patient', 'Disease', 'Treatment', 'Inflammation', 'Site', 'sample_id', 'coarse_identity', 'final_analysis']
pool = pd.read_parquet(root / '03_reference/GSE282122/audit_v2/output/04_identity_coverage_v1/mapped_candidate_cells_v1.parquet', columns=['cell_id'] + fields)
assert cells[fields].astype('string').fillna('NA').to_numpy().tolist() == pool.set_index('cell_id').loc[cells.cell_id, fields].astype('string').fillna('NA').to_numpy().tolist()

# Step2. 读取进程和旧结果stat，不发送任何进程控制指令
workers = {}
for pid in [99113, 99125, 99137, 99149, 99161, 99173, 99185, 99197]:
    values = Path(f'/proc/{pid}/stat').read_text().split(') ', 1)[1].split()
    workers[str(pid)] = {'state': values[0], 'utime_ticks': int(values[11]), 'stime_ticks': int(values[12])}
old_results = []
for name in ['02_bayesprism_joint20_MSCCR120/BayesPrism_result.rds', '05_existing_pseudobulk_benchmark/BayesPrism_natural28_result.rds']:
    path = root / '06_validation/deconv_fraction_sensitivity_joint20_20261013' / name
    old_results.append({'path': str(path), 'bytes': path.stat().st_size, 'mtime': datetime.fromtimestamp(path.stat().st_mtime).astimezone().isoformat()})
script_checks = {}
for line in (run / 'logs/script_SHA256.txt').read_text().splitlines():
    expected, filename = line.split(None, 1)
    path = Path(filename.strip())
    script_checks[path.name] = hashlib.file_digest(path.open('rb'), 'sha256').hexdigest() == expected
assert all(script_checks.values())

# Step3. 输出最后时间点和可供本地解码的快照
result = {'timestamp': datetime.now().astimezone().isoformat(), 'reference_labels_equal_candidate_metadata': True, 'workers': workers, 'main_PID_exists': Path('/proc/99062').exists(), 'completed': (run / 'logs/completed_at.txt').exists(), 'formal_exit_exists': (run / 'logs/formal_exit_status.txt').exists(), 'current_result_exists': (run / '02_bayesprism/BayesPrism_result.rds').exists(), 'old_result_stats': old_results, 'retained_integerization_bound_fraction': 0.002438101, 'script_hash_checks': script_checks}
payload = base64.b64encode(json.dumps(result).encode()).decode()
print('FINAL_PAYLOAD_BEGIN')
for i in range(0, len(payload), 64):
    print(payload[i:i + 64])
print('FINAL_PAYLOAD_END')
