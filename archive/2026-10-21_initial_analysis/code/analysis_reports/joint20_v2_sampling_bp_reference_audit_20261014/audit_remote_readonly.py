# 在已建立的单一SSH会话以stdin执行，只读实际输入与metadata，结果打印回本地。
# 不写服务器文件、不读取全量单细胞counts、不调用BayesPrism计算或管理进程。
from pathlib import Path
import base64
import gzip
import hashlib
import json
import subprocess
from datetime import datetime
import numpy as np
import pandas as pd
from scipy import sparse

# Step1. 核验冻结输入哈希、实际细胞及mixture的SampleID
root = Path('/home/mazekai/IBD_EcoTyper')
run = root / '06_validation/joint20_bayesprism_ecotyper_matched_sensitivity_20261014'
training = root / '03_reference/joint20_donor_site_v2/output/run_20261007_v2'
benchmark = root / '06_validation/deconv_fraction_sensitivity_joint20_20261013'
cells = pd.read_csv(training / 'selected_cell_manifest.tsv', sep='\t')
actual = pd.read_csv(run / '01_reference/reference_cells.tsv', sep='\t')
pd.testing.assert_frame_equal(cells, actual)
annotation = pd.read_csv(run / '00_audit/MSCCR_CD_UC_discovery_manifest.tsv', sep='\t')
main = pd.read_csv(root / '03_reference/web_reference_three_pilots/msccr_ibd_joint_web_20261006/MSCCR_CD_UC_2025_annotation.tsv', sep='\t')
legacy = pd.read_csv(root / '05_EcoTyper/input/MSCCR_CORE12/annotation/MSCCR_CORE12_annotation.tsv', sep='\t')
with gzip.open(run / '02_bayesprism/mixture_raw_counts.tsv.gz', 'rt') as handle:
    header = handle.readline().rstrip('\n').split('\t')[1:]
assert header == annotation.ID.tolist() == main.ID.tolist()
assert set(legacy.loc[legacy.Disease.isin(['CD', 'UC']), 'ID']) == set(header)
assert annotation.ID.is_unique and annotation.Disease.isin(['CD', 'UC']).all()
hash_records = []
for row in pd.read_csv(run / '00_audit/input_SHA256.tsv', sep='\t').itertuples():
    path = Path(row.path)
    digest = hashlib.file_digest(path.open('rb'), 'sha256').hexdigest()
    hash_records.append({'path': str(path), 'sha256': digest, 'frozen_sha256': row.sha256, 'match': digest == row.sha256})
    assert digest == row.sha256, str(path)
for name in ['clean_reference.rds', 'gene_cleanup_audit.tsv', 'reference_cells.tsv', 'reference_genes.tsv']:
    digest = hashlib.file_digest((run / '01_reference' / name).open('rb'), 'sha256').hexdigest()
    old = hashlib.file_digest((benchmark / '01_reference' / name).open('rb'), 'sha256').hexdigest()
    assert digest == old, name
    hash_records.append({'path': str(run / '01_reference' / name), 'sha256': digest, 'frozen_sha256': old, 'match': True})
counts = sparse.load_npz(training / 'selected_training_counts_csr.npz')
assert counts.shape[0] == 9294 and np.isfinite(counts.data).all() and (counts.data >= 0).all()
assert np.array_equal(counts.data, np.rint(counts.data))
assert np.array_equal(counts.sum(axis=1).A1, cells.extracted_total_counts)
del counts

# Step2. 对metadata执行原始抽样算法，全部取消导出以防修改
construction = json.loads((training / 'construction_summary.json').read_text())
candidate = pd.read_parquet(construction['candidate_file'])
roles = pd.read_csv(training / 'frozen_donor_roles.tsv', sep='\t')
candidate['split'] = candidate.Patient.map(roles.set_index('Patient').split)
candidate['region'] = np.where(candidate.Site.eq('Terminal_Ileum'), 'Ileum', 'ColonRectum')
types = construction['types']
action = candidate.primary_reference_action.eq('KEEP') | candidate.coarse_identity.isin(['Glial', 'Tuft'])
identity = action & candidate.identity_status.eq('CONFIDENT') & candidate.coarse_identity.isin(types)
training_mask = identity & candidate.split.eq('training')
doublet_mask = training_mask & ~candidate.predicted_doublets
valid = doublet_mask & np.isfinite(candidate.total_counts) & candidate.total_counts.gt(0) & np.isfinite(candidate.n_genes_by_counts) & candidate.n_genes_by_counts.gt(0)
pool = candidate.loc[valid].copy().sort_values('cell_id')
script = (training / 'build_joint20_donor_site_reference_v2.py').read_text()
segment = script.split('# Step2.')[1].split('# Step4.')[0]
segment = '# Step2.' + segment
segment = '\n'.join(line for line in segment.splitlines() if '.to_csv(' not in line and not line.startswith('print('))
namespace = {'pool': pool, 'types': types, 'pd': pd, 'np': np, 'donor_roles': roles}
exec(segment, namespace)
replayed = namespace['selected']
assert set(replayed.cell_id) == set(cells.cell_id)
replayed_by_id = replayed.set_index('cell_id').reindex(cells.cell_id)
assert np.array_equal(replayed_by_id.profile_id.to_numpy(), cells.profile_id.to_numpy())
units = namespace['units']
assert len(namespace['profiles']) == 212
candidate_hashes = []
for path, frozen in construction['input_sha256'].items():
    digest = hashlib.file_digest(Path(path).open('rb'), 'sha256').hexdigest()
    candidate_hashes.append({'path': path, 'sha256': digest, 'frozen_sha256': frozen, 'match': digest == frozen})
    assert digest == frozen, path
source_metadata = root / '03_reference/GSE282122/audit_v2/output/00_schema_audit/author_obs_metadata.parquet'
hierarchy = pd.read_parquet(source_metadata, columns=['bucket', 'sub_bucket', 'major', 'minor', 'final_analysis'])
neutrophil_labels = {field: [str(x) for x in hierarchy[field].dropna().unique() if 'neutro' in str(x).lower()] for field in hierarchy.columns}
candidate_neutrophil = {field: int(candidate[field].astype(str).str.contains('neutro', case=False).sum()) for field in hierarchy.columns}

# Step3. 固定供者划分与开发清单交叉核查
development_path = root / '03_reference/web_reference_three_pilots/validation_joint20_20261006/selected_pool_manifest.tsv'
development = pd.read_csv(development_path, sep='\t')
donor_overlap = sorted(set(cells.Patient) & set(development.Patient))
cell_overlap = sorted(set(cells.cell_id) & set(development.cell_id))
assert not donor_overlap and not cell_overlap
missing_train_donors = sorted(set(roles.loc[roles.split.eq('training'), 'Patient']) - set(cells.Patient))

# Step4. 保存只读运行快照，检查进程CPU累计值及实际输出路径
pid = int((run / 'logs/formal_R.pid').read_text())
processes = subprocess.check_output(['ps', '-ww', '-eo', 'pid,ppid,stat,args'], text=True)
process_lines = [line for line in processes.splitlines() if 'run_full_joint20.R' in line or 'RSOCKnode.R' in line]
proc_ticks = {}
for line in process_lines:
    process_pid = int(line.split()[0])
    stat = Path(f'/proc/{process_pid}/stat').read_text().split(') ', 1)[1].split()
    proc_ticks[str(process_pid)] = {'utime_ticks': int(stat[11]), 'stime_ticks': int(stat[12]), 'state': stat[0], 'starttime_ticks': int(stat[19])}
tmux = subprocess.check_output(['tmux', 'list-panes', '-t', 'joint20_bp_ibd_20261014', '-F', '#{session_name} #{pane_pid} #{pane_current_command} #{pane_dead}'], text=True).strip()
logs = {}
for name in ['formal.stdout.log', 'formal.stderr.log', 'formal_command.txt', 'current_stage.txt', 'started_at.txt', 'formal_exit_status.txt', 'exit_status.txt', 'completed_at.txt', 'existing_result_inventory.txt', 'script_SHA256.txt', 'input_preparation.stdout.log']:
    path = run / 'logs' / name
    logs[name] = path.read_text() if path.exists() else None
build_log = (training / 'build_20261007.log').read_text()
files = [{'path': str(path), 'bytes': path.stat().st_size, 'mtime': datetime.fromtimestamp(path.stat().st_mtime).isoformat()} for path in (run / '02_bayesprism').iterdir() if path.is_file()]
result = {'timestamp': datetime.now().astimezone().isoformat(), 'reference_cells_exact': True, 'integer_raw_counts': True, 'source_counts_shape': [9294, 33075], 'samples': len(header), 'mixture_header_equals_mainline_ID_order': True, 'legacy_mainline_IBD_ID_set_exact': True, 'nonIBD_included': 0, 'sha256': hash_records, 'construction_input_sha256': candidate_hashes, 'sampling_replay_cells_and_profiles_exact': True, 'source_metadata_cells': len(hierarchy), 'candidate_cells': len(candidate), 'identity_keep_or_glial_tuft': int(identity.sum()), 'training_identity_cells': int(training_mask.sum()), 'after_doublet_exclusion': int(doublet_mask.sum()), 'final_training_pool_cells': len(pool), 'training_pool_donors': pool.Patient.nunique(), 'unit_count': len(units), 'eligible_units_at_least_20': int(units.available_cells.ge(20).sum()), 'missing_train_donors': missing_train_donors, 'candidate_neutrophil_counts': candidate_neutrophil, 'full_source_neutrophil_labels': neutrophil_labels, 'development_cells': len(development), 'development_donors': sorted(development.Patient.unique().tolist()), 'development_cell_overlap': len(cell_overlap), 'development_donor_overlap': len(donor_overlap), 'processes': process_lines, 'proc_ticks': proc_ticks, 'tmux': tmux, 'formal_PID': pid, 'logs': logs, 'build_log': build_log, 'outputs': files, 'sampling_script_sha256': hashlib.sha256(script.encode()).hexdigest()}
payload = base64.b64encode(gzip.compress(json.dumps(result, ensure_ascii=False, default=int).encode())).decode()
print('AUDIT_PAYLOAD_BEGIN')
for i in range(0, len(payload), 64):
    print(payload[i:i + 64])
print('AUDIT_PAYLOAD_END')
