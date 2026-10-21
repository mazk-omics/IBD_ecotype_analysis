# 核对已有官方比例、训练清单与模拟真值，不运行或覆盖 CIBERSORTx。
# 在本机科学 Python 环境运行；输出独立敏感性目录，并保存同一 SSH 会话需要的小型清单。
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

# Step1. 新建输出，保护历史文件
project = Path(__file__).resolve().parents[2]
out = project / '06_validation/deconv_fraction_sensitivity_joint20_20261013'
if out.exists():
    out = out / pd.Timestamp.now().strftime('run_%H%M%S')
out.mkdir(parents=True, exist_ok=False)
for folder in ['00_audit', '01_reference', '02_bayesprism_joint20_MSCCR120', '03_cibersortx_existing', '04_MSCCR120_concordance', '05_existing_pseudobulk_benchmark', '06_figures', '07_reports', 'logs', 'scripts']:
    (out / folder).mkdir()
(project / '06_validation/fraction_sensitivity_current_path.txt').write_text(str(out), encoding='utf-8')
ready = project / 'CIBERSORTx_ready_20261007'
training_dir = ready / '06_Joint20v2_DonorSite_reference/核对记录_勿上传'
job20_dir = project / 'CIBERSORTx_signature_audit/MSCCR_Joint20v2_Job20_review_20261009'
signature_path = ready / '09_MSCCR_Joint20v2_Fractions/GJ20v2_Signature20.txt'
signature = pd.read_csv(signature_path, sep='\t', index_col=0)
types = signature.columns.tolist()
assert signature.shape == (4474, 20) and signature.index.is_unique and signature.columns.is_unique
pd.DataFrame({'CIBERSORTx_Job13': types, 'BayesPrism_type': types, 'mapping': 'exact frozen label'}).to_csv(out / '01_reference/cell_type_mapping.tsv', sep='\t', index=False)

# Step2. 核对正式样本和官方输出，不重新归一化
annotation = pd.read_csv(project / 'CIBERSORTx_ready_20261008/01_MSCCR_platform120_pilot/核对记录_勿上传/selected_annotation.tsv', sep='\t', dtype=str)
official_paths = {20: job20_dir / 'original_outputs/CIBERSORTx_Job20_Adjusted.txt', 17: project / 'CIBERSORTx_signature_audit/MSCCR_Joint20v2_Job17_review_20261008/original_outputs/CIBERSORTx_Job17_Adjusted.txt'}
official_paths[14] = project / 'CIBERSORTx_signature_audit/Joint20v2_Job14_review_20261007/original_outputs/CIBERSORTx_Job14_Results.txt'
official_paths[15] = project / 'CIBERSORTx_signature_audit/Joint20v2_Job15_review_20261007/original_outputs/CIBERSORTx_Job15_Results.txt'
inventory = []
for job, path in official_paths.items():
    frame = pd.read_csv(path, sep='\t', index_col=0)
    values = frame[types].to_numpy()
    assert frame.index.is_unique and np.isfinite(values).all() and (values >= 0).all()
    assert np.allclose(values.sum(axis=1), 1, atol=1e-8)
    inventory.append({'job': job, 'official_path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest(), 'bytes': path.stat().st_size, 'samples': len(frame), 'types': len(types), 'relative_fractions_valid': True, 'software_version': 'official web build not recorded'})
    if job in [20, 17]:
        assert set(frame.index) == set(annotation.ID)
assert len(annotation) == 120 and annotation.ID.is_unique and annotation.PatientID.nunique() == 120
annotation.to_csv(out / '00_audit/MSCCR120_annotation.tsv', sep='\t', index=False)
pd.DataFrame(inventory).to_csv(out / '03_cibersortx_existing/official_job_inventory.tsv', sep='\t', index=False)
log_path = job20_dir / 'original_outputs/Job20_runtime_excerpt_user_provided.txt'
log = log_path.read_text(encoding='utf-8')
assert all(x in log for x in ['rmbatchSmode: TRUE', 'QN: FALSE', 'relative', 'perm: 100', 'All done.'])
(out / '03_cibersortx_existing/job20_QC.md').write_text('# Job20 QC\n\nCOMPLETE: existing official Adjusted fractions; 120 unique patients, full 20 types, finite nonnegative rows summing to one. No renormalization.\n\nSettings from a user-provided runtime excerpt: relative, S-mode with SourceSubset1060, QN disabled, 100 permutations. Full runtime log, web build and uploaded-file hashes are unavailable. Signature originated from 212 donor/site pooled CPM profiles, while the S-mode source contains 1060 individual training cells. This is a documented exploratory correction pairing. Job17 is retained as repeat-run sensitivity.\n', encoding='utf-8')

# Step3. 区分训练、供者汇总和 S-mode 子集
cells = pd.read_csv(training_dir / 'selected_cell_manifest.tsv', sep='\t')
profiles = pd.read_csv(training_dir / 'profile_manifest.tsv', sep='\t')
roles = pd.read_csv(training_dir / 'frozen_donor_roles.tsv', sep='\t')
subset_path = project / 'CIBERSORTx_ready_20261008/01_MSCCR_platform120_pilot/核对记录_勿上传/source_subset_cell_manifest.tsv'
subset = pd.read_csv(subset_path, sep='\t')
assert len(cells) == 9294 and cells.Patient.nunique() == 32 and len(profiles) == 212
assert cells.cell_id.is_unique and cells.obs_name.is_unique and cells.split.eq('training').all()
assert set(cells.coarse_identity) == set(types) and set(subset.cell_id) <= set(cells.cell_id) and len(subset) == 1060
summary = json.loads((training_dir / 'construction_summary.json').read_text(encoding='utf-8'))
cells['Dataset'] = 'GSE282122'
cells['Raw_counts_source'] = summary['source'] + ':/X; selected_training_counts_csr.npz checkpoint'
cells['used_signature'] = True
cells['used_Smode'] = cells.cell_id.isin(subset.cell_id)
cells['used_Job14_15'] = False
cells.to_csv(out / '01_reference/reference_provenance.tsv', sep='\t', index=False)
provenance = cells.groupby(['Dataset', 'sample_id', 'Patient', 'Site', 'coarse_identity', 'split', 'used_signature', 'used_Smode', 'used_Job14_15', 'Raw_counts_source'], observed=True).size().rename('Number_of_cells').reset_index()
provenance.to_csv(out / '00_audit/signature_reference_provenance.tsv', sep='\t', index=False)
pd.DataFrame({'cell_type': types, 'cells': [int(cells.coarse_identity.eq(t).sum()) for t in types], 'donors': [int(cells.loc[cells.coarse_identity.eq(t), 'Patient'].nunique()) for t in types], 'profiles': [int(profiles.coarse_identity.eq(t).sum()) for t in types]}).to_csv(out / '01_reference/reference_type_capacity.tsv', sep='\t', index=False)
roles.to_csv(out / '00_audit/frozen_donor_roles.tsv', sep='\t', index=False)

# Step4. 复用已有模拟审计数字和真值，新的准确性结果需要 BayesPrism 实际输出
audit_records = []
for job, folder, truth_file in [(14, '07_Joint20v2_EqualRNA8_validation', 'ground_truth_equalRNA8_fractions.tsv'), (15, '08_Joint20v2_Natural28_validation', 'ground_truth_export_UMI_fractions.tsv')]:
    truth_path = ready / folder / '核对记录_勿上传' / truth_file
    truth = pd.read_csv(truth_path, sep='\t', index_col=0)
    prediction = pd.read_csv(official_paths[job], sep='\t', index_col=0).loc[truth.index, types]
    truth = truth[types]
    assert np.allclose(truth.sum(axis=1), 1)
    audit_records.append({'job': job, 'samples': len(truth), 'historical_official_MAE_pp': float(np.abs(prediction.to_numpy() - truth.to_numpy()).mean() * 100), 'truth_path': str(truth_path), 'truth_sha256': hashlib.sha256(truth_path.read_bytes()).hexdigest(), 'truth_definition': 'equal RNA/cell weights' if job == 14 else 'exported full-gene UMI contribution'})
pd.DataFrame(audit_records).to_csv(out / '00_audit/existing_pseudobulk_audit.tsv', sep='\t', index=False)
input_paths = [signature_path, log_path, subset_path, training_dir / 'selected_cell_manifest.tsv', training_dir / 'source_gene_manifest.tsv', training_dir / 'construction_summary.json']
pd.DataFrame([{'path': str(p), 'bytes': p.stat().st_size, 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()} for p in input_paths]).to_csv(out / '00_audit/local_input_manifest.tsv', sep='\t', index=False)
(out / '00_audit/frozen_types.txt').write_text('\n'.join(types) + '\n', encoding='utf-8')
(out / 'execution_status.json').write_text(json.dumps({'phase0_local_audit': 'COMPLETE', 'phase1_reference': 'NOT_RUN', 'phase2_MSCCR120': 'NOT_RUN', 'phase3_CIBERSORTx': 'COMPLETE', 'phase4_concordance': 'NOT_RUN', 'phase5_pseudobulk': 'NOT_RUN', 'phase6_reports': 'NOT_RUN', 'SSH_session_id': 69993, 'CIBERSORTx_rerun': False}, indent=2), encoding='utf-8')
print('OUTPUT:', out)
print(pd.DataFrame(audit_records).to_string(index=False))
print('Reference:', len(cells), 'cells;', len(profiles), 'profiles;', cells.Patient.nunique(), 'donors')
