# 通过已建立的同一个 SSH 隧道取回 Joint20 v2，核对哈希、网页格式和供者隔离。
# 在本机 .tools/cibersortx_review Python 环境运行；SSH 的本地转发端口为 18527，不创建任何 SSH/SCP 连接。
# 网页上传文件直接放在新目录，清单和诊断分开放在“核对记录_勿上传”；不会改动旧候选或验证结果。
from pathlib import Path
import hashlib
import json
import urllib.request
import numpy as np
import pandas as pd

# Step1. 读入服务器文件清单
out = Path('D:/IBD_ecotype/IBD_ecotype_analysis/CIBERSORTx_ready_20261007/06_Joint20v2_DonorSite_reference')
audit = out / '核对记录_勿上传'
audit.mkdir(parents=True, exist_ok=True)
base_url = 'http://127.0.0.1:18527/'
with urllib.request.urlopen(base_url + 'transfer_checksums.tsv', timeout=30) as response:
    (audit / 'transfer_checksums.tsv').write_bytes(response.read())
files = pd.read_csv(audit / 'transfer_checksums.tsv', sep='\t')
assert files.file.is_unique and files.file.map(lambda x: Path(x).name == x).all()
web_names = ['GJ20v2_DonorSite_CPM.txt', 'GJ20v2_PhenotypeClasses.txt']
assert set(web_names).issubset(files.file)

# Step2. 沿同一隧道下载松散文件，每个文件立即核对 SHA256
local_checks = []
for row in files.itertuples():
    destination = (out if row.file in web_names else audit) / row.file
    with urllib.request.urlopen(base_url + row.file, timeout=180) as response:
        data = response.read()
    assert len(data) == row.bytes and hashlib.sha256(data).hexdigest() == row.sha256, row.file
    destination.write_bytes(data)
    local_checks.append({'file': row.file, 'bytes': len(data), 'sha256_verified': True})
    print('Retrieved and verified:', row.file, len(data), flush=True)
pd.DataFrame(local_checks).to_csv(audit / 'local_transfer_verification.tsv', sep='\t', index=False)

# Step3. 核对两份网页文件的数值和列顺序
reference = pd.read_csv(out / web_names[0], sep='\t', index_col=0)
classes = pd.read_csv(out / web_names[1], sep='\t', header=None, index_col=0)
profiles = pd.read_csv(audit / 'profile_manifest.tsv', sep='\t')
summary = json.loads((audit / 'construction_summary.json').read_text(encoding='utf-8'))
assert reference.shape == (summary['uploaded_genes'], summary['profiles'])
assert classes.shape == (len(summary['types']), summary['profiles'])
assert reference.columns.tolist() == profiles.profile_id.tolist()
classes.columns = reference.columns  # 类别文件按网页要求没有标题行，数值列由参考表达标题定位。
assert classes.index.tolist() == summary['types']
assert reference.index.is_unique and reference.index.notna().all()
assert reference.columns.is_unique and reference.size < 50_000_000
assert np.isfinite(reference.to_numpy()).all() and reference.ge(0).all().all()
assert reference.to_numpy().max() > 50
assert np.isin(classes.to_numpy(), [1, 2]).all() and classes.eq(1).sum(axis=0).eq(1).all()
for cell_type in summary['types']:
    assert classes.columns[classes.loc[cell_type].eq(1)].tolist() == profiles.loc[profiles.coarse_identity.eq(cell_type), 'profile_id'].tolist()
assert np.allclose(reference.sum(axis=0), 1e6, rtol=0, atol=1)

# 逐行确认字段数、无双引号和缺失项，避免小文件看起来正确但中间有坏行。
for name in web_names:
    rows = 0
    with (out / name).open(encoding='utf-8') as handle:
        for line in handle:
            fields = line.rstrip('\r\n').split('\t')
            assert len(fields) == len(profiles) + 1 and all(fields) and '"' not in line, (name, rows)
            rows += 1
    assert rows == (len(reference) + 1 if name == web_names[0] else len(summary['types']))

# Step4. 核对每列的供者、部位、样本和类型，留出供者不允许进入训练
selected = pd.read_csv(audit / 'selected_cell_manifest.tsv', sep='\t')
roles = pd.read_csv(audit / 'frozen_donor_roles.tsv', sep='\t')
assert selected.cell_id.is_unique and selected.obs_name.is_unique and selected.split.eq('training').all()
assert not selected.Patient.isin(roles.loc[roles.split.ne('training'), 'Patient']).any()
checks = selected.groupby('profile_id')[['Patient', 'Site', 'sample_id', 'coarse_identity']].nunique()
assert checks.eq(1).all().all()
observed_counts = selected.groupby('profile_id').size().reindex(profiles.profile_id)
assert np.array_equal(observed_counts.to_numpy(), profiles.selected_cells.to_numpy())
assert observed_counts.between(summary['minimum_cells_per_profile'], summary['target_cells_per_profile']).all()
for field in ['Patient', 'Site', 'sample_id', 'coarse_identity']:
    actual = selected.groupby('profile_id')[field].first().reindex(profiles.profile_id)
    assert np.array_equal(actual.to_numpy(), profiles[field].to_numpy())
composition = pd.read_csv(audit / 'reference_type_region_composition.tsv', sep='\t')
shared = composition.loc[~composition.coarse_identity.isin(['Ileal_absorptive', 'Colonic_absorptive'])]
regional_counts = shared.pivot(index='coarse_identity', columns='region', values='profiles')
assert regional_counts.Ileum.eq(regional_counts.ColonRectum).all()
assert not profiles.duplicated(['coarse_identity', 'region', 'Patient']).any()

# Step5. 输出简短核对记录，区分输入准备完成与网页生成/验证完成
verification = {'web_files': web_names, 'all_transfer_hashes_match': True, 'reference_shape_genes_by_profiles': list(reference.shape), 'phenotype_shape_types_by_profiles': list(classes.shape), 'expression_data_points': reference.size, 'selected_cells': len(selected), 'selected_donors': selected.Patient.nunique(), 'profiles_single_donor_site_sample_type': True, 'shared_types_balanced_profile_counts_by_region': True, 'holdout_cells_used': 0, 'CPM_column_sum_min': float(reference.sum().min()), 'CPM_column_sum_max': float(reference.sum().max()), 'CPM_max_expression': float(reference.to_numpy().max()), 'web_signature_generated': False, 'signature_holdout_validation_done': False}
(audit / 'local_format_and_isolation_verification.json').write_text(json.dumps(verification, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
print(json.dumps(verification, ensure_ascii=False, indent=2), flush=True)
