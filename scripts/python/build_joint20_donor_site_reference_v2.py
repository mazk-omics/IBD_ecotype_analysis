# 从 GSE282122 原始 counts 抽样，按单一供者、部位、样本和类型汇总，供 CIBERSORTx 网页选择 signature genes。
# 在服务器 ibd_refaudit 环境运行；读取既有候选池与冻结供者划分，输出新的 joint20_donor_site_v2 目录，不修改源数据。
# 下方参数为本轮项目约定；输出列是供者伪 bulk，网页应选 RNA-Seq 参考模式，不能作为单细胞列上传。
from pathlib import Path
import hashlib
import json
import sys
import time
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

# 仅复用已验证的 H5AD 注释解码和 CSR 行读取；抽样、汇总和导出在下面直接展开。
root = Path('/home/mazekai/IBD_EcoTyper')
helper_dir = root / 'scripts/python/three_reference_pilots_v1'
sys.path.insert(0, str(helper_dir))
from inspect_ibdverse_h5ad import read_column
from export_ibdverse_cibersortx_pilot import selected_csr

# Step1. 读入已审计候选池，保留冻结训练供者
started = time.time()
source = Path('/home/mazekai/enteric_glia/data_scRNA/GSE282122/TAURUS_raw_counts_annotated_final.h5ad')
candidate_file = root / '03_reference/GSE282122/audit_v2/output/04_identity_coverage_v1/mapped_candidate_cells_v1.parquet'
roles_file = root / '03_reference/web_reference_three_pilots/validation_joint20_20261006/donor_roles.tsv'
out = root / '03_reference/joint20_donor_site_v2/output/run_20261007_v2'
out.mkdir(parents=True, exist_ok=False)  # 防止覆盖本轮或历史输出。
candidate = pd.read_parquet(candidate_file)
donor_roles = pd.read_csv(roles_file, sep='\t')
assert donor_roles.Patient.is_unique and candidate.obs_name.is_unique
assert candidate.reference_eligible.eq(True).all()
assert (candidate.Treatment.eq('Pre') | candidate.Disease.eq('Healthy')).all()
candidate['split'] = candidate.Patient.map(donor_roles.set_index('Patient').split)
assert candidate.split.notna().all()
types = ['B', 'Plasma', 'CD4_T', 'CD8_T', 'NK', 'ILC', 'Monocyte_Macrophage', 'DC', 'Mast', 'Fibroblast', 'Pericyte', 'Endothelial', 'Ileal_absorptive', 'Colonic_absorptive', 'Goblet', 'Stem_TA_progenitor', 'Paneth', 'Enteroendocrine', 'Glial', 'Tuft']
candidate['region'] = np.where(candidate.Site.eq('Terminal_Ileum'), 'Ileum', 'ColonRectum')
action = candidate.primary_reference_action.eq('KEEP') | candidate.coarse_identity.isin(['Glial', 'Tuft'])
valid = action & candidate.identity_status.eq('CONFIDENT') & candidate.coarse_identity.isin(types) & candidate.split.eq('training') & ~candidate.predicted_doublets
valid &= np.isfinite(candidate.total_counts) & candidate.total_counts.gt(0) & np.isfinite(candidate.n_genes_by_counts) & candidate.n_genes_by_counts.gt(0)
pool = candidate.loc[valid].copy().sort_values('cell_id')
donor_roles.to_csv(out / 'frozen_donor_roles.tsv', sep='\t', index=False)
print('Training candidate cells:', len(pool), '; donors:', pool.Patient.nunique(), flush=True)

# Step2. 先查看真实供者容量，再选择参考列
# 每列至少 20、最多 50 细胞；共享类型按两区域相同供者数，上限各 6。稀有类型不复制供者凑数。
min_cells = 20
target_cells = 50
max_donors_per_region = 6
seed = 20261007
rng = np.random.default_rng(seed)
keys = ['coarse_identity', 'region', 'Patient', 'Site', 'sample_id']
units = pool.groupby(keys, observed=True).agg(available_cells=('cell_id', 'size'), Disease=('Disease', 'first'), Inflammation=('Inflammation', 'first')).reset_index()
units.to_csv(out / 'training_unit_capacity.tsv', sep='\t', index=False)
eligible_units = units.loc[units.available_cells.ge(min_cells)].copy()
capacity = eligible_units.groupby(['coarse_identity', 'region'], observed=True).agg(eligible_donors=('Patient', 'nunique'), eligible_units=('sample_id', 'size')).reset_index()
capacity.to_csv(out / 'training_type_region_capacity.tsv', sep='\t', index=False)
print(capacity.to_string(index=False), flush=True)

# 每个类型/区域中每位供者最多贡献一列；结直肠尽量覆盖 Rectum 和至少一个结肠部位。
chosen_units = []
for cell_type in types:
    regions = ['Ileum'] if cell_type == 'Ileal_absorptive' else ['ColonRectum'] if cell_type == 'Colonic_absorptive' else ['Ileum', 'ColonRectum']
    donor_counts = [eligible_units.loc[eligible_units.coarse_identity.eq(cell_type) & eligible_units.region.eq(region), 'Patient'].nunique() for region in regions]
    target_donors = min(max_donors_per_region, min(donor_counts))
    assert target_donors >= 3, f'{cell_type}: fewer than three eligible donors per required region'
    for region in regions:
        q = eligible_units.loc[eligible_units.coarse_identity.eq(cell_type) & eligible_units.region.eq(region)].copy()
        q['random_order'] = rng.random(len(q))
        q = q.sort_values('random_order')
        chosen_indices = []
        used_donors = set()
        if region == 'ColonRectum':
            for rectum in [True, False]:
                anchor = q.loc[q.Site.eq('Rectum').eq(rectum) & ~q.Patient.isin(used_donors)]
                if len(anchor):
                    i = anchor.index[0]
                    chosen_indices.append(i)
                    used_donors.add(q.at[i, 'Patient'])
        for i, row in q.iterrows():
            if len(chosen_indices) >= target_donors:
                break
            if row.Patient not in used_donors:
                chosen_indices.append(i)
                used_donors.add(row.Patient)
        chosen_units.append(q.loc[chosen_indices].drop(columns='random_order'))
profiles = pd.concat(chosen_units, ignore_index=True)
profiles['profile_id'] = [f'GJ2_P{i + 1:03d}' for i in range(len(profiles))]
assert not profiles.duplicated(['coarse_identity', 'region', 'Patient']).any()
profiles.to_csv(out / 'planned_profiles.tsv', sep='\t', index=False)

# Step3. 在每个真实参考单元内随机抽细胞，保存可以追溯的清单
pool_groups = pool.groupby(keys, observed=True).indices
selected_parts = []
for row in profiles.itertuples():
    key = (row.coarse_identity, row.region, row.Patient, row.Site, row.sample_id)
    group_positions = pool_groups[key]
    positions = rng.choice(group_positions, size=min(target_cells, len(group_positions)), replace=False)
    part = pool.iloc[positions].copy()
    part['profile_id'] = row.profile_id
    selected_parts.append(part)
selected = pd.concat(selected_parts, ignore_index=True)
assert selected.cell_id.is_unique and selected.obs_name.is_unique
assert selected.split.eq('training').all()
assert not selected.Patient.isin(donor_roles.loc[donor_roles.split.ne('training'), 'Patient']).any()
manifest_columns = ['profile_id', 'cell_id', 'obs_name', 'Patient', 'Disease', 'Treatment', 'Inflammation', 'Site', 'region', 'sample_id', 'coarse_identity', 'final_analysis', 'pct_counts_mt', 'total_counts', 'n_genes_by_counts', 'split']
selected[manifest_columns].to_csv(out / 'selected_cell_manifest.tsv', sep='\t', index=False)
profile_checks = selected.groupby('profile_id', observed=True)[['Patient', 'Site', 'sample_id', 'coarse_identity']].nunique()
assert profile_checks.eq(1).all().all()
subtypes = selected.groupby(['coarse_identity', 'region', 'final_analysis'], observed=True).size().rename('selected_cells').reset_index()
subtypes.to_csv(out / 'selected_author_subtype_composition.tsv', sep='\t', index=False)
print('Selected profiles:', len(profiles), '; cells:', len(selected), '; donors:', selected.Patient.nunique(), flush=True)

# Step4. 只读取选中细胞的原始 counts，不读取或改写全矩阵
source_before = {'bytes': source.stat().st_size, 'mtime_ns': source.stat().st_mtime_ns}
with h5py.File(source, 'r') as h:
    obs_ids = read_column(h['obs'][h['obs'].attrs.get('_index', '_index')]).astype(str)
    rows = pd.Index(obs_ids).get_indexer(selected.obs_name.astype(str))
    assert obs_ids.is_unique and (rows >= 0).all()
    for field in ['Patient', 'Site', 'sample_id', 'final_analysis']:
        actual = read_column(h['obs'][field]).iloc[rows].astype(str).to_numpy()
        assert np.array_equal(actual, selected[field].astype(str).to_numpy()), field
    order = np.argsort(rows)
    rows = rows[order]
    selected = selected.iloc[order].reset_index(drop=True)
    selected['source_cell_row'] = rows
    genes = pd.DataFrame({'source_gene_id': read_column(h['var']['gene_id']).astype('string'), 'gene_symbol': read_column(h['var']['gene_symbol']).astype('string')})
    indices, values, pointers, cell_totals, detected = selected_csr(h['X'], rows, len(genes))
counts = sparse.csr_matrix((values, indices, pointers), shape=(len(selected), len(genes)))
assert np.allclose(cell_totals, selected.total_counts.to_numpy(), rtol=0, atol=1)
assert np.array_equal(detected, selected.n_genes_by_counts.to_numpy())
sparse.save_npz(out / 'selected_training_counts_csr.npz', counts)  # 服务器检查点：细胞 × 源基因，不作为网页文件。
selected['extracted_total_counts'] = cell_totals
selected[manifest_columns + ['source_cell_row', 'extracted_total_counts']].to_csv(out / 'selected_cell_manifest.tsv', sep='\t', index=False)
genes.to_csv(out / 'source_gene_manifest.tsv', sep='\t', index=False)

# Step5. 同一参考单元累加 counts，使用完整基因库总量计算非 log CPM
profile_index = pd.Index(profiles.profile_id)
profile_codes = profile_index.get_indexer(selected.profile_id)
membership = sparse.csr_matrix((np.ones(len(selected), dtype=np.int64), (profile_codes, np.arange(len(selected)))), shape=(len(profiles), len(selected)))
pooled_counts = (membership @ counts).toarray().astype(np.int64)
full_library_totals = pooled_counts.sum(axis=1)
assert (full_library_totals > 0).all()
np.save(out / 'pooled_raw_counts_profiles_by_source_gene.npy', pooled_counts)
profiles['selected_cells'] = membership.sum(axis=1).A1.astype(int)
profiles['full_library_counts'] = full_library_totals
profiles['median_cell_mito_pct'] = profiles.profile_id.map(selected.groupby('profile_id', observed=True).pct_counts_mt.median())
profiles.to_csv(out / 'profile_manifest.tsv', sep='\t', index=False)

# 同名 gene symbol 累加，空名称不上传；分母始终包含完整源基因库，不能删基因后再重新归一化。
genes['gene_symbol'] = genes.gene_symbol.str.strip()
valid_symbols = genes.gene_symbol.notna() & genes.gene_symbol.ne('') & ~genes.gene_symbol.str.contains(r'[\t\r\n\"]', na=True)
raw_profiles = pd.DataFrame(pooled_counts.T, index=genes.gene_symbol, columns=profile_index)
symbol_counts = raw_profiles.loc[valid_symbols.to_numpy()].groupby(level=0, sort=True).sum()
symbol_counts = symbol_counts.loc[symbol_counts.sum(axis=1).gt(0)]
cpm = symbol_counts.div(pd.Series(full_library_totals, index=profile_index), axis=1) * 1e6
assert cpm.index.is_unique and cpm.columns.is_unique
assert np.isfinite(cpm.to_numpy()).all() and cpm.ge(0).all().all() and cpm.to_numpy().max() > 50
assert cpm.size < 50_000_000
genes.to_csv(out / 'source_gene_manifest.tsv', sep='\t', index=False)
cpm.to_csv(out / 'GJ20v2_DonorSite_CPM.txt', sep='\t', index_label='GeneSymbol', float_format='%.8g')

# Step6. 生成与参考列严格同序的 phenotype classes，网页做每一类对其他类的比较
# 1=该类型参考；2=其他类型比较参考。不能将这些伪 bulk 列填入 scRNA-Seq 模式。
classes = pd.DataFrame(2, index=types, columns=profile_index)
for cell_type in types:
    classes.loc[cell_type, profiles.loc[profiles.coarse_identity.eq(cell_type), 'profile_id']] = 1
assert classes.columns.equals(cpm.columns) and classes.eq(1).sum(axis=0).eq(1).all()
# phenotype classes 与网页 inferred_phenoclasses 一致：首列类型名，没有列标题行；数值列顺序见 profile_manifest。
classes.to_csv(out / 'GJ20v2_PhenotypeClasses.txt', sep='\t', header=False)
composition = profiles.groupby(['coarse_identity', 'region'], observed=True).agg(profiles=('profile_id', 'size'), donors=('Patient', 'nunique'), cells=('selected_cells', 'sum'), min_cells_per_profile=('selected_cells', 'min'), max_cells_per_profile=('selected_cells', 'max')).reset_index()
composition.to_csv(out / 'reference_type_region_composition.tsv', sep='\t', index=False)
profiles.groupby(['Disease', 'region'], observed=True).agg(profiles=('profile_id', 'size'), donors=('Patient', 'nunique'), cells=('selected_cells', 'sum')).to_csv(out / 'reference_disease_composition.tsv', sep='\t')

# Step7. 就地核对关键上皮表达，记录实际构成而不按留出结果调表达系数
markers = ['DEFA5', 'DEFA6', 'LYZ', 'REG3A', 'MUC2', 'TFF3', 'SPINK4', 'LGR5', 'OLFM4', 'MKI67', 'TOP2A', 'CHGA', 'CHGB', 'NEUROG3', 'ALPI', 'APOA4', 'FABP6', 'KRT20', 'BEST4', 'OTOP2', 'SOX10', 'S100B', 'POU2F3', 'TRPM5']
marker_expression = cpm.reindex(markers).T
marker_expression.index.name = 'profile_id'
marker_expression.to_csv(out / 'key_marker_CPM_by_profile.tsv', sep='\t', na_rep='NA', float_format='%.8g')
print('Marker availability:', marker_expression.notna().any(axis=0).to_dict(), flush=True)
print(composition.to_string(index=False), flush=True)

# Step8. 保存构建依据和文件校验；源文件大小与修改时间必须保持不变
source_after = {'bytes': source.stat().st_size, 'mtime_ns': source.stat().st_mtime_ns}
assert source_before == source_after
summary = {'source': str(source), 'candidate_file': str(candidate_file), 'donor_roles_file': str(roles_file), 'counts_layer': '/X', 'source_opened_read_only': True, 'source_before': source_before, 'source_after': source_after, 'environment_python': sys.executable, 'seed': seed, 'max_donors_per_region': max_donors_per_region, 'minimum_cells_per_profile': min_cells, 'target_cells_per_profile': target_cells, 'types': types, 'profiles': len(profiles), 'selected_cells': len(selected), 'selected_donors': selected.Patient.nunique(), 'source_genes': len(genes), 'uploaded_genes': len(cpm), 'expression_data_points': cpm.size, 'holdout_cells_used': 0, 'CPM_column_sum_min': float(cpm.sum().min()), 'CPM_column_sum_max': float(cpm.sum().max()), 'duplicate_symbol_rows': int(genes.gene_symbol.duplicated(keep=False).sum()), 'invalid_symbol_rows': int((~valid_symbols).sum()), 'source_single_cell_web_data_points': len(selected) * len(genes), 'elapsed_seconds': round(time.time() - started, 1), 'web_input_mode': 'RNA-Seq reference samples + phenotype classes; NOT scRNA-Seq', 'unselected_datasets': {'IBDverse': 'No explicit Paneth/Glial/NK labels; retained for later reference/recovery assessment', 'SCP259': 'Cell-resolved colon/rectum location unavailable in current audited metadata; not needed for this one main candidate', 'SCP1884': 'Retained as recovery candidate; no new training cells taken in this round; historical baseline use still requires a donor overlap audit'}}
input_hashes = {}
for path in [candidate_file, roles_file, helper_dir / 'inspect_ibdverse_h5ad.py', helper_dir / 'export_ibdverse_cibersortx_pilot.py']:
    input_hashes[str(path)] = hashlib.file_digest(path.open('rb'), 'sha256').hexdigest()
summary['input_sha256'] = input_hashes
(out / 'construction_summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2, default=int) + '\n', encoding='utf-8')
checksums = []
for path in sorted(out.iterdir()):
    if path.is_file() and path.suffix != '.npz' and path.suffix != '.npy':
        with path.open('rb') as handle:
            digest = hashlib.file_digest(handle, 'sha256').hexdigest()
        checksums.append({'file': path.name, 'bytes': path.stat().st_size, 'sha256': digest})
pd.DataFrame(checksums).to_csv(out / 'transfer_checksums.tsv', sep='\t', index=False)
print('READY:', out, '; reference shape:', cpm.shape, '; data points:', cpm.size, flush=True)
