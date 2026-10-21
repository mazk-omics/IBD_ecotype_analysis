# 审计冻结模型及两轮单细胞恢复，提取同一批细胞的完整计数；原件只读，不重新抽样。
# 服务器 ibd_refaudit 环境运行；所有新表和全基因稀疏矩阵写入本分析目录。
# counts 为基因×细胞；同名基因相加，保留完整基因库和 Unassigned 供后续核对。
from pathlib import Path
import gzip
import hashlib
import json
import shutil
import subprocess
import sys
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread, mmwrite

# Step1. 现场核对模型，而不是把既有解释表当成模型定义
root = Path('/home/mazekai/IBD_EcoTyper')
out = root / '06_validation/paper_aligned_recovery_characterization_20261012'
model = root / 'software/ecotyper/EcoTyper/MSCCR_CORE12_CTSPEC_REGION_IBD/Cell_type_specific_genes'
discovery = root / '05_EcoTyper/output/MSCCR_CORE12/ctspec_region_ibd'
states = pd.read_csv(model / 'Ecotypes/discovery/ecotypes.txt', sep='\t').rename(columns={'ID': 'CellState'})
assert states.CellState.is_unique and len(states) == 53 and states.Ecotype.nunique() == 7 and states.CellType.nunique() == 12
states.to_csv(out / '00_audit/frozen_states.tsv', sep='\t', index=False)
prior = root / '06_validation/ecotyper_original_method_annotation_20261010'
prior_states = pd.read_csv(prior / 'inputs/frozen_states.tsv', sep='\t')
pd.testing.assert_frame_equal(states[['CellState', 'CellType', 'State', 'Ecotype']].sort_values('CellState').reset_index(drop=True), prior_states[['CellState', 'CellType', 'State', 'Ecotype']].sort_values('CellState').reset_index(drop=True))
print('冻结模型：', len(states), 'states;', states.Ecotype.nunique(), 'ecotypes', flush=True)
manifest = []
model_manifest = []
for path in list((model / 'Ecotypes/discovery').glob('*.txt')) + list(discovery.glob('*/gene_info.txt')) + list(discovery.glob('*/state_assignment.txt')) + [discovery / 'rank_data.txt']:
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    model_manifest.append({'File': str(path), 'Bytes': path.stat().st_size, 'SHA256': digest})
pd.DataFrame(model_manifest).to_csv(out / '00_audit/frozen_model_manifest.tsv', sep='\t', index=False)
genes = []
for cell_type in states.CellType.unique():
    info = pd.read_csv(discovery / cell_type / 'gene_info.txt', sep='\t')
    info['CellType'] = cell_type
    info['CellState'] = cell_type + '_' + info.State
    genes.append(info[['Gene', 'CellType', 'State', 'CellState', 'MaxFC']])
pd.concat(genes).to_csv(out / '00_audit/discovery_gene_assignment.tsv', sep='\t', index=False)
for name in ['frozen_states.tsv', 'frozen_ecotypes.tsv', 'literature_marker_panels.tsv', 'discovery_associated_genes.tsv']:
    shutil.copyfile(prior / 'inputs' / name, out / '00_audit' / ('prior_' + name))
subprocess.run(['git', '-C', str(root / 'software/ecotyper'), 'rev-parse', 'HEAD'], stdout=(out / '00_audit/official_git_commit.txt').open('w'), check=True)
subprocess.run(['git', '-C', str(root / 'software/ecotyper'), 'status', '--porcelain'], stdout=(out / '00_audit/official_git_status_before.txt').open('w'), check=True)

# Step2. 固定恢复细胞和真实疾病/肠段；HC 与 Healthy 只在派生标签中合并为 Control
cohorts = [('Colon_GSE214695', 'GSE214695_scRNA_recovery_20261009'), ('Ileum_IBDverse', 'IBDverse_scRNA_recovery_20261009')]
eligibility = []
coverage = []
for dataset, folder in cohorts:
    path = root / '06_validation' / folder / 'audit/cell_assignments_all.tsv'
    cells = pd.read_csv(path, sep='\t')
    assert cells.ID.is_unique and cells.FormalRecovery.all()
    cells['CellState'] = np.where(cells.State.eq('Unassigned'), 'Unassigned', cells.CellType + '_' + cells.State)
    assert set(cells.CellState) <= set(states.CellState) | {'Unassigned'}
    cells['DiseaseGroup'] = cells.Disease.replace({'HC': 'Control', 'Healthy': 'Control'})
    assert cells.DiseaseGroup.isin(['CD', 'UC', 'Control']).all()
    if dataset == 'Ileum_IBDverse':
        assert cells.biopsy_type.eq('terminal ileum').all()
        cells['Region3'] = 'Ileum'
        cells['InflammationGroup'] = cells.Inflammation.replace({'uninflamed': 'Non-inflamed', 'inflamed': 'Inflamed'})
    else:
        cells['InflammationGroup'] = cells.Stratum.map({'CD_inflamed': 'Inflamed', 'UC_inflamed': 'Inflamed', 'HC_nonIBD': 'Non-inflamed'})
        cells['InflammationProvenance'] = '既有作者来源分层；HC 非肿瘤手术组织，非临床炎症评分'
    assert not cells[['Sample', 'Donor']].drop_duplicates().Sample.duplicated().any()
    cells.to_csv(out / '00_audit' / (dataset + '_cells.tsv'), sep='\t', index=False)
    manifest.append({'Role': 'Frozen recovery assignments', 'File': str(path), 'Bytes': path.stat().st_size, 'SHA256': hashlib.sha256(path.read_bytes()).hexdigest()})
    counts = cells.groupby(['CellType', 'CellState', 'DiseaseGroup', 'Region3'], observed=True).agg(Cells=('ID', 'size'), Donors=('Donor', 'nunique'), Samples=('Sample', 'nunique')).reset_index()
    counts.insert(0, 'Dataset', dataset)
    coverage.append(counts)
    for row in states.itertuples():
        parent = cells[cells.CellType.eq(row.CellType)]
        assigned = parent[parent.CellState.eq(row.CellState)]
        eligibility.append({'Dataset': dataset, 'CellState': row.CellState, 'CellType': row.CellType, 'InputCells': len(parent), 'AssignedCells': len(assigned), 'InputDonors': parent.Donor.nunique(), 'AssignedDonors': assigned.Donor.nunique(), 'InputSamples': parent.Sample.nunique(), 'AssignedSamples': assigned.Sample.nunique(), 'UnassignedCells': int(parent.State.eq('Unassigned').sum()), 'MarkerEligible': len(assigned) >= 5, 'LineagePresent': len(parent) > 0})
pd.concat(coverage).to_csv(out / '00_audit/disease_region_state_coverage.tsv', sep='\t', index=False)
pd.DataFrame(eligibility).to_csv(out / '00_audit/dataset_eligibility.tsv', sep='\t', index=False)
pd.DataFrame(manifest).to_csv(out / '00_audit/input_manifest.tsv', sep='\t', index=False)
(out / '00_audit/preflight_identity_pass.json').write_text(json.dumps({'FrozenModel': 'MSCCR_CORE12_CTSPEC_REGION_IBD', 'States': 53, 'Ecotypes': 7, 'PriorMembershipIdentical': True, 'RawFilesReadOnly': True, 'n4Variants': ['DiseaseGroup', 'Region3'], 'Upstream': 'BayesPrism G-like; sensitivity analysis; not CIBERSORTx-equivalent'}, indent=2))
print('Preflight identity PASS; disease/site coverage saved', flush=True)

# Step3. 结直肠逐样本读取作者完整计数，精确匹配 barcode，禁止使用投影基因子集做 DEG
dataset = 'Colon_GSE214695'
source_root = root / '06_validation/GSE214695_scRNA_recovery_20261009'
cells = pd.read_csv(out / '00_audit' / (dataset + '_cells.tsv'), sep='\t')
downloads = pd.read_csv(source_root / 'download_urls.tsv', sep='\t')
blocks = []
ordered_cells = []
symbols_ref = None
for sample, selected in cells.groupby('Sample', sort=True):
    files = downloads[downloads.Sample.eq(sample)]
    matrix_path = source_root / 'source' / files.loc[files.Filename.str.endswith('_matrix.mtx.gz'), 'Filename'].item()
    feature_path = source_root / 'source' / files.loc[files.Filename.str.endswith('_features.tsv.gz'), 'Filename'].item()
    barcode_path = source_root / 'source' / files.loc[files.Filename.str.endswith('_barcodes.tsv.gz'), 'Filename'].item()
    features = pd.read_csv(feature_path, sep='\t', header=None)
    barcodes = pd.read_csv(barcode_path, sep='\t', header=None).iloc[:, 0].astype(str)
    symbols = features.iloc[:, 1].astype(str).to_numpy()
    assert barcodes.is_unique and (pd.Index(barcodes).get_indexer(selected.Barcode) >= 0).all()
    if symbols_ref is None:
        symbols_ref = symbols
    assert np.array_equal(symbols_ref, symbols), '不同样本基因顺序需显式匹配，不能直接拼接'
    with gzip.open(matrix_path, 'rb') as handle:
        full_counts = mmread(handle).tocsc()
    assert full_counts.shape == (len(features), len(barcodes))
    selected_counts = full_counts[:, pd.Index(barcodes).get_indexer(selected.Barcode)]
    blocks.append(selected_counts)
    ordered_cells.append(selected)
    for path in [matrix_path, feature_path, barcode_path]:
        manifest.append({'Role': 'Raw count source', 'File': str(path), 'Bytes': path.stat().st_size, 'SHA256': hashlib.sha256(path.read_bytes()).hexdigest()})
    print('Full count extraction:', sample, selected_counts.shape, flush=True)
full_counts = sparse.hstack(blocks, format='csr')
ordered = pd.concat(ordered_cells, ignore_index=True)
unique_symbols, inverse = np.unique(symbols_ref, return_inverse=True)
collapse = sparse.csr_matrix((np.ones(len(symbols_ref)), (inverse, np.arange(len(symbols_ref)))), shape=(len(unique_symbols), len(symbols_ref)))
counts = (collapse @ full_counts).tocsr()
assert np.isfinite(counts.data).all() and (counts.data >= 0).all() and np.allclose(counts.data, counts.data.round())
for cell_type, selected in ordered.groupby('CellType', sort=True):
    directory = out / '02_external_markers/full_counts' / dataset / cell_type
    directory.mkdir(parents=True, exist_ok=False)
    subset = counts[:, selected.index.to_numpy()]
    with gzip.open(directory / 'counts.mtx.gz', 'wb') as handle:
        mmwrite(handle, subset)
    sparse.save_npz(directory / 'counts.npz', subset)
    pd.Series(unique_symbols).to_csv(directory / 'genes.tsv', index=False, header=False, sep='\t')
    selected.to_csv(directory / 'cells.tsv', index=False, sep='\t')

# Step4. 回肠仅提取冻结细胞的原始 CSR 行，使用既有分类字段解码函数
sys.path.insert(0, str(root / 'scripts/python'))
from inspect_ibdverse_h5ad import read_column
dataset = 'Ileum_IBDverse'
cells = pd.read_csv(out / '00_audit' / (dataset + '_cells.tsv'), sep='\t').sort_values('SourceRow').reset_index(drop=True)
source = Path('/home/mazekai/enteric_glia/data_scRNA/scIBD_external/IBDverse_E-MTAB-16998/E-MTAB-16999/Full_cohort.h5ad')
with h5py.File(source, 'r') as h5:
    index_key = h5['obs'].attrs['_index']
    ids = read_column(h5['obs'][index_key]).astype(str).to_numpy()
    assert np.array_equal(ids[cells.SourceRow.to_numpy()], cells.ID.to_numpy()), 'SourceRow 必须与冻结 ID 完全一致'
    layer = h5['layers/counts']
    assert layer.attrs['encoding-type'] == 'csr_matrix'
    shape = layer.attrs['shape']
    ptr = layer['indptr'][:]
    blocks = []
    for row in cells.SourceRow:
        start, stop = ptr[row:row + 2]
        block = sparse.csr_matrix((layer['data'][start:stop], layer['indices'][start:stop], [0, stop-start]), shape=(1, int(shape[1])))
        blocks.append(block)
    raw_counts = sparse.vstack(blocks, format='csr')
    symbols = read_column(h5['var']['gene_symbols']).astype(str).to_numpy()
unique_symbols, inverse = np.unique(symbols, return_inverse=True)
collapse = sparse.csr_matrix((np.ones(len(symbols)), (np.arange(len(symbols)), inverse)), shape=(len(symbols), len(unique_symbols)))
counts = (raw_counts @ collapse).T.tocsr()
assert np.isfinite(counts.data).all() and (counts.data >= 0).all() and np.allclose(counts.data, counts.data.round())
for cell_type, selected in cells.groupby('CellType', sort=True):
    directory = out / '02_external_markers/full_counts' / dataset / cell_type
    directory.mkdir(parents=True, exist_ok=False)
    subset = counts[:, selected.index.to_numpy()]
    with gzip.open(directory / 'counts.mtx.gz', 'wb') as handle:
        mmwrite(handle, subset)
    sparse.save_npz(directory / 'counts.npz', subset)
    pd.Series(unique_symbols).to_csv(directory / 'genes.tsv', index=False, header=False, sep='\t')
    selected.to_csv(directory / 'cells.tsv', index=False, sep='\t')
manifest.append({'Role': 'Raw full h5ad; counts layer subset by frozen ID', 'File': str(source), 'Bytes': source.stat().st_size, 'SHA256': ''})
pd.DataFrame(manifest).to_csv(out / '00_audit/input_manifest.tsv', sep='\t', index=False)
# 大型原件的整文件 SHA256 单列保存；避免读矩阵时重复扫描25GB。
subprocess.run(['sha256sum', str(source)], stdout=(out / '00_audit/IBDverse_full_file.sha256').open('w'), check=True)
(out / '02_external_markers/full_counts_ready.txt').write_text('PASS: exact frozen cells, full transcriptome, raw counts, gene symbols summed\n')
print('Full transcriptome export PASS', flush=True)
