# 按Cell2021高扰动比例思路，从外部作者注释细胞生成真正的count-based混合物。
# 服务器ibd_refaudit运行；smoke先生成3份，通过后另用冻结种子生成100份，每份1000细胞、有放回。
from pathlib import Path
import gzip
import hashlib
import json
import sys
import time
import h5py
import numpy as np
import pandas as pd
from scipy import io, sparse

# Step1. 读取冻结作者真值池和比例参数，不使用recovery结果
started = time.monotonic()
parent = Path(sys.argv[1])
dataset = sys.argv[2]
mode = sys.argv[3]
folder = parent / '08_external_pseudobulk' / dataset
out = folder / ('01_smoke' if mode == 'smoke' else '02_simulation')
assert not (out / 'generation_summary.json').exists()
if mode == 'formal':
    assert (folder / '01_smoke/smoke_complete.txt').exists(), '正式模拟只能在smoke通过后生成'
cells = pd.read_csv(folder / '00_audit/eligible_author_cells.tsv.gz', sep='\t', keep_default_na=False)
parameters = pd.read_csv(folder / '00_audit/gaussian_sampling_parameters.tsv', sep='\t')
settings = json.loads((folder / '00_audit/endpoint_prespecification.json').read_text())
types = (parent / '00_audit/frozen_types.txt').read_text().splitlines()
active = parameters.cell_type.tolist()
seed = settings['seed_smoke'] if mode == 'smoke' else settings['seed_formal']
n_mixtures = 3 if mode == 'smoke' else 100
rng = np.random.default_rng(seed)
members = []
target_weights = []
for number in range(n_mixtures):
    weights = rng.normal(parameters.mu, parameters.gaussian_SD)
    weights = np.clip(weights, 0, None)
    while weights.sum() == 0:
        weights = np.clip(rng.normal(parameters.mu, parameters.gaussian_SD), 0, None)
    weights /= weights.sum()
    allocations = rng.multinomial(1000, weights)
    sample = f'{dataset}_{mode}_{number + 1:03d}'
    target_weights.append(dict(zip(['sample'] + active, [sample] + weights.tolist())))
    draw_offset = 0
    for cell_type, allocation in zip(active, allocations):
        if allocation == 0:
            continue
        available = cells.index[cells.CellType.eq(cell_type)].to_numpy()
        indices = rng.choice(available, size=allocation, replace=True)
        part = cells.loc[indices, ['CellID', 'Sample', 'Donor', 'AuthorLabel', 'CellType', 'SourceRow']].copy()
        part['sample'] = sample
        part['draw_number'] = np.arange(draw_offset, draw_offset + allocation)
        members.append(part)
        draw_offset += allocation
    assert draw_offset == 1000
membership = pd.concat(members, ignore_index=True)
membership.to_csv(out / 'cell_sampling_manifest.tsv.gz', sep='\t', index=False)
pd.DataFrame(target_weights).to_csv(out / 'target_sampling_fractions.tsv', sep='\t', index=False)
pool = cells.loc[cells.CellID.isin(membership.CellID)].copy()
assert pool.CellID.is_unique and len(membership) == n_mixtures * 1000
print(dataset, mode, 'unique source cells:', len(pool), 'draws:', len(membership), flush=True)

# Step2. 读取实际被抽中的原始稀疏counts；完整基因库保留，不裁剪为模型或signature基因
source_records = []
if dataset == 'IBDverse':
    source = Path('/home/mazekai/enteric_glia/data_scRNA/scIBD_external/IBDverse_E-MTAB-16998/E-MTAB-16999/Full_cohort.h5ad')
    sys.path.insert(0, '/home/mazekai/IBD_EcoTyper/scripts/python/three_reference_pilots_v1')
    from inspect_ibdverse_h5ad import read_column
    from export_ibdverse_cibersortx_pilot import selected_csr
    pool = pool.sort_values('SourceRow').reset_index(drop=True)
    before = {'bytes': source.stat().st_size, 'mtime_ns': source.stat().st_mtime_ns}
    with h5py.File(source, 'r') as h:
        symbols = read_column(h['var']['gene_symbols']).astype(str)
        layer = h['layers/counts']
        positions = pool.SourceRow.to_numpy(dtype=int)
        ids = read_column(h['obs'][h['obs'].attrs.get('_index', '_index')]).astype(str)
        assert np.array_equal(ids.iloc[positions], pool.CellID)
        indices, values, pointers, totals, detected = selected_csr(layer, positions, len(symbols))
    counts = sparse.csr_matrix((values, indices, pointers), shape=(len(pool), len(symbols)))
    assert before == {'bytes': source.stat().st_size, 'mtime_ns': source.stat().st_mtime_ns}
    source_records.append({'source': str(source), 'layer': '/layers/counts', **before, 'source_opened_read_only': True})
else:
    raw = Path('/home/mazekai/IBD_EcoTyper/06_validation/GSE214695_scRNA_recovery_20261009/source')
    downloads = pd.read_csv(raw.parent / 'download_urls.tsv', sep='\t')
    parts = []
    pool_parts = []
    source_symbols = None
    for sample, block in pool.groupby('Sample', observed=True, sort=True):
        filenames = downloads.loc[downloads.Sample.eq(sample), 'Filename']
        matrix_file = raw / filenames.loc[filenames.str.endswith('_matrix.mtx.gz')].item()
        barcode_file = raw / filenames.loc[filenames.str.endswith('_barcodes.tsv.gz')].item()
        feature_file = raw / filenames.loc[filenames.str.endswith('_features.tsv.gz')].item()
        barcodes = pd.read_csv(barcode_file, sep='\t', header=None).iloc[:, 0].astype(str)
        features = pd.read_csv(feature_file, sep='\t', header=None)
        assert barcodes.is_unique
        current_symbols = features.iloc[:, 1].astype(str)
        if source_symbols is None:
            source_symbols = current_symbols
        assert np.array_equal(source_symbols, current_symbols), '不同样本features必须精确对齐'
        positions = pd.Index(barcodes).get_indexer(block.Barcode)
        assert (positions >= 0).all()
        with gzip.open(matrix_file, 'rb') as handle:
            matrix = io.mmread(handle).tocsc()
        assert matrix.shape == (len(features), len(barcodes))
        selected = matrix[:, positions].T.tocsr()
        parts.append(selected)
        pool_parts.append(block)
        with matrix_file.open('rb') as handle:
            digest = hashlib.file_digest(handle, 'sha256').hexdigest()
        source_records.append({'source': str(matrix_file), 'sha256': digest, 'bytes': matrix_file.stat().st_size, 'source_opened_read_only': True})
        print('Read', sample, selected.shape, flush=True)
        del matrix
    counts = sparse.vstack(parts, format='csr')
    pool = pd.concat(pool_parts, ignore_index=True)
    symbols = source_symbols
    totals = counts.sum(axis=1).A1
    detected = counts.getnnz(axis=1)
assert np.isfinite(counts.data).all() and (counts.data >= 0).all() and np.allclose(counts.data, np.rint(counts.data))
assert (totals > 0).all()
valid_symbols = symbols.str.strip().ne('') & ~symbols.str.contains(r'[\t\r\n\"]')
unique = pd.Index(sorted(symbols.loc[valid_symbols].unique()))
collapse = sparse.csr_matrix((np.ones(valid_symbols.sum()), (np.where(valid_symbols)[0], unique.get_indexer(symbols.loc[valid_symbols]))), shape=(len(symbols), len(unique)))
counts = (counts @ collapse).tocsr()
pool['full_raw_counts'] = totals
pool['detected_source_features'] = detected
mito_genes = np.flatnonzero(unique.str.startswith('MT-'))
pool['mitochondrial_percent'] = counts[:, mito_genes].sum(axis=1).A1 / totals * 100
pool.groupby('CellType', observed=True).agg(cells=('CellID', 'size'), donors=('Donor', 'nunique'), min_depth=('full_raw_counts', 'min'), median_depth=('full_raw_counts', 'median'), median_detected_features=('detected_source_features', 'median'), median_mito_percent=('mitochondrial_percent', 'median')).to_csv(out / 'sampled_type_counts_QC.tsv', sep='\t')
markers = ['MS4A1', 'CD79A', 'MZB1', 'JCHAIN', 'CD3D', 'CD4', 'CD8A', 'NKG7', 'GNLY', 'KLRD1', 'IL7R', 'CSF1R', 'CD68', 'FCER1A', 'CD1C', 'TPSAB1', 'TPSB2', 'KIT', 'CPA3', 'CSF3R', 'S100A8', 'COL1A1', 'DCN', 'RGS5', 'PECAM1', 'VWF', 'EPCAM', 'ALPI', 'FABP6', 'KRT20', 'BEST4', 'OTOP2', 'MUC2', 'TFF3', 'LGR5', 'OLFM4', 'MKI67', 'DEFA5', 'DEFA6', 'CHGA', 'SOX10', 'S100B', 'POU2F3', 'TRPM5']
marker_positions = unique.get_indexer(markers)
marker_present = marker_positions >= 0
cell_cpm_markers = counts[:, marker_positions[marker_present]].multiply((1e6 / totals)[:, None]).toarray()
marker_frame = pd.DataFrame(cell_cpm_markers, columns=np.array(markers)[marker_present])
marker_frame['AuthorLabel'] = pool.AuthorLabel.to_numpy()
marker_frame.groupby('AuthorLabel', observed=True).mean().to_csv(out / 'independent_author_marker_mean_CPM.tsv', sep='\t')
pool.to_csv(out / 'unique_source_cell_manifest.tsv.gz', sep='\t', index=False)
sparse.save_npz(out / 'unique_source_cell_counts_csr.npz', counts)
pd.Series(unique).to_csv(out / 'source_genes.tsv', sep='\t', header=False, index=False)
pd.DataFrame(source_records).to_csv(out / 'input_source_manifest.tsv', sep='\t', index=False)

# Step3. 含重复抽样的counts累加，保存实测细胞数量及全部基因RNA贡献
sample_names = pd.Index(pd.unique(membership['sample']))
cell_position = pd.Index(pool.CellID).get_indexer(membership.CellID)
sample_position = sample_names.get_indexer(membership['sample'])
assert (cell_position >= 0).all()
assignment = sparse.csr_matrix((np.ones(len(membership), dtype=np.int64), (sample_position, cell_position)), shape=(n_mixtures, len(pool)))
mixtures = (assignment @ counts).toarray()
assert np.allclose(mixtures, np.rint(mixtures)) and np.all(mixtures.sum(axis=1) > 0)
expressed = mixtures.sum(axis=0) > 0
mixture_frame = pd.DataFrame(mixtures[:, expressed].T, index=unique[expressed], columns=sample_names)
mixture_frame.to_csv(out / 'pseudo_bulk_raw_counts.tsv.gz', sep='\t', index_label='GeneSymbol')
truth_cells = membership.groupby(['sample', 'CellType'], observed=True).size().unstack(fill_value=0).reindex(index=sample_names, columns=types, fill_value=0) / 1000
truth_cells.to_csv(out / 'truth_cell_count_fraction.tsv', sep='\t', index_label='sample')
cell_library = counts.sum(axis=1).A1
membership['full_gene_RNA_counts'] = cell_library[cell_position]
truth_rna = membership.groupby(['sample', 'CellType'], observed=True).full_gene_RNA_counts.sum().unstack(fill_value=0).reindex(index=sample_names, columns=types, fill_value=0)
truth_rna = truth_rna.div(truth_rna.sum(axis=1), axis=0)
truth_rna.to_csv(out / 'truth_full_gene_RNA_fraction.tsv', sep='\t', index_label='sample')
membership.groupby(['sample', 'Donor'], observed=True).size().rename('draws').reset_index().to_csv(out / 'donor_draw_composition.tsv', sep='\t', index=False)
assert np.allclose(truth_cells.sum(axis=1), 1) and np.allclose(truth_rna.sum(axis=1), 1)
if mode == 'formal':
    cpm = mixture_frame.div(mixtures.sum(axis=1), axis=1) * 1e6
    cpm.to_csv(folder / '05_cibersortx_official_pending/pseudo_bulk100_CPM.txt.gz', sep='\t', index_label='GeneSymbol', float_format='%.8g')
    pending = {'status': 'BLOCKED', 'reason': 'No official Docker/runtime authorization available; no official prediction fabricated', 'signature': 'Unchanged Joint20 v2 Job13', 'mixture': 'pseudo_bulk100_CPM.txt (decompress supplied .gz); same summed counts as BayesPrism before method-specific gene filtering', 'run_mode': 'relative', 'QN': False, 'batch_correction': False, 'permutations': 100, 'scope': 'Impute Cell Fractions only', 'source_sample_ids': sample_names.tolist(), 'paired_predictions_available': False, 'research_SVR_used': False}
    (folder / '05_cibersortx_official_pending/run_parameters.json').write_text(json.dumps(pending, indent=2), encoding='utf-8')
summary = {'dataset': dataset, 'mode': mode, 'seed': seed, 'mixtures': n_mixtures, 'draws_per_mixture': 1000, 'unique_cells': len(pool), 'source_donors': pool.Donor.nunique(), 'genes': int(expressed.sum()), 'counts_nonnegative_integer': True, 'author_labels_used': True, 'EcoTyper_state_labels_used': False, 'source_raw_cells_read': True, 'cell_sampling_with_replacement': True, 'counts_not_normalized_for_BP': True, 'invalid_source_symbols': int((~valid_symbols).sum()), 'duplicate_source_symbol_rows': int(symbols.duplicated(keep=False).sum()), 'full_type_axis': types, 'evaluable_types': active, 'elapsed_seconds': time.monotonic() - started, 'python': sys.executable}
(out / 'generation_summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
print('GENERATION_COMPLETE', json.dumps(summary), flush=True)
