# 复用 Joint20 v2 的原始训练 counts，提取真实 MSCCR120 并精确恢复自然模拟。
# 在服务器 ibd_refaudit 环境运行；从既有只读检查点取数，输出新的敏感性目录供 R 运行。
from pathlib import Path
import csv
import gzip
import hashlib
import json
import sys
import time
import h5py
import numpy as np
import pandas as pd
from scipy import sparse, io

# Step1. 对齐冻结训练细胞、标签及供者划分
started = time.monotonic()
project = Path('/home/mazekai/IBD_EcoTyper')
out = Path(sys.argv[1])
training = project / '03_reference/joint20_donor_site_v2/output/run_20261007_v2'
validation = project / '03_reference/web_reference_three_pilots/validation_joint20_20261006'
types = (out / '00_audit/frozen_types.txt').read_text().splitlines()
cells = pd.read_csv(training / 'selected_cell_manifest.tsv', sep='\t')
genes = pd.read_csv(training / 'source_gene_manifest.tsv', sep='\t', keep_default_na=False)
roles = pd.read_csv(training / 'frozen_donor_roles.tsv', sep='\t')
counts = sparse.load_npz(training / 'selected_training_counts_csr.npz')
assert counts.shape == (9294, len(genes)) and cells.cell_id.is_unique and cells.obs_name.is_unique
assert set(cells.coarse_identity) == set(types) and len(types) == 20 and cells.Patient.nunique() == 32
assert cells.split.eq('training').all() and not cells.Patient.isin(roles.loc[roles.split.ne('training'), 'Patient']).any()
assert np.isfinite(counts.data).all() and (counts.data >= 0).all() and np.array_equal(counts.data, np.rint(counts.data))
totals = counts.sum(axis=1).A1
assert (totals > 0).all() and np.array_equal(totals, cells.extracted_total_counts)
source = Path(json.loads((training / 'construction_summary.json').read_text())['source'])
# H5AD 源行和条码同时核验；只读元数据与选中稀疏行，避免重新读取全矩阵。
sys.path.insert(0, str(project / 'scripts/python/three_reference_pilots_v1'))
from inspect_ibdverse_h5ad import read_column
from export_ibdverse_cibersortx_pilot import selected_csr
with h5py.File(source, 'r') as h:
    source_ids = read_column(h['obs'][h['obs'].attrs.get('_index', '_index')]).astype(str)
    rows = cells.source_cell_row.to_numpy()
    assert np.array_equal(source_ids.iloc[rows], cells.obs_name)
    for field in ['Patient', 'Site', 'sample_id', 'final_analysis']:
        assert np.array_equal(read_column(h['obs'][field]).iloc[rows].astype(str), cells[field].astype(str)), field
    idx, values, ptr, _, _ = selected_csr(h['X'], rows, len(genes))
    recovered = sparse.csr_matrix((values, idx, ptr), shape=counts.shape)
    assert (recovered != counts).nnz == 0
symbols = genes.gene_symbol.str.strip()
valid = symbols.ne('') & ~symbols.str.contains(r'[\t\r\n\"]')
unique = pd.Index(sorted(symbols[valid].unique()))
mapping = sparse.csr_matrix((np.ones(valid.sum()), (np.where(valid)[0], unique.get_indexer(symbols[valid]))), shape=(len(genes), len(unique)))
reference = (counts @ mapping).tocsr()
keep = reference.sum(axis=0).A1 > 0
reference = reference[:, keep]
reference_genes = unique[keep]
io.mmwrite(out / '01_reference/reference_counts_cells_by_genes.mtx', reference)
pd.Series(reference_genes).to_csv(out / '01_reference/reference_genes.tsv', sep='\t', index=False, header=False)
cells.to_csv(out / '01_reference/reference_cells.tsv', sep='\t', index=False)
qc = cells.groupby('coarse_identity', observed=True).agg(cells=('cell_id', 'size'), donors=('Patient', 'nunique'), min_raw_depth=('extracted_total_counts', 'min'), median_raw_depth=('extracted_total_counts', 'median')).reindex(types)
assert qc.cells.gt(20).all()
qc['counts_nonnegative_integer'] = True
qc['barcode_and_source_counts_exact'] = True
qc['training_only'] = True
qc.to_csv(out / '01_reference/reference_QC.tsv', sep='\t', index_label='cell_type')
print('Reference verified:', reference.shape, 'types:', len(types), flush=True)

# Step2. 从原始 Ensembl counts 提取同120例，沿用已审计的一对一基因身份桥
annotation = pd.read_csv(out / '00_audit/MSCCR120_annotation.tsv', sep='\t', dtype=str)
assert len(annotation) == 120 and annotation.ID.is_unique and annotation.PatientID.nunique() == 120
bridge_path = project / '03_reference/web_reference_three_pilots/msccr_gene_identity_20261006/joint_reference_to_MSCCR_gene_bridge.tsv'
bridge = pd.read_csv(bridge_path, sep='\t', keep_default_na=False)
accepted = bridge.loc[bridge.accepted.eq(True) & bridge.gene_symbol.isin(reference_genes)].copy()
assert accepted.gene_symbol.is_unique and accepted.mapped_msccr_ensembl_id.is_unique
raw_path = project / '01_raw_data/bulk/MSCCR_GSE193677/GSE193677_MSCCR_Biopsy_counts.txt.gz'
raw_position = pd.Series(np.arange(len(accepted)), index=accepted.mapped_msccr_ensembl_id).to_dict()
bulk = np.zeros((len(accepted), 120))
seen = set()
library = np.zeros(120)
with gzip.open(raw_path, 'rt') as handle:
    samples = next(csv.reader([handle.readline().strip()], delimiter=' ', skipinitialspace=True))
    positions = pd.Index(samples).get_indexer(annotation.ID)
    assert (positions >= 0).all() and len(set(samples)) == len(samples)
    for number, line in enumerate(handle, 1):
        fields = line.split(None, 1)
        gene_id = fields[0].strip('"')
        vals = np.fromstring(fields[1], sep=' ')[positions]
        assert len(vals) == 120 and np.isfinite(vals).all() and (vals >= 0).all()
        library += vals
        if gene_id in raw_position:
            assert gene_id not in seen
            bulk[raw_position[gene_id]] = vals
            seen.add(gene_id)
        if number % 10000 == 0:
            print('Raw bulk genes scanned:', number, flush=True)
assert set(raw_position) == seen and (library > 0).all()
pd.DataFrame(bulk, index=accepted.gene_symbol, columns=annotation.ID).to_csv(out / '02_bayesprism_joint20_MSCCR120/mixture_raw_counts.tsv.gz', sep='\t', index_label='gene', compression='gzip')
accepted.to_csv(out / '02_bayesprism_joint20_MSCCR120/bulk_gene_identity.tsv', sep='\t', index=False)
sample_qc = annotation.copy()
sample_qc['full_raw_library_counts'] = library
sample_qc['mapped_raw_counts'] = bulk.sum(axis=0)
sample_qc['fractional_raw_entries'] = (bulk != np.rint(bulk)).sum(axis=0)
sample_qc.to_csv(out / '02_bayesprism_joint20_MSCCR120/sample_QC.tsv', sep='\t', index=False)
# raw 文件存在非整数值，保留来源实值；既往 BayesPrism 也未舍入。
print('Bulk raw counts:', bulk.shape, 'fractional entries:', int((bulk != np.rint(bulk)).sum()), flush=True)

# Step3. 从保存的抽样清单逐一恢复28个自然混合物，校验历史CPM文本精度
membership = pd.read_csv(validation / 'mixture_cell_membership.tsv', sep='\t')
pool = pd.read_csv(validation / 'selected_pool_manifest.tsv', sep='\t')
with h5py.File(validation / 'holdout_pool_counts.h5', 'r') as h:
    pool_counts = sparse.csr_matrix((h['data'][:], h['indices'][:], h['indptr'][:]), shape=tuple(h.attrs['shape']))
    pool_ids = pd.Index(h['cell_id'].asstr()[:])
    pool_genes = pd.Index(h['gene_symbol'].asstr()[:])
assert pool_ids.is_unique and np.array_equal(pool_ids, pool.cell_id)
cell_overlap = len(set(cells.cell_id) & set(membership.cell_id))
donor_overlap = len(set(cells.Patient) & set(membership.Patient))
assert cell_overlap == 0 and donor_overlap == 0
historical = pd.read_csv(validation / 'GJ20_Holdout28_Mixture_CPM.txt', sep='\t', index_col=0)
raw_mixtures = []
filtered_truth_records = []
for sample in historical.columns:
    part = membership.loc[membership['sample'].eq(sample)]
    positions = pool_ids.get_indexer(part.cell_id)
    assert len(part) == 200 and part.cell_id.is_unique and (positions >= 0).all()
    mixture = pool_counts[positions].sum(axis=0).A1
    raw_mixtures.append(mixture)
    reconstructed = mixture[pool_genes.get_indexer(historical.index)] / mixture.sum() * 1e6
    assert np.allclose(reconstructed, historical[sample], rtol=5e-8, atol=5e-5), sample
raw_mixtures = np.array(raw_mixtures)
assert np.array_equal(raw_mixtures, np.rint(raw_mixtures))
pd.DataFrame(raw_mixtures.T, index=pool_genes, columns=historical.columns).to_csv(out / '05_existing_pseudobulk_benchmark/natural28_raw_counts.tsv.gz', sep='\t', index_label='gene', compression='gzip')
for filename in ['ground_truth_export_UMI_fractions.tsv', 'ground_truth_cell_fractions.tsv', 'ground_truth_equalRNA_fractions.tsv', 'mixture_design.tsv', 'paired_equalRNA_samples.tsv']:
    (out / '05_existing_pseudobulk_benchmark' / filename).write_bytes((validation / filename).read_bytes())
membership.to_csv(out / '05_existing_pseudobulk_benchmark/mixture_cell_membership.tsv', sep='\t', index=False)
pool.groupby(['sample_id', 'Patient', 'Site', 'pilot_type'], observed=True).size().rename('cells').reset_index().to_csv(out / '00_audit/validation_pool_provenance.tsv', sep='\t', index=False)
leakage = {'cell_overlap': cell_overlap, 'donor_overlap': donor_overlap, 'dataset_overlap': 'same GSE282122 dataset', 'classification': 'Cell-disjoint; Donor-disjoint; development benchmark; not Dataset-disjoint', 'final_reserved_donors': ['CD12', 'UC17'], 'reserved_expression_used': False, 'seed': 20261006, 'cells_per_mixture': 200, 'natural28_exact_CPM_reconstruction': True, 'equalRNA8_raw_counts_status': 'BLOCKED', 'equalRNA8_reason': 'Historical mixture is mean per-cell CPM with noninteger equal-RNA weights, not a sum of raw counts. No artificial read-depth scaling or rounding was used.'}
assert not membership.Patient.isin(leakage['final_reserved_donors']).any()
(out / '05_existing_pseudobulk_benchmark/reconstruction_audit.json').write_text(json.dumps(leakage, indent=2), encoding='utf-8')

# Step4. 保存输入证据，说明基因桥和两算法的输入差别
inputs = [training / 'selected_training_counts_csr.npz', training / 'selected_cell_manifest.tsv', training / 'source_gene_manifest.tsv', bridge_path, raw_path, validation / 'holdout_pool_counts.h5', validation / 'mixture_cell_membership.tsv']
records = []
for path in inputs:
    with path.open('rb') as handle:
        digest = hashlib.file_digest(handle, 'sha256').hexdigest()
    records.append({'path': str(path), 'bytes': path.stat().st_size, 'sha256': digest})
pd.DataFrame(records).to_csv(out / '00_audit/server_input_manifest.tsv', sep='\t', index=False)
report = {'reference_cells': len(cells), 'reference_genes': len(reference_genes), 'source_genes': len(genes), 'duplicate_symbols': int(symbols.duplicated(keep=False).sum()), 'invalid_symbols': int((~valid).sum()), 'bulk_mapped_genes': len(accepted), 'bulk_samples': 120, 'raw_fractional_entries_preserved': int((bulk != np.rint(bulk)).sum()), 'training_counts_source_verified_exact': True, 'reference_match': 'EXACT_TRAINING_CELLS_AND_TYPES', 'elapsed_seconds': time.monotonic() - started, 'python': sys.executable}
(out / '01_reference/reference_build_audit.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
(out / '01_reference/reference_build_report.md').write_text('# Joint20 v2 reference build\n\nAll 9294 original training cells and all 20 frozen labels recovered exactly from raw /X and the unchanged training checkpoint. 32 donors, 212 pooled signature profiles. No holdout cells/donors or label merging. BayesPrism receives individual-cell counts, whereas Job13 used donor/site pooled CPM profiles. Bulk uses original numeric counts on the existing validated one-to-one gene identity bridge; counts are neither CPM nor rounded. Full source library QC is retained; the bridge is an identity filter, not the 4474-gene signature restriction.\n', encoding='utf-8')
print('PREPARATION_COMPLETE', json.dumps(report), flush=True)
