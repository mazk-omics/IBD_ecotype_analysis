# 复用已审计Joint20参考并从原始MSCCR counts提取完全匹配的CD+UC样本。
# 服务器ibd_refaudit环境运行；所有输出仅写本次独立目录，审计失败即停止。
from pathlib import Path
import csv
import gzip
import hashlib
import json
import shutil
import sys
import time
import numpy as np
import pandas as pd
from scipy import sparse, io

# Step1. 核验原始输入哈希和已审计9294训练细胞
started = time.monotonic()
project = Path('/home/mazekai/IBD_EcoTyper')
out = Path(sys.argv[1])
benchmark = project / '06_validation/deconv_fraction_sensitivity_joint20_20261013'
training = project / '03_reference/joint20_donor_site_v2/output/run_20261007_v2'
web = project / '03_reference/web_reference_three_pilots/msccr_ibd_joint_web_20261006'
local = json.loads((out / '00_audit/local_mainline_evidence.json').read_text())
manifest = pd.read_csv(benchmark / '00_audit/server_input_manifest.tsv', sep='\t')
records = []
for row in manifest.iloc[:5].itertuples():
    path = Path(row.path)
    with path.open('rb') as handle:
        digest = hashlib.file_digest(handle, 'sha256').hexdigest()
    assert digest == row.sha256, f'BLOCKER: source changed {path}'
    records.append({'path': str(path), 'bytes': path.stat().st_size, 'sha256': digest, 'role': 'unchanged_benchmark_source'})
cells = pd.read_csv(benchmark / '01_reference/reference_cells.tsv', sep='\t')
original_cells = pd.read_csv(training / 'selected_cell_manifest.tsv', sep='\t')
pd.testing.assert_frame_equal(cells, original_cells)
types = (benchmark / '00_audit/frozen_types.txt').read_text().splitlines()
assert types == local['types'] and len(types) == 20
assert len(cells) == 9294 and cells.cell_id.is_unique and cells.obs_name.is_unique
assert cells.Patient.nunique() == 32 and cells.split.eq('training').all()
assert set(cells.coarse_identity) == set(types)
genes = pd.Index((benchmark / '01_reference/reference_genes.tsv').read_text().splitlines())
assert genes.is_unique and len(genes) == 25703
source_genes = pd.read_csv(training / 'source_gene_manifest.tsv', sep='\t', keep_default_na=False)
symbols = source_genes.gene_symbol.str.strip()
assert symbols.is_unique and symbols.ne('').all()
counts = sparse.load_npz(training / 'selected_training_counts_csr.npz')
reference = io.mmread(benchmark / '01_reference/reference_counts_cells_by_genes.mtx').tocsr()
assert reference.shape == (9294, len(genes)) and counts.shape == (9294, len(symbols))
positions = pd.Index(symbols).get_indexer(genes)
assert (positions >= 0).all() and (reference != counts[:, positions]).nnz == 0
assert np.isfinite(reference.data).all() and (reference.data >= 0).all() and np.array_equal(reference.data, np.rint(reference.data))
assert np.array_equal(counts.sum(axis=1).A1, cells.extracted_total_counts)
for filename in ['reference_cells.tsv', 'reference_genes.tsv', 'reference_QC.tsv', 'reference_build_audit.json', 'gene_cleanup_audit.tsv', 'clean_reference.rds']:
    source = benchmark / '01_reference' / filename
    shutil.copyfile(source, out / '01_reference' / filename)
    with source.open('rb') as handle:
        digest = hashlib.file_digest(handle, 'sha256').hexdigest()
    records.append({'path': str(source), 'bytes': source.stat().st_size, 'sha256': digest, 'role': 'reused_audited_reference'})
(out / '00_audit/frozen_types.txt').write_text('\n'.join(types) + '\n')
print('REFERENCE_PASS: 9294 exact cells, 32 donors, 20 exact types, raw counts exact', flush=True)

# Step2. 核对主线实际注释、完整Job16样本轴和CPM轴
annotation_path = web / 'MSCCR_CD_UC_2025_annotation.tsv'
assert hashlib.sha256(annotation_path.read_bytes()).hexdigest() == local['metadata_sha256']
annotation = pd.read_csv(annotation_path, sep='\t', dtype=str, keep_default_na=False)
ids = annotation.ID.tolist()
assert hashlib.sha256(('\n'.join(ids) + '\n').encode()).hexdigest() == local['ordered_sample_sha256']
assert len(ids) == local['samples'] and annotation.ID.is_unique
assert set(annotation.Disease) == {'CD', 'UC'}
assert annotation[['ID', 'PatientID', 'Disease', 'Inflammation', 'Region3']].ne('').all().all()
with (web / 'MSCCR_CD_UC_2025_GJ20_CPM.txt').open() as handle:
    cpm_header = handle.readline().rstrip('\n').split('\t')[1:]
assert ids == cpm_header
legacy = pd.read_csv(project / '05_EcoTyper/input/MSCCR_CORE12/annotation/MSCCR_CORE12_annotation.tsv', sep='\t', dtype=str, keep_default_na=False)
assert legacy.ID.is_unique
aligned = legacy.set_index('ID').loc[ids]
assert np.array_equal(aligned[['PatientID', 'Disease', 'Inflammation', 'Region3']].to_numpy(), annotation[['PatientID', 'Disease', 'Inflammation', 'Region3']].to_numpy())
assert set(legacy.loc[legacy.Disease.isin(['CD', 'UC']), 'ID']) == set(ids)
full_metadata_path = project / '03_reference/combined_reference/output/21_bayesprism_epithelial_qc/GSE193677_metadata_used_for_stage2_qc.tsv'
full_metadata = pd.read_csv(full_metadata_path, sep='\t', dtype=str, keep_default_na=False)
assert full_metadata.sample_id.is_unique
technical_path = project / '03_reference/combined_reference/output/26_panibd_panel_classification/technical_exclusion_samples.tsv'
technical = pd.read_csv(technical_path, sep='\t', dtype=str, keep_default_na=False)
assert technical.sample_id.is_unique
full_ibd = set(full_metadata.loc[full_metadata.disease.isin(['CD', 'UC']), 'sample_id'])
assert full_ibd - set(technical.sample_id) == set(ids)
assert full_ibd - set(ids) == set(technical.sample_id)
annotation.insert(0, 'frozen_order', np.arange(1, len(ids) + 1))
annotation.to_csv(out / '00_audit/MSCCR_CD_UC_discovery_manifest.tsv', sep='\t', index=False)
print('SAMPLES_PASS:', len(ids), annotation.Disease.value_counts().to_dict(), 'patients', annotation.PatientID.nunique(), flush=True)

# Step3. 沿用已审计一对一基因桥，从原始文件保留实数counts和完整库分母
bridge_path = project / '03_reference/web_reference_three_pilots/msccr_gene_identity_20261006/joint_reference_to_MSCCR_gene_bridge.tsv'
bridge = pd.read_csv(bridge_path, sep='\t', keep_default_na=False)
accepted = bridge.loc[bridge.accepted.eq(True) & bridge.gene_symbol.isin(genes)].copy()
assert accepted.gene_symbol.is_unique and accepted.mapped_msccr_ensembl_id.is_unique
raw_path = project / '01_raw_data/bulk/MSCCR_GSE193677/GSE193677_MSCCR_Biopsy_counts.txt.gz'
raw_position = dict(zip(accepted.mapped_msccr_ensembl_id, range(len(accepted))))
bulk = np.zeros((len(accepted), len(ids)))
library = np.zeros(len(ids))
full_fractional = np.zeros(len(ids), dtype=int)
seen = set()
with gzip.open(raw_path, 'rt') as handle:
    samples = next(csv.reader([handle.readline().strip()], delimiter=' ', skipinitialspace=True))
    assert len(samples) == len(set(samples))
    positions = pd.Index(samples).get_indexer(ids)
    assert (positions >= 0).all()
    for number, line in enumerate(handle, 1):
        fields = line.split(None, 1)
        gene_id = fields[0].strip('"')
        all_values = np.fromstring(fields[1], sep=' ')
        assert len(all_values) == len(samples)
        vals = all_values[positions]
        assert np.isfinite(vals).all() and (vals >= 0).all()
        library += vals
        full_fractional += vals != np.rint(vals)
        if gene_id in raw_position:
            assert gene_id not in seen
            bulk[raw_position[gene_id]] = vals
            seen.add(gene_id)
        if number % 10000 == 0:
            print('Raw bulk genes scanned:', number, flush=True)
assert set(raw_position) == seen and (library > 0).all()
excluded = pd.DataFrame({'ID': [sample for sample in samples if sample not in set(ids)]})
excluded = excluded.merge(full_metadata[['sample_id', 'disease']].rename(columns={'sample_id': 'ID', 'disease': 'Disease'}), on='ID', how='left', validate='one_to_one')
assert excluded.Disease.notna().all()
excluded['reason'] = np.where(excluded.Disease.isin(['CD', 'UC']), 'Frozen mainline technical exclusion; see technical_exclusion_samples.tsv', 'Non-IBD Control excluded from matched discovery cohort')
assert set(excluded.loc[excluded.Disease.isin(['CD', 'UC']), 'ID']) == set(technical.sample_id)
excluded['also_excluded_from_legacy_2486_technical_panel'] = ~excluded.ID.isin(legacy.ID)
excluded.to_csv(out / '00_audit/excluded_bulk_samples.tsv', sep='\t', index=False)
shutil.copyfile(technical_path, out / '00_audit/technical_exclusion_samples.tsv')
assert set(samples) == set(full_metadata.sample_id)
pd.DataFrame(bulk, index=accepted.gene_symbol, columns=ids).to_csv(out / '02_bayesprism/mixture_raw_counts.tsv.gz', sep='\t', index_label='gene', compression='gzip')
accepted.to_csv(out / '02_bayesprism/bulk_gene_identity.tsv', sep='\t', index=False)
qc = annotation.copy()
qc['full_raw_library_counts'] = library
qc['mapped_raw_counts'] = bulk.sum(axis=0)
qc['full_raw_fractional_entries'] = full_fractional
qc['mapped_fractional_entries'] = (bulk != np.rint(bulk)).sum(axis=0)
qc['pre_prism_truncation_loss_upper_bound'] = (bulk - np.floor(bulk)).sum(axis=0) / bulk.sum(axis=0)
qc.to_csv(out / '00_audit/bulk_sample_QC.tsv', sep='\t', index=False)

# Step4. 与主线CPM独立回算核对，同源benchmark120的counts逐值核对
max_cpm_error = 0.0
for block in pd.read_csv(web / 'MSCCR_CD_UC_2025_GJ20_CPM.txt', sep='\t', index_col=0, chunksize=500):
    shared = block.index.intersection(accepted.gene_symbol)
    idx = pd.Index(accepted.gene_symbol).get_indexer(shared)
    reconstructed = bulk[idx] / library * 1e6
    max_cpm_error = max(max_cpm_error, float(np.abs(reconstructed - block.loc[shared].to_numpy()).max(initial=0)))
assert max_cpm_error <= 0.000501
old = pd.read_csv(benchmark / '02_bayesprism_joint20_MSCCR120/mixture_raw_counts.tsv.gz', sep='\t', index_col=0)
assert list(old.index) == accepted.gene_symbol.tolist()
assert np.array_equal(old.to_numpy(), bulk[:, pd.Index(ids).get_indexer(old.columns)])
for path in [annotation_path, full_metadata_path, technical_path, web / 'MSCCR_CD_UC_2025_GJ20_CPM.txt', out / '02_bayesprism/mixture_raw_counts.tsv.gz']:
    with path.open('rb') as handle:
        digest = hashlib.file_digest(handle, 'sha256').hexdigest()
    records.append({'path': str(path), 'bytes': path.stat().st_size, 'sha256': digest, 'role': 'matched_bulk_or_manifest'})
pd.DataFrame(records).to_csv(out / '00_audit/input_SHA256.tsv', sep='\t', index=False)
report = {'status': 'PASS_WITH_KNOWN_FRACTIONAL_COUNTS_LIMITATION', 'samples': len(ids), 'disease_counts': annotation.Disease.value_counts().to_dict(), 'inflammation_counts': annotation.Inflammation.value_counts().to_dict(), 'region_counts': annotation.Region3.value_counts().to_dict(), 'patients': int(annotation.PatientID.nunique()), 'excluded_nonIBD': int(excluded.Disease.eq('Control').sum()), 'frozen_technical_IBD_exclusions': len(technical), 'reference_cells': len(cells), 'reference_donors': int(cells.Patient.nunique()), 'types': types, 'reference_genes': len(genes), 'bulk_mapped_genes': len(accepted), 'ordered_sample_sha256': local['ordered_sample_sha256'], 'metadata_local_server_byte_identical': True, 'Job16_order_exact': True, 'mainline_CPM_roundtrip_max_abs_error': max_cpm_error, 'benchmark120_counts_exact': True, 'mapped_fractional_entries': int(qc.mapped_fractional_entries.sum()), 'fractional_samples': int(qc.mapped_fractional_entries.gt(0).sum()), 'max_preprism_truncation_loss_bound_fraction': float(qc.pre_prism_truncation_loss_upper_bound.max()), 'seed': 20261013, 'ncores': 8, 'elapsed_seconds': time.monotonic() - started}
(out / '00_audit/input_audit.json').write_text(json.dumps(report, indent=2) + '\n')
print('INPUT_AUDIT_PASS', json.dumps(report), flush=True)
