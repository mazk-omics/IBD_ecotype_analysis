# 审计两个外部单细胞数据的完整计数、作者标签和训练独立性，冻结模拟设计。
# 服务器ibd_refaudit运行；不使用EcoTyper恢复状态作真值，不修改原始数据或正在运行的MSCCR任务。
from pathlib import Path
import gzip
import hashlib
import json
import sys
import h5py
import numpy as np
import pandas as pd

# Step1. 新建外部验证目录，与真实bulk及开发验证分开
project = Path('/home/mazekai/IBD_EcoTyper')
parent = Path(sys.argv[1])
types = (parent / '00_audit/frozen_types.txt').read_text().splitlines()
training = pd.read_csv(parent / '01_reference/reference_cells.tsv', sep='\t')
external = parent / '08_external_pseudobulk'
external.mkdir(exist_ok=False)
for dataset in ['GSE214695', 'IBDverse']:
    folder = external / dataset
    folder.mkdir()
    for subdir in ['00_audit', '01_smoke', '02_simulation', '03_bayesprism', '04_accuracy', '05_cibersortx_official_pending', 'logs']:
        (folder / subdir).mkdir()

# Step2. GSE214695全作者注释；明确标签映射，隔离不能定义Joint20真值的细胞
gse = project / '06_validation/GSE214695_scRNA_recovery_20261009'
raw = gse / 'source'
cells = pd.read_csv(raw / 'GSE214695_cell_annotation.csv.gz').rename(columns={'sample': 'Sample', 'cell_id': 'Barcode', 'annotation': 'AuthorLabel'})
conflicts = cells.duplicated(['Sample', 'Barcode'], keep=False)
cells.loc[conflicts].to_csv(external / 'GSE214695/00_audit/annotation_conflicts.tsv', sep='\t', index=False)
cells = cells.loc[~conflicts].copy()
old_map = pd.read_csv(gse / 'author_label_map.tsv', sep='\t', keep_default_na=False)
mapping = old_map[['AuthorLabel', 'CellType']].copy()
mapping.loc[~mapping.CellType.isin(types), 'CellType'] = ''
mapping.loc[mapping.AuthorLabel.str.startswith(('PC ', 'Plasmablast ')), 'CellType'] = 'Plasma'
mapping.loc[mapping.AuthorLabel.str.contains('colonocyte', case=False) | mapping.AuthorLabel.eq('BEST4 OTOP2'), 'CellType'] = 'Colonic_absorptive'
mapping.loc[mapping.AuthorLabel.eq('Cycling TA'), 'CellType'] = 'Stem_TA_progenitor'
mapping.loc[mapping.AuthorLabel.str.startswith('DCs '), 'CellType'] = 'DC'
mapping.loc[mapping.AuthorLabel.eq('NK'), 'CellType'] = 'NK'
mapping.loc[mapping.AuthorLabel.eq('Tuft cells'), 'CellType'] = 'Tuft'
mapping['mapping_basis'] = 'Author annotation, coarse biological lineage; no EcoTyper state labels'
mapping['limitation'] = ''
mapping.loc[mapping.AuthorLabel.eq('Paneth-like'), 'limitation'] = 'Paneth-like is not a confirmed Paneth identity; excluded from primary truth pool'
mapping.loc[mapping.AuthorLabel.eq('Secretory progenitor'), 'limitation'] = 'No independent allocation to one frozen parent type; excluded'
mapping.loc[mapping.CellType.eq(''), 'status'] = 'EXCLUDED_AMBIGUOUS_OR_OUTSIDE_JOINT20'
mapping.loc[mapping.CellType.ne(''), 'status'] = 'AUTHOR_MAPPED'
cells = cells.merge(mapping, on='AuthorLabel', validate='many_to_one')
metadata = pd.read_csv(gse / 'sample_metadata.tsv', sep='\t')
cells = cells.merge(metadata, on='Sample', validate='many_to_one')
cells['CellID'] = 'GSE214695|' + cells.Sample + '|' + cells.Barcode
cells['Donor'] = 'GSE214695|' + cells.Sample
cells['SourceRow'] = -1
mapping.to_csv(external / 'GSE214695/00_audit/author_type_mapping.tsv', sep='\t', index=False)
cells.to_csv(external / 'GSE214695/00_audit/all_author_cells.tsv.gz', sep='\t', index=False)
cells.loc[cells.CellType.ne('')].to_csv(external / 'GSE214695/00_audit/eligible_author_cells.tsv.gz', sep='\t', index=False)
downloads = pd.read_csv(gse / 'download_urls.tsv', sep='\t')
checks = []
for row in downloads.itertuples(index=False):
    path = raw / row.Filename
    assert path.exists() and path.stat().st_size > 0
    if path.name.endswith('_matrix.mtx.gz'):
        with gzip.open(path, 'rt') as handle:
            line = handle.readline()
            while line.startswith('%'):
                line = handle.readline()
            shape = list(map(int, line.split()))
        checks.append({'Sample': row.Sample, 'matrix_path': str(path), 'features': shape[0], 'barcodes': shape[1], 'nnz': shape[2], 'bytes': path.stat().st_size})
pd.DataFrame(checks).to_csv(external / 'GSE214695/00_audit/raw_matrix_inventory.tsv', sep='\t', index=False)
assert len(checks) == 18 and cells.Sample.nunique() == 18

# Step3. IBDverse完整回肠cohort，关联已审计供者表，不读取既往recovery预测
source = Path('/home/mazekai/enteric_glia/data_scRNA/scIBD_external/IBDverse_E-MTAB-16998/E-MTAB-16999/Full_cohort.h5ad')
sys.path.insert(0, str(project / 'scripts/python/ibdverse_pilot_v1'))
from inspect_ibdverse_h5ad import read_column
coverage = pd.read_csv(project / 'logs/ibdverse_reference_inputs_20261006_161527/label_category_sample_donor_coverage.tsv', sep='\t')
meta_columns = ['sanger_sample_id', 'individual_id', 'disease_status', 'ti_ses_inflamed', 'biopsy_type', 'reference_metadata_complete']
metadata = coverage[meta_columns].drop_duplicates()
assert metadata.sanger_sample_id.is_unique
with h5py.File(source, 'r') as h:
    obs = h['obs']
    ibd = pd.DataFrame({'CellID': read_column(obs[obs.attrs.get('_index', '_index')]).astype(str), 'Sample': read_column(obs['sanger_sample_id']).astype(str), 'AuthorLabel': read_column(obs['label']).astype(str), 'AuthorCategory': read_column(obs['category']).astype(str), 'TotalCountsMetadata': read_column(obs['total_counts']), 'DetectedGenesMetadata': read_column(obs['n_genes_by_counts']), 'MitoPercentMetadata': read_column(obs['pct_counts_gene_group__mito_transcript'])})
    ibd['SourceRow'] = np.arange(len(ibd))
    assert tuple(h['layers/counts'].attrs['shape']) == (len(ibd), len(read_column(h['var']['gene_symbols'])))
    assert h['layers/counts'].attrs['encoding-type'] == 'csr_matrix'
    first_data = h['layers/counts/data'][:100000]
    assert np.isfinite(first_data).all() and (first_data >= 0).all() and np.allclose(first_data, np.rint(first_data))
ibd = ibd.merge(metadata, left_on='Sample', right_on='sanger_sample_id', how='left', validate='many_to_one')
ibd['Donor'] = ibd.individual_id
label_map = pd.DataFrame({'AuthorLabel': sorted(ibd.AuthorLabel.unique()), 'CellType': ''})
label_map.loc[label_map.AuthorLabel.str.startswith('B cells'), 'CellType'] = 'B'
label_map.loc[label_map.AuthorLabel.str.startswith('Plasma B cell '), 'CellType'] = 'Plasma'
label_map.loc[label_map.AuthorLabel.str.startswith('T cells CD4+'), 'CellType'] = 'CD4_T'
label_map.loc[label_map.AuthorLabel.str.startswith('T cells CD8+'), 'CellType'] = 'CD8_T'
label_map.loc[label_map.AuthorLabel.str.startswith(('Macrophages ', 'Monocytes ')), 'CellType'] = 'Monocyte_Macrophage'
label_map.loc[label_map.AuthorLabel.str.startswith(('Stem cells ', 'Enterocytes progenitor crypt ')), 'CellType'] = 'Stem_TA_progenitor'
label_map.loc[label_map.AuthorLabel.str.startswith('Enterocytes ') & label_map.CellType.eq(''), 'CellType'] = 'Ileal_absorptive'
label_map.loc[label_map.AuthorLabel.str.startswith('Goblet cells '), 'CellType'] = 'Goblet'
for label, cell_type in {'ILCs': 'ILC', 'Mast cells': 'Mast', 'cDC1': 'DC', 'pDCs/cDC2': 'DC', 'Fibroblasts': 'Fibroblast', 'Pericyte cells': 'Pericyte', 'Endothelial cells': 'Endothelial', 'Tuft cells': 'Tuft', 'Enteroendocrine cells': 'Enteroendocrine', 'Enterochromaffin cells': 'Enteroendocrine'}.items():
    label_map.loc[label_map.AuthorLabel.eq(label), 'CellType'] = cell_type
label_map.loc[label_map.AuthorLabel.eq('B cells germinal centre/plasmablasts'), 'CellType'] = ''
label_map['mapping_basis'] = 'Author annotation; explicit broad type mapping, independent of inferred fractions/states'
label_map['status'] = np.where(label_map.CellType.ne(''), 'AUTHOR_MAPPED', 'EXCLUDED_AMBIGUOUS_OR_OUTSIDE_JOINT20')
label_map['limitation'] = ''
label_map.loc[label_map.AuthorLabel.eq('T cells CD8+ TRGC2+'), 'limitation'] = 'Author CD8 label retained; gamma-delta-associated subtype noted'
label_map.loc[label_map.AuthorLabel.eq('Enterocytes BEST4+'), 'limitation'] = 'Coarse regional absorptive mapping; subtype differs from source reference composition'
ibd = ibd.merge(label_map, on='AuthorLabel', validate='many_to_one')
ibd['metadata_eligible'] = ibd.reference_metadata_complete.astype(str).str.lower().eq('true') & ibd.biopsy_type.eq('terminal ileum') & ibd.Donor.notna()
ibd['counts_QC_eligible'] = np.isfinite(ibd.TotalCountsMetadata) & ibd.TotalCountsMetadata.gt(0) & ibd.DetectedGenesMetadata.gt(0)
ibd['eligible'] = ibd.metadata_eligible & ibd.counts_QC_eligible & ibd.CellType.ne('')
assert ibd.CellID.is_unique
label_map.to_csv(external / 'IBDverse/00_audit/author_type_mapping.tsv', sep='\t', index=False)
ibd.to_csv(external / 'IBDverse/00_audit/all_author_cells.tsv.gz', sep='\t', index=False)
ibd.loc[ibd.eligible].to_csv(external / 'IBDverse/00_audit/eligible_author_cells.tsv.gz', sep='\t', index=False)

# Step4. 覆盖和独立性分层；不以不同供者命名空间代替患者级身份核验
for dataset, frame in [('GSE214695', cells.loc[cells.CellType.ne('')]), ('IBDverse', ibd.loc[ibd.eligible])]:
    folder = external / dataset
    capacity = frame.groupby('CellType', observed=True).agg(cells=('CellID', 'size'), donors=('Donor', 'nunique'), samples=('Sample', 'nunique')).reindex(types, fill_value=0)
    capacity['evaluable'] = capacity.cells.gt(0)
    capacity['unsupported_reason'] = np.where(capacity.evaluable, '', 'No unambiguous author-mapped cells in the simulated pool; predicted mass retained as absent-type false-positive control')
    capacity.to_csv(folder / '00_audit/Full20_type_coverage.tsv', sep='\t', index_label='cell_type')
    frame.groupby(['Donor', 'CellType'], observed=True).size().rename('cells').reset_index().to_csv(folder / '00_audit/donor_type_counts.tsv', sep='\t', index=False)
    assert not set(frame.CellID) & set(training.cell_id)
    source_note = 'GSE214695 public original feature-barcode matrices and author annotation; GEO describes 6HC+6UC+6CD independent subjects' if dataset == 'GSE214695' else 'IBDverse E-MTAB-16999 terminal-ileal cohort, separate accession/source file from TAURUS/GSE282122'
    audit = {'dataset': dataset, 'source': source_note, 'training_source': 'GSE282122 only, exactly 9294 frozen training cells', 'source_dataset_disjoint': True, 'exact_cell_ID_intersection': 0, 'sample_ID_intersection': len(set(frame.Sample) & set(training.sample_id)), 'donor_ID_intersection_same_namespace_not_available': True, 'cross_study_patient_linkage': 'UNKNOWN: patient pseudonyms cannot establish a biological donor intersection; distinct recruitment/source context supports external-cohort use', 'independence_claim': 'Dataset-external; no training counts or inferred state labels used. Biological cross-study donor identity not independently verifiable.', 'previous_project_use': 'Previously inspected for reference pilots/recovery; not untouched prospective validation. Not used in frozen Joint20 Job13 or this BayesPrism training reference.', 'eligible_cells': len(frame), 'eligible_donors': frame.Donor.nunique(), 'evaluable_types': capacity.index[capacity.evaluable].tolist(), 'unsupported_types': capacity.index[~capacity.evaluable].tolist(), 'official_CIBERSORTx': 'BLOCKED: server docker executable absent; official runtime token not provided; no official external-mixture predictions exist; research SVR prohibited', 'author_QC': 'Published author-filtered annotations and original count layer; no new reclustering, doublet caller or recovery-based truth assignment', 'new_QC': 'finite, nonnegative integer counts; positive cell depth; barcode-feature alignment; source totals consistency where available'}
    (folder / '00_audit/source_independence_audit.json').write_text(json.dumps(audit, indent=2), encoding='utf-8')
    # 源数据供者组成替代原文癌种bulk均值/SD，这是预先定义的IBD适配。
    donor_counts = frame.groupby(['Donor', 'CellType'], observed=True).size().unstack(fill_value=0)
    donor_fraction = donor_counts.div(donor_counts.sum(axis=1), axis=0)
    parameters = pd.DataFrame({'cell_type': donor_fraction.columns, 'mu': donor_fraction.mean(), 'observed_donor_SD': donor_fraction.std(ddof=1)})
    parameters['gaussian_SD'] = np.maximum(0.25, 3 * parameters.observed_donor_SD)
    parameters.to_csv(folder / '00_audit/gaussian_sampling_parameters.tsv', sep='\t', index=False)
    prespecified = {'main_mixtures': 100, 'cells_per_mixture': 1000, 'smoke_mixtures': 3, 'seed_smoke': 20261013, 'seed_formal': 20261014 if dataset == 'GSE214695' else 20261015, 'fraction_sampling': 'Gaussian per covered type using donor-composition mean and SD=max(0.25,3*observed donor SD); clip negatives to zero and renormalize; all-zero vector redrawn', 'cell_sampling': 'multinomial integer allocations, then original annotated cells sampled within type with replacement; 1000 draws exactly', 'aggregation': 'sum original molecule counts in linear space; duplicate symbols summed, no CPM used as BayesPrism counts', 'primary_accuracy_endpoint': 'Final theta MAE on author-evaluable types against RNA-contribution truth on exact BayesPrism retained genes, without renormalizing predictions', 'secondary_endpoints': ['cell-count truth, reproducing the core known-composition goal', 'full-count-gene RNA truth', 'Full20 error and unsupported-type predicted mass', 'initial theta'], 'independent_replicates': '100 designed mixtures reuse source cells and donors; not 100 independent biological donors', 'reference': str(parent / '01_reference/clean_reference.rds'), 'reference_types_changed': False}
    (folder / '00_audit/endpoint_prespecification.json').write_text(json.dumps(prespecified, indent=2), encoding='utf-8')
    (folder / 'execution_status.json').write_text(json.dumps({'audit': 'COMPLETE', 'smoke': 'NOT_RUN', 'generation': 'NOT_RUN', 'BayesPrism': 'NOT_RUN', 'accuracy': 'NOT_RUN', 'official_CIBERSORTx_paired': 'BLOCKED'}, indent=2), encoding='utf-8')
    print(dataset, len(frame), 'cells;', frame.Donor.nunique(), 'donors;', int(capacity.evaluable.sum()), 'types', flush=True)
print('EXTERNAL_AUDIT_COMPLETE')
