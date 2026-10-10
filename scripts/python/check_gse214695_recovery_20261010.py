# recovery完成后检查原件哈希、输入顺序、状态映射与软丰度；不改变官方输出。
# 在服务器ibd_refaudit环境运行；读本轮source/input/formal/audit，检查记录另存audit。
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

# Step1. 原始下载文件逐个核对；冲突注释的隔离不能改变作者原件
root = Path('/home/mazekai/IBD_EcoTyper/06_validation/GSE214695_scRNA_recovery_20261009')
checksums = pd.read_csv(root / 'source_checksums.tsv', sep='\t')
assert len(checksums) == 55
for row in checksums.itertuples(index=False):
    path = root / 'source' / row.Filename
    assert path.stat().st_size == row.Bytes
    with path.open('rb') as handle:
        assert hashlib.file_digest(handle, 'sha256').hexdigest() == row.SHA256
print('原始55文件哈希一致', flush=True)

# Step2. 基因×细胞CPM、元数据顺序及库大小必须一致
annotation = pd.read_csv(root / 'input/GSE214695_Pilot18_20261009_annotation.tsv', sep='\t')
expression = pd.read_csv(root / 'input/GSE214695_Pilot18_20261009_CPM.txt', sep='\t', index_col=0)
assert expression.index.is_unique and expression.columns.is_unique and annotation.ID.is_unique
assert expression.columns.tolist() == annotation.ID.tolist()
values = expression.to_numpy()
assert np.isfinite(values).all() and (values >= 0).all()
assert (values.sum(axis=0) > 0).all()
summary = json.loads((root / 'input/preparation_summary.json').read_text())
assert expression.shape == (summary['genes'], summary['formal_cells'])
assert annotation.Donor.nunique() == summary['donors']
conflicts = pd.read_csv(root / 'input/quarantined_annotation_conflicts.tsv', sep='\t')
assert not annotation.merge(conflicts[['Sample', 'Barcode']].drop_duplicates(), on=['Sample', 'Barcode'], how='inner').shape[0]
print('输入矩阵和注释对齐', expression.shape, flush=True)
del expression, values

# 诊断marker中哪些基因也用于该类型投影，逐项标记，不能笼统声称独立验证。
diagnostic_columns = pd.read_csv(root / 'input/diagnostic_marker_with_annotation.tsv', sep='\t', index_col=0, nrows=0).columns
diagnostic_genes = [g for g in diagnostic_columns if g not in ['CellType', 'AuthorLabel', 'Donor', 'Stratum']]
model = Path('/home/mazekai/IBD_EcoTyper/software/ecotyper/EcoTyper/MSCCR_CORE12_CTSPEC_REGION_IBD/Cell_type_specific_genes/Cell_States/discovery')
projection_overlap = []
for cell_type in annotation.CellType.unique():
    info_path = list((model / cell_type).glob('*/initial_gene_info.txt'))
    assert len(info_path) == 1
    info = pd.read_csv(info_path[0], sep='\t', index_col=0)
    projection_genes = set(info.Gene)
    projection_overlap.extend([{'CellType': cell_type, 'DiagnosticGene': g, 'UsedInInitialRecoveryGeneList': g in projection_genes} for g in diagnostic_genes])
pd.DataFrame(projection_overlap).to_csv(root / 'audit/diagnostic_marker_projection_overlap.tsv', sep='\t', index=False)

# Step3. 全部细胞及53状态记录必须完整；未分配单独保留，不能补派状态
cells = pd.read_csv(root / 'audit/cell_assignments_all.tsv', sep='\t')
states = pd.read_csv(root / 'audit/state_recovery_summary.tsv', sep='\t')
members = pd.read_csv(root / 'audit/frozen_ecotype_members.tsv', sep='\t')
assert len(cells) == len(annotation) and cells.ID.is_unique and set(cells.ID) == set(annotation.ID)
assert len(states) == 53 and states.CellState.is_unique and set(states.CellState) == set(members.ID)
counted = cells[cells.State.ne('Unassigned')].groupby(['CellType', 'State']).size()
for row in states.itertuples(index=False):
    assert row.AssignedCells == counted.get((row.CellType, row.State), 0)
    assert row.AssignedDonors == cells.loc[cells.CellType.eq(row.CellType) & cells.State.eq(row.State), 'Donor'].nunique()
    assert 0 <= row.AssignedCells <= row.InputCells and 0 <= row.AssignedDonors <= row.InputDonors

# Step4. 官方保留状态软丰度不重归一化；零值/NA按实际输出记录
formal = Path('/home/mazekai/IBD_EcoTyper/software/ecotyper/EcoTyper/MSCCR_CORE12_CTSPEC_REGION_IBD/Cell_type_specific_genes/Cell_States/recovery/scRNA/GSE214695_Pilot18_20261009_CPM')
h_checks = []
for cell_type in annotation.CellType.unique():
    h_path = list((formal / cell_type).glob('*/state_abundances.txt'))
    assert len(h_path) == 1
    h = pd.read_csv(h_path[0], sep='\t', index_col=0)
    assert set(h.columns) == set(annotation.loc[annotation.CellType.eq(cell_type), 'ID'])
    assert set(h.index) == set(states.loc[states.CellType.eq(cell_type), 'State'])
    v = h.to_numpy()
    assert np.isfinite(v).all() and (v >= 0).all() and (v <= 1 + 1e-7).all()
    assert (h.sum(axis=0) <= 1 + 1e-7).all()
    h_checks.append({'CellType': cell_type, 'States': len(h), 'Cells': h.shape[1], 'MedianRetainedHSum': float(h.sum(axis=0).median())})
pd.DataFrame(h_checks).to_csv(root / 'audit/soft_abundance_check.tsv', sep='\t', index=False)
eco = pd.read_csv(root / 'audit/ecotype_abundance.txt', sep='\t', index_col=0)
assert set(eco.columns) == set(annotation.Sample) and set(eco.index) == {f'E{i}' for i in range(1, 8)}
finite = np.isfinite(eco.to_numpy())
assert (eco.to_numpy()[finite] >= 0).all() and (eco.to_numpy()[finite] <= 1 + 1e-7).all()
eco_assignment = pd.read_csv(root / 'audit/ecotype_assignment.txt', sep='\t', index_col=0)
assert eco_assignment.ID.is_unique and set(eco_assignment.ID) <= set(annotation.Sample)
assert (root / 'source_before.sha256').read_bytes() == (root / 'source_after.sha256').read_bytes()
model_check = (root / 'model_integrity.log').read_text()
assert len(model_check.splitlines()) == 650 and all(line.endswith(': OK') for line in model_check.splitlines())
check = {'raw_files_sha256_unchanged': 55, 'input_genes': summary['genes'], 'input_cells': len(annotation), 'donors': annotation.Donor.nunique(), 'types': annotation.CellType.nunique(), 'states': len(states), 'quarantined_rows': len(conflicts), 'quarantined_unique_cells': len(conflicts[['Sample', 'Barcode']].drop_duplicates()), 'matrix_annotation_order_passed': True, 'assignments_counts_passed': True, 'soft_abundance_passed': True, 'ecotype_missing_entries': int((~finite).sum()), 'official_source_files_unchanged': 4, 'frozen_model_files_unchanged': 650, 'scientific_accuracy_claim': False}
(root / 'audit/final_data_check.json').write_text(json.dumps(check, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps(check, ensure_ascii=False), flush=True)
