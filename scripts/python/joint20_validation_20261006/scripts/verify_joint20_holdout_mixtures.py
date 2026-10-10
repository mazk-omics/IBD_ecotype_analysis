# 独立重算留出混合矩阵和真实比例，核对供者隔离、行列名称及网页格式。
# 在 ibd_refaudit 环境运行；读入生成目录中的 CSR 检查点、细胞清单及导出文件。
# 本脚本只审计，不运行 CIBERSORTx，也不把模拟样本当成独立供者。
import argparse
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

# Step1. 读入生成结果和原始抽样身份，按名称匹配
parser = argparse.ArgumentParser()
parser.add_argument('--input-dir', required=True)
args = parser.parse_args()
p = Path(args.input_dir)
summary = json.loads((p / 'summary.json').read_text())
pool = pd.read_csv(p / 'selected_pool_manifest.tsv', sep='\t')
members = pd.read_csv(p / 'mixture_cell_membership.tsv', sep='\t')
roles = pd.read_csv(p / 'donor_roles.tsv', sep='\t')
paired = pd.read_csv(p / 'paired_equalRNA_samples.tsv', sep='\t')
main = pd.read_csv(p / 'GJ20_Holdout28_Mixture_CPM.txt', sep='\t', index_col=0)
equal = pd.read_csv(p / 'GJ20_Holdout8_EqualRNA_Mixture_CPM.txt', sep='\t', index_col=0)
truth = pd.read_csv(p / 'ground_truth_cell_fractions.tsv', sep='\t', index_col=0)
rna = pd.read_csv(p / 'ground_truth_export_UMI_fractions.tsv', sep='\t', index_col=0)
truth_equal = pd.read_csv(p / 'ground_truth_equalRNA_fractions.tsv', sep='\t', index_col=0)
assert pool.cell_id.is_unique and main.index.is_unique and equal.index.is_unique
assert not members.duplicated(['sample', 'cell_id']).any()
assert set(pool.Patient) <= set(roles.loc[roles.validation_role.eq('pilot_validation'), 'Patient'])
training = pd.read_csv('/home/mazekai/IBD_EcoTyper/03_reference/web_reference_three_pilots/output/run_20261006_183211/gse_joint/selected_cell_manifest.tsv', sep='\t')
assert set(pool.Patient).isdisjoint(set(training.individual_id))
assert set(pool.cell_id).isdisjoint(set(training.cell_id))
assert set(pool.Patient).isdisjoint(set(summary['reserved_final_test_donors']))
with h5py.File(p / 'holdout_pool_counts.h5', 'r') as h:
    counts = csr_matrix((h['data'][:], h['indices'][:], h['indptr'][:]), shape=tuple(h.attrs['shape']))
    genes = h['gene_symbol'].asstr()[:]
    cell_ids = h['cell_id'].asstr()[:]
assert np.array_equal(cell_ids, pool.cell_id.astype(str).to_numpy())
library = np.asarray(counts.sum(axis=1)).ravel()
assert np.allclose(library, pool.export_library_counts)
assert np.array_equal(main.columns, truth.index) and np.array_equal(equal.columns, truth_equal.index)
gene_positions = pd.Index(genes).get_indexer(main.index)
assert (gene_positions >= 0).all() and np.array_equal(main.index, equal.index)

# Step2. 由细胞成员重新汇总表达，不使用生成脚本中的混合算法
maximum_error = 0.0
for sample in main.columns:
    selected_ids = members.loc[members['sample'].eq(sample), 'cell_id']
    rows = pd.Index(cell_ids).get_indexer(selected_ids)
    assert len(rows) == 200 and (rows >= 0).all() and len(np.unique(rows)) == 200
    expected_counts = np.asarray(counts[rows].sum(axis=0)).ravel()
    expected_cpm = expected_counts / expected_counts.sum() * 1e6
    difference = np.abs(expected_cpm[gene_positions] - main[sample].to_numpy())
    maximum_error = max(maximum_error, float(difference.max()))
    assert np.allclose(expected_cpm[gene_positions], main[sample], rtol=6e-8, atol=1e-7)
    expected_truth = pool.loc[rows, 'pilot_type'].value_counts().reindex(truth.columns, fill_value=0) / 200
    assert np.allclose(expected_truth, truth.loc[sample])
    expected_rna = pool.loc[rows].groupby('pilot_type', observed=True).export_library_counts.sum().reindex(truth.columns, fill_value=0)
    assert np.allclose(expected_rna / expected_rna.sum(), rna.loc[sample])
for record in paired.itertuples(index=False):
    selected_ids = members.loc[members['sample'].eq(record.paired_main_sample), 'cell_id']
    rows = pd.Index(cell_ids).get_indexer(selected_ids)
    expected = np.asarray(counts[rows].multiply((1e6 / library[rows])[:, None]).mean(axis=0)).ravel()
    assert np.allclose(expected[gene_positions], equal[record.sample], rtol=6e-8, atol=1e-7)
    assert np.allclose(truth.loc[record.paired_main_sample], truth_equal.loc[record.sample])

# Step3. 核对网页格式、真实比例总和及源文件摘要
for name in ['GJ20_Holdout28_Mixture_CPM.txt', 'GJ20_Holdout8_EqualRNA_Mixture_CPM.txt']:
    text = (p / name).read_text()
    assert '"' not in text and '\t' in text.splitlines()[0]
    assert hashlib.sha256((p / name).read_bytes()).hexdigest() == summary['files'][name]['sha256']
for matrix in [main, equal]:
    assert np.isfinite(matrix.to_numpy()).all() and (matrix.to_numpy() >= 0).all()
    assert matrix.size < 50_000_000 and matrix.to_numpy().max() >= 50
    assert np.allclose(matrix.sum(axis=0), 1e6, rtol=1e-7)
for proportions in [truth, rna, truth_equal]:
    assert np.allclose(proportions.sum(axis=1), 1) and (proportions.to_numpy() >= 0).all()
report = {'status': 'passed', 'main_shape': list(main.shape), 'equalRNA_shape': list(equal.shape), 'all_36_profiles_recomputed': True,
    'all_ground_truth_fractions_recomputed': True, 'no_duplicate_cells_within_sample': True, 'no_training_cell_or_donor_overlap': True,
    'no_reserved_final_test_donors_used': True, 'web_format_valid': True, 'max_absolute_CPM_roundtrip_error': maximum_error}
(p / 'independent_validation.json').write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps(report, ensure_ascii=False, indent=2))
