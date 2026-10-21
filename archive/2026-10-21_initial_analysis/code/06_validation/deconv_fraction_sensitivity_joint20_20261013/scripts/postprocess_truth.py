# 计算BayesPrism实际保留基因集上的UMI真值，保留历史全基因真值作为主要配对评价。
# 服务器ibd_refaudit环境运行；只读冻结counts与成员清单，不改变模拟表达或算法输出。
from pathlib import Path
import json
import sys
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

# Step1. 精确对齐既有抽样成员和本轮过滤基因
out = Path(sys.argv[1])
validation = Path('/home/mazekai/IBD_EcoTyper/03_reference/web_reference_three_pilots/validation_joint20_20261006')
folder = out / '05_existing_pseudobulk_benchmark'
types = (out / '00_audit/frozen_types.txt').read_text().splitlines()
genes = (folder / 'BayesPrism_retained_genes.txt').read_text().splitlines()
members = pd.read_csv(validation / 'mixture_cell_membership.tsv', sep='\t')
pool_metadata = pd.read_csv(validation / 'selected_pool_manifest.tsv', sep='\t')
equal_pairs = pd.read_csv(validation / 'paired_equalRNA_samples.tsv', sep='\t')
pool_metadata['used_Job15'] = pool_metadata.cell_id.isin(members.cell_id)
equal_cells = members.loc[members['sample'].isin(equal_pairs.paired_main_sample), 'cell_id']
pool_metadata['used_Job14'] = pool_metadata.cell_id.isin(equal_cells)
pool_metadata['Dataset'] = 'GSE282122'
pool_metadata['Assignment'] = 'development_validation; donor and cell disjoint from Joint20 training'
pool_metadata['Raw_counts_source'] = str(validation / 'holdout_pool_counts.h5')
pool_metadata['used_signature'] = False
pool_metadata['used_Smode'] = False
pool_metadata.to_csv(out / '00_audit/development_validation_cell_provenance.tsv', sep='\t', index=False)
truth = pd.read_csv(validation / 'ground_truth_export_UMI_fractions.tsv', sep='\t', index_col=0)
with h5py.File(validation / 'holdout_pool_counts.h5', 'r') as h:
    matrix = sparse.csr_matrix((h['data'][:], h['indices'][:], h['indptr'][:]), shape=tuple(h.attrs['shape']))
    source_genes = pd.Index(h['gene_symbol'].asstr()[:])
    ids = pd.Index(h['cell_id'].asstr()[:])
positions = source_genes.get_indexer(genes)
assert (positions >= 0).all()
retained_library = matrix[:, positions].sum(axis=1).A1
records = []
for sample in truth.index:
    part = members.loc[members['sample'].eq(sample)].copy()
    rows = ids.get_indexer(part.cell_id)
    assert (rows >= 0).all()
    part['retained_library'] = retained_library[rows]
    contribution = part.groupby('pilot_type', observed=True).retained_library.sum().reindex(types, fill_value=0)
    records.append(contribution.to_numpy() / contribution.sum())
pd.DataFrame(records, index=truth.index, columns=types).to_csv(folder / 'ground_truth_BP_retained_gene_UMI_fractions.tsv', sep='\t', index_label='sample')

# Step2. 输出完整服务器文件清单，主对象留服务器避免不必要大文件传输
records = []
for path in sorted(out.rglob('*')):
    if path.is_file() and path.name != 'server_output_manifest.tsv':
        records.append({'relative_path': str(path.relative_to(out)), 'bytes': path.stat().st_size, 'server_path': str(path)})
pd.DataFrame(records).to_csv(out / '00_audit/server_output_manifest.tsv', sep='\t', index=False)
print('RETAINED_GENE_TRUTH_COMPLETE', len(genes))
