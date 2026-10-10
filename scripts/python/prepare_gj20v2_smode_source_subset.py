# 从v2已保存的训练细胞counts检查点导出网页可用的同来源单细胞校正诊断参考，不重新抽取其他队列。
# 服务器ibd_refaudit环境运行；读v2目录的稀疏counts、细胞/基因清单及原CPM，写1060细胞左右的非log CPM与来源记录。
# v2 signature由212个汇总参考构建，本文件是这些参考的来源子集；不能宣称它是作者指南要求的完全相同构建输入。
from pathlib import Path
import argparse
import hashlib
import json
import sys
import numpy as np
import pandas as pd
from scipy import sparse

# Step1. 读取原v2检查点，细胞和基因顺序沿用当时保存的清单
parser = argparse.ArgumentParser()
parser.add_argument('--v2-dir', required=True)
parser.add_argument('--outdir', required=True)
args = parser.parse_args()
source = Path(args.v2_dir)
out = Path(args.outdir)
out.mkdir(parents=True, exist_ok=False)
cells = pd.read_csv(source / 'selected_cell_manifest.tsv', sep='\t')
genes = pd.read_csv(source / 'source_gene_manifest.tsv', sep='\t', dtype=str)
profiles = pd.read_csv(source / 'profile_manifest.tsv', sep='\t')
roles = pd.read_csv(source / 'frozen_donor_roles.tsv', sep='\t')
construction = json.loads((source / 'construction_summary.json').read_text())
counts = sparse.load_npz(source / 'selected_training_counts_csr.npz')
assert counts.shape == (len(cells), len(genes)) and cells.cell_id.is_unique
assert cells.split.eq('training').all()
assert not cells.Patient.isin(roles.loc[roles.split.ne('training'), 'Patient']).any()
full_library = counts.sum(axis=1).A1
assert (full_library > 0).all() and np.array_equal(full_library, cells.extracted_total_counts.to_numpy())
print('v2 counts checkpoint (cells, genes):', counts.shape, '; training donors:', cells.Patient.nunique(), flush=True)

# Step2. 每个原参考单元固定抽5个细胞，只为控制网页体积，不再按亚型、marker或预测结果筛选
cells_per_profile = 5  # 本次诊断约定；保留原212个参考单元，而不是声称每种细胞只有5个
seed = 20261008
chosen = cells.groupby('profile_id', sort=True, observed=True).sample(n=cells_per_profile, random_state=seed).sort_index()
assert chosen.profile_id.nunique() == len(profiles) and chosen.Patient.nunique() == cells.Patient.nunique()
assert set(chosen.coarse_identity) == set(construction['types'])
chosen.to_csv(out / 'source_subset_cell_manifest.tsv', sep='\t', index=False)
composition = chosen.groupby(['coarse_identity', 'region'], observed=True).agg(cells=('cell_id', 'size'), profiles=('profile_id', 'nunique'), donors=('Patient', 'nunique')).reset_index()
composition.to_csv(out / 'source_subset_type_region_composition.tsv', sep='\t', index=False)
print(composition.to_string(index=False), flush=True)

# Step3. 保留原v2的完整上传基因背景；分母为每细胞全部33075源基因counts，不对删选后的基因再归一化
background = pd.read_csv(source / 'GJ20v2_DonorSite_CPM.txt', sep='\t', usecols=[0]).iloc[:, 0]
assert genes.gene_symbol.is_unique and genes.gene_symbol.notna().all()
gene_positions = pd.Index(genes.gene_symbol).get_indexer(background)
assert (gene_positions >= 0).all() and background.is_unique
selected_counts = counts[chosen.index.to_numpy(), :][:, gene_positions]
cell_cpm = selected_counts.multiply((1e6 / full_library[chosen.index.to_numpy()])[:, None]).toarray()
assert np.isfinite(cell_cpm).all() and (cell_cpm >= 0).all() and cell_cpm.max() > 50
assert cell_cpm.size < 50_000_000
expression = pd.DataFrame(cell_cpm.T, index=background.to_numpy(), columns=chosen.coarse_identity.to_numpy())
# scRNA输入首行重复细胞类型标签是网页要求；唯一cell ID及列顺序保存在manifest，而不是把供者汇总列假装成单细胞。
expression.to_csv(out / 'GJ20v2_SourceSubset1060_scRNA_CPM.txt', sep='\t', index_label='GeneSymbol', float_format='%.8g')

# Step4. 保存限制与检查点哈希，不修改v2 signature、不使用留出/recovery细胞、不提交网页任务
input_paths = ['selected_training_counts_csr.npz', 'selected_cell_manifest.tsv', 'source_gene_manifest.tsv', 'profile_manifest.tsv', 'frozen_donor_roles.tsv', 'GJ20v2_DonorSite_CPM.txt']
input_hashes = {}
for name in input_paths:  # 多个输入文件的校验所需
    with (source / name).open('rb') as handle:
        input_hashes[name] = hashlib.file_digest(handle, 'sha256').hexdigest()
summary = {'environment_python': sys.executable, 'source_v2_dir': str(source), 'counts_origin': 'original /X integer counts saved in selected_training_counts_csr.npz', 'seed': seed, 'cells_per_original_profile': cells_per_profile, 'cells': len(chosen), 'genes': len(background), 'data_points': int(expression.size), 'original_training_cells': len(cells), 'original_profiles_retained': len(profiles), 'original_training_donors_retained': chosen.Patient.nunique(), 'cell_types': construction['types'], 'holdout_cells_used': 0, 'full_source_library_genes': len(genes), 'log_transformed': False, 'signature_modified': False, 'exact_signature_construction_input': False, 'limitation': 'Job13 v2 used 212 donor-site pooled CPM profiles; this is a 5-cell-per-profile subset of their source cells, not the literal original signature-building input. A web S-mode diagnostic using this file is exploratory, not a confirmed standard author-protocol pairing.', 'input_sha256': input_hashes}
(out / 'source_subset_summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
checksums = []
for path in sorted(out.iterdir()):
    if path.is_file():
        with path.open('rb') as handle:
            digest = hashlib.file_digest(handle, 'sha256').hexdigest()
        checksums.append({'file': path.name, 'bytes': path.stat().st_size, 'sha256': digest})
pd.DataFrame(checksums).to_csv(out / 'transfer_checksums.tsv', sep='\t', index=False)
print('SOURCE_SUBSET_READY:', out, '; gene x cell:', expression.shape, '; data points:', expression.size, flush=True)
