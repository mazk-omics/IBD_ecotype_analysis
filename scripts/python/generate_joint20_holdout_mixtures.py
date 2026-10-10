# 用未参与 signature 构建的 GSE282122 供者生成模拟 bulk、真实比例及细胞清单。
# 在服务器 ibd_refaudit 环境运行；读取原始 X 计数和既有供者划分，不修改来源数据。
# 输出到一个新的目录：28 个原始计数混合样本，以及其中 8 个样本的等 RNA 贡献配对对照。
import argparse
import gzip
import hashlib
import json
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

# 复用已检查的 H5AD 注释解码和稀疏行读取；抽样、混合及保存步骤在本脚本直接展开。
from inspect_ibdverse_h5ad import read_column
from export_ibdverse_cibersortx_pilot import selected_csr

# Step1. 明确输入输出，固定这轮验证的数据和参数
parser = argparse.ArgumentParser()
parser.add_argument('--outdir', required=True)
args = parser.parse_args()
out = Path(args.outdir)
out.mkdir(parents=True, exist_ok=False)  # 不覆盖已经生成的测试数据。
started = time.monotonic()
project = Path('/home/mazekai/IBD_EcoTyper')
source = Path('/home/mazekai/enteric_glia/data_scRNA/GSE282122/TAURUS_raw_counts_annotated_final.h5ad')
candidate_path = project / '03_reference/GSE282122/audit_v2/output/04_identity_coverage_v1/mapped_candidate_cells_v1.parquet'
pilot = project / '03_reference/web_reference_three_pilots/output/run_20261006_183211/gse_joint'
types = json.loads((pilot / 'summary.json').read_text())['web_reference_types']
candidate = pd.read_parquet(candidate_path)
donors = pd.read_csv(pilot / 'donor_split.tsv', sep='\t')
training_cells = pd.read_csv(pilot / 'selected_cell_manifest.tsv', sep='\t')
assert donors.Patient.is_unique and candidate.cell_id.is_unique

# 本项目约定：两位最终检验供者现在只保留身份，不读取其表达；按注释预先固定，不看性能挑选。
reserved = ['CD12', 'UC17']
assert set(reserved) <= set(donors.loc[donors.split.eq('holdout'), 'Patient'])
donors['validation_role'] = donors['split'].replace({'training': 'signature_training', 'holdout': 'pilot_validation'})
donors.loc[donors.Patient.isin(reserved), 'validation_role'] = 'reserved_final_test'
donors.to_csv(out / 'donor_roles.tsv', sep='\t', index=False)
development = donors.loc[donors.validation_role.eq('pilot_validation'), 'Patient'].tolist()
print('验证供者:', development, '未读取表达的最终检验供者:', reserved, flush=True)

# Step2. 沿用训练参考的可用细胞规则，核对部位和类型容量
# 不新增线粒体或疾病阈值；既有注释中的 Glial/Tuft 与主候选一致纳入。
valid = candidate.Patient.isin(development) & candidate.reference_eligible.eq(True)
valid &= candidate.identity_status.eq('CONFIDENT') & candidate.coarse_identity.isin(types)
valid &= candidate.primary_reference_action.eq('KEEP') | candidate.coarse_identity.isin(['Glial', 'Tuft'])
valid &= ~candidate.predicted_doublets & candidate.total_counts.gt(0) & candidate.n_genes_by_counts.gt(0)
eligible = candidate.loc[valid].copy()
eligible['pilot_type'] = eligible['coarse_identity']
eligible['region'] = np.where(eligible.Site.eq('Terminal_Ileum'), 'Ileum', 'ColonRectum')
assert not eligible.Patient.isin(donors.loc[donors.split.eq('training'), 'Patient']).any()
assert not eligible.cell_id.isin(training_cells.cell_id).any()
assert not eligible.Patient.isin(reserved).any()
capacity = eligible.groupby(['region', 'pilot_type'], observed=True).agg(n_cells=('cell_id', 'size'), n_donors=('Patient', 'nunique')).reset_index()
capacity.to_csv(out / 'eligible_type_region_capacity.tsv', sep='\t', index=False)
print('合格留出细胞:', len(eligible), '类型/部位容量:', capacity.shape, flush=True)

# 每种类型、部位最多读取 200 个细胞。供者轮流取样，避免大样本供者垄断池子。
# 该循环是分组抽样所需结构；未达到 200 的稀有类型保留全部可用细胞，不补训练细胞。
rng = np.random.default_rng(20261006)
pool_indices = []
for (region, cell_type), group in eligible.groupby(['region', 'pilot_type'], observed=True, sort=True):
    donor_groups = [rng.permutation(x.index.to_numpy()).tolist() for _, x in group.groupby('Patient', observed=True, sort=True)]
    target = min(200, len(group))
    chosen = []
    while len(chosen) < target:
        for donor_group in donor_groups:
            if donor_group and len(chosen) < target:
                chosen.append(donor_group.pop())
    pool_indices.extend(chosen)
pool = eligible.loc[pool_indices].sort_values('cell_id').reset_index(drop=True)
pool.to_csv(out / 'selected_pool_manifest.tsv', sep='\t', index=False)
print('待读取表达的细胞池:', len(pool), '供者:', pool.Patient.nunique(), flush=True)

# Step3. 按真实细胞 ID 对齐原始 H5AD，仅提取所选稀疏行
with h5py.File(source, 'r') as h:
    ids = read_column(h['obs'][h['obs'].attrs.get('_index', '_index')]).astype(str)
    assert ids.is_unique
    source_rows = pd.Index(ids).get_indexer(pool.obs_name.astype(str))
    assert (source_rows >= 0).all()
    for column in ['Patient', 'Site', 'Disease', 'final_analysis']:
        observed = read_column(h['obs'][column]).iloc[source_rows].astype(str).to_numpy()
        assert np.array_equal(observed, pool[column].astype(str).to_numpy()), column
    order = np.argsort(source_rows)
    source_rows = source_rows[order]
    pool = pool.iloc[order].reset_index(drop=True)
    pool['source_cell_row'] = source_rows
    genes = pd.DataFrame({'source_gene_id': read_column(h['var']['gene_id']).astype(str), 'gene_symbol': read_column(h['var']['gene_symbol']).astype(str)})
    indices, values, indptr, totals, detected = selected_csr(h['X'], source_rows, len(genes))

assert np.allclose(totals, pool.total_counts, rtol=0, atol=.01)
assert np.array_equal(detected, pool.n_genes_by_counts.to_numpy())
pool['extracted_total_counts'] = totals
pool['extracted_n_genes'] = detected
counts_original = csr_matrix((values, indices, indptr), shape=(len(pool), len(genes)))

# 与当前参考一致：剔除空白及重复 symbol 的所有行，不任意选择一个重复基因位点。
# 不用 signature 基因列表筛选测试表达；保留完整可用表达背景。
genes['duplicate_symbol'] = genes.gene_symbol.duplicated(keep=False)
keep = genes.gene_symbol.ne('') & ~genes.duplicate_symbol & ~genes.gene_symbol.str.contains(r'[\t\r\n]')
counts = counts_original[:, keep.to_numpy()].astype(np.float64)
export_genes = genes.loc[keep, 'gene_symbol'].to_numpy()
library = np.asarray(counts.sum(axis=1)).ravel()
assert (library > 0).all()
pool['export_library_counts'] = library
pool.to_csv(out / 'selected_pool_manifest.tsv', sep='\t', index=False)
genes['valid_unique_symbol'] = keep
genes.to_csv(out / 'gene_export_audit.tsv', sep='\t', index=False)
print('稀疏计数:', counts.shape, '与原始细胞计数核对全部通过', flush=True)

# 提取检查点用于复核或复现；它是自定义 CSR H5 文件，不是 AnnData 对象。
with h5py.File(out / 'holdout_pool_counts.h5', 'w') as h:
    h.attrs['format'] = 'Custom CSR, cells x unique genes; not AnnData'
    h.attrs['shape'] = counts.shape
    h.create_dataset('data', data=counts.data, compression='gzip')
    h.create_dataset('indices', data=counts.indices, compression='gzip')
    h.create_dataset('indptr', data=counts.indptr, compression='gzip')
    strings = h5py.string_dtype('utf-8')
    h.create_dataset('gene_symbol', data=np.asarray(export_genes, dtype=object), dtype=strings)
    h.create_dataset('cell_id', data=np.asarray(pool.cell_id, dtype=object), dtype=strings)

# Step4. 固定 28 个主测试场景；每个样本 200 个不重复细胞
# 这些是受控场景，比例不代表 MSCCR 真值。跨模拟样本可复用细胞，不能按 28 位独立供者统计。
immune = ['B', 'Plasma', 'CD4_T', 'CD8_T', 'NK', 'ILC', 'Monocyte_Macrophage', 'DC', 'Mast']
stromal = ['Fibroblast', 'Pericyte', 'Endothelial', 'Glial']
epithelial = [x for x in types if x not in immune + stromal]
specifications = []
for region in ['Ileum', 'ColonRectum']:
    region_types = [x for x in types if x != 'Colonic_absorptive'] if region == 'Ileum' else [x for x in types if x not in ['Ileal_absorptive', 'Paneth']]
    for scenario in ['balanced', 'epithelial_high', 'immune_high']:
        for replicate in range(1, 5):
            specifications.append({'sample': f'{region}_{scenario}_{replicate}', 'region': region, 'scenario': scenario, 'types': region_types, 'absent': ''})
for absent in ['DC', 'Glial', 'Tuft', 'Paneth']:
    specifications.append({'sample': f'Joint_absent_{absent}', 'region': 'Joint', 'scenario': 'absent_and_rare', 'types': types, 'absent': absent})

mixtures = []
controlled = []
cell_truth = []
rna_truth = []
original_rna_truth = []
sample_records = []
membership = []
paired_samples = []
for specification in specifications:
    sample = specification['sample']
    region = specification['region']
    active = [x for x in specification['types'] if x != specification['absent']]
    scenario = specification['scenario']
    base = np.ones(len(active))
    if scenario in ['epithelial_high', 'immune_high', 'absent_and_rare']:
        group_mass = {'epithelial_high': [.70, .20, .10], 'immune_high': [.20, .70, .10], 'absent_and_rare': [.50, .35, .15]}[scenario]
        for group, mass in zip([epithelial, immune, stromal], group_mass):
            positions = [i for i, x in enumerate(active) if x in group]
            base[positions] = mass / len(positions)
    weights = rng.dirichlet(base / base.sum() * 100)
    n_cells = 200  # 本轮项目约定，改变后会改变稀有类型最小分辨率。
    allocated = np.ones(len(active), dtype=int)
    if scenario == 'absent_and_rare':
        fixed = {'Paneth': 1, 'Glial': 2, 'Tuft': 2, 'DC': 6}
        for cell_type, number in fixed.items():
            if cell_type in active:
                allocated[active.index(cell_type)] = number
        unfixed = np.array([x not in fixed for x in active])
        scaled = weights[unfixed] / weights[unfixed].sum() * (n_cells - allocated.sum())
        add = np.floor(scaled).astype(int)
        add[np.argsort(-(scaled - add))[:n_cells - allocated.sum() - add.sum()]] += 1
        allocated[unfixed] += add
    else:
        scaled = weights * (n_cells - allocated.sum())
        add = np.floor(scaled).astype(int)
        add[np.argsort(-(scaled - add))[:n_cells - allocated.sum() - add.sum()]] += 1
        allocated += add
    selected = []
    for cell_type, number in zip(active, allocated):
        available = pool.index[pool.pilot_type.eq(cell_type)]
        if region != 'Joint':
            available = pool.index[pool.pilot_type.eq(cell_type) & pool.region.eq(region)]
        assert len(available) >= number, f'{sample}/{cell_type}: need {number}, available {len(available)}'
        selected.extend(rng.choice(available.to_numpy(), int(number), replace=False).tolist())
    selected = np.sort(selected)
    assert len(selected) == n_cells and len(np.unique(selected)) == n_cells
    mixture = np.asarray(counts[selected].sum(axis=0)).ravel()
    mixture = mixture / mixture.sum() * 1e6
    mixtures.append(mixture)
    type_counts = pool.loc[selected, 'pilot_type'].value_counts().reindex(types, fill_value=0)
    type_rna = pool.loc[selected].groupby('pilot_type', observed=True).export_library_counts.sum().reindex(types, fill_value=0)
    type_original_rna = pool.loc[selected].groupby('pilot_type', observed=True).extracted_total_counts.sum().reindex(types, fill_value=0)
    cell_truth.append(type_counts.to_numpy() / n_cells)
    rna_truth.append(type_rna.to_numpy() / type_rna.sum())
    original_rna_truth.append(type_original_rna.to_numpy() / type_original_rna.sum())
    selected_manifest = pool.loc[selected, ['cell_id', 'Patient', 'Disease', 'Site', 'pilot_type', 'export_library_counts']].copy()
    selected_manifest.insert(0, 'sample', sample)
    membership.append(selected_manifest)
    sample_records.append({'sample': sample, 'region': region, 'scenario': scenario, 'absent_type': specification['absent'], 'n_cells': n_cells, 'n_donors': pool.loc[selected, 'Patient'].nunique(), 'unit': 'sum_raw_counts_then_CPM'})
    # 8 个配对对照使用完全相同的细胞，先把每个细胞归一化再平均，使 cell/RNA 比例一致。
    if region != 'Joint' and (sample.endswith('balanced_1') or sample.endswith('balanced_2') or sample.endswith('epithelial_high_1') or sample.endswith('immune_high_1')):
        weighted = counts[selected].multiply((1e6 / library[selected])[:, None])
        controlled.append(np.asarray(weighted.mean(axis=0)).ravel())
        paired_samples.append({'sample': sample + '_EqualRNA', 'paired_main_sample': sample})

main_samples = [x['sample'] for x in specifications]
main = np.column_stack(mixtures)
equal_rna = np.column_stack(controlled)
expressed = (main.sum(axis=1) + equal_rna.sum(axis=1)) > 0
main = main[expressed]
equal_rna = equal_rna[expressed]
export_genes = export_genes[expressed]
genes['mixture_expressed'] = genes.gene_symbol.isin(export_genes)
genes.to_csv(out / 'gene_export_audit.tsv', sep='\t', index=False)
assert np.isfinite(main).all() and np.isfinite(equal_rna).all() and (main >= 0).all() and (equal_rna >= 0).all()
assert np.allclose(main.sum(axis=0), 1e6) and np.allclose(equal_rna.sum(axis=0), 1e6)
assert main.size < 50_000_000 and equal_rna.size < 50_000_000 and main.max() >= 50 and equal_rna.max() >= 50

# Step5. 导出网页输入和真实比例，顺序与 signature 列顺序分开保留，评估时按名称匹配
pd.DataFrame(main, index=export_genes, columns=main_samples).to_csv(out / 'GJ20_Holdout28_Mixture_CPM.txt', sep='\t', index_label='GeneSymbol', float_format='%.8g')
equal_samples = [x['sample'] for x in paired_samples]
pd.DataFrame(equal_rna, index=export_genes, columns=equal_samples).to_csv(out / 'GJ20_Holdout8_EqualRNA_Mixture_CPM.txt', sep='\t', index_label='GeneSymbol', float_format='%.8g')
truth_cells = pd.DataFrame(cell_truth, index=main_samples, columns=types)
truth_rna = pd.DataFrame(rna_truth, index=main_samples, columns=types)
truth_original = pd.DataFrame(original_rna_truth, index=main_samples, columns=types)
truth_cells.to_csv(out / 'ground_truth_cell_fractions.tsv', sep='\t', index_label='sample')
truth_rna.to_csv(out / 'ground_truth_export_UMI_fractions.tsv', sep='\t', index_label='sample')
truth_original.to_csv(out / 'ground_truth_original_UMI_fractions.tsv', sep='\t', index_label='sample')
equal_truth = truth_cells.loc[[x['paired_main_sample'] for x in paired_samples]].copy()
equal_truth.index = equal_samples
equal_truth.to_csv(out / 'ground_truth_equalRNA_fractions.tsv', sep='\t', index_label='sample')
pd.DataFrame(sample_records).to_csv(out / 'mixture_design.tsv', sep='\t', index=False)
pd.DataFrame(paired_samples).to_csv(out / 'paired_equalRNA_samples.tsv', sep='\t', index=False)
pd.concat(membership, ignore_index=True).to_csv(out / 'mixture_cell_membership.tsv', sep='\t', index=False)
discrepancy = (truth_cells - truth_rna).abs()
pd.DataFrame({'cell_type': types, 'mean_abs_cell_vs_export_UMI_fraction': discrepancy.mean().to_numpy(), 'max_abs_cell_vs_export_UMI_fraction': discrepancy.max().to_numpy()}).to_csv(out / 'cell_vs_UMI_truth_difference.tsv', sep='\t', index=False)

# 精确重现抽样和已输出表达所需的中间对象保存在服务器；本地回传网页输入和审计表。
files = {}
for file in sorted(out.glob('*.txt')):
    content = file.read_bytes()
    files[file.name] = {'bytes': len(content), 'sha256': hashlib.sha256(content).hexdigest()}
summary = {'status': 'completed', 'seed': 20261006, 'source': str(source), 'source_layer': 'X raw counts', 'source_opened_read_only': True,
    'python': sys.executable, 'development_donors': development, 'reserved_final_test_donors': reserved, 'source_pool_cells': len(pool),
    'unique_used_cells': int(pd.concat(membership).cell_id.nunique()), 'types': types, 'main_genes': len(export_genes), 'main_samples': len(main_samples),
    'equalRNA_samples': len(equal_samples), 'cells_per_main_sample': 200, 'within_sample_sampling': 'without replacement',
    'across_sample_cells_reused': True, 'simulation_samples_are_independent_donors': False, 'training_cell_overlap': 0,
    'training_donor_overlap': 0, 'reserved_donor_expression_read': False, 'gene_filter': 'unique valid source symbols, positive across mixtures; no signature preselection',
    'main_unit': 'sum raw counts then CPM on exported unique genes', 'equalRNA_unit': 'mean per-cell CPM on exported unique genes',
    'colorectal_panel_Paneth': 'deliberately excluded in regional baseline; source does have annotated colorectal Paneth cells',
    'main_data_points': int(main.size), 'equalRNA_data_points': int(equal_rna.size), 'files': files,
    'seconds': round(time.monotonic() - started, 2), 'CIBERSORTx_deconvolution_run': False}
(out / 'summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)
