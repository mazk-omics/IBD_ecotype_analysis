# 从作者已质控的 IBDverse 回肠图谱抽取独立供者，准备官方 EcoTyper 单细胞 recovery 输入。
# 在服务器 ibd_refaudit 环境运行；读取原始 counts 和已审计的元数据，原文件只读。
# 输出抽样记录、线性 CPM（基因×细胞）、原作者注释和独立谱系 marker 汇总；不重新聚类或修改发现模型。
import argparse
import json
import sys
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

# Step1. 读入已审计的供者映射和原作者细胞注释
parser = argparse.ArgumentParser()
parser.add_argument('--outdir', required=True)
args = parser.parse_args()
out = Path(args.outdir)
out.mkdir(parents=True, exist_ok=False)
project = Path('/home/mazekai/IBD_EcoTyper')
source = Path('/home/mazekai/enteric_glia/data_scRNA/scIBD_external/IBDverse_E-MTAB-16998/E-MTAB-16999/Full_cohort.h5ad')
sys.path.insert(0, str(project / 'scripts/python/ibdverse_pilot_v1'))
from inspect_ibdverse_h5ad import read_column, frame_columns
coverage = pd.read_csv(project / 'logs/ibdverse_reference_inputs_20261006_161527/label_category_sample_donor_coverage.tsv', sep='\t')
cols = ['sanger_sample_id', 'individual_id', 'disease_status', 'ti_ses_inflamed', 'biopsy_type', 'reference_metadata_complete']
metadata = coverage[cols].drop_duplicates()
assert not metadata.sanger_sample_id.duplicated().any(), '样本映射必须一对一'
complete = metadata.reference_metadata_complete.astype(str).str.lower().eq('true')
known = metadata[complete & metadata.biopsy_type.eq('terminal ileum')].copy()
known['Stratum'] = known.disease_status + '_' + known.ti_ses_inflamed
label_map = pd.read_csv(Path(__file__).parent / 'ibdverse_recovery_label_map.tsv', sep='\t')
types = ['B', 'CD4_T', 'CD8_T', 'Endothelial', 'Enteroendocrine', 'Fibroblast', 'Goblet', 'ILC', 'Mast', 'Monocyte_Macrophage', 'Pericyte']

# Step2. 先选择供者，避免根据表达或恢复结果挑选表现好的细胞
rng = np.random.default_rng(20261009)  # 项目试验种子；每组 10 人，不代表原队列患病率
chosen = []
for stratum in ['Healthy_uninflamed', 'CD_inflamed', 'CD_uninflamed']:
    donors = np.sort(known.loc[known.Stratum.eq(stratum), 'individual_id'].dropna().unique())
    selected = rng.choice(donors, size=min(10, len(donors)), replace=False)
    chosen.append(pd.DataFrame({'individual_id': selected, 'Stratum': stratum}))
donors = pd.concat(chosen, ignore_index=True)
assert not donors.individual_id.duplicated().any(), '供者不能跨分层重复入选'
samples = known.merge(donors, on=['individual_id', 'Stratum'], validate='many_to_one')
samples = samples.sort_values('sanger_sample_id').drop_duplicates('individual_id')  # 多样本时固定取排序首个回肠样本
samples.to_csv(out / 'selected_donors_samples.tsv', sep='\t', index=False)
print('Selected donors:', samples.groupby('Stratum').size().to_dict(), flush=True)

with h5py.File(source, 'r') as h5:
    obs = h5['obs']
    index_key, obs_columns = frame_columns(obs)
    cells = pd.DataFrame({'ID': read_column(obs[index_key]).astype(str), 'sanger_sample_id': read_column(obs['sanger_sample_id']).astype(str), 'AuthorLabel': read_column(obs['label']).astype(str), 'AuthorCategory': read_column(obs['category']).astype(str)})
    cells['SourceRow'] = np.arange(len(cells))
    cells = cells.merge(samples, on='sanger_sample_id', how='inner', validate='many_to_one')
    cells = cells.merge(label_map[['label', 'pilot_type']], left_on='AuthorLabel', right_on='label', how='left', validate='many_to_one')
    cells['CellType'] = cells.pilot_type
    print('Available cells in selected biopsies:', len(cells), flush=True)

    # Step3. 每位供者×作者大类随机保留至多 20 个细胞，不强求亚型等额
    picked = []
    for (donor, cell_type), group in cells.dropna(subset=['CellType']).groupby(['individual_id', 'CellType'], sort=True):
        take = rng.choice(group.index.to_numpy(), size=min(20, len(group)), replace=False)
        picked.append(cells.loc[take])
    selected_cells = pd.concat(picked).sort_values('SourceRow').reset_index(drop=True)
    selected_cells['FormalRecovery'] = selected_cells.CellType.isin(types)
    selected_cells = selected_cells.rename(columns={'individual_id': 'Donor', 'sanger_sample_id': 'Sample', 'disease_status': 'Disease', 'ti_ses_inflamed': 'Inflammation'})
    selected_cells.to_csv(out / 'all_selected_cells.tsv', sep='\t', index=False)
    formal = selected_cells[selected_cells.FormalRecovery].copy()
    formal.to_csv(out / 'IBDverse_Pilot30_20261009_annotation.tsv', sep='\t', index=False)
    selected_cells.groupby(['CellType', 'Stratum'], observed=True).size().unstack(fill_value=0).to_csv(out / 'sampling_counts.tsv', sep='\t')
    print('Formal recovery cells:', len(formal), formal.CellType.value_counts().to_dict(), flush=True)
    assert not formal.ID.duplicated().any()
    assert set(types) <= set(formal.CellType), '未覆盖的谱系应先审阅，而不是虚构输入'

    # Step4. 只提取抽中的 CSR 行；基因符号重复时相加，CPM 分母使用完整基因库
    layer = h5['layers/counts']
    n_genes = int(layer.attrs['shape'][1])  # gene_symbols 是分类编码，不能用其 HDF5 节点个数作为基因数
    ptr = layer['indptr'][:]
    pieces = []
    for row in selected_cells.SourceRow:
        start, stop = ptr[row:row + 2]
        indices = layer['indices'][start:stop].astype(np.int32)
        values = layer['data'][start:stop]
        assert len(indices) == len(values) and (indices >= 0).all() and (indices < n_genes).all()
        pieces.append(sparse.csr_matrix((values, indices, [0, len(values)]), shape=(1, n_genes)))
    counts = sparse.vstack(pieces, format='csr')
    symbols = read_column(h5['var']['gene_symbols']).astype(str).to_numpy()
    unique_symbols, inverse = np.unique(symbols, return_inverse=True)
    collapse = sparse.csr_matrix((np.ones(len(symbols)), (np.arange(len(symbols)), inverse)), shape=(len(symbols), len(unique_symbols)))
    totals = np.asarray(counts.sum(axis=1)).ravel()
    assert (totals > 0).all() and np.isfinite(counts.data).all() and (counts.data >= 0).all()
    count_qc = selected_cells[['ID', 'CellType', 'Donor', 'Stratum', 'FormalRecovery']].copy()
    count_qc['FullLibraryCounts'] = totals
    count_qc['DetectedFeatures'] = counts.getnnz(axis=1)
    count_qc.to_csv(out / 'selected_count_qc.tsv', sep='\t', index=False)
    cpm = sparse.diags(1e6 / totals) @ counts @ collapse
    print('Counts:', counts.shape, 'unique symbols:', len(unique_symbols), 'library min/max:', float(totals.min()), float(totals.max()), flush=True)

# Step5. 保留模型需要的全部基因；独立谱系 marker 仅供检查，不参与拟合或选细胞
model = project / 'software/ecotyper/EcoTyper/MSCCR_CORE12_CTSPEC_REGION_IBD/Cell_type_specific_genes/Cell_States/discovery'
gene_sets = []
overlap = []
for cell_type in types:
    gene_file = list((model / cell_type).glob('*/initial_gene_info.txt'))
    assert len(gene_file) == 1
    info = pd.read_csv(gene_file[0], sep='\t', index_col=0)
    genes = info['Gene'].astype(str).unique()
    gene_sets.extend(genes)
    overlap.append({'CellType': cell_type, 'ModelGenes': len(genes), 'PresentGenes': int(np.isin(genes, unique_symbols).sum()), 'MissingGenes': '|'.join(genes[~np.isin(genes, unique_symbols)])})
pd.DataFrame(overlap).to_csv(out / 'model_gene_overlap.tsv', sep='\t', index=False)
markers = ['TPSAB1', 'TPSB2', 'KIT', 'CPA3', 'CSF3R', 'CXCR2', 'S100A8', 'S100A9', 'CSF1R', 'MS4A7', 'MS4A1', 'CD22', 'PAX5', 'CD3D', 'CD3E', 'TCF7', 'NKG7', 'EPCAM', 'KRT19', 'CHGA', 'CHGB', 'TPH1', 'MUC2', 'TFF3', 'MKI67', 'TOP2A', 'OLFM4', 'NOS2', 'DUOX2', 'LCN2', 'CXCL1', 'CXCL2', 'COL1A1', 'DCN', 'PECAM1', 'VWF', 'RGS5', 'PDGFRB', 'SOX10', 'S100B', 'MPZ', 'CUBN', 'ALDOB', 'DEFA5', 'DEFA6']
marker_index = pd.Index(unique_symbols).get_indexer(markers)
valid_markers = np.array(markers)[marker_index >= 0]
marker_values = cpm[:, marker_index[marker_index >= 0]].toarray()
marker_frame = pd.DataFrame(marker_values, columns=valid_markers, index=selected_cells.ID)
marker_frame.to_csv(out / 'diagnostic_marker_CPM.tsv', sep='\t')
marker_frame.join(selected_cells.set_index('ID')[['CellType', 'AuthorLabel', 'Donor', 'Stratum']]).to_csv(out / 'diagnostic_marker_with_annotation.tsv', sep='\t')

genes = sorted(set(gene_sets) & set(unique_symbols))
gene_index = pd.Index(unique_symbols).get_indexer(genes)
formal_rows = np.flatnonzero(selected_cells.FormalRecovery.to_numpy())
expression = pd.DataFrame(cpm[formal_rows][:, gene_index].toarray().T, index=genes, columns=formal.ID)
expression.index.name = 'GeneSymbol'
expression.to_csv(out / 'IBDverse_Pilot30_20261009_CPM.txt', sep='\t', float_format='%.7g')
assert expression.columns.to_list() == formal.ID.to_list()
summary = {'source': str(source), 'discovery': 'MSCCR_CORE12_CTSPEC_REGION_IBD', 'seed': 20261009, 'donors': len(samples), 'formal_cells': len(formal), 'diagnostic_total_cells': len(selected_cells), 'formal_types': types, 'genes': len(genes), 'normalization': 'linear CPM; full-library denominator; duplicate symbols summed; official recovery applies log2(CPM+1)', 'cap_per_donor_type': 20, 'glial': 'not available in author labels; omitted', 'model_gene_overlap': overlap}
(out / 'preparation_summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
print('READY:', out, expression.shape, flush=True)
