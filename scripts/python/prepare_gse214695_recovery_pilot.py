# 下载 GSE214695 作者计数矩阵，并为冻结的 Region3 CD+UC discovery 准备独立结肠 recovery。
# 在服务器 ibd_refaudit 环境运行；原始文件单独保存，读取作者细胞注释，不重新聚类或修改 EcoTyper。
# 每个样本×大类最多 100 个细胞；完整基因库计算 CPM 后才提取模型基因，输出基因×细胞矩阵。
from pathlib import Path
import gzip
import hashlib
import json
import shutil
import subprocess
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread

# Step1. 读取下载清单、样本信息和显式谱系映射
root = Path('/home/mazekai/IBD_EcoTyper/06_validation/GSE214695_scRNA_recovery_20261009')
raw = root / 'source'
out = root / 'input'
raw.mkdir(exist_ok=True)
out.mkdir(exist_ok=True)
manifest = pd.read_csv(root / 'download_urls.tsv', sep='\t')
samples = pd.read_csv(root / 'sample_metadata.tsv', sep='\t')
label_map = pd.read_csv(root / 'author_label_map.tsv', sep='\t', keep_default_na=False)
assert not label_map.AuthorLabel.duplicated().any() and not samples.Sample.duplicated().any()
types = ['B', 'CD4_T', 'CD8_T', 'Endothelial', 'Enteroendocrine', 'Fibroblast', 'Glial', 'Goblet', 'ILC', 'Mast', 'Monocyte_Macrophage', 'Pericyte']

# Step2. 顺序下载，不建立 SSH 新连接；失败文件不冒充已完成文件
downloads = []
for row in manifest.itertuples(index=False):
    target = raw / row.Filename
    if not target.exists():
        partial = target.with_suffix(target.suffix + '.part')
        print('Downloading:', row.Filename, flush=True)
        # NCBI 下载发生 TLS 中断时仅重试文件传输；这不会新建 SSH 登录。
        subprocess.run(['curl', '--fail', '--location', '--retry', '2', '--retry-delay', '5', '--connect-timeout', '30', '--max-time', '900', '--silent', '--show-error', '--output', str(partial), row.URL], check=True)
        partial.replace(target)
    with target.open('rb') as handle:
        digest = hashlib.file_digest(handle, 'sha256').hexdigest()
    downloads.append({'Filename': row.Filename, 'Bytes': target.stat().st_size, 'SHA256': digest})
    print('Downloaded:', row.Filename, target.stat().st_size, flush=True)
pd.DataFrame(downloads).to_csv(root / 'source_checksums.tsv', sep='\t', index=False)

# Step3. 注释先抽样；不按恢复得分、marker 表达或细亚型挑选细胞
cells = pd.read_csv(raw / 'GSE214695_cell_annotation.csv.gz').rename(columns={'sample': 'Sample', 'cell_id': 'Barcode', 'annotation': 'AuthorLabel', 'Unnamed: 0': 'AuthorID'})
# 作者注释有2个同样本同barcode却标签冲突的细胞；隔离全部冲突记录，不任意保留一行，原件不改。
conflicting_cells = cells[cells.duplicated(['Sample', 'Barcode'], keep=False)].copy()
conflicting_cells.to_csv(out / 'quarantined_annotation_conflicts.tsv', sep='\t', index=False)
cells = cells[~cells.duplicated(['Sample', 'Barcode'], keep=False)].copy()
print('Annotation quarantine:', len(conflicting_cells), 'rows;', len(conflicting_cells[['Sample', 'Barcode']].drop_duplicates()), 'ambiguous cells excluded', flush=True)
assert not cells[['Sample', 'Barcode']].duplicated().any()
assert set(cells.AuthorLabel) == set(label_map.AuthorLabel)
cells = cells.merge(label_map, on='AuthorLabel', validate='many_to_one')
cells = cells.merge(samples, on='Sample', validate='many_to_one')
cells['ID'] = 'GSE214695_' + cells.Sample + '_' + cells.Barcode
cells['Donor'] = 'GSE214695_' + cells.Sample
cells['FormalRecovery'] = cells.CellType.isin(types)
cells.groupby(['CellType', 'Sample']).size().to_csv(out / 'available_cell_counts.tsv', sep='\t')
rng = np.random.default_rng(20261009)  # 项目试验约定；每个样本×大类上限 100，不强制亚型等额
picked = []
for (sample, cell_type), group in cells.groupby(['Sample', 'CellType'], sort=True):
    take = rng.choice(group.index.to_numpy(), size=min(100, len(group)), replace=False)
    picked.append(cells.loc[take])
selected = pd.concat(picked).sort_values(['Sample', 'Barcode']).reset_index(drop=True)
selected.to_csv(out / 'all_selected_cells.tsv', sep='\t', index=False)
formal = selected[selected.FormalRecovery].copy()
formal.to_csv(out / 'GSE214695_Pilot18_20261009_annotation.tsv', sep='\t', index=False)
selected.groupby(['CellType', 'Sample']).size().to_csv(out / 'sampling_counts.tsv', sep='\t')
assert formal.ID.is_unique and set(types) <= set(formal.CellType)
print('Formal cells:', len(formal), formal.CellType.value_counts().to_dict(), flush=True)

# Step4. 读取冻结模型的基因并另存诊断 marker；诊断 marker 不参与模型拟合
model = Path('/home/mazekai/IBD_EcoTyper/software/ecotyper/EcoTyper/MSCCR_CORE12_CTSPEC_REGION_IBD/Cell_type_specific_genes/Cell_States/discovery')
model_genes = {}
for cell_type in types:
    gene_files = list((model / cell_type).glob('*/initial_gene_info.txt'))
    assert len(gene_files) == 1
    info = pd.read_csv(gene_files[0], sep='\t', index_col=0)
    model_genes[cell_type] = info.Gene.astype(str).unique()
genes = sorted(set(np.concatenate(list(model_genes.values()))))
source_symbols = []
for filename in manifest.loc[manifest.Filename.str.endswith('_features.tsv.gz'), 'Filename']:
    source_symbols.append(set(pd.read_csv(raw / filename, sep='\t', header=None).iloc[:, 1].astype(str)))
common_symbols = set.intersection(*source_symbols)
missing_genes = sorted(set(genes) - common_symbols)
pd.Series(missing_genes, name='MissingModelGene').to_csv(out / 'missing_model_genes.tsv', sep='\t', index=False)
genes = sorted(set(genes) & common_symbols)  # 官方入口负责缺失基因；这里只输出实际有 feature 的行
print('Model gene union present:', len(genes), 'missing:', len(missing_genes), flush=True)
markers = ['TPSAB1', 'TPSB2', 'KIT', 'CPA3', 'CSF3R', 'CXCR2', 'S100A8', 'S100A9', 'CSF1R', 'MS4A7', 'MS4A1', 'CD22', 'PAX5', 'CD3D', 'CD3E', 'TCF7', 'NKG7', 'EPCAM', 'KRT19', 'CHGA', 'CHGB', 'TPH1', 'MUC2', 'TFF3', 'MKI67', 'TOP2A', 'OLFM4', 'NOS2', 'DUOX2', 'LCN2', 'CXCL1', 'CXCL2', 'COL1A1', 'DCN', 'PECAM1', 'VWF', 'RGS5', 'PDGFRB', 'SOX10', 'S100B', 'MPZ', 'CUBN', 'ALDOB', 'DEFA5', 'DEFA6']
expression_parts = []
diagnostic_parts = []
qc_parts = []
overlap = []

# Step5. 按样本提取作者已质控的细胞，计数矩阵的列顺序必须通过 barcode 精确匹配
for sample_row in samples.itertuples(index=False):
    sample = sample_row.Sample
    files = manifest[manifest.Sample.eq(sample)]
    barcodes_file = raw / files.loc[files.Filename.str.endswith('_barcodes.tsv.gz'), 'Filename'].item()
    features_file = raw / files.loc[files.Filename.str.endswith('_features.tsv.gz'), 'Filename'].item()
    matrix_file = raw / files.loc[files.Filename.str.endswith('_matrix.mtx.gz'), 'Filename'].item()
    barcodes = pd.read_csv(barcodes_file, sep='\t', header=None).iloc[:, 0].astype(str)
    features = pd.read_csv(features_file, sep='\t', header=None)
    symbols = features.iloc[:, 1].astype(str).to_numpy()
    assert barcodes.is_unique and not pd.isna(symbols).any()
    sample_cells = selected[selected.Sample.eq(sample)].copy()
    columns = pd.Index(barcodes).get_indexer(sample_cells.Barcode)
    assert (columns >= 0).all(), '作者注释与计数 barcode 必须全部匹配'
    with gzip.open(matrix_file, 'rb') as handle:
        matrix = mmread(handle).tocsc()
    assert matrix.shape == (len(features), len(barcodes))
    counts = matrix[:, columns].T.tocsr()
    del matrix
    totals = np.asarray(counts.sum(axis=1)).ravel()
    assert np.isfinite(counts.data).all() and (counts.data >= 0).all() and (totals > 0).all()
    assert np.allclose(counts.data, np.round(counts.data)), '输入必须是原始计数'
    unique_symbols, inverse = np.unique(symbols, return_inverse=True)
    collapse = sparse.csr_matrix((np.ones(len(symbols)), (np.arange(len(symbols)), inverse)), shape=(len(symbols), len(unique_symbols)))
    cpm = (sparse.diags(1e6 / totals) @ counts @ collapse).tocsr()
    qc = sample_cells[['ID', 'CellType', 'Donor', 'Sample', 'Stratum', 'FormalRecovery']].copy()
    qc['FullLibraryCounts'] = totals
    qc['DetectedFeatures'] = counts.getnnz(axis=1)
    qc_parts.append(qc)
    gene_index = pd.Index(unique_symbols).get_indexer(genes)
    assert (gene_index >= 0).all(), '缺失模型基因应先审阅，不能静默错位'
    mask = sample_cells.FormalRecovery.to_numpy()
    expression_parts.append(pd.DataFrame(cpm[mask][:, gene_index].toarray().T, index=genes, columns=sample_cells.loc[mask, 'ID']))
    marker_index = pd.Index(unique_symbols).get_indexer(markers)
    marker_values = pd.DataFrame(cpm[:, marker_index[marker_index >= 0]].toarray(), columns=np.array(markers)[marker_index >= 0], index=sample_cells.ID)
    diagnostic_parts.append(marker_values)
    for cell_type in types:
        overlap.append({'Sample': sample, 'CellType': cell_type, 'ModelGenes': len(model_genes[cell_type]), 'PresentGenes': int(np.isin(model_genes[cell_type], unique_symbols).sum())})
    print(sample, 'counts', counts.shape, 'formal', int(mask.sum()), 'library min/max', float(totals.min()), float(totals.max()), flush=True)

# Step6. 保存可核对的输入和中间摘要；正式 recovery 不纳入模型没有的谱系
expression = pd.concat(expression_parts, axis=1).loc[:, formal.ID]
expression.index.name = 'GeneSymbol'
expression.to_csv(out / 'GSE214695_Pilot18_20261009_CPM.txt', sep='\t', float_format='%.7g')
diagnostic = pd.concat(diagnostic_parts).loc[selected.ID]
diagnostic.to_csv(out / 'diagnostic_marker_CPM.tsv', sep='\t')
diagnostic.join(selected.set_index('ID')[['CellType', 'AuthorLabel', 'Donor', 'Stratum']]).to_csv(out / 'diagnostic_marker_with_annotation.tsv', sep='\t')
pd.concat(qc_parts).to_csv(out / 'selected_count_qc.tsv', sep='\t', index=False)
pd.DataFrame(overlap).to_csv(out / 'model_gene_overlap.tsv', sep='\t', index=False)
summary = {'dataset': 'GSE214695', 'discovery': 'MSCCR_CORE12_CTSPEC_REGION_IBD', 'seed': 20261009, 'donors': len(samples), 'formal_cells': len(formal), 'diagnostic_total_cells': len(selected), 'formal_types': types, 'genes': len(genes), 'cap_per_sample_type': 100, 'normalization': 'linear CPM; full-library denominator; duplicate symbols summed', 'control_limit': 'HC comes from non-tumor areas; disease tissues are inflamed; not a representative healthy population', 'sampling_limit': 'Type fractions after stratified cell sampling are not original tissue fractions; no clinical association analysis'}
(out / 'preparation_summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
print('READY:', expression.shape, flush=True)
