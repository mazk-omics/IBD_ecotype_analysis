# 根据已下载的 GEO 元数据和作者注释生成 GSE214695 recovery 的下载与谱系映射小表。
# 在本机项目根目录运行；只读候选元数据，写入本轮独立结果目录，不读取大表达矩阵。
from pathlib import Path
import json
import pandas as pd

# Step1. 读取本地元数据；样本标题与疾病字段交叉核对
project = Path(__file__).resolve().parents[3]
source = project / 'IBD_ecotype_analysis/06_validation/recovery_candidates_20261007/scRNA_extension'
out = project / 'IBD_ecotype_analysis/06_validation/GSE214695_scRNA_recovery_20261009'
out.mkdir(exist_ok=True)
metadata = json.loads((source / 'GSE214695_metadata.json').read_text(encoding='utf-8'))
annotation = pd.read_csv(source / 'author_annotations/GSE214695_cell_annotation.csv.gz')
sample_rows = []
url_rows = [{'Sample': 'annotation', 'Filename': 'GSE214695_cell_annotation.csv.gz', 'URL': 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE214nnn/GSE214695/suppl/GSE214695_cell_annotation.csv.gz'}]
for sample in metadata['samples']:
    name = sample['title'].split(',')[0].replace(' ', '')
    characteristics = sample['channels'][0]['characteristics']
    region = characteristics['tissue']
    disease = name.rstrip('0123456789')
    sample_rows.append({'Sample': name, 'Accession': sample['accession'], 'Disease': disease, 'Stratum': 'HC_nonIBD' if disease == 'HC' else disease + '_inflamed', 'Region': region, 'Region3': 'Rectum' if region == 'rectum' else 'Colon', 'DiseaseSource': characteristics['disease state']})
    for url in sample['supplementary_data']:
        url_rows.append({'Sample': name, 'Filename': url.split('/')[-1], 'URL': url.replace('ftp://', 'https://')})
samples = pd.DataFrame(sample_rows)
assert len(samples) == 18 and samples.Sample.is_unique and set(samples.Sample) == set(annotation['sample'])
samples.to_csv(out / 'sample_metadata.tsv', sep='\t', index=False, lineterminator='\n')
pd.DataFrame(url_rows).to_csv(out / 'download_urls.tsv', sep='\t', index=False, lineterminator='\n')

# Step2. 明确采用作者大谱系；模糊 T 细胞、增殖髓系和中性粒不强派到当前模型
mapping = {}
mapping.update(dict.fromkeys(['B cell', 'B cell Ribhi', 'GC B cell', 'Memory B cell', 'Naïve B cell'], 'B'))
mapping.update(dict.fromkeys(['CD4 ANXA1', 'CD4 naïve', 'S1PR1 T cells', 'ThF', 'Tregs'], 'CD4_T'))
mapping.update(dict.fromkeys(['CD8 CTL', 'CD8 CTL TRM', 'CD8 FGFBP2'], 'CD8_T'))
mapping.update(dict.fromkeys(['Activated endothelium', 'Endothelium', 'Lymphatic endothelium'], 'Endothelial'))
mapping.update(dict.fromkeys(['IER fibroblasts', 'MT fibroblasts', 'Inflammatory fibroblasts', 'Fibroblasts Ribhi', 'Myofibroblasts', 'S1', 'S1 Ribhi', 'S1.2', 'S2a', 'S2b', 'S3', 'FRCs'], 'Fibroblast'))
mapping.update(dict.fromkeys(['IDA macrophage', 'Inflammatory monocytes', 'M0', 'M0_Ribhi', 'M1 ACOD1', 'M1 CXCL5', 'M2', 'M2.2'], 'Monocyte_Macrophage'))
mapping.update(dict.fromkeys(['Mast 1', 'Mast 2', 'Mast Ribhi'], 'Mast'))
mapping.update(dict.fromkeys(['Goblet', 'Mature goblet'], 'Goblet'))
mapping.update({'Enteroendocrine': 'Enteroendocrine', 'Glia': 'Glial', 'ILC3': 'ILC', 'Pericytes': 'Pericyte'})
rows = []
for label, group in annotation.groupby('annotation', sort=True):
    cell_type = mapping.get(label, 'Diagnostic_' + str(group.nanostring_reference.iloc[0]).replace(' ', '_'))
    note = '作者谱系合并；不作为 EcoTyper state 注释' if label in mapping else '未纳入发现模型或谱系不明确；仅保留诊断抽样'
    if label == 'ILC3':
        note = '按作者 annotation=ILC3 映射；nanostring_reference=ILC4 不一致；不代表覆盖全部 ILC 亚型'
    rows.append({'AuthorLabel': label, 'CellType': cell_type, 'AuthorReference': group.nanostring_reference.iloc[0], 'AvailableCells': len(group), 'Note': note})
pd.DataFrame(rows).to_csv(out / 'author_label_map.tsv', sep='\t', index=False, lineterminator='\n')
print(samples[['Sample', 'Disease', 'Region']].to_string(index=False))
print('Labels:', len(rows), 'formal cells before sampling:', annotation.annotation.isin(mapping).sum())

# Step3. 下载期间先核对抽样预算和覆盖；实际随机细胞 ID 仍由服务器脚本保存
coverage = annotation[['sample', 'annotation']].copy()
coverage['CellType'] = coverage.annotation.map(mapping)
coverage = coverage.dropna(subset=['CellType'])
counts = coverage.groupby(['CellType', 'sample']).size().rename('AvailableCells').reset_index()
counts['ExpectedSampledCells'] = counts.AvailableCells.clip(upper=100)
counts.to_csv(out / 'sampling_budget_preflight.tsv', sep='\t', index=False, lineterminator='\n')
summary = counts.groupby('CellType').agg(AvailableCells=('AvailableCells', 'sum'), ExpectedSampledCells=('ExpectedSampledCells', 'sum'), SamplesWithLineage=('sample', 'nunique'))
summary.to_csv(out / 'sampling_budget_by_type.tsv', sep='\t', lineterminator='\n')
print(summary.to_string())
print('Expected formal cells:', int(counts.ExpectedSampledCells.sum()))
