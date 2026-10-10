# 汇总 SeqOut 单细胞候选的实际元数据及作者注释；在本机独立 Python 3.13 环境运行。
# 读取 scRNA_extension 已保存的 JSON/注释，输出样本、文件、类型覆盖和元数据疑点。
# 本脚本不重新联网、不读表达矩阵、不修改作者标签，也不运行 recovery。
from pathlib import Path
import json
import pandas as pd

# Step1. 将元数据展开成可审阅的表，明确 GEO 样本记录不是供者数。
project = Path(__file__).resolve().parents[3]
outdir = project / 'IBD_ecotype_analysis/06_validation/recovery_candidates_20261007/scRNA_extension'
study_rows = []
sample_rows = []
file_rows = []
reference = json.loads((outdir / 'GSE282122_metadata.json').read_text(encoding='utf-8'))
reference_gsm = {sample['accession'] for sample in reference['samples']}
for path in sorted(outdir.glob('*_metadata.json')):
    data = json.loads(path.read_text(encoding='utf-8'))
    study = data['project']
    gsm = {sample['accession'] for sample in data['samples']}
    study_rows.append({'accession': data['accession'], 'title': study['title'], 'design': study.get('overall_design'), 'summary': study.get('summary'), 'geo_sample_records': len(gsm), 'publication_doi': '|'.join(item['doi'] for item in study.get('publications') or [] if item.get('doi')), 'same_GSM_as_GSE282122': len(gsm & reference_gsm), 'donor_overlap_resolved': False})
    for item in study.get('supplementary_data', []):
        url = item[0] if isinstance(item, list) else item
        if isinstance(url, str) and url not in ['', 'NONE']:
            file_rows.append({'dataset': data['accession'], 'sample': '', 'url': url, 'expression_downloaded': False})
    for sample in data['samples']:
        attributes = sample['channels'][0].get('characteristics', {})
        sample_rows.append({'dataset': data['accession'], 'sample': sample['accession'], 'title': sample['title'], 'source': sample['channels'][0].get('source'), 'attributes': json.dumps(attributes, ensure_ascii=False)})
        for item in sample.get('supplementary_data', []):
            url = item[0] if isinstance(item, list) else item
            if isinstance(url, str) and url not in ['', 'NONE']:
                file_rows.append({'dataset': data['accession'], 'sample': sample['accession'], 'url': url, 'expression_downloaded': False})
studies = pd.DataFrame(study_rows)
studies.to_csv(outdir / '研究设计与GSM交集核对.tsv', sep='\t', index=False)
pd.DataFrame(sample_rows).to_csv(outdir / '重点研究样本属性.tsv', sep='\t', index=False)
pd.DataFrame(file_rows).to_csv(outdir / '重点研究文件链接_尚未下载.tsv', sep='\t', index=False)
print('保存元数据:', len(studies), '个 accession；含 GSE282122 作为交集核对参照。')

# Step2. 直接核对 GSE214695 的作者标签，不把 Paneth-like 当成正常回肠 Paneth。
annotation_dir = outdir / 'author_annotations'
colon = pd.read_csv(annotation_dir / 'GSE214695_cell_annotation.csv.gz')
colon['disease_group_from_sample'] = colon['sample'].str.extract(r'^(HC|CD|UC)', expand=False)
coverage = colon.groupby('annotation').agg(cells=('cell_id', 'size'), sample_records=('sample', 'nunique')).reset_index()
group_counts = pd.crosstab(colon.annotation, colon.disease_group_from_sample).reset_index()
coverage = coverage.merge(group_counts, on='annotation', validate='one_to_one')
coverage.to_csv(outdir / 'GSE214695_作者细胞类型及组别覆盖.tsv', sep='\t', index=False)
colon.groupby(['sample', 'disease_group_from_sample']).size().rename('cells').to_csv(outdir / 'GSE214695_18样本细胞数.tsv', sep='\t')
print('GSE214695 注释:', colon.shape, '样本:', colon['sample'].nunique(), '标签:', colon.annotation.nunique())
print(coverage.loc[coverage.annotation.isin(['Glia', 'Enteroendocrine', 'Pericytes', 'ILC3', 'Paneth-like', 'NK', 'Cycling TA'])].to_string(index=False))

# Step3. 合并 GSE266546 四个作者标签文件，记录文件名与列值不一致，禁止直接覆盖 Source。
frames = []
source_rows = []
for path in sorted(annotation_dir.glob('GSE266546_Annotation*.txt.gz')):
    cells = pd.read_csv(path, sep='\t')
    cells['annotation_file'] = path.name
    frames.append(cells)
    for source, n_cells in cells.Source.value_counts().items():
        source_rows.append({'annotation_file': path.name, 'Source_column': source, 'cells': int(n_cells)})
cd = pd.concat(frames, ignore_index=True)
cd_coverage = cd.groupby('Celltypes').agg(cells=('Barcodes', 'size'), patient_ids=('PatientID', 'nunique'), sample_ids=('SampleID', 'nunique')).reset_index()
cd_coverage.to_csv(outdir / 'GSE266546_作者类型覆盖.tsv', sep='\t', index=False)
pd.DataFrame(source_rows).to_csv(outdir / 'GSE266546_文件名与来源列疑点.tsv', sep='\t', index=False)
cd.groupby(['SampleID', 'PatientID', 'Tissue', 'Histology', 'Source', 'annotation_file']).size().rename('cells').to_csv(outdir / 'GSE266546_样本供者部位原始映射.tsv', sep='\t')
audit = {'downloaded_annotation_rows': len(cd), 'unique_barcodes': cd.Barcodes.nunique(), 'sample_ids': cd.SampleID.nunique(), 'patient_ids_as_written': cd.PatientID.nunique(), 'patient_ids_after_strip': cd.PatientID.str.strip().nunique(), 'paper_cells': 202359, 'paper_donors': 83, 'source_conflict': 'AC_Surgical filename contains only Endoscopy in Source column; unresolved', 'do_not_auto_merge_patient_ids': True, 'expression_matrices_read': False}
(outdir / 'GSE266546_注释与论文数量疑点.json').write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding='utf-8')
print('GSE266546 原始注释核对:', audit)

# Step4. 分开 CD 多部位样本和 Martin 图谱中的 PBMC 分支，避免跨组织混用。
metadata = json.loads((outdir / 'GSE290695_metadata.json').read_text(encoding='utf-8'))
multi = pd.DataFrame([{'sample': sample['accession'], 'title': sample['title'], **sample['channels'][0]['characteristics']} for sample in metadata['samples']])
multi['donor_number_from_title'] = multi.title.str.extract(r'donor\s+(\d+)', expand=False)
multi['donor_identity_verified'] = False
multi.to_csv(outdir / 'GSE290695_部位与标题供者编号.tsv', sep='\t', index=False)
print('GSE290695 文库:', len(multi), '标题供者编号:', multi.donor_number_from_title.nunique(), '部位:', multi.tissue.value_counts().to_dict())
metadata = json.loads((outdir / 'GSE134809_metadata.json').read_text(encoding='utf-8'))
martin = pd.DataFrame([{'sample': sample['accession'], 'title': sample['title'], 'source': sample['channels'][0]['source']} for sample in metadata['samples']])
martin['branch'] = 'ileal_tissue'
martin.loc[martin.title.str.contains('PBMC', case=False), 'branch'] = 'exclude_PBMC'
martin.to_csv(outdir / 'GSE134809_回肠与PBMC分支.tsv', sep='\t', index=False)
print('GSE134809 分支:', martin.branch.value_counts().to_dict())

# Step5. 将人工审阅后的分级补入总候选表；关键词误命中的 bulk 不进入这张表。
additional = pd.read_csv(outdir / '新增候选分级.tsv', sep='\t')
catalog_path = outdir.parent / '候选队列与优先级.tsv'
catalog = pd.read_csv(catalog_path, sep='\t')
catalog = pd.concat([catalog, additional], ignore_index=True).drop_duplicates('accession', keep='last')
catalog.to_csv(catalog_path, sep='\t', index=False)
print('新增单细胞候选:', len(additional), '总候选:', len(catalog))
