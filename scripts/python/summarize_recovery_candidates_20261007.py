# 汇总本机 SeqOut 检索与元数据，生成单细胞/空间 recovery 候选的样本清单。
# 使用 .tools/seqout/environment 的 Python 3.13；输出到 06_validation/recovery_candidates_20261007。
# 这是数据来源与输入准备的盘点，不计算 EcoTyper recovery，也不把候选数当作独立患者数。
from pathlib import Path
from importlib.metadata import version
import json
import pandas as pd
import requests

# Step1. 保留原始检索结果，并生成便于人工筛选的标题索引。
project = Path(__file__).resolve().parents[3]
outdir = project / 'IBD_ecotype_analysis/06_validation/recovery_candidates_20261007'
search_rows = []
for path in sorted(outdir.glob('*_search.json')):
    hits = json.loads(path.read_text(encoding='utf-8'))
    for hit in hits:
        search_rows.append({'query_file': path.name, 'accession': hit.get('accession'), 'title': hit.get('title'), 'source': hit.get('source'), 'doi': hit.get('doi'), 'is_single_cell': hit.get('is_single_cell')})
search_index = pd.DataFrame(search_rows)
# IBDverse 被服务纠正成 adverse；这些原始命中留作诊断，不能纳入候选发现。
search_index['query_status'] = 'title_review_required'
search_index.loc[search_index.query_file.str.startswith('ibdverse'), 'query_status'] = 'excluded_query_corrected_to_adverse'
search_index.to_csv(outdir / '检索标题索引.tsv', sep='\t', index=False)
search_index.loc[search_index.query_status.ne('excluded_query_corrected_to_adverse')].to_csv(outdir / '有效查询标题索引_仍需人工筛选.tsv', sep='\t', index=False)
print('检索索引行数:', len(search_index))

# Step2. 展开 GEO 样本和下载链接；不下载链接指向的表达或图像文件。
sample_rows = []
file_rows = []
for path in sorted(outdir.glob('*_seqout_metadata.json')):
    data = json.loads(path.read_text(encoding='utf-8'))
    for url in data['project'].get('supplementary_data', []):
        if isinstance(url, str) and url not in ['NONE', '']:
            file_rows.append({'dataset': data['accession'], 'sample': '', 'url': url, 'downloaded': False})
    for sample in data['samples']:
        channels = sample.get('channels') or []
        characteristics = channels[0].get('characteristics', {}) if channels else {}
        row = {'dataset': data['accession'], 'accession': sample['accession'], 'title': sample['title'], 'description': sample.get('description'), 'source': channels[0].get('source') if channels else None}
        row.update(characteristics)
        sample_rows.append(row)
        for url in sample.get('supplementary_data', []):
            if isinstance(url, str) and url not in ['NONE', '']:
                file_rows.append({'dataset': data['accession'], 'sample': sample['accession'], 'url': url, 'downloaded': False})
samples = pd.DataFrame(sample_rows)
samples.to_csv(outdir / 'GEO候选样本元数据.tsv', sep='\t', index=False)
pd.DataFrame(file_rows).to_csv(outdir / '待下载文件链接_仅清单.tsv', sep='\t', index=False)
print('已整理样本记录:', len(samples), '文件链接:', len(file_rows))

# Step3. 固定低分辨率候选的疾病分支，不把治疗相关结肠炎混入 UC。
visium = samples.loc[samples.dataset.eq('GSE189184')].copy()
visium['proposed_branch'] = 'separate_ICI_colitis'
visium.loc[visium.disease.eq('Ulcerative Colitis'), 'proposed_branch'] = 'primary_UC'
visium.loc[visium.disease.eq('Healthy Control'), 'proposed_branch'] = 'primary_healthy'
visium['patient_id_verified'] = False
visium.to_csv(outdir / 'GSE189184_16切片候选分组.tsv', sep='\t', index=False)
visium.loc[visium.proposed_branch.isin(['primary_UC', 'primary_healthy'])].to_csv(outdir / 'GSE189184_主候选9切片.tsv', sep='\t', index=False)
assert visium.proposed_branch.value_counts().to_dict() == {'primary_UC': 7, 'separate_ICI_colitis': 7, 'primary_healthy': 2}
print('GSE189184 分组:', visium.proposed_branch.value_counts().to_dict())

ileum = samples.loc[samples.dataset.eq('GSE228360')].copy()
ileum['assay_branch'] = ileum.title.str.extract(r'\[(Sc\d+|ST\d+)\]', expand=False).str.replace(r'\d+', '', regex=True)
ileum['patient_id_verified'] = False
ileum.to_csv(outdir / 'GSE228360_单细胞与空间样本分开.tsv', sep='\t', index=False)
print('GSE228360 测量分支:', ileum.assay_branch.value_counts().to_dict())

cosmx = samples.loc[samples.dataset.eq('GSE312415')].copy()
cosmx['panel_from_title'] = cosmx.title.str.extract(r'(?i)(CosMx1k|CosMx6k)', expand=False)
cosmx['proposed_branch'] = 'technical_or_other_pending_review'
cosmx.loc[cosmx.panel_from_title.str.lower().eq('cosmx1k'), 'proposed_branch'] = 'biological_1k_candidate'
cosmx.loc[cosmx.panel_from_title.str.lower().eq('cosmx6k'), 'proposed_branch'] = 'biological_6k_priority'
cosmx['independent_donor_verified'] = False
cosmx.to_csv(outdir / 'GSE312415_面板与技术重复待审清单.tsv', sep='\t', index=False)
print('CosMx 标题明确的面板:', cosmx.panel_from_title.value_counts().to_dict())

# Step4. 延续已有供者隔离，不对 IBDverse 再随机划分；这里仅保存旧留出名单。
pilot = project / 'IBD_ecotype_analysis/IBDverse_reference/pilot_20261006_172734'
donors = pd.read_csv(pilot / 'donor_split.tsv', sep='\t')
holdout = donors.loc[donors.split.eq('holdout')].copy()
holdout.to_csv(outdir / 'IBDverse_既有68供者留出名单_勿重抽.tsv', sep='\t', index=False)
donor_counts = pd.crosstab(donors.stratum, donors.split)
donor_counts.to_csv(outdir / 'IBDverse_既有供者分层.tsv', sep='\t')
print('IBDverse 既有供者分层:\n', donor_counts)

# Step5. 核对类型名称覆盖；标签缺失不代表生物学缺失，也不证明任何状态恢复成功。
pilot_summary = json.loads((pilot / 'summary.json').read_text(encoding='utf-8'))
ibdverse_types = set(pilot_summary['web_reference_types'])
joint_file = project / 'IBD_ecotype_analysis/CIBERSORTx_ready_20261007/03_MSCCR_Joint20_exploratory/GJ20_Signature20.txt'
joint_types = set(joint_file.open(encoding='utf-8').readline().strip().split('\t')[1:])
core12 = {'B', 'CD4_T', 'CD8_T', 'Endothelial', 'Enteroendocrine', 'Fibroblast', 'Glial', 'Goblet', 'ILC', 'Mast', 'Monocyte_Macrophage', 'Pericyte'}
type_check = {'ibdverse_pilot_types': sorted(ibdverse_types), 'joint20_missing_named_types': sorted(joint_types - ibdverse_types), 'joint20_extra_ibdverse_named_types': sorted(ibdverse_types - joint_types), 'core12_missing_named_types': sorted(core12 - ibdverse_types), 'state_gene_coverage_computed': False, 'comment': '仅按已审阅的试验类型名称比较；实际 recovery 必须回到作者标签及最终 discovery 模型，不据此补造细胞类型。'}
(outdir / 'IBDverse_类型覆盖边界.json').write_text(json.dumps(type_check, ensure_ascii=False, indent=2), encoding='utf-8')
print('Joint20 当前缺少同名注释:', type_check['joint20_missing_named_types'])

# Step6. 在 SeqOut 未正确检出的情况下，保存官方 BioStudies 元数据作为补充来源。
url = 'https://www.ebi.ac.uk/biostudies/api/v1/studies/E-MTAB-16999'
response = requests.get(url, timeout=30)
official_record = {'url': url, 'http_status': response.status_code}
if response.status_code == 200:
    official = response.json()
    (outdir / 'E-MTAB-16999_official_biostudies.json').write_text(json.dumps(official, ensure_ascii=False, indent=2), encoding='utf-8')
official_record['seqout_unquoted_query_corrected_to_adverse'] = True
official_record['seqout_quoted_exact_query_result'] = 'No results; does not establish that dataset is unavailable'
official_record['seqout_version'] = version('seqout')
(outdir / '补充来源与SeqOut检索限制.json').write_text(json.dumps(official_record, ensure_ascii=False, indent=2), encoding='utf-8')
print('BioStudies 官方补充来源 HTTP:', response.status_code)
