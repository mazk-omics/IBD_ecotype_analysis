# 只读盘点旧版证据并取得公开论文XML，不下载单细胞队列。
# 本机单线程运行；每个网络资源仅请求一次，文件限制5 MB。
import os
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
from pathlib import Path
import hashlib
import json
import urllib.request
import urllib.parse
import xml.etree.ElementTree as ET
import pandas as pd

# Step1. 记录冻结来源和未提交文件的初始字节指纹
root = Path('D:/IBD_ecotype/IBD_ecotype_analysis')
out = Path(__file__).resolve().parents[1]
old = root / 'analysis_reports/legacy_bayesprism_annotation/run_20261013_01'
base = root / '06_validation/paper_aligned_recovery_characterization_20261012'
audit = root / '06_validation/all_states_lineage_audit_20261012'
paths = list(old.glob('*.tsv')) + list((audit / 'tables').glob('*.tsv')) + list((audit / 'tables').glob('*.tsv.gz'))
paths += [base / '00_audit/prior_discovery_associated_genes.tsv', base / '00_audit/Colon_GSE214695_cells.tsv', base / '00_audit/Ileum_IBDverse_cells.tsv', base / '02_external_markers/marker_composite_ranking.tsv', base / '02_external_markers/marker_LOOCV_per_dataset.tsv', base / '04_state_phenotypes/state_phenotype_enrichment.tsv', base / '04_state_phenotypes/state_region_fisher.tsv', base / '00_audit/NMF_selected_model_manifest.tsv']
status = (out / 'audit/git_status_before.txt').read_text(encoding='utf-8-sig').splitlines()
paths += [root / line[3:] for line in status if line[:2].strip() == 'M']
paths += [root / 'scripts/R/run_fraction_sensitivity_joint20_20261013.R', root / 'scripts/R/run_external_pseudobulk_joint20_20261013.R']
inventory = []
for path in dict.fromkeys(paths):
    inventory.append({'Path': str(path), 'Bytes': path.stat().st_size, 'SHA256Before': hashlib.sha256(path.read_bytes()).hexdigest(), 'Access': 'read-only'})
pd.DataFrame(inventory).to_csv(out / 'audit/INPUT_IMMUTABILITY_BASELINE.tsv', sep='\t', index=False)
schemas = []
for path in paths:
    if path.suffix == '.tsv' or path.name.endswith('.tsv.gz'):
        table = pd.read_csv(path, sep='\t', nrows=1)
        schemas.append({'Path': str(path), 'Bytes': path.stat().st_size, 'Columns': ';'.join(table.columns)})
pd.DataFrame(schemas).to_csv(out / 'audit/LOCAL_RESOURCE_SCHEMA.tsv', sep='\t', index=False)
counts_inventory = []
for dataset in ['Colon_GSE214695', 'Ileum_IBDverse']:
    for suffix in ['audit_counts.npz', 'state_donor_expression.npz', 'count_order.tsv']:
        path = audit / 'source_data' / (dataset + '_' + suffix)
        counts_inventory.append({'Dataset': dataset, 'File': str(path), 'ExistsLocally': path.exists(), 'BytesIfPresent': path.stat().st_size if path.exists() else None, 'Purpose': 'same-cell coexpression or donor-resolved expression; summary tables do not substitute'})
pd.DataFrame(counts_inventory).to_csv(out / 'audit/CELL_LEVEL_RESOURCE_INVENTORY.tsv', sep='\t', index=False)
print('Local read-only resource inventory:', len(inventory), 'hashed files', flush=True)

# Step2. 按真实DOI查询论文身份；取公开XML中的补充资料链接
papers = [('Martin2019', '10.1016/j.cell.2019.08.008'), ('Smillie2019', '10.1016/j.cell.2019.06.029'), ('Friedrich2021', '10.1038/s41591-021-01520-5'), ('Elmentaite2021', '10.1038/s41586-021-03852-1'), ('GarridoTrigo2023', '10.1038/s41467-023-40156-6'), ('Corridoni2020', '10.1038/s41591-020-1003-4'), ('Parikh2019', '10.1038/s41586-019-0992-y'), ('Kinchen2018', '10.1016/j.cell.2018.08.067')]
paper_rows = []
supplements = []
for publication, doi in papers:
    row = {'Publication': publication, 'DOI': doi, 'SearchStatus': '', 'FullTextStatus': ''}
    query = urllib.parse.urlencode({'query': 'DOI:"' + doi + '"', 'format': 'json'})
    url = 'https://www.ebi.ac.uk/europepmc/webservices/rest/search?' + query
    try:
        with urllib.request.urlopen(url, timeout=20) as response:
            data = response.read(5 * 1024 * 1024 + 1)
        assert len(data) <= 5 * 1024 * 1024
        (out / 'literature_resources' / (publication + '_record.json')).write_bytes(data)
        record = json.loads(data)['resultList']['result'][0]
        assert record.get('doi', '').lower() == doi.lower()
        row.update({'Title': record.get('title'), 'PMID': record.get('id'), 'PMCID': record.get('pmcid', ''), 'IsOpenAccess': record.get('isOpenAccess'), 'SearchStatus': 'VERIFIED_BY_DOI'})
        if not row['PMCID']:
            row['FullTextStatus'] = 'NO_PMCID'
        else:
            xml_url = 'https://www.ebi.ac.uk/europepmc/webservices/rest/' + row['PMCID'] + '/fullTextXML'
            with urllib.request.urlopen(xml_url, timeout=20) as response:
                xml = response.read(5 * 1024 * 1024 + 1)
            assert len(xml) <= 5 * 1024 * 1024
            article = ET.fromstring(xml)
            (out / 'literature_resources' / (publication + '_fulltext.xml')).write_bytes(xml)
            row['FullTextStatus'] = 'XML_SAVED'
            row['XMLSourceURL'] = xml_url
            row['XML_SHA256'] = hashlib.sha256(xml).hexdigest()
            for element in article.iter():
                if element.tag in ['supplementary-material', 'supplementary-materials']:
                    label = ' '.join(' '.join(element.itertext()).split())
                    for media in element.iter():
                        href = media.attrib.get('{http://www.w3.org/1999/xlink}href', '')
                        if href:
                            supplements.append({'Publication': publication, 'DOI': doi, 'PMCID': row['PMCID'], 'Label': label[:1500], 'Href': href, 'Tag': media.tag})
    except Exception as error:
        row['Error'] = str(error)
        if not row['SearchStatus']:
            row['SearchStatus'] = 'FAILED'
        if not row['FullTextStatus']:
            row['FullTextStatus'] = 'NOT_RETRIEVED'
    paper_rows.append(row)
    print(publication, row['SearchStatus'], row['FullTextStatus'], row.get('PMCID', ''), row.get('Error', ''), flush=True)
    pd.DataFrame(paper_rows).to_csv(out / 'audit/PUBLISHED_SIGNATURE_RESOURCE_INVENTORY.tsv', sep='\t', index=False)
    pd.DataFrame(supplements).to_csv(out / 'audit/PUBLISHED_SUPPLEMENT_LINKS.tsv', sep='\t', index=False)
