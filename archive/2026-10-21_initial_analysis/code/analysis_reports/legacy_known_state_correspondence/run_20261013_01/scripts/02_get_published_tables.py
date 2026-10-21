# 只下载已定位的公开小型补充签名表，单个文件上限20 MB。
# 不下载新的单细胞队列，不修改软件或历史结果。
from pathlib import Path
import hashlib
import json
import urllib.request
import urllib.parse
import zipfile
import io
import xml.etree.ElementTree as ET
import pandas as pd

# Step1. 核实过的公开补充表链接；原文marker统计作为完整来源保留
out = Path(__file__).resolve().parents[1]
links = pd.read_csv(out / 'audit/PUBLISHED_SUPPLEMENT_LINKS.tsv', sep='\t').drop_duplicates(['Publication', 'Href'])
chosen = links[(links.Publication.eq('Friedrich2021') & links.Href.isin(['41591_2021_1520_MOESM8_ESM.xlsx', '41591_2021_1520_MOESM11_ESM.xlsx'])) | (links.Publication.eq('Elmentaite2021') & links.Href.eq('41586_2021_3852_MOESM9_ESM.xls')) | (links.Publication.eq('GarridoTrigo2023') & links.Href.eq('41467_2023_40156_MOESM4_ESM.zip'))].copy()
download_rows = []
for table in chosen.itertuples():
    url = 'https://static-content.springer.com/esm/art%3A' + urllib.parse.quote(table.DOI, safe='') + '/MediaObjects/' + table.Href
    row = {'Publication': table.Publication, 'DOI': table.DOI, 'Label': table.Label, 'URL': url, 'File': table.Href}
    try:
        with urllib.request.urlopen(url, timeout=25) as response:
            declared = response.headers.get('Content-Length')
            if declared and int(declared) > 20 * 1024 * 1024:
                raise ValueError('DECLARED_SIZE_EXCEEDS20MB; no download')
            content = response.read(20 * 1024 * 1024 + 1)
        assert len(content) <= 20 * 1024 * 1024
        assert content[:2] == b'PK' or content[:4] == bytes.fromhex('D0CF11E0'), 'Not an Excel/ZIP binary'
        target = out / 'literature_resources' / table.Href
        target.write_bytes(content)
        row.update({'Status': 'DOWNLOADED', 'Bytes': len(content), 'SHA256': hashlib.sha256(content).hexdigest()})
        if target.suffix == '.zip':
            with zipfile.ZipFile(io.BytesIO(content)) as archive:
                members = [{'File': info.filename, 'BytesUncompressed': info.file_size} for info in archive.infolist()]
                pd.DataFrame(members).to_csv(out / 'audit/GarridoTrigo2023_SOURCE_DATA_ZIP_INVENTORY.tsv', sep='\t', index=False)
    except Exception as error:
        row.update({'Status': 'NOT_RETRIEVED', 'Error': str(error)})
    download_rows.append(row)
    print(table.Publication, table.Href, row['Status'], row.get('Bytes', ''), row.get('Error', ''), flush=True)
    pd.DataFrame(download_rows).to_csv(out / 'audit/PUBLISHED_TABLE_DOWNLOAD_LOG.tsv', sep='\t', index=False)

# Step2. 复用本地BioC全文，不重复下载；提取候选相关原文段落
root = Path('D:/IBD_ecotype/IBD_ecotype_analysis')
sources = root / '06_validation/ecotyper_original_method_annotation_20261010/sources'
passages = []
for publication, pmcid in [('Martin2019', 'PMC7060942'), ('Smillie2019', 'PMC6662628')]:
    path = sources / (pmcid + '_bioc.xml')
    tree = ET.parse(path)
    for index, passage in enumerate(tree.findall('.//passage')):
        text = passage.findtext('text') or ''
        if any(word in text.lower() for word in ['table s2', 'table s3', 'pla2g2d', 'il23r', 'goblet', 'dendritic', 'cycling', 'naive', 'resident memory']):
            passages.append({'Publication': publication, 'PMCID': pmcid, 'SourcePath': str(path), 'PassageIndex0': index, 'Text': text})
for publication in ['Friedrich2021', 'Elmentaite2021', 'GarridoTrigo2023', 'Kinchen2018']:
    path = out / 'literature_resources' / (publication + '_fulltext.xml')
    tree = ET.parse(path)
    for index, passage in enumerate(tree.findall('.//body//p')):
        text = ' '.join(' '.join(passage.itertext()).split())
        if any(word in text.lower() for word in ['pla2g2d', 'il23r', 'goblet', 'dendritic', 'cycling', 'marker gene', 'ccl19', 'ccl21', 'pdgfra', 'abca8']):
            passages.append({'Publication': publication, 'SourcePath': str(path), 'PassageIndex0': index, 'Text': text})
pd.DataFrame(passages).to_csv(out / 'audit/VERIFIED_PRIMARY_PASSAGES.tsv', sep='\t', index=False)
print('Relevant primary passages saved:', len(passages), flush=True)

# Step3. 读取.xlsx工作簿结构；不改原工作簿
schemas = []
for path in (out / 'literature_resources').glob('*.xlsx'):
    workbook = pd.ExcelFile(path, engine='openpyxl')
    for sheet in workbook.sheet_names:
        preview = pd.read_excel(workbook, sheet_name=sheet, header=None, nrows=8)
        schemas.append({'File': str(path), 'Sheet': sheet, 'PreviewRowsJSON': preview.fillna('').to_json(orient='values', force_ascii=False)})
pd.DataFrame(schemas).to_csv(out / 'audit/PUBLISHED_WORKBOOK_SCHEMA.tsv', sep='\t', index=False)
print('Excel worksheet previews:', len(schemas), flush=True)
