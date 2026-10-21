# 从公开EuropePMC补充资料归档取得小型表；不取原始测序队列。
# 请求无重试；归档超过20 MB即停止，不解压大文件或路径外内容。
from pathlib import Path
import urllib.request
import zipfile
import io
import hashlib
import pandas as pd
import xml.etree.ElementTree as ET

# Step1. 使用公开归档API，限制传输和解压规模
out = Path(__file__).resolve().parents[1]
logs = []
for publication, pmcid in [('Friedrich2021', 'PMC8604730'), ('Elmentaite2021', 'PMC8426186'), ('GarridoTrigo2023', 'PMC10372067'), ('Kinchen2018', 'PMC6176871')]:
    url = 'https://www.ebi.ac.uk/europepmc/webservices/rest/' + pmcid + '/supplementaryFiles'
    row = {'Publication': publication, 'URL': url}
    try:
        with urllib.request.urlopen(url, timeout=25) as response:
            declared = response.headers.get('Content-Length')
            if declared and int(declared) > 20 * 1024 * 1024:
                raise ValueError('Declared supplementary archive exceeds20MB')
            content = response.read(20 * 1024 * 1024 + 1)
        assert len(content) <= 20 * 1024 * 1024, 'Archive exceeds20MB; not saved'
        assert content[:2] == b'PK', 'Not a ZIP archive'
        (out / 'literature_resources' / (publication + '_public_supplements.zip')).write_bytes(content)
        members = []
        with zipfile.ZipFile(io.BytesIO(content)) as archive:
            for info in archive.infolist():
                members.append({'Publication': publication, 'Name': info.filename, 'Bytes': info.file_size})
                name = Path(info.filename).name
                if name.endswith(('.xlsx', '.xls')) and info.file_size <= 15 * 1024 * 1024:
                    (out / 'literature_resources' / (publication + '__' + name)).write_bytes(archive.read(info))
        pd.DataFrame(members).to_csv(out / 'audit' / (publication + '_PUBLIC_ARCHIVE_MEMBERS.tsv'), sep='\t', index=False)
        row.update({'Status': 'ARCHIVE_SAVED', 'Bytes': len(content), 'SHA256': hashlib.sha256(content).hexdigest()})
    except Exception as error:
        row.update({'Status': 'NOT_RETRIEVED', 'Error': str(error)})
    logs.append(row)
    pd.DataFrame(logs).to_csv(out / 'audit/PUBLIC_ARCHIVE_DOWNLOAD_LOG.tsv', sep='\t', index=False)
    print(publication, row['Status'], row.get('Bytes', ''), row.get('Error', ''), flush=True)

# Step2. 摘录既有BioC补充表登记，不凭空猜文件名
root = Path('D:/IBD_ecotype/IBD_ecotype_analysis/06_validation/ecotyper_original_method_annotation_20261010/sources')
records = []
for publication, pmcid in [('Martin2019', 'PMC7060942'), ('Smillie2019', 'PMC6662628')]:
    tree = ET.parse(root / (pmcid + '_bioc.xml'))
    for passage in tree.findall('.//passage'):
        text = passage.findtext('text') or ''
        metadata = {element.attrib.get('key'): element.text for element in passage.findall('infon')}
        if any(word in text.lower() for word in ['table s2', 'mmc', 'supplementary material', 'supplemental information']):
            records.append({'Publication': publication, 'Metadata': str(metadata), 'Text': text[:1500]})
pd.DataFrame(records).to_csv(out / 'audit/BIOC_SUPPLEMENT_METADATA.tsv', sep='\t', index=False)

# Step3. 保存真实表格结构，后续只按列意义解析
schemas = []
for path in (out / 'literature_resources').glob('*.xlsx'):
    workbook = pd.ExcelFile(path, engine='openpyxl')
    for sheet in workbook.sheet_names:
        preview = pd.read_excel(workbook, sheet_name=sheet, header=None, nrows=8)
        schemas.append({'File': str(path), 'Sheet': sheet, 'PreviewRowsJSON': preview.fillna('').to_json(orient='values', force_ascii=False)})
pd.DataFrame(schemas).to_csv(out / 'audit/PUBLISHED_WORKBOOK_SCHEMA.tsv', sep='\t', index=False)
print('Retrieved Excel worksheet previews:', len(schemas), flush=True)
