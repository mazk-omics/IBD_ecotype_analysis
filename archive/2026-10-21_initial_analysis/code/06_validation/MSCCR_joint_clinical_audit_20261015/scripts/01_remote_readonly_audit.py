# 通过唯一持久 SSH 会话在服务器内存中运行，读取 MSCCR 资产及小型临床文件。
# 不写服务器文件；压缩 JSON 经原会话返回本机，表达矩阵只读取表头及少量行。
from pathlib import Path
import os, json, hashlib, gzip, csv, base64, zlib, stat

# Step1. 在项目和上传目录内按名称限定检索，完整保留清单及权限错误
root = Path('/home/mazekai/IBD_EcoTyper')
inventory = []
errors = []
for base in [root, Path('/home/mazekai/uploads')]:
    for directory, dirs, files in os.walk(base, onerror=lambda e: errors.append(str(e))):
        dirs[:] = [d for d in dirs if d not in ['.git', '__pycache__', 'node_modules']]
        for name in files:
            p = Path(directory) / name
            relative = str(p.relative_to(base)).lower()
            relevant = any(x in relative for x in ['msccr', '193677', 'ctspec_region_ibd', '02_metadata'])
            candidate = any(x in name.lower() for x in ['clinical', 'bmis', 'supplement', 'metadata', 'annotation'])
            if not (relevant or candidate):
                continue
            s = p.stat()
            if not stat.S_ISREG(s.st_mode):
                continue
            row = dict(path=str(p), bytes=s.st_size, modified=s.st_mtime, suffix=p.suffix, candidate_only=not relevant)
            if relevant and s.st_size < 3000000 and (name in ['state_abundances.txt', 'ecotype_abundance.txt', 'ecotypes.txt', 'ecotype_assignment.txt'] or 'annotation.tsv' in name):
                row['sha256'] = hashlib.sha256(p.read_bytes()).hexdigest()
            inventory.append(row)

# Step2. 小型权威输入和 schema 证据返回原文，保留来源及校验值
selected = [
 '05_EcoTyper/input/MSCCR_CORE12/annotation/MSCCR_CORE12_annotation.tsv',
 '05_EcoTyper/input/MSCCR_CORE12/region_ibd/annotation/MSCCR_CORE12_IBD_annotation.tsv',
 '05_EcoTyper/input/MSCCR_CORE12/manifests/core12_Glike_QC.tsv',
 '05_EcoTyper/input/MSCCR_CORE12/region_ibd/input_QC.tsv',
 '03_reference/combined_reference/output/08_msccr_bulk_audit/MSCCR_bulk_schema_audit.json',
 '03_reference/combined_reference/output/21_bayesprism_epithelial_qc/GSE193677_metadata_used_for_stage2_qc.tsv',
 '03_reference/web_reference_three_pilots/msccr_input_audit_20261006/MSCCR_sample_library_audit.tsv']
texts = {}
for relative in selected:
    p = root / relative
    if p.exists():
        texts[str(p)] = p.read_text()
soft = root / '03_reference/combined_reference/output/15_bayesprism_broad14_sanity_audit/GSE193677_family.soft.gz'
if soft.exists():
    with gzip.open(soft, 'rt') as handle:
        texts[str(soft)] = ''.join(line for line in handle if line.startswith(('^SAMPLE', '!Sample_title', '!Sample_characteristics_ch1', '!Series_')))

# Step3. 大矩阵只读表头和两行；维度以既有 schema 与 QC 交叉核对
headers = {}
for row in inventory:
    p = Path(row['path'])
    if row['candidate_only'] or p.stat().st_size < 3000000 or p.suffix not in ['.txt', '.gz']:
        continue
    if 'expression_matrices' in str(p) or 'Biopsy_counts' in p.name or 'Biopsy_adjcounts' in p.name or 'GSE193677_MSCCR' in p.name:
        opener = gzip.open if p.suffix == '.gz' else open
        with opener(p, 'rt') as handle:
            header = handle.readline().strip()
            first = handle.readline().strip()
            second = handle.readline().strip()
        delimiter = '\t' if '\t' in header else ' '
        columns = next(csv.reader([header], delimiter=delimiter, skipinitialspace=True))
        first_fields = next(csv.reader([first], delimiter=delimiter, skipinitialspace=True))
        headers[str(p)] = dict(header_fields=columns, first_row_width=len(first_fields), first_gene=first_fields[0], second_gene=second.split()[0], first_values=first_fields[1:6])
payload = json.dumps(dict(inventory=inventory, texts=texts, headers=headers, errors=errors), ensure_ascii=True).encode()
encoded = base64.b64encode(zlib.compress(payload)).decode()
print('AUDIT_READY', len(encoded), len(inventory), len(texts), len(headers), flush=True)
