# 核对新候选的作者细胞注释，下载小型标签文件，不下载表达矩阵。
# 在本机 SeqOut Python 3.13 环境运行；读 scRNA_extension 元数据，写注释原文件与摘要。
# 这里核定谱系覆盖；作者 cell type 与 EcoTyper CellType 的具体映射尚未确定。
from pathlib import Path
import hashlib
import json
import pandas as pd
import requests

# Step1. 从已读取的 GEO 元数据选取作者注释，不猜造文件下载地址。
project = Path(__file__).resolve().parents[3]
outdir = project / 'IBD_ecotype_analysis/06_validation/recovery_candidates_20261007/scRNA_extension'
annotation_dir = outdir / 'author_annotations'
annotation_dir.mkdir(exist_ok=True)
downloads = []
for accession in ['GSE214695', 'GSE266546']:
    metadata = json.loads((outdir / (accession + '_metadata.json')).read_text(encoding='utf-8'))
    for item in metadata['project'].get('supplementary_data', []):
        url = item[0] if isinstance(item, list) else item
        if 'annotation' not in url.lower():
            continue
        url = url.replace('ftp://', 'https://', 1)
        response = requests.get(url, timeout=60)
        response.raise_for_status()
        assert len(response.content) < 10_000_000, '注释超过本次小文件范围，先核对文件。'
        path = annotation_dir / url.rsplit('/', 1)[1]
        path.write_bytes(response.content)
        downloads.append({'dataset': accession, 'url': url, 'file': str(path), 'bytes': len(response.content), 'sha256': hashlib.sha256(response.content).hexdigest()})
        print('保存:', path.name, len(response.content), flush=True)
(annotation_dir / 'download_record.json').write_text(json.dumps(downloads, ensure_ascii=False, indent=2), encoding='utf-8')

# Step2. 显示作者原始列名和少量行，再决定后续类型、供者及部位统计用哪些列。
for item in downloads:
    path = Path(item['file'])
    separator = ',' if '.csv' in path.name else '\t'
    annotation = pd.read_csv(path, sep=separator)
    print('\n', path.name, annotation.shape, list(annotation.columns), flush=True)
    print(annotation.head(2).to_string(index=False), flush=True)
