# 核验 MSCCR 回传文件哈希，解压得到网页 TXT，并逐行检查格式与样本/基因顺序。
# 在本机项目环境运行；读取本轮回传目录，写入同目录，不改来源或 signature。
# 原始到 CPM 的数值核对已在服务器完成，本脚本确认本机交付与服务器文件一致。
import gzip
import hashlib
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

# Step1. 核对压缩文件，按服务器哈希确认完整回传
root = Path('D:/IBD_ecotype/IBD_ecotype_analysis/CIBERSORTx_signature_audit/MSCCR_web_input_adaptation_20261006')
reports = root / 'server_outputs/audit/msccr_ibd_joint_web_20261006'
summary = json.loads((reports / 'summary.json').read_text(encoding='utf-8'))
compressed = root / 'web_upload/MSCCR_CD_UC_2025_GJ20_CPM.txt.gz'
web = compressed.with_suffix('')
compressed_hash = hashlib.file_digest(compressed.open('rb'), 'sha256').hexdigest()
assert compressed_hash == summary['files'][compressed.name]['sha256']
with gzip.open(compressed, 'rb') as reader, web.open('wb') as writer:
    shutil.copyfileobj(reader, writer, length=4 * 1024 * 1024)
web_hash = hashlib.file_digest(web.open('rb'), 'sha256').hexdigest()
assert web_hash == summary['files'][web.name]['sha256']
print('本机 GZ/TXT SHA256 均与服务器一致:', web.stat().st_size, 'bytes', flush=True)

# Step2. 按已有注释和身份清单检查每列、每行，以及所有数值格式
annotation = pd.read_csv(reports / 'MSCCR_CD_UC_2025_annotation.tsv', sep='\t')
identity = pd.read_csv(reports / 'exported_gene_identity.tsv', sep='\t')
genes = []
maximum = 0.0
with web.open() as handle:
    header = handle.readline().rstrip('\n').split('\t')
    assert header[0] == 'GeneSymbol' and header[1:] == annotation.ID.tolist()
    for row_number, line in enumerate(handle, 1):
        assert '"' not in line
        fields = line.rstrip('\n').split('\t', 1)
        genes.append(fields[0])
        assert fields[1].count('\t') == len(annotation) - 1 and '\t\t' not in fields[1]
        values = np.fromstring(fields[1], sep='\t')
        assert len(values) == len(annotation) and np.isfinite(values).all() and (values >= 0).all(), row_number
        maximum = max(maximum, float(values.max()))
assert genes == identity.gene_symbol.tolist() and len(genes) == 17040 and len(annotation) == 2025
assert len(genes) * len(annotation) == 34506000 and maximum >= 50

# Step3. 保存交付核验状态，不把网页结果或 EcoTyper 写成已完成
report = {'status': 'local_delivery_verified', 'gzip_sha256': compressed_hash, 'txt_sha256': web_hash,
    'genes': len(genes), 'samples': len(annotation), 'data_points': len(genes) * len(annotation), 'max_exported_CPM': maximum,
    'all_values_finite_nonnegative': True, 'no_double_quotes_or_missing_entries': True, 'gene_and_sample_order_match_server': True,
    'MSCCR_CIBERSORTx_submitted': False, 'signature_accuracy_validated': False, 'EcoTyper_run': False}
(root / 'local_delivery_verification.json').write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps(report, ensure_ascii=False, indent=2), flush=True)
