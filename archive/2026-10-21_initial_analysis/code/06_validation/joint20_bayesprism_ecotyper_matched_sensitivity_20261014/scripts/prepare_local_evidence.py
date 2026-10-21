# 核对本机主线冻结注释、完整Job16结果及9294细胞参考对应signature。
# 本机运行，只在本次独立目录保存轴哈希；不改动主线文件。
from pathlib import Path
import csv
import hashlib
import json

# Step1. 读取主线实际文件，逐项验证样本及类型顺序
project = Path('D:/IBD_ecotype/IBD_ecotype_analysis')
out = project / '06_validation/joint20_bayesprism_ecotyper_matched_sensitivity_20261014'
annotation_path = project / 'CIBERSORTx_ready_20261007/09_MSCCR_Joint20v2_Fractions/核对记录_勿上传/MSCCR_CD_UC_2025_annotation.tsv'
result_path = project / 'CIBERSORTx_signature_audit/MSCCR_Joint20v2_Job16_review_20261008/original_outputs/CIBERSORTx_Job16_Results.txt'
signature_path = project / 'CIBERSORTx_ready_20261007/09_MSCCR_Joint20v2_Fractions/GJ20v2_Signature20.txt'
with annotation_path.open(encoding='utf-8-sig') as handle:
    metadata = list(csv.DictReader(handle, delimiter='\t'))
with result_path.open(encoding='utf-8-sig') as handle:
    reader = csv.reader(handle, delimiter='\t')
    result_header = next(reader)
    result_ids = [row[0] for row in reader]
with signature_path.open(encoding='utf-8-sig') as handle:
    reader = csv.reader(handle, delimiter='\t')
    types = next(reader)[1:]
    genes = [row[0] for row in reader]
ids = [row['ID'] for row in metadata]
assert ids == result_ids and len(ids) == len(set(ids))
assert types == result_header[1:21] and len(types) == len(set(types)) == 20
assert len(genes) == len(set(genes)) == 4474
assert set(row['Disease'] for row in metadata) == {'CD', 'UC'}
assert all(row[field] for row in metadata for field in ['ID', 'PatientID', 'Disease', 'Inflammation', 'Region3'])

# Step2. 保存字节和顺序证据，供同一SSH连接上传服务器核验
files = []
for path in [annotation_path, result_path, signature_path]:
    files.append({'path': str(path), 'bytes': path.stat().st_size, 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()})
snapshot = {'samples': len(ids), 'disease_counts': {d: sum(row['Disease'] == d for row in metadata) for d in ['CD', 'UC']}, 'ordered_sample_sha256': hashlib.sha256(('\n'.join(ids) + '\n').encode()).hexdigest(), 'metadata_sha256': files[0]['sha256'], 'types': types, 'signature_genes': len(genes), 'Job16_order_exact': True, 'files': files, 'full_fractions_status': 'Job16 full 2025 NS fractions; corrected Job20 S-mode is 120-sample development subset; full paired HiRes and discovery design not frozen'}
(out / '00_audit').mkdir(exist_ok=True)
(out / '00_audit/local_mainline_evidence.json').write_text(json.dumps(snapshot, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
print(json.dumps(snapshot, ensure_ascii=False, indent=2))
