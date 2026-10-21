# 准备精简交付与官方上传组合，严格检查解压后总字节数不超过1GB。
# 本机运行；只打包小型表、脚本、报告与必要图件，不包括原始单细胞矩阵、大RDS和counts缓存。
from pathlib import Path
import hashlib
import json
import shutil
import zipfile
import pandas as pd

# Step1. 官方上传只需要同批CPM与一个冻结signature
project = Path(__file__).resolve().parents[2]
out = Path((project / '06_validation/fraction_sensitivity_current_path.txt').read_text(encoding='utf-8'))
upload = out / '09_CIBERSORTx_upload_small'
upload.mkdir(exist_ok=True)
signature = project / 'CIBERSORTx_ready_20261007/09_MSCCR_Joint20v2_Fractions/GJ20v2_Signature20.txt'
shutil.copyfile(signature, upload / signature.name)
files = [upload / signature.name]
for dataset in ['GSE214695', 'IBDverse']:
    path = upload / (dataset + '_pseudo_bulk100_CPM.txt')
    if path.exists():
        files.append(path)
assert len(files) == 3, '两个同批官方输入取回后才能验收上传组合'
total_upload = sum(path.stat().st_size for path in files)
assert total_upload <= 1_000_000_000
rows = []
for path in files:
    header = pd.read_csv(path, sep='\t', nrows=0).columns
    rows.append({'file': path.name, 'uncompressed_bytes': path.stat().st_size, 'numeric_columns': len(header) - 1, 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()})
pd.DataFrame(rows).to_csv(upload / 'upload_size_and_sha256.tsv', sep='\t', index=False)
assert (upload / 'README_官方fractions运行.md').exists(), '保留已有逐项网页操作说明'

# Step2. 验证结果压缩包只收必要小文件，选择性排除规则可审计
bundle = out / 'IBD_fraction_validation_compact.zip'
selected = []
for path in out.rglob('*'):
    if not path.is_file() or path == bundle or '09_CIBERSORTx_upload_small' in path.parts:
        continue
    if path.name in ['local_output_manifest.tsv', 'COMPACT_DELIVERY_SIZE_AUDIT.json']:
        continue
    if path.suffix in ['.rds', '.npz', '.npy', '.mtx']:
        continue
    if (path.name.endswith('.tsv.gz') or path.name.endswith('.txt.gz')) and path.name != 'cell_sampling_manifest.tsv.gz':
        continue
    if path.name.endswith('.collision-audit.pdf') or path.name.endswith('.alignment.svg') or path.name.endswith('.svg'):
        continue
    if path.name == 'all_author_cells.tsv' or path.stat().st_size > 20_000_000:
        continue
    selected.append(path)
uncompressed = sum(path.stat().st_size for path in selected)
assert uncompressed + total_upload <= 1_000_000_000, '本地精简结果和原始上传TXT合计不得超过1GB'
with zipfile.ZipFile(bundle, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=6) as archive:
    for path in selected:
        archive.write(path, str(path.relative_to(out)))
manifest = {'official_upload_uncompressed_bytes': total_upload, 'validation_bundle_uncompressed_bytes': uncompressed, 'validation_zip_bytes': bundle.stat().st_size, 'validation_zip_sha256': hashlib.sha256(bundle.read_bytes()).hexdigest(), 'combined_uncompressed_bytes': uncompressed + total_upload, 'limit_bytes': 1_000_000_000, 'passed': True, 'full_source_scRNA_counts_not_uploaded': True, 'large_server_intermediates_excluded': True, 'official_CIBERSORTx_predictions': 'BLOCKED; inputs only', 'size_audit_and_local_manifest': 'sidecar files outside zip to avoid recursive self-hashes'}
(out / '07_reports/COMPACT_DELIVERY_SIZE_AUDIT.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
print(json.dumps(manifest, indent=2))
