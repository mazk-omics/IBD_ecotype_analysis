# 单一SSH连接取回小型结果；大单细胞数据、counts缓存和RDS不进入结果传输包。
# 服务器ibd_refaudit运行；列出全部服务器对象路径，压缩必要数值表、日志和脚本。
from pathlib import Path
import hashlib
import json
import sys
import zipfile
import pandas as pd

# Step1. 全部新输出保存路径和大小，重点结果计算哈希
out = Path(sys.argv[1])
transfer = out / 'transfer'
transfer.mkdir(exist_ok=True)
records = []
selected = []
for path in sorted(out.rglob('*')):
    if not path.is_file() or transfer in path.parents:
        continue
    relative = str(path.relative_to(out))
    small = path.suffix in ['.tsv', '.json', '.txt', '.log', '.md', '.py', '.R', '.sh', '.exit_code', '.sha256'] and path.stat().st_size < 15_000_000
    # 完整取样manifest本身仅数MB，保留精确抽样而不带大细胞表达。
    manifest = path.name == 'cell_sampling_manifest.tsv.gz'
    digest = ''
    if small or manifest:
        with path.open('rb') as handle:
            digest = hashlib.file_digest(handle, 'sha256').hexdigest()
        selected.append(path)
    records.append({'relative_path': relative, 'bytes': path.stat().st_size, 'sha256': digest, 'server_path': str(path), 'included_in_small_transfer': small or manifest})
inventory = out / '00_audit/server_output_manifest.tsv'
pd.DataFrame(records).to_csv(inventory, sep='\t', index=False)
if inventory not in selected:
    selected.append(inventory)

# Step2. 官方CPM只记精确压缩和未压缩大小，后续同连接分块取回
uploads = []
for dataset in ['GSE214695', 'IBDverse']:
    path = out / '08_external_pseudobulk' / dataset / '05_cibersortx_official_pending/pseudo_bulk100_CPM.txt.gz'
    if path.exists():
        import gzip
        raw_bytes = 0
        with gzip.open(path, 'rb') as handle:
            for block in iter(lambda: handle.read(1_048_576), b''):
                raw_bytes += len(block)
        with path.open('rb') as handle:
            digest = hashlib.file_digest(handle, 'sha256').hexdigest()
        uploads.append({'dataset': dataset, 'server_path': str(path), 'compressed_bytes': path.stat().st_size, 'uncompressed_bytes': raw_bytes, 'sha256': digest})
assert sum(x['uncompressed_bytes'] for x in uploads) < 1_000_000_000
upload_manifest = out / '00_audit/official_upload_mixture_sizes.tsv'
pd.DataFrame(uploads).to_csv(upload_manifest, sep='\t', index=False)
selected.append(upload_manifest)

# Step3. 只传小结果，单次验证结果包不超过1GB
bundle = transfer / 'fraction_validation_small_results.zip'
with zipfile.ZipFile(bundle, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=6) as archive:
    for path in dict.fromkeys(selected):
        archive.write(path, str(path.relative_to(out)))
uncompressed = sum(path.stat().st_size for path in dict.fromkeys(selected))
assert uncompressed < 1_000_000_000
with bundle.open('rb') as handle:
    digest = hashlib.file_digest(handle, 'sha256').hexdigest()
print(json.dumps({'bundle': str(bundle), 'bytes': bundle.stat().st_size, 'sha256': digest, 'uncompressed_bytes': uncompressed, 'mixture_files': uploads}, indent=2))
