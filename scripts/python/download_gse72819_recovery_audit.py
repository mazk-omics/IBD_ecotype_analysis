# 通过已建立的同一SSH隧道读取服务器结果，只写本地副本，不修改服务器模型或分析源码。
# 本机 cibersortx_review 环境运行；输入为127.0.0.1:18765文件清单，输出为本轮审阅目录。
from pathlib import Path
from urllib.request import urlopen
from urllib.parse import quote
import hashlib
import json

# Step1. 保存服务器结果清单；只下载本轮明确列出的结果和模型审阅副本。
outdir = Path('D:/IBD_ecotype/IBD_ecotype_analysis/06_validation/GSE72819_recovery_review_20261009/raw_server')
outdir.mkdir(parents=True, exist_ok=True)
base_url = 'http://127.0.0.1:18765/'
with urlopen(base_url + 'audit_download_manifest.json', timeout=30) as response:
    manifest_bytes = response.read()
manifest = json.loads(manifest_bytes)
(outdir / 'audit_download_manifest.json').write_bytes(manifest_bytes)
print('Files:', len(manifest), 'Bytes:', sum(row['bytes'] for row in manifest), flush=True)

# Step2. 每个文件核对大小与SHA256，防止错传；所有请求复用同一SSH转发，不产生新登录。
for i, row in enumerate(manifest, start=1):
    destination = outdir / row['path']
    if not destination.resolve().is_relative_to(outdir.resolve()):
        raise ValueError('Manifest path is outside the audit directory')
    with urlopen(base_url + quote(row['path'], safe='/'), timeout=90) as response:
        data = response.read()
    if len(data) != row['bytes'] or hashlib.sha256(data).hexdigest() != row['sha256']:
        raise ValueError('Transfer checksum failed: ' + row['path'])
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_bytes(data)
    if i % 40 == 0 or i == len(manifest):
        print('Verified', i, '/', len(manifest), flush=True)
print('Local results:', outdir, flush=True)
