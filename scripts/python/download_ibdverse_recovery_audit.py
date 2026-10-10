# 经用户授权的同一 SSH 隧道回传结果；不建立 SSH/SCP/SFTP 新连接。
# 仅在现有 127.0.0.1:18765 转发和服务器临时 HTTP 服务开启时运行，写入本地 raw_server。
# 逐文件核对服务器清单的字节数与 SHA256；大表达矩阵和随机置换缓存留在服务器。
from pathlib import Path
import hashlib
import json
import urllib.request

# Step1. 读取本次隔离结果目录的文件清单
project = Path(__file__).resolve().parents[3]
out = project / 'IBD_ecotype_analysis/06_validation/IBDverse_scRNA_recovery_20261009/raw_server'
out.mkdir(exist_ok=True)
base = 'http://127.0.0.1:18765/'
manifest = json.loads(urllib.request.urlopen(base + 'download_manifest.json', timeout=30).read())
(out / 'download_manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
print('Downloading', len(manifest), 'files', flush=True)

# Step2. 在同一隧道中读取文件并核对，拒绝越过本次输出目录的路径
for item in manifest:
    target = (out / item['path']).resolve()
    assert target.is_relative_to(out.resolve())
    target.parent.mkdir(parents=True, exist_ok=True)
    data = urllib.request.urlopen(base + item['path'], timeout=90).read()
    assert len(data) == item['bytes'] and hashlib.sha256(data).hexdigest() == item['sha256']
    target.write_bytes(data)
    print(item['path'], flush=True)
print('All files verified', flush=True)
