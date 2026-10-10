# 只读核对官方 HiRes 公开应用层的文件清单、源码和许可文件，不运行容器或修改认证。
# 本地 cibersortx_review 环境运行；读已有镜像元数据，保存小型应用层和目录审计。
from pathlib import Path
import json
import urllib.request
import hashlib
import tarfile
import io
import pandas as pd

# Step1. 对齐构建记录和应用层，不下载大型系统依赖
project = Path(r'D:\IBD_ecotype\IBD_ecotype_analysis')
metadata = project / 'CIBERSORTx_server_plan_20261009'
output = project / 'CIBERSORTx_server_plan_20261010'
output.mkdir(exist_ok=True)
manifest = json.loads((metadata / 'hires_registry_manifest.json').read_text())
config = json.loads((metadata / 'hires_image_config.json').read_text())
history = [row for row in config['history'] if not row.get('empty_layer', False)]
assert len(history) == len(manifest['layers'])  # 防止构建步骤与层错配
selected = [(step, layer) for step, layer in zip(history, manifest['layers']) if 'COPY ' in step.get('created_by', '')]
assert sum(layer['size'] for step, layer in selected) < 5 * 1024 * 1024
print('Application COPY layers:', len(selected), 'bytes:', sum(layer['size'] for step, layer in selected))

# Step2. 匿名下载公开应用层；registry 下载凭据不是 CIBERSORTx 运行 token，不打印或保存
auth_url = 'https://auth.docker.io/token?service=registry.docker.io&scope=repository:cibersortx/hires:pull'
with urllib.request.urlopen(auth_url, timeout=30) as response:
    registry_token = json.load(response)['token']
file_rows = []
for step, layer in selected:  # 每层分别校验，仅列目录，不解压或执行归档内容
    digest = layer['digest']
    layer_path = output / (digest.replace(':', '_') + '.tar.gz')
    if not layer_path.exists():
        request = urllib.request.Request('https://registry-1.docker.io/v2/cibersortx/hires/blobs/' + digest, headers={'Authorization': 'Bearer ' + registry_token})
        with urllib.request.urlopen(request, timeout=60) as response:
            layer_path.write_bytes(response.read())
    data = layer_path.read_bytes()
    assert hashlib.sha256(data).hexdigest() == digest.split(':')[1]
    with tarfile.open(fileobj=io.BytesIO(data), mode='r:gz') as archive:
        for member in archive.getmembers():
            if member.isfile():
                file_rows.append({'Layer': digest, 'Path': member.name, 'Bytes': member.size, 'BuildStep': step.get('created_by', '')})
files = pd.DataFrame(file_rows)
files.to_csv(output / 'hires_public_application_files.tsv', sep='\t', index=False)
print(files[['Path', 'Bytes']].to_string(index=False))

# Step3. 文件名只能用于初筛；缺少许可或源码文件不能证明所有其他层都没有
source_names = files.loc[files.Path.str.contains(r'(?i)(?:\.R$|\.cpp$|\.c$|\.py$|license|readme|copying)')].copy()
source_names.to_csv(output / 'source_and_license_filename_candidates.tsv', sep='\t', index=False)
print('Source/document filename candidates:', len(source_names))
