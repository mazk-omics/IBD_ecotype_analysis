# 核对 Docker Hub 官方 CIBERSORTx 镜像的公开标签、清单和启动配置，不下载镜像层、不运行反卷积。
# 在本机 Python 环境运行；输出到服务器重构方案目录。此处的匿名 registry 凭据只用于公开镜像查询，并非 CIBERSORTx 运行 token。
import datetime
import hashlib
import json
import urllib.error
import urllib.request
from pathlib import Path

# Step1. 保存查询位置和日期
outdir = Path('D:/IBD_ecotype/IBD_ecotype_analysis/CIBERSORTx_server_plan_20261009')
outdir.mkdir(parents=True, exist_ok=True)
records = []
print('Query public registry metadata only; no image layers or runtime credentials are saved.')

# Step2. 逐个核对官方三个模块，保留清单及配置以便固定版本
# 模块各有独立清单；一个模块查询失败不应被误写成其他模块也不可下载。
for module in ['hires', 'fractions', 'gep']:
    record = {'repository': f'cibersortx/{module}', 'queried_at_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'image_layers_downloaded': False, 'runtime_tested': False}
    try:
        request = urllib.request.Request(f'https://hub.docker.com/v2/repositories/cibersortx/{module}/tags/latest/', headers={'User-Agent': 'IBD-EcoTyper-image-audit'})
        with urllib.request.urlopen(request, timeout=20) as response:
            tag = json.load(response)
        (outdir / f'{module}_dockerhub_latest_tag.json').write_text(json.dumps(tag, ensure_ascii=False, indent=2), encoding='utf-8')
        record['hub_last_updated'] = tag.get('last_updated')
        record['hub_images'] = [{'os': image.get('os'), 'architecture': image.get('architecture'), 'digest': image.get('digest'), 'compressed_bytes': image.get('size')} for image in tag.get('images', [])]
    except (urllib.error.URLError, TimeoutError, ValueError) as error:
        record['hub_query_error'] = str(error)

    try:
        with urllib.request.urlopen(f'https://auth.docker.io/token?service=registry.docker.io&scope=repository:cibersortx/{module}:pull', timeout=20) as response:
            registry_access = json.load(response)['token']
        request = urllib.request.Request(f'https://registry-1.docker.io/v2/cibersortx/{module}/manifests/latest', headers={'Authorization': f'Bearer {registry_access}', 'Accept': 'application/vnd.docker.distribution.manifest.v2+json'})
        with urllib.request.urlopen(request, timeout=20) as response:
            manifest_bytes = response.read()
            digest = response.headers.get('Docker-Content-Digest')
        manifest = json.loads(manifest_bytes)
        calculated_digest = 'sha256:' + hashlib.sha256(manifest_bytes).hexdigest()
        if digest and digest != calculated_digest:
            raise ValueError('Registry manifest digest does not match received content.')
        (outdir / f'{module}_registry_manifest.json').write_bytes(manifest_bytes)
        record['manifest_digest'] = calculated_digest
        record['anonymous_manifest_download'] = True
        record['compressed_layer_bytes'] = sum(layer['size'] for layer in manifest.get('layers', []))

        config_digest = manifest['config']['digest']
        request = urllib.request.Request(f'https://registry-1.docker.io/v2/cibersortx/{module}/blobs/{config_digest}', headers={'Authorization': f'Bearer {registry_access}'})
        with urllib.request.urlopen(request, timeout=20) as response:
            config_bytes = response.read()
        if 'sha256:' + hashlib.sha256(config_bytes).hexdigest() != config_digest:
            raise ValueError('Registry configuration digest does not match received content.')
        config = json.loads(config_bytes)
        (outdir / f'{module}_image_config.json').write_bytes(config_bytes)
        record['config_digest'] = config_digest
        record['architecture'] = config.get('architecture')
        record['os'] = config.get('os')
        record['entrypoint'] = config.get('config', {}).get('Entrypoint')
        record['command'] = config.get('config', {}).get('Cmd')
        record['image_created'] = config.get('created')
    except (urllib.error.URLError, TimeoutError, ValueError, KeyError) as error:
        record['registry_query_error'] = str(error)

    records.append(record)
    (outdir / 'official_image_metadata_audit.json').write_text(json.dumps(records, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps(record, ensure_ascii=False))

# Step3. 保留结论边界：公开清单可读取，不等于服务器已拉取镜像或获得算法运行授权
print('Metadata audit completed. Server installation, container help, and authenticated execution remain untested.')
