# 本机下载并校验 Ghidra 和配套 JDK，用于静态读取 Linux 程序；不运行 CIBERSORTx。
# cibersortx_review 环境运行，工具放在项目 .tools，版本和摘要写入源码审计目录。
from pathlib import Path
import urllib.request
import json
import hashlib
import zipfile

# Step1. 固定本轮已查询的官方发行版；只安装项目内的便携工具，不改系统 Java
project = Path(__file__).resolve().parents[2]
tool_dir = project.parent / '.tools' / 'cibersortx_static'
tool_dir.mkdir(parents=True, exist_ok=True)
outdir = project / 'CIBERSORTx_server_plan_20261011'
records = []
assets = [
    {'repository': 'NationalSecurityAgency/ghidra', 'release': 'Ghidra_12.1.4_build', 'name': 'ghidra_12.1.4_PUBLIC_20260921.zip', 'size': 569649598, 'digest': 'sha256:ddac49f903da9d5bac833e5cc79395098b9c33cfd3279be5f31bd00387d2d4db', 'browser_download_url': 'https://github.com/NationalSecurityAgency/ghidra/releases/download/Ghidra_12.1.4_build/ghidra_12.1.4_PUBLIC_20260921.zip'},
    {'repository': 'adoptium/temurin25-binaries', 'release': 'jdk-25.0.4.1+1', 'name': 'OpenJDK25U-jdk_x64_windows_hotspot_25.0.4.1_1.zip', 'size': 141167264, 'digest': '', 'browser_download_url': 'https://github.com/adoptium/temurin25-binaries/releases/download/jdk-25.0.4.1%2B1/OpenJDK25U-jdk_x64_windows_hotspot_25.0.4.1_1.zip'}
]

for asset in assets:  # 两个不同官方发行包分别下载、校验
    repository = asset['repository']
    expected = asset['digest']
    if not expected:
        with urllib.request.urlopen(asset['browser_download_url'] + '.sha256.txt', timeout=30) as response:
            expected = 'sha256:' + response.read().decode().split()[0]

    # Step2. 大文件流式下载，保留明确进度；无自动重试或未经校验的解压
    target = tool_dir / asset['name']
    if not target.exists():
        partial = target.with_suffix('.download')
        print('Downloading:', repository, asset['release'], asset['size'], flush=True)
        total = 0
        progress = 0
        with urllib.request.urlopen(asset['browser_download_url'], timeout=60) as response, partial.open('wb') as stream:
            while True:
                chunk = response.read(1024 * 1024)
                if not chunk:
                    break
                stream.write(chunk)
                total += len(chunk)
                if total - progress >= 50 * 1024 * 1024:
                    print('Downloaded MB:', round(total / 1024 / 1024), flush=True)
                    progress = total
        assert total == asset['size'], 'Download size mismatch'
        partial.rename(target)
    digest = hashlib.sha256()
    with target.open('rb') as stream:
        while chunk := stream.read(1024 * 1024):
            digest.update(chunk)
    assert 'sha256:' + digest.hexdigest() == expected, 'Official checksum mismatch'
    print('Verified:', target.name, flush=True)

    # Step3. 解压到独立工具目录，校验归档路径不会越界
    with zipfile.ZipFile(target) as archive:
        for entry in archive.infolist():
            assert (tool_dir / entry.filename).resolve().is_relative_to(tool_dir.resolve())
        archive.extractall(tool_dir)
    records.append({'repository': repository, 'release': asset['release'], 'asset': asset['name'], 'url': asset['browser_download_url'], 'sha256': digest.hexdigest(), 'bytes': asset['size']})
    (outdir / 'static_tool_versions.json').write_text(json.dumps(records, indent=2) + '\n', encoding='utf-8')
    print('Extracted:', repository, flush=True)
