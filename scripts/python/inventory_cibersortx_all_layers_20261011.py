# 只读扫描 Fractions/GEP/HiRes 的所有官方镜像层，寻找应用源码和早期层残留。
# 本机 cibersortx_review 环境运行；保存完整层目录、应用文件版本和无认证内容的 R 辅助模块，不运行容器。
from pathlib import Path, PurePosixPath
import urllib.request
import tarfile
import hashlib
import json
import csv
import re
import shutil

# Step1. 使用已核对的镜像 manifest，公共下载凭据只保留在内存
project = Path(__file__).resolve().parents[2]
metadata = project / 'CIBERSORTx_server_plan_20261009'
previous = project / 'CIBERSORTx_server_plan_20261010'
outdir = project / 'CIBERSORTx_server_plan_20261011'
cache = project.parent / '.tools' / 'cibersortx_image_layers'
reference = outdir / 'vendor_reference'
cache.mkdir(parents=True, exist_ok=True)
reference.mkdir(exist_ok=True)
exclude = project / '.git' / 'info' / 'exclude'
rules = exclude.read_text(encoding='utf-8') if exclude.exists() else ''
for rule in ['/CIBERSORTx_server_plan_20261011/vendor_reference/', '/CIBERSORTx_server_plan_20261011/decompiler_workspace/', '/CIBERSORTx_server_plan_20261011/decompiled_numeric/']:
    if rule not in rules.splitlines():
        rules += '\n' + rule + '\n'
exclude.write_text(rules, encoding='utf-8')  # 供应商文件留在本机，不加入研究代码 Git

layers_checked = []
applications = []
summaries = []
for name in ['hires', 'fractions', 'gep']:  # 三个镜像及共享层需分别核对最终路径
    manifest = json.loads((metadata / (name + '_registry_manifest.json')).read_text())
    config = json.loads((metadata / (name + '_image_config.json')).read_text())
    history = [step for step in config['history'] if not step.get('empty_layer', False)]
    assert len(history) == len(manifest['layers'])
    auth_url = 'https://auth.docker.io/token?service=registry.docker.io&scope=repository:cibersortx/' + name + ':pull'
    with urllib.request.urlopen(auth_url, timeout=30) as response:
        registry_token = json.load(response)['token']
    final_files = {}
    application_versions = 0
    members_count = 0
    filename_sources = 0
    project_cpp = 0

    # Step2. 每层校验 SHA256 并逐条读取目录，避免解压整个系统或执行任何文件
    inventory_path = outdir / (name + '_all_layer_files.tsv')
    with inventory_path.open('w', encoding='utf-8', newline='') as inventory_stream:
        writer = csv.writer(inventory_stream, delimiter='\t')
        writer.writerow(['LayerIndex', 'LayerSHA256', 'Path', 'Bytes', 'MemberType'])
        for index, (layer, step) in enumerate(zip(manifest['layers'], history)):
            digest = layer['digest'].split(':')[1]
            target = cache / ('sha256_' + digest + '.tar.gz')
            old = previous / target.name
            if not target.exists() and old.exists():
                shutil.copyfile(old, target)
            if not target.exists():
                print('Downloading layer:', name, index, round(layer['size'] / 1024 / 1024, 2), 'MB', flush=True)
                request = urllib.request.Request('https://registry-1.docker.io/v2/cibersortx/' + name + '/blobs/' + layer['digest'], headers={'Authorization': 'Bearer ' + registry_token})
                partial = target.with_suffix('.download')
                total = 0
                with urllib.request.urlopen(request, timeout=60) as response, partial.open('wb') as stream:
                    while chunk := response.read(1024 * 1024):
                        stream.write(chunk)
                        total += len(chunk)
                assert total == layer['size'], 'Download size mismatch'
                partial.rename(target)
            checksum = hashlib.sha256()
            with target.open('rb') as stream:
                while chunk := stream.read(1024 * 1024):
                    checksum.update(chunk)
            assert checksum.hexdigest() == digest, 'Layer hash mismatch'
            layer_applications = 0
            with tarfile.open(target, 'r|gz') as archive:
                for member in archive:
                    path = member.name.removeprefix('./').lstrip('/')
                    writer.writerow([index, digest, path, member.size, member.type.decode('ascii', errors='replace')])
                    members_count += 1
                    if member.isfile() and re.search(r'(?i)\.(R|cpp|cc|cxx|c|h|hpp|py|sh)$', path):
                        filename_sources += 1
                    if not path.startswith('src/'):
                        continue
                    parts = PurePosixPath(path)
                    if parts.name.startswith('.wh.'):
                        removed = str(parts.parent / parts.name.removeprefix('.wh.'))
                        if parts.name == '.wh..wh..opq':
                            final_files = {p: row for p, row in final_files.items() if not p.startswith(str(parts.parent) + '/')}
                        else:
                            final_files.pop(removed, None)
                        continue
                    if not member.isfile():
                        continue
                    data = archive.extractfile(member).read()
                    row = {'Image': name, 'LayerIndex': index, 'LayerSHA256': digest, 'Path': path, 'Bytes': len(data), 'FileSHA256': hashlib.sha256(data).hexdigest(), 'Kind': 'ELF' if data.startswith(b'\x7fELF') else 'R' if path.lower().endswith('.r') else 'Other'}
                    if re.search(r'(?i)\.(cpp|cc|cxx|c|h|hpp)$', path):
                        project_cpp += 1
                    applications.append(row)
                    final_files[path] = row
                    application_versions += 1
                    layer_applications += 1

                    # Step3. 仅导出无认证/联网标记的 R 辅助文件；主 GEP 单独进行语法级筛选
                    destination = reference / name / path
                    if row['Kind'] == 'ELF' or row['Kind'] == 'R':
                        destination.parent.mkdir(parents=True, exist_ok=True)
                    if row['Kind'] == 'ELF':
                        destination.write_bytes(data)  # 仅作为反编译输入，绝不执行
                    elif row['Kind'] == 'R':
                        source_text = data.decode('utf-8-sig')
                        sensitive = bool(re.search(r'(?i)token|username|auth|password|credential|https?://|\b(?:POST|GET)\s*\(|[\w.+-]+@[\w.-]+\.[A-Za-z]{2,}', source_text))
                        row['RExportStatus'] = 'withheld_for_expression_level_review' if sensitive else 'readable_helper_exported'
                        if not sensitive:
                            destination.write_bytes(data)
            layers_checked.append({'Image': name, 'Index': index, 'LayerSHA256': digest, 'CompressedBytes': layer['size'], 'ProjectFileVersions': layer_applications})
            print('Checked layer:', name, index, 'application files:', layer_applications, flush=True)

    # Step4. 每完成一个镜像即保存检查点；系统/第三方源码数量与应用源码分开说明
    summaries.append({'Image': name, 'LayersScanned': len(manifest['layers']), 'ArchiveMembers': members_count, 'SourceFilenameEntriesAllDirectories': filename_sources, 'ApplicationFileVersions': application_versions, 'ApplicationCppHeaderVersions': project_cpp, 'FinalApplicationFiles': len(final_files), 'FinalPaths': sorted(final_files)})
    (outdir / 'full_image_scan_summary.json').write_text(json.dumps(summaries, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    (outdir / 'full_image_application_versions.json').write_text(json.dumps(applications, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    (outdir / 'full_image_layers_checked.json').write_text(json.dumps(layers_checked, indent=2) + '\n', encoding='utf-8')
    print('Completed image:', name, 'files:', len(final_files), 'C++ source/header versions:', project_cpp, flush=True)
print('Output:', outdir, flush=True)
