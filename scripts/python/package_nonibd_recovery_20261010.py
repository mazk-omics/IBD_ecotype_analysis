# 打包本轮健康恢复的正式输出、审计表、注释、脚本和环境，不复制大型表达输入。
# 在服务器 ibd_refaudit 环境运行；同一 SSH 隧道回传 ZIP，另存校验值。
# 表达输入的实际路径和指纹保留在 input 的审计记录中，原始输入不动。
from pathlib import Path
import argparse
import hashlib
import zipfile
import json

# Step1. 确认计算和审计完成，明确需要回传的文件
parser = argparse.ArgumentParser()
parser.add_argument('--root', required=True)
args = parser.parse_args()
root = Path(args.root)
assert (root / 'run_exit_code.txt').read_text().strip() == '0'
assert (root / 'coassociation_exit_code.txt').read_text().strip() == '0'
assert json.loads((root / 'audit/integrity_check.json').read_text())['Status'] == 'PASS'
files = []
for folder in ['audit', 'formal', 'scripts', 'environment']:
    files.extend(p for p in (root / folder).rglob('*') if p.is_file())
files.extend(p for p in (root / 'input').glob('*') if p.suffix in ['.tsv', '.json'])
files.extend(p for p in root.glob('*') if p.is_file() and (p.suffix == '.log' or p.name.endswith('_code.txt') or p.name in ['run_started.txt', 'run_finished.txt']))
files = sorted(set(files))
print('Transfer files:', len(files), 'uncompressed bytes:', sum(p.stat().st_size for p in files))

# Step2. 保存逐文件指纹与压缩包指纹，通过隧道下载后核验
manifest = root / 'transfer_manifest.tsv'
with manifest.open('w', encoding='utf-8') as handle:
    handle.write('File\tBytes\tSHA256\n')
    for p in files:
        handle.write(f'{p.relative_to(root).as_posix()}\t{p.stat().st_size}\t{hashlib.sha256(p.read_bytes()).hexdigest()}\n')
archive = root / 'nonIBD_healthy_recovery_review.zip'
# 完整包已生成时保留该快照；再次执行只更新小包，避免重复压缩大图。
if not archive.exists():
    with zipfile.ZipFile(archive, 'w', compression=zipfile.ZIP_DEFLATED) as handle:
        for p in files + [manifest]:
            handle.write(p, p.relative_to(root).as_posix())
digest = hashlib.sha256(archive.read_bytes()).hexdigest()
(root / 'nonIBD_healthy_recovery_review.zip.sha256').write_text(digest + '\n', encoding='ascii')
print('Archive:', archive, 'bytes:', archive.stat().st_size, 'SHA256:', digest)

# Step3. 隧道较慢时先回传表和生态型图，状态热图保留服务器且逐文件列出
retained = [p for p in files if 'formal' in p.relative_to(root).parts and (p.suffix.lower() in ['.pdf', '.png'] or p.name == 'heatmap_data.txt')]
retained_manifest = root / 'server_retained_large_outputs.tsv'
with retained_manifest.open('w', encoding='utf-8') as handle:
    handle.write('ServerPath\tBytes\tSHA256\n')
    for p in retained:
        handle.write(f'{p}\t{p.stat().st_size}\t{hashlib.sha256(p.read_bytes()).hexdigest()}\n')
light_files = [p for p in files if p not in retained]
light_manifest = root / 'light_transfer_manifest.tsv'
with light_manifest.open('w', encoding='utf-8') as handle:
    handle.write('File\tBytes\tSHA256\n')
    for p in light_files + [retained_manifest]:
        handle.write(f'{p.relative_to(root).as_posix()}\t{p.stat().st_size}\t{hashlib.sha256(p.read_bytes()).hexdigest()}\n')
light_archive = root / 'nonIBD_healthy_recovery_tables.zip'
with zipfile.ZipFile(light_archive, 'w', compression=zipfile.ZIP_DEFLATED) as handle:
    for p in light_files + [retained_manifest, light_manifest]:
        handle.write(p, p.relative_to(root).as_posix())
light_digest = hashlib.sha256(light_archive.read_bytes()).hexdigest()
(root / 'nonIBD_healthy_recovery_tables.zip.sha256').write_text(light_digest + '\n', encoding='ascii')
print('Tables archive:', light_archive, 'bytes:', light_archive.stat().st_size, 'SHA256:', light_digest)
