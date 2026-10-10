# 汇集研究目录中尚未集中保存的项目脚本；相同内容只记录来源，差异版本另存历史归档。
# 本地 Python 运行；读取研究目录中的代码文件，输出 scripts/archived 和 docs/provenance，不读取矩阵。
from pathlib import Path
import hashlib
import csv
import shutil
import subprocess

# Step1. 先记录已集中保存的源码指纹，以免重复归档下载包中的同一份脚本
root = Path(__file__).resolve().parents[2]
extensions = {'.py', '.r', '.sh', '.ps1'}
known = {}
tracked = subprocess.check_output(['git', 'ls-files', '-z'], cwd=root).decode('utf-8').split('\0')
existing = list((root / 'scripts').rglob('*')) + [root / name for name in tracked if name]
for path in existing:
    if path.is_file() and path.suffix.lower() in extensions:
        fingerprint = hashlib.sha256(path.read_bytes().replace(b'\r\n', b'\n')).hexdigest()
        known.setdefault(fingerprint, path.relative_to(root).as_posix())

# Step2. 只扫描项目准备/审计目录；不收集第三方官方源码快照或镜像内的程序
folders = ['06_validation', 'IBDverse_reference', 'CIBERSORTx_ready_20261007', 'CIBERSORTx_ready_20261008', 'CIBERSORTx_signature_audit']
rows = []
for folder in folders:
    for path in sorted((root / folder).rglob('*')):
        if not path.is_file() or path.suffix.lower() not in extensions or 'software' in path.parts or '__pycache__' in path.parts:
            continue
        name = path.relative_to(root).as_posix()
        data = path.read_bytes()
        fingerprint = hashlib.sha256(data.replace(b'\r\n', b'\n')).hexdigest()
        if fingerprint in known:
            destination = known[fingerprint]
            decision = 'identical normalized source; existing version retained'
        else:
            target = root / 'scripts' / 'archived' / '20261010_research_folder_sources' / name
            target = target.resolve()
            assert target.is_relative_to(root) and not target.exists(), str(target)
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(path, target)
            destination = target.relative_to(root).as_posix()
            known[fingerprint] = destination
            decision = 'distinct source copied; original retained'
        rows.append([name, hashlib.sha256(data).hexdigest(), decision, destination])

# Step3. 原目录不改动，保存每个脚本的来源与最终代码入口
out = root / 'docs' / 'provenance' / 'research_folder_source_inventory_20261010.tsv'
out.parent.mkdir(parents=True, exist_ok=True)
with out.open('w', encoding='utf-8', newline='') as stream:
    writer = csv.writer(stream, delimiter='\t')
    writer.writerow(['SourcePath', 'SourceSHA256', 'Decision', 'CodePath'])
    writer.writerows(rows)
print('Source snapshots:', len(rows))
print('Distinct copies:', sum('distinct source' in row[2] for row in rows))
print('Source record:', out)
