# 核对服务器研究目录里剩余的未跟踪代码，保存不同版本；原始脚本和研究数据不改动。
# 在项目根目录运行 python3 scripts/python/check_git_source_coverage_20261010.py；输出差异源码与来源清单。
from pathlib import Path
import subprocess
import hashlib
import shutil
import csv

# Step1. 以已提交源码的内容核对遗漏，临时 HTTP 传输服务不属于分析归档
root = Path(__file__).resolve().parents[2]
extensions = {'.py', '.r', '.sh', '.ps1'}
tracked = subprocess.check_output(['git', 'ls-files', '-z'], cwd=root).decode('utf-8').split('\0')
known = {}
for name in tracked:
    path = root / name
    if path.is_file() and path.suffix.lower() in extensions:
        fingerprint = hashlib.sha256(path.read_bytes().replace(b'\r\n', b'\n')).hexdigest()
        known.setdefault(fingerprint, name)
untracked = subprocess.check_output(['git', 'ls-files', '--others', '--exclude-standard', '-z'], cwd=root).decode('utf-8').split('\0')
rows = []

# Step2. 同内容只记录归档入口，不同内容复制到服务器历史版本目录
for name in sorted(untracked):
    path = root / name
    if not name or path.suffix.lower() not in extensions or not name.startswith(('03_reference/', '06_validation/')):
        continue
    if path.name == 'transfer_http.py':
        rows.append([name, '', 'temporary transfer utility excluded', ''])
        continue
    data = path.read_bytes()
    fingerprint = hashlib.sha256(data.replace(b'\r\n', b'\n')).hexdigest()
    if fingerprint in known:
        destination = known[fingerprint]
        decision = 'identical normalized source; already archived'
    else:
        target = (root / 'scripts' / 'archived' / '20261010_server_research_sources' / name).resolve()
        assert target.is_relative_to(root) and not target.exists(), str(target)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, target)
        destination = target.relative_to(root).as_posix()
        known[fingerprint] = destination
        decision = 'distinct server source copied; original retained'
    rows.append([name, hashlib.sha256(data).hexdigest(), decision, destination])

# Step3. 保留来源与版本差异，文件原位保留；后续按显式清单暂存并检查
out = root / 'docs' / 'provenance' / 'server_research_source_coverage_20261010.tsv'
with out.open('w', encoding='utf-8', newline='') as stream:
    writer = csv.writer(stream, delimiter='\t')
    writer.writerow(['SourcePath', 'SourceSHA256', 'Decision', 'CodePath'])
    writer.writerows(rows)
print('Remaining source candidates:', len(rows))
print('Distinct server copies:', sum('distinct server' in row[2] for row in rows))
print('Coverage record:', out)
