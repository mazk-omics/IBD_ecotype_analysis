# 在本地同步服务器提交前保存同名未跟踪代码；不合并分析逻辑，不移动数据矩阵。
# 本地 Python 运行，先 git fetch origin；输入 origin/main 新增路径，输出差异归档和来源清单。
from pathlib import Path
import subprocess
import hashlib
import csv
import shutil

# Step1. 核对服务器新增文件与本地未跟踪代码的重名情况
root = Path(__file__).resolve().parents[2]
incoming = subprocess.check_output(['git', 'diff', '--name-only', '-z', 'HEAD', 'origin/main'], cwd=root).decode('utf-8').split('\0')
untracked = subprocess.check_output(['git', 'ls-files', '--others', '--exclude-standard', '-z'], cwd=root).decode('utf-8').split('\0')
collisions = sorted((set(incoming) & set(untracked)) - {''})
rows = []

# Step2. 不同内容另存历史版本；相同内容保存在忽略的临时备份中，以允许 fast-forward
for name in collisions:
    source = (root / name).resolve()
    assert source.is_relative_to(root) and source.is_file(), name
    local = source.read_bytes()
    remote = subprocess.check_output(['git', 'show', 'origin/main:' + name], cwd=root)
    same = local.replace(b'\r\n', b'\n') == remote.replace(b'\r\n', b'\n')
    if same:
        destination = root / 'tmp' / 'pre_server_sync_20261010' / name
        decision = 'identical_except_line_endings; server version retained'
    else:
        destination = root / 'scripts' / 'archived' / '20261010_local_pre_server_sync' / name
        decision = 'different; local source archived without merging'
    destination = destination.resolve()
    assert destination.is_relative_to(root) and not destination.exists(), str(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.move(str(source), str(destination))
    rows.append([name, hashlib.sha256(local).hexdigest(), hashlib.sha256(remote).hexdigest(), decision, destination.relative_to(root).as_posix()])

# Step3. 保存可核对的来源记录，再单独执行 git merge --ff-only origin/main
out = root / 'docs' / 'provenance' / 'local_server_reconciliation_20261010.tsv'
out.parent.mkdir(parents=True, exist_ok=True)
with out.open('w', encoding='utf-8', newline='') as stream:
    writer = csv.writer(stream, delimiter='\t')
    writer.writerow(['OriginalPath', 'LocalSHA256', 'ServerSHA256', 'Decision', 'PreservedLocalPath'])
    writer.writerows(rows)
print('Overlapping untracked files:', len(rows))
print('Distinct local versions archived:', sum('different;' in row[3] for row in rows))
print('Source record:', out)
