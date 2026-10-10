# 生成本轮 Git 代码归档清单，检查语法、文件大小和明显凭据，不读取或提交研究矩阵。
# 本地 Python 或服务器 python3 运行；参数 --root 指仓库，--label 区分服务器与本地来源。
from pathlib import Path
import argparse
import subprocess
import ast
import hashlib
import re
import csv
import shutil
import os

# Step1. 仅选择现有 Git 规则允许的代码、配置和整理后的方法记录
parser = argparse.ArgumentParser()
parser.add_argument('--root', required=True)
parser.add_argument('--label', required=True)
args = parser.parse_args()
root = Path(args.root).resolve()
untracked = subprocess.check_output(['git', 'ls-files', '--others', '--exclude-standard', '-z'], cwd=root).decode('utf-8').split('\0')
modified = subprocess.check_output(['git', 'diff', '--name-only', '-z'], cwd=root).decode('utf-8').split('\0')
staged = subprocess.check_output(['git', 'diff', '--cached', '--name-only', '-z'], cwd=root).decode('utf-8').split('\0')
selected = []
for name in sorted(set(untracked + modified + staged)):
    path = root / name
    if not name or not path.is_file() or name == 'docs/provenance/code_inventory_' + args.label + '_20261010.tsv':
        continue
    code = (name.startswith('scripts/') or (name.startswith('06_validation/') and '/scripts/' in name and '/server_snapshot/' not in name)) and path.suffix.lower() in {'.py', '.r', '.sh', '.ps1', '.md', '.yml', '.yaml'}
    configuration = name.startswith('05_EcoTyper/config/') and path.suffix.lower() in {'.yml', '.yaml'}
    public_docs = name.startswith('docs/') and path.suffix.lower() in {'.md', '.patch', '.json', '.tsv', '.txt', '.r'}
    environment = name.startswith('environments/') and path.suffix.lower() in {'.yml', '.yaml', '.txt', '.json', '.tsv'}
    policy = name in {'.gitignore', '.gitattributes', 'AGENTS.md', 'Profile/CODE_STYLE.md', 'Profile/SERVER_ACCESS_POLICY.md', 'README.md'}
    if code or configuration or public_docs or environment or policy:
        selected.append(name)

# Step2. 防止大数据和真实凭据混入；只输出文件名与错误类型，不打印匹配内容
issues = []
rows = []
bash = shutil.which('bash') if os.name == 'posix' else str(Path(shutil.which('git')).resolve().parents[1] / 'bin' / 'bash.exe')
secret_pattern = re.compile(r'(?:gh[pousr]_[A-Za-z0-9]{20,}|github_pat_[A-Za-z0-9_]{30,}|-----BEGIN (?:OPENSSH |RSA |EC )?PRIVATE KEY-----|(?:token|password|secret)\s*[:=]\s*[\"\x27][A-Za-z0-9+/=_-]{24,}[\"\x27])', re.I)
for name in selected:
    path = root / name
    data = path.read_bytes()
    text = data.decode('utf-8-sig')
    if len(data) > 2 * 1024 * 1024:
        issues.append((name, 'over_2_MiB'))
    if secret_pattern.search(text):
        issues.append((name, 'credential_candidate'))
    if path.suffix.lower() == '.py':
        try:
            ast.parse(text, filename=name)
        except SyntaxError as error:
            issues.append((name, 'python_syntax_line_' + str(error.lineno)))
    if path.suffix.lower() == '.sh' and bash and Path(bash).exists():
        shell_check = subprocess.run([bash, '-n', str(path)], capture_output=True)
        if shell_check.returncode != 0:
            issues.append((name, 'bash_syntax'))
    rows.append([name, len(data), hashlib.sha256(data).hexdigest()])
print('Selected files:', len(rows), 'bytes:', sum(row[1] for row in rows))
print('Issues:', issues)
if issues:
    raise SystemExit('Resolve publication checks before staging')

# Step3. 清单记录来源和 SHA256；路径列表用于显式 git add，不使用 git add .
provenance = root / 'docs' / 'provenance'
provenance.mkdir(parents=True, exist_ok=True)
manifest = provenance / ('code_inventory_' + args.label + '_20261010.tsv')
with manifest.open('w', encoding='utf-8', newline='') as stream:
    writer = csv.writer(stream, delimiter='\t')
    writer.writerow(['RepositoryPath', 'Bytes', 'SHA256'])
    writer.writerows(rows)
stage_paths = selected + [manifest.relative_to(root).as_posix()]
temporary = root / 'tmp'
temporary.mkdir(exist_ok=True)
(temporary / ('git_closeout_' + args.label + '.paths')).write_bytes(('\0'.join(stage_paths) + '\0').encode('utf-8'))
print('Manifest:', manifest)
print('Pathspec:', temporary / ('git_closeout_' + args.label + '.paths'))
