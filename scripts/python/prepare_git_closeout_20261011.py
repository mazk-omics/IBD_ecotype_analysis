# 归档2026-10-11自写源码、方法报告和环境；供应商代码、矩阵及其他日期文件不纳入。
# 本机科学Python环境运行，生成明确暂存清单，检查语法/大小/凭据；不自动commit或push。
from pathlib import Path
import ast
import csv
import hashlib
import json
import re
import subprocess

# Step1. 明确今天的文件范围，避免git add .带入未整理研究结果
project = Path(__file__).resolve().parents[2]
git = r'D:\Git\cmd\git.exe'
bash = r'D:\Git\bin\bash.exe'
rscript = r'D:\R\R-4.4.1\bin\Rscript.exe'
selected = sorted({p.relative_to(project).as_posix() for p in (project / 'scripts').rglob('*20261011*') if p.is_file() and p.suffix.lower() in {'.py', '.r', '.sh', '.ps1'}})
selected += sorted(p.relative_to(project).as_posix() for p in (project / 'scripts/java').glob('*.java'))
selected += sorted(p.relative_to(project).as_posix() for p in (project / 'CIBERSORTx_server_plan_20261011').glob('*.md'))
selected += sorted(p.relative_to(project).as_posix() for p in (project / 'CIBERSORTx_server_plan_20261011').glob('*_stage_manifest.json'))
selected += sorted(p.relative_to(project).as_posix() for p in (project / 'environments/local_research_20261011').glob('*') if p.is_file())
selected += ['.gitignore', '.gitattributes', 'docs/audits/cibersortx_research_closeout_20261011.md']
selected = sorted(set(selected))
existing_staged = subprocess.check_output([git, 'diff', '--cached', '--name-only', '-z'], cwd=project).decode('utf-8').split('\0')
assert set(filter(None, existing_staged)) <= set(selected + ['docs/provenance/code_inventory_local_20261011.tsv', 'docs/provenance/local_archive_checks_20261011.json']), 'Existing staged work outside this closeout must be kept separate'
temporary = project / 'tmp'
temporary.mkdir(exist_ok=True)

# Step2. 只输出凭据候选的文件名，不输出匹配内容；文本换行归一化前后分别记录SHA
secret_pattern = re.compile(r'(?:gh[pousr]_[A-Za-z0-9]{20,}|github_pat_[A-Za-z0-9_]{30,}|-----BEGIN (?:OPENSSH |RSA |EC )?PRIVATE KEY-----|(?:token|password|secret)\s*[:=]\s*[\"\x27][A-Za-z0-9+/=_-]{24,}[\"\x27])', re.I)
rows = []
issues = []
python_count = bash_count = 0
for name in selected:
    path = project / name
    data = path.read_bytes()
    text = data.decode('utf-8-sig')
    if len(data) > 2 * 1024 * 1024:
        issues.append([name, 'over_2_MiB'])
    if secret_pattern.search(text):
        issues.append([name, 'credential_candidate'])
    if path.suffix.lower() == '.py':
        ast.parse(text, filename=name)
        python_count += 1
    if path.suffix.lower() == '.sh':
        subprocess.run([bash, '-n', str(path)], check=True, capture_output=True)
        bash_count += 1
    if path.suffix.lower() == '.json':
        json.loads(text)
    blob = subprocess.check_output([git, 'hash-object', '--path=' + name, str(path)], cwd=project).decode().strip()
    rows.append([name, len(data), hashlib.sha256(data).hexdigest(), blob])
assert not issues, issues

# Step3. R只parse，不重跑计算；PowerShell文件清单交给原生AST检查
r_files = [str(project / name) for name in selected if Path(name).suffix.lower() == '.r']
r_expression = '\n'.join(["# Git归档临时语法检查，不执行被审阅脚本。", "# Step1. 逐个解析显式文件列表", "args <- commandArgs(trailingOnly = TRUE)", "for (path in args) parse(file = path, encoding = 'UTF-8')", "cat('R parse passed:', length(args), '\\n')"])
r_parse_script = temporary / 'git_closeout_20261011_parse.R'
r_parse_script.write_text(r_expression + '\n', encoding='utf-8')
r_check = subprocess.run([rscript, str(r_parse_script), *r_files], text=True, capture_output=True, encoding='utf-8', errors='replace')
assert r_check.returncode == 0 and f'R parse passed: {len(r_files)}' in r_check.stdout, 'R parse failed or did not confirm file count; inspect stderr locally'
print(r_check.stdout.strip())
ps_files = [name for name in selected if Path(name).suffix.lower() == '.ps1']
(temporary / 'git_closeout_20261011_powershell.paths').write_text('\n'.join(ps_files) + '\n', encoding='utf-8')

# Step4. 生成清单与方法摘要的Git入口，私有证据只保留哈希和所在位置
provenance = project / 'docs/provenance'
provenance.mkdir(exist_ok=True)
manifest = provenance / 'code_inventory_local_20261011.tsv'
with manifest.open('w', encoding='utf-8', newline='') as stream:
    writer = csv.writer(stream, delimiter='\t')
    writer.writerow(['RepositoryPath', 'WorkingTreeBytes', 'WorkingTreeSHA256', 'GitBlobSHA1_after_attributes'])
    writer.writerows(rows)
checks = {'date': '2026-10-11', 'selected_files_before_closeout_records': len(rows), 'total_source_bytes': sum(row[1] for row in rows), 'python_AST_files': python_count, 'Bash_syntax_files': bash_count, 'R_parse_files': len(r_files), 'PowerShell_AST': 'pending native check', 'Java_Ghidra_scripts': 'Archive own scripts; static analyses already recorded, no compilation rerun during closeout', 'publication_scan': 'No credential candidates; selected files below 2 MiB', 'numerical_validation': {'HiRes': 'full_target_stage_manifest.json', 'Group': 'group_reconstruction_stage_manifest.json'}, 'runtime_results': 'Local Windows only; no SSH or Linux integrated rerun', 'excluded': 'Vendor files, decompiled pseudocode, paper full text, original and derived matrices, sample-level files, unrelated untracked research'}
check_path = provenance / 'local_archive_checks_20261011.json'
check_path.write_text(json.dumps(checks, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
stage_paths = selected + [manifest.relative_to(project).as_posix(), check_path.relative_to(project).as_posix()]
(temporary / 'git_closeout_20261011.paths').write_bytes(('\0'.join(stage_paths) + '\0').encode('utf-8'))
print(json.dumps(checks, ensure_ascii=False, indent=2))
