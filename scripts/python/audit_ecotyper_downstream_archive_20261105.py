# 核对上一轮 EcoTyper 下游代码归档的脚本哈希、Python 语法和方法入口。
# 在项目根目录或任意目录运行；可选第1参数为项目根目录，第2参数为JSON输出路径。
# 本脚本不读取患者级结果、不运行统计，也不修改历史分析文件。
from pathlib import Path
import ast
import csv
import hashlib
import json
import sys

project = Path(sys.argv[1]).resolve() if len(sys.argv) >= 2 else Path(__file__).resolve().parents[2]
output = Path(sys.argv[2]).resolve() if len(sys.argv) >= 3 else None

# Step1. 读取两套脚本清单
packages = [
    project / '06_validation/paper_aligned_recovery_characterization_20261012',
    project / '06_validation/paper_figures_20261012',
]
manifests = [
    packages[0] / 'scripts_manifest_repository_20261105.tsv',
    packages[1] / 'scripts_manifest_repository_20261105.tsv',
]

rows = []
for package, manifest in zip(packages, manifests):
    with manifest.open(encoding='utf-8', newline='') as handle:
        records = list(csv.DictReader(handle, delimiter='\t'))
    for record in records:
        path = package / record['File']
        actual = hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None
        rows.append({'package': package.name, 'file': record['File'], 'exists': path.exists(), 'sha256_match': actual == record['SHA256']})

# Step2. 对归档Python脚本做语法检查
python_rows = []
for package in packages:
    for path in sorted((package / 'scripts').rglob('*.py')):
        try:
            ast.parse(path.read_text(encoding='utf-8'))
            syntax_ok = True
            error = None
        except SyntaxError as exc:
            syntax_ok = False
            error = f'{exc.msg} at line {exc.lineno}'
        python_rows.append({'file': path.relative_to(project).as_posix(), 'syntax_ok': syntax_ok, 'error': error})

# Step3. 核对统一阅读入口
required_docs = [
    project / 'docs/methods/ECOTYPER_DOWNSTREAM_ANALYSIS.md',
    project / 'docs/methods/ECOTYPER_DOWNSTREAM_CODE_STYLE.md',
    packages[0] / 'README.md',
    packages[1] / 'README.md',
    packages[0] / '00_audit/scripts_manifest_final.tsv',
    packages[0] / 'scripts_manifest_repository_20261105.tsv',
    packages[1] / 'scripts_manifest_repository_20261105.tsv',
    packages[0] / '00_audit/method_parity_matrix.tsv',
    packages[1] / 'reports/figure_contract.json',
]

summary = {
    'manifest_files': len(rows),
    'manifest_pass': all(row['exists'] and row['sha256_match'] for row in rows),
    'python_files': len(python_rows),
    'python_syntax_pass': all(row['syntax_ok'] for row in python_rows),
    'required_docs': len(required_docs),
    'required_docs_pass': all(path.exists() for path in required_docs),
}
report = {'summary': summary, 'manifest': rows, 'python_syntax': python_rows, 'missing_docs': [path.relative_to(project).as_posix() for path in required_docs if not path.exists()]}

print(json.dumps(summary, ensure_ascii=False, indent=2))
if output is not None:
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')

assert summary['manifest_pass']
assert summary['python_syntax_pass']
assert summary['required_docs_pass']
