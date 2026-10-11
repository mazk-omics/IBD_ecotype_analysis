# 将镜像中的 R 文本交给 R 的 parse() 做语法级审计，不执行供应商程序。
# 本机 cibersortx_review 环境运行；临时原文只写入已被 Git 排除的 vendor_reference 目录。
from pathlib import Path
import json
import tarfile
import hashlib
import subprocess
import re

# Step1. 从已有完整应用文件清单找到所有 R 文件版本，按镜像保留来源
project = Path(__file__).resolve().parents[2]
outdir = project / 'CIBERSORTx_server_plan_20261011'
cache = project.parent / '.tools' / 'cibersortx_image_layers'
reference = outdir / 'vendor_reference'
versions = json.loads((outdir / 'full_image_application_versions.json').read_text(encoding='utf-8'))
inputs = []
for row in versions:
    if row['Kind'] != 'R':
        continue
    with tarfile.open(cache / ('sha256_' + row['LayerSHA256'] + '.tar.gz'), 'r:gz') as archive:
        data = archive.extractfile(row['Path']).read()
    assert hashlib.sha256(data).hexdigest() == row['FileSHA256']
    inputs.append({'image': row['Image'], 'path': row['Path'], 'sha256': row['FileSHA256'], 'text': data.decode('utf-8-sig')})
temporary = reference / 'parse_inputs.json'
temporary.write_text(json.dumps(inputs, ensure_ascii=False), encoding='utf-8')

# Step2. R 只解析函数定义并导出通过筛选的片段；结束后删除这一个明确的临时文件
script = project / 'scripts' / 'R' / 'export_cibersortx_R_functions_20261011.R'
destination = reference / 'R_function_fragments'
result = subprocess.run([r'D:\R\R-4.4.1\bin\Rscript.exe', str(script), str(temporary), str(destination)], text=True, encoding='utf-8', errors='replace', capture_output=True)
temporary.unlink()
print(result.stdout)
if result.returncode:
    errors = [line for line in result.stderr.splitlines() if line.startswith('Error') or line.startswith('错误')]
    for line in errors:
        safe = re.sub(r'([\"\x27]).*?\1', '[literal withheld]', line)
        if not re.search(r'(?i)token|password|credential|https?://|@', safe):
            print(safe[:240])
    print('R parsing failed; inspect the local R script, without dumping raw supplier source.')
    raise SystemExit(result.returncode)
print('Source-level function inventory:', outdir / 'R_function_expression_audit.json')
