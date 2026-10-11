# 固定官方EcoTyper提交，阅读公开接口脚本；只保存来源和接口线索，不执行第三方代码。
# 本机Python运行，读取本轮仓库清单，输出源码缓存与逐行定位，供后续设计正式输入接口。
from pathlib import Path
from urllib.request import Request, urlopen
import hashlib
import json

# Step1. 仅取官方公开的接口和presorted基因过滤脚本，使用已记录的提交
project = Path(__file__).resolve().parents[2]
audit = project / 'CIBERSORTx_server_plan_20261011'
tree = json.loads((audit / 'public_repository_audit/digitalcytometry_ecotyper_file_tree.json').read_text(encoding='utf-8'))
source_paths = ['pipeline/csx_hires_worker.R', 'pipeline/csx_hires_scheduler.R', 'pipeline/csx_hires_aggregate_worker_results.R', 'pipeline/state_discovery_presorted_filter_genes.R', 'pipeline/csx_fractions_S_mode.R']
outdir = audit / 'public_repository_audit/ecotyper_interface_sources'
outdir.mkdir(exist_ok=False)
records = []
for source_path in source_paths:
    assert source_path in tree['files']
    url = f"https://raw.githubusercontent.com/digitalcytometry/ecotyper/{tree['commit']}/{source_path}"
    with urlopen(Request(url, headers={'User-Agent': 'IBD-scientific-method-audit'}), timeout=20) as response:
        contents = response.read()
    (outdir / Path(source_path).name).write_bytes(contents)
    lines = contents.decode('utf-8').splitlines()
    evidence = [{'line_1based': i, 'text': line.strip()} for i, line in enumerate(lines, 1) if any(term in line.lower() for term in ['cibersortx', 'write.table', 'read.table', 'hires', 'na.omit', 'is.na', 'presort'])]
    records.append({'path': source_path, 'commit': tree['commit'], 'url': url, 'sha256': hashlib.sha256(contents).hexdigest(), 'lines': len(lines), 'interface_locations': evidence})
    print(source_path, len(lines), 'lines cached, not executed')
(audit / 'public_repository_audit/ecotyper_interface_evidence.json').write_text(json.dumps(records, ensure_ascii=False, indent=2), encoding='utf-8')
