# 在本机独立的 SeqOut Python 3.13 环境中搜索 recovery 候选，保存检索记录和原始元数据。
# 输入是下方明确列出的查询与 accession；输出在 06_validation/recovery_candidates_20261007。
# 只检索和下载小型元数据，不下载表达矩阵、FASTQ、模型权重或整个数据库。
from pathlib import Path
from importlib.metadata import version
import json
import os
import subprocess
import sys

# Step1. 确认独立工具路径，记录版本以便复查检索结果。
project = Path(__file__).resolve().parents[3]
outdir = project / 'IBD_ecotype_analysis/06_validation/recovery_candidates_20261007'
outdir.mkdir(parents=True, exist_ok=True)
seqout = project / '.tools/seqout/environment/Scripts/seqout.exe'
environment = os.environ.copy()
environment['PYTHONUTF8'] = '1'
record = {'python': sys.version, 'python_executable': sys.executable, 'seqout_version': version('seqout'), 'source_commit': 'aea316e4cc9a778995bb324636d682f37129a54a', 'queries': [], 'metadata': []}
(outdir / 'seqout_execution_record.json').write_text(json.dumps(record, ensure_ascii=False, indent=2), encoding='utf-8')
print('SeqOut:', record['seqout_version'], 'Python:', sys.version.split()[0], flush=True)

# Step2. 分开检索 CD、UC 空间数据和单细胞数据；结果只是候选，需要回到原始档案确认。
queries = [('cd_spatial', 'Crohn spatial transcriptomics', 'geo'), ('uc_visium', 'ulcerative colitis Visium', 'geo'), ('ibdverse', 'IBDverse', 'arrayexpress'), ('uc_scrna', 'ulcerative colitis single cell', 'geo')]
for name, query, database in queries:
    command = [str(seqout), 'search', query, '--organism', 'Homo sapiens', '--db', database, '--max', '30', '-o', str(outdir / (name + '_search.json'))]
    print('查询:', name, query, flush=True)
    result = subprocess.run(command, capture_output=True, text=True, encoding='utf-8', env=environment, timeout=180)
    (outdir / (name + '_search_console.txt')).write_text(result.stdout + result.stderr, encoding='utf-8')
    record['queries'].append({'name': name, 'query': query, 'database': database, 'command': command, 'returncode': result.returncode, 'output_exists': (outdir / (name + '_search.json')).exists()})
    (outdir / 'seqout_execution_record.json').write_text(json.dumps(record, ensure_ascii=False, indent=2), encoding='utf-8')
    print('返回码:', result.returncode, flush=True)

# Step3. 对重点 GEO 候选保存原始元数据；不启用 supplementary/FASTQ 下载选项。
for accession in ['GSE312415', 'GSE189184', 'GSE277348', 'GSE298464']:
    command = [str(seqout), 'download', accession, '-o', str(outdir / (accession + '_seqout_metadata.json'))]
    print('元数据:', accession, flush=True)
    result = subprocess.run(command, capture_output=True, text=True, encoding='utf-8', env=environment, timeout=180)
    (outdir / (accession + '_metadata_console.txt')).write_text(result.stdout + result.stderr, encoding='utf-8')
    record['metadata'].append({'accession': accession, 'command': command, 'returncode': result.returncode, 'output_exists': (outdir / (accession + '_seqout_metadata.json')).exists()})
    (outdir / 'seqout_execution_record.json').write_text(json.dumps(record, ensure_ascii=False, indent=2), encoding='utf-8')
    print('返回码:', result.returncode, flush=True)
