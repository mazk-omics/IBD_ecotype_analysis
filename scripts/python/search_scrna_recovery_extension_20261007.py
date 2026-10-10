# 用已安装的 SeqOut 补搜肠道单细胞 recovery 候选；在本机独立 Python 3.13 环境运行。
# 输入为下方疾病/部位查询，输出到 recovery_candidates_20261007/scRNA_extension。
# 只保存原始检索和标题索引；候选的细胞覆盖、文件可用性与独立性另行核对。
from pathlib import Path
from importlib.metadata import version
import json
import os
import subprocess
import pandas as pd

# Step1. 分开查询不同疾病、部位与实验设计，避免一次搜索被某类研究占满。
project = Path(__file__).resolve().parents[3]
outdir = project / 'IBD_ecotype_analysis/06_validation/recovery_candidates_20261007/scRNA_extension'
outdir.mkdir(parents=True, exist_ok=True)
seqout = project / '.tools/seqout/environment/Scripts/seqout.exe'
environment = os.environ.copy()
environment['PYTHONUTF8'] = '1'
queries = [('cd_ileum', 'Crohn single cell ileum'), ('uc_colon', 'ulcerative colitis single cell colon'), ('ibd_scrna', 'inflammatory bowel single cell'), ('cd_fibrosis', 'Crohn fibrosis single cell'), ('uc_treatment', 'ulcerative colitis remission single cell'), ('healthy_gut', 'healthy intestine single cell'), ('ibd_nuclei', 'inflammatory bowel single nucleus')]
record = {'seqout_version': version('seqout'), 'queries': []}
rows = []

# Step2. 逐条执行真实查询并保存响应；最多 60 条命中是本次检索范围，不代表完整文献集合。
for name, query in queries:
    path = outdir / (name + '_search.json')
    command = [str(seqout), 'search', query, '--organism', 'Homo sapiens', '--db', 'geo', '--max', '60', '-o', str(path)]
    print('查询:', name, query, flush=True)
    result = subprocess.run(command, capture_output=True, text=True, encoding='utf-8', env=environment, timeout=180)
    (outdir / (name + '_console.txt')).write_text(result.stdout + result.stderr, encoding='utf-8')
    hits = json.loads(path.read_text(encoding='utf-8')) if result.returncode == 0 and path.exists() else []
    record['queries'].append({'name': name, 'query': query, 'command': command, 'returncode': result.returncode, 'hits': len(hits)})
    (outdir / 'execution_record.json').write_text(json.dumps(record, ensure_ascii=False, indent=2), encoding='utf-8')
    for hit in hits:
        rows.append({'query': name, 'accession': hit.get('accession'), 'title': hit.get('title'), 'summary': hit.get('summary'), 'doi': hit.get('doi'), 'is_single_cell': hit.get('is_single_cell'), 'n_samples': hit.get('n_samples')})
    print('返回码:', result.returncode, '命中:', len(hits), flush=True)

# Step3. 合并索引便于逐个审阅；同一 accession 的多次命中保留在原表。
index = pd.DataFrame(rows)
index.to_csv(outdir / '所有查询命中.tsv', sep='\t', index=False)
unique = index.drop_duplicates('accession').copy()
unique.to_csv(outdir / '不重复研究标题索引.tsv', sep='\t', index=False)
print('查询命中:', len(index), '不重复 accession:', len(unique), flush=True)
print(unique[['accession', 'title']].to_string(index=False), flush=True)
