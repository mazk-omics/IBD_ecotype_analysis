# 只读检查旧版小表的结构；在独立目录运行，不读取原始表达矩阵。
# 使用本机已有 Python，不安装软件，不连接服务器。
import os
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
from pathlib import Path
import gzip
import json
import pandas as pd

# Step1. 检查候选证据表
root = Path('D:/IBD_ecotype/IBD_ecotype_analysis')
out = Path(__file__).resolve().parents[1]
base = root / '06_validation/paper_aligned_recovery_characterization_20261012'
paths = list((base / '00_audit').glob('*.tsv')) + list((base / '02_external_markers').glob('*.tsv')) + list((base / '06_ecotype_features').glob('*.tsv'))
paths += list((root / '06_validation/all_states_lineage_audit_20261012/tables').glob('*.tsv'))
rows = []
for path in paths:
    table = pd.read_csv(path, sep='\t', nrows=2)
    rows.append({'Path': str(path), 'Bytes': path.stat().st_size, 'Columns': ';'.join(table.columns)})
pd.DataFrame(rows).to_csv(out / 'audit/source_schema_inventory.tsv', sep='\t', index=False)

# Step2. 只看压缩DE首行与状态审计关键列
with gzip.open(base / '02_external_markers/external_DEG_by_dataset.tsv.gz', 'rt', encoding='utf-8') as handle:
    print('DE_HEADER', handle.readline().strip())
    print('DE_ROW', handle.readline().strip())
review = pd.read_csv(root / '06_validation/all_states_lineage_audit_20261012/tables/all_53_states_final_evidence_review.tsv', sep='\t')
print(review[['CellState', 'CandidateProgram', 'Colon_MatchedDonors_Min5CellsEach', 'Ileum_MatchedDonors_Min5CellsEach', 'Colon_PairedStatus', 'Ileum_PairedStatus']].to_string(index=False))
print('Small-table inventory', len(rows))
