# 完成解释和绘图后检查回传原件、状态计数、证据来源与图件；不改变任何结果。
# 本机cibersortx_review环境运行，检查摘要写入独立解析目录，保留未分配及缺失Glial。
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from pypdf import PdfReader

# Step1. 重新核对回传清单，确认本地审阅没有改写服务器原始输出
project = Path(__file__).resolve().parents[3] / 'IBD_ecotype_analysis'
root = project / '06_validation/GSE214695_scRNA_recovery_20261009'
raw = root / 'raw_server'
interpretation = project / '06_validation/recovery_interpretation_20261010'
manifest = json.loads((raw / 'download_manifest.json').read_text())
for item in manifest:
    path = raw / item['path']
    assert path.stat().st_size == item['bytes']
    assert hashlib.sha256(path.read_bytes()).hexdigest() == item['sha256']
print('回传原件哈希通过：', len(manifest))

# Step2. 状态与细胞总量一致；两轮不同队列不能伪作配对比较
cells = pd.read_csv(raw / 'audit/cell_assignments_all.tsv', sep='\t')
states = pd.read_csv(root / '53状态单细胞恢复审阅.tsv', sep='\t')
comparison = pd.read_csv(root / '53状态_结肠回肠描述性比较.tsv', sep='\t')
summary = json.loads((root / 'review_summary.json').read_text())
assert len(cells) == summary['input_cells'] == 9752 and cells.ID.is_unique
assert len(states) == len(comparison) == 53 and states.CellState.is_unique
assert states.AssignedCells.sum() == summary['assigned_cells'] == cells.State.ne('Unassigned').sum()
assert states.AboveNominalZ1_96.sum() == summary['states_above_nominal_z_1_96']
assert len(comparison[comparison.InputCells_Ileum.isna()]) == 6
assert comparison.loc[comparison.CellType.eq('Glial'), 'InputCells_Ileum'].isna().all()
assignments = pd.read_csv(root / '18样本完整生态型分配.tsv', sep='\t')
assert len(assignments) == 18 and assignments.Sample.is_unique
assert sum(summary['ecotype_counts'].values()) == 18

# Step3. 53个程序和7生态型逐项附来源，作者亚型比例不能超过1
defined = pd.read_csv(interpretation / '53状态_初步判断及来源.tsv', sep='\t')
ecotypes = pd.read_csv(interpretation / '7生态型_初步判断及来源.tsv', sep='\t')
sources = pd.read_csv(interpretation / '文献方法与数据库来源.tsv', sep='\t')
assert len(defined) == 53 and set(defined.CellState) == set(states.CellState)
assert len(ecotypes) == 7 and ecotypes.Ecotype.is_unique
for source_ids in pd.concat([defined.SourceIDs, ecotypes.SourceIDs]):
    assert set(source_ids.split(';')) <= set(sources.SourceID)
cross = pd.read_csv(interpretation / '作者亚型与恢复状态_描述性对应.tsv', sep='\t')
totals = cross.groupby(['Dataset', 'CellType', 'State']).DonorBalancedFractionWithinState.sum()
assert np.allclose(totals, 1)
hashes = pd.read_csv(interpretation / '只读输入校验.tsv', sep='\t')
assert hashes.SHA256Before.eq(hashes.SHA256After).all()
for row in hashes.itertuples(index=False):
    assert hashlib.sha256(Path(row.Input).read_bytes()).hexdigest() == row.SHA256Before

# Step4. 全部17图与阅读用PDF一致，自动审计通过；目视记录单独检查
qa = json.loads((root / 'figure_qa/audit_exit_codes.json').read_text())
assert len(qa) == 17 and all(x['text_exit_code'] == x['collision_exit_code'] == 0 for x in qa)
assert len(PdfReader(root / 'recovery图件汇总.pdf').pages) == 17
manual = json.loads((root / 'figure_qa/final_review.json').read_text(encoding='utf-8'))
assert manual['reviewed_figure_count'] == 17 and manual['all_panels_reviewed']
server_check = json.loads((raw / 'audit/final_data_check.json').read_text())
check = {'returned_files_unchanged': len(manifest), 'cells': len(cells), 'states': len(states), 'ecotypes': len(ecotypes), 'authors_subtype_fractions_sum_to_one': True, 'source_records': len(sources), 'interpretation_inputs_unchanged': len(hashes), 'pdf_figures_checked': 17, 'pdf_packet_pages': 17, 'server_check': server_check, 'accuracy_claim': False, 'clinical_analysis_performed': False}
(interpretation / '最终检查摘要.json').write_text(json.dumps(check, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps(check, ensure_ascii=False, indent=2))
