# 就地审阅已回传的官方 recovery，区分状态分配、表达支持和供者覆盖。
# 在本地 .tools/cibersortx_review 环境运行；输入 raw_server/audit，输出逐状态审阅表和简要数值记录。
# 只作探索性描述，不重新命名状态，不删除弱状态，不以单细胞数量代替供者数量做疾病检验。
from pathlib import Path
import json
import numpy as np
import pandas as pd

# Step1. 合并官方统计和既有的初步程序注释
project = Path(__file__).resolve().parents[3]
root = project / 'IBD_ecotype_analysis/06_validation/IBDverse_scRNA_recovery_20261009'
audit = root / 'raw_server/audit'
states = pd.read_csv(audit / 'state_recovery_summary.tsv', sep='\t')
bulk = pd.read_csv(project / 'IBD_ecotype_analysis/06_validation/GSE72819_recovery_review_20261009/53状态注释与恢复审阅.tsv', sep='\t')
states = states.merge(bulk[['CellState', 'WinnerSamples', 'InitialProgram', 'ProgramConfidence', 'ReviewReason', 'Top10Genes']], on='CellState', validate='one_to_one')
cells = pd.read_csv(audit / 'cell_assignments_all.tsv', sep='\t')
programs = pd.read_csv(audit / 'state_program_expression.tsv', sep='\t')

# Step2. 查看冻结 marker 程序是否在自身状态中比其他状态更高；这只是同源表达描述
contrast_rows = []
for marker_state, frame in programs.groupby('MarkerState', sort=False):
    own = frame.loc[frame.TargetState.eq(marker_state), 'RecoveryMeanZ'].iloc[0]
    other = frame.loc[frame.TargetState.ne(marker_state), 'RecoveryMeanZ']
    best_other = other.max()
    contrast_rows.append({'CellState': marker_state, 'OwnMarkerMeanZ': own, 'HighestOtherStateMarkerMeanZ': best_other, 'OwnMinusBestOther_InternalOnly': own - best_other})
states = states.merge(pd.DataFrame(contrast_rows), on='CellState', validate='one_to_one')
states['FiniteOfficialZ'] = np.isfinite(states.Z)
states['AboveNominalZ1_96'] = states.FiniteOfficialZ & states.Z.gt(1.96)
states['RecoveryReview'] = '有分配细胞；官方 Z 未超过名义参考，保留为待核对状态'
states.loc[states.AssignedCells.eq(0), 'RecoveryReview'] = '未分配到细胞；本轮不能验证，不能据此自动删除发现状态'
states.loc[states.AboveNominalZ1_96, 'RecoveryReview'] = '官方 Z 超过名义参考；仍需检查谱系身份、供者覆盖和外部队列'
states.loc[states.AssignedCells.gt(0) & ~states.FiniteOfficialZ, 'RecoveryReview'] = '有分配细胞但 Z 不可计算；支持细胞不足或空零分布需核对'
states.to_csv(root / '47状态单细胞恢复审阅.tsv', sep='\t', index=False)
type_summary = cells.groupby('CellType').agg(InputCells=('ID', 'size'), InputDonors=('Donor', 'nunique'))
type_summary['AssignedCells'] = cells[cells.State.ne('Unassigned')].groupby('CellType').size().reindex(type_summary.index, fill_value=0)
type_summary['AssignedFraction'] = type_summary.AssignedCells / type_summary.InputCells
type_summary['EvaluatedStates'] = states.groupby('CellType').size()
type_summary['StatesWithWinners'] = states[states.AssignedCells.gt(0)].groupby('CellType').size().reindex(type_summary.index, fill_value=0)
type_summary['StatesAboveNominalZ1_96'] = states[states.AboveNominalZ1_96].groupby('CellType').size().reindex(type_summary.index, fill_value=0)
type_summary.to_csv(root / '11类型恢复摘要.tsv', sep='\t')
print(type_summary.to_string())

# Step3. 核对官方生态型分配；未输出的样本在完整样本集合中记为未分配
eco = pd.read_csv(audit / 'ecotype_assignment.txt', sep='\t', index_col=0)
samples = cells[['Sample', 'Donor', 'Stratum']].drop_duplicates()
samples['Ecotype'] = samples.Sample.map(eco.set_index('ID').Ecotype).fillna('Unassigned')
samples = samples.merge(eco.drop(columns='Ecotype'), left_on='Sample', right_on='ID', how='left', validate='one_to_one').drop(columns='ID')
samples.to_csv(root / '30活检完整生态型分配.tsv', sep='\t', index=False)
ecotype_counts = samples.Ecotype.value_counts().reindex([f'E{i}' for i in range(1, 8)] + ['Unassigned'], fill_value=0)
ecotype_counts.to_csv(root / '生态型计数.tsv', sep='\t')
summary = {'input_cells': len(cells), 'assigned_cells': int(cells.State.ne('Unassigned').sum()), 'assigned_cell_fraction': float(cells.State.ne('Unassigned').mean()), 'states': len(states), 'states_with_winners': int(states.AssignedCells.gt(0).sum()), 'states_with_finite_z': int(states.FiniteOfficialZ.sum()), 'states_above_nominal_z_1_96': int(states.AboveNominalZ1_96.sum()), 'states_own_marker_higher_than_best_other_internal_only': int(states.OwnMinusBestOther_InternalOnly.gt(0).sum()), 'ecotype_counts': {str(k): int(v) for k, v in ecotype_counts.items()}, 'official_shuffle_repeats': 30}
summary['assigned_biopsies_with_assignment_q_ge_0_25'] = int((samples.Ecotype.ne('Unassigned') & samples.AssignmentQ.ge(0.25)).sum())
(root / 'review_summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
print(summary)
print(states[~states.AboveNominalZ1_96][['CellState', 'AssignedCells', 'AssignedDonors', 'Z']].to_string(index=False))
