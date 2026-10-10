# 将结直肠与回肠两轮独立 recovery 对齐到冻结状态，保留缺失谱系和失败状态。
# 在本机 .tools/cibersortx_review 环境运行；读取两轮审阅表与诊断 marker，写出描述性比较。
# 队列、细胞预算和检测深度不同；不以跨队列 Z 排序作为准确率，不重新筛选发现 marker。
from pathlib import Path
import json
import numpy as np
import pandas as pd

# Step1. 以全部 53 个发现状态为骨架；IBDverse 缺失的 Glial 保留 NA
project = Path(__file__).resolve().parents[3]
root = project / 'IBD_ecotype_analysis/06_validation/GSE214695_scRNA_recovery_20261009'
ileum_root = project / 'IBD_ecotype_analysis/06_validation/IBDverse_scRNA_recovery_20261009'
colon = pd.read_csv(root / '53状态单细胞恢复审阅.tsv', sep='\t')
ileum = pd.read_csv(ileum_root / '47状态单细胞恢复审阅.tsv', sep='\t')
identity = ['CellState', 'CellType', 'State', 'InitialProgram', 'ProgramConfidence', 'Top10Genes']
metrics = ['InputCells', 'InputDonors', 'AssignedCells', 'AssignedDonors', 'AssignedFraction', 'Z', 'OwnMinusBestOther_InternalOnly', 'MarkerMeanDetectionInAssignedCells']
comparison = colon[identity].merge(colon[['CellState'] + metrics].add_suffix('_Colon').rename(columns={'CellState_Colon': 'CellState'}), on='CellState', validate='one_to_one')
comparison = comparison.merge(ileum[['CellState'] + metrics].add_suffix('_Ileum').rename(columns={'CellState_Ileum': 'CellState'}), on='CellState', how='left', validate='one_to_one')
comparison['ColonDonorSupportFraction'] = comparison.AssignedDonors_Colon / comparison.InputDonors_Colon
comparison['IleumDonorSupportFraction'] = comparison.AssignedDonors_Ileum / comparison.InputDonors_Ileum
comparison['BothOwnMarkerHigher_InternalOnly'] = comparison.OwnMinusBestOther_InternalOnly_Colon.gt(0) & comparison.OwnMinusBestOther_InternalOnly_Ileum.gt(0)
comparison['ColonOnlyLineage'] = comparison.InputCells_Ileum.isna()
comparison.to_csv(root / '53状态_结肠回肠描述性比较.tsv', sep='\t', index=False)

# Step2. 预先关注的四类，按状态保存谱系 marker 检出和平均表达；保留未分配细胞
focus = ['Mast', 'Goblet', 'Enteroendocrine', 'Glial']
marker_genes = ['TPSAB1', 'TPSB2', 'KIT', 'CPA3', 'CSF3R', 'CXCR2', 'S100A8', 'S100A9', 'CSF1R', 'MS4A7', 'CHGA', 'CHGB', 'TPH1', 'MUC2', 'TFF3', 'MKI67', 'TOP2A', 'NOS2', 'DUOX2', 'LCN2', 'SOX10', 'S100B', 'MPZ', 'CUBN', 'ALDOB']
rows = []
for dataset, data_root in [('Colon_GSE214695', root), ('Ileum_IBDverse', ileum_root)]:
    raw = data_root / 'raw_server'
    diagnostic_path = raw / ('input/diagnostic_marker_with_annotation.tsv' if dataset == 'Colon_GSE214695' else 'diagnostic_marker_with_annotation.tsv')
    diagnostic = pd.read_csv(diagnostic_path, sep='\t', index_col=0)
    cells = pd.read_csv(raw / 'audit/cell_assignments_all.tsv', sep='\t').set_index('ID')
    selected = cells[cells.CellType.isin(focus)].join(diagnostic[[g for g in marker_genes if g in diagnostic]], validate='one_to_one')
    for (cell_type, state), frame in selected.groupby(['CellType', 'State'], sort=True):
        for gene in [g for g in marker_genes if g in selected]:
            donor_mean = frame.groupby('Donor')[gene].mean()
            donor_detection = frame[gene].gt(0).to_frame('Detected').join(frame[['Donor']]).groupby('Donor').Detected.mean()
            rows.append({'Dataset': dataset, 'CellType': cell_type, 'State': state, 'Gene': gene, 'Cells': len(frame), 'Donors': frame.Donor.nunique(), 'CellDetectionFraction': float(frame[gene].gt(0).mean()), 'DonorBalancedDetectionFraction': float(donor_detection.mean()), 'DonorBalancedMeanCPM': float(donor_mean.mean())})
pd.DataFrame(rows).to_csv(root / '重点谱系_逐状态诊断marker.tsv', sep='\t', index=False)

# Step3. 描述性计数用于报告，不构成新增统计验收阈值
shared = comparison[~comparison.ColonOnlyLineage]
summary = {'shared_states': len(shared), 'all_discovery_states': len(comparison), 'both_self_marker_higher_internal_only': int(shared.BothOwnMarkerHigher_InternalOnly.sum()), 'colon_self_marker_higher_internal_only': int(shared.OwnMinusBestOther_InternalOnly_Colon.gt(0).sum()), 'ileum_self_marker_higher_internal_only': int(shared.OwnMinusBestOther_InternalOnly_Ileum.gt(0).sum()), 'glial_states_evaluated_only_in_colon': int(comparison.ColonOnlyLineage.sum()), 'comparison_limit': 'Different cohorts and sampling caps; descriptive recovery comparison, not paired patient analysis or comparative accuracy'}
(root / 'comparison_summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
print(summary)
print(comparison[comparison.CellType.isin(focus)][['CellState', 'AssignedCells_Colon', 'AssignedDonors_Colon', 'Z_Colon', 'AssignedCells_Ileum', 'AssignedDonors_Ileum', 'Z_Ileum']].to_string(index=False))
