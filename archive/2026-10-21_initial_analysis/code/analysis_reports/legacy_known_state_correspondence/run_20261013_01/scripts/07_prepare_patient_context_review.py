# 小规模候选状态的MSCCR患者水平背景审阅；不读取benchmark输入或重新拟合模型。
# 用已冻结的53状态丰度审计表，患者×疾病×肠段内汇总，避免把重复活检当独立样本。
import os
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
from pathlib import Path
import hashlib
import pandas as pd

# Step1. 确认旧版2025活检来源与ID，不能把GSE72819恢复快照当MSCCR
out = Path(__file__).resolve().parents[1]
source = Path('D:/IBD_ecotype/Profile/audit/2026-10-06_ctspec_region_ibd/tables')
matrix_path = source / 'combined_H_for_audit.tsv'
meta_path = source / '01_discovery_sample_metadata.tsv'
matrix = pd.read_csv(matrix_path, sep='\t', index_col=0)
meta = pd.read_csv(meta_path, sep='\t')
assert matrix.index.is_unique and meta.ID.is_unique and len(meta) == 2025
assert set(matrix.columns) == set(meta.ID) and len(matrix) == 53
assert meta.groupby('PatientID').Disease.nunique().le(1).all()
pd.DataFrame([{'Path': str(path), 'SHA256': hashlib.sha256(path.read_bytes()).hexdigest(), 'SourceScope': 'frozen legacy MSCCR2025 discovery audit; read-only'} for path in [matrix_path, meta_path]]).to_csv(out / 'audit/MSCCR_PATIENT_REVIEW_INPUT_HASHES.tsv', sep='\t', index=False)

# Step2. 预先限定7状态，不做53状态临床筛选后选择最佳P
candidates = ['Monocyte_Macrophage_S02', 'Monocyte_Macrophage_S03', 'Goblet_S01', 'CD8_T_S03', 'CD4_T_S04', 'Fibroblast_S02', 'Fibroblast_S03']
long = matrix.loc[candidates].rename_axis('CellState').reset_index().melt(id_vars='CellState', var_name='ID', value_name='FrozenWithinParentStateFraction')
long = long.merge(meta[['ID', 'PatientID', 'Disease', 'Region3', 'Inflammation']], on='ID', validate='many_to_one')
assert long.FrozenWithinParentStateFraction.dropna().between(0, 1 + 1e-8).all()
summary = long.groupby(['CellState', 'PatientID', 'Disease', 'Region3', 'Inflammation']).agg(StateFractionMedian=('FrozenWithinParentStateFraction', 'median'), AvailableBiopsies=('FrozenWithinParentStateFraction', 'count'), TotalBiopsies=('ID', 'size')).reset_index()
summary.to_csv(out / 'audit/MSCCR_PATIENT_REGION_INFLAMMATION_MEDIANS.tsv', sep='\t', index=False)
pairs = summary.pivot(index=['CellState', 'PatientID', 'Disease', 'Region3'], columns='Inflammation', values='StateFractionMedian').reset_index()
assert 'Inflamed' in pairs and 'Non-inflamed' in pairs
pairs['WithinPatientDifference'] = pairs.Inflamed - pairs['Non-inflamed']
pairs['EvaluablePair'] = pairs.Inflamed.notna() & pairs['Non-inflamed'].notna()
pairs.to_csv(out / 'audit/MSCCR_WITHIN_PATIENT_INFLAMMATION_PAIRS.tsv', sep='\t', index=False)
disease = long.groupby(['CellState', 'PatientID', 'Disease', 'Region3']).FrozenWithinParentStateFraction.median().rename('StateFractionMedian').reset_index()
disease.to_csv(out / 'audit/MSCCR_PATIENT_REGION_DISEASE_MEDIANS.tsv', sep='\t', index=False)
print('Biopsy IDs:', len(meta), '; unique patients:', meta.PatientID.nunique())
print(pairs.groupby(['CellState', 'Disease', 'Region3']).EvaluablePair.sum().to_string())
