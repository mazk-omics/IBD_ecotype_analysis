# 单线程读现有小型压缩汇总，检查优先候选的实际证据。
# 不读原始表达矩阵，不重跑DE；汇总不当作单细胞共表达。
import os
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
from pathlib import Path
import pandas as pd

# Step1. 读取全状态表达与供者配对方向
root = Path('D:/IBD_ecotype/IBD_ecotype_analysis')
out = Path(__file__).resolve().parents[1]
base = root / '06_validation/all_states_lineage_audit_20261012'
expression = pd.read_csv(base / 'tables/all_state_gene_expression.tsv.gz', sep='\t')
paired = pd.read_csv(base / 'tables/all_associated_genes_paired_donor_replication.tsv.gz', sep='\t')
author = pd.read_csv(base / 'tables/author_lineage_gene_expression.tsv.gz', sep='\t')
subtypes = pd.read_csv(base / 'tables/state_author_subtype_composition.tsv', sep='\t')
print('Existing summary rows:', len(expression), len(paired), len(author), flush=True)

# Step2. 核心线索与替代谱系marker分别检查
panels = {'Monocyte_Macrophage': ['PLA2G2D', 'ABCA7', 'P2RY13', 'MERTK', 'CSF1R', 'CD68', 'MRC1', 'C1QA', 'FLT3', 'CD1C', 'FCER1A', 'CLEC10A', 'CD1E', 'LILRA4', 'CLEC9A', 'BATF3', 'IRF8'], 'Goblet': ['MKI67', 'TOP2A', 'FOXM1', 'E2F2', 'MYBL2', 'PKMYT1', 'MUC2', 'SPDEF', 'TFF3', 'AGR2', 'OLFM4', 'LGR5', 'SOX9', 'PCNA'], 'CD8_T': ['IL23R', 'JUND', 'CD27', 'LY9', 'LAPTM5', 'RORC', 'IL17A', 'IL17F', 'IL22', 'CCR6', 'KLRB1', 'CD8A', 'CD8B', 'TRAC', 'CD3D', 'TCF7', 'LEF1', 'GZMK', 'GZMB'], 'Fibroblast': ['CCL19', 'CCL21', 'FDCSP', 'CXCL13', 'FN1', 'ECM1', 'ECM2', 'ITGA11', 'HMOX1', 'NDRG1', 'FMO1', 'OAT', 'IL11', 'IL13RA2', 'PDPN', 'FAP', 'THY1', 'ABCA8', 'PDGFRA'], 'CD4_T': ['CD69', 'KLRB1', 'TCF7', 'LEF1', 'FOXP3', 'CTLA4', 'IL2RA', 'BATF', 'STAT1', 'GBP5', 'RUNX3', 'CD4', 'CD3D']}
selected = []
for cell_type, genes in panels.items():
    part = expression[expression.CellState.str.startswith(cell_type + '_') & expression.Gene.isin(genes)].copy()
    part['Parent'] = cell_type
    selected.append(part)
evidence = pd.concat(selected, ignore_index=True)
evidence = evidence.merge(paired, on=['Dataset', 'CellState', 'Gene'], how='left', suffixes=('', '_Paired'), validate='one_to_one')
evidence.to_csv(out / 'audit/CANDIDATE_GENE_PROBE.tsv', sep='\t', index=False)
for state, genes in [('Monocyte_Macrophage_S02', ['PLA2G2D', 'ABCA7', 'P2RY13', 'CSF1R', 'CD68', 'MERTK', 'FLT3', 'CD1C', 'FCER1A', 'CLEC10A']), ('Goblet_S01', ['FOXM1', 'MKI67', 'TOP2A', 'MUC2', 'SPDEF', 'TFF3', 'OLFM4']), ('CD8_T_S03', ['IL23R', 'JUND', 'CD27', 'RORC', 'IL17A', 'KLRB1', 'CD8A', 'CD8B'])]:
    print(state, flush=True)
    print(evidence[evidence.CellState.eq(state) & evidence.Gene.isin(genes)][['Dataset', 'Gene', 'DetectionFraction', 'MeanLog2CPM', 'MatchedDonors_Min5CellsEach', 'MeanPairedLog2CPM_Difference', 'FractionMatchedDonors_Positive']].to_string(index=False), flush=True)
    print(subtypes[subtypes.CellState.eq(state)].sort_values('DonorBalancedFraction', ascending=False).head(6).to_string(index=False), flush=True)
for gene in ['PLA2G2D', 'IL23R']:
    print('Author lineage:', gene, flush=True)
    print(author[author.Gene.eq(gene)].sort_values(['Dataset', 'DonorBalancedMeanLog2CPM'], ascending=[True, False]).groupby('Dataset').head(6)[['Dataset', 'ComparisonLineage', 'DonorBalancedDetectionFraction', 'DonorBalancedMeanLog2CPM']].to_string(index=False), flush=True)
