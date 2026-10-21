# 只读核查Joint20冻结清单与主线样本，生成本次审计表；不改原参考或运行目录。
# 在本机项目目录运行，读取现有metadata和212列CPM，不读取完整单细胞矩阵。
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

# Step1. 读入冻结清单，核对profile和供者隔离
root = Path('D:/IBD_ecotype/IBD_ecotype_analysis')
out = root / 'analysis_reports/joint20_v2_sampling_bp_reference_audit_20261014'
source = root / 'CIBERSORTx_ready_20261007/06_Joint20v2_DonorSite_reference/核对记录_勿上传'
cells = pd.read_csv(source / 'selected_cell_manifest.tsv', sep='\t')
profiles = pd.read_csv(source / 'profile_manifest.tsv', sep='\t')
roles = pd.read_csv(source / 'frozen_donor_roles.tsv', sep='\t')
types = json.loads((source / 'construction_summary.json').read_text(encoding='utf-8'))['types']
assert cells.cell_id.is_unique and len(cells) == 9294
assert cells.groupby('profile_id')[['Patient', 'Site', 'sample_id', 'coarse_identity']].nunique().eq(1).all().all()
observed = cells.groupby('profile_id').size().reindex(profiles.profile_id)
assert np.array_equal(observed, profiles.selected_cells) and observed.between(20, 50).all()
assert not cells.Patient.isin(roles.loc[roles.split.ne('training'), 'Patient']).any()
profiles.to_csv(out / 'profile_counts.tsv', sep='\t', index=False)

# Step2. 汇总父类型、疾病、炎症和实际部位
summary = cells.groupby('coarse_identity').agg(cells=('cell_id', 'size'), donors=('Patient', 'nunique'), samples=('sample_id', 'nunique'), profiles=('profile_id', 'nunique'), subtypes=('final_analysis', 'nunique')).reindex(types)
for field in ['Disease', 'Inflammation', 'Site', 'region']:
    values = cells[field].fillna('UNKNOWN')
    table = pd.crosstab(cells.coarse_identity, values).reindex(types, fill_value=0)
    table.to_csv(out / f'type_by_{field}.tsv', sep='\t')
    summary = summary.join(table.add_prefix(field + '_'))
    totals = cells.assign(audit_value=values).groupby('audit_value').agg(cells=('cell_id', 'size'), donors=('Patient', 'nunique'), samples=('sample_id', 'nunique'), profiles=('profile_id', 'nunique'))
    totals.to_csv(out / f'overall_{field}.tsv', sep='\t')
contribution = cells.groupby(['coarse_identity', 'Patient']).agg(cells=('cell_id', 'size'), profiles=('profile_id', 'nunique'), samples=('sample_id', 'nunique')).reset_index()
contribution['within_type_fraction'] = contribution.cells / contribution.coarse_identity.map(summary.cells)
contribution.to_csv(out / 'type_donor_contribution.tsv', sep='\t', index=False)
summary['max_donor_fraction'] = contribution.groupby('coarse_identity').within_type_fraction.max()
summary['top3_donor_fraction'] = contribution.groupby('coarse_identity').within_type_fraction.apply(lambda x: x.nlargest(3).sum())
summary['effective_donors'] = 1 / contribution.groupby('coarse_identity').within_type_fraction.apply(lambda x: (x ** 2).sum())
summary.to_csv(out / 'type_summary.tsv', sep='\t')
cells.groupby('coarse_identity').agg(cells=('cell_id', 'size'), min_counts=('total_counts', 'min'), median_counts=('total_counts', 'median'), median_genes=('n_genes_by_counts', 'median'), median_mito_pct=('pct_counts_mt', 'median'), max_mito_pct=('pct_counts_mt', 'max')).reindex(types).to_csv(out / 'cell_QC_by_type.tsv', sep='\t')
cells.groupby(['coarse_identity', 'Disease', 'Inflammation', 'Site'], dropna=False).agg(cells=('cell_id', 'size'), donors=('Patient', 'nunique'), samples=('sample_id', 'nunique'), profiles=('profile_id', 'nunique')).to_csv(out / 'type_disease_inflammation_site.tsv', sep='\t')
cells.groupby(['coarse_identity', 'region', 'final_analysis']).size().rename('cells').to_csv(out / 'selected_author_subtypes.tsv', sep='\t')
roles.groupby(['split', 'validation_role', 'Disease']).size().rename('donors').to_csv(out / 'donor_roles_summary.tsv', sep='\t')

# Step3. 实际SampleID逐一核对主线结果与当前manifest
run = root / '06_validation/joint20_bayesprism_ecotyper_matched_sensitivity_20261014'
annotation = pd.read_csv(run / '00_audit/MSCCR_CD_UC_discovery_manifest.tsv', sep='\t')
main = pd.read_csv(root / 'CIBERSORTx_ready_20261007/09_MSCCR_Joint20v2_Fractions/核对记录_勿上传/MSCCR_CD_UC_2025_annotation.tsv', sep='\t')
job16 = pd.read_csv(root / 'CIBERSORTx_signature_audit/MSCCR_Joint20v2_Job16_review_20261008/original_outputs/CIBERSORTx_Job16_Results.txt', sep='\t')
assert annotation.ID.tolist() == main.ID.tolist() == job16.iloc[:, 0].tolist()
assert annotation.Disease.isin(['CD', 'UC']).all() and annotation.ID.is_unique
comparison = annotation[['ID', 'Disease', 'Inflammation', 'Region3']].copy()
comparison['mainline_order'] = np.arange(1, len(main) + 1)
comparison['Job16_ID'] = job16.iloc[:, 0].values
comparison['exact_ID_match'] = comparison.ID.eq(comparison.Job16_ID)
comparison.to_csv(out / 'sampleID_match_2025.tsv', sep='\t', index=False)

# Step4. 复核已保存参考表达中的典型marker与表达相近的父类型
matrix = pd.read_csv(source.parent / 'GJ20v2_DonorSite_CPM.txt', sep='\t', index_col=0)
assert matrix.columns.tolist() == profiles.profile_id.tolist()
panels = {'B': ['MS4A1', 'CD79A', 'CD74'], 'Plasma': ['MZB1', 'JCHAIN', 'SDC1'], 'CD4_T': ['CD3D', 'CD3E', 'IL7R'], 'CD8_T': ['CD3D', 'CD8A', 'CD8B'], 'NK': ['NKG7', 'GNLY', 'KLRD1'], 'ILC': ['IL7R', 'KIT', 'KLRB1'], 'Monocyte_Macrophage': ['LST1', 'CD68', 'FCER1G'], 'DC': ['CD1C', 'FCER1A', 'CLEC10A'], 'Mast': ['TPSAB1', 'TPSB2', 'CPA3'], 'Fibroblast': ['COL1A1', 'COL3A1', 'DCN'], 'Pericyte': ['RGS5', 'PDGFRB', 'CSPG4'], 'Endothelial': ['PECAM1', 'VWF', 'KDR'], 'Ileal_absorptive': ['ALPI', 'APOA4', 'FABP6'], 'Colonic_absorptive': ['KRT20', 'CA1', 'SLC26A3'], 'Goblet': ['MUC2', 'TFF3', 'SPINK4'], 'Stem_TA_progenitor': ['LGR5', 'OLFM4', 'MKI67', 'TOP2A'], 'Paneth': ['DEFA5', 'DEFA6', 'LYZ', 'REG3A'], 'Enteroendocrine': ['CHGA', 'CHGB', 'NEUROD1'], 'Glial': ['SOX10', 'S100B', 'PLP1'], 'Tuft': ['POU2F3', 'TRPM5', 'AVIL']}
marker_records = []
for target, markers in panels.items():
    own = profiles.loc[profiles.coarse_identity.eq(target), 'profile_id']
    for gene in markers:
        if gene not in matrix.index:
            marker_records.append({'parent_type': target, 'gene': gene, 'present': False})
            continue
        means = matrix.loc[gene].groupby(profiles.set_index('profile_id').coarse_identity).mean()
        marker_records.append({'parent_type': target, 'gene': gene, 'present': True, 'own_mean_CPM': matrix.loc[gene, own].mean(), 'own_median_CPM': matrix.loc[gene, own].median(), 'own_positive_profiles': int(matrix.loc[gene, own].gt(0).sum()), 'profiles': len(own), 'top_type': means.idxmax(), 'top_type_mean_CPM': means.max()})
pd.DataFrame(marker_records).to_csv(out / 'marker_profile_audit.tsv', sep='\t', index=False)
means = matrix.T.groupby(profiles.set_index('profile_id').coarse_identity).mean().T
corr = np.log1p(means).corr(method='pearson')
pairs = [{'type1': a, 'type2': b, 'log1p_CPM_Pearson_r': corr.loc[a, b]} for i, a in enumerate(corr.index) for b in corr.index[i + 1:]]
pd.DataFrame(pairs).sort_values('log1p_CPM_Pearson_r', ascending=False).to_csv(out / 'parent_expression_overlap.tsv', sep='\t', index=False)

# Step5. 保存核查结论和输入校验，不修改任何来源文件
result = {'cells': len(cells), 'donors': cells.Patient.nunique(), 'samples': cells.sample_id.nunique(), 'profiles': len(profiles), 'profile_min': int(observed.min()), 'profile_max': int(observed.max()), 'profiles_at_50': int(observed.eq(50).sum()), 'profile_size_distribution': observed.value_counts().sort_index().to_dict(), 'Disease': cells.Disease.value_counts().to_dict(), 'Inflammation': cells.Inflammation.value_counts(dropna=False).to_dict(), 'Site': cells.Site.value_counts().to_dict(), 'Treatment': cells.Treatment.value_counts(dropna=False).to_dict(), 'selected_donors': sorted(cells.Patient.unique().tolist()), 'all_training_donors_used': set(cells.Patient) == set(roles.loc[roles.split.eq('training'), 'Patient']), 'mainline_ID_set_exact': True, 'mainline_ID_order_exact': True, 'Job16_ID_order_exact': True, 'ordered_sample_sha256': hashlib.sha256(('\n'.join(annotation.ID) + '\n').encode()).hexdigest(), 'no_frozen_holdout_donor_overlap': True}
(out / 'local_metadata_audit.json').write_text(json.dumps(result, ensure_ascii=False, indent=2, default=int) + '\n', encoding='utf-8')
print(json.dumps(result, ensure_ascii=False, indent=2, default=int))
print(summary.to_string())
