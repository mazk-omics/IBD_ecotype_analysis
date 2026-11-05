# 独立核验已有F1/F2、四套相关和置换Z，复用已完成的共关联；补算疾病、炎症、肠段检验。
# 服务器 ibd_refaudit 环境运行；冻结数据只读，全部新结果写入本次目录。
# 细胞级与样本级置换分开；多活检和供者依赖只记录，不把细胞P当成患者P。
from pathlib import Path
import hashlib
import json
import shutil
import numpy as np
import pandas as pd
from scipy import stats

# Step1. 读入现场核验的模型与历史共关联，不重复拟合NMF或重跑相同1000次置换
root = Path('/home/mazekai/IBD_EcoTyper')
out = root / '06_validation/paper_aligned_recovery_characterization_20261012'
prior = root / '06_validation/ecotyper_original_method_annotation_20261010'
states = pd.read_csv(out / '00_audit/frozen_states.tsv', sep='\t')
coassoc = pd.read_csv(prior / 'ecotype_coassociation_both_datasets.tsv', sep='\t')
recovery_rows = []
coverage_rows = []
qc = []
phenotypes = []
rng = np.random.default_rng(20261012)
for dataset, folder in [('Colon_GSE214695', 'GSE214695_scRNA_recovery_20261009'), ('Ileum_IBDverse', 'IBDverse_scRNA_recovery_20261009')]:
    cells = pd.read_csv(out / '00_audit' / (dataset + '_cells.tsv'), sep='\t')
    official = pd.read_csv(root / '06_validation' / folder / 'audit/state_recovery_summary.tsv', sep='\t').set_index('CellState')
    f1 = pd.read_csv(prior / (dataset + '_F1_assigned_denominator.tsv'), sep='\t', index_col=0)
    f2 = pd.read_csv(prior / (dataset + '_F2_all_cells_denominator.tsv'), sep='\t', index_col=0)
    for row in states.itertuples():
        parent = cells[cells.CellType.eq(row.CellType)]
        own = parent[parent.State.eq(row.State)]
        recovery = {'Dataset': dataset, 'CellState': row.CellState, 'CellType': row.CellType, 'Ecotype': row.Ecotype, 'InputCells': len(parent), 'AssignedCells': len(own), 'InputDonors': parent.Donor.nunique(), 'AssignedDonors': own.Donor.nunique(), 'InputSamples': parent.Sample.nunique(), 'AssignedSamples': own.Sample.nunique(), 'UnassignedCells': int(parent.State.eq('Unassigned').sum()), 'Status': 'EVALUABLE' if len(parent) else 'NOT_EVALUABLE', 'OfficialP': np.nan, 'SignificancePermutations': 30, 'PaperPermutations': 100}
        if row.CellState in official.index:
            value = official.loc[row.CellState]
            assert int(value.AssignedCells) == len(own) and int(value.InputCells) == len(parent)
            recovery.update({'OfficialZ': value.Z, 'PositiveZgt1p96': value.Z > 1.96, 'MarkerGenes': value.MarkerGenes, 'MarkerGenesDetectedInLineage': value.MarkerGenesDetectedInLineage})
        recovery_rows.append(recovery)
        if row.CellState not in f1.index:
            assert len(parent) == 0
            continue
        numerator = own.groupby('Sample').size().reindex(f1.columns, fill_value=0)
        d1 = parent[parent.State.ne('Unassigned')].groupby('Sample').size().reindex(f1.columns, fill_value=0)
        d2 = parent.groupby('Sample').size().reindex(f1.columns, fill_value=0)
        np.testing.assert_allclose(f1.loc[row.CellState], numerator/d1.replace(0, np.nan), atol=1e-12, equal_nan=True)
        np.testing.assert_allclose(f2.loc[row.CellState], numerator/d2.replace(0, np.nan), atol=1e-12, equal_nan=True)
    matrices = [f1, f2, f1.replace(0, np.nan), f2.replace(0, np.nan)]
    expected = sum(x.T.corr(method='spearman', min_periods=2).to_numpy() for x in matrices)/4
    np.fill_diagonal(expected, np.nan)
    c = pd.read_csv(prior / (dataset + '_C_four_matrix_mean.tsv'), sep='\t', index_col=0)
    np.testing.assert_allclose(c, expected, atol=1e-12, equal_nan=True)
    null = pd.read_csv(prior / (dataset + '_permutation_null.tsv'), sep='\t')
    assert null.shape == (1000, 7)
    for row in coassoc[coassoc.Dataset.eq(dataset)].itertuples():
        inside = c.index.intersection(states.loc[states.Ecotype.eq(row.Ecotype), 'CellState'])
        outside = c.columns.difference(inside)
        delta = c.loc[inside, inside].mean(axis=1)-c.loc[inside, outside].mean(axis=1)
        np.testing.assert_allclose(delta.mean(), row.AvgDif, atol=1e-12, equal_nan=True)
        np.testing.assert_allclose((row.AvgDif-null[row.Ecotype].mean())/null[row.Ecotype].std(ddof=1), row.Z, atol=1e-12, equal_nan=True)
        frozen = states.loc[states.Ecotype.eq(row.Ecotype), 'CellState']
        coverage_rows.append({'Dataset': dataset, 'Ecotype': row.Ecotype, 'FrozenMembers': ';'.join(frozen), 'InputLineageMembers': ';'.join(inside), 'MissingMembers': ';'.join(frozen[~frozen.isin(inside)]), 'DefinedReferenceMembers': int(delta.notna().sum())})
    for path in prior.glob(dataset + '*.tsv'):
        shutil.copyfile(path, out / '05_ecotype_coassociation' / path.name)
    qc.append({'Dataset': dataset, 'F1F2Denominators': 'PASS', 'FourSpearmanAverage': 'PASS', 'NullZ': 'PASS', 'DefinedUndirectedEdges': int(np.isfinite(c.to_numpy()).sum()//2)})

    # Step2. 原文式1000次标签打乱；仅在同父类型且成功分配的细胞内检验，不引入跨谱系对照
    for cell_type, parent in cells[cells.State.ne('Unassigned')].groupby('CellType'):
        for unit in ['Cell', 'BiologicalSampleDominant']:
            if unit == 'Cell':
                observations = parent[['State', 'DiseaseGroup', 'InflammationGroup', 'Sample', 'Donor']].copy()
            else:
                counts = parent.groupby(['Sample', 'State']).size().unstack(fill_value=0)
                tied = counts.eq(counts.max(axis=1), axis=0).sum(axis=1) > 1
                counts = counts.loc[~tied]
                observations = parent[['Sample', 'Donor', 'DiseaseGroup', 'InflammationGroup']].drop_duplicates().set_index('Sample').loc[counts.index].reset_index()
                observations['State'] = counts.idxmax(axis=1).to_numpy()
            for feature in ['DiseaseGroup', 'InflammationGroup']:
                targets = sorted(observations[feature].dropna().unique())
                labels = observations.State.to_numpy()
                unique_states = sorted(set(labels))
                reference = np.array([[np.sum((labels==state) & observations[feature].eq(target).to_numpy()) for state in unique_states] for target in targets])
                nulls = np.zeros((1000, len(targets), len(unique_states)))
                for perm in range(1000):
                    shuffled = rng.permutation(labels)
                    nulls[perm] = [[np.sum((shuffled==state) & observations[feature].eq(target).to_numpy()) for state in unique_states] for target in targets]
                means = nulls.mean(axis=0)
                sd = nulls.std(axis=0, ddof=1)
                for a, target in enumerate(targets):
                    for b, state in enumerate(unique_states):
                        z = (reference[a, b]-means[a, b])/sd[a, b] if sd[a, b]>0 else np.nan
                        phenotypes.append({'Dataset': dataset, 'CellState': cell_type+'_'+state, 'Unit': unit, 'Feature': feature, 'Target': target, 'Ns': int(reference[a, b]), 'NullMean': means[a, b], 'NullSD': sd[a, b], 'Z': z, 'PositiveZgt1p96': bool(z>1.96), 'Permutations': 1000, 'Seed': 20261012, 'Method': 'PAPER_ADAPTED_IBD', 'Dependence': 'Cell-level donor dependence' if unit=='Cell' else 'Dominant per sample; ties omitted; exploratory unit adaptation'})
        print('Phenotypes complete:', dataset, cell_type, flush=True)
pd.DataFrame(recovery_rows).to_csv(out / '01_state_recovery/state_recovery_evidence_by_dataset.tsv', sep='\t', index=False)
pd.DataFrame(coverage_rows).to_csv(out / '05_ecotype_coassociation/ecotype_member_coverage.tsv', sep='\t', index=False)
coassoc['PvalueStatus'] = 'Not supplied: original tail ambiguity; Z retained'
coassoc.to_csv(out / '05_ecotype_coassociation/ecotype_coassociation_1000permutations.tsv', sep='\t', index=False)
pd.DataFrame(qc).to_csv(out / '00_audit/existing_coassociation_independent_QC.tsv', sep='\t', index=False)
pd.DataFrame(phenotypes).to_csv(out / '04_state_phenotypes/state_phenotype_enrichment.tsv', sep='\t', index=False)

# Step3. Discovery dominant state 与Region3 Fisher；不存在的样本分配不能补成某一状态
annotation = pd.read_csv(root / '05_EcoTyper/input/MSCCR_CORE12/region_ibd/annotation/MSCCR_CORE12_IBD_annotation.tsv', sep='\t')
assert annotation.ID.is_unique and len(annotation)==2025
region_rows = []
for cell_type in states.CellType.unique():
    assignments = pd.read_csv(root / '05_EcoTyper/output/MSCCR_CORE12/ctspec_region_ibd' / cell_type / 'state_assignment.txt', sep='\t')
    assert assignments.ID.is_unique and set(assignments.ID) <= set(annotation.ID)
    observations = assignments.merge(annotation, on='ID', validate='one_to_one')
    for state in states.loc[states.CellType.eq(cell_type), 'State']:
        for region in sorted(annotation.Region3.unique()):
            own = observations.State.eq(state)
            target = observations.Region3.eq(region)
            table = np.array([[np.sum(own & target), np.sum(own & ~target)], [np.sum(~own & target), np.sum(~own & ~target)]])
            odds, p = stats.fisher_exact(table, alternative='greater')
            region_rows.append({'CellState': cell_type+'_'+state, 'CellType': cell_type, 'Region3': region, 'AssignedSampleDenominator': len(observations), 'PatientCount': observations.PatientID.nunique(), 'Ns': int(table[0, 0]), 'OddsRatio': odds, 'P': p, 'FisherTail': 'greater; explicit IBD enrichment adaptation', 'BHFamily': 'All 53 states x 3 regions'})
region_results = pd.DataFrame(region_rows)
region_results['Q'] = stats.false_discovery_control(region_results.P.to_numpy(), method='bh')
region_results.to_csv(out / '04_state_phenotypes/state_region_fisher.tsv', sep='\t', index=False)

# Step4. Discovery疾病/炎症以连续ecotype abundance Wilcoxon；Region3用三组KW，不混用dominant Fisher
abundance = pd.read_csv(root / 'software/ecotyper/EcoTyper/MSCCR_CORE12_CTSPEC_REGION_IBD/Cell_type_specific_genes/Ecotypes/discovery/ecotype_abundance.txt', sep='\t', index_col=0).T
assert set(abundance.index)==set(annotation.ID) and abundance.shape==(2025, 7)
joined = annotation.set_index('ID').join(abundance, validate='one_to_one')
dominant = abundance.idxmax(axis=1)
tie = abundance.eq(abundance.max(axis=1), axis=0).sum(axis=1)>1
joined['DominantEcotype'] = dominant.where(~tie)
joined.to_csv(out / '06_ecotype_features/ecotype_dominant_assignments.tsv', sep='\t')
feature_rows = []
for feature, group_a, group_b in [('Disease', 'CD', 'UC'), ('Inflammation', 'Inflamed', 'Non-inflamed')]:
    for ecotype in abundance.columns:
        a = joined.loc[joined[feature].eq(group_a), ecotype]
        b = joined.loc[joined[feature].eq(group_b), ecotype]
        test = stats.mannwhitneyu(a, b, alternative='two-sided', method='asymptotic')
        feature_rows.append({'Dataset': 'MSCCR_discovery', 'Feature': feature, 'Ecotype': ecotype, 'GroupA': group_a, 'GroupB': group_b, 'N_A': len(a), 'N_B': len(b), 'MeanDifference': a.mean()-b.mean(), 'P': test.pvalue, 'Unit': 'Biopsy', 'RepeatedPatientCaveat': 'Multiple biopsies per patient; not patient-independent clinical evidence'})
features = pd.DataFrame(feature_rows)
features['Q'] = stats.false_discovery_control(features.P.to_numpy(), method='bh')
features['BHFamily'] = 'All 7 ecotypes x 2 discrete features'
features.to_csv(out / '06_ecotype_features/ecotype_discrete_features.tsv', sep='\t', index=False)
kw = []
for ecotype in abundance.columns:
    test = stats.kruskal(*(joined.loc[joined.Region3.eq(region), ecotype] for region in sorted(joined.Region3.unique())))
    kw.append({'Ecotype': ecotype, 'Feature': 'Region3', 'Statistic': test.statistic, 'P': test.pvalue, 'Method': 'PAPER_ADAPTED_IBD: 3-group Kruskal-Wallis'})
kw = pd.DataFrame(kw)
kw['Q'] = stats.false_discovery_control(kw.P.to_numpy(), method='bh')
kw.to_csv(out / '06_ecotype_features/ecotype_region_KW.tsv', sep='\t', index=False)
(out / '00_audit/reuse_phenotypes_complete.json').write_text(json.dumps({'QC': qc, 'StatePhenotypeRows': len(phenotypes), 'DiscoveryBiopsies': len(annotation), 'DiscoveryPatients': annotation.PatientID.nunique(), 'MaxBiopsiesPerPatient': int(annotation.PatientID.value_counts().max()), 'SourceFilesWritten': False}, indent=2))
print('REUSE_AND_PHENOTYPE_COMPLETE', flush=True)
