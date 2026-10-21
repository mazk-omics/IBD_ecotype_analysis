# 单线程整合发表完整导出表、同谱系表达和供者方向；不拟合或重跑DE。
# 统计匹配是冻结特征空间内的候选对应，不是独立身份验证或新颖性检验。
import os
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
from pathlib import Path
import csv
import hashlib
import math
import re
import xml.etree.ElementTree as ET
import pandas as pd
import numpy as np

# Step1. 读冻结身份与既有全基因汇总，不再读取原始单细胞矩阵
root = Path('D:/IBD_ecotype/IBD_ecotype_analysis')
out = Path(__file__).resolve().parents[1]
old = root / 'analysis_reports/legacy_bayesprism_annotation/run_20261013_01'
base = root / '06_validation/paper_aligned_recovery_characterization_20261012'
audit = root / '06_validation/all_states_lineage_audit_20261012'
master = pd.read_csv(old / 'CELL_STATE_ANNOTATION_MASTER.tsv', sep='\t').fillna('')
discovery = pd.read_csv(base / '00_audit/prior_discovery_associated_genes.tsv', sep='\t')
ranking = pd.read_csv(base / '02_external_markers/marker_composite_ranking.tsv', sep='\t')
expr = pd.read_csv(audit / 'tables/all_state_gene_expression.tsv.gz', sep='\t')
paired = pd.read_csv(audit / 'tables/all_associated_genes_paired_donor_replication.tsv.gz', sep='\t')
de = pd.read_csv(old / 'CELL_STATE_MARKER_EVIDENCE.tsv', sep='\t', low_memory=False)
de = de[de.EvidenceType.eq('Existing_external_DE') & de.Scope.eq('Dataset')].drop_duplicates(['Dataset', 'CellState', 'Gene'])
model = pd.read_csv(base / '00_audit/NMF_selected_model_manifest.tsv', sep='\t')
backgrounds = pd.read_csv(old / 'audit/ORA_GENE_BACKGROUNDS.tsv', sep='\t')
print('Read existing summaries:', len(expr), 'expression;', len(paired), 'paired-gene rows', flush=True)

# Step2. 完整发表表与完整原文程序列表分开登记，短诊断panel不进入统计主分析
signature_parts = []
f = pd.read_excel(out / 'literature_resources/41591_2021_1520_MOESM8_ESM.xlsx', header=1, engine='openpyxl')
f = f.rename(columns={'feature': 'Gene', 'cell_type': 'PublishedState', 'log2FoldChange': 'PublishedEffect', 'padj': 'PublishedAdjustedP'})
f['Publication'] = 'Friedrich2021'
f['DOI'] = '10.1038/s41591-021-01520-5'
f['SignatureScope'] = 'Complete Table6 exported significant upregulated cluster-marker table; not all tested genes'
f['SourceFile'] = '41591_2021_1520_MOESM8_ESM.xlsx'
f['Direction'] = 'up'
f['ExpectedParent'] = f.PublishedState.map({'Endothelial_cells': 'Endothelial', 'Pericytes': 'Pericyte', 'ABCA8_fibroblasts': 'Fibroblast', 'PDGFRA_fibroblasts': 'Fibroblast', 'Myofibroblasts': 'Fibroblast', 'Inflammatory_fibroblasts': 'Fibroblast'})
f = f[f.PublishedEffect.gt(0) & f.PublishedAdjustedP.lt(.05)].copy()
signature_parts.append(f)
il1 = pd.read_excel(out / 'literature_resources/41591_2021_1520_MOESM11_ESM.xlsx', header=1, engine='openpyxl')
il1 = il1.rename(columns={'gene_name': 'Gene', 'log2FoldChange': 'PublishedEffect', 'padj': 'PublishedAdjustedP'})
il1 = il1[il1.PublishedEffect.gt(0) & il1.PublishedAdjustedP.lt(.05)].copy()
il1['PublishedState'] = 'IL1beta_response_3h_fibroblasts'
il1['ExpectedParent'] = 'Fibroblast'
il1['Publication'] = 'Friedrich2021'
il1['DOI'] = '10.1038/s41591-021-01520-5'
il1['SignatureScope'] = 'Complete Table9 DE results, positive BH<0.05 subset; stimulation program, not cell identity'
il1['SourceFile'] = '41591_2021_1520_MOESM11_ESM.xlsx'
il1['Direction'] = 'up'
signature_parts.append(il1)
elm = pd.read_csv(out / 'literature_resources/Elmentaite2021_sheet1_raw.tsv', sep='\t')
elm = elm.rename(columns={'genes': 'Gene', 'cluster': 'PublishedState', 'logfoldchanges': 'PublishedEffect', 'pvals_adj': 'PublishedAdjustedP'})
elm['ExpectedParent'] = ''
for index, row in elm.iterrows():
    label = row.PublishedState
    parent = ''
    if 'Goblet' in label:
        parent = 'Goblet'
    elif label in ['D cells (SST+)', 'EC cells (NPW+)', 'EC cells (TAC1+)', 'I cells (CCK+)', 'K cells (GIP+)', 'L cells (PYY+)', 'M/X cells (MLN/GHRL+)', 'N cells (NTS+)']:
        parent = 'Enteroendocrine'
    elif label in ['Activated CD4 T', 'Tfh', 'Th1', 'Th17', 'Treg']:
        parent = 'CD4_T'
    elif 'CD8' in label:
        parent = 'CD8_T'
    elif 'gdT' in label or 'MAIT' in label:
        parent = 'T_gamma_delta' if 'gdT' in label else 'MAIT_T'
    elif 'ILC' in label or 'LTi-like' in label:
        parent = 'ILC'
    elif label in ['Cycling B cell', 'DZ GC cell', 'FCRL4+ Memory B', 'GC B cell', 'Immature B', 'LZ GC cell', 'Memory B', 'Naive B', 'Pre-B', 'Pro-B', 'STAT1+ Naive B']:
        parent = 'B'
    elif 'Mast cell' in label:
        parent = 'Mast'
    elif 'Macrophage' in label or 'macrophage' in label or label == 'Monocytes':
        parent = 'Monocyte_Macrophage'
    elif label in ['cDC1', 'cDC2', 'pDC', 'Lymphoid DC']:
        parent = 'DC'
    elif 'pericyte' in label.lower():
        parent = 'Pericyte'
    elif label.startswith(('Stromal ', 'Transitional Stromal ', 'myofibroblast')) or label in ['FDC', 'T reticular', 'mLTo', 'mLN Stroma (FMO2+)', 'cycling stromal']:
        parent = 'Fibroblast'
    elif label.startswith(('LEC', 'Fetal arterial EC', 'Fetal venous EC', 'Mature arterial EC', 'Mature venous EC')) or label in ['arterial capillary', 'cycling EC']:
        parent = 'Endothelial'
    elif 'glia' in label.lower():
        parent = 'Glial'
    elm.loc[index, 'ExpectedParent'] = parent
elm['Publication'] = 'Elmentaite2021'
elm['DOI'] = '10.1038/s41586-021-03852-1'
elm['SignatureScope'] = 'Complete published Table7 export of up to50 markers per cluster; mixed developmental/healthy atlas, not exhaustive DE'
elm['SourceFile'] = '41586_2021_3852_MOESM9_ESM.xls'
elm['Direction'] = 'up'
elm = elm[elm.PublishedEffect.gt(0) & elm.PublishedAdjustedP.lt(.05)].copy()
signature_parts.append(elm)
prior_panels = pd.read_csv(base / '00_audit/prior_literature_marker_panels.tsv', sep='\t')
martin = prior_panels[prior_panels.Publication.eq('Martin2019') & prior_panels.PanelStatus.eq('正文该程序完整基因列表')].copy()
martin = martin.rename(columns={'Program': 'PublishedState', 'ExpectedLineage': 'ExpectedParent'})
martin['Publication'] = 'Martin2019'
martin['DOI'] = '10.1016/j.cell.2019.08.008'
martin['SignatureScope'] = 'Complete explicitly enumerated main-text gene program; not complete cell-cluster DE table'
martin['SourceFile'] = 'PMC7060942_bioc.xml; passage61'
martin['Direction'] = 'positive_program_list'
signature_parts.append(martin)
signature = pd.concat(signature_parts, ignore_index=True, sort=False)
signature = signature.drop_duplicates(['Publication', 'PublishedState', 'Gene'])
signature.to_csv(out / 'literature_resources/PUBLISHED_SIGNATURE_GENES.tsv', sep='\t', index=False)
signature.groupby(['Publication', 'PublishedState', 'ExpectedParent', 'SignatureScope'], dropna=False).size().rename('PublishedGeneN').reset_index().to_csv(out / 'audit/PUBLISHED_SIGNATURE_DEFINITIONS.tsv', sep='\t', index=False)
print('Published signature rows:', len(signature), '; definitions:', signature.groupby(['Publication', 'PublishedState']).ngroups, flush=True)

# Step3. 同父类型全部state的NMF基因效应，区别唯一assignment与共享正效应
full_effect_parts = []
gene_info_sources = []
for parent in master.CellType.unique():
    rank = int(model.loc[model.CellType.eq(parent), 'SelectedRank'].iloc[0])
    folder = root / f'06_validation/GSE72819_recovery_review_20261009/raw_server/audit_snapshot/frozen_model/{parent}/{rank}'
    path = folder / 'gene_info.txt'
    gene_info_sources.append({'Path': str(path), 'SHA256': hashlib.sha256(path.read_bytes()).hexdigest()})
    with path.open(encoding='utf-8') as handle:
        reader = csv.reader(handle, delimiter='\t')
        headers = next(reader)
        records = list(reader)
    # R write.table保留了额外行名列，按原列名对齐，不能错移IS效应。
    info = pd.DataFrame([row[1:] if len(row) == len(headers) + 1 else row for row in records], columns=headers)
    mapping = info[['State', 'InitialState']].drop_duplicates()
    for state in master.loc[master.CellType.eq(parent), 'CellState']:
        state_short = state.removeprefix(parent + '_')
        initial = mapping.loc[mapping.State.eq(state_short), 'InitialState'].iloc[0]
        part = pd.DataFrame({'CellType': parent, 'CellState': state, 'Gene': info.Gene, 'DiscoveryScaledMeanDifference': pd.to_numeric(info[initial]), 'FrozenInitialState': initial})
        full_effect_parts.append(part)
full_effect = pd.concat(full_effect_parts, ignore_index=True)
pd.DataFrame(gene_info_sources).to_csv(out / 'audit/GENE_INFO_SOURCE_HASHES.tsv', sep='\t', index=False)
full_effect['PositiveWithinParentStates'] = full_effect.groupby(['CellType', 'Gene']).DiscoveryScaledMeanDifference.transform(lambda x: (x > .1).sum())
full_effect.to_csv(out / 'audit/FROZEN_WITHIN_PARENT_GENE_EFFECTS.tsv', sep='\t', index=False)

# Step4. 以原NMF候选背景做单侧超几何匹配，三种marker输入分别给统计量
match_rows = []
for state in master.itertuples():
    universe = set(backgrounds.loc[backgrounds.CellType.eq(state.CellType), 'Gene'])
    queries = {'DiscoveryTop50': set(discovery.loc[discovery.CellState.eq(state.CellState) & discovery.WithinStateRank.le(50), 'Gene'])}
    for variant in ['DiseaseGroup', 'Region3']:
        queries['ExternalTop50_' + variant] = set(ranking.loc[ranking.CellState.eq(state.CellState) & ranking.HeldOut.eq('All') & ranking.Variant.eq(variant) & ranking.WithinStateOrder.le(50), 'Gene'])
    for (publication, published), table in signature.groupby(['Publication', 'PublishedState'], sort=False):
        expected = table.ExpectedParent.iloc[0]
        same = state.CellType in str(expected).split(';') or expected == 'Any_shared_program'
        alternative = (state.CellType == 'Monocyte_Macrophage' and expected == 'DC') or (state.CellType in ['CD4_T', 'CD8_T'] and expected in ['T_gamma_delta', 'MAIT_T'])
        if not same and not alternative:
            continue
        target = set(table.Gene) & universe
        for query_name, query_genes in queries.items():
            selected = query_genes & universe
            N, n, K = len(universe), len(selected), len(target)
            overlap = selected & target
            k = len(overlap)
            p = np.nan
            evaluable = n >= 5 and K >= 5
            if evaluable:
                denominator = math.lgamma(N + 1) - math.lgamma(n + 1) - math.lgamma(N - n + 1)
                terms = []
                for hit_n in range(k, min(K, n) + 1):
                    if 0 <= n - hit_n <= N - K:
                        log_value = math.lgamma(K + 1) - math.lgamma(hit_n + 1) - math.lgamma(K - hit_n + 1)
                        log_value += math.lgamma(N - K + 1) - math.lgamma(n - hit_n + 1) - math.lgamma(N - K - n + hit_n + 1) - denominator
                        terms.append(math.exp(log_value))
                p = min(1, math.fsum(terms))
            match_rows.append({'CellState': state.CellState, 'CellType': state.CellType, 'Publication': publication, 'DOI': table.DOI.iloc[0], 'PublishedStateOrProgram': published, 'PublishedExpectedParent': expected, 'ComparisonScope': 'same_parent_or_shared_program' if same else 'explicit_alternative_lineage_audit', 'Query': query_name, 'BackgroundN': N, 'QueryInBackgroundN': n, 'PublishedSignatureN': table.Gene.nunique(), 'SignatureInBackgroundN': K, 'OverlapN': k, 'OverlapGenes': ';'.join(sorted(overlap)), 'JaccardInBackground': k / (n + K - k) if n + K - k else np.nan, 'QueryCoverage': k / n if n else np.nan, 'SignatureCoverageInBackground': k / K if K else np.nan, 'P': p, 'Evaluable': evaluable, 'EffectDirection': 'positive markers/program list', 'SignatureSourceScope': table.SignatureScope.iloc[0], 'EvidenceIndependence': 'ExternalTop50 reuse marker-selection cohorts; discovery not independent verification', 'Limit': 'Conditional NMF feature-space overlap, not identity accuracy or novelty test'})
matches = pd.DataFrame(match_rows)
matches['BH_Q_GlobalEvaluableMatches'] = np.nan
evaluable = matches[matches.Evaluable].sort_values('P')
q = evaluable.P.to_numpy() * len(evaluable) / np.arange(1, len(evaluable) + 1)
matches.loc[evaluable.index, 'BH_Q_GlobalEvaluableMatches'] = np.minimum(1, np.minimum.accumulate(q[::-1])[::-1])
matches.to_csv(out / 'audit/ALL_PUBLISHED_SIGNATURE_MATCH_STATISTICS.tsv', sep='\t', index=False)
print('Signature comparisons:', len(matches), '; evaluable:', len(evaluable), flush=True)

# Step5. 每个Top50基因记录同父类型的共享/区分和供者复现；不把assignment独占当成特异性
top = discovery[discovery.WithinStateRank.le(50)][['CellState', 'Gene', 'WithinStateRank', 'MaxFC']].merge(master[['CellState', 'CellType', 'Ecotype']], on='CellState', validate='many_to_one')
within = top.merge(full_effect, on=['CellType', 'CellState', 'Gene'], validate='one_to_one')
state_n = master.groupby('CellType').size()
within['ParentFrozenStateN'] = within.CellType.map(state_n)
for index, gene in within.iterrows():
    others = full_effect[full_effect.CellType.eq(gene.CellType) & full_effect.Gene.eq(gene.Gene) & full_effect.CellState.ne(gene.CellState)]
    within.loc[index, 'DiscoveryBestOtherStateEffect'] = others.DiscoveryScaledMeanDifference.max()
    within.loc[index, 'DiscoveryEffectGapToBestOther'] = gene.DiscoveryScaledMeanDifference - others.DiscoveryScaledMeanDifference.max()
    within.loc[index, 'OtherStatesWithPositiveEffectGt0p1'] = ';'.join(others.loc[others.DiscoveryScaledMeanDifference.gt(.1), 'CellState'])
for dataset in ['Colon_GSE214695', 'Ileum_IBDverse']:
    parent_expression = expr[expr.Dataset.eq(dataset) & expr.CellState.isin(master.CellState)].merge(master[['CellState', 'CellType']], on='CellState', validate='many_to_one')
    prefix = 'Colon' if dataset.startswith('Colon') else 'Ileum'
    for index, gene in within.iterrows():
        rows = parent_expression[parent_expression.CellType.eq(gene.CellType) & parent_expression.Gene.eq(gene.Gene)]
        own = rows[rows.CellState.eq(gene.CellState)]
        other = rows[rows.CellState.ne(gene.CellState) & rows.AssignedCells.gt(0) & rows.SourceMeasured.eq(True)]
        if len(own) and own.SourceMeasured.iloc[0] and own.AssignedCells.iloc[0] > 0:
            within.loc[index, prefix + '_DetectionFraction'] = own.DetectionFraction.iloc[0]
            within.loc[index, prefix + '_MeanLog2CPM'] = own.MeanLog2CPM.iloc[0]
            within.loc[index, prefix + '_AssignedCells'] = own.AssignedCells.iloc[0]
            within.loc[index, prefix + '_AssignedDonors'] = own.AssignedDonors.iloc[0]
            within.loc[index, prefix + '_MeanLog2CPMGapToBestOther'] = own.MeanLog2CPM.iloc[0] - other.MeanLog2CPM.max()
            within.loc[index, prefix + '_DetectionGapToBestOther'] = own.DetectionFraction.iloc[0] - other.DetectionFraction.max()
        rep = paired[paired.Dataset.eq(dataset) & paired.CellState.eq(gene.CellState) & paired.Gene.eq(gene.Gene)]
        if len(rep):
            r = rep.iloc[0]
            within.loc[index, prefix + '_MatchedDonors'] = r.MatchedDonors_Min5CellsEach
            within.loc[index, prefix + '_PairedMeanDifference'] = r.MeanPairedLog2CPM_Difference
            within.loc[index, prefix + '_FractionDonorsPositive'] = r.FractionMatchedDonors_Positive
            within.loc[index, prefix + '_DonorDirectionReplicated'] = r.MatchedDonors_Min5CellsEach >= 3 and r.MeanPairedLog2CPM_Difference > 0 and r.FractionMatchedDonors_Positive >= .7
        d = de[de.Dataset.eq(dataset) & de.CellState.eq(gene.CellState) & de.Gene.eq(gene.Gene)]
        if len(d):
            within.loc[index, prefix + '_ExistingDE_log2FC'] = d.log2FC.iloc[0]
            within.loc[index, prefix + '_ExistingDE_BonferroniP'] = d.p_val_adj.iloc[0]
within['BothDatasetDonorDirectionReplicated'] = within.Colon_DonorDirectionReplicated.fillna(False) & within.Ileum_DonorDirectionReplicated.fillna(False)
within['EvidenceScope'] = 'Same pooled recovery/marker-selection cohorts; descriptive donor checks, no new independent test'
within.to_csv(out / 'WITHIN_LINEAGE_STATE_COMPARISON.tsv', sep='\t', index=False)
state_summary = within.groupby(['CellState', 'CellType', 'Ecotype']).agg(Top50N=('Gene', 'size'), MedianDiscoveryEffectGap=('DiscoveryEffectGapToBestOther', 'median'), SharedPositiveGeneN=('PositiveWithinParentStates', lambda x: (x > 1).sum()), BothDatasetStableGeneN=('BothDatasetDonorDirectionReplicated', 'sum')).reset_index()
state_summary['BothDatasetStableGenes'] = state_summary.CellState.map(within[within.BothDatasetDonorDirectionReplicated].groupby('CellState').Gene.agg(';'.join)).fillna('')
state_summary.to_csv(out / 'audit/ALL53_WITHIN_LINEAGE_SUMMARY.tsv', sep='\t', index=False)
print(state_summary.sort_values('BothDatasetStableGeneN', ascending=False).head(20).to_string(index=False), flush=True)

# Step6. 仅用同群体边际检出率计算共检出界限，不能冒充已观测联合细胞计数
core_panels = [('Monocyte_Macrophage_S02', 'P2RY13+CLEC10A+CSF1R', ['P2RY13', 'CLEC10A', 'CSF1R']), ('Monocyte_Macrophage_S02', 'PLA2G2D+P2RY13', ['PLA2G2D', 'P2RY13']), ('Monocyte_Macrophage_S02', 'FLT3+CSF1R', ['FLT3', 'CSF1R']), ('Goblet_S01', 'MUC2+MKI67', ['MUC2', 'MKI67']), ('Goblet_S01', 'MUC2+MKI67+TOP2A', ['MUC2', 'MKI67', 'TOP2A']), ('CD8_T_S03', 'CD8A+IL23R', ['CD8A', 'IL23R']), ('CD8_T_S03', 'CD8A+JUND+LAPTM5', ['CD8A', 'JUND', 'LAPTM5'])]
bounds = []
for state, label, genes in core_panels:
    for dataset in ['Colon_GSE214695', 'Ileum_IBDverse']:
        rows = expr[expr.Dataset.eq(dataset) & expr.CellState.eq(state) & expr.Gene.isin(genes)]
        evaluable = len(rows) == len(genes) and rows.SourceMeasured.all() and rows.AssignedCells.gt(0).all()
        bounds.append({'CellState': state, 'Dataset': dataset, 'Panel': label, 'MarginalDetectionFractions': ';'.join(f'{r.Gene}:{r.DetectionFraction:.6g}' for r in rows.itertuples()), 'LowerBoundDonorBalancedJointDetection': max(0, rows.DetectionFraction.sum() - (len(genes) - 1)) if evaluable else np.nan, 'UpperBoundDonorBalancedJointDetection': rows.DetectionFraction.min() if evaluable else np.nan, 'Method': 'Frechet bounds from same-group donor-balanced marginal detection fractions', 'ObservedSameCellCoexpression': 'NOT_MEASURED_LOCAL_CHECKPOINT_MISSING', 'Limit': 'Bounds only; no observed joint counts, donor-specific joint positivity, multi-gene score or independent validation'})
pd.DataFrame(bounds).to_csv(out / 'audit/CANDIDATE_COEXPRESSION_BOUNDS.tsv', sep='\t', index=False)

# Step7. 逐谱系全部state两两比较可测表达，保留负差、缺测及不独立属性
pair_rows = []
for parent in master.CellType.unique():
    states = master.loc[master.CellType.eq(parent), 'CellState'].tolist()
    for a_index, state_a in enumerate(states):
        genes_a = set(top.loc[top.CellState.eq(state_a), 'Gene'])
        for state_b in states[a_index + 1:]:
            genes_b = set(top.loc[top.CellState.eq(state_b), 'Gene'])
            union = genes_a | genes_b
            for dataset in ['Colon_GSE214695', 'Ileum_IBDverse']:
                a = expr[expr.Dataset.eq(dataset) & expr.CellState.eq(state_a) & expr.Gene.isin(union)].set_index('Gene')
                b = expr[expr.Dataset.eq(dataset) & expr.CellState.eq(state_b) & expr.Gene.isin(union)].set_index('Gene')
                common = a.index.intersection(b.index)
                delta = a.loc[common, 'MeanLog2CPM'] - b.loc[common, 'MeanLog2CPM']
                status = 'DESCRIPTIVE_POOLED_STATE_MEANS' if len(a) and len(b) and a.AssignedCells.iloc[0] >= 5 and b.AssignedCells.iloc[0] >= 5 else 'LOW_OR_ABSENT_STATE_SUPPORT'
                pair_rows.append({'CellType': parent, 'StateA': state_a, 'StateB': state_b, 'Dataset': dataset, 'Top50SharedN': len(genes_a & genes_b), 'Top50UnionN': len(union), 'CommonMeasuredGeneN': len(common), 'A_FavoringGenes_MeanLog2CPMDeltaGt0p5': ';'.join(delta[delta > .5].sort_values(ascending=False).index[:20]), 'B_FavoringGenes_MeanLog2CPMDeltaLtMinus0p5': ';'.join(delta[delta < -.5].sort_values().index[:20]), 'Status': status, 'Limit': 'Donor-balanced pooled means; not a donor-matched pairwise DE test; threshold0.5 is descriptive only'})
pd.DataFrame(pair_rows).to_csv(out / 'audit/WITHIN_LINEAGE_PAIRWISE_DESCRIPTIVE.tsv', sep='\t', index=False)
best = matches[matches.Query.eq('DiscoveryTop50') & matches.ComparisonScope.eq('same_parent_or_shared_program') & matches.Evaluable].sort_values(['BH_Q_GlobalEvaluableMatches', 'OverlapN'], ascending=[True, False]).groupby('CellState').head(3)
best.to_csv(out / 'audit/BEST_DISCOVERY_PUBLISHED_MATCHES.tsv', sep='\t', index=False)
print(best[['CellState', 'Publication', 'PublishedStateOrProgram', 'OverlapN', 'BH_Q_GlobalEvaluableMatches']].to_string(index=False), flush=True)
