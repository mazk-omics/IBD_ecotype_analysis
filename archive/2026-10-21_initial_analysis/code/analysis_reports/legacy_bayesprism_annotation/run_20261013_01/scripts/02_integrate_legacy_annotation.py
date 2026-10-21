# 整合冻结53状态和7生态型，所有输出仅写入本次独立目录。
# 单线程流式读既有DE，不重跑BayesPrism、NMF、recovery或差异分析。
import os
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['NUMEXPR_NUM_THREADS'] = '1'
from pathlib import Path
import csv
import hashlib
import json
import math
import time
import pandas as pd
import numpy as np

# Step1. 固定只读来源与状态身份
started = time.time()
root = Path('D:/IBD_ecotype/IBD_ecotype_analysis')
out = Path(__file__).resolve().parents[1]
base = root / '06_validation/paper_aligned_recovery_characterization_20261012'
lineage = root / '06_validation/all_states_lineage_audit_20261012'
source_paths = [base / '00_audit/frozen_states.tsv', base / '00_audit/prior_discovery_associated_genes.tsv', base / '00_audit/NMF_selected_model_manifest.tsv', base / '03_state_annotation/53_states_annotation_evidence.tsv', lineage / 'tables/all_53_states_final_evidence_review.tsv', base / '02_external_markers/marker_composite_ranking.tsv', base / '02_external_markers/marker_n1_n7.tsv', base / '02_external_markers/marker_LOOCV_per_dataset.tsv', base / '02_external_markers/n4_disease_vs_region_comparison.tsv', base / '02_external_markers/DE_stratum_eligibility_parallel4.tsv', base / '01_state_recovery/state_recovery_evidence_by_dataset.tsv', base / '00_audit/Colon_GSE214695_cells.tsv', base / '00_audit/Ileum_IBDverse_cells.tsv', base / '00_audit/h.all.v7.4.symbols.gmt', base / '04_state_phenotypes/state_phenotype_enrichment.tsv', base / '04_state_phenotypes/state_region_fisher.tsv', base / '06_ecotype_features/7_ecotypes_evidence_summary.tsv', base / '06_ecotype_features/ecotype_discrete_features.tsv', base / '06_ecotype_features/ecotype_continuous_features.tsv', base / '06_ecotype_features/ecotype_region_KW.tsv', base / '00_audit/prior_literature_marker_panels.tsv', lineage / 'tables/state_author_subtype_composition.tsv']
manifest_rows = []
for path in source_paths:
    manifest_rows.append({'Path': str(path), 'Bytes': path.stat().st_size, 'SHA256': hashlib.sha256(path.read_bytes()).hexdigest()})
pd.DataFrame(manifest_rows).to_csv(out / 'audit/INPUT_SHA256_BEFORE.tsv', sep='\t', index=False)
frozen = pd.read_csv(source_paths[0], sep='\t')
assert len(frozen) == 53 and frozen.CellState.is_unique and frozen.Ecotype.nunique() == 7
discovery = pd.read_csv(source_paths[1], sep='\t')
prior = pd.read_csv(base / '03_state_annotation/53_states_annotation_evidence.tsv', sep='\t').fillna('')
review = pd.read_csv(lineage / 'tables/all_53_states_final_evidence_review.tsv', sep='\t').fillna('')
ranking = pd.read_csv(base / '02_external_markers/marker_composite_ranking.tsv', sep='\t')
loocv = pd.read_csv(base / '02_external_markers/marker_LOOCV_per_dataset.tsv', sep='\t')
recovery = pd.read_csv(base / '01_state_recovery/state_recovery_evidence_by_dataset.tsv', sep='\t')
n4 = pd.read_csv(base / '02_external_markers/n4_disease_vs_region_comparison.tsv', sep='\t')
n4.to_csv(out / 'audit/N4_EXISTING_COMPARISON.tsv', sep='\t', index=False)
print('Frozen identity:', len(frozen), 'states;', frozen.Ecotype.nunique(), 'ecotypes', flush=True)

# Step2. 保存Top50、两版排名和可直接核对的基因层证据
top_discovery = discovery[discovery.WithinStateRank <= 50].copy()
top_discovery['EvidenceType'] = 'Discovery_gene_info'
top_discovery['Rank'] = top_discovery.WithinStateRank
top_discovery['Variant'] = 'FrozenDiscovery'
top_discovery['HeldOut'] = 'NotApplicable'
top_discovery['EffectMetric'] = 'EcoTyper_gene_info_MaxFC; not external log2FC'
top_discovery['Direction'] = np.where(top_discovery.MaxFC > 0, 'positive_state_association', 'nonpositive')
top_rank = ranking[ranking.WithinStateOrder <= 50].copy()
top_rank['EvidenceType'] = 'External_composite_ranking'
top_rank['Rank'] = top_rank.WithinStateOrder
top_rank['EffectMetric'] = 'CompositeScore is rank aggregation; n7 is existing external FC summary'
top_rank['Direction'] = np.where(top_rank.n7 > 0, 'positive_external_mean', 'nonpositive_external_mean')
marker = pd.concat([top_discovery, top_rank], ignore_index=True, sort=False)
candidate_keys = set(zip(marker.CellState, marker.Gene))
marker_de_parts = []
de_top_parts = []
de_counts = []
de_path = base / '02_external_markers/external_DEG_by_dataset.tsv.gz'
de_stat_before = de_path.stat()
for chunk in pd.read_csv(de_path, sep='\t', chunksize=25000):
    chunk['CellState'] = chunk.CellType + '_' + chunk.State
    counts = chunk.groupby(['Dataset', 'CellState', 'Scope', 'Stratum']).size().rename('SavedDEGenes').reset_index()
    de_counts.append(counts)
    retain = [(state, gene) in candidate_keys for state, gene in zip(chunk.CellState, chunk.Gene)]
    marker_de_parts.append(chunk.loc[retain].copy())
    dataset_de = chunk[(chunk.Scope == 'Dataset') & (chunk.p_val_adj < 0.05) & (chunk.log2FC > 0)].copy()
    dataset_de = dataset_de.sort_values(['Dataset', 'CellState', 'p_val_adj', 'log2FC'], ascending=[True, True, True, False])
    de_top_parts.append(dataset_de.groupby(['Dataset', 'CellState'], sort=False).head(50))
de_evidence = pd.concat(marker_de_parts, ignore_index=True)
de_evidence['EvidenceType'] = 'Existing_external_DE'
de_evidence['EffectMetric'] = 'log2FC = Seurat3.1.3 avg_logFC / ln(2); pseudocount=1'
de_evidence['Direction'] = np.where(de_evidence.log2FC > 0, 'up', 'down')
de_evidence['AdjustedPMethod'] = 'Native_Seurat_Bonferroni; not BH_FDR'
de_top = pd.concat(de_top_parts, ignore_index=True).drop_duplicates(['Dataset', 'CellState', 'Gene'])
de_top = de_top.sort_values(['Dataset', 'CellState', 'p_val_adj', 'log2FC'], ascending=[True, True, True, False])
de_top = de_top.groupby(['Dataset', 'CellState'], sort=False).head(50).copy()
de_top['Rank'] = de_top.groupby(['Dataset', 'CellState']).cumcount() + 1
de_top['EvidenceType'] = 'External_transcriptome_DE_top50'
de_top['Direction'] = 'up'
de_top['AdjustedPMethod'] = 'Native_Seurat_Bonferroni; not BH_FDR'
marker = pd.concat([marker, de_evidence, de_top], ignore_index=True, sort=False)
marker['SourcePath'] = np.select([marker.EvidenceType.eq('Discovery_gene_info'), marker.EvidenceType.eq('External_composite_ranking')], [str(source_paths[1]), str(base / '02_external_markers/marker_composite_ranking.tsv')], default=str(de_path))
marker.to_csv(out / 'CELL_STATE_MARKER_EVIDENCE.tsv', sep='\t', index=False)
de_coverage = pd.concat(de_counts).groupby(['Dataset', 'CellState', 'Scope', 'Stratum'], as_index=False).SavedDEGenes.sum()
de_coverage.to_csv(out / 'audit/EXISTING_DE_COVERAGE.tsv', sep='\t', index=False)
print('Existing DE streamed once; marker evidence rows:', len(marker), flush=True)

# Step3. 用冻结NMF候选基因背景补轻量Hallmark ORA
# 初始gene_info包含后来过滤掉的初始状态；不能只用53状态关联基因作背景。
gene_sets = {}
with open(base / '00_audit/h.all.v7.4.symbols.gmt', encoding='utf-8') as handle:
    for line in handle:
        fields = line.rstrip('\n').split('\t')
        gene_sets[fields[0]] = set(fields[2:])
assert len(gene_sets) == 50
universes = {}
model_manifest = pd.read_csv(base / '00_audit/NMF_selected_model_manifest.tsv', sep='\t')
background_rows = []
for cell_type in frozen.CellType.unique():
    rank = int(model_manifest[model_manifest.CellType.eq(cell_type)].SelectedRank.iloc[0])
    path = root / f'06_validation/GSE72819_recovery_review_20261009/raw_server/audit_snapshot/frozen_model/{cell_type}/{rank}/initial_gene_info.txt'
    expected = model_manifest[model_manifest.CellType.eq(cell_type) & model_manifest.File.str.endswith('/initial_gene_info.txt')]
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    assert len(expected) == 1 and digest == expected.SHA256.iloc[0], f'Model snapshot mismatch: {cell_type}'
    with open(path, encoding='utf-8') as handle:
        reader = csv.reader(handle, delimiter='\t')
        next(reader)
        genes = {row[0] for row in reader}
    universes[cell_type] = genes
    background_rows.extend({'CellType': cell_type, 'Gene': gene, 'Background': 'Frozen_initial_gene_info_all_initial_states', 'Path': str(path), 'SHA256': digest} for gene in sorted(genes))
    manifest_rows.append({'Path': str(path), 'Bytes': path.stat().st_size, 'SHA256': digest})
pd.DataFrame(background_rows).to_csv(out / 'audit/ORA_GENE_BACKGROUNDS.tsv', sep='\t', index=False)
ora_rows = []
for state in frozen.itertuples():
    universe = universes[state.CellType]
    selected = set(discovery.loc[discovery.CellState.eq(state.CellState) & discovery.MaxFC.gt(0), 'Gene']) & universe
    N = len(universe)
    n = len(selected)
    for name, pathway in gene_sets.items():
        measured = pathway & universe
        hits = selected & measured
        K = len(measured)
        k = len(hits)
        p = 1.0
        if k > 0 and K > 0 and n > 0:
            log_denominator = math.lgamma(N + 1) - math.lgamma(n + 1) - math.lgamma(N - n + 1)
            terms = []
            for x in range(k, min(K, n) + 1):
                if 0 <= n - x <= N - K:
                    log_numerator = math.lgamma(K + 1) - math.lgamma(x + 1) - math.lgamma(K - x + 1)
                    log_numerator += math.lgamma(N - K + 1) - math.lgamma(n - x + 1) - math.lgamma(N - K - n + x + 1)
                    terms.append(math.exp(log_numerator - log_denominator))
            p = min(1.0, math.fsum(terms))
        ora_rows.append({'CellState': state.CellState, 'CellType': state.CellType, 'Ecotype': state.Ecotype, 'GeneSet': name, 'Database': 'MSigDB Hallmark v7.4 frozen GMT', 'Method': 'one-sided hypergeometric ORA', 'Direction': 'positive discovery-associated genes; not activation score', 'BackgroundN': N, 'SelectedN': n, 'SetInBackgroundN': K, 'SetOriginalN': len(pathway), 'OverlapN': k, 'HitGenes': ';'.join(sorted(hits)), 'P': p, 'OddsRatio': (k * (N - K - n + k) / ((n - k) * (K - k))) if (n - k) * (K - k) > 0 else np.nan, 'Limit': 'Conditional on NMF feature selection; no claim of whole-transcriptome unbiased enrichment'})
enrichment = pd.DataFrame(ora_rows)
order = np.argsort(enrichment.P.to_numpy(), kind='stable')
q_sorted = enrichment.P.to_numpy()[order] * len(enrichment) / np.arange(1, len(enrichment) + 1)
q_sorted = np.minimum(1, np.minimum.accumulate(q_sorted[::-1])[::-1])
enrichment['FDR_Global2650'] = np.nan
enrichment.loc[order, 'FDR_Global2650'] = q_sorted
enrichment['FDR_WithinState50'] = np.nan
for state, index in enrichment.groupby('CellState', sort=False).groups.items():
    index = list(index)
    local = enrichment.loc[index].sort_values('P').index
    q = enrichment.loc[local, 'P'].to_numpy() * 50 / np.arange(1, 51)
    enrichment.loc[local, 'FDR_WithinState50'] = np.minimum(1, np.minimum.accumulate(q[::-1])[::-1])
enrichment['PrimaryFDRFamily'] = 'all 53 states x 50 Hallmarks = 2650 tests'
enrichment.to_csv(out / 'CELL_STATE_FUNCTIONAL_ENRICHMENT.tsv', sep='\t', index=False)
print('Hallmark ORA:', len(enrichment), 'tests;', int((enrichment.FDR_Global2650 < .05).sum()), 'global FDR<0.05', flush=True)

# Step4. 联合疾病-部位仅汇总现有分配，不替代历史n4
cells_parts = []
for dataset in ['Colon_GSE214695', 'Ileum_IBDverse']:
    cells = pd.read_csv(base / f'00_audit/{dataset}_cells.tsv', sep='\t')
    cells['Dataset'] = dataset
    cells_parts.append(cells)
cells = pd.concat(cells_parts, ignore_index=True, sort=False)
cells['Context'] = cells.DiseaseGroup + '_' + cells.Region3
cells[['Dataset', 'Context', 'DiseaseGroup', 'Region3', 'Donor', 'Sample']].drop_duplicates().groupby(['Dataset', 'Context', 'DiseaseGroup', 'Region3']).agg(Donors=('Donor', 'nunique'), Samples=('Sample', 'nunique')).reset_index().to_csv(out / 'audit/CONTEXT_DONOR_METADATA.tsv', sep='\t', index=False)
external_rows = []
for state in frozen.itertuples():
    state_review = review[review.CellState.eq(state.CellState)].iloc[0]
    for dataset, context in [('Ileum_IBDverse', 'CD_Ileum'), ('Colon_GSE214695', 'CD_Colon'), ('Colon_GSE214695', 'UC_Colon')]:
        parent = cells[cells.Dataset.eq(dataset) & cells.Context.eq(context) & cells.CellType.eq(state.CellType)]
        own = parent[parent.CellState.eq(state.CellState)]
        donor_counts = own.groupby('Donor').size()
        assigned_counts = parent[parent.State.ne('Unassigned')].groupby('CellState').size()
        eligible_states = assigned_counts[assigned_counts >= 5]
        contrast_evaluable = len(own) >= 5 and len(eligible_states) >= 2
        rec = recovery[recovery.Dataset.eq(dataset) & recovery.CellState.eq(state.CellState)].iloc[0]
        prefix = 'Colon' if dataset.startswith('Colon') else 'Ileum'
        native_scope = 'DiseaseGroup'
        native_stratum = context.split('_')[0]
        mixed_regions = sorted(cells.loc[cells.Dataset.eq(dataset) & cells.DiseaseGroup.eq(native_stratum), 'Region3'].dropna().unique())
        exact_scope = len(mixed_regions) == 1 and mixed_regions[0] == context.split('_')[1]
        state_de = de_evidence[de_evidence.Dataset.eq(dataset) & de_evidence.CellState.eq(state.CellState) & de_evidence.Scope.eq(native_scope) & de_evidence.Stratum.eq(native_stratum)]
        top20 = set(discovery.loc[discovery.CellState.eq(state.CellState) & discovery.WithinStateRank.le(20), 'Gene'])
        support_genes = sorted(set(state_de.loc[state_de.Gene.isin(top20) & state_de.log2FC.gt(0) & state_de.p_val_adj.lt(.05), 'Gene']))
        folds = loocv[loocv.HeldOut.eq(dataset) & loocv.CellState.eq(state.CellState)]
        row = {'CellState': state.CellState, 'CellType': state.CellType, 'Ecotype': state.Ecotype, 'Dataset': dataset, 'Context': context, 'InputCellsInParentContext': len(parent), 'InputDonorsInParentContext': parent.Donor.nunique(), 'AssignedCellsInContext': len(own), 'AssignedDonorsInContext': own.Donor.nunique(), 'AssignedSamplesInContext': own.Sample.nunique(), 'DonorsWithGE5AssignedCells': int((donor_counts >= 5).sum()), 'ContextStatus': 'NO_INPUT_LINEAGE' if len(parent) == 0 else ('NO_ASSIGNED_CELLS' if len(own) == 0 else ('DESCRIPTIVE_GE5_CELLS' if len(own) >= 5 else 'LOW_CELL_SUPPORT')), 'JointContextDEEvaluabilityFromExistingCells': contrast_evaluable, 'JointContextEligibleStates': ';'.join(eligible_states.index), 'RecoveryOfficialZ_PooledDataset': rec.OfficialZ, 'RecoveryZgt1p96_PooledDataset': rec.PositiveZgt1p96, 'RecoveryInferenceScope': 'Pooled dataset only; no context-specific permutation test', 'MarkerDEAvailableScope': native_scope + ':' + native_stratum, 'MarkerDEMatchesJointContext': exact_scope, 'MarkerDESourceRegions': ';'.join(mixed_regions), 'Top20PositiveBonferroniGenes_ExistingDiseaseScope': ';'.join(support_genes), 'JointContextMarkerStatus': ('SUPPORTED_GE2_GENES_EXACT_SCOPE' if len(support_genes) >= 2 else 'FEWER_THAN2_SIGNIFICANT_TOP20_GENES') if exact_scope and contrast_evaluable else ('NOT_EVALUABLE_EXISTING_CONTRAST' if not contrast_evaluable else 'POOLED_COLON_RECTUM_PROXY_ONLY'), 'PooledLineageAnchorsDetected': state_review[prefix + '_DetectedAnchorGenes_GE10pct'], 'PooledMatchedDonorsMin5CellsEach': state_review[prefix + '_MatchedDonors_Min5CellsEach'], 'PooledPairedTop20PositiveGenes': state_review[prefix + '_PairedTop20PositiveGenes'], 'PooledPairedStatus': state_review[prefix + '_PairedStatus'], 'LOOCV_DiseaseGroup_Status': folds.loc[folds.Variant.eq('DiseaseGroup'), 'Status'].iloc[0], 'LOOCV_DiseaseGroup_Concordant': folds.loc[folds.Variant.eq('DiseaseGroup'), 'Concordant'].iloc[0], 'LOOCV_Region3_Status': folds.loc[folds.Variant.eq('Region3'), 'Status'].iloc[0], 'LOOCV_Region3_Concordant': folds.loc[folds.Variant.eq('Region3'), 'Concordant'].iloc[0], 'Limit': 'Frozen sampled cells; marker-dependent assignment; pooled evidence cannot be relabeled joint-context validation'}
        external_rows.append(row)
external = pd.DataFrame(external_rows)
external.to_csv(out / 'CELL_STATE_EXTERNAL_SUPPORT.tsv', sep='\t', index=False)
print('External context support:', len(external), 'state-context rows', flush=True)

# Step5. 核实过的文献短面板；命中是候选程序对应，不是已发表亚型的身份证明
publications = [
    ('Luca2021', 'Atlas of clinically distinct cell states and ecosystems across human solid tumors', '10.1016/j.cell.2021.09.014', 'https://pmc.ncbi.nlm.nih.gov/articles/PMC8526411/', 'Methods and Table S4 link verified; supplement binary access blocked by browser challenge'),
    ('Martin2019', 'Single-Cell Analysis of Crohn’s Disease Lesions Identifies a Pathogenic Cellular Module Associated with Resistance to Anti-TNF Therapy', '10.1016/j.cell.2019.08.008', 'https://pmc.ncbi.nlm.nih.gov/articles/PMC7060942/', 'Bibliography verified online; detailed gene panels reused from prior passage-indexed extraction'),
    ('Friedrich2021', 'IL-1-driven stromal–neutrophil interactions define a subset of patients with inflammatory bowel disease that does not respond to therapies', '10.1038/s41591-021-01520-5', 'https://www.nature.com/articles/s41591-021-01520-5', 'Primary full-text indexed paragraphs and PubMed verified online'),
    ('Smillie2019', 'Intra- and Inter-cellular Rewiring of the Human Colon during Ulcerative Colitis', '10.1016/j.cell.2019.06.029', 'https://pubmed.ncbi.nlm.nih.gov/31348891/', 'Primary abstract verified; no full supplementary-signature equivalence claimed'),
    ('Elmentaite2021', 'Cells of the human intestinal tract mapped across space and time', '10.1038/s41586-021-03852-1', 'https://www.nature.com/articles/s41586-021-03852-1', 'Primary atlas and lymphoid-organizer paragraphs verified online'),
    ('GarridoTrigo2023', 'Macrophage and neutrophil heterogeneity at single-cell spatial resolution in human inflammatory bowel disease', '10.1038/s41467-023-40156-6', 'https://pmc.ncbi.nlm.nih.gov/articles/PMC10372067/', 'Primary article DOI and cohort confirmed online')]
sources = pd.DataFrame(publications, columns=['Publication', 'Title', 'DOI', 'URL', 'Verification'])
sources['VerifiedDate'] = '2026-10-13'
sources.to_csv(out / 'audit/VERIFIED_LITERATURE_SOURCES.tsv', sep='\t', index=False)
panels = pd.read_csv(base / '00_audit/prior_literature_marker_panels.tsv', sep='\t').fillna('')
panels = panels[['Program', 'Gene', 'ExpectedLineage', 'Publication', 'PanelStatus', 'BioCPassageIndex0']].copy()
additional_panels = [('Friedrich2021', 'Inflammatory_fibroblast_identity', 'Fibroblast', 'PDPN;THY1;FAP', 'Primary Fig4c; short selected panel'), ('Friedrich2021', 'Neutrophil_attracting_stromal_program', 'Fibroblast;Pericyte;Endothelial', 'CXCL1;CXCL2;CXCL3;CXCL5;CXCL6;CXCL8', 'Primary Fig5b; complete chemokine list in indexed paragraph'), ('Friedrich2021', 'Pericyte_identity', 'Pericyte', 'NOTCH3;MCAM', 'Primary Fig4c; short identity panel'), ('Friedrich2021', 'Endothelial_identity', 'Endothelial', 'ACKR1;CD34;PECAM1', 'Primary Fig4c and Extended Data Fig7e; selected panel'), ('Smillie2019', 'Inflammatory_fibroblast', 'Fibroblast', 'IL13RA2;IL11', 'Primary abstract; short selected panel'), ('Elmentaite2021', 'Lymphoid_organizer_stromal', 'Fibroblast', 'CCL19;CCL21;CXCL13', 'Primary Extended Data Fig13b-f; short selected panel')]
new_panel_rows = []
for publication, program, expected, genes, status in additional_panels:
    new_panel_rows.extend({'Program': program, 'Gene': gene, 'ExpectedLineage': expected, 'Publication': publication, 'PanelStatus': status, 'BioCPassageIndex0': ''} for gene in genes.split(';'))
panels = pd.concat([panels, pd.DataFrame(new_panel_rows)], ignore_index=True).drop_duplicates(['Publication', 'Program', 'Gene'])
panels.to_csv(out / 'audit/LITERATURE_GENE_PANELS.tsv', sep='\t', index=False)
literature_rows = []
for state in frozen.itertuples():
    associated = set(discovery.loc[discovery.CellState.eq(state.CellState), 'Gene'])
    top20 = set(discovery.loc[discovery.CellState.eq(state.CellState) & discovery.WithinStateRank.le(20), 'Gene'])
    external_positive = set(de_evidence.loc[de_evidence.CellState.eq(state.CellState) & de_evidence.Scope.eq('Dataset') & de_evidence.log2FC.gt(0) & de_evidence.p_val_adj.lt(.05), 'Gene'])
    for (publication, program), panel in panels.groupby(['Publication', 'Program'], sort=False):
        genes = set(panel.Gene)
        hits = associated & genes
        expected = panel.ExpectedLineage.iloc[0]
        lineage_match = expected == 'Any_shared_program' or state.CellType in expected.split(';')
        src = sources[sources.Publication.eq(publication)]
        if len(src) == 0:
            continue
        literature_rows.append({'CellState': state.CellState, 'FrozenCellType': state.CellType, 'PublishedProgram': program, 'Publication': publication, 'DOI': src.DOI.iloc[0], 'URL': src.URL.iloc[0], 'PanelN': len(genes), 'DiscoveryAssociatedHitN': len(hits), 'DiscoveryTop20HitN': len(top20 & genes), 'HitGenes': ';'.join(sorted(hits)), 'Top20HitGenes': ';'.join(sorted(top20 & genes)), 'ExternalPositiveBonferroniHitGenes': ';'.join(sorted(external_positive & genes)), 'ExpectedLineage': expected, 'FrozenLineageCompatible': lineage_match, 'MatchingLevel': 'MULTIGENE_PROGRAM_CANDIDATE' if len(hits) >= 2 and lineage_match else ('CROSS_LINEAGE_PROGRAM_ONLY' if hits and not lineage_match else ('SINGLE_GENE_INSUFFICIENT' if len(hits) == 1 else 'NO_PANEL_OVERLAP')), 'PanelSourceStatus': panel.PanelStatus.iloc[0], 'PriorPassageIndex0': panel.BioCPassageIndex0.iloc[0], 'Verification': src.Verification.iloc[0], 'Limit': 'Selected panel overlap; no formal atlas classifier, no complete-signature GSEA, no same-state identity assertion'})
literature = pd.DataFrame(literature_rows)
literature.to_csv(out / 'CELL_STATE_LITERATURE_MATCHING.tsv', sep='\t', index=False)

# Step6. 逐状态证据审阅，置信度针对功能候选解释，细亚型身份仍保留验证要求
curation = pd.read_csv(out / 'audit/STATE_CURATION.tsv', sep='\t').fillna('')
assert set(curation.CellState) == set(frozen.CellState) and curation.CellState.is_unique
master = frozen.merge(curation, on='CellState', validate='one_to_one')
master = master.merge(review.drop(columns=['CellType', 'Ecotype']), on='CellState', validate='one_to_one')
master = master.merge(prior[['CellState', 'ExternalTop10_DiseaseGroup', 'ExternalTop10_Region3', 'LOOCV_EvaluableFolds_DiseaseGroup', 'LOOCV_ConcordantFolds_DiseaseGroup', 'LOOCV_EvaluableFolds_Region3', 'LOOCV_ConcordantFolds_Region3']], on='CellState', validate='one_to_one')
phenotype = pd.read_csv(base / '04_state_phenotypes/state_phenotype_enrichment.tsv', sep='\t')
region = pd.read_csv(base / '04_state_phenotypes/state_region_fisher.tsv', sep='\t')
subtypes = pd.read_csv(lineage / 'tables/state_author_subtype_composition.tsv', sep='\t')
master['AnnotationScope'] = 'Frozen-parent functional candidate; no validated fine-subtype identity'
master['NamingPolicy'] = 'Original State ID retained; Unresolved states receive descriptive program only'
master['EvidenceSourceRoot'] = str(base)
for index, state in master.iterrows():
    state_id = state.CellState
    sig = enrichment[enrichment.CellState.eq(state_id) & enrichment.FDR_Global2650.lt(.05)].sort_values('FDR_Global2650')
    master.loc[index, 'ORA_GlobalFDRlt0p05'] = ';'.join(f'{r.GeneSet}(q={r.FDR_Global2650:.3g};hits={r.OverlapN})' for r in sig.itertuples())
    master.loc[index, 'DiscoveryTop50'] = ';'.join(discovery.loc[discovery.CellState.eq(state_id) & discovery.WithinStateRank.le(50)].sort_values('WithinStateRank').Gene)
    context_support = external[external.CellState.eq(state_id)]
    master.loc[index, 'JointContextAssignmentSummary'] = ';'.join(f'{r.Context}:cells={r.AssignedCellsInContext},donors={r.AssignedDonorsInContext},marker={r.JointContextMarkerStatus}' for r in context_support.itertuples())
    ph = phenotype[phenotype.CellState.eq(state_id) & phenotype.Unit.eq('BiologicalSampleDominant') & phenotype.PositiveZgt1p96.eq(True)]
    master.loc[index, 'PhenotypeSampleDominantPositiveZ'] = ';'.join(f'{r.Dataset}:{r.Feature}:{r.Target}:Z={r.Z:.3g}' for r in ph.itertuples())
    reg = region[region.CellState.eq(state_id) & region.Q.lt(.05)]
    master.loc[index, 'DiscoveryRegionFisherQlt0p05'] = ';'.join(f'{r.Region3}:OR={r.OddsRatio:.3g},q={r.Q:.3g}' for r in reg.itertuples())
    state_lit = literature[literature.CellState.eq(state_id) & literature.FrozenLineageCompatible.eq(True) & literature.DiscoveryAssociatedHitN.ge(2)]
    master.loc[index, 'LiteratureCandidatePrograms'] = ';'.join(f'{r.Publication}:{r.PublishedProgram}[{r.HitGenes}]' for r in state_lit.itertuples())
    change = n4[n4.HeldOut.eq('All') & n4.CellState.eq(state_id)].iloc[0]
    master.loc[index, 'N4_Top10Intersection'] = change.Top10Intersection
    master.loc[index, 'N4_DiseaseOnlyTop'] = change.DiseaseOnlyTop
    master.loc[index, 'N4_RegionOnlyTop'] = change.RegionOnlyTop
    author = subtypes[subtypes.CellState.eq(state_id)].sort_values('DonorBalancedFraction', ascending=False).groupby('Dataset').head(3)
    master.loc[index, 'PooledAuthorSubtypeSummary'] = ';'.join(f'{r.Dataset}:{r.AuthorLabel}={r.DonorBalancedFraction:.2f}' for r in author.itertuples())
    genes = set(discovery.loc[discovery.CellState.eq(state_id), 'Gene'])
    master.loc[index, 'CandidateTFObserved'] = ';'.join(g for g in str(state.CandidateTFPanel).split(';') if g in genes)
    master.loc[index, 'CandidateSurfaceObserved'] = ';'.join(g for g in str(state.CandidateSurfacePanel).split(';') if g in genes)
    master.loc[index, 'HistologyProgramGenesObserved'] = ';'.join(g for g in str(state.HistologyProgramPanel).split(';') if g in genes)
master.to_csv(out / 'CELL_STATE_ANNOTATION_MASTER.tsv', sep='\t', index=False)

# Step7. 所有7生态型保持原始成员，逐成员合并功能与限制
eco_prior = pd.read_csv(base / '06_ecotype_features/7_ecotypes_evidence_summary.tsv', sep='\t')
eco_rows = master[['Ecotype', 'CellState', 'CellType', 'CandidateAnnotation', 'Confidence', 'BiologicalRationale', 'Uncertainty', 'DiscoveryTop20', 'ExternalTop10_DiseaseGroup', 'ExternalTop10_Region3', 'ORA_GlobalFDRlt0p05', 'CandidateTFObserved', 'CandidateSurfaceObserved', 'HistologyProgramGenesObserved', 'ValidationAnchorPanel', 'Colon_AssignedCells', 'Colon_AssignedDonors', 'Ileum_AssignedCells', 'Ileum_AssignedDonors']].merge(eco_prior, on='Ecotype', validate='many_to_one')
eco_rows.to_csv(out / 'ECOTYPE_FUNCTIONAL_COMPOSITION.tsv', sep='\t', index=False)
assert len(eco_rows) == 53
for eco in eco_prior.itertuples():
    assert set(eco.FrozenMembers.split(';')) == set(master.loc[master.Ecotype.eq(eco.Ecotype), 'CellState'])
de_audit = frozen[['CellState', 'CellType', 'Ecotype']].merge(de_coverage[de_coverage.Scope.eq('Dataset')].groupby('CellState', as_index=False).agg(DatasetsWithSavedDE=('Dataset', 'nunique'), SavedDatasetDEGeneRows=('SavedDEGenes', 'sum')), on='CellState', how='left')
de_audit['DiscoveryAssociatedGenes'] = de_audit.CellState.map(discovery.groupby('CellState').size())
de_audit['CompositeTop10Variants'] = de_audit.CellState.map(ranking[ranking.HeldOut.eq('All') & ranking.WithinStateOrder.le(10)].groupby('CellState').Variant.nunique()).fillna(0)
de_audit['ExternalEvaluabilityNote'] = 'Missing DE table can reflect <5 cells, absent parent lineage, or no eligible contrast; not negative marker evidence'
de_audit.to_csv(out / 'audit/STATE_EVIDENCE_COMPLETENESS.tsv', sep='\t', index=False)

# Step8. 输入保护复核和阶段产物摘要
manifest = pd.DataFrame(manifest_rows)
manifest.to_csv(out / 'audit/INPUT_SHA256_BEFORE.tsv', sep='\t', index=False)
after = []
for row in manifest.itertuples():
    path = Path(row.Path)
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    after.append({'Path': str(path), 'SHA256Before': row.SHA256, 'SHA256After': digest, 'Unchanged': digest == row.SHA256})
integrity = pd.DataFrame(after)
integrity.to_csv(out / 'audit/INPUT_IMMUTABILITY_AFTER.tsv', sep='\t', index=False)
assert integrity.Unchanged.all()
de_stat_after = de_path.stat()
assert de_stat_before.st_size == de_stat_after.st_size and de_stat_before.st_mtime_ns == de_stat_after.st_mtime_ns
summary = {'FrozenStates': len(master), 'FrozenEcotypes': master.Ecotype.nunique(), 'Confidence': master.Confidence.value_counts().to_dict(), 'MarkerRows': len(marker), 'ORATests': len(enrichment), 'ORASignificantGlobalFDR': int(enrichment.FDR_Global2650.lt(.05).sum()), 'ExternalContextRows': len(external), 'LiteratureComparisonRows': len(literature), 'InputHashChecksPassed': int(integrity.Unchanged.sum()), 'ExternalDEMetadataUnchanged': True, 'WallSeconds': round(time.time() - started, 2), 'ThreadCap': 1, 'RemoteConnectionsOpened': 0, 'HeavyModelsRun': 0}
(out / 'audit/COMPUTATION_SUMMARY.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)
