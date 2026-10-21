# 复用新取得的Garrido-Trigo原始Source Data补充已知巨噬signature。
# 只补充文献匹配；不重复全状态表达比较或DE，不把同队列对应当独立验证。
import os
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
from pathlib import Path
import math
import pandas as pd
import numpy as np

# Step1. 全部发表导出marker和成对DE，方向及来源单独保存
root = Path('D:/IBD_ecotype/IBD_ecotype_analysis')
out = Path(__file__).resolve().parents[1]
file = out / 'literature_resources/GarridoTrigo2023__Source_data_file_supplementary_figure_5.xlsx'
cluster = pd.read_excel(file, sheet_name='e.fig 5d', engine='openpyxl')
contrast = pd.read_excel(file, sheet_name='e.fig 5e', engine='openpyxl')
cluster.to_csv(out / 'literature_resources/GarridoTrigo2023_COMPLETE_EXPORTED_CLUSTER_MARKERS.tsv', sep='\t', index=False)
contrast.to_csv(out / 'literature_resources/GarridoTrigo2023_COMPLETE_EXPORTED_PAIRWISE_DE.tsv', sep='\t', index=False)
print('Cluster marker rows and counts:', len(cluster), cluster.groupby('cluster').size().to_dict(), flush=True)
print('Pairwise DE rows and counts:', len(contrast), contrast.groupby('comp').size().to_dict(), flush=True)
cluster = cluster.rename(columns={'cluster': 'PublishedState', 'gene': 'Gene', 'avg_log2FC': 'PublishedEffect', 'p_val_adj': 'PublishedAdjustedP'})
cluster = cluster[cluster.PublishedEffect.gt(0) & cluster.PublishedAdjustedP.lt(.05)].copy()
cluster['PublishedSignedEffect'] = cluster.PublishedEffect
cluster['SignatureScope'] = 'Complete published FigS5d exported top-marker lists (paper uses top200), not exhaustive transcriptome DE'
cluster['SourceFile'] = file.name + ':e.fig 5d'
cluster['Direction'] = 'up'
contrast = contrast.rename(columns={'gene': 'Gene', 'avg_log2FC': 'PublishedEffect', 'p_val_adj': 'PublishedAdjustedP'})
contrast = contrast[contrast.PublishedAdjustedP.lt(.05)].copy()
contrast['PublishedSignedEffect'] = contrast.PublishedEffect
contrast['PublishedState'] = contrast['comp'] + np.where(contrast.PublishedEffect > 0, '_UP_in_IDA', '_DOWN_in_IDA')
contrast['Direction'] = np.where(contrast.PublishedEffect > 0, 'up_in_IDA', 'down_in_IDA')
contrast['PublishedEffect'] = contrast.PublishedEffect.abs()
contrast['SignatureScope'] = 'Complete published FigS5e pairwise DEG export, supplied adjustedP<0.05; UP/DOWN relative to IDA retained'
contrast['SourceFile'] = file.name + ':e.fig 5e'
added = pd.concat([cluster, contrast], ignore_index=True, sort=False)
added['Publication'] = 'GarridoTrigo2023'
added['DOI'] = '10.1038/s41467-023-40156-6'
added['ExpectedParent'] = 'Monocyte_Macrophage'
old_signature = pd.read_csv(out / 'literature_resources/PUBLISHED_SIGNATURE_GENES.tsv', sep='\t', low_memory=False)
signature = pd.concat([old_signature[old_signature.Publication.ne('GarridoTrigo2023')], added], ignore_index=True, sort=False).drop_duplicates(['Publication', 'PublishedState', 'Gene'])
signature.to_csv(out / 'literature_resources/PUBLISHED_SIGNATURE_GENES.tsv', sep='\t', index=False)
signature.groupby(['Publication', 'PublishedState', 'ExpectedParent', 'SignatureScope'], dropna=False).size().rename('PublishedGeneN').reset_index().to_csv(out / 'audit/PUBLISHED_SIGNATURE_DEFINITIONS.tsv', sep='\t', index=False)

# Step2. 三个髓系状态与全部新signature做背景限制的候选比较
old = root / 'analysis_reports/legacy_bayesprism_annotation/run_20261013_01'
master = pd.read_csv(old / 'CELL_STATE_ANNOTATION_MASTER.tsv', sep='\t')
discovery = pd.read_csv(root / '06_validation/paper_aligned_recovery_characterization_20261012/00_audit/prior_discovery_associated_genes.tsv', sep='\t')
ranking = pd.read_csv(root / '06_validation/paper_aligned_recovery_characterization_20261012/02_external_markers/marker_composite_ranking.tsv', sep='\t')
backgrounds = pd.read_csv(old / 'audit/ORA_GENE_BACKGROUNDS.tsv', sep='\t')
universe = set(backgrounds.loc[backgrounds.CellType.eq('Monocyte_Macrophage'), 'Gene'])
new_rows = []
for state in master[master.CellType.eq('Monocyte_Macrophage')].itertuples():
    queries = {'DiscoveryTop50': set(discovery.loc[discovery.CellState.eq(state.CellState) & discovery.WithinStateRank.le(50), 'Gene'])}
    for variant in ['DiseaseGroup', 'Region3']:
        queries['ExternalTop50_' + variant] = set(ranking.loc[ranking.CellState.eq(state.CellState) & ranking.HeldOut.eq('All') & ranking.Variant.eq(variant) & ranking.WithinStateOrder.le(50), 'Gene'])
    for published, table in added.groupby('PublishedState', sort=False):
        target = set(table.Gene) & universe
        for query_name, query_genes in queries.items():
            selected = query_genes & universe
            N, n, K = len(universe), len(selected), len(target)
            hits = selected & target
            k = len(hits)
            evaluable = n >= 5 and K >= 5
            p = np.nan
            if evaluable:
                denominator = math.lgamma(N + 1) - math.lgamma(n + 1) - math.lgamma(N - n + 1)
                values = []
                for x in range(k, min(K, n) + 1):
                    if 0 <= n - x <= N - K:
                        v = math.lgamma(K + 1) - math.lgamma(x + 1) - math.lgamma(K - x + 1)
                        v += math.lgamma(N - K + 1) - math.lgamma(n - x + 1) - math.lgamma(N - K - n + x + 1) - denominator
                        values.append(math.exp(v))
                p = min(1, math.fsum(values))
            new_rows.append({'CellState': state.CellState, 'CellType': state.CellType, 'Publication': 'GarridoTrigo2023', 'DOI': '10.1038/s41467-023-40156-6', 'PublishedStateOrProgram': published, 'PublishedExpectedParent': 'Monocyte_Macrophage', 'ComparisonScope': 'same_parent_or_shared_program', 'Query': query_name, 'BackgroundN': N, 'QueryInBackgroundN': n, 'PublishedSignatureN': table.Gene.nunique(), 'SignatureInBackgroundN': K, 'OverlapN': k, 'OverlapGenes': ';'.join(sorted(hits)), 'JaccardInBackground': k / (n + K - k) if n + K - k else np.nan, 'QueryCoverage': k / n if n else np.nan, 'SignatureCoverageInBackground': k / K if K else np.nan, 'P': p, 'Evaluable': evaluable, 'EffectDirection': table.Direction.iloc[0], 'SignatureSourceScope': table.SignatureScope.iloc[0], 'EvidenceIndependence': 'GarridoTrigo2023 is GSE214695 recovery source cohort; published matching and recovery support are not independent', 'Limit': 'Same-cohort conditional overlap; pairwise DOWN is relative depletion in IDA, not a new cell identity'})
matches = pd.read_csv(out / 'audit/ALL_PUBLISHED_SIGNATURE_MATCH_STATISTICS.tsv', sep='\t')
matches = pd.concat([matches[matches.Publication.ne('GarridoTrigo2023')], pd.DataFrame(new_rows)], ignore_index=True)
matches['BH_Q_GlobalEvaluableMatches'] = np.nan
eligible = matches[matches.Evaluable].sort_values('P')
q = eligible.P.to_numpy() * len(eligible) / np.arange(1, len(eligible) + 1)
matches.loc[eligible.index, 'BH_Q_GlobalEvaluableMatches'] = np.minimum(1, np.minimum.accumulate(q[::-1])[::-1])
matches.to_csv(out / 'audit/ALL_PUBLISHED_SIGNATURE_MATCH_STATISTICS.tsv', sep='\t', index=False)
best = matches[matches.Query.eq('DiscoveryTop50') & matches.ComparisonScope.eq('same_parent_or_shared_program') & matches.Evaluable].sort_values(['BH_Q_GlobalEvaluableMatches', 'OverlapN'], ascending=[True, False]).groupby('CellState').head(3)
best.to_csv(out / 'audit/BEST_DISCOVERY_PUBLISHED_MATCHES.tsv', sep='\t', index=False)
print('All comparisons/evaluable:', len(matches), len(eligible), flush=True)
print(matches[matches.Publication.eq('GarridoTrigo2023') & matches.BH_Q_GlobalEvaluableMatches.lt(.05)][['CellState', 'PublishedStateOrProgram', 'Query', 'OverlapN', 'BH_Q_GlobalEvaluableMatches']].to_string(index=False), flush=True)
