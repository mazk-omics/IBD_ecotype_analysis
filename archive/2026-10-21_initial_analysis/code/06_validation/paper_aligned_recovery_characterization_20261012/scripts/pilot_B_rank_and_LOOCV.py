# 构建疾病覆盖/肠段覆盖两版七指标，并用严格留数据集法重建排名和评价top10。
# 服务器 ibd_refaudit 环境运行；读取全基因Seurat3.1.3 DEG与冻结gene assignment，不改模型。
# 两版只替换n4；ties用平均rank，未测P不补0，来源/部位混杂写入独立报告。
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy import sparse, stats

# Step1. 汇集逐数据集/父类型/分层的已保存检验；保留全部候选FC用于可评价数据集间求平均
out = Path('/home/mazekai/IBD_EcoTyper/06_validation/paper_aligned_recovery_characterization_20261012')
# 两个B父类型检查点已完成，先做不依赖其他谱系的数值烟雾测试。
for dataset in ['Colon_GSE214695', 'Ileum_IBDverse']:
    assert (out / '02_external_markers/DE_checkpoints_parallel4' / dataset / 'B/completed.txt').read_text().strip() == 'PASS'
pilot_out = out / '02_external_markers/pilot_B'
pilot_out.mkdir(exist_ok=False)
checkpoints = out / '02_external_markers/DE_checkpoints_parallel4'
discovery = pd.read_csv(out / '00_audit/discovery_gene_assignment.tsv', sep='\t').dropna(subset=['State'])
assert not discovery[['CellType', 'Gene']].duplicated().any()
de_parts = []
fc_parts = []
expression_parts = []
for directory in sorted(checkpoints.glob('*/B')):
    dataset, cell_type = directory.parent.name, directory.name
    assert (directory / 'completed.txt').read_text().strip() == 'PASS'
    expression = pd.read_csv(directory / 'gene_expression.tsv', sep='\t')
    expression['Dataset'] = dataset
    expression['CellType'] = cell_type
    expression_parts.append(expression)
    for path in sorted(directory.glob('*_DEG.tsv')):
        scope, stratum = path.name.removesuffix('_DEG.tsv').split('_', 1)
        de = pd.read_csv(path, sep='\t')
        de['Dataset'] = dataset
        de['CellType'] = cell_type
        de['Scope'] = scope
        de['Stratum'] = stratum
        de_parts.append(de)
    for path in sorted(directory.glob('*_candidate_full_FC.tsv')):
        scope, stratum = path.name.removesuffix('_candidate_full_FC.tsv').split('_', 1)
        fc = pd.read_csv(path, sep='\t')
        fc['Dataset'] = dataset
        fc['CellType'] = cell_type
        fc['Scope'] = scope
        fc['Stratum'] = stratum
        fc_parts.append(fc)
de = pd.concat(de_parts, ignore_index=True)
folds = pd.concat(fc_parts, ignore_index=True)
expression = pd.concat(expression_parts, ignore_index=True)
de.to_csv(pilot_out / 'external_DEG_by_dataset.tsv.gz', sep='\t', index=False)
de_keys = ['Dataset', 'CellType', 'Scope', 'Stratum', 'Gene', 'State']
assert not de[de_keys].duplicated().any()
positive = de[de.log2FC.gt(0)].copy()
max_fc = positive.groupby(de_keys[:-1]).log2FC.transform('max')
assigned = positive[np.isclose(positive.log2FC, max_fc, rtol=0, atol=1e-12)].copy()
tie_sizes = assigned.groupby(de_keys[:-1]).State.transform('size')
assigned['UniqueDatasetAssignment'] = tie_sizes.eq(1)
assigned.to_csv(pilot_out / 'external_gene_state_assignment.tsv.gz', sep='\t', index=False)
eligibility = pd.read_csv(out / '02_external_markers/DE_stratum_eligibility_parallel4.tsv', sep='\t')
datasets = sorted(expression.Dataset.unique())
all_rankings = []
all_metrics = []
exclusions = []

# Step2. 所有数据集和每个留出折都重建基因归属、n1–n7、n4，避免训练集泄漏
training_sets = [('All', datasets)] + [(dataset, [x for x in datasets if x!=dataset]) for dataset in datasets]
for held_out, training in training_sets:
    train_de = de[de.Dataset.isin(training)]
    train_folds = folds[folds.Dataset.isin(training)]
    train_assign = assigned[assigned.Dataset.isin(training)]
    train_expression = expression[expression.Dataset.isin(training)]
    overall_de = train_de[train_de.Scope.eq('Dataset') & train_de.log2FC.gt(0)].copy()
    maximum = overall_de.groupby(['CellType', 'Gene']).log2FC.transform('max')
    winners = overall_de[np.isclose(overall_de.log2FC, maximum, atol=1e-12, rtol=0)].copy()
    means = train_folds[train_folds.Scope.eq('Dataset')].groupby(['CellType', 'Gene', 'State']).log2FC.mean().rename('MeanAcrossEvaluableDatasets').reset_index()
    winners = winners.merge(means, on=['CellType', 'Gene', 'State'], validate='many_to_one')
    winner_mean = winners.groupby(['CellType', 'Gene']).MeanAcrossEvaluableDatasets.transform('max')
    winners = winners[np.isclose(winners.MeanAcrossEvaluableDatasets, winner_mean, atol=1e-12, rtol=0)]
    winners = winners[['CellType', 'Gene', 'State']].drop_duplicates()
    unresolved = winners[winners[['CellType', 'Gene']].duplicated(keep=False)]
    # 双重并列仍无法判别时隔离，不任意分配；排名内部的并列另用平均rank。
    winners = winners[~winners[['CellType', 'Gene']].duplicated(keep=False)]
    candidates = winners.merge(discovery[['CellType', 'Gene', 'State']], on=['CellType', 'Gene'], suffixes=('_External', '_Discovery'), validate='one_to_one')
    consistent = candidates[candidates.State_External.eq(candidates.State_Discovery)].copy()
    consistent['State'] = consistent.State_External
    consistent['CellState'] = consistent.CellType + '_' + consistent.State
    exclusions.append({'HeldOut': held_out, 'TrainingDatasets': ';'.join(training), 'UnresolvedGeneTies': len(unresolved), 'DiscoveryMatchedGenes': len(candidates), 'ConsistentCandidates': len(consistent), 'DiscordantExcluded': len(candidates)-len(consistent)})
    # 先按关键列分组查找，避免逐基因扫描整张DE表；统计定义和候选集合不变。
    expression_n = train_expression[train_expression.MeanCPM.gt(0)].groupby(['CellType', 'Gene']).Dataset.nunique().to_dict()
    assigned_n = train_assign[train_assign.Scope.eq('Dataset') & train_assign.UniqueDatasetAssignment].groupby(['CellType', 'Gene', 'State']).Dataset.nunique().to_dict()
    specific_de = train_de[train_de.Scope.eq('Dataset')].merge(consistent[['CellType', 'Gene', 'State']], on=['CellType', 'Gene', 'State'], validate='many_to_one')
    de_lookup = dict(tuple(specific_de.groupby(['CellType', 'Gene', 'State'])))
    fc_lookup = means.set_index(['CellType', 'Gene', 'State']).MeanAcrossEvaluableDatasets.to_dict()
    strata_assignment_n = train_assign[train_assign.UniqueDatasetAssignment].groupby(['CellType', 'Gene', 'State', 'Scope']).Stratum.nunique().to_dict()
    strata_eligible_n = eligibility[eligibility.Dataset.isin(training) & eligibility.Eligible & eligibility.EvaluableContrast].groupby(['CellType', 'State', 'Scope']).Stratum.nunique().to_dict()
    metric_rows = []
    for row in consistent.itertuples():
        key = (row.CellType, row.Gene, row.State)
        specific = de_lookup.get(key, specific_de.iloc[:0])
        mean_fc = fc_lookup[key]
        zs = stats.norm.isf(np.maximum(specific.p_val.to_numpy(), 1e-300)/2) * np.sign(specific.log2FC.to_numpy())
        meta_z = zs.sum()/np.sqrt(len(zs)) if len(zs) else np.nan
        metric = {'HeldOut': held_out, 'TrainingDatasets': ';'.join(training), 'CellState': row.CellState, 'CellType': row.CellType, 'State': row.State, 'Gene': row.Gene, 'n1': expression_n.get((row.CellType, row.Gene), 0), 'n2': assigned_n.get(key, 0), 'n5': int(specific.p_val_adj.lt(0.05).sum()), 'MetaZUnadjusted': meta_z, 'MetaP': 2*stats.norm.sf(abs(meta_z)), 'n7': mean_fc, 'TestedDatasets': len(specific)}
        metric['n3'] = metric['n2']/metric['n1'] if metric['n1']>0 else np.nan
        for scope in ['DiseaseGroup', 'Region3']:
            evaluated_n = strata_eligible_n.get((row.CellType, row.State, scope), 0)
            metric['n4_'+scope] = strata_assignment_n.get((*key, scope), 0) if evaluated_n>0 else np.nan
            metric['Evaluable_'+scope] = evaluated_n
        metric_rows.append(metric)
    metrics = pd.DataFrame(metric_rows)
    assert len(metrics)>0
    metrics['n6'] = np.nan
    for state, idx in metrics.groupby('CellState').groups.items():
        tested_idx = idx[metrics.loc[idx, 'MetaP'].notna()]
        q = stats.false_discovery_control(metrics.loc[tested_idx, 'MetaP'].to_numpy(), method='bh')
        metrics.loc[tested_idx, 'MetaQ'] = q
        metrics.loc[tested_idx, 'n6'] = stats.norm.isf(np.maximum(q, 1e-300)/2) * np.sign(metrics.loc[tested_idx, 'MetaZUnadjusted'])
    all_metrics.append(metrics)
    for variant in ['DiseaseGroup', 'Region3']:
        ranking = metrics.copy()
        ranking['n4'] = ranking['n4_'+variant]
        metrics7 = ['n1', 'n2', 'n3', 'n4', 'n5', 'n6', 'n7']
        ranking['AllSevenEvaluable'] = ranking[metrics7].notna().all(axis=1)
        ranking = ranking[ranking.AllSevenEvaluable].copy()
        # 每个状态内部rank越大越好；平均rank保证七项等权，不使用额外置信度加权。
        ranks = ranking.groupby('CellState')[metrics7].rank(method='average', ascending=True)
        ranking['CompositeScore'] = ranks.mean(axis=1)
        ranking['Variant'] = variant
        ranking = ranking.sort_values(['CellState', 'CompositeScore', 'Gene'], ascending=[True, False, True])
        ranking['WithinStateOrder'] = ranking.groupby('CellState').cumcount()+1
        all_rankings.append(ranking)
metrics = pd.concat(all_metrics, ignore_index=True)
rankings = pd.concat(all_rankings, ignore_index=True)
metrics.to_csv(pilot_out / 'marker_n1_n7.tsv', sep='\t', index=False)
rankings.to_csv(pilot_out / 'marker_composite_ranking.tsv', sep='\t', index=False)
rankings[rankings.WithinStateOrder.le(10)].to_csv(pilot_out / 'state_top_markers.tsv', sep='\t', index=False)
pd.DataFrame(exclusions).to_csv(pilot_out / 'gene_assignment_exclusions.tsv', sep='\t', index=False)

# Step3. 在同一候选基因集合上比较两版rank与top10；不因某版显著更多就选该版
comparison = []
for (held_out, state), group in rankings.groupby(['HeldOut', 'CellState']):
    disease = group[group.Variant.eq('DiseaseGroup')].set_index('Gene')
    site = group[group.Variant.eq('Region3')].set_index('Gene')
    common = disease.index.intersection(site.index)
    top_a = set(disease.index[disease.WithinStateOrder.le(10)])
    top_b = set(site.index[site.WithinStateOrder.le(10)])
    rho = stats.spearmanr(disease.loc[common, 'CompositeScore'], site.loc[common, 'CompositeScore']).statistic if len(common)>=3 else np.nan
    comparison.append({'HeldOut': held_out, 'CellState': state, 'CommonRankableGenes': len(common), 'RankSpearman': rho, 'DiseaseTopN': len(top_a), 'RegionTopN': len(top_b), 'Top10Intersection': len(top_a & top_b), 'Top10Jaccard': len(top_a & top_b)/len(top_a | top_b) if top_a | top_b else np.nan, 'DiseaseOnlyTop': ';'.join(sorted(top_a-top_b)), 'RegionOnlyTop': ';'.join(sorted(top_b-top_a))})
pd.DataFrame(comparison).to_csv(pilot_out / 'n4_disease_vs_region_comparison.tsv', sep='\t', index=False)

# Step4. 留出队列log2(CPM+1)后，按对应父类型全部输入细胞单位方差缩放；不使用留出DE建排名
loocv = []
loocv_expression = []
for held_out in datasets:
    train_ranking = rankings[rankings.HeldOut.eq(held_out) & rankings.WithinStateOrder.le(10)]
    for directory in sorted((out / '02_external_markers/full_counts' / held_out).iterdir()):
        cell_type = directory.name
        raw = sparse.load_npz(directory / 'counts.npz')
        genes = pd.read_csv(directory / 'genes.tsv', header=None).iloc[:, 0]
        cells = pd.read_csv(directory / 'cells.tsv', sep='\t')
        totals = np.asarray(raw.sum(axis=0)).ravel()
        state_counts = cells.loc[cells.State.ne('Unassigned'), 'State'].value_counts()
        allowed = state_counts.index[state_counts.ge(5)]
        candidates = train_ranking[train_ranking.CellType.eq(cell_type)]
        for (variant, state), markers in candidates.groupby(['Variant', 'State']):
            target = cell_type+'_'+state
            available = markers[markers.Gene.isin(genes)]
            if len(markers)!=10 or len(available)!=10 or state not in allowed or len(allowed)<2:
                loocv.append({'HeldOut': held_out, 'Variant': variant, 'CellState': target, 'Status': 'NOT_EVALUABLE', 'Reason': 'Need 10 trained/present markers, >=5 assigned target cells and >=2 eligible held-out states', 'TrainedMarkers': len(markers), 'PresentMarkers': len(available), 'Concordant': np.nan})
                continue
            index = pd.Index(genes).get_indexer(available.Gene)
            values = raw[index, :].toarray() * (1e6/totals)
            values = np.log2(values+1)
            sd = values.std(axis=1, ddof=1)
            if np.any(sd==0):
                loocv.append({'HeldOut': held_out, 'Variant': variant, 'CellState': target, 'Status': 'NOT_EVALUABLE', 'Reason': 'Zero marker variance in held-out parent cells', 'TrainedMarkers': 10, 'PresentMarkers': 10, 'Concordant': np.nan})
                continue
            values = (values-values.mean(axis=1, keepdims=True))/sd[:, None]
            scores = pd.Series({other: values[:, cells.State.eq(other)].mean() for other in allowed})
            winners = scores.index[np.isclose(scores, scores.max(), atol=1e-12, rtol=0)]
            for other, score in scores.items():
                loocv_expression.append({'HeldOut': held_out, 'Variant': variant, 'MarkerState': target, 'ExpressionState': cell_type+'_'+other, 'MeanUnitVarianceExpression': score})
            loocv.append({'HeldOut': held_out, 'Variant': variant, 'CellState': target, 'Status': 'EVALUABLE' if len(winners)==1 else 'TIE', 'PredictedState': ';'.join(cell_type+'_'+x for x in winners), 'Concordant': float(len(winners)==1 and winners[0]==state), 'TrainedMarkers': 10, 'PresentMarkers': 10, 'TopScore': scores.max(), 'TargetScore': scores[state]})
loocv = pd.DataFrame(loocv)
frozen_states = pd.read_csv(out / '00_audit/frozen_states.tsv', sep='\t')
frozen_states = frozen_states[frozen_states.CellType.eq('B')]
full_grid = pd.MultiIndex.from_product([datasets, ['DiseaseGroup', 'Region3'], frozen_states.CellState], names=['HeldOut', 'Variant', 'CellState']).to_frame(index=False)
loocv = full_grid.merge(loocv, on=['HeldOut', 'Variant', 'CellState'], how='left', validate='one_to_one')
missing_folds = loocv.Status.isna()
loocv['Status'] = loocv.Status.fillna('NOT_EVALUABLE')
if 'Reason' not in loocv:
    loocv['Reason'] = ''
loocv['Reason'] = loocv.Reason.fillna('')
loocv.loc[missing_folds, 'Reason'] = 'No trained seven-metric ranking or no input parent lineage'
loocv.to_csv(pilot_out / 'marker_LOOCV_per_dataset.tsv', sep='\t', index=False)
state_summary = loocv.groupby(['CellState', 'Variant']).agg(TotalFolds=('HeldOut', 'size'), EvaluableFolds=('Status', lambda x: int(x.eq('EVALUABLE').sum())), ConcordantFolds=('Concordant', 'sum')).reset_index()
state_summary['Concordance'] = state_summary.ConcordantFolds/state_summary.EvaluableFolds.replace(0, np.nan)
state_summary.to_csv(pilot_out / 'marker_LOOCV_per_state.tsv', sep='\t', index=False)
pd.DataFrame(loocv_expression).to_csv(pilot_out / 'marker_LOOCV_mean_expression.tsv', sep='\t', index=False)
summary = loocv[loocv.Status.eq('EVALUABLE')].groupby(['HeldOut', 'Variant']).agg(EvaluableStates=('CellState', 'size'), ConcordantStates=('Concordant', 'sum'), Concordance=('Concordant', 'mean')).reset_index()
summary.to_csv(pilot_out / 'marker_LOOCV_summary.tsv', sep='\t', index=False)
(pilot_out / 'marker_LOOCV_summary.md').write_text('# 两版n4的实际留出结果\n\n每次整队列留出后重建排名，训练仅另一个队列；完整53状态逐折包含NOT_EVALUABLE。下表分母只含可评价且无并列的状态，不能解释成全部53状态的成功率。\n\n'+summary.to_csv(sep='\t', index=False)+'\n疾病/部位来源混杂见n4比较说明。', encoding='utf8')
(pilot_out / 'ranking_LOOCV_complete.json').write_text(json.dumps({'TrainingFolds': training_sets, 'MetricRows': len(metrics), 'RankingRows': len(rankings), 'StrataDEPerformed': True, 'NoHeldOutDEUsedForRanking': True, 'Scope': 'B only pilot; not all53 states', 'n4Variants': ['DiseaseGroup', 'Region3'], 'LOOCV': summary.to_dict('records')}, indent=2))
print(summary.to_string(index=False), flush=True)
