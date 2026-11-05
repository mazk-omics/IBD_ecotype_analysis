# 最终比较疾病/Region3两版n4；只在两版都可评价的状态上做成对留出比较。
# 本机cibersortx_review或服务器ibd_refaudit运行；不改任何marker score、DE或冻结状态。
# 一版没有完整七指标排名时，重合率记NA，不把缺少排名误写成0个共同marker。
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy import stats

# Step1. 对全部53状态保留明确覆盖，按三组训练输入分别检查两版排名
out = Path(__file__).resolve().parents[1]
states = pd.read_csv(out / '00_audit/frozen_states.tsv', sep='\t')
ranking = pd.read_csv(out / '02_external_markers/marker_composite_ranking.tsv', sep='\t')
rows = []
for held_out in ranking.HeldOut.unique():
    fold = ranking[ranking.HeldOut.eq(held_out)]
    for state in states.CellState:
        a = fold[fold.CellState.eq(state) & fold.Variant.eq('DiseaseGroup')].set_index('Gene')
        b = fold[fold.CellState.eq(state) & fold.Variant.eq('Region3')].set_index('Gene')
        common = a.index.intersection(b.index)
        top_a = set(a.index[a.WithinStateOrder.le(10)])
        top_b = set(b.index[b.WithinStateOrder.le(10)])
        comparable = bool(top_a and top_b)
        rho = stats.spearmanr(a.loc[common, 'CompositeScore'], b.loc[common, 'CompositeScore']).statistic if comparable and len(common)>=3 else np.nan
        rows.append({'HeldOut': held_out, 'CellState': state, 'ComparisonStatus': 'TOP10_COMPARABLE' if len(top_a)==10 and len(top_b)==10 else 'PARTIAL_TOP_N' if comparable else 'NOT_COMPARABLE', 'DiseaseRankingGenes': len(a), 'RegionRankingGenes': len(b), 'CommonRankableGenes': len(common), 'RankSpearman': rho, 'DiseaseTopN': len(top_a), 'RegionTopN': len(top_b), 'Top10Intersection': len(top_a & top_b) if comparable else np.nan, 'Top10Jaccard': len(top_a & top_b)/len(top_a | top_b) if comparable else np.nan, 'DiseaseOnlyTop': ';'.join(sorted(top_a-top_b)) if comparable else '', 'RegionOnlyTop': ';'.join(sorted(top_b-top_a)) if comparable else ''})
comparison = pd.DataFrame(rows)
comparison.to_csv(out / '02_external_markers/n4_disease_vs_region_comparison.tsv', sep='\t', index=False)
assert len(comparison)==159

# Step2. 两版留出结果先取共同可评价状态，不能拿不同分母的成功率直接判断哪版优
loocv = pd.read_csv(out / '02_external_markers/marker_LOOCV_per_dataset.tsv', sep='\t')
paired = []
for dataset, fold in loocv.groupby('HeldOut'):
    a = fold[fold.Variant.eq('DiseaseGroup') & fold.Status.eq('EVALUABLE')].set_index('CellState')
    b = fold[fold.Variant.eq('Region3') & fold.Status.eq('EVALUABLE')].set_index('CellState')
    common = a.index.intersection(b.index)
    paired.append({'HeldOut': dataset, 'CommonEvaluable': len(common), 'DiseaseConcordant': int(a.loc[common, 'Concordant'].sum()), 'RegionConcordant': int(b.loc[common, 'Concordant'].sum()), 'DiseaseOnlyCorrect': ';'.join(common[(a.loc[common, 'Concordant'].eq(1) & b.loc[common, 'Concordant'].eq(0)).to_numpy()]), 'RegionOnlyCorrect': ';'.join(common[(a.loc[common, 'Concordant'].eq(0) & b.loc[common, 'Concordant'].eq(1)).to_numpy()]), 'DiseaseOnlyEvaluable': ';'.join(a.index.difference(b.index)), 'RegionOnlyEvaluable': ';'.join(b.index.difference(a.index))})
paired = pd.DataFrame(paired)
paired.to_csv(out / '02_external_markers/n4_paired_LOOCV_comparison.tsv', sep='\t', index=False)
main = comparison[comparison.HeldOut.eq('All') & comparison.ComparisonStatus.eq('TOP10_COMPARABLE')]
summary = {'Status': 'PASS', 'MainComparableStates': len(main), 'MedianRankSpearman': float(main.RankSpearman.median()), 'MedianTop10Intersection': float(main.Top10Intersection.median()), 'MinTop10Intersection': int(main.Top10Intersection.min()), 'MaxTop10Intersection': int(main.Top10Intersection.max()), 'PairedLOOCV': paired.to_dict('records'), 'MissingRankingOverlapIsNA': True}
(out / '02_external_markers/n4_final_comparison_summary.json').write_text(json.dumps(summary, indent=2))
report = '# n4疾病覆盖与肠段覆盖：完整比较\n\n两版保留同一其余六指标、候选gene assignment和冻结53状态；各折独立重建排名。\n\n'
report += f'全数据下两版都有完整top10的状态：{len(main)}/53。rank Spearman中位数{main.RankSpearman.median():.3f}；top10重合中位数{main.Top10Intersection.median():.0f}个，范围{main.Top10Intersection.min():.0f}–{main.Top10Intersection.max():.0f}个。一版无法建立排名时重合率记NA，而不是0。\n\n'
report += '| 留出队列 | 共同可评价状态 | 疾病版吻合 | 肠段版吻合 |\n|---|---|---|---|\n'
for row in paired.itertuples():
    report += f'| {row.HeldOut} | {row.CommonEvaluable} | {row.DiseaseConcordant} | {row.RegionConcordant} |\n'
report += '\n回肠共同可评价状态中，B_S03和Fibroblast_S01在疾病版吻合、肠段版不吻合；两版结直肠共同状态的结果一致。少量状态差异和只有两个数据集的留出，不足以证明疾病版普遍更好。\n\n完整单版分母：结直肠疾病34/37、肠段35/38；回肠疾病31/32、肠段33/37。前者的覆盖较窄，不能仅比较百分比。\n\nDisease版强调跨CD/UC/Control的marker对应；Region3版强调跨Ileum/Colon/Rectum。真正疾病/部位特异的状态可能只有少数背景可评价，n4低不自动代表marker差。直肠全部UC、回肠没有UC、部位与数据集绑定的混杂仍在。疾病分层还有不足5细胞导致n4缺失的问题，例如回肠B_S06为CD4细胞、健康3细胞。\n\n建议保留两版，不以留出成功率机械选优。若pan-IBD主报告需要单一入口，可继续用疾病版作为概念上的主定义，同时提供Region3版敏感性对照和每状态覆盖；这不是证明疾病版数值更准确。共用frozen recovery assignment及discovery基因一致性限制，检验的是marker程序对应，不是独立细胞身份准确率。\n'
(out / '09_reports/n4完整比较与建议.md').write_text(report, encoding='utf8')
print(json.dumps(summary, indent=2))
