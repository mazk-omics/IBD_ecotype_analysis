# 根据复算证据建立20类父谱系的条件化可靠性矩阵，不用事后数值阈值筛选好看的类型。
# 本机运行；A/B/C/D为有解释依据的收尾证据分级，不是预设统计验收或新训练分类器。
from pathlib import Path
import json
import pandas as pd

# Step1. 读取全部尺度与检出结果，CIBERSORTx主分析和BP敏感性分开记录
out = Path(__file__).resolve().parents[1]
root = out.parent
types = (root / '00_audit/frozen_types.txt').read_text().splitlines()
metrics = pd.read_csv(out / 'EXTERNAL_RECOMPUTED_METRICS.tsv', sep='\t')
zeros = pd.read_csv(out / 'DETECTION_AND_GENE_DENOMINATORS.tsv', sep='\t')
capacity = pd.read_csv(root / '01_reference/reference_type_capacity.tsv', sep='\t').set_index('cell_type')

# Step2. 先明确分级理由与适用条件，不能把BP自己的基因分母作为唯一中立评分
judgments = {
    'B': ('B', 'B', 'B', 'B', 'B', 'GSE214695有22/58个真值阳性混合物精确零预测；两数据集数量bias为负，IBDverse约-5.61 pp。RNA量较低误差不抵消检出和数量校准限制。', '可讨论组成变化，尤其低丰度或绝对数量不可靠；不能把B与Plasma差异都归因于某一种生物机制。'),
    'Plasma': ('B/C', 'C', 'B', 'B', 'A/B', 'GSE214695数量bias +27.50 pp，完整RNA仍+5.21 pp；IBDverse完整RNA误差较低。BP保留基因分母显著改变Plasma truth。', 'C仅指GSE214695按真实细胞数量解读的严重校准失败；不是所有RNA尺度或所有组织失败。'),
    'CD4_T': ('B', 'B', 'A/B', 'B', 'A/B', '两数据集数量bias约-4.07/-4.46 pp，相关仍可追踪；完整RNA MAE约1.37/1.60 pp。', '相对变化有支持，精确细胞数量有系统低估；不据此辨认CD4内部亚型。'),
    'CD8_T': ('A/B', 'B', 'A', 'B', 'A', '两数据集完整RNA MAE约0.78/0.82 pp，数量MAE较高且偏低；阳性漏检0/58及2/54。', 'A限于完整RNA和所测模拟条件；数量百分比与低丰度检出仍需B级警惕。'),
    'NK': ('A/B', 'B', 'A/B', 'D', 'D', 'GSE214695有50个阳性混合物、漏检1个；IBDverse没有明确非零设计真值。', '支持限于GSE214695结肠模拟；IBDverse非零预测只能作为设计中的阴性控制，不能代表回肠检出已验证。'),
    'ILC': ('B/C', 'B', 'B', 'C', 'C', 'GSE214695阳性漏检13/46；IBDverse 100/100预测0，49/49阳性全部漏检，相关不可计算。', 'IBDverse该配置中不可作可靠检出或数量指标；GSE214695也存在低估，不能整体列为已充分验证。'),
    'Monocyte_Macrophage': ('B', 'B', 'A/B', 'B', 'B', 'GSE214695漏检8/48，IBDverse漏检2/54；数量偏低，RNA表现和signature中与DC的相似性需并看。', '父谱系变化有一定支持，绝对数量及与DC的分界需谨慎；共线性不证明reads去向。'),
    'DC': ('B', 'B', 'B', 'B', 'B', '数量r从GSE214695约0.95降至IBDverse约0.50；IBDverse设计无DC样本中的最大预测约10.04%。', '数据集与truth依赖明显，低丰度非零不能自动当作DC浸润证据。'),
    'Mast': ('A/B', 'B', 'A', 'B', 'A', '完整RNA MAE约0.38/0.37 pp，阳性漏检1/48及0/50；细胞数量仍约2.26/2.96 pp误差且偏低。', 'RNA恢复相对更有支持；不将A外推为精确数量计数。'),
    'Fibroblast': ('A/B', 'A/B', 'A/B', 'B', 'B', '两数据集相关可追踪，GSE214695校准相对较好；IBDverse数量与完整RNA存在低估，漏检3/56。', '数量和状态构成可能随数据集变化，不能解释为所有成纤维亚群均准确。'),
    'Pericyte': ('A/B', 'A/B', 'A/B', 'A/B', 'A/B', '两数据集阳性均无精确零漏检；数量MAE约1.60/2.15 pp，仍有条件依赖的bias及无真值样本非零预测。', '在所测范围相对稳定；A/B不是无误计数或普适检出限承诺。'),
    'Endothelial': ('B', 'B', 'B', 'B', 'B', '阳性漏检10/53及15/53；完整RNA bias两数据集均约-2.5 pp。', '存在可重复的数量/检出低估，不宜以近零结果认定缺乏内皮。'),
    'Ileal_absorptive': ('B', 'D', 'D', 'B', 'B', '仅IBDverse有非零验证，数量bias约+9.21 pp、完整RNA bias约+2.90 pp。GSE214695为阴性设计。', '单一回肠队列证据且系统高估，不能声称已在所有组织中校准。'),
    'Colonic_absorptive': ('B', 'B', 'B', 'D', 'D', '仅GSE214695有非零验证；设计无该类时仍有非零预测，IBDverse阴性设计最大约9.64%。', '结肠模拟内有相对变化支持；回肠负对照不代表结肠类型可在回肠作生物数量解释。'),
    'Goblet': ('B', 'B', 'B', 'B', 'B', 'GSE214695表现较好，IBDverse数量bias +6.83 pp、完整RNA约+4.73 pp；阳性无精确零漏检。', '变化趋势有支持，但数据集相关的高估影响定量，不能只看高r。'),
    'Stem_TA_progenitor': ('B/C', 'C', 'C', 'B', 'B', 'GSE214695 97/100预测0、阳性漏检51/54；IBDverse能追踪变化但数量+8.79 pp、完整RNA-3.50 pp。', 'GSE214695配置存在明确检出失败；IBDverse也不能作精确计数。Stem/TA为粗类型，作者标签对应不完全。'),
    'Paneth': ('D', 'D', 'D', 'D', 'D', '两个数据集均无非零真值；IBDverse设计阴性时平均预测约2.62%、最大17.13%。', '外部灵敏度和校准未验证，且有阴性条件假阳性风险；不得以真值全零时低MAE认定通过。'),
    'Enteroendocrine': ('A/B', 'B', 'A', 'A/B', 'A/B', '完整RNA MAE约0.33/1.24 pp，阳性漏检1/56及1/49；数量尺度仍有bias。', 'RNA与组成趋势相对更有支持，稀有类型的精确数量和亚型仍需谨慎。'),
    'Glial': ('A/B', 'B', 'A', 'D', 'D', 'GSE214695完整RNA MAE约0.31 pp、阳性无漏检；IBDverse没有明确非零真值。', 'A仅限所测结肠完整RNA；不外推回肠检出或独立亚型有效性。'),
    'Tuft': ('B', 'B', 'A/B', 'B', 'B', 'CX GSE214695完整RNA误差较低；IBDverse阳性漏检12/54。BP在IBDverse主要RNA与cell-count均严重高估。', 'CX主路线仍有低丰度与尺度限制；BP敏感性中该谱系为C级异常，不以BP失败自动判定CX失败。')
}
bp_notes = {'Tuft': 'C in IBDverse: retained-RNA bias +17.03 pp and cell-count MAE 18.11 pp; mechanism unresolved', 'Stem_TA_progenitor': 'C in IBDverse retained-RNA calibration (bias -13.22 pp), while cell-count MAE 3.36 pp; scope depends on estimand', 'Plasma': 'Gene-set dependence: only 17.79%/33.03% of type-associated RNA retained; BP target differs from full-gene RNA'}
rows = []
for cell_type in types:
    grade, g_cell, g_rna, i_cell, i_rna, basis, limitation = judgments[cell_type]
    row = {'cell_type': cell_type, 'primary_method': 'CIBERSORTx frozen Joint20 v2', 'overall_CIBERSORTx_evidence': grade, 'GSE214695_cell_count_grade': g_cell, 'GSE214695_full_RNA_grade': g_rna, 'IBDverse_cell_count_grade': i_cell, 'IBDverse_full_RNA_grade': i_rna, 'reference_cells': int(capacity.loc[cell_type, 'cells']), 'reference_donors': int(capacity.loc[cell_type, 'donors']), 'reference_profiles': int(capacity.loc[cell_type, 'profiles']), 'basis': basis, 'interpretation_limit': limitation, 'BP_sensitivity_note': bp_notes.get(cell_type, 'All three BP truth metrics retained separately; no inference that all methods are interchangeable'), 'grading_rule': 'Post-hoc descriptive evidence appraisal with stated conditions; no tuned acceptance cutoff or pass/fail score'}
    for dataset in ['GSE214695', 'IBDverse']:
        for scale, label in [('Cell_count_SECONDARY', 'cell'), ('Full_gene_RNA_SECONDARY', 'fullRNA'), ('Retained_gene_RNA_PRIMARY', 'BPgeneRNA')]:
            block = metrics.loc[metrics.dataset.eq(dataset) & metrics.truth_scale.eq(scale) & metrics.cell_type.eq(cell_type)]
            for method in ['CIBERSORTx', 'BayesPrism_final']:
                part = block.loc[block.method.eq(method)].iloc[0]
                for metric in ['Pearson_r', 'Spearman_rho', 'MAE_pp', 'RMSE_pp', 'bias_pp']:
                    row[dataset + '_' + method + '_' + label + '_' + metric] = part[metric]
            row[dataset + '_evaluable'] = bool(block.evaluable.iloc[0])
        z = zeros.loc[zeros.dataset.eq(dataset) & zeros.method.eq('CIBERSORTx') & zeros.cell_type.eq(cell_type)].iloc[0]
        for metric in ['n_true_positive', 'n_true_absent', 'false_negative_exact_zero', 'false_positive_exact_positive', 'false_positive_above_0p01pp', 'mean_prediction_when_absent_percent', 'max_prediction_when_absent_percent', 'retained_fraction_of_type_RNA']:
            row[dataset + '_CX_' + metric] = z[metric]
    rows.append(row)
reliability = pd.DataFrame(rows)
assert len(reliability) == 20 and reliability.cell_type.is_unique
assert reliability.loc[reliability.cell_type.eq('Paneth'), 'overall_CIBERSORTx_evidence'].item() == 'D'
reliability.to_csv(out / 'JOINT20_CELLTYPE_RELIABILITY.tsv', sep='\t', index=False)
print(reliability[['cell_type', 'overall_CIBERSORTx_evidence', 'GSE214695_cell_count_grade', 'IBDverse_cell_count_grade']].to_string(index=False))
(out / 'reliability_grading_context.json').write_text(json.dumps({'primary_method': 'CIBERSORTx', 'grading': 'Qualitative, retrospective, condition-specific evidence appraisal; no newly optimized cutoff', 'A': 'Reasonably supported in explicitly named observed context and estimand', 'B': 'Substantial quantitative caveats', 'C': 'Clear failure in named tested condition, not universal biological invalidity', 'D': 'Insufficient positive external validation; negative-design false positives may still be reported', 'BP_retained_RNA': 'Method-specific; retained numerically but not treated as uniquely neutral grading scale for CIBERSORTx'}, indent=2), encoding='utf-8')
