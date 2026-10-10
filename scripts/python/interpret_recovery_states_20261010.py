# 合并冻结53状态与bulk、回肠、结直肠recovery，给出附来源的初步程序判断。
# 本机cibersortx_review环境运行；只读已有结果，写入独立recovery_interpretation_20261010目录。
# 不改原始注释、表达、W、状态或生态型；重叠与分配统计不是准确率或正式富集检验。
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

# Step1. 读入冻结定义和三轮恢复；校验值用于确认输入未改写
project = Path(__file__).resolve().parents[3] / 'IBD_ecotype_analysis'
validation = project / '06_validation'
root = validation / 'recovery_interpretation_20261010'
root.mkdir(exist_ok=True)
definitions_path = validation / 'bayesprism_recovery_readiness_20261009/53状态初步定义.tsv'
bulk_path = validation / 'GSE72819_recovery_review_20261009/53状态注释与恢复审阅.tsv'
ileum_path = validation / 'IBDverse_scRNA_recovery_20261009/47状态单细胞恢复审阅.tsv'
colon_path = validation / 'GSE214695_scRNA_recovery_20261009/53状态单细胞恢复审阅.tsv'
definitions = pd.read_csv(definitions_path, sep='\t')
bulk = pd.read_csv(bulk_path, sep='\t')
ileum = pd.read_csv(ileum_path, sep='\t')
colon = pd.read_csv(colon_path, sep='\t')
database = validation / 'bayesprism_recovery_readiness_20261009/database'
input_paths = [definitions_path, bulk_path, ileum_path, colon_path] + list(database.glob('*.json'))
input_hashes = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in input_paths}
assert len(definitions) == 53 and definitions.CellState.is_unique
print('输入状态：', len(definitions), len(bulk), len(ileum), len(colon))

# Step2. 区分论文方法、数据库与队列，不把参考思路写成已复现的计算
sources = pd.DataFrame([
    ['M1', '方法文献', 'Luca et al., Cell 2021', 'https://pmc.ncbi.nlm.nih.gov/articles/PMC8526411/', 'reference-guided recovery、marker对照及生态型共现证据层级；实际执行沿用官方代码', '当前代码30次打乱；本轮没有新增生态型共现置换'],
    ['M2', '方法文献', 'Steen et al., Cancer Cell 2021', 'https://doi.org/10.1016/j.ccell.2021.08.011', '借鉴图2–4以外部单细胞marker/已知表型解释状态、图5–6成员组合解释生态型的思路', '未复现fgsea、CytoTRACE、肿瘤克隆或临床预测分析'],
    ['D1', '数据库', 'MSigDB Human Hallmark：预先选定16集', 'https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp', '与冻结状态全部关联基因做集合重叠；保存命中基因、数量、集合大小及SHA256', '2026-10-09保存的JSON未注明发布版本；不能追认文件版本；非全库GSEA/ORA，无P值'],
    ['D2', '数据库', 'Human Protein Atlas：CSF3R单细胞/免疫细胞页面', 'https://www.proteinatlas.org/ENSG00000119535-CSF3R/single%2Bcell', '2026-10-10页面显示CSF3R的中性粒细胞相关表达，辅助核对Mast_S01谱系疑点', '单基因辅助证据；不直接改谱系，未全库扫描HPA'],
    ['Q1', '外部队列', 'GSE72819：UC基线结肠bulk', 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72819', '73例官方bulk recovery；状态胜出及既有成员支持描述', 'EUCALYPTUS而非RISK；胜出不等于谱系验证，不做临床关联'],
    ['Q2', '外部队列', 'IBDverse E-MTAB-16999：回肠scRNA', 'https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-16999', '作者大类label与冻结状态交叉；30供者、每供者×类型至多20细胞', '无Glial注释；限量抽样不代表组织整体组成'],
    ['Q3', '外部队列与文献', 'GSE214695；Macrophage and neutrophil heterogeneity... Nat Commun 2023', 'https://doi.org/10.1038/s41467-023-40156-6', '沿用作者annotation，显式大类映射后官方recovery；18供者、每供者×类型至多100细胞', '不重聚类或复现Milo/轨迹/空间分析；HC为非肿瘤区域；两barcode冲突隔离'],
    ['P1', '本项目描述方法', '预先指定诊断marker与供者等权汇总', '', '供者内计算检出率/平均CPM后等权汇总；比对作者标签，保留Unassigned', '部分marker参与投影；对应不是准确率，少供者保留不确定性']
], columns=['SourceID', 'Category', 'Source', 'URL', 'ActuallyUsed', 'Limit'])
sources.to_csv(root / '文献方法与数据库来源.tsv', sep='\t', index=False)
database_provenance = pd.read_csv(database / 'MSigDB来源记录.tsv', sep='\t')
for row in database_provenance.itertuples(index=False):
    assert hashlib.sha256((database / (row.GeneSet + '.json')).read_bytes()).hexdigest() == row.SHA256
database_provenance.to_csv(root / 'MSigDB原件来源与校验.tsv', sep='\t', index=False)
overlaps = pd.read_csv(validation / 'bayesprism_recovery_readiness_20261009/53状态_MSigDB定向面板重叠.tsv', sep='\t')
overlaps['FractionOfAssociatedGenes'] = overlaps.OverlapCount / overlaps.StateAssociatedGenes
overlaps['SourceID'] = 'D1'
overlaps.to_csv(root / '53状态_数据库基因重叠证据.tsv', sep='\t', index=False)

# Step3. 对齐恢复支持，不用Z大小自动排生物学可信度
metrics = ['CellState', 'InputCells', 'InputDonors', 'AssignedCells', 'AssignedDonors', 'Z', 'OwnMinusBestOther_InternalOnly']
states = definitions.merge(bulk[['CellState', 'WinnerSamples']], on='CellState', validate='one_to_one')
states = states.merge(colon[metrics].add_suffix('_Colon').rename(columns={'CellState_Colon': 'CellState'}), on='CellState', validate='one_to_one')
states = states.merge(ileum[metrics].add_suffix('_Ileum').rename(columns={'CellState_Ileum': 'CellState'}), on='CellState', how='left', validate='one_to_one')
states['Judgment'] = states.InitialProgram
updated_judgments = {
    'Mast_S01': '髓系/粒细胞样候选程序；典型CSF3R/CXCR2在作者Mast中缺乏表达支持，不确认炎症Mast亚型',
    'Mast_S03': '淋巴样混合程序，结直肠恢复支持不足，暂不定义Mast功能亚型',
    'Glial_S02': '吸收上皮/代谢样候选程序，在结直肠胶质中恢复支持不足',
    'Glial_S03': 'CCL21/补体共享候选程序，在结直肠胶质中恢复支持不足',
    'Glial_S05': '胶质身份相关候选程序，但独立结直肠中状态间marker区分不足',
    'Goblet_S01': '细胞周期相关杯状候选程序；仅部分细胞有周期marker，不等同全群增殖',
    'Goblet_S02': '屏障/转运与代谢候选程序，作者成熟杯状标签有所对应，尚非确定成熟亚型',
    'Enteroendocrine_S01': '损伤上皮/趋化炎症候选程序，结直肠仅3细胞且恢复Z不足，功能命名待定',
    'Enteroendocrine_S03': '屏障/离子转运候选程序，结直肠仅3细胞，不足以定义激素亚型'
}
states['Judgment'] = states.CellState.map(updated_judgments).fillna(states.Judgment)
states['SourceIDs'] = 'M1;M2;D1;Q1;Q2;Q3;P1'
states.loc[states.CellState.eq('Mast_S01'), 'SourceIDs'] += ';D2'
states['DefinitionStatus'] = '初步程序解释；保留原谱系/状态，未验证为最终细胞亚型'
states['EvidenceNote'] = '恢复统计、同源marker对照、作者表型对应分开评估；不是准确率'
states.loc[states.AssignedCells_Colon.eq(0), 'EvidenceNote'] += '；结直肠无分配细胞'
states.loc[states.AssignedCells_Ileum.eq(0), 'EvidenceNote'] += '；回肠无分配细胞'
states.loc[states.InputCells_Ileum.isna(), 'EvidenceNote'] += '；回肠无该谱系注释，不可评估'
states.loc[states.OwnMinusBestOther_InternalOnly_Colon.le(0), 'EvidenceNote'] += '；结直肠同源marker未高于最强其他状态'
states.loc[states.OwnMinusBestOther_InternalOnly_Ileum.le(0), 'EvidenceNote'] += '；回肠同源marker未高于最强其他状态'
states.to_csv(root / '53状态_初步判断及来源.tsv', sep='\t', index=False)
print(states[['CellState', 'AssignedCells_Colon', 'AssignedDonors_Colon', 'AssignedCells_Ileum', 'AssignedDonors_Ileum']].head().to_string(index=False))

# Step4. 作者亚型交叉：供者内计算状态组成，再等权汇总，未分配保留
cross_rows = []
for dataset, data_root, label in [('Colon_GSE214695', validation / 'GSE214695_scRNA_recovery_20261009', 'AuthorLabel'), ('Ileum_IBDverse', validation / 'IBDverse_scRNA_recovery_20261009', 'label')]:
    cells = pd.read_csv(data_root / 'raw_server/audit/cell_assignments_all.tsv', sep='\t')
    assert cells.ID.is_unique and label in cells
    counts = cells.groupby(['CellType', 'State', label, 'Donor']).size().rename('Cells').reset_index()
    denominators = cells.groupby(['CellType', 'State', 'Donor']).size().rename('StateCellsInDonor').reset_index()
    counts = counts.merge(denominators, on=['CellType', 'State', 'Donor'], validate='many_to_one')
    counts['FractionWithinStateDonor'] = counts.Cells / counts.StateCellsInDonor
    counts.to_csv(root / (dataset + '_作者亚型_供者中间表.tsv'), sep='\t', index=False)
    for (cell_type, state), group in cells.groupby(['CellType', 'State'], sort=True):
        state_counts = counts[counts.CellType.eq(cell_type) & counts.State.eq(state)]
        donor_number = group.Donor.nunique()
        for author_label, evidence in state_counts.groupby(label, sort=True):
            cross_rows.append({'Dataset': dataset, 'CellType': cell_type, 'State': state, 'AuthorSubtype': author_label, 'Cells': int(evidence.Cells.sum()), 'LabelDonors': evidence.Donor.nunique(), 'StateDonors': donor_number, 'DonorBalancedFractionWithinState': float(evidence.FractionWithinStateDonor.sum() / donor_number), 'SourceID': 'Q3' if dataset.startswith('Colon') else 'Q2'})
cross = pd.DataFrame(cross_rows)
cross.to_csv(root / '作者亚型与恢复状态_描述性对应.tsv', sep='\t', index=False)
top_labels = cross[cross.State.ne('Unassigned')].sort_values('DonorBalancedFractionWithinState', ascending=False).groupby(['Dataset', 'CellType', 'State'], sort=False).head(3)
top_labels.to_csv(root / '各状态前三作者亚型_仅描述.tsv', sep='\t', index=False)

# 附加诊断类型仍按作者标签保留；中性粒细胞用于谱系对照，不冒充已做模型recovery。
diagnostic = pd.read_csv(validation / 'GSE214695_scRNA_recovery_20261009/raw_server/input/diagnostic_marker_with_annotation.tsv', sep='\t', index_col=0)
marker_genes = [g for g in diagnostic.columns if g not in ['CellType', 'AuthorLabel', 'Donor', 'Stratum']]
mean_cpm = diagnostic.groupby(['CellType', 'Donor'])[marker_genes].mean().groupby('CellType').mean()
detection = diagnostic[marker_genes].gt(0).join(diagnostic[['CellType', 'Donor']]).groupby(['CellType', 'Donor']).mean().groupby('CellType').mean()
diagnostic_rows = []
for cell_type, group in diagnostic.groupby('CellType', sort=True):
    for gene in marker_genes:
        diagnostic_rows.append({'CellType': cell_type, 'Gene': gene, 'Cells': len(group), 'Donors': group.Donor.nunique(), 'DonorBalancedMeanCPM': mean_cpm.loc[cell_type, gene], 'DonorBalancedDetectionFraction': detection.loc[cell_type, gene], 'FormalRecoveryLineage': cell_type in set(states.CellType), 'SourceIDs': 'Q3;P1'})
diagnostic_table = pd.DataFrame(diagnostic_rows)
projection_overlap = pd.read_csv(validation / 'GSE214695_scRNA_recovery_20261009/raw_server/audit/diagnostic_marker_projection_overlap.tsv', sep='\t')
diagnostic_table = diagnostic_table.merge(projection_overlap.rename(columns={'DiagnosticGene': 'Gene'}), on=['CellType', 'Gene'], how='left', validate='one_to_one')
diagnostic_table.to_csv(root / '作者谱系_供者等权诊断marker.tsv', sep='\t', index=False)
colon_cells = pd.read_csv(validation / 'GSE214695_scRNA_recovery_20261009/raw_server/audit/cell_assignments_all.tsv', sep='\t')
colon_cells.groupby(['CellType', 'AuthorLabel', 'InitialState', 'State']).size().rename('Cells').reset_index().to_csv(root / '作者亚型_初始状态与保留状态.tsv', sep='\t', index=False)

# Step5. 每样本成员出现数；至少1细胞即记出现，受预算影响，不是共现检验
presence_rows = []
for dataset, data_root in [('Colon_GSE214695', validation / 'GSE214695_scRNA_recovery_20261009'), ('Ileum_IBDverse', validation / 'IBDverse_scRNA_recovery_20261009')]:
    cells = pd.read_csv(data_root / 'raw_server/audit/cell_assignments_all.tsv', sep='\t')
    official = pd.read_csv(data_root / 'raw_server/audit/ecotype_assignment.txt', sep='\t', index_col=0).set_index('ID')
    cells['CellState'] = cells.CellType + '_' + cells.State
    evaluated = set(states.loc[states.CellType.isin(cells.CellType.unique()), 'CellState'])
    for sample, sample_cells in cells.groupby('Sample', sort=True):
        present = set(sample_cells.loc[sample_cells.State.ne('Unassigned'), 'CellState'])
        for ecotype, members in states.groupby('Ecotype', sort=True):
            expected = set(members.CellState) & evaluated
            present_members = sorted(present & expected)
            presence_rows.append({'Dataset': dataset, 'Sample': sample, 'Donor': sample_cells.Donor.iloc[0], 'Ecotype': ecotype, 'EvaluableMembers': len(expected), 'PresentMemberStates': len(present_members), 'PresentMembers': ';'.join(present_members), 'OfficialAssignedEcotype': official.loc[sample, 'Ecotype'] if sample in official.index else 'Unassigned', 'SourceIDs': 'M1;M2;Q3' if dataset.startswith('Colon') else 'M1;M2;Q2'})
presence = pd.DataFrame(presence_rows)
presence.to_csv(root / '生态型_逐样本成员出现_描述性.tsv', sep='\t', index=False)
eco_bulk = pd.read_csv(validation / 'GSE72819_recovery_review_20261009/7生态型恢复与匹配discovery分布.tsv', sep='\t')
eco_rows = []
eco_names = {'E1': '炎症/干扰素与损伤相关候选组合', 'E2': '淋巴分化与CCL19/CCL21基质候选组合', 'E3': '核内/结构调控、基质与增殖相关混合组合', 'E4': '糖代谢/生物合成相关候选组合', 'E5': '应激/热休克与重塑相关候选组合', 'E6': '屏障/转运、基质及免疫调控混合组合', 'E7': '代谢/氧化还原及驻留样混合候选组合'}
for ecotype, members in states.groupby('Ecotype', sort=True):
    row = {'Ecotype': ecotype, 'InitialInterpretation': eco_names[ecotype], 'MemberStates': len(members), 'Members': ';'.join(members.CellState), 'MemberPrograms': ';'.join(members.CellState + ':' + members.InitialProgram), 'SourceIDs': 'M1;M2;D1;Q1;Q2;Q3', 'Limit': '按冻结成员描述程序组合；不证明细胞相互作用、空间共定位或人群生态型比例'}
    b = eco_bulk[eco_bulk.Ecotype.eq(ecotype)].iloc[0]
    row['BulkAssignedSamples73'] = int(b.RecoverySamples)
    row['BulkMedianMemberWinnerCount'] = float(b.MedianMemberWinnerCount)
    for suffix, dataset in [('Colon', 'Colon_GSE214695'), ('Ileum', 'Ileum_IBDverse')]:
        p = presence[presence.Dataset.eq(dataset) & presence.Ecotype.eq(ecotype)]
        assigned = p[p.OfficialAssignedEcotype.eq(ecotype)]
        row[suffix + 'OfficialAssignedSamples'] = len(assigned)
        row[suffix + 'MedianPresentMembersAmongAssigned'] = float(assigned.PresentMemberStates.median()) if len(assigned) else np.nan
        row[suffix + 'MembersWithAnyAssignedCells'] = int(members['AssignedCells_' + suffix].gt(0).sum())
        row[suffix + 'EvaluableMembers'] = int(members['InputCells_' + suffix].notna().sum())
    eco_rows.append(row)
ecotypes = pd.DataFrame(eco_rows)
ecotypes.to_csv(root / '7生态型_初步判断及来源.tsv', sep='\t', index=False)

# Step6. 保存证据卡；数据库重叠、恢复与最终亚型定义分开表述
lines = ['# Cell state与生态型初步解析', '', '2026-10-10；冻结Region3 CD+UC、BayesPrism presort模型。只读原始结果，给出候选程序解释；不修改状态编号、W或生态型成员，不做临床关联。', '', '## 实际用了什么方法与数据库', '', '- Cell 2021：官方reference-guided recovery及marker对照；当前代码30次打乱，名义Z参考1.96。', '- Cancer Cell 2021：借鉴外部单细胞与已知表型对照、成员组合解释生态型的思路。本轮未运行fgsea、CytoTRACE或临床预测模型。', '- MSigDB：复用2026-10-09保存的16个Human Hallmark集合，只做集合重叠；SHA256已核对。无正式GSEA/ORA、NES、富集P值或FDR。', '- HPA：CSF3R页面辅助核对Mast_S01谱系疑点；单基因不能确定完整谱系。', '- 作者表型：沿用IBDverse的label与GSE214695的annotation，不重聚类；供者等权计算状态内作者亚型组成。', '- 独立队列不等于所有指标独立：marker与H可能共享投影基因；表达对照不是准确率。', '', '详见文献方法与数据库来源.tsv，各SourceID对应URL、实际操作和限制。', '', '## 53状态证据卡', '']
for row in states.itertuples(index=False):
    support = f'结直肠 {int(row.AssignedCells_Colon)}细胞/{int(row.AssignedDonors_Colon)}供者；'
    support += '回肠不可评估' if pd.isna(row.InputCells_Ileum) else f'回肠 {int(row.AssignedCells_Ileum)}细胞/{int(row.AssignedDonors_Ileum)}供者'
    gene_sets = overlaps[overlaps.CellState.eq(row.CellState) & overlaps.OverlapCount.gt(0)].sort_values('OverlapCount', ascending=False).head(3)
    hallmark = '；'.join(f'{x.GeneSet} {x.OverlapCount}/{x.PanelGeneSetSize}个集合基因（例：{",".join(str(x.OverlapGenes).split(",")[:5])}）' for x in gene_sets.itertuples(index=False))
    author_rows = top_labels[top_labels.CellType.eq(row.CellType) & top_labels.State.eq(row.State)]
    author_evidence = '；'.join(f'{x.Dataset}: {x.AuthorSubtype} {x.DonorBalancedFractionWithinState:.0%}（{x.LabelDonors}供者）' for x in author_rows.itertuples(index=False))
    lines += [f'### {row.CellState}（{row.Ecotype}）', '', f'- 初步判断：{row.Judgment}；沿用此前程序判断支持级别“{row.ProgramConfidence}”，不等于统计显著性。', '- 主要基因：' + ', '.join(row.Top20Genes.split(',')[:10]) + '。', '- 数据库重叠：' + (hallmark or '定向16集无命中') + '；按命中数列示，非富集检验（D1）。', '- 外部支持：' + support + f'；bulk胜出 {int(row.WinnerSamples)}/73。', '- 作者亚型对应：' + (author_evidence or '没有分配细胞，无法对应') + '；供者等权、状态内比例，非准确率。', '- 限制：' + str(row.ReviewReason) + '；' + row.EvidenceNote + '。', '- 来源：' + row.SourceIDs + '。', '']
lines += ['## 7生态型证据卡', '']
for row in ecotypes.itertuples(index=False):
    lines += [f'### {row.Ecotype}：{row.InitialInterpretation}', '', '- 冻结成员：' + row.Members + '。', f'- bulk：{row.BulkAssignedSamples73}/73例分配；成员胜出数中位数{row.BulkMedianMemberWinnerCount:g}。', f'- 结直肠scRNA：{row.ColonOfficialAssignedSamples}/18样本官方分配；{row.ColonMembersWithAnyAssignedCells}/{row.ColonEvaluableMembers}成员至少一细胞恢复。', f'- 回肠scRNA：{row.IleumOfficialAssignedSamples}/30样本官方分配；{row.IleumMembersWithAnyAssignedCells}/{row.IleumEvaluableMembers}可评估成员至少一细胞恢复。', '- 限制：' + row.Limit + '。出现数以至少1细胞计，不是显著性阈值；少细胞成员另看供者。', '- 来源：' + row.SourceIDs + '。', '']
lines += ['## 下一步', '', '先审阅四个重点谱系及E1/E2。E3–E7保留候选组合名称；E7不能称回肠独占。正式功能富集或临床分析另设步骤，不把初步支持升级成最终细胞亚型。']
(root / '初步解析_逐状态与生态型证据卡.md').write_text('\n'.join(lines) + '\n', encoding='utf-8')
assert all(hashlib.sha256(p.read_bytes()).hexdigest() == input_hashes[str(p)] for p in input_paths)
pd.DataFrame([{'Input': p, 'SHA256Before': h, 'SHA256After': hashlib.sha256(Path(p).read_bytes()).hexdigest()} for p, h in input_hashes.items()]).to_csv(root / '只读输入校验.tsv', sep='\t', index=False)
(root / '解析检查摘要.json').write_text(json.dumps({'states': len(states), 'ecotypes': len(ecotypes), 'source_records': len(sources), 'hallmark_sets': len(database_provenance), 'input_hashes_unchanged': True, 'formal_enrichment_performed': False, 'clinical_analysis_performed': False, 'definition_changes': False}, ensure_ascii=False, indent=2), encoding='utf-8')
print('初步解析输出：', root)
