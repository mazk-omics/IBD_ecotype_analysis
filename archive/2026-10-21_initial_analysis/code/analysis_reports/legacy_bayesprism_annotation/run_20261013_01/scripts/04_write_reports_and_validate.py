# 复用本次已保存的小表，补全逐状态和生态型报告及交付检查。
# 不再读取压缩全转录组DE，不运行任何恢复或模型计算。
import os
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
from pathlib import Path
import hashlib
import json
import subprocess
import pandas as pd

# Step1. 更新本次人工审阅字段与真实表型单位
root = Path('D:/IBD_ecotype/IBD_ecotype_analysis')
out = Path(__file__).resolve().parents[1]
base = root / '06_validation/paper_aligned_recovery_characterization_20261012'
lineage = root / '06_validation/all_states_lineage_audit_20261012'
master = pd.read_csv(out / 'CELL_STATE_ANNOTATION_MASTER.tsv', sep='\t').fillna('')
curation = pd.read_csv(out / 'audit/STATE_CURATION.tsv', sep='\t').fillna('')
master = master.drop(columns=[c for c in curation.columns if c != 'CellState']).merge(curation, on='CellState', validate='one_to_one')
marker = pd.read_csv(out / 'CELL_STATE_MARKER_EVIDENCE.tsv', sep='\t', low_memory=False)
discovery = marker[marker.EvidenceType.eq('Discovery_gene_info')]
all_discovery = pd.read_csv(base / '00_audit/prior_discovery_associated_genes.tsv', sep='\t')
phenotype = pd.read_csv(base / '04_state_phenotypes/state_phenotype_enrichment.tsv', sep='\t')
external = pd.read_csv(out / 'CELL_STATE_EXTERNAL_SUPPORT.tsv', sep='\t').fillna('')
ora = pd.read_csv(out / 'CELL_STATE_FUNCTIONAL_ENRICHMENT.tsv', sep='\t').fillna('')
literature = pd.read_csv(out / 'CELL_STATE_LITERATURE_MATCHING.tsv', sep='\t').fillna('')
for index, match in literature.iterrows():
    literature.loc[index, 'MatchingLevel'] = 'MULTIGENE_TOP20_PROGRAM_CANDIDATE' if match.DiscoveryTop20HitN >= 2 and match.FrozenLineageCompatible else ('CROSS_LINEAGE_PROGRAM_ONLY' if match.DiscoveryAssociatedHitN > 0 and not match.FrozenLineageCompatible else ('WEAK_ASSOCIATED_PANEL_OVERLAP' if match.DiscoveryAssociatedHitN >= 2 else ('SINGLE_GENE_INSUFFICIENT' if match.DiscoveryAssociatedHitN == 1 else 'NO_PANEL_OVERLAP')))
literature.to_csv(out / 'CELL_STATE_LITERATURE_MATCHING.tsv', sep='\t', index=False)
recovery = pd.read_csv(base / '01_state_recovery/state_recovery_evidence_by_dataset.tsv', sep='\t').fillna('')
completeness = pd.read_csv(out / 'audit/STATE_EVIDENCE_COMPLETENESS.tsv', sep='\t').fillna(0)
for index, state in master.iterrows():
    genes = set(all_discovery.loc[all_discovery.CellState.eq(state.CellState), 'Gene'])
    master.loc[index, 'CandidateTFObserved'] = ';'.join(g for g in state.CandidateTFPanel.split(';') if g in genes)
    master.loc[index, 'CandidateSurfaceObserved'] = ';'.join(g for g in state.CandidateSurfacePanel.split(';') if g in genes)
    master.loc[index, 'HistologyProgramGenesObserved'] = ';'.join(g for g in state.HistologyProgramPanel.split(';') if g in genes)
    ph = phenotype[phenotype.CellState.eq(state.CellState) & phenotype.Unit.eq('BiologicalSampleDominant') & phenotype.PositiveZgt1p96.eq(True)]
    master.loc[index, 'PhenotypeSampleDominantPositiveZ'] = ';'.join(f'{r.Dataset}:{r.Feature}:{r.Target}:Z={r.Z:.3g}' for r in ph.itertuples())
    lit = literature[literature.CellState.eq(state.CellState) & literature.MatchingLevel.eq('MULTIGENE_TOP20_PROGRAM_CANDIDATE')]
    master.loc[index, 'LiteratureCandidatePrograms'] = ';'.join(f'{r.Publication}:{r.PublishedProgram}[{r.Top20HitGenes}]' for r in lit.itertuples())
master.to_csv(out / 'CELL_STATE_ANNOTATION_MASTER.tsv', sep='\t', index=False)
eco_prior = pd.read_csv(base / '06_ecotype_features/7_ecotypes_evidence_summary.tsv', sep='\t').fillna('')
composition = master[['Ecotype', 'CellState', 'CellType', 'CandidateAnnotation', 'Confidence', 'BiologicalRationale', 'Uncertainty', 'DiscoveryTop20', 'ExternalTop10_DiseaseGroup', 'ExternalTop10_Region3', 'ORA_GlobalFDRlt0p05', 'CandidateTFObserved', 'CandidateSurfaceObserved', 'HistologyProgramGenesObserved', 'ValidationAnchorPanel', 'Colon_AssignedCells', 'Colon_AssignedDonors', 'Ileum_AssignedCells', 'Ileum_AssignedDonors']].merge(eco_prior, on='Ecotype', validate='many_to_one')
composition.to_csv(out / 'ECOTYPE_FUNCTIONAL_COMPOSITION.tsv', sep='\t', index=False)
discrete = pd.read_csv(base / '06_ecotype_features/ecotype_discrete_features.tsv', sep='\t')
continuous = pd.read_csv(base / '06_ecotype_features/ecotype_continuous_features.tsv', sep='\t')
region = pd.read_csv(base / '06_ecotype_features/ecotype_region_KW.tsv', sep='\t')
coassociation = pd.concat([pd.read_csv(base / f'05_ecotype_coassociation/{dataset}_ecotype_coassociation.tsv', sep='\t') for dataset in ['Colon_GSE214695', 'Ileum_IBDverse']], ignore_index=True)
pd.concat([discrete.assign(EvidenceScope='Discovery_bulk_discrete'), continuous.assign(EvidenceScope='Discovery_bulk_Hallmark'), region.assign(EvidenceScope='Discovery_Region3_KW'), coassociation.assign(EvidenceScope='Pooled_external_coassociation')], ignore_index=True, sort=False).to_csv(out / 'audit/ECOTYPE_CONTEXT_AND_ASSOCIATION_EVIDENCE.tsv', sep='\t', index=False)

# Step2. 逐状态给出数值证据和验证候选，禁止表型反向命名
lines = ['# 旧版 BayesPrism 53 Cell States 与 7 Ecotypes 的探索性注释', '', '日期：2026-10-13。范围：冻结MSCCR_CORE12_CTSPEC_REGION_IBD旧版53状态/7生态型。与Joint20 fraction benchmark独立；不构成未来CIBERSORTx discovery基准。', '', 'High/Moderate/Low/Unresolved衡量候选功能解释的可置信度，不是细胞身份分类准确率。所有原ID和父类型不变。High仅用于CD8_T_S01的干扰素程序和CD8_T_S02的低分化程序；不表示细亚型已验证。', '', '恢复Z、表达检出、marker方向、供者内复现、LOOCV和文献对应分别报告。既定父类型内的投影依赖模型基因，不能作为独立谱系真值。ORA是冻结NMF候选空间内的条件富集，不是通路活性或全转录组GSEA。', '', '## 逐状态证据卡', '']
for state in master.itertuples():
    lines += [f'### {state.CellState}（{state.Ecotype}；父类型 {state.CellType}）', '', f'**候选解释：{state.CandidateAnnotation}。置信度：{state.Confidence}。**', '', f'依据：{state.BiologicalRationale}。限制：{state.Uncertainty}。', '', f'Discovery Top20：{state.DiscoveryTop20}。完整Top50和MaxFC保存在marker证据表；MaxFC不是外部Wilcoxon log2FC。', '', f'外部疾病版Top10：{state.ExternalTop10_DiseaseGroup or "无可复用完整排名"}。Region3版Top10：{state.ExternalTop10_Region3 or "无可复用完整排名"}。', '', f'n4两版Top10交集：{state.N4_Top10Intersection or "不可评价"}；Disease-only：{state.N4_DiseaseOnlyTop or "无/不可评价"}；Region-only：{state.N4_RegionOnlyTop or "无/不可评价"}。']
    de = marker[marker.CellState.eq(state.CellState) & marker.EvidenceType.eq('Existing_external_DE') & marker.Scope.eq('Dataset') & marker.Gene.isin(state.DiscoveryTop20.split(','))].copy()
    de = de.drop_duplicates(['Dataset', 'Gene']).sort_values(['p_val_adj', 'log2FC'], ascending=[True, False]).head(6)
    de_text = '; '.join(f'{r["Dataset"]}:{r["Gene"]} log2FC={r["log2FC"]:.3g},Bonferroni P={r["p_val_adj"]:.3g},检出={r["pct.1"]:.3g}/{r["pct.2"]:.3g}' for _, r in de.iterrows())
    lines += ['', 'Top20基因的外部DE示例（按调整P排序，保留负方向）：' + (de_text or '没有可复用状态间统计结果；不能补成P=1或无marker') + '。', '', '外部恢复与供者支持：']
    for rec in recovery[recovery.CellState.eq(state.CellState)].itertuples():
        z = f'{rec.OfficialZ:.3g}' if isinstance(rec.OfficialZ, (float, int)) else 'NA'
        lines.append(f'- {rec.Dataset}：{rec.AssignedCells}细胞/{rec.AssignedDonors}供者；官方Z={z}；状态={rec.Status}。Z为全队列范围，P未提供，原恢复仅30次置换。')
    lines += ['', f'供者内marker方向：Colon {state.Colon_PairedStatus}，匹配供者={state.Colon_MatchedDonors_Min5CellsEach}；Ileum {state.Ileum_PairedStatus}，匹配供者={state.Ileum_MatchedDonors_Min5CellsEach}。仅描述性，不是新的供者水平DE显著性。', '', f'LOOCV：DiseaseGroup {state.LOOCV_ConcordantFolds_DiseaseGroup}/{state.LOOCV_EvaluableFolds_DiseaseGroup}可评价折吻合；Region3 {state.LOOCV_ConcordantFolds_Region3}/{state.LOOCV_EvaluableFolds_Region3}。分母为可评价折；NA不作失败。', '', '疾病–部位分层：']
    for support in external[external.CellState.eq(state.CellState)].itertuples():
        lines.append(f'- {support.Context}：{support.AssignedCellsInContext}细胞/{support.AssignedDonorsInContext}供者，≥5分配细胞供者={support.DonorsWithGE5AssignedCells}；{support.ContextStatus}；marker={support.JointContextMarkerStatus}；现有DE来源={support.MarkerDEAvailableScope}（{support.MarkerDESourceRegions}）。')
    lines += ['', f'功能富集（全局BH q<0.05）：{state.ORA_GlobalFDRlt0p05 or "无；不能据此否定候选程序"}。', '', f'作者亚型分布（供者等权、全队列）：{state.PooledAuthorSubtypeSummary or "未评价"}。作者标签不是一对一状态真值。', '', f'表型背景：生物样本dominant分析正Z：{state.PhenotypeSampleDominantPositiveZ or "无正向Z>1.96结果"}；discovery肠段Fisher q<0.05：{state.DiscoveryRegionFisherQlt0p05 or "无"}。不以这些关联决定身份。']
    lit = literature[literature.CellState.eq(state.CellState) & literature.MatchingLevel.eq('MULTIGENE_TOP20_PROGRAM_CANDIDATE')].sort_values(['DiscoveryTop20HitN', 'DiscoveryAssociatedHitN'], ascending=False).head(3)
    lines += ['', '文献程序对应（候选面板匹配，非亚型等价）：']
    for match in lit.itertuples():
        lines.append(f'- {match.Publication}：{match.PublishedProgram}，关联基因命中{match.DiscoveryAssociatedHitN}/{match.PanelN}，Top20命中{match.DiscoveryTop20HitN}；命中={match.HitGenes}；[DOI {match.DOI}](https://doi.org/{match.DOI})。')
    if len(lit) == 0:
        lines.append('- 当前核实过的短面板在Discovery Top20无至少2基因同谱系对应；全关联基因弱重叠仍保存在TSV，不用于强行命名或虚构文献亚型。')
    lines += ['', f'候选TF/转录辅因子（仅关联基因中实有）：{state.CandidateTFObserved or "无稳定候选"}；候选表面蛋白基因：{state.CandidateSurfaceObserved or "无稳定候选"}；组织学程序读出：{state.HistologyProgramGenesObserved or "无稳定候选"}。', '', f'组织学联合身份锚点建议：{state.ValidationAnchorPanel}。锚点是待验证panel，不声称均为本状态阳性；需在同一细胞确认父谱系+程序读出及阴性谱系，RNA表达不等同可用抗体/蛋白阳性。', '']

# Step3. 深入E1/E3/E7，同时显式保留其他生态型矛盾
eco_interpretation = {
    'E1': ('炎症/干扰素与损伤相关跨谱系候选组合', 'B/T成员有多个ISG；内皮SELE/ICAM1和周细胞CXCL1/2、PLAU/SERPINE1描述血管黏附与趋化重塑；髓系成员有CXCL1/2/3和S100A9。成员层面的证据与整体炎症Hallmark方向相容，但Enteroendocrine_S01仅少量细胞，Mast_S01有显著粒细胞样矛盾，Glial_S01只在单一可评价队列得到有限支持。', 'Martin2019 GIMATS提供“炎症免疫+基质”比较框架，但此模型没有IgG plasma状态，不能命名为GIMATS。Friedrich2021的IL-1基质-中性粒细胞程序也仅为生物学参照：E1没有冻结Fibroblast成员，Mast_S01不等于已验证中性粒细胞；不能宣称IL-1驱动、药物抵抗或真实通讯。E1的联合共关联与成员身份可信度必须分开。'),
    'E2': ('淋巴分化/淋巴组织基质样候选组合', 'BACH2/EBF1、TCF7/LEF1与Fibroblast_S02的CCL19/CCL21程序支持淋巴组织样解释。', 'Glial_S03、Mast_S03和Pericyte_S04身份/程序支持较弱，外部联合共关联Z均未达到1.96。没有空间聚集和B/T区结构，不能确认TLS或结构成熟的淋巴滤泡。'),
    'E3': ('基质黏附/结构调节与基础免疫程序的共现组合', '结构成员Fibroblast_S03、Endothelial_S03和Pericyte_S03分别具有ECM/ITGA11、ELN/LTBP4和BCAM/NOTCH相关读出；免疫成员显示受体信号、白细胞骨架或核内加工；Goblet_S01另有细胞周期读出，仅在部分细胞出现。', '两队列较高共关联Z支持该冻结组合的复现，但B_S03/ILC_S05的核内加工不能自动升级为独特功能，Enteroendocrine_S02、Glial_S05/S06细身份未确定。不能把E3统一命名为增殖、修复、纤维化或NOTCH生态型。整体Apical/WNT/Myogenesis等Hallmark与结构重塑可相容，仍不能分派给全部成员；共享上游尺度/基础程序是必须保留的替代解释。'),
    'E4': ('代谢/生物合成与AREG共享候选组合', 'B_S04及ILC_S04显示AREG与糖代谢/生物合成读出，Mast_S02有代谢程序。', 'Colon共关联较弱、Ileum负Z明显，且回肠缺少完整可定义成员。AREG不能独立证明ILC2或组织修复；不掩盖跨队列不一致。'),
    'E5': ('应激/氧化还原相关异质候选组合', 'B_S05/Fibroblast_S05有HMOX1和应激读出；Glial_S04为热休克，Mast_S05只呈基础加工。', 'Colon共关联有正向线索，Ileum支持不足；热休克可能有处理因素。整体OXPHOS关联不能推断所有成员氧化磷酸化增加，也不能定义保护性生态型。'),
    'E6': ('上皮屏障/转运与混合免疫调节候选组合', 'Goblet_S02和Enteroendocrine_S03有CLDN8/TRPM6转运读出，Fibroblast_S04显示FN1/受体程序，CD8_T_S03有IL23R/JUND。', 'B_S06、Mast_S04、Mast_S07未解决，两个队列联合共关联接近零。整体负炎症关联只描述discovery背景；不证明生态型为抗炎/稳态或独立修复网络。'),
    'E7': ('代谢/氧化还原背景与CD4活化或驻留样程序的弱支持组合', 'Endothelial_S02/Pericyte_S02共享PTGR1、CD36、NDRG1等代谢应激读出；Fibroblast_S01有OAT/FMO1/ECM2；CD4_T_S01的CD69/KLRB1/IL32具有独立于基质代谢的免疫程序。', 'Glial_S02的HNF4G/CUBN/KHK/SLC5A9偏吸收上皮，结肠仅3细胞且恢复Z为负，回肠无胶质输入；必须保留身份未解决。Fibroblast_S01的n4两版回肠LOOCV亦不同，其他结构成员供者内支持弱。两队列共关联Z均低，不支持把E7命名为已验证的稳态/抗炎生态型。它更适合作为需要谱系与独立联合复现验证的探索性组合；CD4成员的驻留样程序不能由整体负Hallmark强行解释。')}
lines += ['## 逐生态型解释', '', '统计共现不等于物理接触、空间邻近或细胞通讯。所有bulk关联均来自同一路线discovery，存在重复活检及共同表达影响；这不是独立临床验证。', '']
eco_summary_rows = []
for eco in eco_prior.itertuples():
    members = master[master.Ecotype.eq(eco.Ecotype)]
    label, rationale, caveat = eco_interpretation[eco.Ecotype]
    eco_summary_rows.append({'Ecotype': eco.Ecotype, 'CandidateInterpretation': label, 'MemberN': len(members), 'MemberConfidenceCounts': ';'.join(f'{k}:{v}' for k, v in members.Confidence.value_counts().items()), 'UnresolvedMembers': ';'.join(members.loc[members.Confidence.eq('Unresolved'), 'CellState']), 'Rationale': rationale, 'CriticalLimit': caveat})
    lines += [f'### {eco.Ecotype}：{label}', '', rationale, '', caveat, '', f'冻结成员：{eco.FrozenMembers}。', '', '| 父类型/State | 候选程序 | 置信度 |', '|---|---|---|']
    lines.extend(f'| {r.CellType}/{r.CellState} | {r.CandidateAnnotation} | {r.Confidence} |' for r in members.itertuples())
    lines += ['', '现有联合恢复/共关联：']
    for rec in coassociation[coassociation.Ecotype.eq(eco.Ecotype)].itertuples():
        lines.append(f'- {rec.Dataset}：共关联Z={rec.Z:.3g}；输入谱系成员{rec.InputLineageMembers}/{rec.FrozenMembers}；可定义参考成员{rec.DefinedReferenceMembers}/{rec.FrozenMembers}；有效置换={rec.ValidPermutations}。不补写未知方向的P或显著性星号。')
    d = discrete[discrete.Ecotype.eq(eco.Ecotype)]
    lines += ['', '| discovery关联 | 效应：GroupA−GroupB | BH q |', '|---|---:|---:|']
    lines.extend(f'| {r.Feature}: {r.GroupA}−{r.GroupB} | {r.MeanDifference:.4g} | {r.Q:.3g} |' for r in d.itertuples())
    kw = region[region.Ecotype.eq(eco.Ecotype)].iloc[0]
    lines += ['', f'Region3总体Kruskal–Wallis：统计量={kw.Statistic:.3g}，q={kw.Q:.3g}。无组均值/事后对比时不指定偏向哪个肠段。', '', '整体bulk Hallmark关联（仅q<0.05前5；完整数据在audit中；不是成员功能赋值）：']
    h = continuous[continuous.Ecotype.eq(eco.Ecotype) & continuous.Q.lt(.05)].sort_values('Q').head(5)
    for r in h.itertuples():
        lines.append(f'- {r.Feature}：均值差={r.MeanDifference:.4g}，q={r.Q:.3g}；生态型优势活检与其他优势活检比较。')
    lines += ['', '后续优先验证：先验证低置信度/未解决成员的父谱系与程序共表达，再在未参与定义的供者和空间数据中评价该完整组合；不得据共现补写配体–受体机制。', '']
pd.DataFrame(eco_summary_rows).to_csv(out / 'audit/ECOTYPE_INTERPRETATION_SUMMARY.tsv', sep='\t', index=False)
sources = pd.read_csv(out / 'audit/VERIFIED_LITERATURE_SOURCES.tsv', sep='\t')
lines += ['## 可核实文献', '']
lines.extend(f'- {r.Publication}：[ {r.Title} ]({r.URL})；DOI：{r.DOI}。{r.Verification}。' for r in sources.itertuples())
(out / 'CELL_STATE_AND_ECOTYPE_INTERPRETATION.md').write_text('\n'.join(lines), encoding='utf-8')

# Step4. 保存图源，不因绘图依赖缺失安装或改变环境
ora.pivot(index='CellState', columns='GeneSet', values='FDR_Global2650').reindex(master.CellState).to_csv(out / 'figures/FUNCTIONAL_HEATMAP_GLOBAL_FDR_SOURCE.tsv', sep='\t')
composition.groupby(['Ecotype', 'Confidence']).size().unstack(fill_value=0).to_csv(out / 'figures/ECOTYPE_MEMBER_CONFIDENCE_SOURCE.tsv', sep='\t')
composition.groupby(['Ecotype', 'CellType']).size().unstack(fill_value=0).to_csv(out / 'figures/ECOTYPE_PARENT_COMPOSITION_SOURCE.tsv', sep='\t')
marker_plot_parts = []
for cell_type in master.CellType.unique():
    path = lineage / f'source_data/{cell_type}_state_marker_plot_source.tsv'
    frame = pd.read_csv(path, sep='\t')
    frame['FrozenParentCellType'] = cell_type
    frame['SourcePath'] = str(path)
    marker_plot_parts.append(frame)
plot_source = pd.concat(marker_plot_parts, ignore_index=True)
plot_source.to_csv(out / 'figures/MARKER_DETECTION_HEATMAP_SOURCE.tsv', sep='\t', index=False)
(out / 'figures/FIGURE_CONTRACT_AND_STATUS.md').write_text('''# 图件合同与状态

核心问题：哪些冻结状态的功能程序具有外部marker支持，以及这些不同置信度的成员怎样组成7个生态型？

绘图偏好检查返回Python；现有两套本机Python均缺matplotlib，未安装软件。遵循nature-figure的Missing runtime/package rule，绘图状态BLOCKED_DEPENDENCY；图源已保存，渲染/碰撞/文字QA均NOT_RUN，不能作为论文最终图。

预备三个独立图：父类型内marker检出率热图（供者等权、Colon/Ileum分别显示，缺测为NA）；53×50 Hallmark ORA全局FDR热图（显示全部状态，不显著结果保留）；E1–E7父类型与成员置信度构成图（成员数，不是组织细胞比例）。各图是quantitative grid；不把低支持状态隐藏，不将基因检出率替换为状态恢复准确率。

导出目标：可编辑PDF/SVG和PNG预览。所有刻度和图例≥5 pt，供者等权中心值无推断性误差条。最终需要在可用Python环境完成PDF文字与碰撞审计及人工视觉检查；未执行多panel时对齐检查为NOT_APPLICABLE。
''', encoding='utf-8')

# Step5. 明确PASS/PARTIAL/BLOCKED，保存实际覆盖与方法限制
without_de = completeness.loc[completeness.DatasetsWithSavedDE.eq(0), 'CellState'].tolist()
both_de = int(completeness.DatasetsWithSavedDE.eq(2).sum())
one_de = int(completeness.DatasetsWithSavedDE.eq(1).sum())
conf = master.Confidence.value_counts().to_dict()
methods = f'''# 注释方法与隔离审计

日期：2026-10-13；旧版探索性注释，非未来CIBERSORTx主路线标准。

## P0 隔离与资源

根目录D:/IBD_ecotype并非Git仓库；实际仓库在IBD_ecotype_analysis/main。首次只读git status发现.gitignore及4个既有R脚本修改、大量未跟踪文件，保存audit/git_status_before.txt；未执行reset/clean/stash/checkout/commit/push。

本次所有脚本、日志、结果写入analysis_reports/legacy_bayesprism_annotation/run_20261013_01。旧版冻结模型和既有结果只读。本地首次内存约15.7 GiB、可用约4.9 GiB、16逻辑CPU、CPU采样负载14%。发现已有ssh.exe PID18136；不附接、不发送指令、不打开服务器连接。远端实时余量无法确认。

benchmark脚本只读显示MSCCR正式阶段n.cores=8，外部pseudo-bulk正式阶段n.cores=12（smoke=3）；这是配置值，不是实时资源监测。本任务单线程，BLAS/OMP/NumExpr限制1线程，未读取大型表达矩阵、未下载队列、未重新拟合或运行恢复、未修改Conda或安装包。一次流式读取既有68 MB压缩DE；其余仅小表。哈希检查34个本次来源一致，DE文件大小/mtime复核一致。benchmark正在写的日志和输出不作冻结校验，也不声称验证了远端执行状态。

## 原文方法对应

[Luca et al., Cell 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8526411/)，DOI 10.1016/j.cell.2021.09.014：借鉴已知程序、显著marker、已发表scRNA状态、TF/表面候选及外部恢复的证据组合。S4链接指向NIHMS1740510-supplement-10.xlsx；本次二进制内容访问遇浏览器验证，未假称逐列复现S4。当前是IBD适配，不是原文癌症亚型的直接迁移，也不是CIBERSORTx上游数值复现。

## P1 来源与完整性

冻结53状态、7成员关系来自paper_aligned_recovery_characterization_20261012/00_audit/frozen_states.tsv，经与7_ecotypes_evidence_summary的成员集合逐一验证；无增删或改ID。Discovery marker来自prior_discovery_associated_genes及gene_info，保留Top50、MaxFC和初始状态；MaxFC不重标为外部log2FC。每个父类型的初始gene_info哈希匹配NMF_selected_model_manifest，背景包括之后被过滤掉的初始状态。

现有外部全转录组来源DE可复用：{both_de}状态有两队列结果，{one_de}状态仅一队列，{len(without_de)}状态无保存的全队列状态间DE：{';'.join(without_de)}。这是实际可评价性，不是缺省53/53“完整DE”。Seurat3.1.3 Wilcoxon结果经min.pct=0.1、原生lnFC阈值0.05预筛；不足5细胞状态排除且要求至少2可评价状态。avg_logFC为ln，log2FC=avg_logFC/ln(2)；p_val_adj原生Bonferroni，不冒充BH FDR。缺行可能是预筛或不可评价，未补P。没有重新执行DE。

marker证据表是长表：Discovery Top50、两版所有fold的Top50排名、这些候选基因的既有DE方向和效应、以及每状态每队列显著正DE Top50。全折排名及n1–n7原件只读保留；audit内记录来源和版本差异。独立全转录组未预筛完整ranking不足，未做pre-ranked GSEA。

## n4历史定义核查

实际脚本/结果两版为DiseaseGroup=CD/UC/Control与Region3=Ileum/Colon/Rectum；不是CD_Ileum/CD_Colon/UC_Colon联合n4。用户描述的联合背景与保存实现之间的差异明确记录，不重定义或重跑n4。原两版同样保留，48/53都有完整Top10；中位Top10交集9，排名Spearman中位0.962。回肠B_S03和Fibroblast_S01在两版LOOCV重要差异，不能据两折成功率宣称某版更优。全部实际差异保存在N4_EXISTING_COMPARISON。

## P2 功能与文献

既有53_states_pathway_summary仅集合重叠，不含P/NES，本次不冒充显著富集。使用冻结MSigDB Hallmark v7.4的50集做轻量单侧超几何ORA；每父类型背景为哈希确认的initial_gene_info全初始状态基因，选集为相应冻结状态MaxFC>0的关联基因。保留背景N、选集N、原集合/背景内集合大小、命中基因及方向。主FDR全53×50=2650个检验的BH；另给每状态50检验的BH作敏感性列。全局FDR<0.05共63项；这是NMF候选空间条件富集，不能代表全转录组无偏功能激活。R phyper+p.adjust独立复核P和全局Q，最大绝对差均<1e-10；结果见ORA_INDEPENDENT_R_QC。

GO/Reactome正式库未在当前审阅来源发现可追溯状态富集，暂不下载大型库、不以Hallmark替代GO/Reactome。未做GSEA和新的TF regulon分析。TF/表面候选仅为人工审阅且实际关联基因中观察到的候选表达，不推断调控活性。组织学锚点建议另列，不能混同已观察marker。

在线核实Luca2021、Martin2019、Friedrich2021、Smillie2019、Elmentaite2021、Garrido-Trigo2023的真实题名/DOI与可访问主文本。文献短面板来自既有passage-indexed提取及新核对段落；保持人工选择/完整短列表标识。每状态×面板记录0命中、单基因不足、跨谱系共享或多基因候选，非正式亚型分类器、非完整signature GSEA。没有强行把所有状态配到某篇文章的亚型；无法匹配的状态在报告明确保留。

## P3 外部支持

使用既有逐细胞metadata和冻结assignment，为每状态保存3背景共159行：父类型输入细胞/供者，分配细胞/供者/样本，≥5分配细胞供者、状态间DE可评价性。NO_INPUT_LINEAGE、NO_ASSIGNED_CELLS、低细胞和描述性可评价分别记录。官方Z/供者配对/作者亚型/LOOCV均是原全队列范围，不重标为背景专属统计。

CD_Colon仅Colon、CD_Ileum仅Ileum，原DiseaseGroup DE可精确对应背景；UC疾病DE混合Colon/Rectum，仅可作为UC_Colon proxy，不能确认UC_Colon marker支持。只汇总现有Top20基因方向+Bonferroni<0.05，多基因标准≥2，缺少显著基因不等于无生物学程序。8个精确背景行达到该严格标准；其余保留局限。需UC_Colon专属DE时留待有资源且授权的后续任务，不在当前补跑。

## P4 Ecotypes与图件

53行成员功能构成完整覆盖E1–E7。逐状态功能证据先确定，随后结合原共关联Z、有效成员数、bulk疾病/炎症、Region3 KW和Hallmark；不从bulk反推成员功能。E1/E3有两队列较高共关联，E7弱；E4回肠负共关联、E5队列不一致、E2/E6弱支持均保留。共现不是物理接触或通讯，未补LR或全局“生态型检测”新检验。

已保存marker检出率、53×50功能FDR和E1–E7父类型/置信度图源，绘图脚本供后续复用。nature-figure保存偏好为Python，两个可用Python都缺matplotlib；根据该技能Missing runtime/package rule，渲染暂缓，未安装或用其他后端替代。图件QA未执行，不能交付为论文终稿。

## 判断与下一步

候选注释分布：{json.dumps(conf, ensure_ascii=False)}。High只针对2个CD8功能程序；Moderate可作谨慎Results候选描述；Low仅展示探索性线索；Unresolved保留ID。Mast_S01/Mast_S07/Glial_S02优先做双谱系锚点+程序RNA或蛋白同细胞验证。E1优先排除粒细胞/损伤上皮归入错误父类型；E3优先确认基础加工是否具有状态特异性；E7优先独立联合复现。

所有判定是本次审阅的候选草稿，不是旧模型“正确性”的回溯认证，也不限制新版模型的细胞类型、状态数、成员或功能解释。
'''
(out / 'ANNOTATION_METHODS_AND_AUDIT.md').write_text(methods, encoding='utf-8')
status_rows = [
    ('P0', 'PASS', '本任务全部写入隔离目录；未连接/控制服务器；无环境安装或现有脚本修改', '远端实时资源未知，仅审阅配置；不验证benchmark是否成功完成'),
    ('P1', 'PASS', '53状态候选注释、Top50/DE/n4/LOOCV/谱系证据和置信度实际交付', '不是53状态均有可评价DE，也不是最终亚型身份确认'),
    ('P2', 'PARTIAL', '2650 Hallmark条件ORA+独立R数值复核，6篇真实文献及1325状态面板比较', 'GO/Reactome、未预筛GSEA、S4二进制逐列复现未完成'),
    ('P3', 'PARTIAL', '3背景159行实际支持、精确/混合DE范围及不可评价状态分开', 'UC_Colon专属DE缺失；背景专属恢复置换和供者DE未执行'),
    ('P4', 'PARTIAL', 'E1–E7全部成员构成和解释；E1/E3/E7深入审阅；图源已保存', 'Python matplotlib缺失，新图渲染和图件QA未完成'),
    ('Figure rendering', 'BLOCKED', '已保存图源和脚本', '本机既有Python缺matplotlib，按隔离规则不安装'),
    ('Final subtype identity', 'PARTIAL', '23状态具有High/Moderate候选程序，16 Low，14 Unresolved', '无新的独立身份验证；重点疑点不强行定名')]
pd.DataFrame(status_rows, columns=['Module', 'Status', 'Completed', 'RemainingOrLimit']).to_csv(out / 'PASS_PARTIAL_BLOCKED.tsv', sep='\t', index=False)

# Step6. 验证交付维度、分母、成员和哈希，保留Git审阅状态
frozen = pd.read_csv(base / '00_audit/frozen_states.tsv', sep='\t')
assert len(master) == 53 and set(master.CellState) == set(frozen.CellState)
assert len(composition) == 53 and composition.Ecotype.nunique() == 7
assert len(external) == 159 and not external.duplicated(['CellState', 'Context']).any()
assert len(ora) == 2650 and ora.FDR_Global2650.between(0, 1).all()
assert not external[external.Context.eq('UC_Colon')].MarkerDEMatchesJointContext.any()
assert set(master.loc[master.Confidence.eq('High confidence'), 'CellState']) == {'CD8_T_S01', 'CD8_T_S02'}
for eco in eco_prior.itertuples():
    assert set(eco.FrozenMembers.split(';')) == set(composition.loc[composition.Ecotype.eq(eco.Ecotype), 'CellState'])
source_checks = pd.read_csv(out / 'audit/INPUT_SHA256_BEFORE.tsv', sep='\t')
assert all(hashlib.sha256(Path(row.Path).read_bytes()).hexdigest() == row.SHA256 for row in source_checks.itertuples())
git_result = subprocess.run(['git', '-C', str(root), '-c', 'core.quotepath=false', 'status', '--porcelain=v1'], capture_output=True, text=True, encoding='utf-8', check=True)
(out / 'audit/git_status_after.txt').write_text(git_result.stdout, encoding='utf-8')
required = ['CELL_STATE_ANNOTATION_MASTER.tsv', 'CELL_STATE_MARKER_EVIDENCE.tsv', 'CELL_STATE_FUNCTIONAL_ENRICHMENT.tsv', 'CELL_STATE_LITERATURE_MATCHING.tsv', 'CELL_STATE_EXTERNAL_SUPPORT.tsv', 'ECOTYPE_FUNCTIONAL_COMPOSITION.tsv', 'ANNOTATION_METHODS_AND_AUDIT.md', 'CELL_STATE_AND_ECOTYPE_INTERPRETATION.md']
assert all((out / name).stat().st_size > 0 for name in required)
summary = {'RequiredFilesPresent': len(required), 'States': len(master), 'Ecotypes': composition.Ecotype.nunique(), 'ConfidenceCounts': conf, 'ExistingDEBothDatasets': both_de, 'ExistingDEOneDataset': one_de, 'ExistingDENoSavedContrast': without_de, 'StatesWithGlobalSignificantHallmark': ora.loc[ora.FDR_Global2650.lt(.05), 'CellState'].nunique(), 'PhenotypeSampleDominantStatesWithPositiveZ': int(master.PhenotypeSampleDominantPositiveZ.ne('').sum()), 'FiguresRendered': 0, 'FigureStatus': 'BLOCKED_DEPENDENCY_matplotlib', 'SourceHashesUnchanged': len(source_checks)}
(out / 'audit/FINAL_DELIVERY_CHECK.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps(summary, ensure_ascii=False, indent=2))
