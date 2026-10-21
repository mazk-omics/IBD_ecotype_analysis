# 汇总53状态文献对应与候选证据，分类由明确生物学审阅决定，不用旧置信度代替新颖性。
# 所有结论仅针对旧版探索性结果；不能把同队列再分析、未匹配或未解决状态叫作Novel。
import os
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
from pathlib import Path
import json
import hashlib
import pandas as pd

# Step1. 复用已完成的分析与核实文献统计
root = Path('D:/IBD_ecotype/IBD_ecotype_analysis')
out = Path(__file__).resolve().parents[1]
old = root / 'analysis_reports/legacy_bayesprism_annotation/run_20261013_01'
base = root / '06_validation/all_states_lineage_audit_20261012'
master = pd.read_csv(old / 'CELL_STATE_ANNOTATION_MASTER.tsv', sep='\t').fillna('')
within = pd.read_csv(out / 'WITHIN_LINEAGE_STATE_COMPARISON.tsv', sep='\t')
summary = pd.read_csv(out / 'audit/ALL53_WITHIN_LINEAGE_SUMMARY.tsv', sep='\t').fillna('')
matches = pd.read_csv(out / 'audit/ALL_PUBLISHED_SIGNATURE_MATCH_STATISTICS.tsv', sep='\t').fillna('')
subtypes = pd.read_csv(base / 'tables/state_author_subtype_composition.tsv', sep='\t')
probe = pd.read_csv(out / 'audit/CANDIDATE_GENE_PROBE.tsv', sep='\t')
expr = pd.read_csv(base / 'tables/all_state_gene_expression.tsv.gz', sep='\t')
paired = pd.read_csv(base / 'tables/all_associated_genes_paired_donor_replication.tsv.gz', sep='\t')
contexts = pd.read_csv(out / 'audit/CANDIDATE_MSCCR_PATIENT_PAIRED_CONTEXT.tsv', sep='\t')
bounds = pd.read_csv(out / 'audit/CANDIDATE_COEXPRESSION_BOUNDS.tsv', sep='\t')
papers = pd.read_csv(out / 'audit/PUBLISHED_SIGNATURE_RESOURCE_INVENTORY.tsv', sep='\t').fillna('')
definitions = pd.read_csv(out / 'audit/PUBLISHED_SIGNATURE_DEFINITIONS.tsv', sep='\t').fillna('')

# Step2. 类别反映已知对应/区分性/未充分解释/无法解释，不表示细胞身份或新颖性已证实
established = {'B_S01': '已有干扰素应答B程序；STAT1+ Naive B signature对应不等于本状态必为naive B', 'CD4_T_S02': '已有IFN/免疫调节T程序；外部Martin免疫调节完整程序匹配，仍不等同确定Treg身份', 'CD8_T_S01': '已有干扰素应答/效应T程序；不以一个GZMB基因定义新细胞毒亚型', 'CD4_T_S03': '已有LEF1/TCF7低分化或记忆T程序；naive/CM等精细身份仍未分开', 'CD8_T_S02': '已有低分化/SELL相关T程序；不能凭SELL signature把作者CTL/TRGC2细胞改为naive T', 'Fibroblast_S02': 'CCL19/CCL21组织基质程序对应mLTo及ABCA8-fibroblast签名，非新TLS/FDC身份', 'Goblet_S01': '已有细胞周期相关杯状程序，含已知crypt MKI67+作者状态；整个S01不等于全群增殖', 'Monocyte_Macrophage_S01': '已发表M1/单核炎症marker有强匹配；同GSE214695数据源，粒细胞共享仍需谱系审阅'}
related = {'B_S02': 'B分化marker与已有B状态相关；组合区分未证明新的B亚型', 'B_S04': '生物合成/AREG与已有cycling B程序相关，但不能将代谢等同增殖', 'B_S05': 'HMOX1/IRF4候选应激程序，与B活化/普遍应激相关，非独立新亚型', 'CD4_T_S01': '与已知活化/驻留T程序相关，PDLIM2/IL32等有双队列方向支持，但驻留身份不确定', 'Endothelial_S01': '已知黏附/炎症血管程序相关；缺少足够双队列配对供者以定义特殊内皮亚型', 'Monocyte_Macrophage_S02': '偏已知M2/驻留谱系，P2RY13/CLEC10A组合值得研究；PLA2G2D不复现，DC-like替代未完全排除', 'Monocyte_Macrophage_S03': '偏IDA相关已知巨噬程序，RUNX3/RARA/GRN/GBA等区别复现；IDA匹配不能忽略，不定名TGFB驱动新亚型', 'Fibroblast_S03': '基质/ITGA11/ECM程序相关，多作者状态混合，缺足够第二队列供者配对', 'Fibroblast_S04': '外部Myofibroblast signature对应，但无充分双队列供者方向支持，FN1不能定义纤维化因果', 'Fibroblast_S05': '应激/氧化还原相关程序，有处理应激替代解释，非新生物学亚型', 'Goblet_S02': '成熟杯状/转运程序相关，作者成熟标签部分对应，不认为确定成熟亚型', 'Glial_S05': '与发表Glia2/神经胶质marker相关；单队列状态区分弱，不确认特殊胶质亚型', 'Mast_S06': 'Mast基础身份与免疫/骨架程序相关，未确认特殊TGFB功能', 'Pericyte_S04': '与已知Pericyte signature相关，CCL21共享读出有区别线索，父谱系与亚型仍需验证'}
undercharacterized = {'CD8_T_S03': '候选为JUND/LAPTM5/NR1H2伴随免疫骨架程序；13个Top50基因双队列供者方向稳定。IL23R并非核心。可能是已知TRM/效应或TRGC2组成轴，尚不能称新亚型', 'CD4_T_S04': '候选为ZNF831/CD6/SEPT9/FERMT3受体-骨架/黏附程序；6个Top50基因双队列方向稳定，跨记忆/Treg/ANXA1标签；已有状态混合和普遍T信号仍是替代解释'}
known_rows = []
for state in master.itertuples():
    category = 'Unresolved'
    reason = '当前细胞身份、可评价性或功能区分证据不足；文献无明确对应不能视为Novel。' + state.Uncertainty
    if state.CellState in established:
        category, reason = 'Established correspondence', established[state.CellState]
    if state.CellState in related:
        category, reason = 'Related but distinct', related[state.CellState]
    if state.CellState in undercharacterized:
        category, reason = 'Potentially undercharacterized', undercharacterized[state.CellState]
    m = matches[matches.CellState.eq(state.CellState) & matches.ComparisonScope.eq('same_parent_or_shared_program') & matches.Evaluable.eq(True)].copy()
    m = m.sort_values(['BH_Q_GlobalEvaluableMatches', 'OverlapN'], ascending=[True, False])
    top = m.head(3)
    sig = m[m.BH_Q_GlobalEvaluableMatches.lt(.05)]
    s = summary[summary.CellState.eq(state.CellState)].iloc[0]
    alt = matches[matches.CellState.eq(state.CellState) & matches.ComparisonScope.eq('explicit_alternative_lineage_audit') & matches.Evaluable.eq(True)].sort_values('BH_Q_GlobalEvaluableMatches').head(2)
    known_rows.append({'CellState': state.CellState, 'FrozenParent': state.CellType, 'FrozenEcotype': state.Ecotype, 'Category': category, 'CategoryScope': 'Provisional program-level correspondence; not novel subtype determination', 'BiologicalJudgment': reason, 'PriorCandidateAnnotation': state.CandidateAnnotation, 'PriorFunctionalConfidence_NotNoveltyEvidence': state.Confidence, 'TopPublishedCorrespondence': ';'.join(f'{r.Publication}:{r.PublishedStateOrProgram}[{r.Query};hits={r.OverlapN};q={r.BH_Q_GlobalEvaluableMatches:.4g}]' for r in top.itertuples()), 'TopCorrespondenceDOIs': ';'.join(top.DOI.drop_duplicates()), 'SignificantPublishedMatches': ';'.join(f'{r.Publication}:{r.PublishedStateOrProgram}:{r.Query}:q={r.BH_Q_GlobalEvaluableMatches:.4g}' for r in sig.itertuples()), 'AlternativeLineageSignatures': ';'.join(f'{r.PublishedStateOrProgram}:{r.Query}:q={r.BH_Q_GlobalEvaluableMatches:.4g}' for r in alt.itertuples()), 'Top50GeneN': s.Top50N, 'SharedPositiveGeneN': s.SharedPositiveGeneN, 'DoubleDatasetDonorStableGeneN': s.BothDatasetStableGeneN, 'DoubleDatasetDonorStableGenes': s.BothDatasetStableGenes, 'PublishedResourceCoverage': 'Full available exports Friedrich/Elmentaite/Garrido; full enumerated Martin programs; Smillie cluster-marker TableS2 not retrieved; mixed-age atlas not IBD-only', 'Independence': 'Colon GSE214695 and Ileum IBDverse used in marker selection/recovery; Garrido signature source overlaps Colon cohort', 'NovelCellSubtypeClaim': False, 'UnmatchedDoesNotProveNovelty': True})
known = pd.DataFrame(known_rows)
known.to_csv(out / 'KNOWN_STATE_CORRESPONDENCE.tsv', sep='\t', index=False)

# Step3. 优先候选与反向证据，身份审计和已知比较对象与程序候选分别标记
candidate_specs = [
    ('CD8_T_S03', 'P1_program', True, 'JUND/LAPTM5/NR1H2-associated immune-cytoskeletal program', '13 Top50基因在两队列供者方向稳定，多个核心基因广泛检出', 'IL23R仅约2%检出，正向供者比例低且回肠均差负；IL17A为0检出。回肠作者CD8+TRGC2+占88%，不可忽略细胞组成', '已知TRM/效应T连续轴或TRGC2相关细胞组成，而非新IL23-responsive CD8亚型', '13 stable genes; pooled distinction reproduced, though many genes broadly shared', 'Published Tc17/IL23R exists in Smillie; full TableS2 missing; do not claim novelty', '15/11 paired donors; IL23R positive fractions 4/15 and 2/11', 'CD8/T anchors detected, alpha-beta vs gamma-delta/TCR not resolved', 'Program hypothesis only; AP1 or NR1H2 activity not measured', '在相同作者CD8/TRGC2亚型内，该多基因程序仍区别于S01/S02/S04，并在多数独立供者中同细胞表达', '若控制作者亚型或TCR后差异消失，或独立供者核心程序方向不稳定，则为组成/连续状态而非独立程序'),
    ('Monocyte_Macrophage_S02', 'P1_known-state-refinement_and_DC_audit', True, 'P2RY13/CLEC10A-related macrophage program', 'P2RY13供者方向8/8、5/6正；CLEC10A在两队列方向6/8、5/6正且高检出；8个Top50双队列稳定基因', 'PLA2G2D结肠检出7.65%、回肠0；FLT3回肠33.3%、FCER1A50.8%有DC-like读出；核心MERTK/CD68未在本地审计基因表提供', '已知M2/驻留巨噬或cDC2/DC-like分化连续谱、父谱系混合；不称新脂质介质状态', 'P2RY13 and CLEC10A distinguish relative states; PLA2G2D not stable core', 'Garrido M2 externalTop50 q=0.0416 and author M2/M2.2/CD163-hi support known identity', '8/6 paired donors; CLEC10A rank102 included beyond Top50 for purposeful audit', 'Author macrophage annotations support, FLT3/FCER1A keep alternative lineage open', 'Reframe as resident/DC-like heterogeneity; no lipid mediator function assay', 'P2RY13与CLEC10A在具有多种排他性巨噬锚点的细胞内稳定共表达，且不由cDC2比例解释', '若核心程序仅出现于FLT3/CD1C/FCER1A DC细胞，或与已发表M2签名完全重合，则不是额外巨噬程序'),
    ('Monocyte_Macrophage_S03', 'P1_related-known-myeloid-program', True, 'RUNX3/RARA/GRN/GBA-associated myeloid program', '9 Top50基因双队列方向稳定；RUNX3和GRN高检出；相对S02分子轴可复现', 'Colon作者IDA占52.9%；IDA_vs_M1_UP签名匹配q=0.0435。TGFB1结肠仅6/10方向正，不是稳定核心', '已发表IDA/驻留巨噬内部连续谱，溶酶体或普遍免疫调节程序，不是新TGFB驱动亚型', '9 stable genes; state distinctions reproducible', 'Published IDA correspondence materially explains program; keep Related but distinct', '10/4 paired donors; no new independent cohort', 'Macrophage anchors and author labels coherent; no proof of new lineage', 'Study residual differentiation within known IDA/M2 continuum', '在匹配已知IDA/M2状态后，RUNX3/RARA/GRN/GBA组合仍有独立分子和功能差异', '若已知IDA/M2标签与签名充分解释该组合，或TGFB功能读出不存在，则撤回独立新程序假设'),
    ('CD4_T_S04', 'P2_program', True, 'ZNF831/CD6/FERMT3-associated receptor-cytoskeletal program', '6 Top50基因两队列供者方向稳定；CD6/FERMT3/SEPT9构成联合受体-骨架线索', '多个基因是广泛T细胞基础信号；作者ANXA1/Treg/memory等混合；完整Smillie/CD4DE未取得', '已知记忆/活化T或Treg混合及细胞大小/通用骨架差异', '6 stable genes; many DiscoveryTop50 shared positive across CD4 states', 'Current published exports do not exhaustively explain; missing references prevents novelty claim', '14/13 paired donors; borrowed marker cohorts', 'T anchors coherent, fine CD4 subtype mixed', 'Potentially undercharacterized component, not new T subtype', '匹配记忆/Treg作者亚型后，ZNF831/CD6/FERMT3联合程序仍区分S04，且不依赖IFN或细胞周期', '若差异被Treg/记忆比例、UMI深度或通用T信号解释，则不支持独立功能程序'),
    ('Goblet_S01', 'P1_known-program_refinement', True, 'Cell-cycle-positive component within Goblet_S01', 'MUC2近100%、TFF3=100%；MKI67检出8.63%/43.89%；完整Martin周期程序外部排名匹配', '全S01不全群增殖；回肠仅2匹配供者；两个周期基因联合检出下界为0；已知crypt MKI67+作者状态占18.5%', '已知cycling/immature goblet或隐窝位置/区域/采样混合；不视为新增殖杯状亚型', 'Cycle signal distinguishes component, not every assigned cell', 'Known cell-cycle program and existing crypt MKI67+ goblet explain much of signal', '8/2 paired donors; region differences and selection circularity', 'MUC2/TFF3 strongly coherent; progenitor/secretory doublet audit incomplete', 'Quantify mixture of cycling and noncycling goblet components', 'S01包含同一细胞内MUC2与至少2周期读出共表达的子群，该子群比S02多，但整体S01不均一增殖', '若周期阳性仅来自前体/双细胞或与已知crypt MKI67+标签完全重合，则不提出新的杯状程序'),
    ('Fibroblast_S02', 'Known_program_comparator', False, 'CCL19/CCL21 organizer-like stromal program', 'Discovery与两版外部marker均匹配mLTo；两版外部匹配Friedrich ABCA8 fibroblasts', '第二队列供者配对不足；没有空间TLS或FDC结构证据', '已有淋巴组织基质/ABCA8 fibroblast程序与发育读出共享', 'Several shared chemokine/ECM markers; not unique TLS identity', 'Established correspondence at program level', '4/0 matched donors', 'Fibroblast background coherent, fine stromal labels mixed', 'Positive control for known-state recovery', '用已知组织基质signature作为阳性对照检验候选筛选方法', '若不能复现该明确对应，应先修正比较背景和signature覆盖，而非宣称其他状态Novel'),
    ('Fibroblast_S03', 'P3_low-priority_context-dependent', False, 'ECM/ITGA11-associated fibroblast program', '结肠11配对供者及ECM/ITGA11多基因线索', '回肠仅2配对供者、Top50无双队列稳定基因；多个S2a/S2b/myofibroblast作者标签混合', '已有基质/机械/成纤维分化连续谱和抽样差异', 'Some local distinction; cross-cohort power inadequate', 'Related known stromal programs, incomplete full Kinchen signatures', '11/2 matched donors', 'Broad fibroblast coherent; fine subtype heterogeneous', 'Do not promote while stronger candidates exist', '有足够新独立供者时，ECM/ITGA11程序在同亚型内保持区分', '若仅由myofibroblast或S2a/S2b比例解释，则属于已知组成变化')]
identity_audits = ['Enteroendocrine_S01', 'Mast_S01', 'Mast_S07', 'Glial_S02']
candidate_rows = []
for state, priority, selected, focus, positive, negative, alternative, specificity, correspondence, reproducibility, coherence, interpretability, hypothesis, falsifier in candidate_specs:
    k = known[known.CellState.eq(state)].iloc[0]
    candidate_rows.append({'CellState': state, 'FrozenParent': k.FrozenParent, 'Ecotype': k.FrozenEcotype, 'Priority': priority, 'SelectedForProgramFollowup': selected, 'Category': k.Category, 'ProgramFocus': focus, 'ProgramSpecificity': specificity, 'LiteratureCorrespondence': correspondence, 'ExternalReproducibility': reproducibility, 'LineageCoherence': coherence, 'BiologicalInterpretability': interpretability, 'PositiveEvidence': positive, 'CounterEvidence': negative, 'AlternativeExplanations': alternative, 'FalsifiableHypothesis': hypothesis, 'FalsificationCondition': falsifier, 'IndependentValidationStatus': 'NOT_DONE; current external cohorts participated in marker selection/recovery', 'NovelSubtypeClaim': False})
for state in identity_audits:
    k = known[known.CellState.eq(state)].iloc[0]
    m = master[master.CellState.eq(state)].iloc[0]
    candidate_rows.append({'CellState': state, 'FrozenParent': k.FrozenParent, 'Ecotype': k.FrozenEcotype, 'Priority': 'P0_identity_audit_before_biology', 'SelectedForProgramFollowup': False, 'Category': 'Unresolved', 'ProgramFocus': m.CandidateAnnotation, 'ProgramSpecificity': 'Not interpretable before identity audit', 'LiteratureCorrespondence': 'Cross-lineage differences are not novelty evidence', 'ExternalReproducibility': 'Recovery or detection does not settle origin', 'LineageCoherence': m.Uncertainty, 'BiologicalInterpretability': 'Identity-first; no new mechanism claim', 'PositiveEvidence': m.BiologicalRationale, 'CounterEvidence': m.Uncertainty, 'AlternativeExplanations': 'Wrong parent, mixed transcriptional source, doublet/ambient signal, rare-cell instability, reference limitations', 'FalsifiableHypothesis': '先在同一细胞确认父谱系+程序+阴性候选谱系', 'FalsificationCondition': '若程序只定位于另一个已知谱系则不解释为新父谱系状态', 'IndependentValidationStatus': 'NOT_DONE', 'NovelSubtypeClaim': False})
candidates = pd.DataFrame(candidate_rows)
candidates.to_csv(out / 'UNDERCHARACTERIZED_STATE_CANDIDATES.tsv', sep='\t', index=False)

# Step4. 重点候选逐基因证据和缺测均保存；诊断不只看单marker
diagnostic_panels = {'Monocyte_Macrophage_S02': ['P2RY13', 'CLEC10A', 'PLA2G2D', 'CSF1R', 'MAFB', 'MRC1', 'MERTK', 'CD68', 'FLT3', 'CD1C', 'FCER1A', 'CLEC9A', 'LILRA4'], 'Monocyte_Macrophage_S03': ['RUNX3', 'RARA', 'GRN', 'GBA', 'NRROS', 'TGFB1', 'PPARD'], 'Goblet_S01': ['MUC2', 'TFF3', 'MKI67', 'TOP2A', 'FOXM1', 'SPDEF', 'OLFM4', 'LGR5'], 'CD8_T_S03': ['CD8A', 'CD8B', 'CD3D', 'TRAC', 'JUND', 'LAPTM5', 'NR1H2', 'IL23R', 'RORC', 'IL17A', 'IL17F', 'TRDC', 'TRGC2'], 'CD4_T_S04': ['CD4', 'CD3D', 'ZNF831', 'CD6', 'FERMT3', 'SEPT9', 'FOXP3', 'TCF7', 'IL7R'], 'Fibroblast_S03': ['COL1A1', 'DCN', 'ITGA11', 'ECM1', 'ADAMTS7', 'ACTA2']}
diagnostic_rows = []
diagnostic_panels.update({'Enteroendocrine_S01': ['CHGA', 'CHGB', 'SYP', 'INSM1', 'NEUROD1', 'EPCAM', 'DUOX2', 'NOS2', 'LCN2'], 'Mast_S01': ['TPSAB1', 'TPSB2', 'KIT', 'MS4A2', 'CPA3', 'HDC', 'S100A9', 'CSF3R', 'CXCR2', 'FPR1', 'OSM'], 'Mast_S07': ['TPSAB1', 'KIT', 'MS4A2', 'CPA3', 'HDC', 'CSF1R', 'MS4A7', 'STAB1', 'LST1', 'TYROBP', 'FLT3'], 'Glial_S02': ['SOX10', 'PLP1', 'S100B', 'MPZ', 'GFAP', 'HNF4G', 'CUBN', 'KHK', 'SLC5A9', 'EPCAM']})
discovery_ranks = pd.read_csv(root / '06_validation/paper_aligned_recovery_characterization_20261012/00_audit/prior_discovery_associated_genes.tsv', sep='\t')
for state, genes in diagnostic_panels.items():
    for dataset in ['Colon_GSE214695', 'Ileum_IBDverse']:
        for gene in genes:
            e = expr[expr.CellState.eq(state) & expr.Dataset.eq(dataset) & expr.Gene.eq(gene)]
            p = paired[paired.CellState.eq(state) & paired.Dataset.eq(dataset) & paired.Gene.eq(gene)]
            row = {'CellState': state, 'Dataset': dataset, 'Gene': gene, 'MeasurementStatus': 'NOT_IN_LOCAL_AUDIT_UNIVERSE' if len(e) == 0 else ('MEASURED_SOURCE' if e.SourceMeasured.iloc[0] else 'NOT_SOURCE_MEASURED')}
            if len(e):
                row.update(e.iloc[0][['DetectionFraction', 'MeanLog2CPM', 'AssignedCells', 'AssignedDonors']].to_dict())
                row['AssignedStateCoverageStatus'] = e.Status.iloc[0]
                if e.AssignedCells.iloc[0] == 0:
                    row['MeasurementStatus'] = 'NOT_EVALUABLE_NO_ASSIGNED_STATE_CELLS'
            rank = discovery_ranks[discovery_ranks.CellState.eq(state) & discovery_ranks.Gene.eq(gene)]
            row['DiscoveryAssociatedRank'] = rank.WithinStateRank.iloc[0] if len(rank) else None
            if len(p):
                row.update(p.iloc[0][['MatchedDonors_Min5CellsEach', 'MeanPairedLog2CPM_Difference', 'FractionMatchedDonors_Positive']].to_dict())
            diagnostic_rows.append(row)
diagnostics = pd.DataFrame(diagnostic_rows)
diagnostics.to_csv(out / 'audit/CANDIDATE_LINEAGE_AND_CORE_MARKER_REVIEW.tsv', sep='\t', index=False)
cards = ['# Candidate State Evidence Review', '', '仅旧版冻结探索性结果。候选优先级不是新颖性分数；已发表状态能解释的部分必须先承认。当前所有外部表达/DE/recovery支持为marker选择队列再分析，非新独立验证。', '', '同细胞共表达当前未直接计算；下文Fréchet界限来自同一群体供者等权边际检出率，只是约束联合检出范围，不是观测联合细胞数或供者内联合频率。', '']
for candidate in candidates.itertuples():
    cards += [f'## {candidate.CellState}：{candidate.ProgramFocus}', '', f'优先级：{candidate.Priority}；类别：{candidate.Category}。', '', '| 审阅维度 | 判断 |', '|---|---|', f'| Program specificity | {candidate.ProgramSpecificity} |', f'| Literature correspondence | {candidate.LiteratureCorrespondence} |', f'| External reproducibility | {candidate.ExternalReproducibility} |', f'| Lineage coherence | {candidate.LineageCoherence} |', f'| Biological interpretability | {candidate.BiologicalInterpretability} |', '', f'正向证据：{candidate.PositiveEvidence}。', '', f'反向证据：{candidate.CounterEvidence}。', '', f'替代解释：{candidate.AlternativeExplanations}。', '', f'可证伪假设：{candidate.FalsifiableHypothesis}。', '', f'撤回/证伪条件：{candidate.FalsificationCondition}。', '', '核心/诊断基因的外部数据：', '', '| 队列 | 基因 | 检出率 | 匹配供者 | 配对均差 | 正方向供者比例 |', '|---|---|---:|---:|---:|---:|']
    d = diagnostics[diagnostics.CellState.eq(candidate.CellState)]
    for row in d.fillna('').itertuples():
        cards.append(f'| {row.Dataset} | {row.Gene} ({row.MeasurementStatus}) | {row.DetectionFraction} | {row.MatchedDonors_Min5CellsEach} | {row.MeanPairedLog2CPM_Difference} | {row.FractionMatchedDonors_Positive} |')
    cards += ['', '发表signature对应（完整背景与全部非显著结果在audit统计表）：']
    m = matches[matches.CellState.eq(candidate.CellState) & matches.Evaluable.eq(True)].sort_values(['BH_Q_GlobalEvaluableMatches', 'OverlapN'], ascending=[True, False]).head(5)
    for r in m.itertuples():
        cards.append(f'- {r.Publication} / {r.PublishedStateOrProgram}；{r.Query}；命中{r.OverlapN}/{r.QueryInBackgroundN}，signature可测{r.SignatureInBackgroundN}/{r.PublishedSignatureN}；q={r.BH_Q_GlobalEvaluableMatches:.4g}；[DOI](https://doi.org/{r.DOI})。{r.EvidenceIndependence}。')
    cards += ['', '作者亚型组成（供者等权，不是一对一真值）：']
    labels = subtypes[subtypes.CellState.eq(candidate.CellState)].sort_values('DonorBalancedFraction', ascending=False).groupby('Dataset').head(3)
    for r in labels.itertuples():
        cards.append(f'- {r.Dataset}：{r.AuthorLabel}={r.DonorBalancedFraction:.3f}。')
    b = bounds[bounds.CellState.eq(candidate.CellState)]
    if len(b):
        cards += ['', '同细胞联合检出界限（不是直接观测）：']
        for r in b.itertuples():
            cards.append(f'- {r.Dataset}，{r.Panel}：供者等权联合检出下界={r.LowerBoundDonorBalancedJointDetection:.4f}，上界={r.UpperBoundDonorBalancedJointDetection:.4f}。')
    clinical = contexts[contexts.CellState.eq(candidate.CellState)]
    if len(clinical):
        cards += ['', 'MSCCR患者配对背景（保留原H，未重新归一；非独立验证）：']
        for r in clinical.itertuples():
            q = 'NA' if pd.isna(r.BH_Q_AllEvaluablePatientComparisons) else f'{r.BH_Q_AllEvaluablePatientComparisons:.4g}'
            delta = 'NA' if pd.isna(r.MedianInflamedMinusNonInflamed) else f'{r.MedianInflamedMinusNonInflamed:.4g}'
            cards.append(f'- {r.Disease}/{r.Region3}：{r.PairedPatients}对患者；Inflamed−Non-inflamed中位差={delta}，BH q={q}；{r.Status}。')
        cards.append('- 主要控制同患者、疾病和肠段；年龄/治疗/时间等额外混杂未建模。NMF H相对分量有组成和被过滤状态残余影响，下降不等于保护/抗炎。MM_S03的UC/Colon结果对求解器数值零处理敏感，未加容差q≈0.061，不能称稳健显著。')
    cards += ['', '必要后续：在资源允许且独立授权后，只取目标父谱系和核心基因的逐细胞检查点；匹配作者亚型、供者、肠段、检测深度/周期并联合阴性谱系。独立供者或原位联合验证才增加新验证层级。本轮不下载新队列、不重跑全recovery。', '']
(out / 'CANDIDATE_STATE_EVIDENCE_REVIEW.md').write_text('\n'.join(cards), encoding='utf-8')

# Step5. 按Cell2021 Results逻辑写科学解释，先已知，再区别，最后未解决
counts = known.Category.value_counts().to_dict()
result = f'''# Known-State Correspondence and Potentially Undercharacterized Cell-State Programs in IBD

2026-10-13；冻结旧版BayesPrism–EcoTyper 53 states / 7 ecotypes。与Joint20 benchmark和未来CIBERSORTx主线独立。类别分布：{json.dumps(counts, ensure_ascii=False)}。这些是程序层面的探索性分类，不是最终细胞亚型或新颖性判定。

## 重新识别已有IBD相关程序

依据Luca et al., Cell 2021的marker/外部恢复/已知程序联合表征思路，我们比较了完整发表导出资源而非只挑几个marker：Friedrich2021补充表6全部上调cluster marker及补充表9完整IL-1β刺激DE；Elmentaite2021补充表7完整导出（每cluster最多50marker，混合发育/健康背景）；Martin2019原文明确列举的5个完整T细胞程序；Garrido-Trigo2023补充图5d各200个M1/M2/IDA marker及图5e完整成对DE。后者与GSE214695 recovery数据源相同，不提供额外独立验证。Smillie完整cluster marker TableS2、Martin完整cluster DE表、Corridoni/Parikh完整表尚未获得，不能宣称穷尽IBD文献。

B_S01重新呈现ISG程序，并与STAT1+ B参考signature对应，但不因此认定其为naive B。CD4_T_S02同时呈现炎症与免疫调节程序，其外部marker匹配Martin完整免疫调节程序；这不是Treg身份的充分证明。CD8_T_S02与SELL相关低分化T程序对应，但作者细胞标签包含CTL/TRGC2，精细分化身份仍须核对。

Fibroblast_S02的CCL19/CCL21相关组合与mLTo及ABCA8-fibroblast签名对应，属于已有淋巴组织基质程序的重新识别。组织基质程序相似不能推导TLS、FDC或空间结构。Goblet_S01的外部marker与已知细胞周期程序匹配，且回肠有既有crypt MKI67+作者标签；不构成新增殖杯状亚型。MM_S01与发表M1/炎症单核签名显著对应，需保留粒细胞共享与同队列来源限制。

## 同父细胞类型内的可复现区别

Discovery assignment给每基因一个胜出状态，并不说明基因只在该状态表达。读取全部冻结gene_info初始列可见多数Top50在同父类型多个状态呈正效应。例如CD8_T_S03有49/50、CD4_T_S04有49/50、MM_S02有40/50、MM_S03有42/50基因同时在至少一个其他冻结状态为正（描述阈值0.1）。真正值得关注的是方向和幅度的组合差异，而非“独占marker”命名。

现有供者内目标vs其余已分配状态对比显示，CD8_T_S03有13、MM_S03有9、MM_S02有8、CD4_T_S04有6个Top50基因同时满足两队列≥3匹配供者、平均差>0、≥70%供者方向正。该规则是本次透明筛查，不是原论文判定门槛或显著性检验；不足供者不记为不复现。B_S03也有6个方向稳定基因，但主要为泛核内加工/调节，未因此晋升生物学新程序。

髓系S02与S03分别显示P2RY13/CLEC10A方向和RUNX3/RARA/GRN/GBA等方向区别。已获得的Garrido-Trigo完整导出为S02提供M2相关解释，为S03提供IDA-relative-to-M1上调程序解释，不能忽略这些已发表状态。S03的TGFB1在结肠仅60%匹配供者方向正，不足以作为稳定致纤维化核心。谱系来源可能是已知巨噬状态连续谱，而非新亚型。

CD4父类型中，S02的IFN/免疫调节、S03的低分化和S04的受体-骨架读出提供不同分子轴。Fibroblast的组织基质、ECM和应激轴亦可描述，但全部5个状态在本次严格两队列供者方向筛查下均无Top50同时通过，主要受回肠配对供者不足影响；不应把该阴性筛查写成程序错误，也不宜优先宣称新基质亚型。

## 值得深入的候选程序，及必须修正的三个预设线索

**CD8_T_S03：值得研究的是JUND/LAPTM5/NR1H2伴随的免疫骨架组合，而不是已证实的IL23R状态。** IL23R检出率在两队列约1.92%和2.07%；结肠4/15供者方向正，回肠2/11正且平均差为负；IL17A在当前审计中为0检出。Smillie2019已发表IL23R/IL17相关CD8程序，因此IL23R本身也不是新颖性证据。当前更稳的13基因组合值得在同一作者CD8/TRGC2亚型内验证。回肠S03约88%来自作者CD8+TRGC2+标签；这个组成替代解释尚未排除。没有TCR证据不能直接把这些细胞叫作γδ T。

**MM_S02：修正为P2RY13/CLEC10A相关已知巨噬状态细化候选。** PLA2G2D在结肠仅7.65%检出、回肠S02为0，不能维持“跨肠段稳定脂质介质状态”命名。P2RY13方向在8/8与5/6供者为正，CLEC10A在6/8与5/6为正；CLEC10A在discovery关联列表排名102，因明确DC-like排查目的纳入，不把它伪装成Top50。CSF1R支持巨噬背景，回肠FLT3/FCER1A读出则要求进一步排除DC-like/已知cDC2程序。当前MERTK/CD68等部分诊断基因缺本地审计结果，缺测不能视为阴性。

**Goblet_S01：明确区分“含周期阳性杯状细胞”与“全群增殖状态”。** MUC2检出率约99.55%/99.63%、TFF3为100%，MKI67为8.63%/43.89%。同群体边际值对MUC2+MKI67的供者等权联合检出给出Fréchet范围：结肠8.18%–8.63%，回肠43.52%–43.89%。这支持周期阳性成分与杯状身份相容，仍不是直接测得的联合细胞数。MUC2+MKI67+TOP2A下界为0，完整多基因周期程序的逐细胞共表达未证实；回肠状态对照只有2合格供者。最有价值的下一步是量化已知crypt cycling与非cycling成分，而非提出新杯状亚型。

此外，**CD4_T_S04的ZNF831/CD6/FERMT3相关受体-骨架程序**有6个双队列供者稳定基因，值得作为较低优先级程序候选；它跨记忆/Treg/ANXA1作者标签出现，完整既有CD4 signature覆盖不全且通用T信号仍可解释。我们仅将CD8_T_S03和CD4_T_S04列为“Potentially undercharacterized”候选组件，没有证明它们是新细胞亚型。MM_S02/MM_S03仍归Related but distinct，Goblet_S01作为已知程序细化对象。

## 可证伪假设与患者背景

首要假设：CD8_T_S03的JUND/LAPTM5/NR1H2组合在匹配作者亚型/TCR、肠段、供者和检测深度后仍区别于同父类型其他状态；若控制组成后差异消失，则撤回独立程序假设。CD4_T_S04的联合受体-骨架特征需在匹配记忆/Treg背景后成立；若仅为通用T信号/UMI差异则不支持额外生物学程序。MM_S02需在具有联合巨噬身份和DC阴性/阳性对照的同一细胞确认P2RY13/CLEC10A；若已知M2/cDC2足以解释，则不提出新程序。MM_S03需在匹配已知IDA/M2状态后仍有残余分子/功能差异。

利用冻结MSCCR2025 H审计表，7个预定候选/比较对象进行了同患者、同疾病、同Region3的炎症配对；重复活检先取中位数，CD/Colon有41对、UC/Colon有46对，其余背景<10对不作推断。CD8_T_S03及MM_S03在结肠炎症活检中相对分量降低；CD4_T_S04在UC/Colon升高、CD/Colon没有同样支持。上述结果是同discovery cohort背景关联；H为旧NMF相对系数且过滤后未重新归一，不能把减少推断为抗炎或保护。年龄/治疗/时间等额外混杂未建模。MM_S03的UC/Colon q值对数值零处理跨越0.05（容差版≈0.0475、原值≈0.0610），不作稳健显著结论。

## 暂不具备可靠解释的状态

Enteroendocrine_S01、Mast_S01/S07及Glial_S02优先身份审计；损伤上皮、粒细胞/髓系或吸收上皮marker偏离父谱系不构成Novel证据。广泛核内加工状态（如B_S03）、稀少细胞、缺失Glial参考及弱供者支持状态也保留Unresolved。高恢复Z、可测基因或既有功能注释置信度均不能替代新颖性证据。

## 核心科学结论

除已有IBD程序的重新识别，本轮最可信的增量是：跨已知作者标签的可重复多基因差异轴、S02/S03髓系内部区别，以及Goblet_S01中区域依赖的周期阳性成分。CD8_T_S03和CD4_T_S04的程序组件值得独立验证，但目前没有证据足以宣布新的细胞亚型或功能机制。明确撤回IL23R和PLA2G2D作为稳定核心的解释，并承认M2/IDA和crypt cycling等已知替代，较强行扩展新名称更接近当前数据。

## 文献来源

- [Luca et al., Cell 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8526411/)，DOI 10.1016/j.cell.2021.09.014：表征思路参照，非癌症状态直接迁移。
- [Martin et al., Cell 2019](https://doi.org/10.1016/j.cell.2019.08.008)：完整原文程序列表复用，本次没有完整cluster DE表。
- [Smillie et al., Cell 2019](https://doi.org/10.1016/j.cell.2019.06.029)：IL23R相关CD8既有解释；全文BioC复用，完整TableS2未取得。
- [Friedrich et al., Nature Medicine 2021](https://doi.org/10.1038/s41591-021-01520-5)：真实补充表6和9已取得。
- [Elmentaite et al., Nature 2021](https://doi.org/10.1038/s41586-021-03852-1)：真实补充表7已取得，发育/健康背景区别保留。
- [Garrido-Trigo et al., Nature Communications 2023](https://doi.org/10.1038/s41467-023-40156-6)：补充图5d/e完整可用导出；与GSE214695外部队列重合。
- [Corridoni et al., Nature Medicine 2020](https://doi.org/10.1038/s41591-020-1003-4)、[Parikh et al., Nature 2019](https://doi.org/10.1038/s41586-019-0992-y)、[Kinchen et al., Cell 2018](https://doi.org/10.1016/j.cell.2018.08.067)：真实文献身份/资源位置核实，但完整signature未取得，不假称完成全面对照。
'''
(out / 'CELL_STATE_DISCOVERY_INTERPRETATION.md').write_text(result, encoding='utf-8')
print('Correspondence categories:', counts, flush=True)
print('Candidate rows:', len(candidates), '; program followup:', int(candidates.SelectedForProgramFollowup.sum()), flush=True)
