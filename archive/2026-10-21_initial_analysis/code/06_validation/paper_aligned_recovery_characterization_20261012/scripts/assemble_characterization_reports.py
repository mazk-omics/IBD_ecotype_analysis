# 汇总真实计算结果，给全部53状态和7生态型生成初步证据卡；未完成/不明确步骤不假造结果。
# 服务器ibd_refaudit环境或本机cibersortx_review环境运行；仅写本次分析目录。
# 读取中间表逐步合并，报告保留BayesPrism上游差异、供者依赖和方法阻断。
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

# Step1. 读入冻结定义和已独立验收结果；marker任务未完成时报告明确保留RUNNING
out = Path(__file__).resolve().parents[1]
states = pd.read_csv(out / '00_audit/frozen_states.tsv', sep='\t')
prior = pd.read_csv(out / '00_audit/prior_frozen_states.tsv', sep='\t')
genes = pd.read_csv(out / '00_audit/discovery_gene_assignment.tsv', sep='\t')
recovery = pd.read_csv(out / '01_state_recovery/state_recovery_evidence_by_dataset.tsv', sep='\t')
coassoc = pd.read_csv(out / '05_ecotype_coassociation/ecotype_coassociation_1000permutations.tsv', sep='\t')
discrete = pd.read_csv(out / '06_ecotype_features/ecotype_discrete_features.tsv', sep='\t')
continuous = pd.read_csv(out / '06_ecotype_features/ecotype_continuous_features.tsv', sep='\t')
marker_done = (out / '02_external_markers/ranking_LOOCV_complete.json').exists()
rankings = pd.read_csv(out / '02_external_markers/marker_composite_ranking.tsv', sep='\t') if marker_done else pd.DataFrame()
loocv = pd.read_csv(out / '02_external_markers/marker_LOOCV_per_dataset.tsv', sep='\t') if marker_done else pd.DataFrame()
comparison = pd.read_csv(out / '02_external_markers/n4_disease_vs_region_comparison.tsv', sep='\t') if marker_done else pd.DataFrame()
annotation = states.merge(prior[['CellState', 'InitialProgram', 'Judgment', 'ReviewReason', 'SourceIDs', 'Top20Genes']], on='CellState', validate='one_to_one')
annotation['AnnotationStatus'] = '候选程序对应；保留冻结谱系；不是经过验证的最终亚型'
annotation['ExternalMarkerStatus'] = 'COMPLETED' if marker_done else 'RUNNING'
for row in annotation.itertuples():
    idx = annotation.CellState.eq(row.CellState)
    evidence = recovery[recovery.CellState.eq(row.CellState)]
    annotation.loc[idx, 'SupportingDatasets_Zgt1p96'] = ';'.join(evidence.loc[evidence.OfficialZ.gt(1.96), 'Dataset'])
    annotation.loc[idx, 'MissingLineageDatasets'] = ';'.join(evidence.loc[evidence.Status.eq('NOT_EVALUABLE'), 'Dataset'])
    annotation.loc[idx, 'AssignedCells_ByDataset'] = ';'.join(evidence.Dataset+':'+evidence.AssignedCells.astype(str))
    if marker_done:
        for variant in ['DiseaseGroup', 'Region3']:
            selected = rankings[rankings.HeldOut.eq('All') & rankings.Variant.eq(variant) & rankings.CellState.eq(row.CellState) & rankings.WithinStateOrder.le(10)]
            annotation.loc[idx, 'ExternalTop10_'+variant] = ';'.join(selected.Gene)
            evaluated = loocv[loocv.CellState.eq(row.CellState) & loocv.Variant.eq(variant) & loocv.Status.eq('EVALUABLE')]
            annotation.loc[idx, 'LOOCV_EvaluableFolds_'+variant] = len(evaluated)
            annotation.loc[idx, 'LOOCV_ConcordantFolds_'+variant] = evaluated.Concordant.sum()
annotation.to_csv(out / '03_state_annotation/53_states_annotation_evidence.tsv', sep='\t', index=False)
annotation[['CellState', 'CellType', 'Ecotype', 'Top20Genes'] + [x for x in annotation if x.startswith('ExternalTop10')]].to_csv(out / '03_state_annotation/53_states_marker_summary.tsv', sep='\t', index=False)

# Step2. 文献/数据库的集合命中只作为程序提示，不赋予ORA/GSEA的P值或用单一通路命名状态
gene_sets = {}
for line in (out / '00_audit/h.all.v7.4.symbols.gmt').read_text().splitlines():
    fields = line.split('\t')
    gene_sets[fields[0]] = set(fields[2:])
pathways = []
for state in states.CellState:
    associated = set(genes.loc[genes.CellState.eq(state), 'Gene'])
    for name, panel in gene_sets.items():
        pathways.append({'CellState': state, 'GeneSet': name, 'AssociatedGenes': len(associated), 'GMTGenes': len(panel), 'OverlapCount': len(associated & panel), 'HitGenes': ';'.join(sorted(associated & panel)), 'Version': 'MSigDB7.4', 'Method': 'EXPLORATORY_EXTENSION: gene-set overlap only; no enrichment P or NES'})
pd.DataFrame(pathways).to_csv(out / '03_state_annotation/53_states_pathway_summary.tsv', sep='\t', index=False)
ecotypes = []
for ecotype, members in states.groupby('Ecotype'):
    state_cards = annotation[annotation.Ecotype.eq(ecotype)]
    support = coassoc[coassoc.Ecotype.eq(ecotype)]
    clinical = discrete[discrete.Ecotype.eq(ecotype) & discrete.Q.lt(0.05)]
    hallmark = continuous[continuous.Ecotype.eq(ecotype) & continuous.Q.lt(0.05)].sort_values('Q').head(5)
    ecotypes.append({'Ecotype': ecotype, 'FrozenMembers': ';'.join(members.CellState), 'ParentLineages': ';'.join(sorted(members.CellType.unique())), 'CandidateMemberPrograms': ';'.join(state_cards.InitialProgram.astype(str).unique()), 'Coassociation_Z': ';'.join(support.Dataset+':'+support.Z.round(3).astype(str)), 'DefinedMembers': ';'.join(support.Dataset+':'+support.DefinedReferenceMembers.astype(str)+'/'+support.FrozenMembers.astype(str)), 'DiscoveryAssociations_Qlt0p05': ';'.join(clinical.Feature+':'+clinical.GroupA+'-'+clinical.GroupB+':delta='+clinical.MeanDifference.round(4).astype(str)), 'Hallmark_Qlt0p05_Top5': ';'.join(hallmark.Feature+':delta='+hallmark.MeanDifference.round(3).astype(str)), 'LRStatus': 'BLOCKED: original scaled-fold definition / state-pair missing rules not confirmed', 'InterpretationLimit': 'BayesPrism sensitivity route; discovery associations have repeated biopsies; external coassociation has incomplete states; not validated functional communication'})
ecotypes = pd.DataFrame(ecotypes)
ecotypes.to_csv(out / '06_ecotype_features/7_ecotypes_evidence_summary.tsv', sep='\t', index=False)

# Step3. 方法对照逐模块写清公式、原版和适配；不把有歧义的全局/LR计算改成另一套算法
methods = [
    ['Model identity', 'Cell state and ecotype recovery', 'All', '53 unique states / 7 fixed ecotypes', 'No refit', 'Frozen official model + actual recovery', 'PASS', 'PAPER_ADAPTED_IBD', 'BayesPrism G-like upstream; Region3 scaling', 'Reference donor overlap not globally provable from normalized cohort IDs'],
    ['Recovery evidence', 'Significance of cell state recovery', 'Recovery', 'Official marker-consistency permutation Z', 'Existing 30 vs paper100 shuffles', 'Frozen recovery assignment and official Z', 'PARTIAL', 'PAPER_ADAPTED_IBD', 'Reuse actual official results; no refit; source different', 'Official P absent; Z evidence is not cell-lineage accuracy'],
    ['External DEG', 'State-specific marker genes in scRNA-seq data', 'Marker', 'Within-parent state vs other eligible states', 'Seurat3.1.3 min.pct0.1 logfc.threshold0.05 native ln; >=5 cells', 'Full raw transcriptome of same frozen recovered cells', 'PASS' if marker_done else 'RUNNING', 'PAPER_ADAPTED_IBD', 'Disease and Region3 stratified repeats for n4; lnFC separately converted to log2FC', 'Native Bonferroni p_val_adj kept; not relabeled as BH'],
    ['Seven metrics', 'State-specific marker genes in scRNA-seq data', 'Marker', 'n1-n7; signed two-sided Stouffer; BH -> signed MetaZ; equal mean ranks', 'Mean ties; per-state candidate BH family; P floor1e-300', 'DEG + discovery consistency + stratified coverage', 'PASS' if marker_done else 'RUNNING', 'PAPER_ADAPTED_IBD', 'n4 DiseaseGroup and Region3 two versions; same six other metrics', 'Original rank/tie/BH-family code unavailable; two datasets restrict variability'],
    ['Marker LOOCV', 'Validation of scRNA-seq markers', 'Marker', 'Retrain excluding held dataset; top10; unit-SD log2(CPM+1); maximum state mean', '10 complete markers; same parent; no held-out DEG', 'Two full-count external recovery datasets', 'PASS' if marker_done else 'RUNNING', 'PAPER_ADAPTED_IBD', 'Only two folds; each training fold one cohort; undefined folds retained', 'State assignment itself uses model markers; not independent identity truth'],
    ['State phenotype', 'Enrichment of cell states in known phenotypes', 'Phenotype', '(Ns-null mean)/null sample SD', '1000 shuffles; positive Z>1.96; seed20261012', 'Actual CD/UC/Control and inflammatory labels', 'PASS', 'PAPER_ADAPTED_IBD', 'Within-parent cells; biological-sample dominant separate with tie exclusion', 'Cell dependence; sample-dominant unit is extra adaptation; HC surgical controls'],
    ['Region representation', 'Cell state representation by tumor type', 'S2G/H', '2x2 one-sided Fisher enrichment; BH159 hypotheses', 'Dominant discovery state; Region3', 'Actual assigned biopsies, no artificial absent-state assignment', 'PASS', 'PAPER_ADAPTED_IBD', 'Tumor type -> Region3; explicit greater tail / all-state BH family', 'Repeated patients; undefined original exact family'],
    ['Co-association', 'CE recovery in scRNA-seq data', 'S5', 'F1 assigned / F2 all; C1-C4 mean; within-outside AvgDif', 'Pairwise Spearman; all4 finite;1000 row shuffles; original-diagonal swap; seed20261010', 'Existing results independent Python denominator/correlation/null-Z check', 'PASS', 'PAPER_ADAPTED_IBD', 'Missing-edge means and strict four-matrix requirement documented', 'Original P tail not explicit; retain Z, no significance stars'],
    ['Global detection', 'Probability of CE detection', 'S5I', 'Original joint CE statistic not fully decoded', 'Paper1000000 distinct from individual1000', '7 frozen ecotypes', 'BLOCKED', 'BLOCKED', 'No fabricated replacement global test', 'Original math missing in accessible BioC; joint criterion / membership null not confirmed'],
    ['Discrete ecotype features', 'Feature analysis of carcinoma ecotypes', '5A', 'Two-sided rank-sum on continuous CE abundance; mean group difference', 'BH14 disease/inflammation tests', '2025 CD/UC MSCCR biopsies', 'PASS', 'PAPER_ADAPTED_IBD', 'Disease / inflammation actual labels; Region3 separate7 KW tests', 'Discovery-associated; repeated biopsies not independent patient outcomes'],
    ['Hallmark ecotype features', 'Feature analysis of carcinoma ecotypes', '5A', 'Mean per-gene unit-SD log2 bulk expression; dominant CE Wilcoxon', '50 Hallmark sets; BH350; ddof1; log2(CPM+1)', 'MSCCR harmonized bulk counts + frozen dominant CE', 'PASS', 'PAPER_ADAPTED_IBD', 'MSigDB7.4 explicitly archived; harmonized full-library CPM', 'Paper MSigDB release not confirmed; not original CIBERSORTx upstream'],
    ['State biological cards', 'State characterization / marker correspondence', 'Interpretation', 'Separate discovery genes / literature programs / external markers / recovery', '53 fixed states; no forced renaming', 'Prior literature panel and actual tables', 'PARTIAL', 'PAPER_ADAPTED_IBD', 'IBD literature and source-attributed provisional programs', 'Mast/Glial cross-lineage concerns preserved; overlap is not formal GSEA'],
    ['Ligand-receptor', 'Ligand-receptor enrichment analysis', 'S5G', 'Welch DE; q0.05 log2FC0.1;100 gene draws; Stouffer state pairs', 'Need original fold scale / weights / undefined-pair rule / source DB', 'BayesPrism G-like not digitally purified HiRes', 'BLOCKED', 'BLOCKED', 'No z-score fold ratios or substitute CellChat', 'Cannot yet claim original-equivalent LR; FigureG not fabricated'],
    ['Existing bulk/nonIBD', 'Cell state and ecotype recovery', 'Sensitivity', 'Existing results not newly fitted', 'Keep within-MSCCR nonIBD separate from external', 'GSE72819 UC/EUCALYPTUS; prior control recoveries', 'PARTIAL', 'PAPER_ADAPTED_IBD', 'Source independence distinguished; healthy subsets of main scRNA cohorts are not new datasets', 'Full new metadata join / all-run manifest not completed this batch']
]
parity = pd.DataFrame(methods, columns=['Module', 'PaperSection', 'Figure', 'Formula', 'Parameters', 'ActualInput', 'Validation', 'ParityClass', 'Adaptation', 'Unresolved'])
parity['OriginalSoftwareVersion'] = 'Seurat3.1.3 / fgsea1.14.0 where specified; full author code not provided'
parity['ActualSoftwareVersion'] = 'R4.0.3 Seurat3.1.3 fgsea1.14.0; scipy versions recorded separately'
parity.to_csv(out / '00_audit/method_parity_matrix.tsv', sep='\t', index=False)
pd.DataFrame([{'Status': 'BLOCKED', 'Reason': 'Original global joint criterion and random membership constraints not confirmed; do not turn individual1000 into global1000000'}]).to_csv(out / '05_ecotype_coassociation/ecotype_global_detection.tsv', sep='\t', index=False)
pd.DataFrame([{'Status': 'BLOCKED', 'Reason': 'Original LR fold scale and pair meta rules incomplete; no functional communication inferred'}]).to_csv(out / '07_ligand_receptor/LR_module_status.tsv', sep='\t', index=False)

# Step4. 可阅读报告只引用实际结果，完整逐状态/逐生态型表作为审阅入口
marker_text = '两版七指标和两个留出折已计算，具体差异见 n4_disease_vs_region_comparison.tsv 与 marker_LOOCV_summary.tsv。' if marker_done else '两版七指标和LOOCV正在tmux计算，尚不能报告两版实际排名差异。'
method_table = '| Module | Validation | Parity class | Unresolved |\n|---|---|---|---|\n'
for row in parity.itertuples():
    method_table += '| '+ ' | '.join([row.Module, row.Validation, row.ParityClass, row.Unresolved]) + ' |\n'
audit = '# 方法对照审计\n\n本轮针对冻结53状态/7生态型，只读原件，不重新discovery或修改recovery assignment。所有涉及上游表达的结果均为BayesPrism敏感性分析，不宣称CIBERSORTx数值等价。\n\n'+marker_text+'\n\n'+method_table+'\n\n主要依据：[Luca et al., Cell 2021 STAR Methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC8526411/)。全局CE与LR阻断只影响各自模块，没有用替代检验冒充原文。服务器环境规格、函数默认源码、数值烟雾测试及输入校验保存在00_audit。\n'
cell_report = '# 53个cell states的初步证据\n\n恢复分配、官方Z、外部marker复现、作者谱系和文献程序是不同证据。不能仅凭高Z确认细胞身份。\n\n'+marker_text+'\n\n'
for row in annotation.itertuples():
    cell_report += f'## {row.CellState}（{row.Ecotype}）\n\n父类型：{row.CellType}。候选程序：{row.Judgment}。\n\nDiscovery top20：{row.Top20Genes}。\n\n官方正向Z支持：{row.SupportingDatasets_Zgt1p96 or "无"}；缺失谱系：{row.MissingLineageDatasets or "无"}。分配细胞：{row.AssignedCells_ByDataset}。\n\n限制：{row.ReviewReason}。文献/数据库来源编号：{row.SourceIDs}，对应既有文献来源表；这不是本轮重新验证后的细亚型定名。\n\n'
eco_report = '# 7个ecotypes的证据\n\n在两轮单细胞中，成员恢复与共关联不是同一件事。E1和E3的共关联Z在两个队列均较高；其余生态型须结合可定义成员、NA边和队列覆盖逐一看。原文共关联P的方向尚未确认，本轮不补写P或显著性星号。\n\nMSCCR疾病/炎症和50个Hallmark是discovery关联，单位为活检，存在患者重复与上游共同表达程序，不能称为独立临床验证。\n\n'
for row in ecotypes.itertuples():
    eco_report += f'## {row.Ecotype}\n\n冻结组成：{row.FrozenMembers}。\n\n成员候选程序：{row.CandidateMemberPrograms}。\n\n共关联Z：{row.Coassociation_Z}；可定义成员：{row.DefinedMembers}。\n\n疾病/炎症关联：{row.DiscoveryAssociations_Qlt0p05 or "无BH q<0.05"}。\n\nHallmark（按q列前5，完整50集另存）：{row.Hallmark_Qlt0p05_Top5}。\n\nLR：{row.LRStatus}。限制：{row.InterpretationLimit}。\n\n'
(out / '09_reports/METHODS_PARITY_AUDIT.md').write_text(audit, encoding='utf8')
(out / '09_reports/CELL_STATE_CHARACTERIZATION.md').write_text(cell_report, encoding='utf8')
(out / '09_reports/ECOTYPE_CHARACTERIZATION.md').write_text(eco_report, encoding='utf8')
scripts_manifest = [{'File': str(p.relative_to(out)), 'SHA256': hashlib.sha256(p.read_bytes()).hexdigest()} for p in sorted((out / 'scripts').glob('*')) if p.is_file()]
pd.DataFrame(scripts_manifest).to_csv(out / '00_audit/scripts_manifest.tsv', sep='\t', index=False)
(out / '09_reports/execution_status.json').write_text(json.dumps({'FrozenStates': len(states), 'FrozenEcotypes': states.Ecotype.nunique(), 'ExternalMarkersComplete': marker_done, 'StateRecoveryIndependentQC': 'PASS', 'EcotypeCoassociationIndependentQC': 'PASS', 'Phenotypes': 'PASS; adapted', 'HallmarkFeatures': 'PASS; adapted', 'GlobalDetection': 'BLOCKED', 'LigandReceptor': 'BLOCKED', 'AllRequestedModulesComplete': False, 'SourceDataModified': False}, indent=2))
print('REPORTS_WRITTEN; marker complete:', marker_done)
