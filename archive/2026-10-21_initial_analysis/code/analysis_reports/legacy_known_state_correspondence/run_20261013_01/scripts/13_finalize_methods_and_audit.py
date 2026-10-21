# 保存方法、全53筛选理由、验证缺口和实际进展，核对旧文件与benchmark脚本保护。
# 只写本次目录，不修改旧annotation、ranking、n4、recovery或discovery。
from pathlib import Path
from datetime import datetime, timezone, timedelta
import hashlib
import json
import subprocess
import pandas as pd

# Step1. 全53状态的五维筛选公开，功能置信度不参与新颖性评分
out = Path(__file__).resolve().parents[1]
root = Path('D:/IBD_ecotype/IBD_ecotype_analysis')
known = pd.read_csv(out / 'KNOWN_STATE_CORRESPONDENCE.tsv', sep='\t').fillna('')
summary = pd.read_csv(out / 'audit/ALL53_WITHIN_LINEAGE_SUMMARY.tsv', sep='\t').fillna('')
master = pd.read_csv(root / 'analysis_reports/legacy_bayesprism_annotation/run_20261013_01/CELL_STATE_ANNOTATION_MASTER.tsv', sep='\t').fillna('')
candidates = pd.read_csv(out / 'UNDERCHARACTERIZED_STATE_CANDIDATES.tsv', sep='\t')
screen = known.merge(summary.drop(columns=['CellType', 'Ecotype', 'SharedPositiveGeneN']), on='CellState', validate='one_to_one')
screen['ProgramSpecificityReview'] = screen.apply(lambda r: f'Top genes={r.Top50N}; shared positive={r.SharedPositiveGeneN}; discovery median effect gap={r.MedianDiscoveryEffectGap:.4g}; not independent validation', axis=1)
screen['LiteratureCorrespondenceReview'] = screen.BiologicalJudgment
screen['ExternalReproducibilityReview'] = screen.BothDatasetStableGeneN.map(lambda x: f'{int(x)} Top50 genes pass descriptive direction rule in both marker-selection cohorts; no independent validation')
screen['LineageCoherenceReview'] = screen.CellState.map(master.set_index('CellState').Uncertainty)
screen['BiologicalInterpretabilityReview'] = screen.CellState.map(master.set_index('CellState').BiologicalRationale)
screen['SelectedForProgramFollowup'] = screen.CellState.isin(candidates.loc[candidates.SelectedForProgramFollowup.eq(True), 'CellState'])
screen['SelectionOrExclusionReason'] = screen.CellState.map(candidates.set_index('CellState').Priority).fillna('Not promoted: weak specific/reproducible evidence, generic known axis, or identity/technical ambiguity; no novelty claim')
screen.loc[screen.CellState.eq('B_S03'), 'SelectionOrExclusionReason'] = 'Six stable genes but generic nuclear/regulatory program; biological novelty not established'
screen.to_csv(out / 'audit/ALL53_FIVE_DIMENSION_CANDIDATE_SCREEN.tsv', sep='\t', index=False)

# Step2. 核验来源哈希与Git；包括本次初检时已有5个修改文件
checks = []
baseline = pd.read_csv(out / 'audit/INPUT_IMMUTABILITY_BASELINE.tsv', sep='\t')
for row in baseline.itertuples():
    path = Path(row.Path)
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    checks.append({'Path': str(path), 'SHA256Before': row.SHA256Before, 'SHA256After': digest, 'Unchanged': digest == row.SHA256Before, 'BaselineScope': 'before substantive analysis'})
for manifest_name in ['GENE_INFO_SOURCE_HASHES.tsv', 'MSCCR_PATIENT_REVIEW_INPUT_HASHES.tsv']:
    manifest = pd.read_csv(out / 'audit' / manifest_name, sep='\t')
    for row in manifest.itertuples():
        path = Path(row.Path)
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        checks.append({'Path': str(path), 'SHA256Before': row.SHA256, 'SHA256After': digest, 'Unchanged': digest == row.SHA256, 'BaselineScope': 'immediately before reading this input'})
checks = pd.DataFrame(checks)
checks.to_csv(out / 'audit/INPUT_IMMUTABILITY_FINAL.tsv', sep='\t', index=False)
assert checks.Unchanged.all()
current = subprocess.run(['git', '-C', str(root), '-c', 'core.quotepath=false', 'status', '--porcelain=v1'], capture_output=True, text=True, encoding='utf-8', check=True).stdout
(out / 'audit/git_status_after.txt').write_text(current, encoding='utf-8')
before = set((out / 'audit/git_status_before.txt').read_text(encoding='utf-8-sig').splitlines())
after = set(current.splitlines())
git_check = {'BeforeEntriesRemoved': sorted(before - after), 'NewEntries': sorted(after - before), 'ExistingModifiedFileContentHashesVerified': 5, 'TaskWriteScope': str(out), 'RemoteConnectionsOpened': 0, 'NMF_BayesPrism_Recovery_Reruns': 0}
assert before == after, 'Git status changed outside already-untracked analysis_reports; inspect before claiming complete protection'
(out / 'audit/GIT_STATUS_AND_ISOLATION_CHECK.json').write_text(json.dumps(git_check, ensure_ascii=False, indent=2), encoding='utf-8')

# Step3. 核对主结果规模、统计、缺测和不可评价性
within = pd.read_csv(out / 'WITHIN_LINEAGE_STATE_COMPARISON.tsv', sep='\t')
matches = pd.read_csv(out / 'audit/ALL_PUBLISHED_SIGNATURE_MATCH_STATISTICS.tsv', sep='\t')
signatures = pd.read_csv(out / 'literature_resources/PUBLISHED_SIGNATURE_GENES.tsv', sep='\t', low_memory=False)
definitions = pd.read_csv(out / 'audit/PUBLISHED_SIGNATURE_DEFINITIONS.tsv', sep='\t')
clinical = pd.read_csv(out / 'audit/CANDIDATE_MSCCR_PATIENT_PAIRED_CONTEXT.tsv', sep='\t')
qc = pd.read_csv(out / 'audit/SIGNATURE_MATCH_INDEPENDENT_R_QC.tsv', sep='\t')
assert len(known) == 53 and known.CellState.is_unique and known.FrozenEcotype.nunique() == 7
assert set(within.CellState) == set(known.CellState) and not within.duplicated(['CellState', 'Gene']).any()
assert (within.MaxFC-within.DiscoveryScaledMeanDifference).abs().max() == 0
assert set(known.Category) == {'Established correspondence', 'Related but distinct', 'Potentially undercharacterized', 'Unresolved'}
assert not known.NovelCellSubtypeClaim.any() and not candidates.NovelSubtypeClaim.any()
assert known[known.CellState.isin(['Enteroendocrine_S01', 'Mast_S01', 'Mast_S07', 'Glial_S02'])].Category.eq('Unresolved').all()
assert qc.PError.max() < 1e-10 and qc.QError.max() < 1e-10
required = ['KNOWN_STATE_CORRESPONDENCE.tsv', 'WITHIN_LINEAGE_STATE_COMPARISON.tsv', 'UNDERCHARACTERIZED_STATE_CANDIDATES.tsv', 'CANDIDATE_STATE_EVIDENCE_REVIEW.md', 'CELL_STATE_DISCOVERY_INTERPRETATION.md']
assert all((out / name).stat().st_size > 0 for name in required)
totals = {'RequiredFiles': len(required), 'States': len(known), 'FrozenEcotypes': known.FrozenEcotype.nunique(), 'CategoryCounts': known.Category.value_counts().to_dict(), 'WithinLineageGeneRows': len(within), 'PublishedSignatureGeneRows': len(signatures), 'PublishedSignatureDefinitions': len(definitions), 'SignatureComparisons': len(matches), 'EvaluableSignatureTests': int(matches.Evaluable.sum()), 'IndependentRMaxPError': qc.PError.max(), 'IndependentRMaxQError': qc.QError.max(), 'CandidateCards': len(candidates), 'SelectedProgramFollowupStates': candidates.loc[candidates.SelectedForProgramFollowup.eq(True), 'CellState'].tolist(), 'PatientContextRows': len(clinical), 'EvaluablePatientPairedTests': int(clinical.P.notna().sum()), 'SourceHashesVerified': len(checks), 'SameCellJointCountsMeasured': False, 'NovelSubtypeClaim': False, 'NewLargeCohorts': 0, 'BenchmarkModifiedByTask': False}
(out / 'audit/FINAL_QC_SUMMARY.json').write_text(json.dumps(totals, ensure_ascii=False, indent=2), encoding='utf-8')

# Step4. 方法和资源限制明确，未检索到不等于文献中不存在
methods = f'''# Methods and Resource Audit

日期2026-10-13；仅冻结旧版53状态/7生态型。不作为CIBERSORTx未来discovery基准。所有新文件位于本版本化目录；旧版annotation、recovery、marker ranking、n4和discovery未改。

## 隔离和资源

已只读检查项目Git/本机相关进程和内存；初检本机可用约4.9 GiB，没有接入服务器确认实时余量。本轮没有SSH登录、附接会话、进程控制、Conda修改、安装、commit/push或任何Git恢复/清理操作。小表分析顺序执行，分析代码限制OMP/BLAS为1线程；观察到Python工作集约161 MiB。没有读取原始大表达矩阵、没有新增scRNA队列或重跑模型/recovery/DE。阶段末{len(checks)}个来源哈希一致，包含最初已修改的.gitignore与4个R脚本、旧注释成果、2个benchmark脚本、冻结gene_info及MSCCR审计输入；Git status条目保持一致。benchmark输入和活跃输出没有被本任务写入，也没有对活跃输出作不合理冻结哈希。不能由这些检查宣称远端benchmark已经完成。

## 文献资源盘点和真实来源

按真实DOI从EuropePMC核实8篇论文身份；新取得Friedrich/Elmentaite/Garrido-Trigo/Kinchen公开XML；Martin/Smillie全文BioC在既有sources中复用。HTTP访问和失败均保留日志。最初下载脚本把static-content域名映射错为springer-cdn.com，后修正为springer.com，失败与成功日志均保留，不将早期错误当作论文资源缺失。

实际导入Friedrich Table6（全部已导出显著上调cluster marker，保留原统计）与Table9（IL1β刺激3h完整DE，统计选集为正向adjustedP<0.05）、Elmentaite Table7（全部原表导出每cluster最多50 marker）。Elmentaite是发育/健康多肠段图谱，不能称纯IBD状态库。现有R/readxl只读导入.xls，无软件安装；原始.xls保存用于核对。一个原始/导出特殊INS标签编码异常保留未映射，不据此作负匹配。

Garrido Source Data ZIP约281 MB，不整包下载。先以HTTP206请求尾目录，再取两个<10 MB的指定成员（补充图4、5），核对长度和central-directory CRC。第一次CRC实现误用了可能为0的local-header CRC，修正后通过；不涉及服务器登录。实际只用图5d的M1/M2/IDA各200 marker，以及图5e全11898行成对DE的明确方向；图4细胞坐标/作者映射不当作新表达验证队列。成对负方向保留Direction和PublishedSignedEffect，不能把DOWN-in-IDA解释为“上调IDA身份”。Garrido signature和Colon recovery都是GSE214695，不能重复计为独立证据。

Martin完整原文枚举的5程序列表进入定量比较；它们是完整程序列表而非完整cluster DE。手工短marker panel只作诊断，不进入主统计。Smillie完整TableS2、Martin完整cluster DE、Corridoni/Parikh及Kinchen完整marker工作簿未取得，文献覆盖为PARTIAL；未匹配绝不能等同新亚型。真实URL、DOI、文件SHA256、下载/范围传输、原文段落、sheet结构均保存在audit和literature_resources。

## 文献对应统计

共{len(definitions)}个可追溯signature定义、{len(signatures)}行signature基因记录；含未映射和不用于某父类型的定义，不把它们默认为负匹配。对各冻结状态使用Discovery Top50、外部All-fold DiseaseGroup Top50及Region3 Top50。只比较同父类型/明确共享程序；MM另比较DC signature，CD4/CD8另比较MAIT/γδ候选signature作替代谱系审计。父类型背景为上一阶段哈希验证的NMF initial_gene_info所有初始状态特征；标明背景内query和signature覆盖。

共有{len(matches)}个比较，其中{int(matches.Evaluable.sum())}个满足query和背景内signature均≥5基因；对这些全部可评价比较共同BH控制，三种query及替代谱系比较均在同一检验族。统计为单侧超几何上尾；保存P、BH q、Jaccard、query coverage、signature coverage、命中基因、方向和来源。该条件特征空间中的统计不是全转录组无偏匹配、身份分类准确率、显著性差异的检验或新颖性检验。精确基因符号匹配，不做未核实的别名转换，覆盖不足标NA。R phyper/p.adjust独立验证全部最终统计，最大P误差{qc.PError.max():.3g}、Q误差{qc.QError.max():.3g}。

外部Top50本身由相同两队列筛选，发表reference独立性不消除query选择循环；原discovery区分能力不能作为独立验证。部分统计显著的是普遍谱系/周期程序，不据此确认精细亚型。对应类别由逐项人工多证据审阅决定，原功能confidence保留但不参与新颖性评分。Related but distinct表示有候选组合区别，不代表统计证实对已发表状态存在新的差异；Potentially undercharacterized仅是尚需排除已知解释的程序组件。29个Unresolved为本阶段对应/区分结论，不覆盖上阶段冻结功能候选注释。

## 同谱系比较与外部证据

解析gene_info的所有IS效应列，保留正确R行名偏移与S-to-IS映射；每个Top50实际MaxFC与解析后的目标效应完全相等。多数gene_info只包含QC后胜出基因，但对这些基因保存所有冻结状态效应；不将唯一assignment误认为独占表达。描述性正效应阈值0.1，保存共有正状态数及目标相对最佳其他状态差。两状态之间的外部表达差是供者等权均值的描述差，不是新供者配对DE；0.5用于呈现，不是显著性阈值。

重用已完成的全状态表达、作者全谱系表达、状态间DE、全部8812关联的供者配对方向和作者亚型构成。目标/其余已分配状态各≥5细胞的既有供者对比保留，Unassigned不用于对照；筛查规则为每队列≥3匹配供者、平均差>0、≥70%供者方向正。报告两队列同时满足的基因，不把不足供者当失败。该筛查是描述性而非P值；粗粒度作者bin与抽样上限限制解释。CLEC10A因明确排除DC-like目的检查全关联列表，原排名102，未重新marker排名或改n4。

当前逐细胞audit_counts和state_donor_expression检查点未在本地找到。没有用均值代替观测同细胞联合表达。对同一群体、同一供者等权定义的边际检出p_i仅计算Fréchet界限：lower=max(0,sum(p_i)-(m-1))、upper=min(p_i)。界限依赖同细胞群体/等权口径一致，不能替代观测联合计数、供者内共表达或多基因得分。缺测/无输入谱系/无分配细胞都与实际0检出分别标记。

## 患者重复与混杂

采用Profile审计的旧MSCCR2025 H和metadata（2025活检、919患者），不是GSE72819的73例recovery快照，也不是当前MSCCR120 fraction benchmark。7个预定候选/比较对象按患者×疾病×Region3×炎症聚合活检中位数，随后同患者同疾病同肠段Inflamed−Non-inflamed配对。仅≥10患者对推断，CD/Colon=41、UC/Colon=46；Ileum/Rectum不足对数保留描述/NA。

两侧Wilcoxon signed-rank，正常近似与连续性校正；主报告将绝对差<1e-10视为数值零，并保留完全原值的敏感性版本，不回写原H。BH只作用于{int(clinical.P.notna().sum())}个可评价比较。MM_S03 UC/Colon显著性对数值零规则敏感，明确标为不稳健。疾病/肠段固定与患者配对控制主要已知混杂，年龄/治疗/时间未建模。旧H是初始NMF相对分量，QC过滤后保留状态之和可能<1，未重新归一；关联不能解释为绝对细胞数、功能因果或独立临床验证。旧Region3 Fisher仅作来源背景，未重复当患者水平验证。

## 本轮结果与下一步

候选优先级及正/反向证据见UNDERCHARACTERIZED_STATE_CANDIDATES；全部53的五维screen单独保存，B_S03虽有6稳定基因但泛核内程序不晋升。CD8_S03的IL23R、MM_S02的PLA2G2D作为稳定核心的解释不成立；Goblet_S01含周期阳性成分，不能把整个状态命名为新增殖群。MM_S02/MMS03的M2/IDA已知解释得到真实完整发表数据支持，必须先排除。CD8_S03与CD4_S04仅为暂定未充分解释的多基因组件，没有Novel亚型结论。

后续任务仅提出有明确目的的验证，不执行：取本来存在的目标谱系计数检查点以测量核心联合表达；以完整Smillie/Martin及已知DC/巨噬/Trm等signature作匹配对照；匹配作者亚型/TCR/深度/周期后评价残余区分；用独立供者和原位联合锚点验证。这些缺口不通过重新NMF、BayesPrism或高成本recovery弥补。任何远端取数/计算需另行明确授权并遵循单persistent SSH规则。
'''
(out / 'METHODS_RESOURCE_AND_SAFETY_AUDIT.md').write_text(methods, encoding='utf-8')
status_rows = [('Isolation and frozen inputs', 'PASS', f'{len(checks)} source hashes unchanged; original modified files preserved; no SSH/model rerun/environment change'), ('Resource inventory', 'PASS', 'Local evidence and actual published supplement resources inventoried; failures and scope preserved'), ('All53 correspondence and within-parent screen', 'PASS', '53 categories, all-state gene rows, complete available signature statistics and all53 five-dimensional screen'), ('Published signature coverage', 'PARTIAL', 'Full available Friedrich/Elmentaite/Garrido exports and Martin enumerated programs used; other full cluster tables missing'), ('Candidates and falsifiable hypotheses', 'PASS', '11 cards, 5 program followup objects, 4 identity-first audits; no new subtype claim'), ('Same-cell joint expression', 'BLOCKED_LOCAL_INPUT', 'Count checkpoints absent locally; only transparent Frechet bounds, not observed joint counts'), ('Patient context', 'PARTIAL', '14 paired inference tests; 28 sparse strata NA; extra clinical covariates not modeled'), ('Independent validation and novelty', 'NOT_ESTABLISHED', 'Same marker-selection cohorts reused; known M2/IDA and crypt-cycling explanations retained')]
pd.DataFrame(status_rows, columns=['Module', 'Status', 'EvidenceOrGap']).to_csv(out / 'PASS_PARTIAL_BLOCKED.tsv', sep='\t', index=False)
future_rows = [
    ('Same-cell joint program', 'CD8_T_S03;CD4_T_S04;MM_S02;MM_S03;Goblet_S01', 'Existing target-lineage audit_counts/count_order checkpoints', 'Current pooled marginals cannot resolve observed joint detection or donor-specific scores', 'Target genes only; one process; no new recovery; remote access not executed'),
    ('Known-state exhaustive comparison', 'All; priority myeloid/T/Goblet', 'Smillie TableS2, Martin complete cluster DE, Corridoni/Parikh/Kinchen published signature tables', 'Current complete exports do not exhaust known-cell programs', 'Small supplementary tables only; no new cohort counts required'),
    ('Alternative lineage and composition', 'MM_S02;CD8_T_S03;identity audits', 'Joint macrophage/DC markers; TCR/author subtype; negative-lineage anchors', 'CSF1R/CLEC10A or CD8/TRGC2 label alone cannot exclude DC/gamma-delta/mixture', 'Subset existing cells, within-donor and same-author-subtype comparisons'),
    ('Independent replication', 'Two tentative undercharacterized components', 'Donors/cells not used in current discovery/marker-selection chain', 'Two current external cohorts do not provide a fresh independent layer', 'New large cohort not downloaded; require purposeful independent donor design')]
pd.DataFrame(future_rows, columns=['Purpose', 'States', 'NeededResource', 'Reason', 'ScopeAndCurrentStatus']).to_csv(out / 'audit/DEFERRED_TARGETED_VALIDATION.tsv', sep='\t', index=False)

# Step5. 实际进展和交付入口
zone = timezone(timedelta(hours=8))
progress = ['# Actual progress log', '', '本次全部在独立目录执行。所有失败日志保留，未恢复/重启任何benchmark任务。', '']
for name in ['01_inventory_resources', '02_get_published_tables_corrected', '03_public_archive_tables', '04_probe_existing_evidence', '05_import_published_xls', '06_compare_states_and_signatures', '07_prepare_patient_context_review', '08_numeric_and_patient_checks_final', '09_inventory_large_supplement_by_range_corrected', '10_get_targeted_source_tables_final', '11_add_verified_myeloid_signatures_final', '12_classify_and_write_scientific_reports_final']:
    path = out / 'logs' / (name + '.log')
    progress.append(f'- {datetime.fromtimestamp(path.stat().st_mtime, zone).isoformat(timespec="seconds")} | {name} | 实际日志：`logs/{path.name}`。')
progress += ['', '最终五项核心文件完整；来源/数值/身份分类检查通过。未补做逐细胞联合检测，也没有独立验证/新亚型结论。']
(out / 'PROGRESS_LOG.md').write_text('\n'.join(progress), encoding='utf-8')
readme = f'''# 交付入口：已知状态对应与待深入程序

冻结旧版53状态/7生态型，2026-10-13；benchmark、旧annotation/ranking/n4等未修改。

- [KNOWN_STATE_CORRESPONDENCE.tsv](KNOWN_STATE_CORRESPONDENCE.tsv)：53状态类别、真实文献对应、统计摘要与独立性限制。
- [WITHIN_LINEAGE_STATE_COMPARISON.tsv](WITHIN_LINEAGE_STATE_COMPARISON.tsv)：{len(within)}行同父类型基因效应、共享性、表达区分、供者方向和已有DE。
- [UNDERCHARACTERIZED_STATE_CANDIDATES.tsv](UNDERCHARACTERIZED_STATE_CANDIDATES.tsv)：11项程序候选/比较对象/身份审计，五维评价与可证伪假设。
- [CANDIDATE_STATE_EVIDENCE_REVIEW.md](CANDIDATE_STATE_EVIDENCE_REVIEW.md)：重点状态完整证据卡，包括反向证据、作者亚型、联合检出界限和患者配对结果。
- [CELL_STATE_DISCOVERY_INTERPRETATION.md](CELL_STATE_DISCOVERY_INTERPRETATION.md)：按Cell2021 Results逻辑组织科学解释。

分类：8 Established、14 Related but distinct、2 Potentially undercharacterized、29 Unresolved。类别针对本轮程序对应/区分，不改变原功能confidence或任何State ID。两个待深入组件为CD8_T_S03（JUND/LAPTM5/NR1H2相关，IL23R非稳定核心）与CD4_T_S04（ZNF831/CD6/FERMT3相关）；不是新亚型认定。

MM_S02更接近已知M2相关程序，P2RY13/CLEC10A值得细化、PLA2G2D跨队列不复现；MM_S03有已知IDA相关对应；Goblet_S01值得量化已知周期阳性成分而非整体定名增殖。四个指定谱系疑点保持Unresolved。

完整已取得signature导出和统计见literature_resources与audit；最终{len(matches)}匹配比较、{int(matches.Evaluable.sum())}个可评价检验有独立R数值核对。不得把Garrido与GSE214695对应重复计为独立验证。

本轮没有足够证据宣布新细胞亚型。逐细胞计数检查点本地缺失、部分完整发表表未取得、稀少患者配对与额外混杂仍限制结论。见[方法与资源审计](METHODS_RESOURCE_AND_SAFETY_AUDIT.md)、[模块状态](PASS_PARTIAL_BLOCKED.tsv)、[进展日志](PROGRESS_LOG.md)。
'''
(out / 'README.md').write_text(readme, encoding='utf-8')
print(json.dumps(totals, ensure_ascii=False, indent=2))
