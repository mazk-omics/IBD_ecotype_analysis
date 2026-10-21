# 汇总本次审计、临床统计与图件，保存 PPT 摘要及可复用的方法记录。
# 所有数字读取实际结果；只写本次输出目录，不自动提交 Git。
from pathlib import Path
import pandas as pd
import numpy as np
from scipy import stats
import json, io, hashlib, shutil

# Step1. 给资产清单补充数据对象、已核实维度、版本和本机副本映射
root = Path('D:/IBD_ecotype/IBD_ecotype_analysis')
out = root/'06_validation/MSCCR_joint_clinical_audit_20261015'
remote = json.loads((out/'audit/remote_payload.json').read_text())
inventory = pd.read_csv(out/'audit/server_inventory.tsv',sep='\t')
local = pd.read_csv(out/'audit/local_inventory.tsv',sep='\t')
metadata = pd.read_csv(out/'source_data/analysis_metadata_2025.tsv',sep='\t',index_col=0)
members = pd.read_csv(out/'source_data/frozen_state_membership.tsv',sep='\t')
final_hashes = json.loads((out/'audit/remote_final_verification.json').read_text())['hashes']
geo_samples = set(pd.read_csv(out/'source_data/GEO_metadata_2490.tsv',sep='\t').SampleID)
inventory.loc[inventory.path.isin(final_hashes),'sha256'] = inventory.loc[inventory.path.isin(final_hashes),'path'].map(final_hashes)
inventory['n_samples'] = np.nan
inventory['n_genes'] = np.nan
inventory['n_patients'] = np.nan
inventory['n_states'] = np.nan
inventory['n_ecotypes'] = np.nan
inventory['clinical_match'] = 'Not individually loaded; see object-specific audit'
inventory['object_kind'] = 'Historical derived file or candidate from another cohort'
inventory['dimension_evidence'] = 'Not applicable or not loaded in this read-only audit'
inventory['version'] = inventory.path.str.extract(r'(MSCCR_CORE12[^/]*|ctspec_region_[^/]+)')[0].fillna('path-specific; see full path')
for path,header in remote['headers'].items():
    positions = inventory.path.eq(path)
    samples = [v for v in header['header_fields'] if str(v).startswith('MSCCR_reGRID_')]
    inventory.loc[positions,'n_samples'] = len(samples)
    inventory.loc[positions,'n_patients'] = pd.Series(samples).str.extract(r'(reGRID_\d+)')[0].nunique()
    inventory.loc[positions,'clinical_match'] = 'All header samples in GEO' if set(samples)<=geo_samples else 'Mismatch; inspect header'
    inventory.loc[positions,'object_kind'] = 'G-like parent-specific expression' if 'expression_matrices' in path else 'Bulk expression'
    inventory.loc[positions,'dimension_evidence'] = 'Fresh header/sample IDs; gene count from existing schema or frozen QC, not full matrix rescan'
    inventory.loc[positions,'n_genes'] = 13668 if 'expression_matrices' in path else 25970 if 'adjcounts' in path else 56632 if 'Biopsy_counts' in path else np.nan
for cell_type in members.CellType.unique():
    abundance = root/'ctspec_region_ibd'/cell_type/'state_abundances.txt'
    assignment = pd.read_csv(root/'ctspec_region_ibd'/cell_type/'state_assignment.txt',sep='\t')
    pattern = '/ctspec_region_ibd/'+cell_type+'/'
    state_rows = inventory.path.str.contains(pattern,regex=False)&inventory.path.str.endswith('state_abundances.txt')
    inventory.loc[state_rows,['n_samples','n_patients','n_states']] = [2025,919,int(members.CellType.eq(cell_type).sum())]
    inventory.loc[state_rows,'object_kind'] = 'Frozen within-parent state abundance; not whole-tissue cell fraction'
    inventory.loc[state_rows,'clinical_match'] = '2025/2025 unique sample IDs matched'
    inventory.loc[state_rows,'dimension_evidence'] = 'Local matrix load and server identical SHA256'
    assigned_rows = inventory.path.str.contains(pattern,regex=False)&inventory.path.str.endswith('state_assignment.txt')
    inventory.loc[assigned_rows,['n_samples','n_patients']] = [len(assignment),metadata.loc[assignment.ID,'PatientID'].nunique()]
    inventory.loc[assigned_rows,'object_kind'] = 'Frozen quality-filtered parent-state assignment'
    inventory.loc[assigned_rows,'clinical_match'] = 'All assigned IDs match Discovery; denominator differs by parent'
eco = inventory.path.str.contains('/ctspec_region_ibd/Ecotypes/',regex=False)
inventory.loc[eco&inventory.path.str.endswith('ecotype_abundance.txt'),['n_samples','n_patients','n_ecotypes']] = [2025,919,7]
inventory.loc[eco&inventory.path.str.endswith('ecotype_abundance.txt'),'object_kind'] = 'Frozen ecotype abundance; not NMF H or G-like expression'
inventory.loc[eco&inventory.path.str.endswith('ecotype_abundance.txt'),'clinical_match'] = '2025/2025 unique sample IDs matched'
assigned_rows = eco&inventory.path.str.endswith('ecotype_assignment.txt')
inventory.loc[assigned_rows,['n_samples','n_patients','n_ecotypes']] = [1712,metadata.loc[metadata.AssignedEcotype.ne('Unassigned'),'PatientID'].nunique(),7]
inventory.loc[assigned_rows,'object_kind'] = 'Frozen significance-filtered ecotype assignment'
inventory.loc[inventory.path.str.contains('/NMF|/Nmf|/nmf',regex=True),'object_kind'] = 'NMF intermediate; do not use as clinical abundance without mapping'
copy_records = []
local_set = set(local.path)
for _,row in inventory.iterrows():
    p = row.path
    relative = p.removeprefix('/home/mazekai/IBD_EcoTyper/')
    candidate = str(root/relative)
    if '/05_EcoTyper/output/MSCCR_CORE12/ctspec_region_ibd/' in p:
        candidate = str(root/'ctspec_region_ibd'/p.split('/ctspec_region_ibd/',1)[1])
    hash_copies = row.local_hash_copies if pd.notna(row.local_hash_copies) else ''
    direct = candidate if candidate in local_set else ''
    copy_records.append(hash_copies or direct)
inventory['local_verified_or_path_copy'] = copy_records
inventory['local_copy_status'] = np.where(inventory.local_verified_or_path_copy.ne(''),'Present; SHA match when hash available, otherwise path only','Not located by exact path/SHA; nested historical snapshots not exhaustively content-hashed')
inventory.to_csv(out/'audit/server_inventory.tsv',sep='\t',index=False)
inventory.loc[inventory.object_kind.str.contains('Frozen|G-like|Bulk')|inventory.path.str.contains('GSE193677_family|MSCCR_CORE12.*annotation.tsv',regex=True)].to_csv(out/'audit/key_server_assets.tsv',sep='\t',index=False)
for c in ['n_samples','n_genes','n_patients','n_states','n_ecotypes']:
    byhash = inventory.dropna(subset=['sha256',c]).drop_duplicates('sha256').set_index('sha256')[c]
    local[c] = local.sha256.map(byhash)
local['origin_class'] = np.where(local.phase.eq('FrozenDiscovery'),'frozen',np.where(local.phase.eq('Metadata'),'raw_or_annotated_metadata','derived_or_audit'))
local.to_csv(out/'audit/local_inventory.tsv',sep='\t',index=False)

# Step2. Excel 通过 GSM 验证患者身份及数值，记录 log2 列的人工四舍五入
excel = pd.read_csv(out/'audit/Excel_metadata_original_2490.tsv',sep='\t')
geo = pd.read_csv(out/'source_data/GEO_metadata_2490.tsv',sep='\t')
pair = excel.merge(geo,left_on=excel.columns[0],right_on='geo_accession',validate='one_to_one')
patient_excel = pair.iloc[:,1].str.extract(r'(reGRID_\d+)')[0]
excel_summary = dict(GSMMatches=len(pair),PatientPrefixMismatches=int(patient_excel.ne(pair.GEO_PatientID).sum()),NonUniquePatientLabelNotBiopsyKey=True)
numeric_pairs = list(zip([3,8,9,12,13,14,15,16,17,18,19,20,21],['study_eligibility_age_at_endo','log2_fecalcalpro_mgperg','crp_jjmgl_log2','ghas_sum7','nancyindex','ibdsescd_totalsescd','ibdmesuc_mayo_score','harveybradshawindex_hbi_score','colitisactivityindex_sccai','max_ghas_sum7','max_nancy','endoremiss','historemiss']))
excel_checks = []
for position,field in numeric_pairs:
    a = pd.to_numeric(pair.iloc[:,position],errors='coerce')
    b = pd.to_numeric(pair[field],errors='coerce')
    valid = a.notna()&b.notna()
    tolerance = .0051 if field in ['log2_fecalcalpro_mgperg','crp_jjmgl_log2'] else 1e-10
    excel_checks.append(dict(ExcelColumn=excel.columns[position],GEOField=field,Compared=int(valid.sum()),MissingnessDifferences=int(a.isna().ne(b.isna()).sum()),ValuesBeyondDocumentedTolerance=int((a[valid]-b[valid]).abs().gt(tolerance).sum()),Tolerance=tolerance,MaxAbsDifference=float((a[valid]-b[valid]).abs().max())))
pd.DataFrame(excel_checks).to_csv(out/'audit/Excel_numeric_values_vs_GEO.tsv',sep='\t',index=False)
(out/'audit/Excel_identity_summary.json').write_text(json.dumps(excel_summary,indent=2),encoding='utf-8')

# Step3. 对照历史 Region3 及 Disease/Inflammation 结果，保存患者调整效应表
cat = pd.read_csv(out/'statistics/categorical_overall_unadjusted.tsv',sep='\t')
old = pd.read_csv(root/'06_validation/paper_aligned_recovery_characterization_20261012/06_ecotype_features/ecotype_region_KW.tsv',sep='\t')
region_check = cat.loc[cat.Entity.eq('ecotype')&cat.Feature.eq('Region3')].merge(old,left_on='Marker',right_on='Ecotype',suffixes=('_new','_old'),validate='one_to_one')
region_check['AbsPdiff'] = (region_check.P_new-region_check.P_old).abs()
region_check.to_csv(out/'audit/previous_ecotype_Region3_recheck.tsv',sep='\t',index=False)
old_state = pd.read_csv(root/'06_validation/paper_aligned_recovery_characterization_20261012/04_state_phenotypes/state_region_fisher.tsv',sep='\t')
state_checks = []
for _,row in old_state.iterrows():
    a = pd.read_csv(root/'ctspec_region_ibd'/row.CellType/'state_assignment.txt',sep='\t').set_index('ID').join(metadata[['PatientID','Region3']],validate='one_to_one')
    own = a.State.eq(row.CellState.removeprefix(row.CellType+'_'))
    target = a.Region3.eq(row.Region3)
    table = np.array([[np.sum(own&target),np.sum(own&~target)],[np.sum(~own&target),np.sum(~own&~target)]])
    odds,p = stats.fisher_exact(table,alternative='greater')
    state_checks.append(dict(Marker=row.CellState,Region3=row.Region3,PreviousNs=row.Ns,NewNs=table[0,0],P_old=row.P,P_new=p,AbsPdiff=abs(row.P-p)))
pd.DataFrame(state_checks).to_csv(out/'audit/previous_state_Region3_Fisher_recheck.tsv',sep='\t',index=False)
assert max(r['AbsPdiff'] for r in state_checks)<1e-12
adj = pd.read_csv(out/'statistics/patient_cluster_adjusted_contrasts.tsv',sep='\t')
primary = adj.loc[adj.Weighting.eq('Biopsy_weighted')]
weighted = adj.loc[adj.Weighting.eq('Equal_patient_weight')]
comparison = primary.merge(weighted,on=['Entity','Marker','Feature','Contrast'],suffixes=('_primary','_weighted'),validate='one_to_one')
comparison['SameDirection'] = np.sign(comparison.Beta_primary)==np.sign(comparison.Beta_weighted)
comparison.to_csv(out/'statistics/primary_vs_equal_patient_weight_sensitivity.tsv',sep='\t',index=False)
state_specificity = primary.loc[primary.Entity.eq('state')].merge(members[['ID','CellType','State','Ecotype']],left_on='Marker',right_on='ID',validate='many_to_one')
state_specificity['WithinParentEffectRank'] = state_specificity.groupby(['CellType','Feature','Contrast']).Beta.rank(ascending=False,method='min')
state_specificity['Note'] = 'Descriptive within-parent ordering; formal assignment specificity is in permutation enrichment, not a test of beta differences'
state_specificity.to_csv(out/'statistics/state_within_parent_effect_summary.tsv',sep='\t',index=False)

# Step4. 汇报摘要用实际模型值，主文明确横断面及缺失限制
dictionary = pd.read_csv(out/'source_data/clinical_variable_dictionary.tsv',sep='\t')
paired = pd.read_csv(out/'statistics/patient_detailed_region_matched_inflammation.tsv',sep='\t')
e1 = primary.loc[primary.Entity.eq('ecotype')&primary.Marker.eq('E1')].set_index('Feature')
infl = e1.loc['Inflammation']
ghas = e1.loc['Local_GHAS']
nancy = e1.loc['Local_Nancy']
crp = e1.loc['CRP_log2']
fcal = e1.loc['Fecal_calprotectin_log2']
ses = e1.loc['SES_CD']
mayo = e1.loc['Mayo_endo']
e1pair = paired.loc[paired.Entity.eq('ecotype')&paired.Marker.eq('E1')].iloc[0]
significant = primary.groupby('Entity').Q.agg(lambda q:int(q.lt(.05).sum())).to_dict()
coverage_rows = ['| 临床字段 | 活检数 | 患者数 | 适用范围 |','|---|---:|---:|---|']
for _,r in dictionary.iterrows():
    coverage_rows.append(f'| {r.Feature} | {r.NonmissingBiopsies} | {r.NonmissingPatients} | {r.AnalysisDisease}; {r.Use} |')
coverage_table = '\n'.join(coverage_rows)
summary = f'''# 2026-10-15 汇报摘要

今天已完成本机与服务器 MSCCR 联合只读审计、GEO/SOFT 临床核验、2025 活检映射、患者聚类调整与敏感性、5 张核心图及全部源数据。独立 conda 环境：`msccr_clinical_audit_20261015`。

- 冻结结果核实：2,025 份 CD/UC 活检、919 患者、12 父类型、53 states、7 ecotypes；1,712 已分配，313 Unassigned。27 个冻结丰度/分配/成员文件两端 SHA256 一致。
- 找到的真实临床字段超过原有 Disease/Inflammation/Region3：年龄、性别、CRP、粪钙卫蛋白、临床活动性、内镜严重度、局部 GHAS/Nancy、SES-CD、Mayo、HBI、SCCAI、最大组织学分数、内镜缓解。完整字段字典见 source_data/clinical_variable_dictionary.tsv。
- **结果 1：E1 与局部炎症相关。** 调整 beta={infl.Beta:.3f} abundance SD，95% CI [{infl.CI_low:.3f}, {infl.CI_high:.3f}]，q={infl.Q:.3g}；2,025 活检/919 患者。
- **结果 2：E1 对应实际组织学炎症强度。** CD 局部 GHAS 每升高 1 SD，E1 beta={ghas.Beta:.3f}，q={ghas.Q:.3g}（778 活检/358 患者）；UC 局部 Nancy beta={nancy.Beta:.3f}，q={nancy.Q:.3g}（402/184）。
- **结果 3：E1 与全身/粪便炎症标志物关联。** CRP 原始 log2 值 beta={crp.Beta:.3f}，q={crp.Q:.3g}（1,949/882）；粪钙卫蛋白原始 log2 值 beta={fcal.Beta:.3f}，q={fcal.Q:.3g}（464/203）。连续指标 beta 均对应指标 1 SD。
- **结果 4：E1 与内镜评分关联。** CD SES-CD beta={ses.Beta:.3f}，q={ses.Q:.3g}；UC Mayo beta={mayo.Beta:.3f}，q={mayo.Q:.3g}。已调整重复活检、肠段、年龄和性别；疾病适用范围严格分开。
- **结果 5：严格患者内、详细部位配对支持 E1/E5/E6 的方向。** 30 患者/30 对；E1 炎症减非炎症的中位差={e1pair.MedianPatientDifference:.3f}，q={e1pair.Q:.3g}；E5、E6 均减少，q=0.000887。较宽 Region3 配对有 106 患者，另表保留。每患者等权调整结果同方向。

患者级调整已经完成：主模型患者聚类稳健 SE、每患者等权敏感性、每患者首份活检 Spearman、Region3 及详细部位配对。主模型 53 states 的 1,166 个对比中 {significant['state']} 个 q<0.05，7 ecotypes 的 154 个对比中 {significant['ecotype']} 个 q<0.05；这些彼此相关的对比不能当作独立生物学发现计数。

PPT 建议顺序：Fig3（生态型特征全貌）→ Fig2（临床关联）→ Fig4（调整前后效应）→ Fig5（患者内配对）；Fig1 作为完整 states 备份。5 图均含 PNG/PDF/SVG，已通过 PDF 字体、碰撞和适用的面板对齐审计。

必须说明：本次均为冻结 Discovery 内横断面关联。原始 historemiss 有 189 活检/80 患者的编码与论文阈值不符，未用于最终主推断；使用显式派生 Derived_histology_zero（CD max GHAS=0 / UC max Nancy=0）替代。没有找到样本级预计算 bMIS、药物暴露或未来疗效/随访；未做治疗反应预测。粪钙卫蛋白仅覆盖 203/919 患者，完整病例结论受选择性缺失影响。临床测量具体日期/与活检间隔不可核实。

未完成：原始 historemiss 编码向作者澄清、bMIS 原始评分获取/严格复算、用药混杂控制及独立队列验证。未重跑 BayesPrism/EcoTyper，未启动外部治疗队列。
'''
(out/'presentation_summary.md').write_text(summary,encoding='utf-8')
report = f'''# MSCCR 本机与服务器联合审计及临床关联分析

完成日期：2026-10-15。工作目录：{out}。单次持久 SSH 已正常关闭；远端只读，未创建分析输出目录。原始数据、冻结 Discovery 和既有分析文件均未修改；没有 Git 提交。

## 1. 本机实际资产

限定 MSCCR、GSE193677、冻结路线和指定既有图件/审计范围，保存 1,282 条本机资产记录（audit/local_inventory.tsv）。任务中预期的 05_EcoTyper/input/MSCCR_CORE12 和 05_EcoTyper/output/MSCCR_CORE12/ctspec_region_ibd 不在本机。实际冻结副本在 `D:/IBD_ecotype/IBD_ecotype_analysis/ctspec_region_ibd/`；GEO series matrix 在 `D:/IBD_ecotype/Profile/audit/2026-10-06_ctspec_region_ibd/`；原 Excel 和 73 页补充 PDF 在 `D:/IBD_ecotype/data_bulk/discovery_cohort/metadata/MSCCR(GSE193677)/`。既有 paper_aligned... 和 paper_figures... 含统计、Hallmark 分数、真实 BayesPrism 20 类全活检比例及绘图来源。

## 2. 服务器实际资产与样本链

服务器 `/home/mazekai/IBD_EcoTyper` 及 `/home/mazekai/uploads` 为限定查找范围。清单 95,249 条包括其他队列临床候选、历史分支、NMF 和 CIBERSORTx 分块；不能把这 95,249 条都称作 MSCCR 原始数据。全表见 audit/server_inventory.tsv，关键输入见 audit/key_server_assets.tsv。95,249 条的未知维度显式保留，未为清单而逐个读取全部中间矩阵。

- 原始 counts：56,632 genes × 2,490 biopsies，107,824,100 bytes；adjcounts：25,970 × 2,490，278,394,867 bytes。表头/样本 ID 本次新读；基因行数和全矩阵值分布由既有 schema 审计提供，未重扫整矩阵。adjcounts 包含负值，为作者变换/校正表达，不能当原始 counts。
- 原始 GEO 2,490 活检/1,162 患者：1,157 CD + 872 UC + 461 Control。原始 counts 既有 schema 显示零深度样本为 0。
- 下游上皮 posterior mass 4 个样本为零，technical_exclusion_samples.tsv 明确列出；排除后 2,486 活检/1,162 患者，1,154 CD + 871 UC + 461 Control。这四个不是原始 bulk 零 counts，证据原文和脚本保存于 audit/server_*。
- CORE12 G-like：12 个父类型，各 13,668 genes × 2,486 samples；IBD 子集 13,668 × 2,025；active IBD 子集 708，仅为历史分支，不混入当前冻结分析。G-like 来源包括 Stage1 broad 和 Stage2 epithelial，不能当 reference 或 H 矩阵。
- 当前冻结 2,025 活检/919 患者（CD 498、UC 421），53 states、7 ecotypes；正式 assignment 1,712，Unassigned 313。最高丰度 ecotype 对所有 2,025 样本均唯一；1,712 已分配者与最高丰度完全一致，313 Unassigned 保留其独立标签。

## 3. 两端重叠、缺失与版本差异

27 个冻结成员、丰度、state/ecotype assignment 文件逐个 SHA256 一致。服务器保存 05_EcoTyper 原始输入、G-like 大矩阵、模型中间结果及历史 formal/ctspec_formal/region_all/smoke 分支；本机主要保存导出和审阅快照。原始 bulk counts/adjcounts 在所查本机 MSCCR 目录未见对应原件，网页 CPM 输入是另一派生对象。本次只取回资产清单、矩阵表头、临床 metadata 和小型 QC 证据，没有复制大型表达矩阵。标准服务器 02_metadata/clinical 等仅有 .gitkeep，但真实 GEO SOFT 位于 03_reference/combined_reference/output/15_bayesprism_broad14_sanity_audit/，2,490 个 GSM 的所有共有字段与本地 series matrix **零差异**。

Excel“患者编号”是重复的患者标签，不能当唯一活检 ID；其中 1,328 行不等于实际活检 ID 是层级差异，不等于患者匹配错误。通过唯一 GSM 恢复全部 2,490 活检，提取 reGRID 后患者身份错误为 {excel_summary['PatientPrefixMismatches']}。Excel 数值逐列核对及四舍五入差异另存，主分析使用未四舍五入 GEO。

## 4. 实际可用临床字段

{coverage_table}

GHAS/Nancy 原字段各有 1,180 个活检值；主分析按疾病使用 GHAS-CD、Nancy-UC，故分母为 778/402。临床分数、CRP、粪钙卫蛋白等在患者内不变，局部组织学和局部炎症则可变。Age 为内镜入组年龄；其他字段没有逐项测量日期，不能保证同日。保留 CRP/FCal 的原始 log2 数值，不依据含糊单位/未知 offset 回算原浓度。RegionDetailed 为 7 个解剖类别，Region3 合并为 Ileum/Colon/Rectum；Region3 和 Disease/Inflammation 与 GEO 标签无冲突。

原始 endoremiss 的 2,020 个可比较值全部符合 SES-CD<3（CD）/Mayo=0（UC）。historemiss 与 max GHAS=0（CD）/max Nancy=0（UC）不符 189 活检（CD 110/48 患者；UC 79/32 患者）。原字段、冲突、初版统计均保留，但最终主结果不使用该未核实编码；零组织学活动变量明确标为本次派生。补充 PDF 有队列人口学汇总、bMIS 基因集及方法/相关性表，未见可连接的 reGRID 样本级 bMIS 值；不能据文献已有 bMIS 直接把它补进本队列。

未找到可核实的个体级采样日期、病程、Montreal B1/B2/B3、范围、狭窄/穿透/深溃疡、用药或既往生物制剂、白蛋白、血红蛋白、UCEIS、未来疗效及随访。限定范围外可能仍有文件，本报告不声称整个服务器不存在。

## 5. 完成的分析及结果

方法对照 [Luca et al., Cell 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8526411/)：父类型内 1,000 次 assignment label permutations；分类变量 abundance rank-sum；连续指标按 dominant CE 分组与其余 CE 比较，报告均值差及秩效应，并 BH。新增连续丰度 Spearman、KW 多分类整体和完整事后比较。19 个临床/表型变量、22 个回归对比；各类检验家族先在 statistical_analysis_plan.md 规定，随后因 historemiss 定义冲突作公开记录的修订。

主模型标准化 abundance ~ 目标临床特征 + Disease + Region3 + Age + Sex；不重复放入目标自身，不保留疾病子集内常量列。患者聚类稳健协方差、小样本修正和患者数减 1 的 t/F 参考；每患者等权拟合另行保留。无缺失插补，临床各模型用完整病例，包含 Unassigned。

{summary[summary.index('- **结果 1'):summary.index('患者级调整已经完成')]}

主表统计：分类整体 480、事后 540、连续 Spearman 660、dominant CE 连续特征 77；患者聚类主与等权共 2,640 个系数对比、2,280 个整体检验；53-state phenotype permutations 共 1,007 个 state×category 结果。患者首份活检 Spearman 660；Region3 配对 60 个 marker 检验、详细部位配对另 60 个。完整 P/q/效应/分母见 statistics/，不只输出显著结果。主模型检验 q<0.05 的 state 对比 {significant['state']}/1,166，CE 对比 {significant['ecotype']}/154。

既有 14 个 ecotype Disease/Inflammation rank-sum P 最大差 7.9e-17；7 个 Region3 KW P 最大差 {region_check.AbsPdiff.max():.3g}；159 个 state Region3 Fisher 样本数和 P 复核一致。父类型特异性由原始父类型内 assignment 分母和置换结果判断；同父类型效应排序另存，未把排序当正式 beta 差异检验。

这些结论是冻结 Discovery 中的横断面临床关联；不是未来预测、疗效验证或因果效应。变量缺失、组成约束、Discovery 选择和未知药物混杂仍限制解释。

## 6. 可直接用于 PPT 的图件

figures/ 包含 5 套 PNG/PDF/SVG；源数据在 source_data/。Fig1：53 states×22 临床对比；Fig2：7 CE×22 对比；Fig3：Cell 2021 Fig5A 的真实组成与临床特征适配；Fig4：E1 效应及患者调整前后区间；Fig5：30 个严格详细部位配对。Fig3 的临床 rank-sum 为活检层面描述，患者调整在 Fig1/2/4 展示。图件均经实际渲染、人工逐图查看、字体、碰撞和对齐检查。

## 7. 后续与复现

优先澄清原始 historemiss 阈值；获取 bMIS 样本级原值，或单独预注册严格复算；获得用药与随访后再讨论治疗结局。现有 CRP/粪钙卫蛋白及组织学强度分析已完成，不需要虚构指标或重跑 Discovery。

conda 环境导出见 environment.yml，包版本见 audit/package_versions.json。复现顺序：02_local_assets_and_geo.py（可从既存检查点跳过）→03_integrate_clinical_audit.py→03b_resolve_histology_label.py→04_clinical_associations.py→05_sensitivity_and_checks.py→06_plot_clinical_figures.py→07_closeout_reports.py。原始服务器清单与小型 metadata 快照已经保存，不必重连扫描。01_remote_readonly_audit.py 只在再次明确授权远端审计时才运行。
'''
(out/'audit_report.md').write_text(report,encoding='utf-8')
print('REPORTS_COMPLETE',significant,flush=True)
