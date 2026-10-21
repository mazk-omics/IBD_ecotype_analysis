# 根据实际统计表生成正式科研报告与中文汇报摘要，不手填显著性。
# 在服务器ibd_outcomes_20261015环境或含相同表格的本机目录运行。
from pathlib import Path
import pandas as pd
import numpy as np
import json, shutil

# Step1. 整合全部效应、患者分母和既有MSCCR横断面背景。
out = Path(__file__).resolve().parents[1]
stats = pd.read_csv(out / '04_cross_cohort_summary/all_associations.tsv', sep='\t')
denom = pd.read_csv(out / '04_cross_cohort_summary/endpoint_denominators.tsv', sep='\t')
eco = stats.loc[stats.Family.eq('Ecotype')].copy()
eco['Direction'] = np.where(eco.OR.gt(1), 'Higher score: higher endpoint odds', 'Higher score: lower endpoint odds')
context = pd.read_csv(out / 'audit/MSCCR_cross_sectional_context.tsv', sep='\t')
inflammation = context.loc[context.Feature.eq('Inflammation'), ['Marker', 'Beta', 'CI_low', 'CI_high', 'Q']].copy()
inflammation.columns = ['Feature', 'MSCCR_Inflammation_Beta', 'MSCCR_Inflammation_CI_low', 'MSCCR_Inflammation_CI_high', 'MSCCR_Inflammation_Q']
cross = eco.merge(inflammation, on='Feature', how='left', validate='many_to_one')
cross.to_csv(out / '04_cross_cohort_summary/ecotype_effects_across_endpoints.tsv', sep='\t', index=False)
hits = stats.loc[stats.Q_Wilcoxon.lt(.05) | stats.Q_Firth.lt(.05)].copy()
hits.to_csv(out / '04_cross_cohort_summary/all_family_adjusted_candidates.tsv', sep='\t', index=False)
stats.loc[~stats.Status.eq('OK')].to_csv(out / '04_cross_cohort_summary/model_failures_or_nonevaluable.tsv', sep='\t', index=False)
folders = {'GSE72819': '01_GSE72819_treatment_outcome', 'GSE93624': '02_RISK_disease_progression', 'GSE73661': '03_GSE73661_mucosal_healing'}
for accession, folder in folders.items():
    expression_qc = pd.read_csv(out / 'audit/expression_QC.tsv', sep='\t')
    if accession in expression_qc.Accession.values:
        expression_qc.loc[expression_qc.Accession.eq(accession)].to_csv(out / folder / 'expression_QC.tsv', sep='\t', index=False)
    m = pd.read_csv(out / folder / 'cohort_endpoint_manifest.tsv', sep='\t')
    m['ExclusionReason'] = m.ExclusionReason.fillna('')
    if accession == 'GSE73661':
        m['RecoveryStandardizationPool'] = '67 W0 biopsies only; all post-treatment expression excluded before scaling'
    for endpoint in m.Endpoint.unique():
        d = denom.loc[denom.Endpoint.eq(endpoint)].iloc[0]
        ix = m.Endpoint.eq(endpoint)
        m.loc[ix, 'EndpointEffectivePatientN'] = int(d.N)
        m.loc[ix, 'EndpointEvents'] = int(d.Events)
        m.loc[ix, 'EndpointNonevents'] = int(d.Nonevents)
    m.to_csv(out / folder / 'cohort_endpoint_manifest.tsv', sep='\t', index=False)
long_path = out / '03_GSE73661_mucosal_healing/all_samples_longitudinal_mapping.tsv'
long = pd.read_csv(long_path, sep='\t')
long['ExpressionAnalysisUse'] = np.select([long.Timepoint.eq('W0'), long.Timepoint.eq('CONTROL')], ['W0 baseline candidate; see drug-visit eligibility manifest', 'Non-IBD control; excluded from treatment association'], default='Post-treatment expression excluded from baseline regression; Mayo may define a matched follow-up endpoint')
long.to_csv(long_path, sep='\t', index=False)

# Step2. 记录预先固定的检验家庭、拟合参数与科学解释边界。
methods = '''# 统计方法与分析边界

分析对象为冻结MSCCR_CORE12_CTSPEC_REGION_IBD的12父类型、53 states和7 ecotypes。GSE72819复用既有完整Recovery。GSE73661与GSE93624使用同一冻结W、状态映射及成员表进行官方bulk Recovery；仅求新样本H，不重训W、不进行Discovery、不重跑BayesPrism。

连续state评分是官方state_abundances.txt中的保留H：H先在该父类型的全部初始状态内归一化，再按冻结mapping保留。保留状态之和可能小于1，不二次归一化。连续ecotype评分是官方ecotype_abundance.txt：先对冻结成员state H取均值，再在7生态型间归一化。已逐数值复核该公式。二者是相对表达程序系数，不是绝对细胞比例，不是Recovery Z，也不是dominant assignment。离散未分配不构成连续评分分析的自动排除理由。

全部结局只使用治疗前/初诊活检评分。GSE73661在Recovery开始前仅保留67份W0活检，剔除所有治疗后及非IBD对照表达；跨样本基因标准化、固定W投影与回归均不使用治疗后表达。原始GEO log2 RMA的作者预处理可能曾联合处理样本，缺乏其独立基线预处理结果，这是存档数据的限制；本轮没有把治疗后数据再次带入Recovery标准化。该设计是回顾性冻结程序关联，不是可对单个未来患者直接应用的预测器。

每一临床终点分别对7 ecotypes、53 states分析，构成相互独立报告的检验家庭。两侧Mann–Whitney/Wilcoxon采用渐近、ties校正、连续性校正；报告中位差和rank-biserial效应。每SD的Firth Logistic采用完整设计矩阵的Jeffreys惩罚：l(beta)+0.5 log|X'WX|。BFGS最大1000次迭代、gtol=1e-8，接受阈值为最大调整得分绝对值<=1e-5。95%CI由惩罚似然剖面与chi-square(1) 3.841459阈值求根；双侧P为同一完整模型惩罚下的似然比检验。未用Wald区间冒充剖面区间。SD来自本终点有效完整病例，ddof=1。

Wilcoxon与Firth的P分别按固定家庭大小7或53进行BH。失败/常数特征保留为不可评价，不缩小家庭大小。没有跨药物、时间窗或队列合并P值，也没有控制所有6终点整体的家族错误率。主要评分比较以Wilcoxon q为依据，Firth提供效应与支持性推断；所有结果并列，不以最小P选模型。未优化cutoff、未联合拟合7生态型、未做特征筛选、未计算训练AUC。

预定小样本规则：任一结局组<3时仅描述；3–5例时仅作非常不稳定的探索，保留全部宽区间。ETR有12例缓解、IFX有8例愈合、VDZ W6有6例愈合，未据显著性挑选协变量。其未调整结果不能排除疾病活动度、既往治疗和剂量混杂。RISK对全部7生态型另作年龄+性别敏感性，作为单独7检验家庭；不使用随访治疗或行为作为协变量。RISK没有逐人事件日期、删失日期及末次随访，采用论文定义的三年固定终点，不拟合Cox。非IBD控制不进入进展回归。

本轮RISK及GSE73661 W0-only新Recovery子进程均使用R_PROFILE_USER固定种子20261015。完整H与最终评分均已冻结保存，统计从这些文件复现。GSE73661最初含随访样本的Recovery已明确作废并隔离，不参与任何最终统计、图件或候选结论。正式报告以W0-only连续评分及最终统计表为准。
'''
(out / 'statistical_analysis_plan.md').write_text(methods, encoding='utf-8')

# Step3. 各队列报告保留完整7生态型表及全部53状态表入口。
notes = {
'GSE72819': '''73份治疗前UC结肠活检对应73个不同GEO inventory sample name（HF编号）。原始SOFT字段是“remission at week 10”，不是泛指clinical response。原文试验层面定义为总Mayo<=2且任何分项<=1；本队列缺少逐人Mayo供复算。70人有结局：12 Remitter、58 Non-remitter，3缺失。不能将3缺失视作未缓解。

全部73样本有53-state/7-ecotype连续评分，65/73有离散ecotype分配。70个有效结局中63个有离散分配，另7个虽离散Unassigned仍有有效连续评分，因此主分析N=70。历史“52/53状态恢复”涉及离散胜出覆盖，并不表示缺少第53个连续H。原始RPKM及既有冻结Recovery直接复用。

GEO总体设计指向etrolizumab治疗患者；没有可核验的逐样本placebo或剂量臂字段，DoseArm全部缺失，不推测100 mg/300 mg。既往anti-TNF资料完整：naive 27人（9缓解、18未缓解），既往暴露46人（3缓解、40未缓解、3结局缺失）。既往暴露与终点差异显著的分布是潜在混杂；本轮不作高维调整，不能称生态型具有独立疗效预测作用。组织诊断No Dx/NOS等不当作可靠基线炎症评分。

7生态型和53 states均无BH q<0.05信号。E5只有未经校正的候选趋势，不能写成已验证缓解标志物。E1横断面炎症关联不能推导为本队列不缓解预测。''',
'GSE93624': '''GEO及原文一致报告210位初诊、未经治疗的CD患者与35位非IBD控制；三年终点为27例B2/B3进展和183例未进展。原文选择的是诊断后至少前90天保持B1者，随后17例B2和10例B3。逐人基线行为、事件日期、删失日期和末次随访日期未公开，不能独立复核每个人的完整观察时长；本轮固定终点依赖该发表队列定义，不构造生存时间。Sample1..245是公开匿名研究个体标识，不是经核实可跨研究连接的患者ID。

原始GEO TAR逐样本就是有小数和少量负值的已处理表达，不是raw counts。服务器旧审计同样判为“exact transformation unresolved”。原文描述raw counts经edgeR TMM，另有逐基因inverse-rank-normal和eQTL/SNM分支；GEO存档没有唯一说明所处步骤。存档的基因均值/方差不符合已经逐基因标准正态的矩阵，不能宣称它就是inverse-normal或精确log2 CPM。标准化步骤/逆变换不可凭数值范围猜测。

本轮逐个原始样本核对存档数值，并保留现有权威基因映射后的13323×245矩阵；648个负值全部保留。不截断、不平移、不取绝对值、不指数逆变换，也不把已处理小数喂给counts模型再次TMM。官方Recovery对带负值或max<=50的输入不再log，随后逐基因标准化并按正负分量展开，固定W投影在数值上兼容该矩阵。该运行仅形成“存档处理尺度下的条件性关联”，不能替代来源完整、非负表达/原始counts的独立验证。

全部245样本得到完整连续评分，210个CD样本进入进展关联。GSE57945属于RISK相关数据，但公开CCFA_Risk编号与GSE93624 Sample编号缺少患者级对照表。既有表达指纹有152个一致候选，且不是零误差的同一表达实验复现，不能当作身份证明、不能转移27个进展标签，也不能当独立训练/验证队列。两者可能来自同一研究人群的不同实验或处理版本，精确重叠数及是否同一活检尚未证实。

E2和Glial_S03出现队列内校正信号。E2按年龄和性别调整后仍为正向，但7检验家庭的q>0.05，证据减弱。没有足够公开基线局部炎症、深溃疡及疾病范围，不能证明关联超越炎症。Glial是冻结父类型标签，本轮bulk程序投影不能单独证明真实胶质细胞数量或细胞特异性机制。''',
'GSE73661': '''原始178份活检包括：IFX分支46份（23个W0、23个W4–6配对），VDZ试验分支120份和12份非IBD控制。VDZ试验44个独立W0，41人基线诱导方案以vdz开头，3人为plac_plac，不纳入VDZ主要分析；vdz4w属于VDZ诱导，vdz_plac属于VDZ诱导后placebo维持，均不以随访治疗重新定义基线暴露。

主要结局为相应随访活检的Mayo内镜亚评分0或1。IFX W4–6为23人（8愈合、15未愈合）；VDZ W6为27人（6/21）、W12为13人（3/10）、W52为19人（13/6）。同一患者可出现在多个分别定义的VDZ终点中，但每个终点仅有一个W0评分。不能把27+13+19当59个独立患者。各随访的可评价病例集不同，缺失随访活检也可能包含选择性取样或退出，不能将缺失者计为未愈合。

表达为GEO明确记录的Affymetrix Human Gene 1.0 ST log2 RMA，既有权威映射18407×178，范围2.374–13.597。直接输入官方max<=50分支，不进行任意指数逆变换。新建Recovery在任何基因标准化前即仅保留67份W0，最终评分完整；所有治疗后表达均未参与本轮评分。后治疗表达没有用来预测已经发生的结局；本轮没有额外进行纵向变化检验。

IFX的E3与较高愈合几率相关；VDZ W6的E2与较低愈合几率相关，二者Wilcoxon和Firth分别BH均支持。VDZ W6的E6仅Wilcoxon家庭校正通过，而Firth的q不支持，须保留方法间不一致。两者样本和事件数很小，区间宽，不是可用的临床预测器。VDZ W12/W52无7生态型家庭校正后信号。E2在W12的未经校正方向与W6相反，不能据此证明时间效应翻转，因为两个终点来自不同且很小的病例集。基线Mayo实际可用，但6–8个愈合事件不足以稳健区分生态型、初始炎症与治疗方案的独立影响，未按显著性追加调整。'''
}
for accession, folder in folders.items():
    lines = [f'# {accession}：冻结生态型与未来结局关联结果', '', '分析日期：2026-10-15。完整病例、单特征、未经临床调整的关联；无预测AUC。', '', '## 队列及输入审计', '', notes[accession], '', '## 全部7生态型的实际结果', '', '|终点|生态型|N/事件|OR / SD|95%CI|Firth P|Firth q|Wilcoxon q|', '|---|---|---:|---:|---|---:|---:|---:|']
    for _, row in eco.loc[eco.Accession.eq(accession)].iterrows():
        lines.append(f'|{row.Endpoint}|{row.Feature}|{int(row.N)}/{int(row.Events)}|{row.OR:.3g}|[{row.CI_low:.3g}, {row.CI_high:.3g}]|{row.P_Firth:.4g}|{row.Q_Firth:.4g}|{row.Q_Wilcoxon:.4g}|')
    lines += ['', '## 状态探索性结果', '', '全部53 states及父谱系、效应方向、零值/常数审计、OR、剖面95%CI、P和q见 `cell_state_associations_complete.tsv`。各终点单独按53检验BH。']
    selected = hits.loc[hits.Accession.eq(accession) & hits.Family.eq('State')]
    if selected.empty:
        lines += ['', '没有state在任一单独报告的检验家庭内达到q<0.05。']
    else:
        for _, row in selected.iterrows():
            lines += ['', f'- {row.Endpoint}：{row.Feature}（父谱系{row.ParentLineage}；{row.Ecotype}），OR={row.OR:.3g}，95%CI [{row.CI_low:.3g}, {row.CI_high:.3g}]；Firth q={row.Q_Firth:.4g}，Wilcoxon q={row.Q_Wilcoxon:.4g}。']
    if accession == 'GSE93624':
        sensitivity = pd.read_csv(out / folder / 'ecotype_age_sex_sensitivity.tsv', sep='\t').set_index('Feature').loc['E2']
        lines += ['', f'E2年龄/性别敏感性：OR={sensitivity.OR:.3g}，95%CI [{sensitivity.CI_low:.3g}, {sensitivity.CI_high:.3g}]，P={sensitivity.P_Firth:.4g}，q={sensitivity.Q_Firth:.4g}；保持q=0.05阈值，不称调整后显著。']
    lines += ['', '## 方法、文件与限制', '', '统计定义与参数见根目录 `statistical_analysis_plan.md`。所有主模型保留独立失败记录；当前360个模型均完成。完整清单见 `cohort_endpoint_manifest.tsv`，临床字段见 `clinical_covariate_missingness.tsv`，连续评分与Recovery QC在本目录。图件在本目录figures及跨队列目录；全部53-state×6终点图见跨队列Fig3。所有缺失分母明确报告，未优化阈值，也未把不显著解读为不存在生物学效应。', '', f'原始来源：[GEO {accession}](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession})。']
    if accession == 'GSE93624':
        lines += ['', '原文：[Marigorta等，2017，PMCID PMC5745037](https://pmc.ncbi.nlm.nih.gov/articles/PMC5745037/)。GSE57945来源：[GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE57945)。']
    if accession == 'GSE72819':
        lines += ['', '缓解定义来源：[Tew等，2016，PMID 26522261](https://pubmed.ncbi.nlm.nih.gov/26522261/)。']
    (out / folder / 'Results_report.md').write_text('\n'.join(lines) + '\n', encoding='utf-8')

# Step4. 跨队列不合并P值，重点报告方向不一致和未通过校正的候选。
lines = ['# 冻结IBD生态型与外部未来结局：研究报告', '', '2026-10-15；服务器完成统计与Recovery。MSCCR横断面结果仅作生物学背景，不重新训练或用于外部筛选。', '', '## 有效分母', '', '|队列/终点|患者N|事件|非事件|', '|---|---:|---:|---:|']
for _, row in denom.iterrows():
    lines.append(f'|{row.Accession} / {row.Endpoint}|{int(row.N)}|{int(row.Events)}|{int(row.Nonevents)}|')
lines += ['', '## 主要结果', '', 'GSE72819的70个有效基线病例中，7个ecotypes和53 states均没有BH校正后显著关联。GSE73661的早期IFX黏膜愈合与E3、VDZ第6周黏膜愈合与E2相关。RISK三年进展与E2、其成员Glial_S03相关，但属于处理尺度溯源未完全解决的条件性分析；年龄/性别敏感性进一步减弱了E2的统计支持。']
for endpoint, feature in [('IFX_W4_6_endoscopic_healing', 'E3'), ('VDZ_W6_endoscopic_healing', 'E2'), ('RISK_3yr_B2_B3', 'E2')]:
    row = eco.loc[eco.Endpoint.eq(endpoint) & eco.Feature.eq(feature)].iloc[0]
    lines += ['', f'- {endpoint}，{feature}：OR/SD={row.OR:.3g}，95%CI [{row.CI_low:.3g}, {row.CI_high:.3g}]；Firth P={row.P_Firth:.4g}、q={row.Q_Firth:.4g}；Wilcoxon q={row.Q_Wilcoxon:.4g}。']
lines += ['', '## 与MSCCR炎症关联的区别', '', 'MSCCR中E1与局部炎症显著相关，而本轮所有未来结局中E1均未通过7检验家庭的BH。MSCCR中E5/E6随炎症降低；本轮E5无校正后信号，E6在VDZ W6的Wilcoxon家庭有支持而Firth家庭不支持，不能把两种方法说成一致验证。因此横断面炎症轴不能自动改写为不应答预测轴。', '', 'E2在MSCCR二分类炎症中缺乏支持（q约0.408），本轮在VDZ W6及RISK三年进展出现关联。虽然OR分别小于1和大于1，较高E2对应更少愈合及更多进展，在不利结局方向上相容；两种终点不能当作同一疗效终点的独立复现。E3在IFX愈合中有信号，而VDZ W52的名义方向相反但未通过校正；不能把它解释成跨药物通用疗效预测。E2在VDZ W12与W6方向不一致，同时W12仅3例愈合且病例集不同。E5在ETR缓解的名义趋势也未达到7生态型家庭q<0.05。', '', '完整42个生态型-终点效应、区间、q以及既有MSCCR炎症beta列于 `ecotype_effects_across_endpoints.tsv`；不同指标、疾病与统计单位不直接合并，也不对原始P作跨队列合并。', '', '## 输入、混杂及外推限制', '', 'RISK的公开表达不是counts，负值已在GEO逐样本存档中。官方带符号标准化投影数值上兼容，但确切存档变换仍未获完整来源证明。其信号须在有原始counts/明确非负表达及患者独立性的队列复核。GSE57945的152个指纹候选不能证明独立性或患者身份，未进行标签转移。', '', 'GSE73661分别按药物与时间窗处理；基线VDZ病例只有41人，各时间点的有效分母不同，不能将重复患者相加。其6–8例早期愈合事件及宽CI不足以支持部署性结论。ETR既往anti-TNF分布与缓解差异明显；剂量臂和可复算疾病活动度缺失，未作可靠混杂调整。RISK年龄/性别调整后的E2 q约0.061，不满足原定0.05。', '', '7 ecotypes及53 states是相对、彼此相关的程序评分，不是细胞绝对丰度。不同终点没有共同的临床含义：UC缓解/愈合为有利事件，CD进展为不利事件。无AUC、无自动筛选、无优化cutoff、无共同7特征模型。', '', '## 下一步独立验证需求', '', '1. 获取RISK存档生成代码或对应raw counts/可信RPKM，固定基因注释及样本ID后重新核对输入；不能猜测逆变换。', '2. 获取GSE57945与GSE93624的作者级患者/活检对照表、实验批次和测序对应关系；在此之前不宣称独立验证。', '3. 针对IFX E3及VDZ W6 E2预先登记同药物、同黏膜愈合定义、同时间窗的独立患者验证，并收集基线Mayo、既往生物制剂与诱导剂量。', '4. 对RISK E2/Glial_S03收集逐人90天基线行为、B2/B3事件时间、删失与治疗资料，先确认独立性和有效随访；有可靠时间后再考虑Cox。', '5. 用足量事件的前瞻性队列评价增量价值与校准；本轮不能替代预测模型开发、外部性能验证或机制实验。', '', '## 文件与图件', '', '四张主图均有PNG/PDF/SVG及源数据。Fig1：UC主要终点全部7生态型；Fig2：RISK与VDZ后续时间点；Fig3：全部53 states×6终点；Fig4：候选全部患者连续评分。原始GEO/论文来源及队列细节见各队列Results_report.md。']
(out / '04_cross_cohort_summary/Research_report.md').write_text('\n'.join(lines) + '\n', encoding='utf-8')

ifx_e3 = eco.loc[eco.Endpoint.eq('IFX_W4_6_endoscopic_healing') & eco.Feature.eq('E3')].iloc[0]
vdz_e2 = eco.loc[eco.Endpoint.eq('VDZ_W6_endoscopic_healing') & eco.Feature.eq('E2')].iloc[0]
vdz_e6 = eco.loc[eco.Endpoint.eq('VDZ_W6_endoscopic_healing') & eco.Feature.eq('E6')].iloc[0]
risk_e2 = eco.loc[eco.Endpoint.eq('RISK_3yr_B2_B3') & eco.Feature.eq('E2')].iloc[0]
ifx_states = stats.loc[stats.Endpoint.eq('IFX_W4_6_endoscopic_healing') & stats.Family.eq('State') & stats.Q_Firth.lt(.05)]
summary = f'''# 课题汇报：外部治疗结局与CD疾病进展

已实际完成三个队列、6个分别定义的结局、360个单特征关联模型（42 ecotype + 318 state），另有7个RISK年龄/性别敏感性模型。所有未来结局只使用基线/初诊评分。

|问题|实际结果|汇报边界|
|---|---|---|
|ETR第10周临床缓解|70人，12缓解；7 ecotypes和53 states均无BH显著信号|不能把“炎症相关E1”说成“不应答预测E1”|
|IFX W4–6内镜愈合|23人、8愈合；E3 OR={ifx_e3.OR:.3g}，95%CI [{ifx_e3.CI_low:.3g}, {ifx_e3.CI_high:.3g}]，Firth q={ifx_e3.Q_Firth:.4g}、Wilcoxon q={ifx_e3.Q_Wilcoxon:.4g}|小样本、宽CI、未可靠调整基线Mayo与既往治疗|
|VDZ W6内镜愈合|27人、6愈合；E2 OR={vdz_e2.OR:.3g}，95%CI [{vdz_e2.CI_low:.3g}, {vdz_e2.CI_high:.3g}]，Firth q={vdz_e2.Q_Firth:.4g}、Wilcoxon q={vdz_e2.Q_Wilcoxon:.4g}|6事件且CI很宽；W12/W52无校正后信号|
|VDZ W6的E6|OR={vdz_e6.OR:.3g}，95%CI [{vdz_e6.CI_low:.3g}, {vdz_e6.CI_high:.3g}]；Wilcoxon q={vdz_e6.Q_Wilcoxon:.4g}，Firth q={vdz_e6.Q_Firth:.4g}|仅秩检验家庭支持，Firth不支持；不能省略检验方法差异|
|RISK三年B2/B3进展|210位CD、27进展；E2 OR={risk_e2.OR:.3g}，95%CI [{risk_e2.CI_low:.3g}, {risk_e2.CI_high:.3g}]，Firth q={risk_e2.Q_Firth:.4g}；Glial_S03 q约0.0422|存档表达尺度溯源有保留；E2年龄/性别调整后q约0.0613|

状态候选：IFX的Enteroendocrine_S02、Glial_S06、Pericyte_S03（均属于冻结E3；53状态家庭Firth q范围{ifx_states.Q_Firth.min():.4g}–{ifx_states.Q_Firth.max():.4g}），以及RISK的Glial_S03（E2）。IFX的CD8_T_S04（E6）仅Wilcoxon q<0.05而Firth q>0.05，单独标记。它们是bulk冻结程序评分，不能直接写成这些细胞数量改变或机制已确立。

最适合汇报的结论：**未来结局关联表现出药物和终点特异性；E2/E3提供了与E1横断面炎症轴不同的候选信息，但证据仍需独立验证。** E5没有校正后显著信号；E6只在VDZ W6的Wilcoxon家庭中有支持，Firth不支持。

必须明确：最终GSE73661在任何Recovery标准化前仅保留67份W0，首次包含治疗后样本的评分/统计/图件已作废隔离，不进入本交付。不同药物/终点/时间窗没有合并P；GSE57945与GSE93624未证明患者独立；随访缺失没有当作未愈合。没有训练AUC、没有优化cutoff、没有联合7特征预测模型。RISK是已运行的条件性探索结果，不能隐藏其输入溯源限制。

PPT建议：Fig1→Fig2→Fig4；Fig3作为全部53状态备份。详细准确值以完整统计TSV和正式Research_report.md为准。
'''
(out / 'presentation_clinical_outcomes_summary.md').write_text(summary, encoding='utf-8')
(out / 'README.md').write_text('''# 外部临床结局分析交付

入口：`presentation_clinical_outcomes_summary.md`；正式整合报告：`04_cross_cohort_summary/Research_report.md`。

三个队列各有患者-样本-终点manifest、临床缺失审计、连续评分、Recovery QC、完整7生态型/53状态统计与Results_report.md。跨队列表保留全部360个模型，RISK另有年龄/性别敏感性。原始大型数据/模型不在本交付中复制，服务器源路径及哈希见audit。

主图位于`04_cross_cohort_summary/figures/`；PNG用于汇报，PDF/SVG用于编辑。RISK为有存档处理尺度保留的条件性结果。方法、复现与安全边界见statistical_analysis_plan.md。服务器分析环境为ibd_outcomes_20261015，完整environment.yml锁定依赖。临床统计可从服务器保留的Recovery输出及manifest复算，无需再次执行Recovery；本地包同时提供完整连续评分供独立核对。
''', encoding='utf-8')
print('REPORTS_COMPLETE', len(stats), flush=True)
