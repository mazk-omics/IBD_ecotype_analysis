from pathlib import Path
import json, hashlib
import pandas as pd
import numpy as np

OUT=Path(__file__).resolve().parents[1]
ext=pd.read_csv(OUT/'statistics/external_all60_rankings.tsv',sep='\t')
disc=pd.read_csv(OUT/'statistics/discovery_one_biopsy_all60_rankings.tsv',sep='\t')
allres=pd.concat([ext,disc],ignore_index=True)
denom=pd.read_csv(OUT/'audit/endpoint_comparability.tsv',sep='\t')
groups=pd.read_csv(OUT/'statistics/ecotype_state_rank_distribution_descriptive.tsv',sep='\t')
boot=pd.read_csv(OUT/'statistics/discovery_E1_member_effect_bootstrap.tsv',sep='\t')
sens=pd.read_csv(OUT/'statistics/discovery_100_selection_sensitivity.tsv',sep='\t')
short={'ETR_W10_clinical_remission':'ETR W10','IFX_W4_6_endoscopic_healing':'IFX W4–6','VDZ_W6_endoscopic_healing':'VDZ W6','RISK_3yr_B2_B3':'RISK 3年 B2/B3','VDZ_W12_endoscopic_healing':'VDZ W12（补充）','VDZ_W52_endoscopic_healing':'VDZ W52（补充）','MSCCR_Inflammation_one_biopsy':'MSCCR Discovery 炎症（一人一活检）'}

def fmt(value):
    if pd.isna(value):return 'NA'
    if isinstance(value,(int,np.integer)):return str(value)
    if isinstance(value,(float,np.floating)):
        return f'{value:.3g}' if abs(value)<.001 else f'{value:.4f}'.rstrip('0').rstrip('.')
    return str(value)

def md(frame):
    columns=list(frame.columns)
    lines=['| '+' | '.join(columns)+' |','| '+' | '.join(['---']*len(columns))+' |']
    lines += ['| '+' | '.join(fmt(v) for v in row)+' |' for row in frame.itertuples(index=False,name=None)]
    return '\n'.join(lines)

def write(path,text):
    (OUT/path).write_text(text,encoding='utf-8')

signed=ext.copy()
signed['OriginalDirectionRank']=signed.Rank_signed_favorable
signed.sort_values(['Endpoint','OriginalDirectionRank','Feature']).to_csv(OUT/'statistics/paper_direction_rankings_separate_endpoints.tsv',sep='\t',index=False,float_format='%.17g')
signed[signed.OriginalDirectionRank<=10].to_csv(OUT/'statistics/paper_direction_top_ranks.tsv',sep='\t',index=False,float_format='%.17g')

primary=ext[ext.Role.str.startswith('Primary')]
e1=disc[disc.Feature.eq('E1')].iloc[0]
members=disc[disc.E1_member]
s_e1=sens[sens.Feature.eq('E1')]
firstn=int((s_e1.Rank_absZ==1).sum())
best=members.sort_values('Rank_absZ').iloc[0]
bestboot=boot[boot.Feature.eq(best.Feature)].iloc[0]
clustered=pd.read_csv(OUT/'statistics/discovery_historical_clustered_beta_separate_model.tsv',sep='\t')
clustered=clustered[clustered.Marker.isin(['E1']+members.Feature.tolist())]
sections=['# 冻结53 states / 7 Ecotypes：统一临床关联排名（2026-10-16）',
'本轮已经实际重算：6个外部终点×60项=360个双侧检验；Discovery一人一活检60项；100轮患者内活检抽样共6000项敏感性结果；E1与10成员的2000轮配对特征bootstrap。原始冻结文件只读，所有新结果独立输出。',
'## 结论',
'外部治疗/进展结局中，没有看到Ecotype在整体排名或效应分布上持续优于单一state。ETR两类中位排名接近；IFX、VDZ W6、RISK的Ecotype中位排名略靠后。这里保留了state胜出的结果。',
f'Discovery内部炎症比较中，E1排第1，rank-biserial={e1.RankBiserial:.4f}；十个成员排第2–11，效应范围{members.RankBiserial.min():.4f}–{members.RankBiserial.max():.4f}。与最强成员{best.Feature}的效应差只有{e1.RankBiserial-best.RankBiserial:.4f}，配对bootstrap 95%区间[{bestboot.Delta_CI_low:.4f}, {bestboot.Delta_CI_high:.4f}]。因此可说E1聚合略靠前，但不能据此声称显著优于每个成员。',
f'100轮一人一份敏感性中，E1有{firstn}/100轮排第1，其余轮次排第2。这支持排名对活检选择较稳定，同时效应增益仍小。该证据完全来自同一Discovery，不构成外部验证或增量预测验证。',
'## 公平比较与实际分母',
md(denom[['Endpoint','N','Events','Nonevents','EffectiveN_min','EffectiveN_max','NA_cells','ExactZeroCells','TiedObservationExcess']].rename(columns={'Endpoint':'终点','Events':'事件','Nonevents':'非事件'})),
'每个终点的60项使用完全相同患者集合；每项有效N均为该终点N。六矩阵无NA、非有限值、精确0及样本评分ties。部分评分小于等于1e-8但保持原始小正数，不取整、不置零。评分ties为零不代表P值/名次无并列：小样本U统计量离散，等P保留平均名次；用于展示的同名次次序按Feature字典序，不能解释为优劣。逐特征near-zero、unique数及P并列人数均在TSV。',
'外部state为同一冻结W与mapping下的Recovery H；Ecotype为成员H均值后在七组间归一。三个队列公式误差均约5e-16。Discovery读取同一冻结原生H及其Ecotype聚合，全部107325个state值与14175个Ecotype值另与原生文件核对；不是把Discovery称作外部W投影。',
'## Ecotype与State群组分布（描述性）',
md(groups.assign(Endpoint=groups.Endpoint.map(short))[['Endpoint','Family','FeatureN','RankMedian','RankQ25','RankQ75','MedianAbsRB','Rank_le10_N','BH60_lt005']]),
'排名越小关联越强；MedianAbsRB为绝对rank-biserial中位数。不同特征相关、两类特征数7和53不同，因此不把这些特征当独立重复做组间显著性检验。Rank≤10计数保留所有并列；VDZ W6存在跨第10位的并列，所以可能有11项入列。未跨不同药物或终点计算总体Liptak，也未提供虚构的全局优势P值。',
'## 每个终点Top10、7个E及E1成员',
'下列Top10按|Z|排序，同时给出原文方向排名。RB正值表示event=1患者评分较高；ETR/IFX/VDZ的event=1是缓解/愈合，RISK的event=1是B2/B3进展，Discovery的event=1是Inflamed。AUC-like=U/(n1×n0)，是同批样本的两组随机配对排序概率（ties计半），并非样本外AUC。完整方向和两种名次均在全60项TSV。']
for ep,d in allres.groupby('Endpoint',sort=False):
    first=d.iloc[0]
    sections.extend(['### '+short[ep],f'N={int(first.N)}，事件={int(first.Events)}，非事件={int(first.Nonevents)}。'+('RISK只作表达处理尺度条件下的关联；E2年龄/性别Firth调整BH7 q=0.0613。' if ep.startswith('RISK') else '补充时间窗；不进入主终点结论。' if 'W12' in ep or 'W52' in ep else '同一患者与缺失口径。'),
       md(d.sort_values('DisplayOrder_absZ').head(10)[['Feature','Family','Rank_absZ','Rank_signed_event' if ep.startswith('MSCCR') else 'Rank_signed_favorable','RankBiserial','AUC_like_event','P_two_sided','Q_BH60']].rename(columns={'Feature':'Top10特征','Rank_absZ':'强度名次','Rank_signed_favorable':'有利方向名次','Rank_signed_event':'炎症正向名次','RankBiserial':'RB','P_two_sided':'P','Q_BH60':'BH60 q'})),
       '7个Ecotype（按强度名次）：'+ '；'.join(f'{r.Feature} #{r.Rank_absZ:g} (RB={r.RankBiserial:+.3f})' for r in d[d.Family.eq('Ecotype')].sort_values('Rank_absZ').itertuples())+'。',
       'E1及其全部10成员：',
       md(d[d.Feature.eq('E1')|d.E1_member].sort_values('Rank_absZ')[['Feature','Rank_absZ','RankBiserial','P_two_sided','Q_BH60']].rename(columns={'Feature':'特征','Rank_absZ':'强度名次','RankBiserial':'RB','P_two_sided':'P','Q_BH60':'BH60 q'}))])
sections.extend(['## Figure 5C原文核对及本轮适配',
'原文先在各免疫治疗数据集中独立估计基线特征，OS用单变量Cox，二分类反应用双侧Wilcoxon并将P转成Z；同治疗类型内按样本数加权Liptak；再对3治疗类型×2结局的6项有方向Z分别排名、取平均名次。Table S6J定义负号为较优结局。原文未披露P转signed Z的明确代数式、组输入顺序、exact/ties选项或最终并列规则。',
'本轮按预设要求比较任意方向关联强度，主名次为降序|Z|；外部终点另提供正值=有利结局的SignedZ_favorable与原文负值=有利的SignedZ_paper_style，以及对应有利方向名次。Discovery炎症不能定义为有利临床结局：这些临床方向列留NA，只提供event=Inflamed的正向Z和Rank_signed_event。主图以预设首终点ETR W10的|Z|排序统一行；其余列保留相同特征行，未生成异质终点综合分。所有终点各自提供完整排序图。',
'原文方法与图注：[Luca et al., Cell 2021全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC8526411/)，[作者稿PDF](https://pmc.ncbi.nlm.nih.gov/articles/PMC8526411/pdf/nihms-1740510.pdf)，[Table S6](https://pmc.ncbi.nlm.nih.gov/articles/instance/8526411/bin/NIHMS1740510-supplement-12.xlsx)。STAR Methods “Prediction of immunotherapy response”为PDF内页31–32，Figure 5图注为内页49。',
'## 统计口径和多重检验',
'所有60项统一采用渐近双侧Mann–Whitney，平均秩处理ties，方差校正ties，连续性校正开启。Z_abs=norm.isf(P/2)，方向=sign(U−n1n0/2)，SignedZ_event=方向×Z_abs；RB=2U/(n1n0)−1。基准P与历史Wilcoxon逐项一致，并用合并平均秩及校正方差独立复算核验。',
'每终点BH60为本轮新增探索性校正；历史Wilcoxon/Firth BH7、BH53及本轮重算家族q分别保留，主排名只取原始P/Z。Discovery历史F062/F063的BH分别为1166和154个contrast家族，另表保留，不能拿q跨模型排序。本轮没有跨6外部终点总体校正，也没有使用选择显著模型、优化cutoff或联合预测器。',
'## Discovery补充比较与限制',
f'主分析从全部2025活检/919患者中，按SHA256(20261016|PatientID|SampleID)最小值固定一份；不使用炎症或评分来选择。入组{int(e1.N)}人，Inflamed {int(e1.Events)}，Non-inflamed {int(e1.Nonevents)}。100轮敏感性使用numpy default_rng(20261016)，每轮同919人、同60项，各患者随机一份。所有91900个迭代sample ID均交付。',
'E1与十成员配对特征bootstrap保留同一患者的所有评分，按炎症组分层抽患者，2000轮，seed=20261017（随机种子数值，不是运行日期），提供效应及效应差的percentile 95%区间。这些区间未校正十成员比较，也不能消除Discovery共享数据、组成相关性和临床混杂。',
'患者聚类调整的历史F062/F063炎症beta另存statistics/discovery_historical_clustered_beta_separate_model.tsv；其Disease+Region3+Age+Sex调整及clustered t推断与本轮未调整MW是不同模型，不作统一名次。E1 beta约1.098，成员约1.064–1.103，支持共享炎症信号接近的解释。',
md(clustered[['Marker','Beta','CI_low','CI_high','P','Q','HistoricalBHFamilySize','N','Patients']].rename(columns={'Marker':'历史模型特征','Beta':'调整beta','Q':'原家族q','HistoricalBHFamilySize':'原BH家族大小'})),
'## 数据来源与审计边界',
'GSE73661先从既有178份deposited log2 RMA筛67份W0，随后才进入Recovery gene-wise标准化；不是重新只用67份做RMA。基线表达与随访Outcome对接。已存档作者RMA是否联合处理随访表达仍是限制。W6/W52与W12/W52分别共享10和8位患者，不能视为独立meta证据。VDZ manifest EndpointDefinition误附IFX时间窗的文字已在本轮新清单重写，原文件保留。',
'外部45个不可变输入的历史核验清单覆盖12个W及mapping等；三队列projection QC支持同一冻结模型。GSE72819本机保留原始Recovery并可逐值核对。GSE73661和RISK的原始RecoveryOutput本机缺失：本轮用冻结全精度派生矩阵、成员公式、交付hash和审计链；未声称在本机从原始Recovery再重建这两队列。RISK的deposited精确转换未闭环；E2调整后BH7未显著。',
'没有重跑BayesPrism/EcoTyper/NMF/Recovery，未改冻结结果，未git commit/push。input_SHA256_before_after.tsv核对本轮所有实际读取来源；另将本机存档12个原始W文件SHA与外部冻结模型清单逐一精确比对，12/12相符。统计运行环境和78项数值QA见audit。Nature-figure规范用于可编辑导出、面板对齐、5pt字号及PDF碰撞检查；仅继承用户参考图的热图/引线风格。',
'## 文件入口',
'主图：figures/Fig1_reference_style_external_heatmap（PNG/PDF/SVG）；群组比较：Fig2_rank_and_effect_distributions；Discovery：Fig3_discovery_E1_and_members。FigS1为晚期VDZ，FigS2_1至_6为每终点完整60项图，FigS3为Discovery完整60项。',
'完整表：statistics/external_all60_rankings.tsv、discovery_one_biopsy_all60_rankings.tsv、all_7_ecotype_rankings.tsv、E1_and_10_members_all_endpoints.tsv；方向排名另见paper_direction_rankings_separate_endpoints.tsv。manifests含真实sample/patient ID及逐患者60项评分。scripts/analyze.py、plot_figures.py、write_report.py和run_qa.py可复跑；建议使用本轮统计环境或同版本依赖。'])
write('中文分析汇报.md','\n\n'.join(sections)+'\n')
write('README.md','''# Cell 2021 Fig5C-inspired IBD association ranking

先读 `中文分析汇报.md`。本目录为独立输出，原始输入只读。

运行顺序（Python 3.13.5 / NumPy 2.1.3 / pandas 2.2.3 / SciPy 1.15.3 / matplotlib 3.10.0）：

```powershell
& 'D:\\anaconda3\\python.exe' scripts/analyze.py
& 'D:\\anaconda3\\python.exe' scripts/plot_figures.py
& 'D:\\anaconda3\\python.exe' scripts/write_report.py
& 'D:\\anaconda3\\python.exe' scripts/run_qa.py
```

分析脚本按绝对冻结项目位置读取输入。绘图QA使用已安装nature-figure脚本（复制到scripts/figure_qa后可离线复用）。无任何恢复、拟合NMF或重跑BayesPrism步骤。

方向：event=1为ETR缓解、IFX/VDZ愈合、RISK进展、Discovery炎症。`SignedZ_event`正号=评分在事件组更高。`SignedZ_favorable`对RISK反号；`SignedZ_paper_style`再反号以与原表负=较优结局一致。主强度排名是absZ，原文方向排名单列。AUC_like为同批两组排序概率；不解释为预测性能。

每终点独立BH60；historical q与family q分列。没有全局Liptak、跨终点合并P或样本外预测声明。RISK尺度未闭环，VDZ W12小样本，Discovery为内部关联。
''')
write('reference/原文规则与适配说明.md','''# Luca et al., Cell 2021: Figure 5C

证据来源为正文、STAR Methods “Prediction of immunotherapy response”、Figure 5图注及Table S6J。链接与页码见中文分析汇报。

1. 仅治疗前RNA-seq，特征在每个数据集独立估计。
2. OS采用单变量Cox的Z；二分类疗效用双侧Wilcoxon的P转Z。
3. 同治疗类型内使用Liptak加权Z，样本量为权重。标准重建式为sum(n_i*z_i)/sqrt(sum(n_i^2))；原文没有打印公式。
4. 三种治疗×两种结局的六项有方向Z各自排名，六名次平均。S6J负号为较优结局。
5. 未公开组次序、P转signed Z代数式、Wilcoxon exact/ties与最终名次并列规则。本轮明确声明这些实现选择。
6. 原文没有要求该排名做BH，也不是AUC或样本外判别能力评估。

本轮IBD终点药物、疾病和事件方向异质，晚期VDZ也存在患者重叠，因此不合并成Liptak/meta-Z或平均全局名次。分别提供|Z|强度排名与有利方向排名。热图共同行按预设首终点ETR W10强度排序，目的只是让相同特征对齐。

原文中Liptak括号的aPD1重复为排版笔误；正文/Fig5C六列对应anti-PDL1、anti-PD1、anti-CTLA4。本报告没有将该笔误当成统计定义。
''')
reference_paths=[Path(r'C:/Windows/Temp/nihms-1740510.pdf'),Path(r'C:/Windows/Temp/NIHMS1740510-supplement-12.xlsx')]
refs=[]
for p in reference_paths:
    if p.exists():refs.append(dict(Path=str(p),SHA256=hashlib.sha256(p.read_bytes()).hexdigest(),Role='Official paper reference, no biomarker values imported'))
pd.DataFrame(refs).to_csv(OUT/'reference/original_reference_SHA256.tsv',sep='\t',index=False)
print('Report written; Discovery inflamed',int(e1.Events),'E1 rank1 sensitivity',firstn,'/100')
