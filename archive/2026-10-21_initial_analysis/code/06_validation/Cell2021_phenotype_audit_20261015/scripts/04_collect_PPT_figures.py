# 将核验后的原图复制到独立汇报批次，保留来源、SHA256、中文图注和逐文件质量记录。
# 本机msccr_clinical_audit_20261015环境运行；不复制原始表达或患者级临床明细。
from pathlib import Path
from datetime import datetime
import hashlib,shutil,json,re,sys,subprocess
import numpy as np
import pandas as pd
import pymupdf
from PIL import Image,ImageDraw,ImageFont

# Step1. 创建全新批次，避免覆盖；明确每一组图的证据对象
root = Path('D:/IBD_ecotype/IBD_ecotype_analysis')
out = Path(__file__).resolve().parents[1]
batch = root.parent/'初步结果汇报/相关结果/20261015_BayesPrism_CORE12_53states_7ecotypes'
if batch.exists():
    batch = batch.with_name(batch.name+'_'+datetime.now().strftime('%H%M%S'))
batch.mkdir(parents=True,exist_ok=False)
folders = ['00_结果导航与清单','01_Reference_and_Discovery','02_CellState_Official_EcoTyper','03_Ecotype_Official_EcoTyper','04_CellState_Molecular_Characterization','05_Ecotype_Functional_Characterization','06_External_Recovery','07_MSCCR_Clinical_Associations','08_Treatment_Outcomes_and_RISK','09_Phenotype_Enrichment','10_Method_QC_and_Backup']
for folder in folders:
    (batch/folder).mkdir()
nav = batch/folders[0]
prior = root/'06_validation/paper_aligned_recovery_characterization_20261012'
paper = root/'06_validation/paper_figures_20261012'
clinical = root/'06_validation/MSCCR_joint_clinical_audit_20261015'
outcomes = root/'06_validation/external_clinical_outcomes_20261015'
frozen = root/'ctspec_region_ibd'
states = pd.read_csv(prior/'00_audit/frozen_states.tsv',sep='\t')
parent_coverage = pd.read_csv(out/'audit/discovery_parent_assignment_coverage.tsv',sep='\t').set_index('CellType')
registry = []

# Step2. 官方Discovery与最终7生态型，所有父类型原生热图保持原文件格式
for ct,g in states.groupby('CellType',sort=True):
    registry.append(dict(Source=frozen/ct/'state_assignment_heatmap',Module=folders[2],Order=3,Title=ct+'：冻结状态原生表达与分配热图',English=ct+' frozen cell-state expression and assignment',Main='Main PPT',Dataset='MSCCR',N=str(int(parent_coverage.loc[ct,'AssignedBiopsies']))+' assigned biopsies / 2025',Patients=str(int(parent_coverage.loc[ct,'Patients']))+' assigned patients / 919',States=str(len(g)),Ecotypes='NA',Method='Official EcoTyper expression/assignment heatmap; descriptive; no BH',Evidence=frozen/ct/'heatmap_top_ann.txt',Result='显示该父类型通过QC状态的marker表达及真实分配结构。',TakeHome='该父类型具有可追溯的冻结状态分配。',Limit='仅显示该父类型成功分配活检；不是疾病富集检验或细胞身份确证。'))
registry.append(dict(Source=frozen/'Ecotypes/heatmap_assigned_samples_viridis',Module=folders[3],Order=5,Title='7生态型与53状态的官方Discovery热图',English='Official discovery heatmap of seven ecotypes',Main='Main PPT',Dataset='MSCCR',N='1712 assigned / 2025 total biopsies',Patients='919 total; assigned subset in metadata',States='53',Ecotypes='7',Method='Official ecotype significant assignment; no additional hypothesis test',Evidence=frozen/'Ecotypes/ecotypes.txt',Result='最终模型为7生态型、53状态；图中显示1712份显著分配活检。',TakeHome='冻结模型的生态型组成和分配可以直接追溯。',Limit='313份Unassigned不在原生已分配热图中；使用P4同时展示其分母。'))
for stem,title in [('rank_plot','官方父类型秩选择QC'),('Ecotypes/jaccard_matrix','官方状态共现Jaccard矩阵'),('Ecotypes/nclusters_jaccard','官方Ecotype簇数QC')]:
    registry.append(dict(Source=frozen/stem,Module=folders[10],Order=12,Title=title,English=title,Main='Backup',Dataset='MSCCR',N='2025 input biopsies',Patients='919',States='96 initial / 53 retained; QC context',Ecotypes='7 final',Method='Native QC diagnostic; descriptive',Evidence=frozen/'Ecotypes/ecotypes.txt',Result='保留原生模型选择及共现诊断。',TakeHome='QC选择与最终53/7分配分别解释。',Limit='秩图涉及初始模型，不把96初始状态称为53个QC后状态。'))

# Step3. 分子、功能与既有方法图；旧版未聚类临床结果只留备份
for stem,title,order,module,status,evidence,claim,limit in [
 ('Fig2B_cell_state_UMAP','53状态的原文式UMAP展示',3,4,'Main PPT',paper/'source_data/Fig2_input_coverage.tsv','父类型内冻结状态的表达异质性可视化。','UMAP距离及分区不是新的状态验证或聚类结果。'),
 ('Fig4B_ecotype_member_networks','7生态型的状态成员网络',5,5,'Main PPT',paper/'source_data/Fig4B_all_within_ecotype_edges.tsv','展示既定成员及状态共现边。','共现边不表示配体受体通讯或因果。'),
 ('Fig5A_ecotype_feature_overview','历史385项Ecotype功能与表型概览',6,10,'Backup',paper/'source_data/Fig5A_all_385_feature_tests.tsv','350 Hallmark、14疾病/炎症及21肠段比较的历史概览。','未聚类调整，不作最终临床证据；不同模块统计对象不同。')]:
    registry.append(dict(Source=paper/'figures'/stem,Module=folders[module],Order=order,Title=title,English=stem.replace('_',' '),Main=status,Dataset='MSCCR',N='2025',Patients='919',States='53',Ecotypes='7',Method='UMAP / frozen network' if stem!='Fig5A_ecotype_feature_overview' else 'Mann-Whitney two-sided; BH385; biopsy dependence unadjusted',Evidence=evidence,Result=claim,TakeHome=claim,Limit=limit))
for stem,title,module,order,status,evidence,claim,limit in [
 ('Figure_A_recovery_Z','两队列53状态Recovery证据',6,7,'Main PPT','01_state_recovery/state_recovery_evidence_by_dataset.tsv','分队列保留Recovery Z和覆盖。','Recovery Z是名义程序恢复证据，不是表型富集一致性。'),
 ('Figure_B_LOOCV','外部状态marker的LOOCV',4,6,'Main PPT','02_external_markers/marker_LOOCV_per_dataset.tsv','共同可评价留出中结肠34/37、回肠疾病版31/32。','只按实际可评价分母报告；两数据集不能视为广泛独立验证。'),
 ('Figure_B_top10_overlap','两种分层策略top10 marker重合',10,12,'Backup','02_external_markers/n4_paired_LOOCV_comparison.tsv','48/53共同完整top10状态的重合中位数9。','对同一数据的替代分层比较不是外部样本复制。'),
 ('Figure_C_region_Fisher','53状态×3肠段Fisher富集',9,4,'Backup','04_state_phenotypes/state_region_fisher.tsv','159项数值与分母独立重算通过。','活检重复患者未控制；优先展示本轮患者敏感性P2。'),
 ('Figure_E_ecotype_Z','外部Ecotype共现Z',6,7,'Main PPT','05_ecotype_coassociation/ecotype_coassociation_1000permutations.tsv','E1和E3在两个队列的共关联Z较高；须结合成员覆盖。','不是完整成员无条件复现；原始共现Z的tail有保留，不添加P星号。'),
 ('Figure_F_Hallmark','7生态型的50个Hallmark',5,6,'Main PPT','06_ecotype_features/ecotype_continuous_features.tsv','展示生态型对应的分子功能程序和全部阴性项。','Hallmark来自同一bulk；不是独立功能实验，未控制重复患者。')]:
    registry.append(dict(Source=prior/'08_figures'/stem,Module=folders[module],Order=order,Title=title,English=stem.replace('_',' '),Main=status,Dataset='MSCCR and external scRNA' if module==6 else 'MSCCR / external markers',N='2025 discovery; external 9752 / 4818 cells',Patients='919 discovery; external 18 / 30 donors',States='53; dataset-evaluable coverage varies',Ecotypes='7',Method='Recovery nominal Z / marker LOOCV / coassociation Z' if module in [4,6] else ('Fisher greater BH159' if stem=='Figure_C_region_Fisher' else 'Mann-Whitney; BH350; descriptive'),Evidence=prior/evidence,Result=claim,TakeHome=claim,Limit=limit))
for dataset in ['Colon_GSE214695','Ileum_IBDverse']:
    registry.append(dict(Source=prior/'08_figures'/('Figure_D_'+dataset),Module=folders[9],Order=12,Title='旧外部表型图：混合对照定义，排除主汇报',English='Superseded phenotype display',Main='Excluded',Dataset=dataset,N='See old806 coverage',Patients='18 / 30',States='Variable',Ecotypes='NA',Method='Historical cell / sample dominant Z; no BH',Evidence=prior/'04_state_phenotypes/state_phenotype_enrichment.tsv',Result='旧计算可复现，但不能用于本轮严格可比表型结论。',TakeHome='以本轮P3替代。',Limit='健康对照与IBD非炎症混合；不复制进本批次。'))

# Step4. 两套scRNA原有正式图及marker备份；只引用外部明细，不将供者明细复制到PPT
for folder,dataset,n,donors,evaluable in [('GSE214695_scRNA_recovery_20261009','GSE214695',9752,18,53),('IBDverse_scRNA_recovery_20261009','IBDverse',4818,30,47)]:
    p = root/'06_validation'/folder
    for image in sorted((p/'figures').glob('*.png')):
        stem = image.stem
        status = 'Main PPT' if stem in ['Fig1_recovery_overview','Fig3_donor_state_support'] else 'Backup'
        titles = {'Fig1_recovery_overview':'状态Recovery与覆盖概览','Fig3_donor_state_support':'供者层状态支持','Fig4_lineage_marker_checks':'父谱系marker核查','Fig5_ecotype_projection':'Ecotype投影','Fig6_colon_ileum_comparison':'结肠与回肠Recovery比较'}
        title = dataset+'：'+titles.get(stem,stem.replace('Fig2_','').replace('_marker_recovery','')+' marker恢复')
        registry.append(dict(Source=image.with_suffix(''),Module=folders[6] if 'Fig2_' not in stem else folders[4],Order=7,Title=title,English=dataset+' '+stem.replace('_',' '),Main=status,Dataset=dataset,N=str(n)+' cells',Patients=str(donors)+' donors',States=str(evaluable)+' input-lineage states / 53 frozen',Ecotypes='7 projected',Method='Frozen scRNA recovery / donor coverage; nominal Z threshold; not phenotype test',Evidence=p/'figures/Fig1_source_data.tsv',Result='保留真实Recovery、marker与供者覆盖结果。',TakeHome='程序恢复和供者支持分开评价。',Limit='Pilot及参考来源独立性受限；恢复Z不是细胞身份或外部表型复现。'))

# Step5. bulk Recovery、最新临床和仅基线治疗结局
bulk = root/'06_validation/GSE72819_recovery_review_20261009'
for stem in ['official_ecotype_heatmap','recovery_assignment_summary','state_mean_by_ecotype']:
    registry.append(dict(Source=bulk/'figures'/stem,Module=folders[6],Order=7,Title='GSE72819：'+stem,English='GSE72819 '+stem.replace('_',' '),Main='Main PPT' if stem=='recovery_assignment_summary' else 'Backup',Dataset='GSE72819',N='73 baseline biopsies; 70 outcome-evaluable',Patients='73 baseline / 70 endpoint',States='53',Ecotypes='7',Method='Frozen bulk Recovery; descriptive',Evidence=outcomes/'01_GSE72819_treatment_outcome/recovery_QC.tsv',Result='73份基线输入与70份结局完整病例区分。',TakeHome='Recovery成功不意味着治疗预测成功。',Limit='ETR W10无BH显著信号；该图不能单独作为疗效关联。'))
for stem,title,order,status,claim,evidence in [
 ('Fig1_state_clinical_adjusted','全53状态患者聚类调整临床关联',8,'Backup','完整保留全部状态的阳性与阴性对比。','source_data/Fig1_state_clinical_adjusted_full_statistics.tsv'),
 ('Fig2_ecotype_clinical_adjusted','7生态型患者聚类调整临床关联',8,'Main PPT','E1与局部炎症、组织学、内镜和炎症标志物相关。','source_data/Fig2_ecotype_clinical_adjusted_full_statistics.tsv'),
 ('Fig3_Cell2021_Fig5A_adaptation','最新Ecotype功能与临床全貌',6,'Main PPT','最新功能概览与正式临床审计一致。','source_data/Fig3_dominant_CE_clinical_statistics.tsv'),
 ('Fig4_E1_effects_adjustment','E1调整前后效应及患者敏感性',8,'Main PPT','E1局部炎症beta=1.098 SD，q约3.87e-85；多类临床指标同向。','source_data/representative_effects_before_after_adjustment.tsv'),
 ('Fig5_exact_region_paired','30位患者严格详细部位配对炎症比较',9,'Main PPT','E1炎症减非炎症中位差0.223；E1/E5/E6配对q约0.000887。','statistics/patient_detailed_region_matched_inflammation.tsv')]:
    registry.append(dict(Source=clinical/'figures'/stem,Module=folders[7],Order=order,Title=title,English=stem.replace('_',' '),Main=status,Dataset='MSCCR',N='30 matched pairs' if order==9 else '2025 total; complete cases vary by feature',Patients='30' if order==9 else '919 total; complete cases vary by feature',States='53',Ecotypes='7',Method='Paired Wilcoxon; BH within entity' if order==9 else 'Patient-cluster robust models; BH by entity/weighting; Fig3 module-specific tests',Evidence=clinical/evidence,Result=claim,TakeHome=claim,Limit='横断面评分关联；日期、用药混杂未充分核实；粪钙卫蛋白仅203患者。'))
for image in sorted(outcomes.glob('*/figures/*.png')):
    stem = image.stem
    status = 'Main PPT' if stem in ['Fig1_UC_primary_ecotype_associations','Fig2_progression_and_later_VDZ','Fig4_candidate_patient_distributions'] else 'Backup'
    registry.append(dict(Source=image.with_suffix(''),Module=folders[8],Order=11 if stem in ['Fig2_progression_and_later_VDZ','Fig_S2_RISK_ecotypes'] else 10,Title='基线外部结局：'+stem,English=stem.replace('_',' '),Main=status,Dataset='ETR / IFX / VDZ / RISK; see endpoint table',N='ETR70/12; IFX23/8; VDZ W6 27/6; RISK210/27 (N/events)',Patients='Endpoint patients, one baseline score each; denominators vary',States='53',Ecotypes='7',Method='Firth logistic per SD and Wilcoxon; BH separately per endpoint and 7 CE / 53 state family',Evidence=outcomes/'04_cross_cohort_summary/all_associations.tsv',Result='IFX E3 q=0.01805；VDZ W6 E2 q=0.02129；ETR及VDZ W12/W52阴性；RISK E2调整后q约0.0613。',TakeHome='结局关联依赖药物与终点，需独立验证。',Limit='小样本宽CI；RISK表达尺度溯源受限；仅67份GSE73661 W0重新Recovery；不是已验证预测器。'))

# Step6. 本轮表型图与最新身份审计备份
titles = ['53状态活检水平表型富集','53状态患者及分层敏感性','Discovery与外部scRNA表型富集一致性','7生态型的官方显著分配表型分布']
claims = ['主要炎症和部位关联在患者敏感性后仍有支持；Disease信号较少。','全局打乱活检标签不等于控制重复患者。','名义细胞富集未形成BH显著的整体供者级一致性。','313 Unassigned保留在各表型分母中，避免分布解释偏差。']
for image,title,claim in zip(sorted((out/'figures').glob('P*.png')),titles,claims):
    stem = image.stem
    registry.append(dict(Source=image.with_suffix(''),Module=folders[9],Order=5 if stem.startswith('P4') else 4,Title=title,English=stem.replace('_',' '),Main='Main PPT',Dataset='MSCCR / GSE214695 / IBDverse' if stem.startswith('P3') else 'MSCCR',N='2025 discovery; external comparisons use 12 / 20 donors' if stem.startswith('P3') else '2025 biopsies',Patients='919 discovery; eligible units vary by parent/scope',States='53',Ecotypes='7',Method='Official assignment descriptive percentages' if stem.startswith('P4') else '1000 permutations; empirical BH over1431 prespecified tests/mode; external106 tests/dataset/mode; concordance Fisher BH8',Evidence=out/'statistics'/('ecotype_phenotype_distribution.tsv' if stem.startswith('P4') else ('discovery_external_concordance_all_states.tsv' if stem.startswith('P3') else 'discovery_state_phenotype_enrichment.tsv')),Result=claim,TakeHome=claim,Limit='BayesPrism敏感性；条件分配与样本选择限制；相关性不证明因果；外部稀有父类型覆盖不足。'))
identity = root/'06_validation/cross_lineage_specificity_identity_audit_20261015'
for image in sorted((identity/'figures').glob('*.png')):
    registry.append(dict(Source=image.with_suffix(''),Module=folders[10],Order=12,Title='细胞身份与共享程序审计：'+image.stem,English=image.stem.replace('_',' '),Main='Backup',Dataset='Frozen 53 states and literature signatures',N='Not a new patient cohort',Patients='NA',States='53',Ecotypes='7',Method='Shared-gene / competing lineage evidence; source audit-specific BH families',Evidence=identity/'final_scientific_report.md',Result='最新证据分层A17/B27/C1/D0/E8；没有达到D门槛。',TakeHome='共享炎症程序不能简单当成父谱系错误。',Limit='未执行完整G-like Phase4和同细胞共表达；不宣称53个状态身份全部确证。'))

# Step7. COPY、逐文件SHA和格式/尺寸检查；Excluded只记清单不复制
records = []
quality = []
unique = []
for i,item in enumerate(sorted(registry,key=lambda r:(r['Order'],str(r['Source']))),1):
    figure_id = 'F'+str(i).zfill(3)
    item['FigureID'] = figure_id
    item['Files'] = []
    files = [item['Source'].with_suffix('.'+ext) for ext in ['png','pdf','svg'] if item['Source'].with_suffix('.'+ext).exists()]
    assert files, str(item['Source'])
    for source in files:
        source_hash = hashlib.sha256(source.read_bytes()).hexdigest()
        destination = batch/item['Module']/(figure_id+'_'+source.name)
        dimension = ''
        status = 'EXCLUDED'
        visual = 'Not reviewed: deliberately excluded'
        page_count = 0
        if item['Main']!='Excluded':
            assert not destination.exists()
            shutil.copy2(source,destination)
            destination_hash = hashlib.sha256(destination.read_bytes()).hexdigest()
            assert source_hash==destination_hash
            if source.suffix=='.png':
                with Image.open(destination) as im:
                    im.load()
                    dimension = str(im.width)+' x '+str(im.height)+' px'
                    spread = float(np.std(np.asarray(im.convert('L').resize((160,100)))))
                assert spread>1, str(source)
                status = 'PASS decode/nonblank/SHA'
            elif source.suffix=='.pdf':
                doc = pymupdf.open(destination)
                page_count = len(doc)
                assert page_count>0
                page = doc[0]
                dimension = str(round(page.rect.width,1))+' x '+str(round(page.rect.height,1))+' pt; '+str(page_count)+' pages'
                pix = page.get_pixmap(matrix=pymupdf.Matrix(.3,.3),alpha=False)
                assert np.std(np.frombuffer(pix.samples,dtype=np.uint8))>1
                status = 'PASS PDF render/nonblank/SHA'
                doc.close()
            else:
                import xml.etree.ElementTree as ET
                tree = ET.parse(destination)
                dimension = str(tree.getroot().get('width'))+' x '+str(tree.getroot().get('height'))+'; viewBox '+str(tree.getroot().get('viewBox'))
                status = 'PASS XML parse/SHA; raster preview reviewed separately'
            visual = 'Contact-sheet review required; source QA retained where present'
        else:
            destination_hash = ''
        item['Files'].append(dict(Source=str(source),Destination=str(destination) if item['Main']!='Excluded' else '',Format=source.suffix[1:],SHA256=source_hash))
        records.append(dict(FigureID=figure_id,PPT候选顺序=item['Order'],文件名=source.name,图件标题=item['Title'],所属分析模块=item['Module'],BayesPrism敏感性路线='MSCCR_CORE12_CTSPEC_REGION_IBD; sensitivity only',Original_Source_Path=str(source),Server_Original_Path='/home/mazekai/IBD_EcoTyper/05_EcoTyper/output/MSCCR_CORE12/ctspec_region_ibd/'+str(source.relative_to(frozen)).replace('\\','/') if source.is_relative_to(frozen) else 'See source archive / script provenance',Local_Destination_Path=str(destination) if item['Main']!='Excluded' else '',Source_SHA256=source_hash,Destination_SHA256=destination_hash,图件版本日期=datetime.fromtimestamp(source.stat().st_mtime).isoformat(timespec='seconds'),关联脚本日志统计表=str(item['Evidence']),数据集=item['Dataset'],样本数=item['N'],独立患者或供者数=item['Patients'],Cell_States数量=item['States'],Ecotypes数量=item['Ecotypes'],统计检验与BH家族=item['Method'],文件格式=source.suffix[1:],分辨率或页面尺寸=dimension,科学解释限制=item['Limit'],Disposition=item['Main']))
        quality.append(dict(FigureID=figure_id,File=source.name,Disposition=item['Main'],Integrity=status,Dimensions=dimension,VisualReview=visual,StatisticalProvenance=str(item['Evidence']),InterpretationLimit=item['Limit']))
    unique.append(item)
manifest = pd.DataFrame(records)
manifest.to_csv(nav/'figure_manifest.tsv',sep='\t',index=False,encoding='utf-8-sig')
pd.DataFrame(quality).to_csv(nav/'figure_quality_checklist.tsv',sep='\t',index=False,encoding='utf-8-sig')

# Step8. 每张主图给出中文图注；导航明确主图、备份、阴性结果及缺项
captions = ['# PPT主要图件中文图注','', '全篇统一标注：BayesPrism敏感性路线；不能表述为CIBERSORTx正式主分析。n及患者/供者分母以原图源表为准。','']
index = ['# PPT图件导航','', 'PNG适合直接插入PPT；PDF/SVG保留原有矢量文本。仅COPY，原件和冻结模型未变。Backup仍为真实科研图；Excluded未复制。','']
sections = {1:'Discovery cohort',2:'BayesPrism expression purification',3:'53 Cell states',4:'State phenotype enrichment',5:'7 Ecotypes',6:'Functional characterization',7:'External recovery',8:'MSCCR clinical associations',9:'Patient-matched inflammation',10:'Treatment outcomes',11:'RISK progression',12:'Methods / QC backup'}
for order,title in sections.items():
    index += ['## '+str(order)+'. '+title,'']
    for item in unique:
        if item['Order']!=order:
            continue
        candidates = [f for f in item['Files'] if f['Destination'] and f['Format']=='png']
        candidates += [f for f in item['Files'] if f['Destination'] and f['Format']=='pdf']
        link = str(Path(candidates[0]['Destination']).relative_to(batch)).replace('\\','/') if candidates else ''
        index.append('- '+item['FigureID']+' | '+item['Main']+' | '+('[ '+item['Title']+' ](../'+link+')' if link else item['Title'])+'。'+item['Limit'])
    if order in [1,2]:
        index.append('- 本轮未确认冻结路线的独立现成cohort/纯化科学图；参见cohort_flow汇总表及方法报告。不以示意图冒充分析结果。')
    if order==12:
        index.append('- 未定位冻结53/7的患者级Ecotype成员bootstrap图，NOT VERIFIABLE；比例benchmark的bootstrap不用于替代。')
    index.append('')
for item in unique:
    if item['Main']!='Main PPT':
        continue
    captions += ['## '+item['FigureID']+' | '+item['Title'],'', '**推荐中文标题：** '+item['Title'],'', '**English title:** '+item['English'],'', '**数据与方法：** '+item['Dataset']+'；'+item['N']+'；'+item['Patients']+'。'+item['Method']+'。','', '**重点结果：** '+item['Result'],'', '**Take-home message：** '+item['TakeHome'],'', '**不应过度解释：** '+item['Limit'],'', '**方法及来源：** [Luca et al., Cell 2021](https://doi.org/10.1016/j.cell.2021.09.014)；具体源表：`'+str(item['Evidence'])+'`。','']
(nav/'figure_captions_for_PPT.md').write_text('\n'.join(captions),encoding='utf-8')
(nav/'PPT_figure_index.md').write_text('\n'.join(index),encoding='utf-8')

# Step9. 只复制汇总级源表与说明，禁止把患者ID、矩阵和供者明细放进PPT目录
support = [
 (frozen/'Ecotypes/ecotypes.txt',batch/folders[3]),
 (prior/'00_audit/frozen_states.tsv',batch/folders[3]),
 (out/'Cell2021_phenotype_enrichment_methods_audit.md',nav),
 (out/'statistics/discovery_state_phenotype_enrichment.tsv',batch/folders[9]),
 (out/'statistics/external_scRNA_phenotype_enrichment.tsv',batch/folders[9]),
 (out/'statistics/discovery_external_concordance_all_states.tsv',batch/folders[9]),
 (out/'statistics/concordance_summary.tsv',batch/folders[9]),
 (out/'statistics/ecotype_phenotype_distribution.tsv',batch/folders[9]),
 (out/'statistics/discovery_summary.tsv',batch/folders[9]),
 (out/'audit/discovery_parent_assignment_coverage.tsv',batch/folders[9]),
 (out/'audit/external_phenotype_coverage.tsv',batch/folders[9]),
 (out/'audit/historical_recalculation_summary.json',nav),
 (out/'audit/server_local_semantic_verification.tsv',nav),
 (out/'server_local_frozen_comparison.tsv',nav),
 (clinical/'audit/cohort_flow.tsv',batch/folders[1]),
 (clinical/'source_data/Fig1_state_clinical_adjusted_full_statistics.tsv',batch/folders[7]),
 (clinical/'source_data/Fig2_ecotype_clinical_adjusted_full_statistics.tsv',batch/folders[7]),
 (clinical/'statistics/patient_detailed_region_matched_inflammation.tsv',batch/folders[7]),
 (clinical/'presentation_summary.md',batch/folders[7]),
 (outcomes/'04_cross_cohort_summary/all_associations.tsv',batch/folders[8]),
 (outcomes/'04_cross_cohort_summary/endpoint_denominators.tsv',batch/folders[8]),
 (outcomes/'02_RISK_disease_progression/ecotype_age_sex_sensitivity.tsv',batch/folders[8]),
 (outcomes/'presentation_clinical_outcomes_summary.md',batch/folders[8]),
 (prior/'01_state_recovery/state_recovery_evidence_by_dataset.tsv',batch/folders[6]),
 (prior/'05_ecotype_coassociation/ecotype_coassociation_1000permutations.tsv',batch/folders[6]),
 (prior/'02_external_markers/marker_LOOCV_summary.tsv',batch/folders[4]),
 (prior/'06_ecotype_features/ecotype_continuous_features.tsv',batch/folders[5]),
 (identity/'final_scientific_report.md',batch/folders[10]),
 (root/'06_validation/nonIBD_healthy_recovery_20261010/恢复结果概览.tsv',batch/folders[6]),
 (root/'06_validation/nonIBD_healthy_recovery_20261010/上游独立性与解释修订.md',batch/folders[6])]
support_records = []
for source,destination_folder in support:
    assert source.exists(), str(source)
    destination = destination_folder/source.name
    assert not destination.exists()
    shutil.copy2(source,destination)
    sha = hashlib.sha256(source.read_bytes()).hexdigest()
    assert sha==hashlib.sha256(destination.read_bytes()).hexdigest()
    support_records.append(dict(Source=str(source),Destination=str(destination),SHA256=sha,Kind='Aggregate support / report; no patient rows copied'))
pd.DataFrame(support_records).to_csv(nav/'support_file_manifest.tsv',sep='\t',index=False,encoding='utf-8-sig')

# Step10. 主图和全部图缩略图，只生成派生浏览件，原图不改变
font_path = Path('C:/Windows/Fonts/arial.ttf')
font = ImageFont.truetype(str(font_path),16)
included = [item for item in unique if item['Main']!='Excluded']
for start in range(0,len(included),12):
    sheet = Image.new('RGB',(1800,1400),'white')
    draw = ImageDraw.Draw(sheet)
    for j,item in enumerate(included[start:start+12]):
        png = next((f['Destination'] for f in item['Files'] if f['Format']=='png' and f['Destination']),None)
        if png:
            preview = Image.open(png).convert('RGB')
        else:
            pdf = next(f['Destination'] for f in item['Files'] if f['Format']=='pdf')
            doc = pymupdf.open(pdf)
            pix = doc[0].get_pixmap(matrix=pymupdf.Matrix(1,1),alpha=False)
            preview = Image.frombytes('RGB',[pix.width,pix.height],pix.samples)
            doc.close()
        preview.thumbnail((580,300))
        x,y=(j%3)*600,(j//3)*350
        sheet.paste(preview,(x+(600-preview.width)//2,y+34))
        draw.text((x+10,y+7),item['FigureID']+' '+item['Source'].name[:43],font=font,fill='black')
    sheet.save(nav/('contact_sheet_'+str(start//12+1).zfill(2)+'.png'))
summary = dict(Batch=str(batch),IncludedScientificFigures=len(included),MainFigures=sum(x['Main']=='Main PPT' for x in unique),BackupFigures=sum(x['Main']=='Backup' for x in unique),ExcludedFigures=sum(x['Main']=='Excluded' for x in unique),CopiedFigureFiles=int(manifest.Disposition.ne('Excluded').sum()),CopiedSupportFiles=len(support_records),AllCopyHashesMatch=True,PatientLevelTablesCopied=False)
(nav/'delivery_summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False),encoding='utf-8')
(out/'PPT_batch_path.txt').write_text(str(batch),encoding='utf-8')
print(json.dumps(summary,indent=2,ensure_ascii=False),flush=True)
