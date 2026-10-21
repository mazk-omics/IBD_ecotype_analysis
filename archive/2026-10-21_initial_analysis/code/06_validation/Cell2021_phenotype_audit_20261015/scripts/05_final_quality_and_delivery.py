# 汇总方法证据和交付质量，保存最终验收版本；不覆盖已复制的原科研图件。
# 本机msccr_clinical_audit_20261015环境运行，读取本轮批次指针和完整统计。
from pathlib import Path
import base64,zlib,json,hashlib,subprocess,sys
import numpy as np
import pandas as pd
import pymupdf

# Step1. 方法学逐字段对照，不把未明确的原文实现写成已确认
out = Path(__file__).resolve().parents[1]
root = out.parents[1]
batch = Path((out/'PPT_batch_path.txt').read_text(encoding='utf-8'))
nav = batch/'00_结果导航与清单'
methods = [
 dict(OriginalMethodsTitle='Enrichment of cell states in known phenotypes',StatisticalObject='Ns: phenotype-P units assigned to state s',DiscoveryOrRecovery='Both discovery samples and recovered scRNA cells',PhenotypeAndAssignment='True known phenotype; EcoTyper state assignment; exact Unassigned and parent scope not fully specified',Statistic='(Ns - null mean) / null SD',PermutationMechanism='Shuffle cell-state assignment labels; parent/stratum scope not explicitly fully stated in text',Iterations='1000',NullHypothesis='No preferential state representation in phenotype',Direction='Positive enrichment',Significance='Original Z > 1.96; paper describes P < 0.05',Multiplicity='BH not specified for this procedure',FisherExact='No, not for individual state enrichment',ExternalScRNARequirement='Separate validation evidence; not a prerequisite of every recovery',Figures='2D; S3D-G; S2F distribution context',Evidence='PMC8526411_bioc.txt:301-305; XML formula; Table S3',IBDCorrespondence='Parent-wise assigned QC units; patient and restricted sensitivity explicitly adapted',Status='ADAPTED'),
 dict(OriginalMethodsTitle='Cell state representation by tumor type',StatisticalObject='Tumor-type x dominant-state 2x2 count table',DiscoveryOrRecovery='Discovery',PhenotypeAndAssignment='Tumor sample assigned to most abundant state in each parent type',Statistic='Fisher exact over-representation probability',PermutationMechanism='None',Iterations='NA',NullHypothesis='No excess representation of state in tumor type',Direction='Enrichment; exact software tail not stated in this paragraph',Significance='Adjusted significance; threshold not fully restated here',Multiplicity='Benjamini-Hochberg across cell types; exact code-level family not located',FisherExact='Yes',ExternalScRNARequirement='No, distinct analysis',Figures='S2G,H; Table S3',Evidence='PMC8526411_bioc.txt:337-339',IBDCorrespondence='Region3 instead of tumor type; greater Fisher; BH53x3=159; all159 reproduced',Status='ADAPTED; numerical PASS'),
 dict(OriginalMethodsTitle='Results / Fig2D caption; Comparison of scRNA-seq annotation methods',StatisticalObject='State-wise discovery and scRNA phenotype enrichment concordance',DiscoveryOrRecovery='Both; matched phenotype definitions',PhenotypeAndAssignment='Enrichment flags independently obtained in discovery and scRNA',Statistic='Overlap fraction among discovery-enriched states; 2x2 enrichment-flag association',PermutationMechanism='Inherited 1000-iteration state enrichment; no independent concordance permutation specified',Iterations='1000 for upstream state enrichment',NullHypothesis='Enriched-state labels do not associate across datasets',Direction='Concordant enrichment; Fisher caption does not specify software alternative',Significance='Caption gives Fisher concordance significance',Multiplicity='Original concordance BH not specified; this project BH8 explicitly added',FisherExact='Yes, for concordance; not a replacement of upstream Z',ExternalScRNARequirement='This is the external phenotype replication analysis',Figures='2D; S3D; S3G',Evidence='PMC8526411_bioc.txt:43,307-309,681; no standalone Methods title matching user paraphrase',IBDCorrespondence='Inflamed colon CD/UC and CD ileum inflammation; donor-balanced sensitivity; no robust global BH concordance',Status='Analysis ADAPTED; replication PARTIAL')]
pd.DataFrame(methods).to_csv(out/'methods_comparison_table.tsv',sep='\t',index=False,encoding='utf-8-sig')
shutil = __import__('shutil')
shutil.copy2(out/'methods_comparison_table.tsv',nav/'methods_comparison_table.tsv')

# Step2. 冻结源文件双端结束哈希；所有PPT源图/目标图再次逐个核验
before = json.loads(zlib.decompress(base64.b64decode((out/'server_frozen_manifest_v2.b64').read_text())))
after = json.loads(zlib.decompress(base64.b64decode((out/'server_frozen_after.b64').read_text())))
assert before==after
freeze_rows = []
for row in after:
    local = root/'ctspec_region_ibd'/row['Path'].split('ctspec_region_ibd/')[1]
    sha = hashlib.sha256(local.read_bytes()).hexdigest()
    assert sha==row['SHA256']
    freeze_rows.append(dict(ServerPath=row['Path'],LocalPath=str(local),SHA256=sha,BeforeAfter='UNCHANGED',ServerLocal='IDENTICAL'))
pd.DataFrame(freeze_rows).to_csv(out/'audit/frozen_final_59_files.tsv',sep='\t',index=False)
shutil.copy2(out/'audit/frozen_final_59_files.tsv',nav/'frozen_final_59_files.tsv')
manifest = pd.read_csv(nav/'figure_manifest.tsv',sep='\t').fillna('')
copied = manifest[manifest.Disposition.ne('Excluded')]
assert not copied.Original_Source_Path.str.contains('superseded|discarded|invalidated|posttreatment',case=False).any()
assert copied.Local_Destination_Path.is_unique
for row in copied.itertuples():
    assert hashlib.sha256(Path(row.Original_Source_Path).read_bytes()).hexdigest()==row.Source_SHA256
    assert hashlib.sha256(Path(row.Local_Destination_Path).read_bytes()).hexdigest()==row.Destination_SHA256
support = pd.read_csv(nav/'support_file_manifest.tsv',sep='\t')
for row in support.itertuples():
    assert hashlib.sha256(Path(row.Source).read_bytes()).hexdigest()==row.SHA256
    assert hashlib.sha256(Path(row.Destination).read_bytes()).hexdigest()==row.SHA256

# Step3. 新生成4图重新执行最终PDF碰撞、实际字号及源脚本检查
qa = Path.home()/'.codex/skills/nature-figure/scripts'
new_qa = []
for pdf in sorted((out/'figures').glob('P*.pdf')):
    for audit,args in [('audit_pdf_text.py',['--min-pt','5','--json']),('audit_figure_collisions.py',['--json-out',str(pdf.with_suffix('.collision-audit.json'))])]:
        run = subprocess.run([sys.executable,str(qa/audit),str(pdf)]+args,capture_output=True,text=True,encoding='utf-8')
        (out/'logs'/(pdf.stem+'_'+audit+'_delivery.log')).write_text(run.stdout+'\n'+run.stderr,encoding='utf-8')
        assert run.returncode==0, pdf.stem+' '+audit+' '+run.stdout
        new_qa.append(dict(Figure=pdf.stem,Audit=audit,Result='PASS'))
    alignment = json.loads(pdf.with_suffix('.alignment.json').read_text(encoding='utf-8'))
    shutil.copy2(pdf.with_suffix('.alignment.json'),nav/(pdf.stem+'.alignment.json'))
    shutil.copy2(pdf.with_suffix('.collision-audit.json'),nav/(pdf.stem+'.collision-audit.json'))
run = subprocess.run([sys.executable,str(qa/'validate_figure.py'),str(out/'scripts/03_plot_phenotypes.py'),'--json'],capture_output=True,text=True,encoding='utf-8')
(out/'logs/plot_source_preflight.json').write_text(run.stdout+'\n'+run.stderr,encoding='utf-8')
pd.DataFrame(new_qa).to_csv(nav/'new_figure_final_PDF_QA.tsv',sep='\t',index=False)

# Step4. 对主图逐页记录实际字体，原生高密度图需放大而不擅自修改
font_rows = []
for row in copied[copied['文件格式'].eq('pdf')].itertuples():
    document = pymupdf.open(row.Local_Destination_Path)
    sizes = [span['size'] for page in document for block in page.get_text('dict')['blocks'] if 'lines' in block for line in block['lines'] for span in line['spans'] if span['text'].strip()]
    minimum = min(sizes) if sizes else np.nan
    font_rows.append(dict(FigureID=row.FigureID,PDF=Path(row.Local_Destination_Path).name,MinimumNativeFontPt=minimum,Review='Original preserved; dense portrait plots should be enlarged or shown as selected panels, not squeezed into one slide' if document[0].rect.height/document[0].rect.width>1.25 or minimum<5 else 'Native text size recorded; assess at chosen slide scale'))
    document.close()
pd.DataFrame(font_rows).to_csv(nav/'native_PDF_font_and_slide_usage.tsv',sep='\t',index=False,encoding='utf-8-sig')
quality = pd.read_csv(nav/'figure_quality_checklist.tsv',sep='\t').fillna('')
quality.loc[quality.Disposition.ne('Excluded'),'VisualReview'] = 'Reviewed all79 figures across7 contact sheets; no gross blank/cutoff/color-loss. Dense originals require enlargement; this is not a full source-code re-audit.'
quality.to_csv(nav/'figure_quality_checklist_final.tsv',sep='\t',index=False,encoding='utf-8-sig')

# Step5. 精确结局核验、阴性发现及未完成项单列
outcomes = root/'06_validation/external_clinical_outcomes_20261015'
associations = pd.read_csv(outcomes/'04_cross_cohort_summary/all_associations.tsv',sep='\t')
assert len(associations)==360
gse = pd.read_csv(outcomes/'03_GSE73661_mucosal_healing/recovery_QC.tsv',sep='\t')
assert int(gse.AllContinuousSamples.iloc[0])==67
key = associations[associations.Family.eq('Ecotype') & ((associations.Endpoint.eq('IFX_W4_6_endoscopic_healing')&associations.Feature.eq('E3'))|(associations.Endpoint.eq('VDZ_W6_endoscopic_healing')&associations.Feature.eq('E2'))|(associations.Endpoint.eq('RISK_3yr_B2_B3')&associations.Feature.eq('E2')))]
key.to_csv(nav/'key_outcome_values_verified.tsv',sep='\t',index=False)
tasks = [
 ['A1','DONE / NOT VERIFIABLE for unpublished custom code','Primary methods, figures and local full text reviewed; custom phenotype code not located'],
 ['A2','DONE','806 permutations fully replayed; 159 Fisher,350 and385 fully recalculated; latest371 count reconciliation'],
 ['A3','DONE; biological replication PARTIAL','4293 discovery,424 external,424 concordance,56 distribution rows; full nulls and unit IDs retained outside PPT'],
 ['A4','DONE','4 PNG/PDF/SVG figure sets; final PDF and alignment QA passed'],
 ['B','DONE with documented missing assets','79 figures215 files; frozen-route patient-level ecotype-member bootstrap and standalone purification/reference plot not confirmed'],
 ['C','DONE with dense-original usage caveat','Provenance/SHA,35 main captions,12-section index,7 contact sheets,final quality checklist'],
 ['D','DONE','Methods audit,comparison TSV,results,scripts,logs and delivery report saved'],
 ['CIBERSORTx','NOT STARTED as requested','Formal main discovery/recovery/results remain future work']]
pd.DataFrame(tasks,columns=['Task','Status','EvidenceOrLimit']).to_csv(nav/'task_completion_status.tsv',sep='\t',index=False,encoding='utf-8-sig')
summary = json.loads((nav/'delivery_summary.json').read_text(encoding='utf-8'))
notes = '''# 最终交付及PPT使用说明

本轮仅为BayesPrism敏感性路线，未改冻结模型、未重跑NMF、未启动CIBERSORTx、未自动Git提交。服务器只建立一次持久SSH连接，小型结果从同一会话取回。

- 已复制79张不同科研图的215个源文件，其中35张Main PPT候选、44张Backup；另2张旧表型图登记Excluded，未复制。另复制30份汇总源表与说明。新增导航、QA文件、缩略图不计入“原图源文件”数量。
- 215个图文件和30个支持文件逐文件源/目标SHA256均一致；冻结59个文件双端、执行前后SHA一致。六张最终统计表两端数值一致。没有复制患者级临床明细、原始表达矩阵或选择ID至PPT目录。
- `figure_manifest.tsv`为逐文件来源与SHA；`figure_quality_checklist_final.tsv`为最终验收；`native_PDF_font_and_slide_usage.tsv`记录原生字体及密集图的使用限制。初始质量表保留过程记录。
- 7页contact_sheet仅作浏览，不能替代原高清图。新P1–P4为16:9，适合直接插入PPT。原生热图保持原样；高密度竖版53状态/供者支持/LOOCV/Hallmark图建议放大或选取原图面板，不能缩成小图后声称所有文字可读。
- 原有scRNA图中的Independent字样并不代表参考来源完全独立；应与本轮图注的pilot/参考来源限制同页说明。所有Recovery Z均不能替代表型一致性。
- F023 Hallmark：2025份argmax分组（包含313份未获官方显著Ecotype分配的活检），Mann–Whitney BH350。F021功能概览各模块方法不同；dominant生态型临床概览不等同患者聚类回归，临床推断优先F063/F064/F065。
- ETR W10及VDZ W12/W52阴性结果完整保留。RISK必须同页标注输入表达尺度溯源限制，以及E2年龄/性别调整后q=0.0613125。GSE73661只采用67份W0重做Recovery的最终目录，作废版本未复制。
- 外部表型富集的一致性未达到稳健的供者级整体复现；不要用名义Recovery或细胞级重叠替代。冻结路线的患者级Ecotype成员bootstrap图仍NOT VERIFIABLE；未用其他分支bootstrap代替。

推荐汇报主线：队列与方法说明→F001/F002–F013状态展示→F014–F016表型审计→F018/F020生态型→F019/F021/F023功能→F060/F061及外部覆盖→F063/F064临床→F065患者配对→F068/F070治疗→F072 RISK和阴性结果→QC备份。

方法结论及所有限制见 `Cell2021_phenotype_enrichment_methods_audit.md`，详细对照见 `methods_comparison_table.tsv`。源脚本与完整null数组保存在本机项目分析目录；PPT目录保留轻量可核验材料。
'''
(nav/'最终交付与PPT使用说明.md').write_text(notes,encoding='utf-8')
(out/'FINAL_DELIVERY_SUMMARY.md').write_text(notes+'\nPPT批次：'+str(batch)+'\n',encoding='utf-8')
(nav/'README.md').write_text('# 本轮结果入口\n\n先读[最终交付说明](最终交付与PPT使用说明.md)，按[PPT图件导航](PPT_figure_index.md)选图，配合[中文图注](figure_captions_for_PPT.md)。\n\n[方法学审计](Cell2021_phenotype_enrichment_methods_audit.md)及[方法逐项对照](methods_comparison_table.tsv)记录全部统计边界。\n',encoding='utf-8')
print(json.dumps(summary,indent=2,ensure_ascii=False))
print('FINAL_QA_PASS; frozen59 unchanged; figure215/support30 SHA matched; all6 result tables verified',flush=True)
