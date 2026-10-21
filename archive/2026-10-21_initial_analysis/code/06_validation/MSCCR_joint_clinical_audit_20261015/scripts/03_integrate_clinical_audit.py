# 本机整合服务器只读快照、GEO 临床字段和冻结结果，保存全部差异与缺失。
# 运行环境 msccr_clinical_audit_20261015；只向本次输出目录写入。
from pathlib import Path
import pandas as pd
import numpy as np
import json, io, csv, hashlib

# Step1. 保存服务器资产清单，并按路径及校验值确认本机副本
root = Path('D:/IBD_ecotype/IBD_ecotype_analysis')
out = root/'06_validation/MSCCR_joint_clinical_audit_20261015'
remote = json.loads((out/'audit/remote_payload.json').read_text())
inventory = pd.DataFrame(remote['inventory'])
local = pd.read_csv(out/'audit/local_inventory.tsv',sep='\t')
local_hash = local.dropna(subset=['sha256']).groupby('sha256').path.agg(list).to_dict()
inventory['local_hash_copies'] = inventory.sha256.map(local_hash).map(lambda v: ';'.join(v) if isinstance(v,list) else '')
inventory['phase'] = np.select([inventory.path.str.contains('01_raw_data'),inventory.path.str.contains('expression_matrices'),inventory.path.str.contains('/ctspec_region_ibd/'),inventory.path.str.contains('/Cell_States/|/Cell_type_specific_genes/|/Ecotypes/')],['Raw','G_like_expression','Frozen_output','NMF_or_model_intermediate'],default='Derived_or_candidate')
inventory['classification'] = inventory.phase.map({'Raw':'raw','Frozen_output':'frozen','G_like_expression':'derived','NMF_or_model_intermediate':'derived','Derived_or_candidate':'derived_or_other_cohort'})
inventory['mtime_utc'] = pd.to_datetime(inventory.modified,unit='s',utc=True)
inventory.to_csv(out/'audit/server_inventory.tsv',sep='\t',index=False)
verified = []
for file in (root/'ctspec_region_ibd').rglob('*'):
    if file.name not in ['state_abundances.txt','state_assignment.txt','ecotype_abundance.txt','ecotypes.txt','ecotype_assignment.txt']:
        continue
    relative = file.relative_to(root/'ctspec_region_ibd').as_posix()
    server_path = '/home/mazekai/IBD_EcoTyper/05_EcoTyper/output/MSCCR_CORE12/ctspec_region_ibd/'+relative
    records = inventory.loc[inventory.path.eq(server_path)]
    digest = hashlib.sha256(file.read_bytes()).hexdigest()
    verified.append(dict(local_path=str(file),server_path=server_path,local_sha256=digest,server_sha256=records.sha256.iloc[0] if len(records) else '',identical=bool(len(records) and records.sha256.iloc[0]==digest)))
pd.DataFrame(verified).to_csv(out/'audit/frozen_copy_hash_comparison.tsv',sep='\t',index=False)
final_verification = out/'audit/remote_final_verification.json'
if final_verification.exists():
    final_hashes = json.loads(final_verification.read_text())['hashes']
    for row in verified:
        row['server_sha256'] = final_hashes[row['server_path']]
        row['identical'] = row['local_sha256']==row['server_sha256']
    assert all(row['identical'] for row in verified)
    pd.DataFrame(verified).to_csv(out/'audit/frozen_copy_hash_comparison.tsv',sep='\t',index=False)
for path,text in remote['texts'].items():
    if not path.endswith('soft.gz'):
        (out/'source_data'/('server_'+Path(path).name)).write_text(text,encoding='utf-8')

# Step2. 服务器 SOFT 原始字段与本地 GEO series matrix 用 GSM 逐值核对
soft = next(v for k,v in remote['texts'].items() if k.endswith('soft.gz'))
records = []
for line in soft.splitlines():
    if line.startswith('^SAMPLE = '):
        records.append({'geo_accession':line.split(' = ',1)[1]})
    elif line.startswith('!Sample_title = '):
        records[-1]['title'] = line.split(' = ',1)[1]
    elif line.startswith('!Sample_characteristics_ch1 = '):
        label,value = line.split(' = ',1)[1].split(': ',1)
        records[-1][label] = value
server_geo = pd.DataFrame(records).replace({'NA':np.nan,'':np.nan})
server_geo.to_csv(out/'source_data/server_GEO_SOFT_metadata_2490.tsv',sep='\t',index=False)
geo = pd.read_csv(out/'source_data/GEO_metadata_2490.tsv',sep='\t')
geo = geo.set_index('geo_accession')
server_geo = server_geo.set_index('geo_accession').reindex(geo.index)
differences = []
for col in server_geo.columns.intersection(geo.columns):
    a = geo[col].fillna('NA').astype(str)
    b = server_geo[col].fillna('NA').astype(str)
    numeric_a = pd.to_numeric(geo[col],errors='coerce')
    numeric_b = pd.to_numeric(server_geo[col],errors='coerce')
    equal = a.eq(b) | np.isclose(numeric_a,numeric_b,equal_nan=False,atol=1e-12)
    for gsm in a.index[~equal]:
        differences.append(dict(GSM=gsm,Field=col,Local=a.loc[gsm],Server=b.loc[gsm]))
pd.DataFrame(differences,columns=['GSM','Field','Local','Server']).to_csv(out/'audit/local_GEO_vs_server_SOFT_differences.tsv',sep='\t',index=False)
assert not differences

# Step3. 样本集合和患者对应逐层核对；原始 counts 的零样本与下游技术零分开
data = pd.read_csv(out/'source_data/Discovery_GEO_full_mapping.tsv',sep='\t').set_index('ID')
assert data.PatientID.eq(data.GEO_PatientID).all()
assert data.Disease.eq(data.ibd_disease).all()
assert data.Inflammation.eq(data.typere.map({'I':'Inflamed','NonI':'Non-inflamed'})).all()
assert data.RegionDetailed.eq(data.regionre).all()
full = pd.read_csv(out/'source_data/server_MSCCR_CORE12_annotation.tsv',sep='\t')
metadata2486 = pd.read_csv(out/'source_data/server_discovery_metadata_2486.tsv',sep='\t')
excluded = geo.loc[~geo.SampleID.isin(full.ID)].copy()
excluded.to_csv(out/'audit/excluded_4_GEO_records.tsv',sep='\t')
assert len(full)==2486 and len(excluded)==4 and set(full.ID)==set(metadata2486.sample_id)
flow = []
for label,frame,idcol,patientcol in [('GEO_raw',geo,'SampleID','GEO_PatientID'),('G_like_eligible',full,'ID','PatientID'),('Frozen_IBD',data.reset_index(),'ID','PatientID')]:
    flow.append(dict(Stage=label,Biopsies=len(frame),Patients=frame[patientcol].nunique(),CD=int(frame['ibd_disease' if label=='GEO_raw' else 'Disease'].eq('CD').sum()),UC=int(frame['ibd_disease' if label=='GEO_raw' else 'Disease'].eq('UC').sum())))
pd.DataFrame(flow).to_csv(out/'audit/cohort_flow.tsv',sep='\t',index=False)
assigned = pd.read_csv(root/'ctspec_region_ibd/Ecotypes/ecotype_assignment.txt',sep='\t')
data['AssignedEcotype'] = assigned.set_index('ID').Ecotype.reindex(data.index).fillna('Unassigned')
ecotypes = pd.read_csv(out/'source_data/ecotype_abundance_2025x7.tsv',sep='\t',index_col=0).reindex(data.index)
tie = ecotypes.eq(ecotypes.max(axis=1),axis=0).sum(axis=1).gt(1)
data['DominantEcotype'] = ecotypes.idxmax(axis=1).mask(tie)
assert data.AssignedEcotype.ne('Unassigned').sum()==1712 and data.AssignedEcotype.eq('Unassigned').sum()==313
pd.crosstab(data.AssignedEcotype,data.DominantEcotype,dropna=False).to_csv(out/'audit/assigned_vs_dominant.tsv',sep='\t')

# Step4. 临床字典采用 GEO 字段原名和原量纲；单位存在歧义时保留原始 log2 数值
spec = [
 ('study_eligibility_age_at_endo','Age','continuous','All','patient','Age at study endoscopy, years'),
 ('demographics_gender','Sex','categorical','All','patient','GEO Male/Female'),
 ('log2_fecalcalpro_mgperg','Fecal_calprotectin_log2','continuous','All','patient','GEO log2_fecalcalpro_mgperg; native log2 scale; do not back-transform without offset/unit clarification'),
 ('crp_jjmgl_log2','CRP_log2','continuous','All','patient','GEO crp_jjmgl_log2; CRP mg/L label confirmed by supplement CRP_JJmgL_log2; native log2; offset not documented'),
 ('ibd_clinicianmeasure_inactive_active','Clinical_activity','categorical','All','patient','Clinician measure Active/Inactive; cross-sectional'),
 ('ibd_endoseverity_4levels','Endoscopic_severity','categorical','All','patient','Inactive/Mild/Moderate/Severe; disease-specific source thresholds not reconstructed'),
 ('ghas_sum7','Local_GHAS','continuous','CD','biopsy','Local Global Histological Activity Score, source ghas_sum7'),
 ('nancyindex','Local_Nancy','continuous','UC','biopsy','Local Nancy histological index'),
 ('ibdsescd_totalsescd','SES_CD','continuous','CD','patient','Total Simple Endoscopic Score for Crohn disease'),
 ('ibdmesuc_mayo_score','Mayo_endo','continuous','UC','patient','Endoscopic Mayo score'),
 ('harveybradshawindex_hbi_score','HBI','continuous','CD','patient','Harvey-Bradshaw Index'),
 ('colitisactivityindex_sccai','SCCAI','continuous','UC','patient','Simple Clinical Colitis Activity Index'),
 ('max_ghas_sum7','Max_GHAS','continuous','CD','patient','Maximum local GHAS over available biopsies'),
 ('max_nancy','Max_Nancy','continuous','UC','patient','Maximum local Nancy over available biopsies'),
 ('endoremiss','Endoscopic_remission','categorical','All','patient','1 remission, 0 no remission; paper SES-CD<3 in CD / Mayo=0 in UC'),
 ('historemiss','Histological_remission','categorical','All','patient','1 remission, 0 no remission; paper max GHAS=0 in CD / max Nancy=0 in UC')]
dictionary = []
coverage = []
for field,label,kind,disease,level,definition in spec:
    data[label] = pd.to_numeric(data[field],errors='coerce') if kind=='continuous' or field in ['endoremiss','historemiss'] else data[field]
    eligible = data.Disease.eq(disease) if disease!='All' else pd.Series(True,index=data.index)
    sub = data.loc[eligible]
    varying = sub.groupby('PatientID')[label].nunique().gt(1).sum()
    dictionary.append(dict(Field=field,Feature=label,Kind=kind,AnalysisDisease=disease,Level=level,Definition=definition,Source='GEO series matrix; server SOFT exact value agreement',SameEncounter='Age is at endoscopy; other clinical fields cohort-linked, exact measurement interval/date unavailable',ObservedMin=sub[label].min() if kind=='continuous' else '',ObservedMax=sub[label].max() if kind=='continuous' else '',NonmissingBiopsies=int(sub[label].notna().sum()),NonmissingPatients=sub.loc[sub[label].notna(),'PatientID'].nunique(),PatientsWithVaryingValues=int(varying),Use='cross-sectional association'))
    for (dis,region),group in data.groupby(['Disease','Region3']):
        observed = group[group[label].notna()]
        coverage.append(dict(Field=field,Feature=label,Disease=dis,Region3=region,NonmissingBiopsies=len(observed),NonmissingPatients=observed.PatientID.nunique(),UsedForDiseaseSpecificAnalysis=disease in ['All',dis]))
data['Disease_feature'] = data.Disease
data['Inflammation_feature'] = data.Inflammation
data['Region3_feature'] = data.Region3
for field,level,definition in [('ID','biopsy','Unique biopsy ID extracted from original GEO title'),('PatientID','patient','reGRID patient identifier extracted from original title; 919 unique patients'),('Disease','patient','CD or UC; original ibd_disease'),('Inflammation','biopsy','Local endoscopically defined inflamed/non-inflamed tissue; GEO typere I/NonI'),('Region3','biopsy','Derived Ileum/Colon/Rectum category; colon subsites merged'),('RegionDetailed','biopsy','Original GEO regionre anatomical category'),('RegionBroad','biopsy','Ileum or Colon_rectum; original established mapping')]:
    dictionary.append(dict(Field=field,Feature=field,Kind='identifier' if field in ['ID','PatientID'] else 'categorical',AnalysisDisease='All',Level=level,Definition=definition,Source='Frozen annotation reconciled with original GEO title/characteristics',SameEncounter='Biopsy-linked identity/anatomy; sampling date unavailable',NonmissingBiopsies=2025,NonmissingPatients=919,Use='ID mapping only' if field in ['ID','PatientID'] else 'cross-sectional association or anatomy annotation'))
pd.DataFrame(dictionary).to_csv(out/'source_data/clinical_variable_dictionary.tsv',sep='\t',index=False)
pd.DataFrame(coverage).to_csv(out/'source_data/clinical_coverage_by_disease_region.tsv',sep='\t',index=False)
data.to_csv(out/'source_data/analysis_metadata_2025.tsv',sep='\t')
missing_requested = ['bMIS precomputed per biopsy','Sampling date','Disease duration','Montreal behavior B1/B2/B3','Disease extent','Stricture','Penetration','Deep ulcer','Prior biologics','Medication exposure','Prospective treatment outcome','Follow-up','Albumin','Hemoglobin','UCEIS']
pd.DataFrame({'RequestedFeature':missing_requested,'Status':'No sample-linked field located in audited MSCCR sources; supplement cohort aggregates/gene lists are not per-sample measurements'}).to_csv(out/'audit/requested_unavailable_features.tsv',sep='\t',index=False)
summary = dict(LocalAssets=len(local),ServerAssets=len(inventory),LocalSOFTFieldDifferences=len(differences),Biopsies=2025,Patients=919,Assigned=1712,Unassigned=313,DominantTies=int(tie.sum()),CopyHashMatched=int(pd.DataFrame(verified).identical.sum()),CopyHashComparisons=len(verified),Excluded4=excluded.SampleID.tolist(),RawCountsZeroSamples=0)
(out/'audit/phase2_3_checkpoint.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
(out/'progress_summary.md').write_text('# 2026-10-15 进度\n\nPhase 1–3 完成：资产清单、SOFT/GEO逐值核对、冻结维度、临床字典及样本映射已保存。统计尚未运行。\n',encoding='utf-8')
print(json.dumps(summary,indent=2),flush=True)
print(pd.DataFrame(dictionary)[['Feature','NonmissingBiopsies','NonmissingPatients']].to_string(index=False),flush=True)
