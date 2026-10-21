# 本机补充严格活检部位配对、同一病例调整前后对比和临床标签核验。
# 固定使用已保存完整输入，不重新扫描资产或重算已完成的主检验。
from pathlib import Path
import pandas as pd
import numpy as np
from scipy import stats
import statsmodels.api as sm
import json, hashlib, importlib.metadata

# Step1. 严格 RegionDetailed 配对，检验 Region3 合并结肠部位是否影响方向
root = Path('D:/IBD_ecotype/IBD_ecotype_analysis')
out = root/'06_validation/MSCCR_joint_clinical_audit_20261015'
metadata = pd.read_csv(out/'source_data/analysis_metadata_2025.tsv',sep='\t',index_col=0)
matrices = {'state':pd.read_csv(out/'source_data/state_abundance_2025x53.tsv',sep='\t',index_col=0),'ecotype':pd.read_csv(out/'source_data/ecotype_abundance_2025x7.tsv',sep='\t',index_col=0)}
paired = []
for entity,abundance in matrices.items():
    merged = metadata[['PatientID','RegionDetailed','Inflammation']].join(abundance)
    medians = merged.groupby(['PatientID','RegionDetailed','Inflammation'])[abundance.columns].median()
    inflamed = medians.xs('Inflamed',level='Inflammation')
    quiet = medians.xs('Non-inflamed',level='Inflammation')
    both = inflamed.index.intersection(quiet.index)
    delta = (inflamed.loc[both]-quiet.loc[both]).groupby('PatientID').mean()
    delta.to_csv(out/('source_data/paired_detailed_region_deltas_'+entity+'.tsv'),sep='\t')
    for marker in abundance.columns:
        test = stats.wilcoxon(delta[marker],alternative='two-sided',zero_method='wilcox',method='auto')
        paired.append(dict(Entity=entity,Marker=marker,Patients=len(delta),PatientDetailedRegionPairs=len(both),MedianPatientDifference=delta[marker].median(),MeanPatientDifference=delta[marker].mean(),P=test.pvalue,Method='Exact RegionDetailed within-patient median contrast; equal region mean per patient; signed-rank'))
paired = pd.DataFrame(paired)
paired['Q'] = paired.groupby('Entity').P.transform(lambda p:stats.false_discovery_control(p.to_numpy(),method='bh'))
paired.to_csv(out/'statistics/patient_detailed_region_matched_inflammation.tsv',sep='\t',index=False)

# Step2. 同一临床完整病例集的未调整和已调整模型，使用一致标准化和患者聚类
comparisons = []
adjusted = pd.read_csv(out/'statistics/patient_cluster_adjusted_contrasts.tsv',sep='\t')
selections = [('E1','Inflammation'),('E1','Local_GHAS'),('E1','Local_Nancy'),('E1','CRP_log2'),('E1','Fecal_calprotectin_log2'),('E1','SES_CD'),('E1','Mayo_endo')]
for marker,feature in selections:
    rows = adjusted.loc[adjusted.Entity.eq('ecotype') & adjusted.Marker.eq(marker) & adjusted.Feature.eq(feature)]
    primary = rows.loc[rows.Weighting.eq('Biopsy_weighted')].iloc[0]
    disease = primary.Disease
    data = metadata.loc[metadata.Disease.eq(disease) if disease!='All' else metadata.index.notna()].dropna(subset=[feature,'Disease','Region3','Age','Sex'])
    original = matrices['ecotype'].loc[data.index,marker]
    y = (original-original.mean())/original.std(ddof=1)
    x = data[feature].eq('Inflamed').astype(float) if feature=='Inflammation' else (data[feature]-data[feature].mean())/data[feature].std(ddof=1)
    fit = sm.OLS(y,sm.add_constant(x.rename('Exposure'))).fit(cov_type='cluster',cov_kwds={'groups':data.PatientID,'use_correction':True,'df_correction':True},use_t=True)
    ci = fit.conf_int().loc['Exposure']
    comparisons.append(dict(Marker=marker,Feature=feature,Model='Unadjusted_cluster_SE',Beta=fit.params.Exposure,CI_low=ci.iloc[0],CI_high=ci.iloc[1],P=fit.pvalues.Exposure,Q=np.nan,N=len(data),Patients=data.PatientID.nunique()))
    for _,row in rows.iterrows():
        comparisons.append(dict(Marker=marker,Feature=feature,Model='Adjusted' if row.Weighting=='Biopsy_weighted' else 'Adjusted_equal_patient_weight',Beta=row.Beta,CI_low=row.CI_low,CI_high=row.CI_high,P=row.P,Q=row.Q,N=row.N,Patients=row.Patients))
pd.DataFrame(comparisons).to_csv(out/'source_data/representative_effects_before_after_adjustment.tsv',sep='\t',index=False)

# Step3. 缓解标签核验；已有标签的冲突保留而不覆盖
checks = []
for disease,score,threshold in [('CD','SES_CD',3),('UC','Mayo_endo',0)]:
    sub = metadata.loc[metadata.Disease.eq(disease)].dropna(subset=[score,'Endoscopic_remission'])
    expected = sub[score].lt(threshold) if disease=='CD' else sub[score].eq(0)
    for id,row in sub.iterrows():
        checks.append(dict(ID=id,PatientID=row.PatientID,Disease=disease,Field='Endoscopic_remission',Score=score,Value=row[score],SourceLabel=row.Endoscopic_remission,Expected=int(expected.loc[id]),Agrees=bool(row.Endoscopic_remission==expected.loc[id])))
for disease,score in [('CD','Max_GHAS'),('UC','Max_Nancy')]:
    sub = metadata.loc[metadata.Disease.eq(disease)].dropna(subset=[score,'Histological_remission'])
    for id,row in sub.iterrows():
        checks.append(dict(ID=id,PatientID=row.PatientID,Disease=disease,Field='Histological_remission',Score=score,Value=row[score],SourceLabel=row.Histological_remission,Expected=int(row[score]==0),Agrees=bool(row.Histological_remission==(row[score]==0))))
checks = pd.DataFrame(checks)
checks.to_csv(out/'audit/remission_label_reconciliation.tsv',sep='\t',index=False)

# Step4. 保存环境、输入校验和临床缺失概况，避免将当前关联当治疗结局
versions = {p:importlib.metadata.version(p) for p in ['numpy','pandas','scipy','statsmodels','matplotlib','openpyxl','pymupdf']}
(out/'audit/package_versions.json').write_text(json.dumps(versions,indent=2),encoding='utf-8')
inputs = [root/'ctspec_region_ibd/Ecotypes/ecotypes.txt', root/'ctspec_region_ibd/Ecotypes/ecotype_abundance.txt', Path('D:/IBD_ecotype/Profile/audit/2026-10-06_ctspec_region_ibd/GSE193677_series_matrix.txt.gz')]
hashes = [dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in inputs]
pd.DataFrame(hashes).to_csv(out/'audit/key_input_sha256.tsv',sep='\t',index=False)
summary = dict(ExactRegionMatchedPatients=int(paired.Patients.max()),ExactRegionPairs=int(paired.PatientDetailedRegionPairs.max()),RemissionCompared=len(checks),RemissionDiscrepancies=int((~checks.Agrees).sum()),PrimaryContrastSignificant=adjusted.loc[adjusted.Weighting.eq('Biopsy_weighted')].groupby('Entity').Q.agg(lambda x:int(x.lt(.05).sum())).to_dict())
(out/'audit/phase4_sensitivity_checkpoint.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
print(json.dumps(summary,indent=2),flush=True)
print(paired.loc[paired.Entity.eq('ecotype')].to_string(index=False),flush=True)
