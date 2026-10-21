# 本机按已保存 statistical_analysis_plan.md 分析冻结 abundance 与实际临床字段。
# conda 环境 msccr_clinical_audit_20261015；输出完整检验表和阶段检查点。
from pathlib import Path
import pandas as pd
import numpy as np
from scipy import stats
import statsmodels.api as sm
import itertools, json

# Step1. 读取并验证唯一样本映射，固定检验家族和分类参照
root = Path('D:/IBD_ecotype/IBD_ecotype_analysis')
out = root/'06_validation/MSCCR_joint_clinical_audit_20261015'
metadata = pd.read_csv(out/'source_data/analysis_metadata_2025.tsv',sep='\t',index_col=0)
dictionary = pd.read_csv(out/'source_data/clinical_variable_dictionary.tsv',sep='\t')
continuous = dictionary.loc[dictionary.Kind.eq('continuous'),['Feature','AnalysisDisease']].values.tolist()
categorical = [('Disease','UC'),('Inflammation','Non-inflamed'),('Region3','Colon'),('Sex','Female'),('Clinical_activity','Inactive'),('Endoscopic_severity','Inactive'),('Endoscopic_remission',0.0),('Derived_histology_zero',0.0)]
matrices = {'state':pd.read_csv(out/'source_data/state_abundance_2025x53.tsv',sep='\t',index_col=0), 'ecotype':pd.read_csv(out/'source_data/ecotype_abundance_2025x7.tsv',sep='\t',index_col=0)}
assert metadata.index.is_unique and metadata.PatientID.nunique()==919
for frame in matrices.values():
    assert frame.index.equals(metadata.index)
members = pd.read_csv(out/'source_data/frozen_state_membership.tsv',sep='\t')

# Step2. 分类变量组间比较及全部两两对比，效应量保留 abundance 原单位
categorical_results = []
pairwise = []
for entity,abundance in matrices.items():
    for feature,reference in categorical:
        valid = metadata[feature].notna()
        values = metadata.loc[valid,feature]
        levels = sorted(values.unique())
        for marker in abundance.columns:
            groups = [abundance.loc[values.index[values.eq(level)],marker].values for level in levels]
            if len(levels)==2:
                a = abundance.loc[values.index[values.ne(reference)],marker]
                b = abundance.loc[values.index[values.eq(reference)],marker]
                test = stats.mannwhitneyu(a,b,alternative='two-sided',method='asymptotic')
                effect = a.mean()-b.mean()
                size = 2*test.statistic/(len(a)*len(b))-1
                method = 'Mann-Whitney two-sided'
                contrast = str(next(v for v in levels if v!=reference))+' vs '+str(reference)
            else:
                test = stats.kruskal(*groups)
                effect = max(np.mean(g) for g in groups)-min(np.mean(g) for g in groups)
                size = max(0,(test.statistic-len(levels)+1)/(valid.sum()-len(levels)))
                method = 'Kruskal-Wallis'
                contrast = 'Overall; range of group mean abundance'
            categorical_results.append(dict(Entity=entity,Marker=marker,Feature=feature,Method=method,Contrast=contrast,Effect=effect,RankBiserial_or_EpsilonSquared=size,Statistic=test.statistic,P=test.pvalue,N=int(valid.sum()),Patients=metadata.loc[valid,'PatientID'].nunique(),GroupCounts=json.dumps(values.value_counts().to_dict())))
            if len(levels)>2:
                for high,low in itertools.combinations(levels,2):
                    a = abundance.loc[values.index[values.eq(high)],marker]
                    b = abundance.loc[values.index[values.eq(low)],marker]
                    test = stats.mannwhitneyu(a,b,alternative='two-sided',method='asymptotic')
                    pairwise.append(dict(Entity=entity,Marker=marker,Feature=feature,Contrast=str(high)+' vs '+str(low),MeanDifference=a.mean()-b.mean(),RankBiserial=2*test.statistic/(len(a)*len(b))-1,P=test.pvalue,N_A=len(a),N_B=len(b),Patients=metadata.loc[a.index.union(b.index),'PatientID'].nunique()))
cat = pd.DataFrame(categorical_results)
posthoc = pd.DataFrame(pairwise)
for frame in [cat,posthoc]:
    frame['Q'] = frame.groupby('Entity').P.transform(lambda p: stats.false_discovery_control(p.to_numpy(),method='bh'))
cat.to_csv(out/'statistics/categorical_overall_unadjusted.tsv',sep='\t',index=False)
posthoc.to_csv(out/'statistics/categorical_posthoc_unadjusted.tsv',sep='\t',index=False)
print('Categorical tests',len(cat),'posthoc',len(posthoc),flush=True)

# Step3. 连续变量 Spearman 及每患者首份活检敏感性，按疾病保留适用分数
correlations = []
one_biopsy = []
dominant = []
for feature,disease in continuous:
    valid = metadata[feature].notna() & (metadata.Disease.eq(disease) if disease!='All' else True)
    data = metadata.loc[valid].copy()
    ids = data.sort_index().drop_duplicates('PatientID').index
    for entity,abundance in matrices.items():
        for marker in abundance.columns:
            a = abundance.loc[data.index,marker]
            test = stats.spearmanr(a,data[feature])
            correlations.append(dict(Entity=entity,Marker=marker,Feature=feature,Disease=disease,Rho=test.statistic,P=test.pvalue,N=len(data),Patients=data.PatientID.nunique()))
            test = stats.spearmanr(abundance.loc[ids,marker],data.loc[ids,feature])
            one_biopsy.append(dict(Entity=entity,Marker=marker,Feature=feature,Disease=disease,Rho=test.statistic,P=test.pvalue,N=len(ids),Selection='Lexicographic first eligible biopsy per patient'))
    for ecotype in matrices['ecotype'].columns:
        own = data.DominantEcotype.eq(ecotype)
        a = data.loc[own,feature]
        b = data.loc[~own,feature]
        test = stats.mannwhitneyu(a,b,alternative='two-sided',method='asymptotic')
        dominant.append(dict(Ecotype=ecotype,Feature=feature,Disease=disease,MeanDifference=a.mean()-b.mean(),RankBiserial=2*test.statistic/(len(a)*len(b))-1,P=test.pvalue,N_Ecotype=len(a),N_Other=len(b),Patients_Ecotype=data.loc[own,'PatientID'].nunique(),Patients_Other=data.loc[~own,'PatientID'].nunique(),Method='Cell2021-guided dominant CE vs rest Wilcoxon; biopsy unit'))
corr = pd.DataFrame(correlations)
one = pd.DataFrame(one_biopsy)
dom = pd.DataFrame(dominant)
for frame in [corr,one]:
    frame['Q'] = frame.groupby('Entity').P.transform(lambda p: stats.false_discovery_control(p.to_numpy(),method='bh'))
dom['Q'] = stats.false_discovery_control(dom.P.to_numpy(),method='bh')
corr.to_csv(out/'statistics/continuous_spearman_unadjusted.tsv',sep='\t',index=False)
one.to_csv(out/'statistics/continuous_spearman_one_biopsy.tsv',sep='\t',index=False)
dom.to_csv(out/'statistics/dominant_ecotype_continuous_features.tsv',sep='\t',index=False)
print('Spearman',len(corr),'dominant',len(dom),flush=True)

# Step4. 聚类稳健线性模型和每患者等权敏感性，全部模型用匹配后的完整病例
adjusted = []
overall = []
model_failures = []
features = [(f,d,'continuous',None) for f,d in continuous]+[(f,'All','categorical',r) for f,r in categorical]
for feature,disease,kind,reference in features:
    eligible = metadata.Disease.eq(disease) if disease!='All' else pd.Series(True,index=metadata.index)
    covariates = [c for c in ['Disease','Region3','Age','Sex'] if c!=feature]
    data = metadata.loc[eligible].dropna(subset=[feature]+covariates).copy()
    covariates = [c for c in covariates if data[c].nunique()>1]
    base = pd.DataFrame({'Intercept':1.0},index=data.index)
    for covariate in covariates:
        if covariate=='Age':
            base['Age_per10'] = data.Age/10
        else:
            base = pd.concat([base,pd.get_dummies(data[covariate],prefix=covariate,drop_first=True,dtype=float)],axis=1)
    if kind=='continuous':
        feature_sd = data[feature].std(ddof=1)
        exposure = ((data[feature]-data[feature].mean())/feature_sd).to_frame('Exposure')
        labels = {'Exposure':'per 1 SD '+feature}
    else:
        levels = [reference]+[v for v in sorted(data[feature].unique()) if v!=reference]
        exposure = pd.get_dummies(pd.Categorical(data[feature],categories=levels),drop_first=True,dtype=float)
        exposure.index = data.index
        exposure.columns = ['Exposure_'+str(v) for v in levels[1:]]
        labels = dict(zip(exposure.columns,[str(v)+' vs '+str(reference) for v in levels[1:]]))
        feature_sd = np.nan
    design = pd.concat([base,exposure],axis=1).astype(float)
    assert np.linalg.matrix_rank(design)==design.shape[1]
    group = data.PatientID
    patients = group.nunique()
    counts = group.map(group.value_counts())
    restriction = np.zeros((exposure.shape[1],design.shape[1]))
    restriction[:, -exposure.shape[1]:] = np.eye(exposure.shape[1])
    for entity,abundance in matrices.items():
        for marker in abundance.columns:
            original = abundance.loc[data.index,marker]
            outcome_sd = original.std(ddof=1)
            y = (original-original.mean())/outcome_sd
            for weighting in ['Biopsy_weighted','Equal_patient_weight']:
                model = sm.OLS(y,design) if weighting=='Biopsy_weighted' else sm.WLS(y,design,weights=1/counts)
                fit = model.fit(cov_type='cluster',cov_kwds={'groups':group,'use_correction':True,'df_correction':True},use_t=True)
                test = fit.f_test(restriction)
                overall.append(dict(Entity=entity,Marker=marker,Feature=feature,Disease=disease,Weighting=weighting,F=float(np.asarray(test.fvalue).item()),P=float(test.pvalue),DfNumerator=exposure.shape[1],DfDenominator=patients-1,N=len(data),Patients=patients,Covariates=';'.join(covariates)))
                for column in exposure.columns:
                    ci = fit.conf_int().loc[column]
                    adjusted.append(dict(Entity=entity,Marker=marker,Feature=feature,Disease=disease,Contrast=labels[column],Weighting=weighting,Beta=fit.params[column],SE=fit.bse[column],CI_low=ci.iloc[0],CI_high=ci.iloc[1],P=fit.pvalues[column],Beta_abundance_units=fit.params[column]*outcome_sd,OutcomeSD=outcome_sd,FeatureSD=feature_sd,N=len(data),Patients=patients,Covariates=';'.join(covariates),ModelRank=design.shape[1],R_squared=fit.rsquared))
    print('Cluster adjusted',feature,'N',len(data),'patients',patients,flush=True)
adj = pd.DataFrame(adjusted)
omnibus = pd.DataFrame(overall)
for frame in [adj,omnibus]:
    frame['Q'] = frame.groupby(['Entity','Weighting']).P.transform(lambda p: stats.false_discovery_control(p.to_numpy(),method='bh'))
adj.to_csv(out/'statistics/patient_cluster_adjusted_contrasts.tsv',sep='\t',index=False)
omnibus.to_csv(out/'statistics/patient_cluster_adjusted_overall.tsv',sep='\t',index=False)

# Step5. 原文式 state phenotype 富集：父类型内逐样本置换，不混入丰度模型
rng = np.random.default_rng(20261015)
enrichments = []
for cell_type in members.CellType.unique():
    assignments = pd.read_csv(root/'ctspec_region_ibd'/cell_type/'state_assignment.txt',sep='\t')
    assert assignments.ID.is_unique
    data = assignments.set_index('ID').join(metadata,validate='one_to_one')
    states = members.loc[members.CellType.eq(cell_type),'State'].tolist()
    for feature,reference in categorical:
        sub = data.dropna(subset=[feature]).copy()
        state_codes = pd.Categorical(sub.State,categories=states).codes
        levels = sorted(sub[feature].unique())
        level_codes = pd.Categorical(sub[feature],categories=levels).codes
        observed = np.bincount(state_codes*len(levels)+level_codes,minlength=len(states)*len(levels)).reshape(len(states),len(levels))
        null = np.zeros((1000,len(states),len(levels)))
        for iteration in range(1000):
            shuffled = rng.permutation(state_codes)
            null[iteration] = np.bincount(shuffled*len(levels)+level_codes,minlength=len(states)*len(levels)).reshape(len(states),len(levels))
        mean = null.mean(axis=0)
        sd = null.std(axis=0,ddof=1)
        z = np.divide(observed-mean,sd,out=np.full_like(mean,np.nan),where=sd>0)
        for i,state in enumerate(states):
            for j,level in enumerate(levels):
                enrichments.append(dict(Marker=cell_type+'_'+state,CellType=cell_type,Feature=feature,Level=level,Observed=int(observed[i,j]),NullMean=mean[i,j],NullSD=sd[i,j],Z=z[i,j],P=stats.norm.sf(z[i,j]),EmpiricalP=(1+np.sum(null[:,i,j]>=observed[i,j]))/1001,N=len(sub),Patients=sub.PatientID.nunique(),Permutations=1000,Seed=20261015,Unit='Biopsy; patient dependence unadjusted'))
enrich = pd.DataFrame(enrichments)
valid = enrich.P.notna()
enrich.loc[valid,'Q'] = stats.false_discovery_control(enrich.loc[valid,'P'].to_numpy(),method='bh')
enrich.loc[valid,'EmpiricalQ'] = stats.false_discovery_control(enrich.loc[valid,'EmpiricalP'].to_numpy(),method='bh')
enrich.to_csv(out/'statistics/state_phenotype_permutation_enrichment.tsv',sep='\t',index=False)

# Step6. 同患者同肠段炎症配对，先患者内等权汇总，不把肠段当独立患者
paired_results = []
for entity,abundance in matrices.items():
    merged = metadata[['PatientID','Region3','Inflammation']].join(abundance)
    medians = merged.groupby(['PatientID','Region3','Inflammation'])[abundance.columns].median()
    inflamed = medians.xs('Inflamed',level='Inflammation')
    quiet = medians.xs('Non-inflamed',level='Inflammation')
    both = inflamed.index.intersection(quiet.index)
    delta = (inflamed.loc[both]-quiet.loc[both]).groupby('PatientID').mean()
    delta.to_csv(out/('source_data/paired_inflammation_patient_deltas_'+entity+'.tsv'),sep='\t')
    for marker in abundance.columns:
        test = stats.wilcoxon(delta[marker],alternative='two-sided',zero_method='wilcox',method='auto')
        paired_results.append(dict(Entity=entity,Marker=marker,Feature='Inflamed_minus_Noninflamed',Patients=len(delta),PatientRegionPairs=len(both),MedianPatientDifference=delta[marker].median(),MeanPatientDifference=delta[marker].mean(),P=test.pvalue,Method='Patient-region medians, equal region mean per patient, paired Wilcoxon'))
paired = pd.DataFrame(paired_results)
paired['Q'] = paired.groupby('Entity').P.transform(lambda p: stats.false_discovery_control(p.to_numpy(),method='bh'))
paired.to_csv(out/'statistics/patient_region_matched_inflammation.tsv',sep='\t',index=False)

# Step7. 与历史主要结果核对，保存统计完成检查点
previous = root/'06_validation/paper_aligned_recovery_characterization_20261012/06_ecotype_features/ecotype_discrete_features.tsv'
old = pd.read_csv(previous,sep='\t')
check = cat.loc[cat.Entity.eq('ecotype') & cat.Feature.isin(['Disease','Inflammation'])].merge(old,left_on=['Marker','Feature'],right_on=['Ecotype','Feature'],suffixes=('_new','_old'),validate='one_to_one')
check['AbsPdiff'] = (check.P_new-check.P_old).abs()
check.to_csv(out/'audit/previous_ecotype_disease_inflammation_recheck.tsv',sep='\t',index=False)
summary = dict(CategoricalOverall=len(cat),Pairwise=len(posthoc),Spearman=len(corr),DominantContinuous=len(dom),ClusterContrasts=len(adj),ClusterOverall=len(omnibus),PermutationRows=len(enrich),PairedTests=len(paired),PairedPatients=int(paired.Patients.max()),PreviousMaxAbsPdiff=float(check.AbsPdiff.max()))
(out/'audit/phase4_statistics_checkpoint.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
(out/'progress_summary.md').write_text('# 2026-10-15 进度\n\nPhase 1–4 完成。临床字典、完整统计表、患者聚类调整、患者等权与配对敏感性已保存；图件和最终汇报摘要制作中。\n',encoding='utf-8')
print(json.dumps(summary,indent=2),flush=True)
