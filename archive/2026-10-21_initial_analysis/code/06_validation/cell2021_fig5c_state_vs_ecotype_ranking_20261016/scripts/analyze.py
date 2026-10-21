"""Read frozen inputs only. Write all new statistics under this script's parent output.
Run: D:/anaconda3/python.exe scripts/analyze.py
No Recovery, NMF, BayesPrism, or existing statistical script is executed.
"""
from pathlib import Path
import hashlib, json, platform, sys, importlib.metadata
import numpy as np
import pandas as pd
from scipy import stats

OUT = Path(__file__).resolve().parents[1]
ROOT = Path('D:/IBD_ecotype/IBD_ecotype_analysis')
EXT = ROOT / '06_validation/external_clinical_outcomes_20261015'
DISC = ROOT / '06_validation/MSCCR_joint_clinical_audit_20261015/source_data'
CORE = ROOT / 'ctspec_region_ibd'
SEED = 20261016
for folder in ['statistics', 'manifests', 'audit', 'figures', 'reference']:
    (OUT / folder).mkdir(exist_ok=True, parents=True)
sources = {}
checks = []

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def register(path, role):
    path = Path(path).resolve()
    assert path.exists(), str(path)
    if str(path) not in sources:
        sources[str(path)] = dict(Path=str(path), Role=role, Bytes=path.stat().st_size, SHA256_before=sha(path))
    return path

def read(path, role, index_col=None):
    return pd.read_csv(register(path, role), sep='\t', index_col=index_col, float_precision='round_trip')

def save(frame, relative):
    path = OUT / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep='\t', index=False, float_format='%.17g')

def check(name, condition, observed, scope):
    checks.append(dict(Check=name, Pass=bool(condition), Observed=str(observed), Scope=scope))
    assert condition, f'{name}: {observed}'

members = read(EXT/'audit/frozen_state_ecotype_membership.tsv', 'Frozen 53-state / 7-ecotype membership')
features = ['E'+str(i) for i in range(1,8)] + members.ID.tolist()
e1_members = members.loc[members.Ecotype.eq('E1'),'ID'].tolist()
check('60 unique features',len(set(features))==60,len(set(features)),'7 ecotypes and 53 states')

def bh(p, size):
    p = np.asarray(p, float)
    valid = np.flatnonzero(np.isfinite(p))
    q = np.full(p.shape, np.nan)
    order = valid[np.argsort(p[valid], kind='stable')]
    q[order] = np.minimum(1, np.minimum.accumulate((p[order]*size/np.arange(1,len(order)+1))[::-1])[::-1])
    return q

def test_features(scores, manifest, endpoint, scope, favorable_event=1):
    ids = manifest.SampleID.tolist()
    x = scores.loc[ids, features]
    y = manifest.Outcome.to_numpy(dtype=int)
    n1, n0 = int(y.sum()), int((y==0).sum())
    check(endpoint+' unique IDs',manifest.SampleID.is_unique and manifest.PatientID.is_unique,len(manifest),scope)
    check(endpoint+' binary and finite',set(y)=={0,1} and np.isfinite(x.to_numpy()).all(),(len(y),n1),scope)
    rows=[]
    for f in features:
        a=x.loc[y==1,f].to_numpy(); b=x.loc[y==0,f].to_numpy()
        u,p=stats.mannwhitneyu(a,b,alternative='two-sided',method='asymptotic',use_continuity=True)
        auc=u/(n1*n0); rb=2*auc-1
        direction=int(np.sign(rb))
        z_abs=float(stats.norm.isf(max(p,np.finfo(float).tiny)/2))
        row=dict(Endpoint=endpoint,Scope=scope,Feature=f,Family='Ecotype' if f.startswith('E') and f[1:].isdigit() else 'State',
                 Ecotype=members.set_index('ID').Ecotype.get(f,f),E1_member=f in e1_members,
                 N=len(y),Events=n1,Nonevents=n0,EffectiveN=len(y),N_event=n1,N_nonevent=n0,
                 P_two_sided=p,Z_two_sided_abs=z_abs,SignedZ_event=direction*z_abs,
                 SignedZ_favorable=direction*z_abs*(1 if favorable_event else -1),
                 EventDirection=direction,Direction='Higher in event=1' if direction>0 else 'Lower in event=1' if direction<0 else 'No U shift',
                 U_event=u,RankBiserial=rb,AUC_like_event=auc,
                 Median_event=float(np.median(a)),Median_nonevent=float(np.median(b)),
                 MedianDifference=float(np.median(a)-np.median(b)),MeanDifference=float(a.mean()-b.mean()),
                 NA_N=int(x[f].isna().sum()),ExactZeroN=int((x[f]==0).sum()),NearZeroN_1e8=int((x[f]<=1e-8).sum()),
                 UniqueScores=x[f].nunique(),TiedObservationExcess=len(x)-x[f].nunique(),
                 FavorableEvent=favorable_event,InputStatus='Conditional_processed_scale' if endpoint.startswith('RISK') else 'Frozen_scores')
        rows.append(row)
    frame=pd.DataFrame(rows)
    frame['Q_BH60']=bh(frame.P_two_sided,60)
    frame['Q_within_family_recomputed']=np.nan
    for fam,size in [('Ecotype',7),('State',53)]:
        mask=frame.Family.eq(fam); frame.loc[mask,'Q_within_family_recomputed']=bh(frame.loc[mask,'P_two_sided'],size)
    frame['Rank_absZ']=frame.Z_two_sided_abs.rank(method='average',ascending=False)
    frame['Rank_signed_favorable']=frame.SignedZ_favorable.rank(method='average',ascending=False)
    frame['SignedZ_paper_style']=-frame.SignedZ_favorable
    frame['Rank_signed_event']=frame.SignedZ_event.rank(method='average',ascending=False)
    frame['DirectionReference']='Favorable clinical outcome'
    if scope.startswith('Discovery'):
        frame['DirectionReference']='Inflamed-high; clinical favorable-outcome ranking not applicable'
        frame[['SignedZ_favorable','Rank_signed_favorable','SignedZ_paper_style','FavorableEvent']]=np.nan
    frame['P_rank_tie_size']=frame.groupby('P_two_sided').P_two_sided.transform('size')
    frame=frame.sort_values(['Z_two_sided_abs','Feature'],ascending=[False,True]).reset_index(drop=True)
    frame['DisplayOrder_absZ']=np.arange(1,61)
    check(endpoint+' all 60 same N',frame.EffectiveN.nunique()==1 and len(frame)==60,len(frame),scope)
    check(endpoint+' BH verification',np.allclose(frame.Q_BH60,stats.false_discovery_control(frame.P_two_sided.to_numpy(),method='bh')),float(frame.Q_BH60.min()),scope)
    # Independently reconstruct U from pooled average ranks and asymptotic p.
    max_u=max_p=0.0
    for row in frame.itertuples():
        values=x[row.Feature].to_numpy(); rank=stats.rankdata(values,method='average')
        u0=rank[y==1].sum()-n1*(n1+1)/2
        _,counts=np.unique(values,return_counts=True)
        variance=n1*n0/12*((len(y)+1)-np.sum(counts**3-counts)/(len(y)*(len(y)-1)))
        z=(abs(u0-n1*n0/2)-0.5)/np.sqrt(variance)
        p0=min(1.,2*stats.norm.sf(z))
        max_u=max(max_u,abs(u0-row.U_event)); max_p=max(max_p,abs(p0-row.P_two_sided))
    check(endpoint+' independent U/P',max_u==0 and max_p<1e-12,(max_u,max_p),'Pooled-rank and tie-corrected variance reconstruction')
    return frame

definitions=[
 ('ETR_W10_clinical_remission','01_GSE72819_treatment_outcome',70,12,'Primary',1,'Week 10 clinical remission: 1=remitter, 0=non-remitter'),
 ('IFX_W4_6_endoscopic_healing','03_GSE73661_mucosal_healing',23,8,'Primary',1,'Follow-up week 4-6 Mayo endoscopic subscore: 1=0/1, 0=2/3'),
 ('VDZ_W6_endoscopic_healing','03_GSE73661_mucosal_healing',27,6,'Primary',1,'Week 6 Mayo endoscopic subscore: 1=0/1, 0=2/3'),
 ('RISK_3yr_B2_B3','02_RISK_disease_progression',210,27,'Primary_conditional',0,'3-year progression: 1=B2/B3, 0=no B2/B3; conditional deposited scale'),
 ('VDZ_W12_endoscopic_healing','03_GSE73661_mucosal_healing',13,3,'Supplementary',1,'Week 12 Mayo endoscopic subscore: 1=0/1, 0=2/3; very small event group'),
 ('VDZ_W52_endoscopic_healing','03_GSE73661_mucosal_healing',19,13,'Supplementary',1,'Week 52 Mayo endoscopic subscore: 1=0/1, 0=2/3')]
tables=[]; denominators=[]; score_qa=[]; patient_manifests=[]; formula=[]
history=read(EXT/'04_cross_cohort_summary/all_associations.tsv','Historical models / P and separate BH7/BH53')
for endpoint,folder,n_expected,event_expected,role,favorable,definition in definitions:
    eco=read(EXT/folder/'continuous_ecotype_scores.tsv','Frozen baseline Ecotype scores',0)
    state=read(EXT/folder/'continuous_state_scores.tsv','Frozen baseline State scores',0)
    m=read(EXT/folder/(endpoint+'_analysis_patients.tsv'),'Endpoint-specific frozen patients and outcome')
    full=read(EXT/folder/'cohort_endpoint_manifest.tsv','Full included/excluded baseline endpoint manifest')
    check(endpoint+' matrix IDs identical',eco.index.equals(state.index),(eco.shape,state.shape),'Full cohort matrices')
    check(endpoint+' expected N/events',len(m)==n_expected and m.Outcome.sum()==event_expected,(len(m),m.Outcome.sum()),definition)
    check(endpoint+' baseline only',m.Timepoint.isin(['Baseline','W0','Diagnosis; treatment naive']).all(),m.Timepoint.unique(),role)
    if endpoint.startswith(('IFX','VDZ')):
        check(endpoint+' actual Mayo coding',np.array_equal(m.Outcome.to_numpy(),(m.FollowupMayo<=1).astype(int).to_numpy()),m.FollowupMayo.value_counts().to_dict(),'Outcome column, never baseline HealingBinary')
    expected=pd.DataFrame({e:state.loc[:,members.loc[members.Ecotype.eq(e),'ID']].mean(axis=1) for e in eco.columns})
    expected=expected.div(expected.sum(axis=1),axis=0)
    err=float(np.max(np.abs(eco-expected)))
    check(endpoint+' ecotype formula',err<1e-12,err,'Member mean H then normalization across 7 ecotypes')
    formula.append(dict(Endpoint=endpoint,MaxAbsError=err,Formula='Member mean H / sum of seven member means'))
    scores=pd.concat([eco,state],axis=1)
    t=test_features(scores,m,endpoint,'External baseline association',favorable)
    old=history.loc[history.Endpoint.eq(endpoint),['Feature','P_Wilcoxon','Q_Wilcoxon','P_Firth','Q_Firth']].rename(columns={'Q_Wilcoxon':'Historical_Q_Wilcoxon_BH7_or_BH53','Q_Firth':'Historical_Q_Firth_BH7_or_BH53','P_Firth':'Historical_P_Firth','P_Wilcoxon':'Historical_P_Wilcoxon'})
    t=t.merge(old,on='Feature',how='left',validate='one_to_one')
    check(endpoint+' historical P match',np.allclose(t.P_two_sided,t.Historical_P_Wilcoxon,rtol=1e-10,atol=1e-13),float((t.P_two_sided-t.Historical_P_Wilcoxon).abs().max()),'Same input, asymptotic ties/continuity policy')
    t['Role']=role; t['OutcomeDefinition']=definition; tables.append(t)
    m=m.copy();m['RankingRole']=role;m['CorrectedOutcomeDefinition']=definition;m['All60ScoresComplete']=True
    patient_manifests.append(m)
    save(m,'manifests/'+endpoint+'_patients.tsv')
    export=m[['SampleID','PatientID','Outcome','Timepoint']].merge(scores.reset_index(),on='SampleID',validate='one_to_one')
    save(export,'manifests/'+endpoint+'_patient_scores60.tsv')
    denominators.append(dict(Endpoint=endpoint,Role=role,N=len(m),Events=int(m.Outcome.sum()),Nonevents=int((m.Outcome==0).sum()),Candidates=int(full.Endpoint.eq(endpoint).sum()),EffectiveN_min=t.EffectiveN.min(),EffectiveN_max=t.EffectiveN.max(),NA_cells=t.NA_N.sum(),ExactZeroCells=t.ExactZeroN.sum(),TiedObservationExcess=t.TiedObservationExcess.sum(),Definition=definition))
    score_qa.append(t[['Endpoint','Feature','Family','N','Events','EffectiveN','NA_N','ExactZeroN','NearZeroN_1e8','UniqueScores','TiedObservationExcess']])

external=pd.concat(tables,ignore_index=True)
save(external,'statistics/external_all60_rankings.tsv')
save(pd.DataFrame(denominators),'audit/endpoint_comparability.tsv')
save(pd.concat(score_qa),'audit/external_per_feature_score_QA.tsv')
save(pd.DataFrame(formula),'audit/ecotype_formula_reconstruction.tsv')
save(pd.concat(patient_manifests,ignore_index=True),'manifests/external_all_endpoint_patient_manifest.tsv')
overlaps=[]
for i,a in enumerate(patient_manifests):
    for b in patient_manifests[i+1:]:
        overlaps.append(dict(EndpointA=a.Endpoint.iloc[0],EndpointB=b.Endpoint.iloc[0],SharedSampleIDs=len(set(a.SampleID)&set(b.SampleID)),SharedPatientIDs=len(set(a.PatientID)&set(b.PatientID))))
save(pd.DataFrame(overlaps),'audit/endpoint_patient_overlap.tsv')

# Discovery: a single biopsy per patient chosen without observing outcome or scores.
meta=read(DISC/'analysis_metadata_2025.tsv','Discovery sample-patient-inflammation metadata')
ds=read(DISC/'state_abundance_2025x53.tsv','Full precision frozen Discovery state H',0)
de=read(DISC/'ecotype_abundance_2025x7.tsv','Full precision frozen Discovery Ecotype abundance',0)
dm=read(DISC/'frozen_state_membership.tsv','Discovery frozen membership')
check('Discovery exact membership',members.equals(dm),len(dm),'All five fields and row order')
check('Discovery matrix IDs',ds.index.equals(de.index) and ds.index.tolist()==meta.ID.tolist(),len(ds),'2025 biopsies / 919 patients')
check('Discovery dimensions',len(meta)==2025 and meta.PatientID.nunique()==919,(len(meta),meta.PatientID.nunique()),'Frozen Discovery')
native=[]
for parent in members.CellType.unique():
    h=read(CORE/parent/'state_abundances.txt','Native full precision Discovery H',0).T
    h.columns=parent+'_'+h.columns
    native.append(h)
native=pd.concat(native,axis=1).loc[ds.index,ds.columns]
maxerr=float(np.abs(native-ds).to_numpy().max())
check('Discovery native H precision',np.allclose(native,ds,atol=1e-12,rtol=1e-10),maxerr,'All 107325 cells')
check('Discovery native H ranks',all(stats.rankdata(native[f]).tolist()==stats.rankdata(ds[f]).tolist() for f in ds),maxerr,'All 53 per-feature rank orders')
ne=read(CORE/'Ecotypes/ecotype_abundance.txt','Native frozen Discovery Ecotype abundance',0).T.loc[de.index,de.columns]
check('Discovery native Ecotype',np.allclose(ne,de,atol=1e-12,rtol=1e-10),float(np.abs(ne-de).to_numpy().max()),'All 14175 cells')
native_members=read(CORE/'Ecotypes/ecotypes.txt','Native frozen Discovery membership')
check('Discovery native members',native_members.equals(members),len(native_members),'Same 53 states / 7 groups')
scores=pd.concat([de,ds],axis=1)
discovery_qc=pd.DataFrame({'Feature':scores.columns,'N':len(scores),'NA_N':scores.isna().sum().to_numpy(),'NonfiniteN':(~np.isfinite(scores)).sum().to_numpy(),'ExactZeroN':scores.eq(0).sum().to_numpy(),'NearZeroN_1e8':scores.le(1e-8).sum().to_numpy(),'UniqueScores':scores.nunique().to_numpy()})
save(discovery_qc,'audit/discovery_all2025_score_QA.tsv')
meta['SampleID']=meta.ID
meta['Outcome']=meta.Inflammation.map({'Inflamed':1,'Non-inflamed':0})
meta['SelectionHash']=meta.apply(lambda r:hashlib.sha256(f'{SEED}|{r.PatientID}|{r.SampleID}'.encode()).hexdigest(),axis=1)
selection=meta.sort_values(['PatientID','SelectionHash','SampleID']).drop_duplicates('PatientID').copy()
selected=test_features(scores,selection,'MSCCR_Inflammation_one_biopsy','Discovery internal association',1)
selected['Role']='Discovery_internal'
save(selected,'statistics/discovery_one_biopsy_all60_rankings.tsv')
manifest_cols=['SampleID','PatientID','Disease','Region3','RegionDetailed','Inflammation','Outcome','SelectionHash']
save(selection[manifest_cols],'manifests/discovery_one_biopsy_919_patient_manifest.tsv')
save(selection[manifest_cols].merge(scores.rename_axis('SampleID').reset_index(),on='SampleID',validate='one_to_one'),'manifests/discovery_one_biopsy_patient_scores60.tsv')
meta['SelectedPrimary']=meta.SampleID.isin(selection.SampleID)
save(meta[manifest_cols+['SelectedPrimary']],'manifests/discovery_all2025_biopsy_manifest.tsv')

# Repeated biopsy selection sensitivity: every feature uses the same 919 patients.
rng=np.random.default_rng(SEED)
groups=[g.SampleID.to_numpy() for _,g in meta.groupby('PatientID',sort=True)]
sens=[]; sensitivity_ids=[]
for iteration in range(100):
    ids=[g[int(rng.integers(len(g)))] for g in groups]
    sub=meta.set_index('SampleID').loc[ids]
    yy=sub.Outcome.to_numpy(); xx=scores.loc[ids,features].to_numpy();n1=int(yy.sum());n0=len(yy)-n1
    u,p=stats.mannwhitneyu(xx[yy==1],xx[yy==0],axis=0,alternative='two-sided',method='asymptotic',use_continuity=True)
    ranks=stats.rankdata(p,method='average')
    for j,f in enumerate(features):
        sens.append(dict(Iteration=iteration+1,Feature=f,N=len(yy),Events=n1,P_two_sided=p[j],Rank_absZ=ranks[j],RankBiserial=2*u[j]/(n1*n0)-1))
    sensitivity_ids.extend(dict(Iteration=iteration+1,PatientID=p,SampleID=s) for p,s in zip(sub.PatientID,ids))
sens=pd.DataFrame(sens)
save(sens,'statistics/discovery_100_selection_sensitivity.tsv')
save(pd.DataFrame(sensitivity_ids),'manifests/discovery_100_selection_sampleIDs.tsv')
summary=sens.groupby('Feature').agg(RankMin=('Rank_absZ','min'),RankP10=('Rank_absZ',lambda x:x.quantile(.1)),RankMedian=('Rank_absZ','median'),RankP90=('Rank_absZ',lambda x:x.quantile(.9)),RankMax=('Rank_absZ','max'),RB_min=('RankBiserial','min'),RB_median=('RankBiserial','median'),RB_max=('RankBiserial','max')).reset_index()
e1r=sens.loc[sens.Feature.eq('E1'),['Iteration','Rank_absZ','RankBiserial']]
cmp=sens[sens.Feature.isin(e1_members)].merge(e1r,on='Iteration',suffixes=('_member','_E1'))
pairsummary=cmp.groupby('Feature').agg(E1_stronger_absRank_N=('Rank_absZ_member',lambda r:0)).reset_index().drop(columns='E1_stronger_absRank_N')
pairsummary['E1_stronger_absRank_N']=[int(((cmp.Feature==f)&(cmp.Rank_absZ_E1<cmp.Rank_absZ_member)).sum()) for f in pairsummary.Feature]
pairsummary['E1_stronger_absRB_N']=[int(((cmp.Feature==f)&(cmp.RankBiserial_E1.abs()>cmp.RankBiserial_member.abs())).sum()) for f in pairsummary.Feature]
save(summary,'statistics/discovery_selection_rank_sensitivity_summary.tsv')
save(pairsummary,'statistics/discovery_E1_member_selection_comparison.tsv')

# Paired-feature bootstrap on the same primary selected patients, stratified by outcome.
# Descriptive intervals for rank-biserial differences, not external validation.
brng=np.random.default_rng(SEED+1)
sf=['E1']+e1_members
sx=scores.loc[selection.SampleID,sf].to_numpy()
sy=selection.Outcome.to_numpy()
a=sx[sy==1];b=sx[sy==0]
bootstrap=[]
for i in range(2000):
    aa=a[brng.integers(len(a),size=len(a))];bb=b[brng.integers(len(b),size=len(b))]
    uu=stats.mannwhitneyu(aa,bb,axis=0,alternative='two-sided',method='asymptotic',use_continuity=True).statistic
    bootstrap.append(2*uu/(len(a)*len(b))-1)
bootstrap=np.array(bootstrap)
effect=selected.set_index('Feature').RankBiserial
bootrows=[]
for j,f in enumerate(sf):
    delta=bootstrap[:,0]-bootstrap[:,j]
    bootrows.append(dict(Feature=f,RankBiserial=effect[f],RB_CI_low=np.quantile(bootstrap[:,j],.025),RB_CI_high=np.quantile(bootstrap[:,j],.975),E1_minus_feature_RB=effect['E1']-effect[f],Delta_CI_low=np.quantile(delta,.025),Delta_CI_high=np.quantile(delta,.975),BootstrapN=2000,Seed=SEED+1,Interval='Percentile 95%; same subjects for all features; stratified event/nonevent resampling; exploratory unadjusted intervals'))
save(pd.DataFrame(bootrows),'statistics/discovery_E1_member_effect_bootstrap.tsv')

# Historic clustered adjusted models: preserve the original family definitions.
adjusted=[]
for filename,sourcefamily in [('Fig1_state_clinical_adjusted_full_statistics.tsv',1166),('Fig2_ecotype_clinical_adjusted_full_statistics.tsv',154)]:
    d=read(DISC/filename,'Historical F062/F063 clustered adjusted model')
    d=d[d.Feature.eq('Inflammation')].copy();d['HistoricalBHFamilySize']=sourcefamily;adjusted.append(d)
save(pd.concat(adjusted),'statistics/discovery_historical_clustered_beta_separate_model.tsv')
register(EXT/'02_RISK_disease_progression/ecotype_age_sex_sensitivity.tsv','RISK prespecified age/sex sensitivity')
for f in ['scripts/03_fit_outcomes.py','scripts/02_prepare_GSE73661.R','scripts/00_prepare_runtime.py','audit/immutable_input_hashes_verified.tsv','audit/model_gene_coverage_and_projection_QC.tsv','audit/delivery_file_manifest.tsv','statistical_analysis_plan.md']:
    register(EXT/f,'External frozen provenance / historical method')
register(ROOT/'06_validation/Cell2021_full_pipeline_audit_20261016/FULL_PIPELINE_METHOD_AUDIT.md','2026-10-16 methodological audit')
for parent in members.CellType.unique():
    for filename in ['W.txt','mapping_to_initial_states.txt']:
        candidate=CORE/parent/filename
        if candidate.exists():register(candidate,'Local frozen Discovery model reference; external hash chain documented separately')

immutable=read(EXT/'audit/immutable_input_hashes_verified.tsv','Archived 45 immutable external model input hashes')
snapshot=ROOT/'06_validation/GSE72819_recovery_review_20261009/raw_server/audit_snapshot/frozen_model'
wqa=[]
for w in sorted(snapshot.glob('*/*/W.txt')):
    register(w,'Archived actual frozen W used by Recovery')
    suffix='/'+w.parent.parent.name+'/'+w.parent.name+'/W.txt'
    matching=immutable[immutable.Path.str.endswith(suffix)]
    matched=len(matching)==1 and sha(w)==matching.SHA256.iloc[0]
    wqa.append(dict(Parent=w.parent.parent.name,LocalPath=str(w),SHA256=sha(w),ArchivedExternalSHA256=matching.SHA256.iloc[0] if len(matching)==1 else '',ExactMatch=matched))
check('Same 12 frozen W hashes',len(wqa)==12 and all(r['ExactMatch'] for r in wqa),len(wqa),'Local archived W vs external immutable model evidence')
save(pd.DataFrame(wqa),'audit/frozen_W_exact_hash_matches.tsv')

combined=pd.concat([external,selected],ignore_index=True)
save(combined[combined.DisplayOrder_absZ<=10],'statistics/top10_each_endpoint.tsv')
save(combined[combined.Rank_absZ<=10],'statistics/top_ranks_including_ties.tsv')
save(combined[combined.Family.eq('Ecotype')],'statistics/all_7_ecotype_rankings.tsv')
save(combined[combined.Feature.isin(['E1']+e1_members)],'statistics/E1_and_10_members_all_endpoints.tsv')
rankgroups=[]
for endpoint,g in combined.groupby('Endpoint',sort=False):
    for family,part in g.groupby('Family'):
        rankgroups.append(dict(Endpoint=endpoint,Family=family,FeatureN=len(part),RankMin=part.Rank_absZ.min(),RankQ25=part.Rank_absZ.quantile(.25),RankMedian=part.Rank_absZ.median(),RankQ75=part.Rank_absZ.quantile(.75),RankMax=part.Rank_absZ.max(),Rank_le10_N=int((part.Rank_absZ<=10).sum()),Rank_le10_Fraction=float((part.Rank_absZ<=10).mean()),MedianAbsRB=part.RankBiserial.abs().median(),MedianAbsZ=part.Z_two_sided_abs.median(),BH60_lt005=int((part.Q_BH60<.05).sum())))
save(pd.DataFrame(rankgroups),'statistics/ecotype_state_rank_distribution_descriptive.tsv')
save(pd.DataFrame(checks),'audit/statistical_QA.tsv')
sourceframe=pd.DataFrame(sources.values())
sourceframe['SHA256_after']=sourceframe.Path.map(sha)
sourceframe['Unchanged']=sourceframe.SHA256_before.eq(sourceframe.SHA256_after)
assert sourceframe.Unchanged.all()
save(sourceframe,'audit/input_SHA256_before_after.tsv')
env=dict(Python=sys.version,Executable=sys.executable,Platform=platform.platform(),Seed=SEED,Packages={p:importlib.metadata.version(p) for p in ['numpy','pandas','scipy','matplotlib','statsmodels','PyMuPDF']},Test='scipy.stats.mannwhitneyu; two-sided; asymptotic; ties corrected; continuity=True',Z='sign(U-n1*n0/2)*norm.isf(P_two_sided/2)',Rank='Descending abs(Z), average ranks for equal P; lexical feature name only resolves display order',MultipleTesting='BH60 per endpoint; separate recomputed BH7/BH53 and historical q retained; no pooled endpoint inference',Selection='Primary: minimum SHA256(seed|PatientID|SampleID); 100 sensitivity draws: numpy default_rng(seed), one random biopsy per patient',NoGlobalLiptak='Heterogeneous drugs/endpoints/diseases and overlapping VDZ follow-up samples; no combined score')
(OUT/'audit/statistical_environment_and_methods.json').write_text(json.dumps(env,ensure_ascii=False,indent=2),encoding='utf-8')
print('COMPLETED',len(external),'external tests;',len(selected),'Discovery tests;',len(sens),'sensitivity rows; QA',len(checks),'PASS')
print(pd.DataFrame(rankgroups).to_string(index=False))
print(selected[selected.Feature.isin(['E1']+e1_members)][['Feature','Rank_absZ','P_two_sided','RankBiserial']].to_string(index=False))
