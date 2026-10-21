# 核验旧表并补充冻结状态的表型富集，不重新拟合 Discovery 或 Recovery。
# 服务器 ibd_refaudit 环境运行，参数为项目根目录；所有新结果保存在脚本上一级目录。
from pathlib import Path
import sys, json, hashlib, platform
import numpy as np
import pandas as pd
from scipy import stats

# Step1. 核对冻结模型、分配与表型，只使用通过QC的父类型内状态
root = Path(sys.argv[1])
out = Path(__file__).resolve().parents[1]
for folder in ['statistics', 'audit', 'logs']:
    (out/folder).mkdir(exist_ok=True)
prior = root/'06_validation/paper_aligned_recovery_characterization_20261012'
frozen = root/'05_EcoTyper/output/MSCCR_CORE12/ctspec_region_ibd'
if not frozen.exists():
    frozen = root/'ctspec_region_ibd'
members = pd.read_csv(prior/'00_audit/frozen_states.tsv', sep='\t')
annotation_path = root/'05_EcoTyper/input/MSCCR_CORE12/region_ibd/annotation/MSCCR_CORE12_IBD_annotation.tsv'
if not annotation_path.exists():
    annotation_path = root/'06_validation/MSCCR_joint_clinical_audit_20261015/source_data/server_MSCCR_CORE12_IBD_annotation.tsv'
annotation = pd.read_csv(annotation_path, sep='\t')
assert len(annotation)==2025 and annotation.ID.is_unique and annotation.PatientID.nunique()==919
assert len(members)==53 and members.Ecotype.nunique()==7 and members.CellType.nunique()==12
input_paths = [annotation_path, prior/'00_audit/frozen_states.tsv']
assignments = []
coverage = []
for ct, states in members.groupby('CellType', sort=False):
    p = frozen/ct/'state_assignment.txt'
    a = pd.read_csv(p, sep='\t')[['ID', 'State']]
    assert a.ID.is_unique and set(a.State)<=set(states.State) and set(a.ID)<=set(annotation.ID)
    a = a.merge(annotation[['ID','PatientID','Disease','Inflammation','Region3']], on='ID', validate='one_to_one')
    a['CellType'] = ct
    assignments.append(a)
    coverage.append(dict(CellType=ct, AssignedBiopsies=len(a), Patients=a.PatientID.nunique(), UnassignedBiopsies=2025-len(a), States=len(states)))
    input_paths.append(p)
assigned = pd.concat(assignments, ignore_index=True)
pd.DataFrame(coverage).to_csv(out/'audit/discovery_parent_assignment_coverage.tsv', sep='\t', index=False)
members.to_csv(out/'audit/frozen_states.tsv', sep='\t', index=False)
print('Frozen:', len(annotation), annotation.PatientID.nunique(), len(members), members.Ecotype.nunique(), flush=True)

# Step2. 完整复核159项Fisher，独立重算Z公式及806项结构，不替换旧表
old_fisher = pd.read_csv(prior/'04_state_phenotypes/state_region_fisher.tsv', sep='\t')
checks = []
for row in old_fisher.itertuples():
    a = assigned[assigned.CellType.eq(row.CellType)]
    own = a.State.eq(row.CellState.removeprefix(row.CellType+'_'))
    target = a.Region3.eq(row.Region3)
    table = [[int((own&target).sum()),int((own&~target).sum())],[int((~own&target).sum()),int((~own&~target).sum())]]
    odds, p = stats.fisher_exact(table, alternative='greater')
    checks.append(dict(CellState=row.CellState, Region3=row.Region3, Observed=table[0][0], N=len(a), Patients=a.PatientID.nunique(), P=p, OldP=row.P, AbsPdiff=abs(p-row.P), CountMatch=table[0][0]==row.Ns, DenominatorMatch=len(a)==row.AssignedSampleDenominator))
checks = pd.DataFrame(checks)
checks['Q'] = stats.false_discovery_control(checks.P.to_numpy())
checks['OldQ'] = old_fisher.Q.to_numpy()
checks['AbsQdiff'] = (checks.Q-checks.OldQ).abs()
assert len(checks)==159 and checks.CountMatch.all() and checks.DenominatorMatch.all() and checks.AbsPdiff.max()<1e-12 and checks.AbsQdiff.max()<1e-12
checks.to_csv(out/'audit/fisher_159_independent_recalculation.tsv', sep='\t', index=False)
old = pd.read_csv(prior/'04_state_phenotypes/state_phenotype_enrichment.tsv', sep='\t')
old['RecalculatedZ'] = (old.Ns-old.NullMean)/old.NullSD.replace(0,np.nan)
old['AbsZdiff'] = (old.Z-old.RecalculatedZ).abs()
old['NullStatus'] = np.where(old.NullSD.gt(0),'VALID','ZERO_SD')
old.to_csv(out/'audit/old_806_formula_check.tsv', sep='\t', index=False)
old.groupby(['Dataset','Unit','Feature','Target'], dropna=False).size().rename('Rows').reset_index().to_csv(out/'audit/old_806_test_families.tsv',sep='\t',index=False)
print('Fisher159 verified; old806 max Z difference',old.AbsZdiff.max(),flush=True)

# Step3. 补充Discovery主要与预先指定分层比较：每种方法均1000次置换
# 原文式活检置换未控制重复患者；患者敏感性每个比较每父类型随机选一份合格活检。
# Restricted额外保持疾病/肠段边际结构；不把全局标签打乱称为患者控制。
tasks = [('All','Disease',{}),('All','Inflammation',{}),('All','Region3',{})]
tasks += [(d,'Inflammation',{'Disease':d}) for d in ['CD','UC']]
tasks += [(r,f,{'Region3':r}) for r in ['Colon','Ileum','Rectum'] for f in ['Disease','Inflammation']]
tasks += [('CD_Ileum','Inflammation',{'Disease':'CD','Region3':'Ileum'})]
tasks += [('Colon_Inflamed','Disease',{'Region3':'Colon','Inflammation':'Inflamed'})]
rows = []
selected = []
null_arrays = {}
for ct, parent in assigned.groupby('CellType', sort=False):
    states = members.loc[members.CellType.eq(ct),'State'].tolist()
    for scope, feature, filters in tasks:
        sub = parent.copy()
        for field, level in filters.items():
            sub = sub[sub[field].eq(level)]
        levels = sorted(sub[feature].dropna().unique())
        for mode in ['BiopsyNominal','OneBiopsyPerPatient','OneBiopsyPerPatientRestricted']:
            data = sub.copy()
            if mode!='BiopsyNominal':
                data = data.sort_values('ID').sample(frac=1,random_state=20261015).drop_duplicates('PatientID').sort_values('ID')
                for r in data.itertuples():
                    selected.append(dict(CellType=ct, Scope=scope, Feature=feature, Mode=mode, ID=r.ID, PatientID=r.PatientID))
                assert data.PatientID.is_unique
            codes = pd.Categorical(data.State,categories=states).codes
            y = pd.Categorical(data[feature],categories=levels).codes
            k = len(levels)
            obs = np.bincount(codes*k+y,minlength=len(states)*k).reshape(len(states),k)
            null = np.zeros((1000,len(states),k))
            block_fields = {'Disease':['Region3','Inflammation'],'Inflammation':['Disease','Region3'],'Region3':['Disease']}[feature]
            blocks = list(data.groupby(block_fields,sort=True).indices.values()) if mode.endswith('Restricted') else [np.arange(len(data))]
            rng = np.random.default_rng(20261015)
            for b in range(1000):
                shuffled = codes.copy()
                for idx in blocks:
                    shuffled[idx] = rng.permutation(codes[idx])
                null[b] = np.bincount(shuffled*k+y,minlength=len(states)*k).reshape(len(states),k)
            mu = null.mean(axis=0)
            sd = null.std(axis=0,ddof=1)
            key = '|'.join([ct,scope,feature,mode])
            null_arrays[key] = null.astype('int32')
            for i,state in enumerate(states):
                for j,level in enumerate(levels):
                    z = (obs[i,j]-mu[i,j])/sd[i,j] if sd[i,j]>0 and k>1 else np.nan
                    status = 'VALID' if np.isfinite(z) else ('SINGLE_PHENOTYPE_LEVEL' if k<2 else 'ZERO_NULL_SD')
                    rows.append(dict(Dataset='MSCCR',CellType=ct,CellState=ct+'_'+state,Scope=scope,Feature=feature,Target=level,Mode=mode,Observed=int(obs[i,j]),NullMean=mu[i,j],NullSD=sd[i,j],Z=z,P=stats.norm.sf(z),EmpiricalP=(1+int((null[:,i,j]>=obs[i,j]).sum()))/1001 if np.isfinite(z) else np.nan,PaperNominal=bool(z>1.96),Permutations=1000,EffectivePermutations=1000 if sd[i,j]>0 else 0,N=len(data),Patients=data.PatientID.nunique(),TargetN=int((y==j).sum()),TargetPatients=data.loc[data[feature].eq(level),'PatientID'].nunique(),Status=status,Shuffle='State labels within parent'+(' and '+','.join(block_fields) if mode.endswith('Restricted') else ''),BHFamily='All prespecified state-scope-target tests within mode'))
    print('Discovery completed',ct,flush=True)
result = pd.DataFrame(rows)
for mode, group in result.groupby('Mode'):
    valid = group.index[group.P.notna()]
    result.loc[valid,'Q'] = stats.false_discovery_control(result.loc[valid,'P'].to_numpy())
    result.loc[valid,'EmpiricalQ'] = stats.false_discovery_control(result.loc[valid,'EmpiricalP'].to_numpy())
result.to_csv(out/'statistics/discovery_state_phenotype_enrichment.tsv',sep='\t',index=False)
pd.DataFrame(selected).to_csv(out/'audit/one_biopsy_selected_units.tsv',sep='\t',index=False)
np.savez_compressed(out/'audit/discovery_permutation_nulls.npz',**null_arrays)

# Step4. 外部表型：真实IBD内部比较，细胞名义置换与供者等权比例分别统计
# 供者比例保留所有父类型细胞作为分母；Unassigned不改派给任何状态。
ext_rows = []
ext_coverage = []
ext_null = {}
for dataset in ['Colon_GSE214695','Ileum_IBDverse']:
    path = prior/'00_audit'/f'{dataset}_cells.tsv'
    cells = pd.read_csv(path,sep='\t')
    input_paths.append(path)
    ext_coverage.append(cells.groupby(['DiseaseGroup','InflammationGroup']).agg(Cells=('ID','size'),Donors=('Donor','nunique'),Samples=('Sample','nunique')).reset_index().assign(Dataset=dataset))
    contrasts = [('Disease','DiseaseGroup',['CD','UC'],cells.DiseaseGroup.isin(['CD','UC']))] if dataset.startswith('Colon') else [('Inflammation','InflammationGroup',['Inflamed','Non-inflamed'],cells.DiseaseGroup.eq('CD'))]
    for feature, field, levels, eligible in contrasts:
        for ct, st in members.groupby('CellType',sort=False):
            parent = cells[eligible & cells.CellType.eq(ct)].copy()
            states = st.State.tolist()
            for mode in ['CellNominal','DonorEqualFraction']:
                if mode=='CellNominal':
                    data = parent[parent.State.ne('Unassigned')].copy()
                    codes = pd.Categorical(data.State,categories=states).codes.copy()
                    assert (codes>=0).all()
                    y = pd.Categorical(data[field],categories=levels).codes
                    obs = np.bincount(codes*2+y,minlength=len(states)*2).reshape(len(states),2)
                    null = np.zeros((1000,len(states),2))
                    rng = np.random.default_rng(20261015)
                    for b in range(1000):
                        null[b] = np.bincount(rng.permutation(codes)*2+y,minlength=len(states)*2).reshape(len(states),2)
                    donors = data.Donor.nunique()
                    n = len(data)
                    target_n = [int((y==j).sum()) for j in range(2)]
                    target_donors = [data.loc[data[field].eq(v),'Donor'].nunique() for v in levels]
                else:
                    meta = parent[['Donor',field]].drop_duplicates()
                    assert not meta.Donor.duplicated().any()
                    denom = parent.groupby('Donor').size()
                    num = parent[parent.State.ne('Unassigned')].groupby(['Donor','State']).size().unstack(fill_value=0).reindex(index=denom.index,columns=states,fill_value=0)
                    fractions = num.div(denom,axis=0).to_numpy()
                    y = pd.Categorical(meta.set_index('Donor').reindex(denom.index)[field],categories=levels).codes.copy()
                    obs = np.column_stack([fractions[y==j].mean(axis=0) if (y==j).sum() else np.full(len(states),np.nan) for j in range(2)])
                    null = np.zeros((1000,len(states),2))
                    rng = np.random.default_rng(20261015)
                    for b in range(1000):
                        shuffled = rng.permutation(y)
                        for j in range(2):
                            null[b,:,j] = fractions[shuffled==j].mean(axis=0) if (shuffled==j).sum() else np.nan
                    donors = len(denom)
                    n = len(parent)
                    target_n = [int(parent[field].eq(v).sum()) for v in levels]
                    target_donors = [int((y==j).sum()) for j in range(2)]
                mu = null.mean(axis=0)
                sd = null.std(axis=0,ddof=1)
                ext_null['|'.join([dataset,ct,feature,mode])] = null
                for i,state in enumerate(states):
                    for j,target in enumerate(levels):
                        status = 'VALID'
                        if n==0:
                            status = 'NO_PARENT_CELLS'
                        elif min(target_donors)<2:
                            status = 'INSUFFICIENT_DONORS'
                        elif not np.isfinite(sd[i,j]) or sd[i,j]==0:
                            status = 'ZERO_NULL_SD_OR_ABSENT_STATE'
                        z = (obs[i,j]-mu[i,j])/sd[i,j] if status=='VALID' else np.nan
                        ext_rows.append(dict(Dataset=dataset,CellType=ct,CellState=ct+'_'+state,Feature=feature,Target=target,Mode=mode,Observed=obs[i,j],NullMean=mu[i,j],NullSD=sd[i,j],Z=z,P=stats.norm.sf(z),EmpiricalP=(1+int((null[:,i,j]>=obs[i,j]-1e-15).sum()))/1001 if np.isfinite(z) else np.nan,PaperNominal=bool(z>1.96),Permutations=1000,EffectivePermutations=1000 if status=='VALID' else 0,N=n,Donors=donors,TargetCells=target_n[j],TargetDonors=target_donors[j],MinGroupDonors=min(target_donors),Status=status,Denominator='Assigned parent cells' if mode=='CellNominal' else 'All parent cells per donor; equal donor weights',BHFamily='All 53 states x 2 targets within dataset and mode'))
    print('External completed',dataset,len(cells),cells.Donor.nunique(),flush=True)
external = pd.DataFrame(ext_rows)
for (dataset,mode), group in external.groupby(['Dataset','Mode']):
    valid = group.index[group.P.notna()]
    external.loc[valid,'Q'] = stats.false_discovery_control(external.loc[valid,'P'].to_numpy())
    external.loc[valid,'EmpiricalQ'] = stats.false_discovery_control(external.loc[valid,'EmpiricalP'].to_numpy())
external.to_csv(out/'statistics/external_scRNA_phenotype_enrichment.tsv',sep='\t',index=False)
pd.concat(ext_coverage).to_csv(out/'audit/external_phenotype_coverage.tsv',sep='\t',index=False)
np.savez_compressed(out/'audit/external_permutation_nulls.npz',**ext_null)

# Step5. 按完全相同对比连接Discovery与外部，未评价/阴性完整保留
concordance = []
for dataset,scope,feature in [('Colon_GSE214695','Colon_Inflamed','Disease'),('Ileum_IBDverse','CD_Ileum','Inflammation')]:
    for discovery_mode, external_mode in [('BiopsyNominal','CellNominal'),('OneBiopsyPerPatient','DonorEqualFraction')]:
        a = result[result.Scope.eq(scope)&result.Feature.eq(feature)&result.Mode.eq(discovery_mode)]
        b = external[external.Dataset.eq(dataset)&external.Mode.eq(external_mode)]
        merged = a.merge(b,on=['CellState','Feature','Target'],suffixes=('_Discovery','_External'),validate='one_to_one')
        merged['Rule'] = discovery_mode+' vs '+external_mode
        merged['Evaluable'] = merged.Status_Discovery.eq('VALID')&merged.Status_External.eq('VALID')
        merged['DirectionAgree'] = (np.sign(merged.Z_Discovery)==np.sign(merged.Z_External))&merged.Evaluable
        merged['BothPaperNominal'] = merged.PaperNominal_Discovery&merged.PaperNominal_External&merged.Evaluable
        merged['BothEmpiricalBH'] = merged.EmpiricalQ_Discovery.lt(.05)&merged.EmpiricalQ_External.lt(.05)&merged.Evaluable
        concordance.append(merged)
concordance = pd.concat(concordance,ignore_index=True)
concordance.to_csv(out/'statistics/discovery_external_concordance_all_states.tsv',sep='\t',index=False)
summaries = []
for (dataset,rule,target),g in concordance.groupby(['Dataset_External','Rule','Target']):
    v = g[g.Evaluable]
    table = pd.crosstab(v.PaperNominal_Discovery,v.PaperNominal_External).reindex(index=[False,True],columns=[False,True],fill_value=0)
    odds,p = stats.fisher_exact(table.to_numpy(),alternative='greater')
    summaries.append(dict(Dataset=dataset,Rule=rule,Target=target,FrozenStates=53,EvaluableStates=len(v),NotEvaluableStates=53-len(v),DiscoveryNominal=int(v.PaperNominal_Discovery.sum()),ExternalNominal=int(v.PaperNominal_External.sum()),BothNominal=int(v.BothPaperNominal.sum()),DirectionAgree=int(v.DirectionAgree.sum()),DirectionDisagree=int((~v.DirectionAgree).sum()),BothEmpiricalBH=int(v.BothEmpiricalBH.sum()),FisherConcordanceP=p,FisherConcordanceOdds=odds,ConcordanceDenominator='Discovery nominal enriched among mutually evaluable states',ConcordanceFraction=float(v.BothPaperNominal.sum()/v.PaperNominal_Discovery.sum()) if v.PaperNominal_Discovery.any() else np.nan))
summary = pd.DataFrame(summaries)
summary['FisherConcordanceQ'] = stats.false_discovery_control(summary.FisherConcordanceP.to_numpy())
summary.to_csv(out/'statistics/concordance_summary.tsv',sep='\t',index=False)

# Step6. 官方显著分配与argmax分开，显示Unassigned及两种分母
ep = frozen/'Ecotypes/ecotype_assignment.txt'
ea = pd.read_csv(ep,sep='\t')
assert len(ea)==1712 and ea.ID.is_unique and ea.Ecotype.nunique()==7
cohort = annotation[['ID','PatientID','Disease','Inflammation','Region3']].merge(ea[['ID','Ecotype']],on='ID',how='left',validate='one_to_one')
cohort['Assignment'] = cohort.Ecotype.fillna('Unassigned')
distribution = []
for feature in ['Disease','Inflammation','Region3']:
    for target,g in cohort.groupby(feature):
        for ecotype in ['E1','E2','E3','E4','E5','E6','E7','Unassigned']:
            own = g[g.Assignment.eq(ecotype)]
            distribution.append(dict(Feature=feature,Target=target,Assignment=ecotype,Count=len(own),Patients=own.PatientID.nunique(),PhenotypeBiopsies=len(g),PhenotypePatients=g.PatientID.nunique(),PercentWithinPhenotype=100*len(own)/len(g),AssignmentBiopsies=int(cohort.Assignment.eq(ecotype).sum()),PercentPhenotypeWithinAssignment=100*len(own)/cohort.Assignment.eq(ecotype).sum(),Definition='Frozen official significant assignment; 313 unassigned retained'))
pd.DataFrame(distribution).to_csv(out/'statistics/ecotype_phenotype_distribution.tsv',sep='\t',index=False)
input_paths += [ep, frozen/'Ecotypes/ecotypes.txt']
hashes = [dict(Path=str(p),Bytes=p.stat().st_size,SHA256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in input_paths]
pd.DataFrame(hashes).to_csv(out/'audit/input_sha256_after.tsv',sep='\t',index=False)
result.groupby(['Mode','Scope','Feature']).agg(Tests=('Z','size'),Evaluable=('Z','count'),NominalPositive=('PaperNominal','sum'),EmpiricalBH=('EmpiricalQ',lambda x:int((x<.05).sum()))).reset_index().to_csv(out/'statistics/discovery_summary.tsv',sep='\t',index=False)
print(summary.to_string(index=False),flush=True)
(out/'audit/completion.json').write_text(json.dumps(dict(DiscoveryRows=len(result),ExternalRows=len(external),ConcordanceRows=len(concordance),FisherRows=len(checks),Python=sys.version,Platform=platform.platform(),SourceChanged=False),indent=2))
print('SUPPLEMENT_COMPLETE',flush=True)
