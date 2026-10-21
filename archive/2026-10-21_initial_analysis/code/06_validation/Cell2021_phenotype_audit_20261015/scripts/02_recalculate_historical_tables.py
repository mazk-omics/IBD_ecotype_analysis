# 本机独立重算806行旧置换及350/385项统计，记录差异，不修改任何历史文件。
# 在msccr_clinical_audit_20261015环境运行；输入均为本机可信归档。
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy import stats

# Step1. 按原始随机数顺序重建806项零分布，检查并列剔除及缺失状态
root = Path('D:/IBD_ecotype/IBD_ecotype_analysis')
out = Path(__file__).resolve().parents[1]
prior = root/'06_validation/paper_aligned_recovery_characterization_20261012'
old = pd.read_csv(prior/'04_state_phenotypes/state_phenotype_enrichment.tsv',sep='\t')
rng = np.random.default_rng(20261012)
rows = []
coverage = []
for dataset in ['Colon_GSE214695','Ileum_IBDverse']:
    cells = pd.read_csv(prior/'00_audit'/f'{dataset}_cells.tsv',sep='\t')
    for ct,parent in cells[cells.State.ne('Unassigned')].groupby('CellType'):
        for unit in ['Cell','BiologicalSampleDominant']:
            if unit=='Cell':
                data = parent[['State','DiseaseGroup','InflammationGroup','Sample','Donor']].copy()
                ties = 0
            else:
                counts = parent.groupby(['Sample','State']).size().unstack(fill_value=0)
                tied = counts.eq(counts.max(axis=1),axis=0).sum(axis=1)>1
                ties = int(tied.sum())
                counts = counts.loc[~tied]
                data = parent[['Sample','Donor','DiseaseGroup','InflammationGroup']].drop_duplicates().set_index('Sample').loc[counts.index].reset_index()
                data['State'] = counts.idxmax(axis=1).to_numpy()
            coverage.append(dict(Dataset=dataset,CellType=ct,Unit=unit,N=len(data),Donors=data.Donor.nunique(),TiedSamplesExcluded=ties,RepresentedStates=data.State.nunique()))
            labels = data.State.to_numpy(copy=True)
            states = sorted(set(labels))
            for feature in ['DiseaseGroup','InflammationGroup']:
                levels = sorted(data[feature].dropna().unique())
                y = pd.Categorical(data[feature],categories=levels).codes
                k = len(levels)
                null = np.zeros((1000,len(levels),len(states)))
                for b in range(1000):
                    shuffled = rng.permutation(labels)
                    sc = pd.Categorical(shuffled,categories=states).codes
                    null[b] = np.bincount(y*len(states)+sc,minlength=len(levels)*len(states)).reshape(len(levels),len(states))
                for j,target in enumerate(levels):
                    for i,state in enumerate(states):
                        ns = int(((labels==state)&data[feature].eq(target).to_numpy()).sum())
                        mu = null[:,j,i].mean()
                        sd = null[:,j,i].std(ddof=1)
                        rows.append(dict(Dataset=dataset,CellState=ct+'_'+state,Unit=unit,Feature=feature,Target=target,Ns_new=ns,NullMean_new=mu,NullSD_new=sd,Z_new=(ns-mu)/sd if sd>0 else np.nan))
    print('806 replay',dataset,flush=True)
replay = pd.DataFrame(rows).merge(old,on=['Dataset','CellState','Unit','Feature','Target'],validate='one_to_one')
for field in ['Ns','NullMean','NullSD','Z']:
    replay[field+'_absdiff'] = (replay[field]-replay[field+'_new']).abs()
    np.testing.assert_allclose(replay[field],replay[field+'_new'],rtol=1e-10,atol=1e-10,equal_nan=True)
assert len(replay)==806
replay.to_csv(out/'audit/old_806_full_permutation_replay.tsv',sep='\t',index=False)
pd.DataFrame(coverage).to_csv(out/'audit/old_806_unit_coverage.tsv',sep='\t',index=False)

# Step2. 重算350项Hallmark组间比较、分母和BH；dominant是2025样本argmax
scores = pd.read_csv(prior/'06_ecotype_features/hallmark_sample_scores.tsv',sep='\t',index_col=0)
metadata = pd.read_csv(prior/'06_ecotype_features/ecotype_dominant_assignments.tsv',sep='\t',index_col=0)
old350 = pd.read_csv(prior/'06_ecotype_features/ecotype_continuous_features.tsv',sep='\t')
checks = []
for row in old350.itertuples():
    own = metadata.DominantEcotype.eq(row.Ecotype)
    other = metadata.DominantEcotype.notna() & ~own
    a,b = scores.loc[own,row.Feature],scores.loc[other,row.Feature]
    p = stats.mannwhitneyu(a,b,alternative='two-sided',method='asymptotic').pvalue
    checks.append(dict(Ecotype=row.Ecotype,Feature=row.Feature,P_new=p,P_old=row.P,MeanDifference_new=a.mean()-b.mean(),MeanDifference_old=row.MeanDifference,N_Ecotype=len(a),N_Other=len(b),N_match=len(a)==row.N_Ecotype and len(b)==row.N_OtherEcotypes))
checks = pd.DataFrame(checks)
checks['Q_new'] = stats.false_discovery_control(checks.P_new.to_numpy())
checks['Q_old'] = old350.Q.to_numpy()
assert checks.N_match.all()
np.testing.assert_allclose(checks.P_new,checks.P_old,atol=1e-12)
np.testing.assert_allclose(checks.Q_new,checks.Q_old,atol=1e-12)
checks.to_csv(out/'audit/old_350_hallmark_recalculation.tsv',sep='\t',index=False)

# Step3. 385项确认为350 Hallmark、14疾病/炎症、21肠段连续丰度比较
old385 = pd.read_csv(root/'06_validation/paper_figures_20261012/source_data/Fig5A_all_385_feature_tests.tsv',sep='\t')
extra = []
for row in old385.itertuples():
    if row.Feature.startswith('HALLMARK_'):
        p = float(checks.loc[checks.Ecotype.eq(row.Ecotype)&checks.Feature.eq(row.Feature),'P_new'].iloc[0])
    elif row.Feature.startswith('Region3_'):
        own = metadata.Region3.eq(row.Feature.removeprefix('Region3_'))
        p = stats.mannwhitneyu(metadata.loc[own,row.Ecotype],metadata.loc[~own,row.Ecotype],alternative='two-sided',method='asymptotic').pvalue
    else:
        target = 'CD' if row.Feature=='Disease' else 'Inflamed'
        other = 'UC' if row.Feature=='Disease' else 'Non-inflamed'
        p = stats.mannwhitneyu(metadata.loc[metadata[row.Feature].eq(target),row.Ecotype],metadata.loc[metadata[row.Feature].eq(other),row.Ecotype],alternative='two-sided',method='asymptotic').pvalue
    extra.append(dict(Ecotype=row.Ecotype,Feature=row.Feature,P_new=p,P_old=row.P,Q_old=row.Q))
extra = pd.DataFrame(extra)
extra['Q_new'] = stats.false_discovery_control(extra.P_new.to_numpy())
np.testing.assert_allclose(extra.P_new,extra.P_old,atol=1e-12)
np.testing.assert_allclose(extra.Q_new,extra.Q_old,atol=1e-12)
extra.to_csv(out/'audit/old_385_complete_recalculation.tsv',sep='\t',index=False)
summary = dict(Old806Rows=len(replay),Old806MaxNullMeanDiff=float(replay.NullMean_absdiff.max()),Old806MaxNullSDDiff=float(replay.NullSD_absdiff.max()),Old806MaxZDiff=float(replay.Z_absdiff.max()),Old350MaxPdiff=float((checks.P_new-checks.P_old).abs().max()),Old385MaxQdiff=float((extra.Q_new-extra.Q_old).abs().max()))
(out/'audit/historical_recalculation_summary.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2),flush=True)
