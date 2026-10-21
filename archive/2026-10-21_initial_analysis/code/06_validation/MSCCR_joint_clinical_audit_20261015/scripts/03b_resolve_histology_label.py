# 原始 historemiss 与论文阈值冲突，保留原字段并显式派生零组织学活动变量。
# 仅修改本次整合副本；初版统计归档，不覆盖原始 metadata 和冻结结果。
from pathlib import Path
import pandas as pd
import numpy as np
import shutil

# Step1. 保留初版统计以审计，主分析排除定义未核实的 historemiss
out = Path('D:/IBD_ecotype/IBD_ecotype_analysis/06_validation/MSCCR_joint_clinical_audit_20261015')
archive = out/'statistics/superseded_unverified_historemiss'
archive.mkdir(exist_ok=True)
for path in (out/'statistics').glob('*.tsv'):
    if not (archive/path.name).exists():
        shutil.copy2(path,archive/path.name)
data = pd.read_csv(out/'source_data/analysis_metadata_2025.tsv',sep='\t',index_col=0)
score = data.Max_GHAS.where(data.Disease.eq('CD'),data.Max_Nancy)
data['Derived_histology_zero'] = np.where(score.notna(),score.eq(0).astype(float),np.nan)
data.to_csv(out/'source_data/analysis_metadata_2025.tsv',sep='\t')

# Step2. 字典记录派生规则、缺失、冲突和报告边界
dictionary = pd.read_csv(out/'source_data/clinical_variable_dictionary.tsv',sep='\t')
dictionary.loc[dictionary.Feature.eq('Histological_remission'),'Definition'] = 'Original binary GEO historemiss; 189 biopsy labels disagree with paper-defined max GHAS=0 CD / max Nancy=0 UC; exact source coding unresolved'
dictionary.loc[dictionary.Feature.eq('Histological_remission'),'Use'] = 'WITHHELD from primary inference; audit only'
derived = dictionary.loc[dictionary.Feature.eq('Histological_remission')].copy()
derived['Field'] = 'Derived from max_ghas_sum7 CD / max_nancy UC'
derived['Feature'] = 'Derived_histology_zero'
derived['Definition'] = 'Explicit derived indicator: maximum GHAS=0 in CD, maximum Nancy=0 in UC; 1 zero activity, 0 nonzero; not original historemiss'
derived['Source'] = 'Derived in this audit from original GEO histology scores; paper threshold https://pmc.ncbi.nlm.nih.gov/articles/PMC10014487/'
derived['Use'] = 'cross-sectional association; explicit derived variable'
dictionary = pd.concat([dictionary.loc[~dictionary.Feature.eq('Derived_histology_zero')],derived],ignore_index=True)
dictionary.to_csv(out/'source_data/clinical_variable_dictionary.tsv',sep='\t',index=False)
coverage = []
for (disease,region),group in data.groupby(['Disease','Region3']):
    observed = group.loc[group.Derived_histology_zero.notna()]
    coverage.append(dict(Field='Derived score zero',Feature='Derived_histology_zero',Disease=disease,Region3=region,NonmissingBiopsies=len(observed),NonmissingPatients=observed.PatientID.nunique(),UsedForDiseaseSpecificAnalysis=True))
prior = pd.read_csv(out/'source_data/clinical_coverage_by_disease_region.tsv',sep='\t')
pd.concat([prior.loc[~prior.Feature.eq('Derived_histology_zero')],pd.DataFrame(coverage)]).to_csv(out/'source_data/clinical_coverage_by_disease_region.tsv',sep='\t',index=False)
print('Derived zero histology',data.Derived_histology_zero.value_counts(dropna=False).to_dict(),flush=True)
