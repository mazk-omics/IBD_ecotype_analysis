# 在服务器 ibd_outcomes_20261015 环境运行，核对外部结局与冻结评分输入。
# 读取既有 GEO、metadata 和表达文件；仅写本轮目录，原始标签保留。
from pathlib import Path
import gzip, hashlib, json, urllib.request, re
import numpy as np
import pandas as pd

# Step1. 固定目录与模型，保存所有实际使用的输入哈希。
root = Path('/home/mazekai/IBD_EcoTyper')
out = root / '06_validation/external_clinical_outcomes_20261015'
model = root / 'software/ecotyper/EcoTyper/MSCCR_CORE12_CTSPEC_REGION_IBD/Cell_type_specific_genes'
ecotypes = pd.read_csv(model / 'Ecotypes/discovery/ecotypes.txt', sep='\t')
ranks = pd.read_csv(model / 'Analysis/rank_selection/rank_data.txt', sep='\t')
assert len(ecotypes) == 53 and ecotypes.ID.is_unique and ecotypes.Ecotype.nunique() == 7 and len(ranks) == 12
ecotypes.to_csv(out / 'audit/frozen_state_ecotype_membership.tsv', sep='\t', index=False)
inputs = [model / 'Ecotypes/discovery/ecotypes.txt', model / 'Analysis/rank_selection/rank_data.txt']
for _, row in ranks.iterrows():
    d = model / 'Cell_States/discovery' / row.iloc[0] / str(row.iloc[1])
    inputs += [d / 'W.txt', d / 'mapping_to_initial_states.txt', d / 'initial_gene_info.txt']

# Step2. 从原始 SOFT 直接恢复字段，避免沿用 response 的错误泛称。
soft = root / '06_validation/02_bulk_treatment/GSE72819/00_provenance/geo_downloads/GSE72819_family.soft.gz'
records = []
for line in gzip.open(soft, 'rt'):
    if line.startswith('^SAMPLE = '):
        records.append({'GSM': line.strip().split(' = ')[1]})
    if line.startswith('!Sample_characteristics_ch1 = '):
        k, v = line.strip().split(' = ', 1)[1].split(': ', 1)
        records[-1][k] = v
raw = pd.DataFrame(records)
raw.to_csv(out / '01_GSE72819_treatment_outcome/raw_GEO_characteristics.tsv', sep='\t', index=False)
path = root / '06_validation/02_bulk_treatment/GSE72819/prepared/04_sample_metadata/GSE72819_canonical_sample_metadata_73.tsv'
uc = pd.read_csv(path, sep='\t')
uc = uc.merge(raw[['GSM', 'remission at week 10', 'inventory sample name']], on='GSM', validate='one_to_one')
assert uc.PatientID.is_unique and (uc.PatientID == uc['inventory sample name']).all()
assert uc.RawResponseWeek10.fillna('NA').equals(uc['remission at week 10'].replace(['NA', 'N/A', 'Unknown', 'Not Available', '-'], np.nan).fillna('NA'))
uc['Outcome'] = uc['remission at week 10'].map({'Remitter': 1, 'Non-remitter': 0})
uc['Study'] = 'EUCALYPTUS'
uc['GEOAccession'] = 'GSE72819'
uc['Endpoint'] = 'ETR_W10_clinical_remission'
uc['EndpointWeek'] = 10
uc['EndpointDefinition'] = 'GEO remission at week 10; trial definition total Mayo <=2, no subscore >1 (trial-level definition; sample-level Mayo unavailable)'
uc['TreatmentStart'] = 'After baseline biopsy; exact date unavailable'
uc['ExpressionScale'] = 'Deposited RPKM; official Recovery log2(x+1), then gene-wise standardization'
uc['PriorAntiTNFBinary'] = uc.PriorAntiTNF.map({'Prior anti-TNF': 1, 'Naive': 0})
uc['ExclusionReason'] = np.where(uc.Outcome.isna(), 'Missing original remission label', '')
uc['PatientIDProvenance'] = 'GEO inventory sample name; unique across 73 biopsies'
uc.to_csv(out / '01_GSE72819_treatment_outcome/cohort_endpoint_manifest.tsv', sep='\t', index=False)
inputs += [path, soft]
print('GSE72819', len(uc), 'patients', uc.PatientID.nunique(), 'outcomes', uc.Outcome.value_counts(dropna=False).to_dict(), flush=True)

# Step3. RISK 采用公开三年固定终点，不制造生存时间或患者跨队列标识。
path = root / '06_validation/01_bulk_progression/GSE93624/prepared/04_sample_metadata/GSE93624_canonical_sample_metadata.tsv'
risk = pd.read_csv(path, sep='\t')
assert risk.ExpressionSampleID.is_unique and risk.GSM.is_unique
risk['PatientID'] = risk.ExpressionSampleID
risk['PatientIDProvenance'] = 'Anonymous study individual per paper; public ID is SampleNumber, not a verified cross-study patient ID'
risk['SampleID'] = risk.ExpressionSampleID
risk['Disease'] = risk.Diagnosis
risk['Outcome'] = risk.Progression3yr.map({'Yes': 1, 'No': 0})
risk['SexBinary'] = risk.Sex.map({'Male': 1, 'Female': 0})
risk['Study'] = 'RISK'
risk['GEOAccession'] = 'GSE93624'
risk['Timepoint'] = 'Diagnosis; treatment naive'
risk['Endpoint'] = 'RISK_3yr_B2_B3'
risk['EndpointWeek'] = 156
risk['EndpointDefinition'] = 'B2/B3 within 3 years; paper reports B1 through first 90 days; 17 B2 plus 10 B3; individual dates/censoring unavailable'
risk['Treatment'] = 'Treatment naive at biopsy; later treatment unavailable and not modeled'
risk['ExpressionScale'] = 'GEO deposited TMM-normalized processed expression; exact deposited transformation unresolved'
risk['InputGate'] = 'Conditional processed-scale projection: deposited values preserved; exact prior transformation unresolved'
risk['ExclusionReason'] = np.where(risk.Diagnosis.eq('CD'), '', 'Non-IBD control')
risk.to_csv(out / '02_RISK_disease_progression/cohort_endpoint_manifest.tsv', sep='\t', index=False)
inputs += [path]
print('RISK', risk.Diagnosis.value_counts().to_dict(), risk.Outcome.value_counts(dropna=False).to_dict(), flush=True)

# Step4. VDZ 以基线诱导方案界定暴露，避免用后续治疗决定基线分组。
path = root / '06_validation/02_bulk_treatment/GSE73661/prepared/04_canonical_metadata/GSE73661_canonical_metadata_178_samples.tsv'
long = pd.read_csv(path, sep='\t')
long['MayoEndoscopicSubscore'] = pd.to_numeric(long.MayoEndoscopicSubscore, errors='coerce')
long['HealingBinary'] = np.where(long.MayoEndoscopicSubscore.notna(), long.MayoEndoscopicSubscore.le(1).astype(float), np.nan)
assert long.GSM.is_unique
baseline = long.loc[long.Timepoint.eq('W0')].copy()
assert baseline.PatientKey.is_unique
print('Baseline schedules', baseline.groupby(['Branch', 'TreatmentSchedule'], dropna=False).size().to_dict(), flush=True)
rows = []
for drug, visit in [('IFX', 'W4_6'), ('VDZ_TRIAL', 'W6'), ('VDZ_TRIAL', 'W12'), ('VDZ_TRIAL', 'W52')]:
    b = baseline.loc[baseline.Branch.eq(drug)].copy()
    follow = long.loc[long.Branch.eq(drug) & long.Timepoint.eq(visit), ['PatientKey', 'GSM', 'HealingBinary', 'MayoEndoscopicSubscore', 'TreatmentSchedule']].copy()
    assert follow.PatientKey.is_unique
    follow.columns = ['PatientKey', 'FollowupGSM', 'Outcome', 'FollowupMayo', 'FollowupSchedule']
    b = b.merge(follow, on='PatientKey', how='left', validate='one_to_one')
    b['BaselineInduction'] = b.TreatmentSchedule.str.split('_').str[0]
    b['BaselineDrugEligible'] = b.BaselineInduction.eq('ifx') if drug == 'IFX' else b.BaselineInduction.str.startswith('vdz')
    b['PatientID'] = b.PatientKey
    b['SampleID'] = b.GSM
    b['Study'] = 'GSE73661'
    b['GEOAccession'] = 'GSE73661'
    b['Endpoint'] = ('IFX' if drug == 'IFX' else 'VDZ') + '_' + visit + '_endoscopic_healing'
    b['EndpointWeek'] = {'W4_6': 6, 'W6': 6, 'W12': 12, 'W52': 52}[visit]
    b['EndpointDefinition'] = 'Follow-up Mayo endoscopic subscore 0 or 1; W4_6 window retained for IFX'
    b['Treatment'] = np.where(b.BaselineDrugEligible, 'IFX' if drug == 'IFX' else 'VDZ baseline induction', 'Placebo baseline induction; excluded from VDZ primary')
    b['Tissue'] = 'Colonic mucosal biopsy'
    b['TreatmentStart'] = 'After W0; date unavailable'
    b['ExpressionScale'] = 'GEO log2 RMA; used as-is under official Recovery max<=50 branch; no inverse transform'
    b['ExclusionReason'] = np.select([~b.BaselineDrugEligible, b.Outcome.isna()], ['Not active-drug induction at W0', 'No biopsy-derived outcome at this visit'], default='')
    rows.append(b)
manifest = pd.concat(rows, ignore_index=True)
manifest.to_csv(out / '03_GSE73661_mucosal_healing/cohort_endpoint_manifest.tsv', sep='\t', index=False)
long.to_csv(out / '03_GSE73661_mucosal_healing/all_samples_longitudinal_mapping.tsv', sep='\t', index=False)
print(manifest.groupby('Endpoint').apply(lambda d: pd.Series({'Baseline': len(d), 'DrugBaseline': int(d.BaselineDrugEligible.sum()), 'Evaluable': int(d.ExclusionReason.eq('').sum()), 'Events': int(d.loc[d.ExclusionReason.eq(''), 'Outcome'].sum())})), flush=True)
inputs += [path]

# Step5. 表达逐矩阵审计，明确识别样本列，不将基因注释列混入数值矩阵。
qc = []
for accession, rel in [('GSE72819', '02_bulk_treatment/GSE72819/prepared/06_ecotyper_recovery_input/GSE72819_bulk_recovery_expression.txt'), ('GSE93624', '01_bulk_progression/GSE93624/prepared/06_ecotyper_recovery_input/GSE93624_EcoTyper_recovery_expression.txt'), ('GSE57945', '01_bulk_progression/GSE57945/prepared/06_ecotyper_recovery_input/GSE57945_EcoTyper_recovery_expression.txt')]:
    path = root / '06_validation' / rel
    x = pd.read_csv(path, sep='\t', index_col=0)
    assert x.index.is_unique and x.columns.is_unique and np.isfinite(x.to_numpy()).all()
    qc.append({'Accession': accession, 'Genes': len(x), 'Samples': x.shape[1], 'Min': x.min().min(), 'Max': x.max().max(), 'NegativeN': int(x.lt(0).sum().sum()), 'DuplicateExpressionProfiles': int(x.T.duplicated().sum()), 'OfficialAutoLog2Plus1': bool(x.ge(0).all().all() and x.max().max() > 50)})
    inputs.append(path)
pd.DataFrame(qc).to_csv(out / 'audit/expression_QC.tsv', sep='\t', index=False)
cross = root / '06_validation/01_bulk_progression/99_cross_cohort/GSE57945_vs_GSE93624'
for f in cross.glob('*.tsv'):
    if 'crosswalk' in f.name or 'summary' in f.name or 'candidate' in f.name:
        (out / '02_RISK_disease_progression' / f.name).write_bytes(f.read_bytes())
cross57945 = pd.read_csv(root / '06_validation/01_bulk_progression/GSE57945/prepared/04_sample_metadata/GSE57945_canonical_sample_metadata_322.tsv', sep='\t')
pd.DataFrame([{'Metric': 'Public patient ID intersection', 'Value': 'Not identifiable: GSE93624 Sample1..245 vs GSE57945 CCFA_Risk_001..322'}, {'Metric': 'Expression fingerprint candidates', 'Value': '152 existing candidates; not proof of individual identity; no outcome-label transfer'}, {'Metric': 'Independent validation status', 'Value': 'Not established; never split these studies into training and validation'}]).to_csv(out / '02_RISK_disease_progression/cohort_overlap_audit.tsv', sep='\t', index=False)

# Step6. 保存可用协变量与缺失分母，网络核验失败也留痕。
for d, folder, fields in [(uc, '01_GSE72819_treatment_outcome', ['PriorAntiTNF', 'DoseArm', 'TissueDiagnosis']), (risk.loc[risk.Diagnosis.eq('CD')], '02_RISK_disease_progression', ['AgeAtDiagnosis', 'Sex', 'Ancestry', 'ParisAge']), (baseline, '03_GSE73661_mucosal_healing', ['MayoEndoscopicSubscore', 'TreatmentSchedule'])]:
    pd.DataFrame([{'Field': c, 'N': len(d), 'Observed': int(d[c].notna().sum()), 'Missing': int(d[c].isna().sum()), 'Levels': json.dumps(d[c].value_counts(dropna=False).to_dict(), ensure_ascii=False)} for c in fields]).to_csv(out / folder / 'clinical_covariate_missingness.tsv', sep='\t', index=False)
network = []
for acc in ['GSE72819', 'GSE93624', 'GSE57945', 'GSE73661']:
    url = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=' + acc
    try:
        b = urllib.request.urlopen(url, timeout=25).read()
        (out / 'audit' / (acc + '_live_GEO.html')).write_bytes(b)
        network.append({'URL': url, 'Status': 'Retrieved', 'Bytes': len(b)})
    except Exception as e:
        network.append({'URL': url, 'Status': repr(e), 'Bytes': 0})
pd.DataFrame(network).to_csv(out / 'audit/source_retrieval.tsv', sep='\t', index=False)
pd.DataFrame([{'Path': str(p), 'Bytes': p.stat().st_size, 'SHA256': hashlib.sha256(p.read_bytes()).hexdigest()} for p in inputs]).to_csv(out / 'audit/immutable_input_hashes_before.tsv', sep='\t', index=False)
print('AUDIT_COMPLETE', flush=True)
