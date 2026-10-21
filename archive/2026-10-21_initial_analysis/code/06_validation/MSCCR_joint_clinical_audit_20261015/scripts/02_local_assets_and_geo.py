# 在本机核查 MSCCR 文件、GEO 原始字段与冻结结果；不修改输入。
# 输出到本次独立目录，GSM 作为 Excel 与 GEO 的连接键，避免患者编号列误作样本 ID。
from pathlib import Path
import pandas as pd
import numpy as np
import gzip, csv, hashlib, json, os

# Step1. 保存本地资产清单，目录范围固定为项目、既有审计及 MSCCR 数据
workspace = Path('D:/IBD_ecotype')
root = workspace / 'IBD_ecotype_analysis'
out = root / '06_validation/MSCCR_joint_clinical_audit_20261015'
rows = []
for base in [root, workspace/'data_bulk', workspace/'Profile']:
    for directory, dirs, files in os.walk(base):
        dirs[:] = [x for x in dirs if x not in ['.git', '__pycache__', out.name]]
        for name in files:
            p = Path(directory)/name
            text = str(p).lower()
            if not any(x in text for x in ['msccr','193677','ctspec_region_ibd','paper_aligned_recovery_characterization_20261012','paper_figures_20261012']):
                continue
            s = p.stat()
            phase = 'FrozenDiscovery' if '/ctspec_region_ibd/' in p.as_posix() else 'Characterization' if 'paper_' in text else 'Metadata' if any(x in text for x in ['metadata','annotation','series_matrix']) else 'DerivedOrAudit'
            row = dict(path=str(p), bytes=s.st_size, modified=s.st_mtime, suffix=p.suffix, phase=phase)
            if s.st_size < 3000000 and (name in ['state_abundances.txt','ecotype_abundance.txt','ecotypes.txt','ecotype_assignment.txt'] or 'annotation.tsv' in name or p.suffix=='.xlsx'):
                row['sha256'] = hashlib.sha256(p.read_bytes()).hexdigest()
            rows.append(row)
pd.DataFrame(rows).to_csv(out/'audit/local_inventory.tsv', sep='\t', index=False)
print('Local assets', len(rows), flush=True)

# Step2. GEO series matrix 保留每一实际字段，characteristics 同名字段冲突立即停止
geo_path = workspace/'Profile/audit/2026-10-06_ctspec_region_ibd/GSE193677_series_matrix.txt.gz'
fields = {}
with gzip.open(geo_path,'rt') as handle:
    for line in handle:
        if not line.startswith('!Sample_'):
            continue
        values = next(csv.reader([line.rstrip('\n')], delimiter='\t'))
        key = values[0].removeprefix('!Sample_')
        if key=='characteristics_ch1':
            labels = [v.split(': ',1)[0] for v in values[1:]]
            assert len(set(labels))==1
            key = labels[0]
            values = [key]+[v.split(': ',1)[1] if ': ' in v else '' for v in values[1:]]
        if key in fields:
            if fields[key]==values[1:]:
                continue
            key = key+'__'+str(sum(k.startswith(key) for k in fields)+1)
        fields[key]=values[1:]
geo = pd.DataFrame(fields).replace({'NA':np.nan,'':np.nan,'NaN':np.nan,'N/A':np.nan})
geo['SampleID'] = geo.title.str.extract(r'(MSCCR_reGRID_\d+_Biopsy_\d+)')[0]
geo['GEO_PatientID'] = geo.SampleID.str.extract(r'MSCCR_(reGRID_\d+)_')[0]
geo.to_csv(out/'source_data/GEO_metadata_2490.tsv', sep='\t', index=False)
assert len(geo)==2490 and geo.geo_accession.is_unique and geo.title.is_unique
print('GEO fields',geo.columns.tolist(), flush=True)
print(geo.head(2).to_string(index=False), flush=True)

# Step3. Excel 只通过唯一 GSM 匹配；逐字段记录而不采用错误的患者编号作 biopsy ID
xlsx = workspace/'data_bulk/discovery_cohort/metadata/MSCCR(GSE193677)/MSCCR_metadata_from_GEO.xlsx'
excel = pd.read_excel(xlsx,sheet_name=1)
excel.to_csv(out/'audit/Excel_metadata_original_2490.tsv',sep='\t',index=False)
notes = pd.read_excel(xlsx,sheet_name=0)
notes.to_csv(out/'audit/Excel_field_notes.tsv',sep='\t',index=False)
comp = excel.iloc[:,:2].copy()
comp.columns=['geo_accession','Excel_patient_label']
comp = comp.merge(geo[['geo_accession','title','SampleID','GEO_PatientID']],on='geo_accession',how='outer',validate='one_to_one',indicator=True)
comp['exact_title_match'] = comp.Excel_patient_label.eq(comp.SampleID)
comp.to_csv(out/'audit/Excel_GEO_ID_audit.tsv',sep='\t',index=False)
print('Excel title mismatch', (~comp.exact_title_match).sum(), flush=True)

# Step4. 检查冻结矩阵的方向、状态全集和 annotation 样本关系
frozen = root/'ctspec_region_ibd'
members = pd.read_csv(frozen/'Ecotypes/ecotypes.txt',sep='\t')
assert len(members)==53 and members.ID.is_unique and members.CellType.nunique()==12 and members.Ecotype.nunique()==7
annotation = pd.read_csv(workspace/'Profile/audit/2026-10-06_ctspec_region_ibd/tables/01_discovery_sample_metadata.tsv',sep='\t')
assert len(annotation)==2025 and annotation.ID.is_unique and annotation.PatientID.nunique()==919
states = []
for cell_type in members.CellType.unique():
    matrix = pd.read_csv(frozen/cell_type/'state_abundances.txt',sep='\t',index_col=0)
    matrix.index = [cell_type+'_'+x for x in matrix.index]
    matrix = matrix.loc[members.loc[members.CellType.eq(cell_type),'ID']]
    assert set(matrix.columns)==set(annotation.ID)
    states.append(matrix.T)
states = pd.concat(states,axis=1).loc[annotation.ID]
ecotypes = pd.read_csv(frozen/'Ecotypes/ecotype_abundance.txt',sep='\t',index_col=0).T.loc[annotation.ID]
assert states.shape==(2025,53) and ecotypes.shape==(2025,7)
assert np.isfinite(states).all().all() and np.isfinite(ecotypes).all().all()
states.to_csv(out/'source_data/state_abundance_2025x53.tsv',sep='\t')
ecotypes.to_csv(out/'source_data/ecotype_abundance_2025x7.tsv',sep='\t')
members.to_csv(out/'source_data/frozen_state_membership.tsv',sep='\t',index=False)
joined = annotation.merge(geo,left_on='ID',right_on='SampleID',how='left',validate='one_to_one',indicator=True)
assert joined._merge.eq('both').all()
joined.to_csv(out/'source_data/Discovery_GEO_full_mapping.tsv',sep='\t',index=False)
(out/'audit/phase1_local_checkpoint.json').write_text(json.dumps(dict(LocalFiles=len(rows),DiscoveryBiopsies=len(annotation),Patients=annotation.PatientID.nunique(),States=53,Ecotypes=7,GEOFields=geo.columns.tolist(),ExcelTitleMismatches=int((~comp.exact_title_match).sum())),indent=2),encoding='utf-8')
print('Frozen dimensions',states.shape,ecotypes.shape,'patients',annotation.PatientID.nunique(),flush=True)
