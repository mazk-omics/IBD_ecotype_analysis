# 准备冻结的 Region3 CD+UC 模型在 nonIBD/健康肠道中的 recovery；不改变发现模型或原始输入。
# 在服务器 ibd_refaudit 环境运行，输出至独立目录；MSCCR 按部位分开，单细胞沿用此前抽样。
# 本轮只评价恢复、基因覆盖和供者支持，不做状态命名或临床关联。
from pathlib import Path
import argparse
import csv
import gzip
import hashlib
import json
import numpy as np
import pandas as pd

# Step1. 固定发现模型，保存源码与原输入指纹
parser = argparse.ArgumentParser()
parser.add_argument('--outdir', required=True)
args = parser.parse_args()
out = Path(args.outdir)
out.mkdir(parents=True, exist_ok=False)
project = Path('/home/mazekai/IBD_EcoTyper')
software = project / 'software/ecotyper'
model = software / 'EcoTyper/MSCCR_CORE12_CTSPEC_REGION_IBD/Cell_type_specific_genes'
rank = pd.read_csv(model / 'Analysis/rank_selection/rank_data.txt', sep='\t')
model_files = [model / 'Analysis/rank_selection/rank_data.txt', model.parent / 'config_used.yml']
# 只校验恢复实际使用的已选 rank 和生态型文件，不读取未选 rank 的大型 NMF 重启归档。
for cell_type, n_states in rank.iloc[:, :2].itertuples(index=False, name=None):
    state_dir = model / 'Cell_States/discovery' / cell_type / str(n_states)
    model_files.extend(state_dir / f for f in ['W.txt', 'mapping_to_initial_states.txt', 'initial_gene_info.txt', 'gene_info.txt', 'heatmap_data.txt', 'heatmap_top_ann.txt'])
model_files.extend(p for p in (model / 'Ecotypes/discovery').glob('*.txt'))
source_files = [software / p for p in ['EcoTyper_recovery_bulk.R', 'EcoTyper_recovery_scRNA.R', 'pipeline/state_recovery_bulk.R', 'pipeline/state_recovery_scRNA.R', 'pipeline/ecotypes_recovery.R', 'pipeline/ecotypes_recovery_scRNA.R', 'pipeline/lib/ecotyper.R']]
fingerprints = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in model_files + source_files if p.is_file()}
(out / 'frozen_model_source_before.json').write_text(json.dumps(fingerprints, indent=2), encoding='utf-8')
gene_info = []
for cell_type, n_states in rank.iloc[:, :2].itertuples(index=False, name=None):
    genes = pd.read_csv(model / 'Cell_States/discovery' / cell_type / str(n_states) / 'initial_gene_info.txt', sep='\t')
    gene_info.append(pd.DataFrame({'CellType': cell_type, 'Gene': genes.Gene}))
model_genes = pd.concat(gene_info, ignore_index=True)
gene_union = np.sort(model_genes.Gene.unique())
print('Frozen model:', rank.shape, 'gene union:', len(gene_union), flush=True)
manifest = []
coverage = []
input_hashes = {}

# Step2. MSCCR 只取明确的 Control、非炎症活检；检查供者不进入 CD/UC discovery
ann_file = project / '05_EcoTyper/input/MSCCR_CORE12/annotation/MSCCR_CORE12_annotation.tsv'
annotation = pd.read_csv(ann_file, sep='\t', dtype=str)
discovery = pd.read_csv(project / '05_EcoTyper/input/MSCCR_CORE12/region_ibd/annotation/MSCCR_CORE12_IBD_annotation.tsv', sep='\t', dtype=str)
control = annotation.loc[annotation.Disease.eq('Control') & annotation.Inflammation.eq('Non-inflamed')].copy()
assert annotation.ID.is_unique and control.ID.is_unique
assert not set(control.ID) & set(discovery.ID)
assert not set(control.PatientID) & set(discovery.PatientID)
control.to_csv(out / 'MSCCR_nonIBD_selected_annotation.tsv', sep='\t', index=False)
print('MSCCR controls:', len(control), 'donors:', control.PatientID.nunique(), flush=True)

# CPM 分母包含全部原始 Ensembl 行；只导出模型需要的基因，重复 symbol 的计数相加。
raw_file = project / '01_raw_data/bulk/MSCCR_GSE193677/GSE193677_MSCCR_Biopsy_counts.txt.gz'
identity_file = project / '03_reference/web_reference_three_pilots/msccr_gene_identity_20261006/MSCCR_Ensembl75_HGNC_identity.tsv'
identity = pd.read_csv(identity_file, sep='\t', keep_default_na=False)
identity = identity.loc[identity.mapping_status.eq('unique_gene_symbol') & identity.ensembl75_gene_name.isin(gene_union)].copy()
assert identity.msccr_ensembl_id.is_unique
symbol_by_id = identity.set_index('msccr_ensembl_id').ensembl75_gene_name.to_dict()
gene_position = dict(zip(gene_union, range(len(gene_union))))
counts = np.zeros((len(gene_union), len(control)))
mapped_rows = np.zeros(len(gene_union), dtype=int)
with gzip.open(raw_file, 'rt') as handle:
    samples = next(csv.reader([handle.readline().strip()], delimiter=' ', skipinitialspace=True))
    positions = pd.Index(samples).get_indexer(control.ID)
    assert min(positions) >= 0 and len(set(samples)) == len(samples)
    library = np.zeros(len(control))
    seen = set()
    for row_number, line in enumerate(handle, 1):
        raw_gene, expression = line.split(None, 1)
        raw_gene = raw_gene.strip('"')
        values = np.fromstring(expression, sep=' ')
        assert len(values) == len(samples) and np.isfinite(values).all() and (values >= 0).all() and raw_gene not in seen
        seen.add(raw_gene)
        library += values[positions]
        if raw_gene in symbol_by_id:
            row = gene_position[symbol_by_id[raw_gene]]
            counts[row] += values[positions]
            mapped_rows[row] += 1
    print('MSCCR raw rows scanned:', row_number, flush=True)
assert (library > 0).all()
present = mapped_rows > 0
cpm = pd.DataFrame(counts[present] / library * 1e6, index=gene_union[present], columns=control.ID)
assert np.isfinite(cpm.to_numpy()).all() and cpm.to_numpy().max() > 50
pd.DataFrame({'Gene': gene_union, 'MappedRawRows': mapped_rows}).to_csv(out / 'MSCCR_model_gene_mapping.tsv', sep='\t', index=False)
control.assign(FullRawLibrary=library).to_csv(out / 'MSCCR_library_audit.tsv', sep='\t', index=False)
for region in ['Ileum', 'Colon', 'Rectum']:
    name = 'MSCCR_nonIBD_' + region + '_20261010'
    selected = control.loc[control.Region3.eq(region)].copy()
    selected.to_csv(out / (name + '_annotation.tsv'), sep='\t', index=False)
    cpm.loc[:, selected.ID].to_csv(out / (name + '.txt'), sep='\t', index_label='GeneSymbol', float_format='%.9g')
    manifest.append({'Dataset': name, 'Assay': 'bulk', 'SamplesOrCells': len(selected), 'Donors': selected.PatientID.nunique(), 'Region': region, 'Scale': 'linear CPM, all raw rows denominator', 'Independence': 'same cohort; sample and donor excluded from discovery'})
    for cell_type, g in model_genes.groupby('CellType', sort=True):
        coverage.append({'Dataset': name, 'CellType': cell_type, 'ModelGenes': len(g), 'PresentGenes': g.Gene.isin(cpm.index).sum()})

# Step3. 外部 bulk 对照：GSE57945 已有 RPKM，不重新归一化，不合并同 RISK 的另一批结果
risk = project / '06_validation/01_bulk_progression/GSE57945/prepared'
risk_ann_file = risk / '04_sample_metadata/GSE57945_canonical_sample_metadata_322.tsv'
risk_expr_file = risk / '06_ecotyper_recovery_input/GSE57945_EcoTyper_recovery_expression.txt'
risk_ann = pd.read_csv(risk_ann_file, sep='\t', dtype=str)
risk_selected = risk_ann.loc[risk_ann.Disease.eq('Non-IBD')].copy().rename(columns={'SampleID': 'ID'})
assert risk_selected.ID.is_unique and risk_selected.Tissue.eq('Ileum').all()
risk_expr = pd.read_csv(risk_expr_file, sep='\t', index_col=0)
assert risk_expr.index.is_unique and risk_expr.columns.is_unique and set(risk_selected.ID) <= set(risk_expr.columns)
risk_subset = risk_expr.loc[:, risk_selected.ID]
assert np.isfinite(risk_subset.to_numpy()).all() and (risk_subset.to_numpy() >= 0).all() and risk_subset.to_numpy().max() > 50
name = 'GSE57945_nonIBD_Ileum_20261010'
risk_subset.to_csv(out / (name + '.txt'), sep='\t', index_label='GeneSymbol', float_format='%.9g')
risk_selected.to_csv(out / (name + '_annotation.tsv'), sep='\t', index=False)
manifest.append({'Dataset': name, 'Assay': 'bulk', 'SamplesOrCells': len(risk_selected), 'Donors': len(risk_selected), 'Region': 'Ileum', 'Scale': 'deposited RPKM unchanged', 'Independence': 'external paediatric RISK; distinct samples; donor identity follows existing sample metadata'})
for cell_type, g in model_genes.groupby('CellType', sort=True):
    coverage.append({'Dataset': name, 'CellType': cell_type, 'ModelGenes': len(g), 'PresentGenes': g.Gene.isin(risk_subset.index).sum()})

# Step4. 单细胞健康子集：仅筛选既有细胞，保留 CPM 数值和既有供者，不再挑选恢复效果好的细胞
cohorts = [('IBDverse_Healthy_Ileum_20261010', 'IBDverse_scRNA_recovery_20261009/input_v2', 'IBDverse_Pilot30_20261009', 'Healthy', 'Ileum'), ('GSE214695_HC_Colon_20261010', 'GSE214695_scRNA_recovery_20261009/input', 'GSE214695_Pilot18_20261009', 'HC', 'Colon')]
for name, folder, prefix, disease, region in cohorts:
    a_file = project / '06_validation' / folder / (prefix + '_annotation.tsv')
    e_file = project / '06_validation' / folder / (prefix + '_CPM.txt')
    a = pd.read_csv(a_file, sep='\t', dtype=str)
    selected = a.loc[a.Disease.eq(disease)].copy()
    expr = pd.read_csv(e_file, sep='\t', index_col=0)
    assert selected.ID.is_unique and expr.index.is_unique and expr.columns.is_unique and set(selected.ID) <= set(expr.columns)
    subset = expr.loc[:, selected.ID]
    assert np.isfinite(subset.to_numpy()).all() and (subset.to_numpy() >= 0).all() and subset.to_numpy().max() > 50
    selected.to_csv(out / (name + '_annotation.tsv'), sep='\t', index=False)
    subset.to_csv(out / (name + '.txt'), sep='\t', index_label='GeneSymbol', float_format='%.9g')
    selected.groupby(['Donor', 'CellType']).size().rename('Cells').reset_index().to_csv(out / (name + '_donor_type_counts.tsv'), sep='\t', index=False)
    manifest.append({'Dataset': name, 'Assay': 'scRNA', 'SamplesOrCells': len(selected), 'Donors': selected.Donor.nunique(), 'Region': region, 'Scale': 'existing full-library CPM unchanged', 'Independence': 'existing recovery donors; subgroup re-analysis, not a new independent cohort'})
    for cell_type, g in model_genes.groupby('CellType', sort=True):
        coverage.append({'Dataset': name, 'CellType': cell_type, 'ModelGenes': len(g), 'PresentGenes': g.Gene.isin(subset.index).sum(), 'Cells': selected.CellType.eq(cell_type).sum()})
    for p in [a_file, e_file]:
        input_hashes[str(p)] = hashlib.sha256(p.read_bytes()).hexdigest()

# Step5. 保存可复核清单；官方 bulk 默认无状态置换检验，单细胞另用官方 -z TRUE
for p in [ann_file, raw_file, identity_file, risk_ann_file, risk_expr_file]:
    input_hashes[str(p)] = hashlib.sha256(p.read_bytes()).hexdigest()
(out / 'source_inputs_before.json').write_text(json.dumps(input_hashes, indent=2), encoding='utf-8')
pd.DataFrame(manifest).to_csv(out / 'run_manifest.tsv', sep='\t', index=False)
pd.DataFrame(coverage).to_csv(out / 'model_gene_coverage.tsv', sep='\t', index=False)
print(pd.DataFrame(manifest).to_string(index=False), flush=True)
print(pd.DataFrame(coverage).groupby('Dataset').PresentGenes.agg(['min', 'max']).to_string(), flush=True)
