# 将 MSCCR 原始 Ensembl 行与联合参考的来源基因逐一对应，解释旧名称及缺失基因。
# 在服务器 ibd_refaudit 环境运行；读取已有 Ensembl75 映射和固定 HGNC 快照，不读取表达、不合并位点。
# 输出逐基因对应表和待人工核对列表；只有身份清楚的一对一关系允许进入后续网页导出。
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

# Step1. 固定来源与输出位置
parser = argparse.ArgumentParser()
parser.add_argument('--outdir', required=True)
parser.add_argument('--signature-genes', required=True)
args = parser.parse_args()
out = Path(args.outdir)
out.mkdir(parents=True, exist_ok=False)
project = Path('/home/mazekai/IBD_EcoTyper')
mapping_path = project / '03_reference/combined_reference/output/09_bayesprism_input_harmonization/MSCCR_Ensembl75_gene_mapping.tsv'
reference_path = project / '03_reference/web_reference_three_pilots/output/run_20261006_183211/gse_joint/gene_export_manifest.tsv'
hgnc_path = project / '01_raw_data/annotation/HGNC_20261006/hgnc_complete_set.txt'
bulk = pd.read_csv(mapping_path, sep='\t', dtype=str, keep_default_na=False)
reference_all = pd.read_csv(reference_path, sep='\t', keep_default_na=False)
reference = reference_all.loc[reference_all.web_included.eq(True), ['source_gene_id', 'gene_symbol']].copy()
hgnc = pd.read_csv(hgnc_path, sep='\t', dtype=str, keep_default_na=False)
hgnc = hgnc.loc[hgnc.status.eq('Approved')].copy()
signature_genes = set(Path(args.signature_genes).read_text().splitlines())
assert bulk.msccr_ensembl_id.is_unique and reference.source_gene_id.is_unique and reference.gene_symbol.is_unique
assert hgnc.hgnc_id.is_unique and hgnc.symbol.is_unique and signature_genes <= set(reference.gene_symbol)
print('来源维度:', len(bulk), 'MSCCR Ensembl 行;', len(reference), '联合参考背景基因;', len(signature_genes), 'signature 基因', flush=True)

# Step2. 只把唯一 HGNC 旧名称/别名用于身份核对，批准名称优先，不自动大写或去标点
approved = hgnc.set_index('symbol').hgnc_id.to_dict()
alias_tables = []
for field in ['prev_symbol', 'alias_symbol']:
    table = hgnc[['hgnc_id', field]].copy()
    table[field] = table[field].str.split('|')
    table = table.explode(field).rename(columns={field: 'old_symbol'})
    table = table.loc[table.old_symbol.ne('')].drop_duplicates()
    alias_tables.append(table)
aliases = pd.concat(alias_tables, ignore_index=True).drop_duplicates()
alias_number = aliases.groupby('old_symbol').hgnc_id.nunique()
unique_aliases = aliases.loc[aliases.old_symbol.isin(alias_number[alias_number.eq(1)].index)]
alias_lookup = unique_aliases.set_index('old_symbol').hgnc_id.to_dict()
symbol_lookup = {**alias_lookup, **approved}
bulk['old_hgnc_id'] = bulk.ensembl75_gene_name.map(symbol_lookup).fillna('')
reference['target_hgnc_id'] = reference.gene_symbol.map(symbol_lookup).fillna('')
bulk.to_csv(out / 'MSCCR_Ensembl75_HGNC_identity.tsv', sep='\t', index=False)

# Step3. 先核对相同 Ensembl ID 的基因身份；相同 ID 但 HGNC 身份不同不能强行改名
raw_by_id = bulk.set_index('msccr_ensembl_id')
reference['mapped_msccr_ensembl_id'] = ''
reference['old_gene_symbol'] = ''
reference['method'] = 'unresolved'
reference['identity_note'] = ''
reference['signature_gene'] = reference.gene_symbol.isin(signature_genes)
for index, row in reference.iterrows():
    source_id = row.source_gene_id
    if source_id not in raw_by_id.index:
        continue
    old = raw_by_id.loc[source_id]
    same_symbol = old.ensembl75_gene_name == row.gene_symbol
    same_hgnc = row.target_hgnc_id != '' and old.old_hgnc_id == row.target_hgnc_id
    reference.at[index, 'old_gene_symbol'] = old.ensembl75_gene_name
    if same_symbol or same_hgnc:
        reference.at[index, 'mapped_msccr_ensembl_id'] = source_id
        reference.at[index, 'method'] = 'exact_ID_same_symbol' if same_symbol else 'exact_ID_HGNC_verified_name_change'
    else:
        reference.at[index, 'method'] = 'quarantine_ID_symbol_disagreement'
        reference.at[index, 'identity_note'] = 'Same Ensembl ID but symbols are not the same uniquely resolved HGNC gene; do not relabel counts'

# Step4. 对来源 ID 已改变的基因，只接受唯一同名或唯一同一 HGNC 基因的旧行
# 不按表达高低挑位点，不合并多个 Ensembl 行；有 ID 身份冲突的行保留隔离状态。
symbol_groups = bulk.groupby('ensembl75_gene_name', sort=False).groups
hgnc_groups = bulk.loc[bulk.old_hgnc_id.ne('')].groupby('old_hgnc_id', sort=False).groups
target_hgnc_number = reference.loc[reference.target_hgnc_id.ne('')].groupby('target_hgnc_id').source_gene_id.nunique()
for index, row in reference.loc[reference.method.eq('unresolved')].iterrows():
    candidates = list(symbol_groups.get(row.gene_symbol, []))
    method = 'unique_same_symbol_changed_ID'
    if not candidates and row.target_hgnc_id != '' and target_hgnc_number.get(row.target_hgnc_id, 0) == 1:
        candidates = list(hgnc_groups.get(row.target_hgnc_id, []))
        method = 'unique_HGNC_gene_changed_ID'
    if len(candidates) == 1:
        old = bulk.loc[candidates[0]]
        reference.at[index, 'mapped_msccr_ensembl_id'] = old.msccr_ensembl_id
        reference.at[index, 'old_gene_symbol'] = old.ensembl75_gene_name
        reference.at[index, 'method'] = method
    elif len(candidates) > 1:
        reference.at[index, 'method'] = 'quarantine_multiple_MSCCR_loci'
        reference.at[index, 'identity_note'] = '|'.join(bulk.loc[candidates, 'msccr_ensembl_id'])

# 一个 bulk 行若被多个来源基因争用，全部隔离，避免复制同一计数向量。
mapped = reference.mapped_msccr_ensembl_id.ne('')
duplicate_rows = reference.loc[mapped, 'mapped_msccr_ensembl_id'].duplicated(keep=False)
conflicting_indices = reference.loc[mapped].index[duplicate_rows]
reference.loc[conflicting_indices, 'method'] = 'quarantine_one_raw_row_multiple_targets'
reference.loc[conflicting_indices, 'mapped_msccr_ensembl_id'] = ''
reference['accepted'] = reference.mapped_msccr_ensembl_id.ne('')
reference.to_csv(out / 'joint_reference_to_MSCCR_gene_bridge.tsv', sep='\t', index=False)
reference.loc[~reference.accepted].to_csv(out / 'unresolved_or_quarantined_reference_genes.tsv', sep='\t', index=False)
reference.loc[reference.signature_gene & ~reference.accepted].to_csv(out / 'unresolved_or_quarantined_signature_genes.tsv', sep='\t', index=False)
reference.loc[reference.accepted & reference.gene_symbol.ne(reference.old_gene_symbol)].to_csv(out / 'verified_gene_name_changes.tsv', sep='\t', index=False)
print('映射方法:', reference.method.value_counts().to_dict(), flush=True)
print('关键基因:', reference.loc[reference.gene_symbol.isin(['JCHAIN', 'PECAM1', 'CD24', 'CXCL8', 'SELENOP', 'FYB1', 'ACKR1', 'PRSS2', 'TCIM', 'TRBC1']), ['gene_symbol', 'source_gene_id', 'old_gene_symbol', 'method', 'accepted']].to_string(index=False), flush=True)

# Step5. 冻结来源版本，报告覆盖；这里只完成身份验证，不宣称已经完成反卷积
sig = reference.loc[reference.signature_gene]
summary = {'status': 'completed', 'background_genes': len(reference), 'mapped_background_genes': int(reference.accepted.sum()),
    'signature_genes': len(sig), 'mapped_signature_genes': int(sig.accepted.sum()), 'mapped_signature_fraction': float(sig.accepted.mean()),
    'method_counts': reference.method.value_counts().to_dict(), 'hgnc_source_url': 'https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt',
    'hgnc_snapshot_sha256': hashlib.sha256(hgnc_path.read_bytes()).hexdigest(), 'ensembl75_mapping_sha256': hashlib.sha256(mapping_path.read_bytes()).hexdigest(),
    'reference_gene_manifest_sha256': hashlib.sha256(reference_path.read_bytes()).hexdigest(), 'reference_source_gene_symbols_preserved': True,
    'ambiguous_loci_aggregated': False, 'counts_read_or_modified': False, 'clinical_analysis_performed': False}
(out / 'summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)
