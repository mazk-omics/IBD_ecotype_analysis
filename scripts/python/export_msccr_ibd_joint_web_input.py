# 导出 MSCCR 的 CD+UC discovery 网页输入，按已核对身份对应联合参考背景。
# 在服务器 ibd_refaudit 环境运行；读取原始 Ensembl counts 与 2,025 样本注释，不修改来源。
# 所有原始基因先参与库总量计算，再取一对一对应行转非 log CPM；不合并、填补或舍入 counts。
import argparse
import csv
import gzip
import hashlib
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

# Step1. 固定用户选择的样本集合，按注释 ID 顺序导出
parser = argparse.ArgumentParser()
parser.add_argument('--outdir', required=True)
parser.add_argument('--identity-dir', required=True)
args = parser.parse_args()
out = Path(args.outdir)
out.mkdir(parents=True, exist_ok=False)
started = time.monotonic()
project = Path('/home/mazekai/IBD_EcoTyper')
source = project / '01_raw_data/bulk/MSCCR_GSE193677/GSE193677_MSCCR_Biopsy_counts.txt.gz'
annotation_path = project / '05_EcoTyper/input/MSCCR_CORE12/region_ibd/annotation/MSCCR_CORE12_IBD_annotation.tsv'
identity = pd.read_csv(Path(args.identity_dir) / 'joint_reference_to_MSCCR_gene_bridge.tsv', sep='\t', keep_default_na=False)
annotation = pd.read_csv(annotation_path, sep='\t', dtype=str, keep_default_na=False)
assert len(annotation) == 2025 and annotation.ID.is_unique and set(annotation.Disease) == {'CD', 'UC'}
accepted = identity.loc[identity.accepted.eq(True)].copy()
assert accepted.gene_symbol.is_unique and accepted.mapped_msccr_ensembl_id.is_unique
assert len(accepted) * len(annotation) <= 50_000_000
annotation.to_csv(out / 'MSCCR_CD_UC_2025_annotation.tsv', sep='\t', index=False)
annotation.groupby(['Disease', 'Inflammation', 'Region3']).size().rename('n_samples').reset_index().to_csv(out / 'sample_group_counts.tsv', sep='\t', index=False)
accepted.to_csv(out / 'exported_gene_identity.tsv', sep='\t', index=False)
print('用户选定 CD+UC:', len(annotation), '样本;', len(accepted), '已核对基因', flush=True)
source_sha256_before = hashlib.sha256(source.read_bytes()).hexdigest()

# Step2. 逐行读取原始 counts，只为已核对的基因保留表达；库总量包括全部原始行
gene_row = pd.Series(np.arange(len(accepted)), index=accepted.mapped_msccr_ensembl_id).to_dict()
counts = np.zeros((len(accepted), len(annotation)), dtype=np.float64)
seen = set()
raw_seen = set()
fractional_values = 0
with gzip.open(source, 'rt') as handle:
    header = handle.readline().strip()
    samples = next(csv.reader([header], delimiter=' ', skipinitialspace=True))
    assert len(samples) == 2490 and len(set(samples)) == len(samples)
    positions = pd.Index(samples).get_indexer(annotation.ID)
    assert (positions >= 0).all()
    library = np.zeros(len(samples))
    for number, line in enumerate(handle, 1):
        fields = line.split(None, 1)
        gene_id = fields[0].strip('"')
        values = np.fromstring(fields[1], sep=' ')
        assert len(values) == len(samples) and np.isfinite(values).all() and (values >= 0).all(), number
        assert gene_id not in raw_seen
        raw_seen.add(gene_id)
        library += values
        fractional_values += int((values != np.rint(values)).sum())
        if gene_id in gene_row:
            counts[gene_row[gene_id]] = values[positions]
            seen.add(gene_id)
        if number % 10000 == 0:
            print('已扫描原始基因:', number, flush=True)
assert set(gene_row) == seen and len(raw_seen) == 56632 and (library > 0).all()
selected_library = library[positions]
retained_counts = counts.sum(axis=0)
sample_audit = annotation[['ID', 'PatientID', 'Disease', 'Region3']].copy()
sample_audit['raw_source_column_index'] = positions
sample_audit['full_raw_library_counts'] = selected_library
sample_audit['exported_background_counts'] = retained_counts
sample_audit['exported_background_library_fraction'] = retained_counts / selected_library
sample_audit.to_csv(out / 'sample_library_and_order_audit.tsv', sep='\t', index=False)
print('原始矩阵行数:', len(raw_seen), '选中计数:', counts.shape, '原始 counts 非整数条目:', fractional_values, flush=True)

# Step3. 保存原始选中计数检查点，再用完整原始库计算 CPM
# 留下检查点用于重复审计；它仍以 Ensembl 行对应表解释，不作为网页输入。
np.savez_compressed(out / 'selected_raw_counts_checkpoint.npz', counts=counts, library=selected_library)
cpm = counts / selected_library[None, :] * 1e6
assert np.isfinite(cpm).all() and (cpm >= 0).all() and cpm.max() >= 50
gene_audit = accepted.copy()
gene_audit['all_zero_in_IBD'] = counts.sum(axis=1) == 0
gene_audit['total_IBD_raw_counts'] = counts.sum(axis=1)
gene_audit.to_csv(out / 'exported_gene_identity.tsv', sep='\t', index=False)
web = out / 'MSCCR_CD_UC_2025_GJ20_CPM.txt'
pd.DataFrame(cpm, index=accepted.gene_symbol, columns=annotation.ID).to_csv(web, sep='\t', index_label='GeneSymbol', float_format='%.8g')
with web.open('rb') as reader, gzip.open(str(web) + '.gz', 'wb', compresslevel=6) as writer:
    while True:
        chunk = reader.read(4 * 1024 * 1024)
        if not chunk:
            break
        writer.write(chunk)
print('已写网页 TXT:', web.stat().st_size, 'bytes;', cpm.size, '数据点', flush=True)

# Step4. 直接回读整个网页 TXT 核对数字和顺序，未把输出存在当作验证通过
roundtrip = pd.read_csv(web, sep='\t', index_col=0)
assert np.array_equal(roundtrip.index, accepted.gene_symbol) and np.array_equal(roundtrip.columns, annotation.ID)
assert np.isfinite(roundtrip.to_numpy()).all() and (roundtrip.to_numpy() >= 0).all()
assert np.allclose(roundtrip.to_numpy(), cpm, rtol=6e-8, atol=1e-7)
max_roundtrip_error = float(np.abs(roundtrip.to_numpy() - cpm).max())
del roundtrip
source_sha256_after = hashlib.sha256(source.read_bytes()).hexdigest()
assert source_sha256_before == source_sha256_after
files = {}
for file in [web, Path(str(web) + '.gz')]:
    files[file.name] = {'bytes': file.stat().st_size, 'sha256': hashlib.sha256(file.read_bytes()).hexdigest()}
signature_rows = gene_audit.loc[gene_audit.signature_gene.eq(True)]
summary = {'status': 'completed_and_roundtrip_verified', 'source': str(source), 'source_sha256_before': source_sha256_before, 'source_sha256_after': source_sha256_after,
    'source_unchanged': True, 'raw_genes': len(raw_seen), 'raw_samples': len(samples), 'user_selected_sample_set': 'CD+UC discovery, existing region_ibd annotation',
    'exported_samples': len(annotation), 'exported_genes': len(accepted), 'data_points': int(cpm.size), 'web_data_point_limit_passed': True,
    'scale': 'nonlog CPM; full original 56632-gene library denominator before background selection', 'fractional_counts_preserved_before_CPM': True,
    'gene_identity_policy': 'audited one-to-one Ensembl/HGNC correspondence; source signature symbols preserved; unresolved and conflicting rows excluded',
    'reference_signature_genes': 5816, 'mapped_signature_genes': len(signature_rows), 'signature_genes_nonzero_in_IBD': int((~signature_rows.all_zero_in_IBD).sum()),
    'sample_order': 'exact existing IBD annotation ID order', 'source_samples_not_in_selected_annotation': len(samples) - len(annotation),
    'background_CPM_column_sums_are_not_one_million': True, 'no_renormalization_after_background_selection': True, 'max_CPM': float(cpm.max()),
    'max_absolute_roundtrip_CPM_error': max_roundtrip_error, 'files': files, 'seconds': round(time.monotonic() - started, 2),
    'CIBERSORTx_run': False, 'cell_type_expression_recovery_run': False, 'EcoTyper_run': False, 'clinical_analysis_performed': False}
(out / 'summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)
