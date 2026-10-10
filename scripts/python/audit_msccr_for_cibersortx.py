# 核对服务器 MSCCR 已注释 counts 的格式、数值、样本注释及三套 signature 基因覆盖。
# 在 ibd_refaudit 环境运行；逐基因读取以避免加载完整矩阵，来源文件保持只读。
# 输出审计表，不把现有 R 格式文件当作已经可上传的网页输入。
import argparse
import csv
import gzip
import hashlib
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

# Step1. 读取已知来源及 signature，输出必须是新目录
parser = argparse.ArgumentParser()
parser.add_argument('--outdir', required=True)
parser.add_argument('--signature-dir', required=True)
args = parser.parse_args()
out = Path(args.outdir)
out.mkdir(parents=True, exist_ok=False)
project = Path('/home/mazekai/IBD_EcoTyper')
harmonized = project / '03_reference/combined_reference/output/09_bayesprism_input_harmonization'
source = harmonized / 'GSE193677_MSCCR_Biopsy_counts_GRCh37_Ensembl75_gene_symbol_unique.txt.gz'
annotation_audit = json.loads((harmonized / 'MSCCR_Ensembl75_annotation_audit.json').read_text())
mapping = pd.read_csv(harmonized / 'MSCCR_Ensembl75_gene_mapping.tsv', sep='\t')
signatures = {}
for file in sorted(Path(args.signature_dir).glob('*signature_gene_symbols_restored.txt')):
    signatures[file.name] = pd.read_csv(file, sep='\t', index_col=0)
print('审计来源:', source, '已恢复名称的 signature 数量:', len(signatures), flush=True)

# Step2. 保留样本身份，核对 discovery / CD+UC 子集注释
annotation_paths = [project / '05_EcoTyper/input/MSCCR_CORE12/annotation/MSCCR_CORE12_annotation.tsv', project / '05_EcoTyper/input/MSCCR_CORE12/region_ibd/annotation/MSCCR_CORE12_IBD_annotation.tsv']
annotation_records = []
annotation_tables = []
for path in annotation_paths:
    table = pd.read_csv(path, sep='\t', dtype=str)
    annotation_tables.append(table)
    annotation_records.append({'path': str(path), 'rows': len(table), 'columns': table.columns.tolist()})
    print('样本注释:', path.name, table.shape, '列名:', table.columns.tolist(), flush=True)

# Step3. 扫描整个矩阵：R write.table 的空格分隔/引号，与网站 TSV 格式不同
genes = []
positive_genes = []
noninteger_values = 0
zero_rows = 0
matrix_max = 0.0
with gzip.open(source, 'rt') as handle:
    header = handle.readline().rstrip('\n')
    samples = next(csv.reader([header], delimiter=' ', skipinitialspace=True))
    assert len(samples) == annotation_audit['msccr_samples'] and len(set(samples)) == len(samples)
    library = np.zeros(len(samples))
    detected = np.zeros(len(samples), dtype=int)
    for line_number, line in enumerate(handle, 2):
        match = re.match(r'^"([^"]+)"\s+(.*)$', line.rstrip('\n'))
        assert match is not None, f'Unrecognized gene row {line_number}'
        gene = match.group(1)
        values = np.fromstring(match.group(2), sep=' ')
        assert len(values) == len(samples) and np.isfinite(values).all() and (values >= 0).all(), line_number
        genes.append(gene)
        library += values
        detected += values > 0
        noninteger_values += int((values != np.rint(values)).sum())
        matrix_max = max(matrix_max, float(values.max()))
        if (values > 0).any():
            positive_genes.append(gene)
        else:
            zero_rows += 1
        if len(genes) % 10000 == 0:
            print('已扫描基因:', len(genes), flush=True)
assert len(set(genes)) == len(genes) and len(genes) == annotation_audit['converted_matrix_rows']
assert (library > 0).all()
pd.DataFrame({'sample': samples, 'total_exported_counts': library, 'detected_genes': detected}).to_csv(out / 'MSCCR_sample_library_audit.tsv', sep='\t', index=False)
pd.DataFrame({'gene_symbol': genes}).to_csv(out / 'MSCCR_available_gene_symbols.tsv', sep='\t', index=False)
pd.DataFrame({'gene_symbol': positive_genes}).to_csv(out / 'MSCCR_nonzero_gene_symbols.tsv', sep='\t', index=False)

# Step4. 逐套、逐细胞类型核对覆盖；数字覆盖和 signature 信号覆盖是不同指标
available = set(genes)
nonzero = set(positive_genes)
coverage = []
signal_coverage = []
missing_rows = []
for name, signature in signatures.items():
    present = signature.index.isin(available)
    detectable = signature.index.isin(nonzero)
    coverage.append({'signature': name, 'signature_genes': len(signature), 'available_genes': int(present.sum()), 'available_fraction': float(present.mean()), 'nonzero_somewhere_genes': int(detectable.sum()), 'nonzero_somewhere_fraction': float(detectable.mean())})
    for cell_type in signature.columns:
        values = signature[cell_type].to_numpy()
        signal_coverage.append({'signature': name, 'cell_type': cell_type, 'available_signature_signal_fraction': float(values[present].sum() / values.sum()), 'nonzero_somewhere_signature_signal_fraction': float(values[detectable].sum() / values.sum())})
    missing = signature.index[~present]
    for gene in missing:
        match = mapping.loc[mapping.ensembl75_gene_name.eq(gene)]
        missing_rows.append({'signature': name, 'gene_symbol': gene, 'ensembl75_mapping_status': '|'.join(match.mapping_status.unique()) if len(match) else 'not_in_this_Ensembl75_gene_map', 'n_msccr_ensembl_rows': len(match)})
pd.DataFrame(coverage).to_csv(out / 'signature_MSCCR_gene_coverage.tsv', sep='\t', index=False)
pd.DataFrame(signal_coverage).to_csv(out / 'signature_MSCCR_signal_coverage.tsv', sep='\t', index=False)
pd.DataFrame(missing_rows).to_csv(out / 'signature_MSCCR_missing_genes.tsv', sep='\t', index=False)

# 样本表的首列在当前文件中是样本 ID；先验证匹配，不推测其它列就是身份字段。
for record, table in zip(annotation_records, annotation_tables):
    id_column = table.columns[0]
    ids = table[id_column].tolist()
    record['sample_id_column'] = id_column
    record['unique_sample_ids'] = len(set(ids))
    record['sample_ids_in_expression'] = sum(x in set(samples) for x in ids)
    record['missing_sample_id_examples'] = [x for x in ids if x not in set(samples)][:10]
    assert len(set(ids)) == len(ids)
    pd.DataFrame({'sample': ids, 'in_expression': [x in set(samples) for x in ids]}).to_csv(out / (Path(record['path']).stem + '_sample_match.tsv'), sep='\t', index=False)

# Step5. 输出结论所需事实；不把 raw counts、批次调整表达或临床结果混为一谈
summary = {'source': str(source), 'source_bytes': source.stat().st_size, 'source_genes': len(genes), 'source_samples': len(samples), 'data_points': len(genes) * len(samples),
    'format': 'space-separated R-style quoted header and gene labels', 'web_upload_ready': False, 'max_expression': matrix_max,
    'all_values_finite_nonnegative': True, 'fractional_count_entries': noninteger_values, 'all_zero_rows': zero_rows, 'nonzero_somewhere_genes': len(positive_genes),
    'scale': 'unlogged count values; annotation conversion performed no normalization', 'annotation_audit': annotation_audit,
    'annotations': annotation_records, 'signature_coverage': coverage, 'duplicate_symbol_strategy': 'existing unique-symbol rows only, no aggregation; excludes all ambiguous and unmapped rows',
    'original_annotation_sources_mutated': False, 'real_MSCCR_deconvolution_run': False, 'clinical_analysis_performed': False}
(out / 'summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps({k: v for k, v in summary.items() if k != 'annotation_audit'}, ensure_ascii=False, indent=2), flush=True)
