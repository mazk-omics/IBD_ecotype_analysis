# 独立抽查 MSCCR 网页输入：从原始 counts 重新算完整库总量，核对改名 marker 和列顺序。
# 在服务器 ibd_refaudit 环境运行；读取已完成的导出目录，不修改原始输入或网页矩阵。
# 只保存抽查报告；它补充全矩阵回读检查，不替代所有基因身份的人工审阅。
import csv
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd

# Step1. 选取可解释的基因和分散的样本，包含改名及全零行
project = Path('/home/mazekai/IBD_EcoTyper')
out = project / '03_reference/web_reference_three_pilots/msccr_ibd_joint_web_20261006'
identity = pd.read_csv(out / 'exported_gene_identity.tsv', sep='\t')
annotation = pd.read_csv(out / 'MSCCR_CD_UC_2025_annotation.tsv', sep='\t')
markers = ['JCHAIN', 'CXCL8', 'SELENOP', 'FYB1', 'ACKR1', 'TCIM', 'EPCAM', 'CD3D', 'CD79A', 'SFT2D3']
chosen = identity.set_index('gene_symbol').loc[markers]
sample_positions = np.linspace(0, len(annotation) - 1, 12, dtype=int)
sample_ids = annotation.ID.iloc[sample_positions].tolist()
print('抽查基因:', markers, '抽查样本:', len(sample_ids), flush=True)

# Step2. 重新扫描全部原始基因，只计算抽查样本的完整库总量
source = project / '01_raw_data/bulk/MSCCR_GSE193677/GSE193677_MSCCR_Biopsy_counts.txt.gz'
raw_library = np.zeros(len(sample_ids))
raw_selected = {}
with gzip.open(source, 'rt') as handle:
    raw_header = next(csv.reader([handle.readline().strip()], delimiter=' ', skipinitialspace=True))
    raw_positions = pd.Index(raw_header).get_indexer(sample_ids)
    assert (raw_positions >= 0).all()
    for raw_rows, line in enumerate(handle, 1):
        fields = line.split(None, 1)
        gene_id = fields[0].strip('"')
        values = np.fromstring(fields[1], sep=' ')[raw_positions]
        raw_library += values
        if gene_id in set(chosen.mapped_msccr_ensembl_id):
            raw_selected[gene_id] = values
assert raw_rows == 56632 and len(raw_selected) == len(markers)
expected_counts = np.stack([raw_selected[gene_id] for gene_id in chosen.mapped_msccr_ensembl_id])
expected_cpm = expected_counts / raw_library * 1e6

# Step3. 直接读取最终 TXT，抽查列顺序、改名映射和 CPM 分母
web_values = {}
with (out / 'MSCCR_CD_UC_2025_GJ20_CPM.txt').open() as handle:
    web_header = handle.readline().rstrip('\n').split('\t')
    assert web_header[1:] == annotation.ID.tolist()
    for line in handle:
        fields = line.rstrip('\n').split('\t')
        if fields[0] in markers:
            web_values[fields[0]] = np.array([float(fields[position + 1]) for position in sample_positions])
observed_cpm = np.stack([web_values[gene] for gene in markers])
assert np.allclose(observed_cpm, expected_cpm, rtol=6e-8, atol=1e-7)
stored_library = pd.read_csv(out / 'sample_library_and_order_audit.tsv', sep='\t').full_raw_library_counts.iloc[sample_positions].to_numpy()
assert np.allclose(stored_library, raw_library, rtol=1e-12)
report = {'status': 'passed', 'source_raw_rows_recomputed': raw_rows, 'samples': sample_ids, 'genes': markers,
    'independent_raw_to_web_values_checked': int(observed_cpm.size), 'column_order_matches_existing_annotation': True,
    'full_raw_library_recomputed_matches_saved_audit': True, 'max_absolute_CPM_difference': float(np.abs(observed_cpm - expected_cpm).max())}
(out / 'independent_raw_spotcheck.json').write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps(report, ensure_ascii=False, indent=2), flush=True)
