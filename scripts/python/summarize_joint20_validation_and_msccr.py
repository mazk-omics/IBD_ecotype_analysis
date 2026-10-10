# 本地整理留出测试包和 MSCCR 输入审计，计算每套 signature 的基因及信号覆盖。
# 使用 Codex 已有 Python/numpy/pandas；读取服务器回传结果和本地已审阅的矩阵。
# 输出网页上传包及覆盖表，不运行反卷积；矩阵数值不作修改。
import hashlib
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

# Step1. 读入回传审计和已恢复基因名称的 signature
project = Path(r'D:\IBD_ecotype\IBD_ecotype_analysis')
root = project / 'CIBERSORTx_signature_audit/main_candidate_validation_20261006'
audit = root / 'server_outputs/msccr_input_audit_20261006'
validation = root / 'server_outputs/validation_joint20_20261006'
comparison = project / 'CIBERSORTx_signature_audit/three_candidate_comparison_20261006'
upload = root / 'web_upload'
upload.mkdir(exist_ok=True)
bulk_summary = json.loads((audit / 'summary.json').read_text(encoding='utf-8'))
validation_summary = json.loads((validation / 'summary.json').read_text(encoding='utf-8'))
available = set(pd.read_csv(audit / 'MSCCR_available_gene_symbols.tsv', sep='\t').gene_symbol)
nonzero = set(pd.read_csv(audit / 'MSCCR_nonzero_gene_symbols.tsv', sep='\t').gene_symbol)
main_mixture = pd.read_csv(validation / 'GJ20_Holdout28_Mixture_CPM.txt', sep='\t', index_col=0)
mixture_genes = set(main_mixture.index)

# Step2. 逐套核对 MSCCR 和测试表达的覆盖，保留缺失基因的类型及表达权重
coverage = []
type_coverage = []
missing = []
for mode in ['gse_joint', 'gse_ileum', 'scp_rectum']:
    path = next((comparison / mode).glob('*signature_gene_symbols_restored.txt'))
    signature = pd.read_csv(path, sep='\t', index_col=0)
    bulk_present = signature.index.isin(available)
    bulk_detected = signature.index.isin(nonzero)
    test_present = signature.index.isin(mixture_genes)
    coverage.append({'candidate': mode, 'signature_genes': len(signature), 'MSCCR_available_genes': int(bulk_present.sum()), 'MSCCR_available_fraction': float(bulk_present.mean()),
        'MSCCR_nonzero_somewhere_genes': int(bulk_detected.sum()), 'test_available_genes': int(test_present.sum()), 'test_available_fraction': float(test_present.mean())})
    for cell_type in signature.columns:
        values = signature[cell_type].to_numpy()
        type_coverage.append({'candidate': mode, 'cell_type': cell_type, 'MSCCR_available_signal_fraction': float(values[bulk_present].sum() / values.sum()),
            'MSCCR_nonzero_somewhere_signal_fraction': float(values[bulk_detected].sum() / values.sum()), 'test_available_signal_fraction': float(values[test_present].sum() / values.sum())})
    for gene in signature.index[~bulk_present]:
        row = signature.loc[gene]
        missing.append({'candidate': mode, 'gene_symbol': gene, 'highest_signature_type': row.idxmax(), 'highest_signature_expression': float(row.max())})
    if mode == 'gse_joint':
        shutil.copyfile(path, upload / 'GJ20_Signature20.txt')
        joint_signature = signature
pd.DataFrame(coverage).to_csv(audit / 'signature_MSCCR_gene_coverage.tsv', sep='\t', index=False)
pd.DataFrame(type_coverage).to_csv(audit / 'signature_MSCCR_signal_coverage.tsv', sep='\t', index=False)
pd.DataFrame(missing).to_csv(audit / 'signature_MSCCR_missing_genes.tsv', sep='\t', index=False)
bulk_summary['signature_coverage'] = coverage
print('三套 signature 覆盖:', pd.DataFrame(coverage).to_string(index=False))

# Step3. 固定完整联合参考背景与 MSCCR 的交集，不按 bulk 高表达排名截取基因
source_reference = project / 'IBDverse_reference/three_reference_pilots_20261006/run_20261006_183211/gse_joint/GSE282122_Joint20_scRNA_reference.txt'
with source_reference.open(encoding='utf-8') as handle:
    next(handle)
    background = [line.split('\t', 1)[0] for line in handle]
common = [gene for gene in background if gene in available]
assert len(common) == len(set(common)) and set(joint_signature.index) & available <= set(common)
pd.DataFrame({'gene_symbol': common}).to_csv(audit / 'GJ20_MSCCR_common_background_gene_list.tsv', sep='\t', index=False)
bulk_summary['joint_reference_common_background_genes'] = len(common)
bulk_summary['joint_reference_common_background_all_samples_data_points'] = len(common) * bulk_summary['source_samples']
bulk_summary['common_background_all_samples_fits_web_data_point_limit'] = len(common) * bulk_summary['source_samples'] <= 50_000_000
bulk_summary['gene_overlap_interpretation'] = 'exact source symbols; no automatic aliases or arbitrary duplicate locus choice'
bulk_summary['downstream_expression_scope_note'] = 'using this background preserves fraction markers present in MSCCR, but restricts later recovered expression to the input genes; EcoTyper gene coverage requires separate audit'
(audit / 'summary.json').write_text(json.dumps(bulk_summary, ensure_ascii=False, indent=2), encoding='utf-8')
print('完整样本共同背景:', len(common), '基因 x', bulk_summary['source_samples'], '样本 =', len(common) * bulk_summary['source_samples'], '数据点')

# Step4. 拷贝可直接上传的两份 mixture；真实比例在本地保留，不是 expression 输入
for name in ['GJ20_Holdout28_Mixture_CPM.txt', 'GJ20_Holdout8_EqualRNA_Mixture_CPM.txt']:
    shutil.copyfile(validation / name, upload / name)
    assert hashlib.sha256((upload / name).read_bytes()).hexdigest() == validation_summary['files'][name]['sha256']
files = []
for file in sorted(upload.glob('*.txt')):
    content = file.read_bytes()
    files.append({'file': file.name, 'bytes': len(content), 'sha256': hashlib.sha256(content).hexdigest()})
(upload / 'upload_file_integrity.json').write_text(json.dumps(files, indent=2), encoding='utf-8')

# 有针对性地列出标志基因覆盖，缺失不等于生物学不存在，也不能据此直接改标签。
markers = ['DEFA5', 'DEFA6', 'PLA2G2A', 'ITLN2', 'SOX10', 'PLP1', 'S100B', 'SOX2', 'POU2F3', 'TRPM5', 'AVIL', 'SH2D6',
    'SI', 'ALPI', 'FABP2', 'APOA4', 'CA1', 'CA2', 'AQP8', 'SLC26A3', 'CLEC9A', 'CD1C', 'CLEC10A', 'CD40LG', 'MT-CO1', 'MT-CO3']
marker_rows = []
for gene in markers:
    marker_rows.append({'gene_symbol': gene, 'in_joint_signature': gene in joint_signature.index, 'in_MSCCR_unique_symbols': gene in available,
        'nonzero_somewhere_MSCCR': gene in nonzero, 'in_holdout_mixture': gene in mixture_genes})
pd.DataFrame(marker_rows).to_csv(audit / 'key_marker_coverage.tsv', sep='\t', index=False)
print('关键 marker 覆盖:', pd.DataFrame(marker_rows).to_string(index=False))
