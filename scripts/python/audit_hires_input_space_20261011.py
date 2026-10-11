# 核对 Job24 原始 mixture 与 Mixture-Adjusted 的数值空间，避免把上游导出误当成 HiRes 实际输入。
# 本机 cibersortx_review Python 运行；读取原件并保存小基准及差异记录，不改变研究数据。
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from scipy.optimize import nnls

# Step1. 同时保留两个输入候选，比较前冻结基因、样本及类型
project = Path(__file__).resolve().parents[2]
audit = project / 'CIBERSORTx_server_plan_20261011'
prior = audit / 'core_benchmark_inputs'
originals = project / 'CIBERSORTx_signature_audit/MSCCR_Joint20v2_Job24_HiRes_review_20261009/original_outputs'
source = project / 'CIBERSORTx_ready_20261008/04_MSCCR120_HiRes_gene_batches/MSCCR_Platform120_CPM.txt'
outdir = audit / 'input_space_audit'
outdir.mkdir(exist_ok=True)
raw = pd.read_csv(source, sep='\t', index_col=0)
adjusted = pd.read_csv(originals / 'CIBERSORTxGEP_Job24_Mixture-Adjusted.txt', sep='\t', index_col=0)
fractions = pd.read_csv(prior / 'fractions_fixed.tsv', sep='\t', index_col=0)
group = pd.read_csv(prior / 'official_group_24genes.tsv', sep='\t', index_col=0)
genes = group.index.tolist()
assert raw.index.is_unique and raw.columns.tolist() == fractions.index.tolist()
assert set(genes).issubset(raw.index) and adjusted.columns.tolist() == raw.columns.tolist()
raw.loc[genes].to_csv(outdir / 'mixture_original_24genes.tsv', sep='\t', index_label='GeneSymbol')

# Step2. 只做可追踪的候选比较，不按最小误差偷偷选择输入或调参
candidates = {'original_linear': raw.loc[genes].to_numpy(), 'original_floor_one': np.maximum(raw.loc[genes].to_numpy(), 1), 'adjusted_export': adjusted.loc[genes].to_numpy()}
rows = []
predictions = {}
for name, values in candidates.items():
    predicted = np.vstack([nnls(fractions.to_numpy(), y)[0] for y in values])
    predictions[name] = predicted
    error = np.abs(np.log2(predicted + 1) - np.log2(group.to_numpy() + 1))
    rows.append({'input': name, 'genes': len(genes), 'median_type_log2_mae_to_official_group': float(np.median(error.mean(axis=0))), 'all_value_log2_mae_to_official_group': float(error.mean()), 'max_input_value': float(values.max()), 'min_input_value': float(values.min())})
comparison = pd.DataFrame(rows)
comparison.to_csv(outdir / 'input_candidates_vs_official_group.tsv', sep='\t', index=False)
common = raw.index.intersection(adjusted.index)
missing = raw.index.difference(adjusted.index)
pd.DataFrame({'GeneSymbol': missing}).to_csv(outdir / 'genes_absent_from_adjusted_export.tsv', sep='\t', index=False)
input_rows = []
for gene in genes:
    original = raw.loc[gene].to_numpy()
    corrected = adjusted.loc[gene].to_numpy()
    input_rows.append({'GeneSymbol': gene, 'original_mean': float(original.mean()), 'adjusted_export_mean': float(corrected.mean()), 'max_absolute_difference': float(np.max(np.abs(original - corrected))), 'log2_mae_between_inputs': float(np.mean(np.abs(np.log2(original + 1) - np.log2(corrected + 1))))})
pd.DataFrame(input_rows).to_csv(outdir / 'gene_input_difference.tsv', sep='\t', index=False)
provenance = {'original_shape': list(raw.shape), 'adjusted_export_shape': list(adjusted.shape), 'common_genes': len(common), 'original_genes_absent_from_export': len(missing), 'original_source': str(source), 'original_sha256': hashlib.sha256(source.read_bytes()).hexdigest(), 'adjusted_source_sha256': hashlib.sha256((originals / 'CIBERSORTxGEP_Job24_Mixture-Adjusted.txt').read_bytes()).hexdigest(), 'conclusion_limit': 'Candidate numerical comparison does not by itself prove current website input routing; use scientific dispatcher and original runtime log evidence together'}
(outdir / 'input_space_provenance.json').write_text(json.dumps(provenance, ensure_ascii=False, indent=2), encoding='utf-8')
print(comparison.to_string(index=False))
print(json.dumps(provenance, ensure_ascii=False, indent=2))
