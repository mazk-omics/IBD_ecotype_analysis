# 修正旧24基因面板的种子索引来源，并与1249基因运行使用同一原始mixture行号。
# 本机科学Python环境运行；历史清单和结果保留，仅在新目录写修正对照与索引差异。
from pathlib import Path
import hashlib
import json
import shutil
import pandas as pd

# Step1. 原始表达不变，单独纠正继承自Adjusted表的旧行号
project = Path(__file__).resolve().parents[2]
audit = project / 'CIBERSORTx_server_plan_20261011'
outdir = audit / 'canonical_24_benchmark_inputs'
assert not outdir.exists()
legacy = pd.read_csv(audit / 'core_benchmark_inputs/gene_selection.tsv', sep='\t')
canonical = pd.read_csv(audit / 'full_target_benchmark_inputs/gene_selection.tsv', sep='\t').set_index('GeneSymbol').loc[legacy['GeneSymbol']].reset_index()
expression = pd.read_csv(audit / 'full_target_benchmark_inputs/mixture_original_target_genes.tsv', sep='\t', index_col=0).loc[legacy['GeneSymbol']]
previous_expression = pd.read_csv(audit / 'input_space_audit/mixture_original_24genes.tsv', sep='\t', index_col=0)
assert expression.equals(previous_expression)
comparison = legacy[['GeneSymbol', 'full_mixture_row_1based']].rename(columns={'full_mixture_row_1based': 'legacy_adjusted_export_row_1based'})
comparison['canonical_original_mixture_row_1based'] = canonical['full_mixture_row_1based']
comparison['seed_index_changed'] = comparison['legacy_adjusted_export_row_1based'] != comparison['canonical_original_mixture_row_1based']

# Step2. 两种面板读取同一完整入口统计；只改变旧worker/refinement的种子索引
outdir.mkdir()
canonical.to_csv(outdir / 'gene_selection.tsv', sep='\t', index=False)
expression.to_csv(outdir / 'mixture_original_24genes.tsv', sep='\t', index_label='GeneSymbol')
shutil.copyfile(audit / 'full_target_benchmark_inputs/fractions_fixed.tsv', outdir / 'fractions_fixed.tsv')
mapping = pd.read_csv(audit / 'full_target_benchmark_inputs/gene_name_mapping.tsv', sep='\t').set_index('GeneSymbol').loc[legacy['GeneSymbol']].reset_index()
mapping.to_csv(outdir / 'gene_name_mapping.tsv', sep='\t', index=False)
comparison.to_csv(outdir / 'seed_index_correction.tsv', sep='\t', index=False)
manifest = {'genes': 24, 'expression_unchanged': True, 'seed_indices_changed': int(comparison['seed_index_changed'].sum()), 'correction': 'Legacy benchmark indices came from Mixture-Adjusted export; canonical indices come from original full mixture. Seed constants and all scientific thresholds unchanged.', 'entry_statistics': 'Same complete 17040-gene entry_statistics_original_run1; no entry recomputation.', 'selection_sha256': hashlib.sha256((outdir / 'gene_selection.tsv').read_bytes()).hexdigest(), 'limits': 'Independent R gene-specific seed convention, not official threaded RNG reproduction.'}
(outdir / 'input_provenance.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps(manifest, ensure_ascii=False, indent=2))
