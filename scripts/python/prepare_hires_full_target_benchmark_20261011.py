# 将冻结比较扩展到Job24全部目标基因，保留原始CPM与固定fraction，不按误差挑基因。
# 本机科学Python环境运行；读取原件，仅在新目录保存研究输入和来源，不修改正式分析。
from pathlib import Path
import hashlib
import json
import shutil
import pandas as pd

# Step1. 固定官方任务的目标基因集合，表达仍从原始mixture提取
project = Path(__file__).resolve().parents[2]
audit = project / 'CIBERSORTx_server_plan_20261011'
outdir = audit / 'full_target_benchmark_inputs'
assert not outdir.exists()
official_dir = project / 'CIBERSORTx_signature_audit/MSCCR_Joint20v2_Job24_HiRes_review_20261009/original_outputs'
mixture_path = project / 'CIBERSORTx_ready_20261008/04_MSCCR120_HiRes_gene_batches/MSCCR_Platform120_CPM.txt'
fraction_path = audit / 'core_benchmark_inputs/fractions_fixed.tsv'
fractions = pd.read_csv(fraction_path, sep='\t', index_col=0)
mixture = pd.read_csv(mixture_path, sep='\t', index_col=0)
assert mixture.index.is_unique and mixture.columns.tolist() == fractions.index.tolist()
first_path = official_dir / f'CIBERSORTxHiRes_Job24_{fractions.columns[0]}_Window60.txt'
first = pd.read_csv(first_path, sep='\t', index_col=0)
official_genes = first.index.tolist()
assert len(official_genes) == 1249 and len(set(official_genes)) == len(official_genes)
# 网页输出把部分基因名中的连字符改成点；只在无歧义时建立对应，原始输入保持官方gene symbol。
name_map = pd.DataFrame({'GeneSymbol': mixture.index, 'official_output_gene': mixture.index.str.replace('-', '.', regex=False)})
target_map = pd.DataFrame({'official_output_gene': official_genes}).merge(name_map, on='official_output_gene', how='left', sort=False, validate='one_to_one')
assert target_map['GeneSymbol'].notna().all() and target_map['GeneSymbol'].is_unique
genes = target_map['GeneSymbol'].tolist()
assert first.columns.tolist() == mixture.columns.tolist()
for celltype in fractions.columns:
    output = pd.read_csv(official_dir / f'CIBERSORTxHiRes_Job24_{celltype}_Window60.txt', sep='\t', index_col=0)
    assert output.index.tolist() == official_genes and output.columns.tolist() == mixture.columns.tolist()

# Step2. 全基因入口统计已固定；种子仍依赖原始mixture行号，扩大面板不改变原24基因随机流
outdir.mkdir()
panel = mixture.loc[genes].copy()
assert panel.notna().all().all() and (panel >= 0).all().all()
selection = pd.DataFrame({'GeneSymbol': genes, 'full_mixture_row_1based': mixture.index.get_indexer(genes) + 1})
selection.to_csv(outdir / 'gene_selection.tsv', sep='\t', index=False)
target_map.to_csv(outdir / 'gene_name_mapping.tsv', sep='\t', index=False)
panel.to_csv(outdir / 'mixture_original_target_genes.tsv', sep='\t', index_label='GeneSymbol')
shutil.copyfile(fraction_path, outdir / 'fractions_fixed.tsv')
source_paths = [mixture_path, fraction_path, first_path, audit / 'entry_statistics_original_run1/entry_statistics.rds']
sources = [{'path': str(path.resolve()), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()} for path in source_paths]
manifest = {'genes': len(genes), 'samples': len(mixture.columns), 'types': len(fractions.columns), 'selection': 'All 1249 target genes found in the existing Job24 output; no selection by agreement.', 'gene_name_aliases': int((target_map['GeneSymbol'] != target_map['official_output_gene']).sum()), 'expression_input': 'Original linear CPM, not official HiRes or Mixture-Adjusted values.', 'entry_statistics': 'Reuse complete 17040-gene entry_statistics_original_run1; no panel-wise BH recalculation.', 'sources': sources, 'upstream_dependence': 'Fractions remain official fixed fractions; independent signature, fractions and S-mode not implemented.'}
(outdir / 'input_provenance.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding='utf-8')
print(manifest['genes'], 'genes;', manifest['samples'], 'samples;', manifest['types'], 'types')
print('Research benchmark inputs:', outdir)
