# 为独立 R 数值原型准备固定的 Job24 小基准；只读已有官方输出，不改变研究输入。
# 本机 cibersortx_review Python 运行；输出 24 基因、全部 120 样本和 20 类型的小表及原件哈希。
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

# Step1. 读取同一官方任务的校正后输入，隔离上游平台校正与下游数值实现的差异
project = Path(__file__).resolve().parents[2]
originals = project / 'CIBERSORTx_signature_audit/MSCCR_Joint20v2_Job24_HiRes_review_20261009/original_outputs'
anchor_path = project / 'CIBERSORTx_ready_20261008/04_MSCCR120_HiRes_gene_batches/核对记录_勿上传/共同100基因_选择依据.tsv'
outdir = project / 'CIBERSORTx_server_plan_20261011/core_benchmark_inputs'
outdir.mkdir(exist_ok=True)
fractions = pd.read_csv(originals / 'CIBERSORTxGEP_Job24_Weights.txt', sep='\t', index_col=0)
mixture = pd.read_csv(originals / 'CIBERSORTxGEP_Job24_Mixture-Adjusted.txt', sep='\t', index_col=0)
group = pd.read_csv(originals / 'CIBERSORTxGEP_Job24_GEPs.txt', sep='\t', index_col=0)
fraction_output = pd.read_csv(originals / 'CIBERSORTxGEP_Job24_Fractions-Adjusted.txt', sep='\t', index_col=0)
assert fractions.index.is_unique and fractions.columns.is_unique and mixture.index.is_unique
assert mixture.columns.tolist() == fractions.index.tolist() and group.columns.tolist() == fractions.columns.tolist()
assert np.isfinite(fractions.to_numpy()).all() and fractions.ge(0).all().all()
assert np.allclose(fractions.sum(axis=1), 1, atol=1e-8)
print('完整调整后输入、比例及 group 输出：', mixture.shape, fractions.shape, group.shape)

# Step2. 在既有100个共同基因中等距取24个，不根据官方预测好坏筛选
anchors = pd.read_csv(anchor_path, sep='\t')
positions = np.linspace(0, len(anchors) - 1, 24, dtype=int)
selection = anchors.iloc[positions].copy()
genes = selection['GeneSymbol'].tolist()
assert len(genes) == len(set(genes)) == 24 and set(genes).issubset(mixture.index)
selection['full_mixture_row_1based'] = mixture.index.get_indexer(genes) + 1
selection.to_csv(outdir / 'gene_selection.tsv', sep='\t', index=False)
mixture.loc[genes].to_csv(outdir / 'mixture_adjusted_24genes.tsv', sep='\t', index_label='GeneSymbol')
fractions.to_csv(outdir / 'fractions_fixed.tsv', sep='\t', index_label='Sample')
group.loc[genes].to_csv(outdir / 'official_group_24genes.tsv', sep='\t', index_label='GeneSymbol')
print(selection[['GeneSymbol', 'selection_type']].to_string(index=False))

# Step3. 保存相同基因的官方 HiRes，作为差异对照而非用于拟合原型
source_paths = [anchor_path, originals / 'CIBERSORTxGEP_Job24_Weights.txt', originals / 'CIBERSORTxGEP_Job24_Mixture-Adjusted.txt', originals / 'CIBERSORTxGEP_Job24_GEPs.txt', originals / 'CIBERSORTxGEP_Job24_Fractions-Adjusted.txt']
for cell_type in fractions.columns:
    path = originals / f'CIBERSORTxHiRes_Job24_{cell_type}_Window60.txt'
    official = pd.read_csv(path, sep='\t', index_col=0)
    assert official.columns.tolist() == fractions.index.tolist() and set(genes).issubset(official.index)
    official.loc[genes].to_csv(outdir / f'official_hires_{cell_type}.tsv', sep='\t', index_label='GeneSymbol', na_rep='NA')
    source_paths.append(path)
fraction_difference = float(np.max(np.abs(fractions.to_numpy() - fraction_output.loc[fractions.index, fractions.columns].to_numpy())))
provenance = {'source_job': 'Job24', 'purpose': 'Independent numerical core benchmark; not a full CIBERSORTx replacement', 'genes': len(genes), 'samples': len(fractions), 'cell_types': len(fractions.columns), 'selection': '24 equally spaced positions in the previously frozen 100-gene anchor list', 'upstream': 'Official adjusted mixture and fixed group-mode weights; no independent fractions or S-mode implemented', 'max_weight_vs_fraction_difference': fraction_difference, 'adjusted_mixture_min': float(mixture.to_numpy().min()), 'sources': [{'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()} for path in source_paths]}
(outdir / 'input_provenance.json').write_text(json.dumps(provenance, ensure_ascii=False, indent=2), encoding='utf-8')
print('Weights 与官方 fraction 最大差：', fraction_difference)
print('输出：', outdir)
