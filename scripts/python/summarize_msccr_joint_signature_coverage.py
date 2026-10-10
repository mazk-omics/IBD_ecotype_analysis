# 汇总 MSCCR 网页输入对联合 signature 的基因覆盖，不据此宣称反卷积准确。
# 在本机项目环境运行；读取已回传身份报告和未改数值的 signature，写入本轮审计目录。
# 分开记录可匹配行、在 IBD 样本有表达的行，以及各类型的参考表达信号覆盖。
from pathlib import Path

import numpy as np
import pandas as pd

# Step1. 读取冻结的 signature 和本轮导出身份表
root = Path('D:/IBD_ecotype/IBD_ecotype_analysis/CIBERSORTx_signature_audit')
out = root / 'MSCCR_web_input_adaptation_20261006'
signature_path = root / 'main_candidate_validation_20261006/web_upload/GJ20_Signature20.txt'
signature = pd.read_csv(signature_path, sep='\t', index_col=0)
identity = pd.read_csv(out / 'server_outputs/audit/msccr_ibd_joint_web_20261006/exported_gene_identity.tsv', sep='\t')
bridge = pd.read_csv(out / 'server_outputs/audit/msccr_gene_identity_20261006/joint_reference_to_MSCCR_gene_bridge.tsv', sep='\t')
assignment = pd.read_csv(root / 'three_candidate_comparison_20261006/gse_joint/signature_gene_diagnostic_assignment.tsv', sep='\t')
assert signature.shape == (5816, 20) and signature.index.is_unique
mapped = signature.index.isin(identity.gene_symbol)
nonzero = signature.index.isin(identity.loc[~identity.all_zero_in_IBD, 'gene_symbol'])
print('signature:', signature.shape, '映射:', int(mapped.sum()), 'IBD 非零:', int(nonzero.sum()))

# Step2. 计算每列参考表达信号的覆盖，辅助发现少量关键基因缺失
# 这里的表达信号来自 signature 列求和，不是 MSCCR 细胞比例或统计准确率。
coverage = pd.DataFrame({'cell_type': signature.columns, 'reference_signal_sum': signature.sum(axis=0).to_numpy(),
    'mapped_signal_fraction': (signature.loc[mapped].sum(axis=0) / signature.sum(axis=0)).to_numpy(),
    'nonzero_IBD_signal_fraction': (signature.loc[nonzero].sum(axis=0) / signature.sum(axis=0)).to_numpy()})
diagnostic = assignment.merge(bridge[['gene_symbol', 'accepted', 'method']], left_on='gene', right_on='gene_symbol', how='left', validate='one_to_one')
assert len(diagnostic) == len(signature) and diagnostic.accepted.notna().all()
diagnostic['IBD_nonzero'] = diagnostic.gene.isin(identity.loc[~identity.all_zero_in_IBD, 'gene_symbol'])
diagnostic_counts = diagnostic.groupby('highest_expression_type').agg(diagnostic_genes=('gene', 'size'), mapped_diagnostic_genes=('accepted', 'sum'), nonzero_diagnostic_genes=('IBD_nonzero', 'sum')).reset_index()
coverage = coverage.merge(diagnostic_counts, left_on='cell_type', right_on='highest_expression_type', validate='one_to_one').drop(columns='highest_expression_type')
coverage.to_csv(out / 'joint_signature_MSCCR_coverage_by_type.tsv', sep='\t', index=False)
print(coverage[['cell_type', 'mapped_signal_fraction', 'nonzero_IBD_signal_fraction', 'diagnostic_genes', 'mapped_diagnostic_genes']].sort_values('nonzero_IBD_signal_fraction').to_string(index=False))

# Step3. 保留未映射和全零 signature 清单，便于后续按实际预测问题回查
missing = diagnostic.loc[~diagnostic.accepted].copy()
missing.to_csv(out / 'unmapped_signature_diagnostic_genes.tsv', sep='\t', index=False)
zero_signature = identity.loc[identity.signature_gene & identity.all_zero_in_IBD].copy()
zero_signature.to_csv(out / 'mapped_but_zero_IBD_signature_genes.tsv', sep='\t', index=False)
print('未映射 signature:', len(missing), '映射但全零 signature:', len(zero_signature))
