# 将冻结的120例HiRes基因请求拆分为网页可接受的小清单，不改变表达输入。
# 本机cibersortx_review环境运行；读取03目录，写入04目录。前两批用于运行及共同基因一致性检查，未提交网页任务。
from pathlib import Path
import hashlib
import json
import shutil
import numpy as np
import pandas as pd

# Step1. 读入冻结输入，确认样本与完整基因背景没有改变
project = Path('D:/IBD_ecotype/IBD_ecotype_analysis')
source = project / 'CIBERSORTx_ready_20261008/03_MSCCR120_HiRes_next'
outdir = project / 'CIBERSORTx_ready_20261008/04_MSCCR120_HiRes_gene_batches'
first_dir = outdir / '先提交这两批'
later_dir = outdir / '其余批次_先勿提交'
records = outdir / '核对记录_勿上传'
first_dir.mkdir(parents=True, exist_ok=True)
later_dir.mkdir(parents=True, exist_ok=True)
records.mkdir(parents=True, exist_ok=True)
checks = pd.read_csv(source / '核对记录_勿上传/上传文件维度与哈希.tsv', sep='\t')
annotation = pd.read_csv(source / '核对记录_勿上传/MSCCR120_annotation.tsv', sep='\t')
mixture = pd.read_csv(source / 'MSCCR_Platform120_CPM.txt', sep='\t', index_col=0)
signature = pd.read_csv(source / 'GJ20v2_Signature20.txt', sep='\t', index_col=0)
assert mixture.shape == (17040, 120) and signature.shape == (4474, 20)
assert mixture.columns.tolist() == annotation['ID'].tolist() and annotation['PatientID'].nunique() == 120
assert mixture.index.is_unique and signature.index.is_unique
print('冻结mixture、signature维度：', mixture.shape, signature.shape)

# 同版多文件复制必须逐一核对哈希；不对三个矩阵删行、缩放或改列顺序。
for filename in ['GJ20v2_Signature20.txt', 'MSCCR_Platform120_CPM.txt', 'GJ20v2_SourceSubset1060_scRNA_CPM.txt']:
    input_path = source / filename
    expected_hash = checks.loc[checks['file'].eq(filename), 'sha256'].iloc[0]
    assert hashlib.sha256(input_path.read_bytes()).hexdigest() == expected_hash
    output_path = outdir / filename
    if output_path.exists():
        assert hashlib.sha256(output_path.read_bytes()).hexdigest() == expected_hash
    else:
        shutil.copy2(input_path, output_path)
    assert hashlib.sha256(output_path.read_bytes()).hexdigest() == expected_hash
annotation.to_csv(records / 'MSCCR120_annotation.tsv', sep='\t', index=False)

# Step2. 选择100个共同检查基因，每类5个，不作为新signature或验证真值
# 从现有signature与bulk的交集中，按该类型相对其他类型平均表达的比值排序。
# 只用于跨任务比较同一基因的表达；本项目选择5/类，不是CIBERSORTx默认或显著DEG判定。
shared_signature = signature.loc[signature.index.intersection(mixture.index)].copy()
anchor_rows = []
anchor_genes = []
for cell_type in signature.columns:
    other_mean = shared_signature.drop(columns=cell_type).mean(axis=1)
    ratio = shared_signature[cell_type] / (other_mean + 1)
    candidates = pd.DataFrame({'type_expression': shared_signature[cell_type], 'ratio': ratio})
    candidates = candidates.loc[candidates['type_expression'].gt(10) & ~candidates.index.isin(anchor_genes)]
    candidates = candidates.sort_values(['ratio', 'type_expression'], ascending=False, kind='stable').head(5)
    assert len(candidates) == 5, f'共同检查基因不足5个：{cell_type}'
    anchor_genes.extend(candidates.index.tolist())
    for gene, row in candidates.iterrows():
        anchor_rows.append({'GeneSymbol': gene, 'selection_type': cell_type, 'signature_expression': row['type_expression'], 'type_over_other_mean_plus_1': row['ratio']})
assert len(anchor_genes) == 100 and len(set(anchor_genes)) == 100
pd.DataFrame(anchor_rows).to_csv(records / '共同100基因_选择依据.tsv', sep='\t', index=False)
print('共同检查基因数：', len(anchor_genes))

# Step3. 拆分全部17040基因，只改变请求顺序，不筛掉低变异基因
# 较高bulk变异的基因排在前面，便于前两批查看样本级表达；所有剩余基因均进入后续批次。
gene_variance = np.log1p(mixture).var(axis=1, ddof=1).sort_values(ascending=False, kind='stable')
remaining_genes = gene_variance.index[~gene_variance.index.isin(anchor_genes)].tolist()
max_genes_per_request = 1249  # 网页报错硬上限1250；另留1行表头，避免行数计数歧义。
new_genes_per_batch = max_genes_per_request - len(anchor_genes)
batch_rows = []
membership_rows = []
covered_genes = set()
for start in range(0, len(remaining_genes), new_genes_per_batch):
    batch_number = start // new_genes_per_batch + 1
    new_genes = remaining_genes[start:start + new_genes_per_batch]
    request_genes = anchor_genes + new_genes
    assert len(request_genes) <= max_genes_per_request and len(set(request_genes)) == len(request_genes)
    batch_dir = first_dir if batch_number <= 2 else later_dir
    gene_path = batch_dir / f'MSCCR120_HiRes_B{batch_number:02d}_Genes{len(request_genes)}.txt'
    gene_text = 'GeneSymbol\n' + '\n'.join(request_genes) + '\n'
    if gene_path.exists():
        assert gene_path.read_text(encoding='utf-8') == gene_text, f'已有批次清单不同：{gene_path}'
    else:
        gene_path.write_text(gene_text, encoding='utf-8')
    parsed_genes = pd.read_csv(gene_path)['GeneSymbol'].tolist()
    assert parsed_genes == request_genes and set(parsed_genes).issubset(mixture.index)
    covered_genes.update(request_genes)
    batch_rows.append({'batch': batch_number, 'file': gene_path.relative_to(outdir).as_posix(), 'requested_genes': len(request_genes), 'new_genes': len(new_genes), 'shared_check_genes': len(anchor_genes), 'lines_including_header': len(request_genes) + 1, 'sha256': hashlib.sha256(gene_path.read_bytes()).hexdigest(), 'stage': 'first_two_checks' if batch_number <= 2 else 'not_for_submission_yet'})
    membership_rows.extend({'GeneSymbol': gene, 'batch': batch_number, 'role': 'shared_check' if gene in anchor_genes else 'new_gene'} for gene in request_genes)
assert covered_genes == set(mixture.index) and sum(row['new_genes'] for row in batch_rows) == len(remaining_genes)
assert len(batch_rows) == 15
pd.DataFrame(batch_rows).to_csv(records / '批次清单与哈希.tsv', sep='\t', index=False)
pd.DataFrame(membership_rows).to_csv(records / '基因批次对应.tsv', sep='\t', index=False)
gene_variance.rename('variance_log1p_CPM').to_csv(records / '基因请求排序依据.tsv', sep='\t', index_label='GeneSymbol')

# Step4. 保存限制证据与待验证事项，不能把本地分批当成已成功运行
screenshot = Path('C:/Users/马泽凯/AppData/Local/Temp/codex-clipboard-2c5b3ef8-8a87-4c0f-b157-478e26d85879.png')
shutil.copy2(screenshot, records / '网页1250基因限制_用户截图.png')
summary = {'status': 'prepared_not_submitted', 'date': '2026-10-09', 'fixed_samples': 120, 'types': 20, 'web_reported_gene_limit': 1250, 'screenshot_detected_previous_list_genes': 17039, 'previous_list_actual_genes': 17040, 'header_policy': 'GeneSymbol header based on apparent first-line consumption; confirm actual web parsing/output', 'first_line_consumption_confirmed': False, 'max_requested_genes_per_batch': max_genes_per_request, 'shared_check_genes': 100, 'batches': len(batch_rows), 'unique_gene_coverage': len(covered_genes), 'batching_equivalent_to_full_transcriptome_verified': False, 'cross_batch_fractions_adjustment_and_shared_expression_verified': False, 'first_two_batches_only_for_initial_submission': True}
(records / '分批准备状态.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
print(pd.DataFrame(batch_rows)[['batch', 'requested_genes', 'new_genes', 'stage']].to_string(index=False))
print('独立目录：', outdir)
print('全部17040基因已覆盖。先顺序提交01、02批；未经一致性检查，不直接拼接表达用于EcoTyper。')
