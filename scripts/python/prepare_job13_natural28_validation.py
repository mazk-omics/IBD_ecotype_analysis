# 整理 Job13 新版 signature 与原 28 个自然混合物，供 CIBERSORTx 网页与 Job6 配对验证。
# 在本机 cibersortx_review 环境运行；读取已冻结文件，输出两个松散 TXT 和本地核对记录，不访问服务器。
# 自然混合物由细胞原始 counts 求和后转 CPM；分别保留 RNA 贡献真值与细胞数真值，不用预测值定义真值。
from pathlib import Path
import hashlib
import json
import shutil
import numpy as np
import pandas as pd

# Step1. 读入冻结输入，新目录与等 RNA 任务分开
root = Path('D:/IBD_ecotype/IBD_ecotype_analysis')
base = root / 'CIBERSORTx_signature_audit/main_candidate_validation_20261006'
frozen = base / 'server_outputs/validation_joint20_20261006'
out = root / 'CIBERSORTx_ready_20261007/08_Joint20v2_Natural28_validation'
audit = out / '核对记录_勿上传'
audit.mkdir(parents=True, exist_ok=True)
(out / '网页结果').mkdir(exist_ok=True)
signature_source = root / 'signature_matrix/gse_joint_v2/CIBERSORTx_Job13_GJ20v2_PhenotypeClasses.CIBERSORTx_Job13_GJ20v2_DonorSite_CPM.bm.K999.txt'
mixture_source = frozen / 'GJ20_Holdout28_Mixture_CPM.txt'
signature = pd.read_csv(signature_source, sep='\t', index_col=0)
mixture = pd.read_csv(mixture_source, sep='\t', index_col=0)
truth_rna = pd.read_csv(frozen / 'ground_truth_export_UMI_fractions.tsv', sep='\t', index_col=0)
truth_cells = pd.read_csv(frozen / 'ground_truth_cell_fractions.tsv', sep='\t', index_col=0)
baseline = pd.read_csv(base / 'Job6_review/predicted_fractions_aligned.tsv', sep='\t', index_col=0)
assert signature.shape == (4474, 20) and mixture.shape == (23633, 28)
assert signature.index.is_unique and mixture.index.is_unique and mixture.columns.is_unique
assert set(signature.columns) == set(truth_rna.columns) == set(truth_cells.columns) == set(baseline.columns)
assert set(mixture.columns) == set(truth_rna.index) == set(truth_cells.index) == set(baseline.index)
truth_rna = truth_rna.loc[mixture.columns, signature.columns]
truth_cells = truth_cells.loc[mixture.columns, signature.columns]
baseline = baseline.loc[mixture.columns, signature.columns]
assert np.allclose(truth_rna.sum(axis=1), 1) and np.allclose(truth_cells.sum(axis=1), 1) and np.allclose(baseline.sum(axis=1), 1)
assert np.isfinite(signature.to_numpy()).all() and signature.ge(0).all().all()
assert np.isfinite(mixture.to_numpy()).all() and mixture.ge(0).all().all() and mixture.to_numpy().max() > 50
assert np.allclose(mixture.sum(axis=0), 1e6, rtol=0, atol=1) and mixture.size < 50_000_000
frozen_summary = json.loads((frozen / 'summary.json').read_text(encoding='utf-8'))
assert hashlib.sha256(mixture_source.read_bytes()).hexdigest() == frozen_summary['files'][mixture_source.name]['sha256']
print('Signature genes × types:', signature.shape, '; mixture genes × samples:', mixture.shape, flush=True)

# Step2. 保留所有原 mixture 基因；新版 signature 的覆盖单独核对
coverage = pd.DataFrame({'signature_gene': signature.index})
coverage['in_original_mixture'] = coverage.signature_gene.isin(mixture.index)
coverage.to_csv(audit / 'signature_gene_overlap.tsv', sep='\t', index=False)
overlap = int(coverage.in_original_mixture.sum())
assert overlap == len(signature)
print('Signature gene overlap:', overlap, '/', len(signature), flush=True)

# Step3. 从原细胞清单重算两种真值，核对 v2 训练与最终留出隔离
membership = pd.read_csv(frozen / 'mixture_cell_membership.tsv', sep='\t')
test_cells = membership.loc[membership['sample'].isin(mixture.columns)].copy()
assert not test_cells.duplicated(['sample', 'cell_id']).any()
assert test_cells.groupby('sample').size().eq(200).all() and test_cells['sample'].nunique() == 28
recounted_cells = pd.crosstab(test_cells['sample'], test_cells.pilot_type).reindex(index=mixture.columns, columns=signature.columns, fill_value=0) / 200
assert np.allclose(recounted_cells.to_numpy(), truth_cells.to_numpy())
rna_counts = test_cells.groupby(['sample', 'pilot_type']).export_library_counts.sum().unstack(fill_value=0)
rna_counts = rna_counts.reindex(index=mixture.columns, columns=signature.columns, fill_value=0)
recounted_rna = rna_counts.div(rna_counts.sum(axis=1), axis=0)
assert np.allclose(recounted_rna.to_numpy(), truth_rna.to_numpy())
training = pd.read_csv(root / 'CIBERSORTx_ready_20261007/06_Joint20v2_DonorSite_reference/核对记录_勿上传/selected_cell_manifest.tsv', sep='\t')
assert not test_cells.Patient.isin(training.Patient).any() and not test_cells.cell_id.isin(training.cell_id).any()
assert not test_cells.Patient.isin(['CD12', 'UC17']).any()
design = pd.read_csv(frozen / 'mixture_design.tsv', sep='\t').set_index('sample').loc[mixture.columns].reset_index()
print(design.groupby(['region', 'scenario']).size().to_string(), flush=True)
test_cells.to_csv(audit / 'frozen_natural28_cell_membership.tsv', sep='\t', index=False)
truth_rna.to_csv(audit / 'ground_truth_export_UMI_fractions.tsv', sep='\t', index_label='sample')
truth_cells.to_csv(audit / 'ground_truth_cell_fractions.tsv', sep='\t', index_label='sample')
baseline.to_csv(audit / 'Job6_same_samples_predictions.tsv', sep='\t', index_label='sample')
design.to_csv(audit / 'natural28_design.tsv', sep='\t', index=False)

# Step4. 原样复制两个网页 TXT；循环只用于两文件的同版导出及逐行格式检查
checks = []
for source, name, n_columns in [(signature_source, 'GJ20v2_Signature20.txt', 21), (mixture_source, mixture_source.name, 29)]:
    destination = out / name
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    if destination.exists():
        assert hashlib.sha256(destination.read_bytes()).hexdigest() == digest, '不覆盖已有不同内容的网页输入'
    else:
        shutil.copyfile(source, destination)
    assert hashlib.sha256(destination.read_bytes()).hexdigest() == digest
    with destination.open(encoding='utf-8') as handle:
        for line_number, line in enumerate(handle):
            fields = line.rstrip('\r\n').split('\t')
            assert len(fields) == n_columns and all(fields) and '"' not in line, (name, line_number)
    checks.append({'source': str(source), 'file': name, 'bytes': destination.stat().st_size, 'sha256': digest, 'identical_to_source': True})
pd.DataFrame(checks).to_csv(audit / 'web_input_checksums.tsv', sep='\t', index=False)

# Step5. 保存提交设置与准备状态；参数沿用 Job6 未校正对照，不在本机提交任务
settings = {'module': 'Impute Cell Fractions', 'analysis_mode': 'Custom', 'signature': 'GJ20v2_Signature20.txt', 'mixture': mixture_source.name, 'disable_quantile_normalization': True, 'batch_correction': False, 'run_mode': 'relative', 'permutations': 100, 'suggested_job_label': 'GJ20v2_N28_NS'}
summary = {'preparation_complete': True, 'web_job_submitted': False, 'signature_shape_genes_by_types': list(signature.shape), 'mixture_shape_genes_by_samples': list(mixture.shape), 'mixture_data_points': mixture.size, 'signature_gene_overlap': overlap, 'signature_gene_overlap_fraction': overlap / len(signature), 'mixture_identical_to_frozen_original': True, 'RNA_and_cell_truth_recount_matches': True, 'training_donor_overlap': 0, 'training_cell_overlap': 0, 'reserved_final_test_donors_used': False, 'unique_test_cells': int(test_cells.cell_id.nunique()), 'test_donors': sorted(test_cells.Patient.unique().tolist()), 'comparison_baseline': 'Job6, identical frozen natural28 input', 'recommended_web_settings': settings}
(audit / 'preparation_summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)
