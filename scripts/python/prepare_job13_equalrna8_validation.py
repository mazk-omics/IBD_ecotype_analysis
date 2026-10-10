# 整理 Job13 与原 8 个等 RNA 开发混合物，供用户在 CIBERSORTx 网页进行同样本验证。
# 在本机 cibersortx_review 环境运行；只复制已经冻结的输入，核对类型、基因覆盖和训练/留出隔离，不访问服务器。
# 本轮只改变 signature；等 RNA 真值来自原始细胞计数，不用新矩阵或预测值重新定义真值。
from pathlib import Path
import hashlib
import json
import shutil
import numpy as np
import pandas as pd

# Step1. 定位新版 signature 和原开发集，输出松散文件
root = Path('D:/IBD_ecotype/IBD_ecotype_analysis')
base = root / 'CIBERSORTx_signature_audit/main_candidate_validation_20261006'
frozen = base / 'server_outputs/validation_joint20_20261006'
out = root / 'CIBERSORTx_ready_20261007/07_Joint20v2_EqualRNA8_validation'
audit = out / '核对记录_勿上传'
audit.mkdir(parents=True, exist_ok=True)
signature_source = root / 'signature_matrix/gse_joint_v2/CIBERSORTx_Job13_GJ20v2_PhenotypeClasses.CIBERSORTx_Job13_GJ20v2_DonorSite_CPM.bm.K999.txt'
mixture_source = frozen / 'GJ20_Holdout8_EqualRNA_Mixture_CPM.txt'
signature = pd.read_csv(signature_source, sep='\t', index_col=0)
mixture = pd.read_csv(mixture_source, sep='\t', index_col=0)
truth = pd.read_csv(frozen / 'ground_truth_equalRNA_fractions.tsv', sep='\t', index_col=0)
baseline = pd.read_csv(base / 'Job7_review/predicted_fractions_aligned.tsv', sep='\t', index_col=0)
assert signature.shape == (4474, 20) and mixture.shape == (23633, 8)
assert signature.index.is_unique and mixture.index.is_unique and mixture.columns.is_unique
assert set(signature.columns) == set(truth.columns) == set(baseline.columns)
assert set(mixture.columns) == set(truth.index) == set(baseline.index)
truth = truth.loc[mixture.columns, signature.columns]
baseline = baseline.loc[mixture.columns, signature.columns]
assert np.allclose(truth.sum(axis=1), 1) and np.allclose(baseline.sum(axis=1), 1)
assert np.isfinite(signature.to_numpy()).all() and signature.ge(0).all().all()
assert np.isfinite(mixture.to_numpy()).all() and mixture.ge(0).all().all() and mixture.to_numpy().max() > 50
assert np.allclose(mixture.sum(axis=0), 1e6, rtol=0, atol=1)
assert mixture.size < 50_000_000
frozen_summary = json.loads((frozen / 'summary.json').read_text(encoding='utf-8'))
mixture_hash = hashlib.sha256(mixture_source.read_bytes()).hexdigest()
assert mixture_hash == frozen_summary['files'][mixture_source.name]['sha256']
print('Signature:', signature.shape, '; frozen mixture:', mixture.shape, flush=True)

# Step2. 保存基因覆盖，保持 mixture 的数值和完整基因背景不变
coverage = pd.DataFrame({'signature_gene': signature.index})
coverage['in_original_mixture'] = coverage.signature_gene.isin(mixture.index)
coverage.to_csv(audit / 'signature_gene_overlap.tsv', sep='\t', index=False)
overlap = int(coverage.in_original_mixture.sum())
print('Signature gene overlap:', overlap, '/', len(signature), flush=True)
assert overlap / len(signature) >= 0.95

# Step3. 重算原细胞计数真值，确认 v2 训练供者未进入这些混合物
membership = pd.read_csv(frozen / 'mixture_cell_membership.tsv', sep='\t')
sample_names = pd.Series(mixture.columns)
base_samples = sample_names.str.removesuffix('_EqualRNA')
sample_map = dict(zip(base_samples, sample_names))
test_cells = membership.loc[membership['sample'].isin(base_samples)].copy()
assert not test_cells.duplicated(['sample', 'cell_id']).any()
assert test_cells.groupby('sample').size().eq(200).all() and test_cells['sample'].nunique() == 8
test_cells['sample'] = test_cells['sample'].map(sample_map)
recounted = pd.crosstab(test_cells['sample'], test_cells.pilot_type).reindex(index=mixture.columns, columns=signature.columns, fill_value=0) / 200
assert np.allclose(recounted.to_numpy(), truth.to_numpy())
training = pd.read_csv(root / 'CIBERSORTx_ready_20261007/06_Joint20v2_DonorSite_reference/核对记录_勿上传/selected_cell_manifest.tsv', sep='\t')
assert not test_cells.Patient.isin(training.Patient).any() and not test_cells.cell_id.isin(training.cell_id).any()
assert not test_cells.Patient.isin(['CD12', 'UC17']).any()
test_cells.to_csv(audit / 'frozen_equalRNA8_cell_membership.tsv', sep='\t', index=False)
truth.to_csv(audit / 'ground_truth_equalRNA8_fractions.tsv', sep='\t', index_label='sample')
baseline.to_csv(audit / 'Job7_same_samples_predictions.tsv', sep='\t', index_label='sample')
design = pd.read_csv(frozen / 'mixture_design.tsv', sep='\t').set_index('sample').loc[base_samples].reset_index()
design['sample'] = design['sample'].map(sample_map)
design['unit'] = 'equal per-cell CPM contribution; original frozen mixture unchanged'
design.to_csv(audit / 'equalRNA8_design.tsv', sep='\t', index=False)

# Step4. 复制网页文件，不改 gene symbol、表达或样本顺序
checks = []
for source, name in [(signature_source, 'GJ20v2_Signature20.txt'), (mixture_source, 'GJ20_Holdout8_EqualRNA_Mixture_CPM.txt')]:
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
            expected_columns = len(signature.columns) + 1 if source == signature_source else len(mixture.columns) + 1
            assert len(fields) == expected_columns and all(fields) and '"' not in line, (name, line_number)
    checks.append({'source': str(source), 'file': name, 'bytes': destination.stat().st_size, 'sha256': digest, 'identical_to_source': True})
pd.DataFrame(checks).to_csv(audit / 'web_input_checksums.tsv', sep='\t', index=False)

# Step5. 保存当前准备状态；Job7 实际参数缺少日志，不能宣称已经全部核验
settings = {'module': 'Impute Cell Fractions', 'analysis_mode': 'Custom', 'signature': 'GJ20v2_Signature20.txt', 'mixture': mixture_source.name, 'disable_quantile_normalization': True, 'batch_correction': False, 'run_mode': 'relative', 'permutations': 100, 'suggested_job_label': 'GJ20v2_EQ8_NS'}
summary = {'preparation_complete': True, 'web_job_submitted': False, 'signature_shape_genes_by_types': list(signature.shape), 'mixture_shape_genes_by_samples': list(mixture.shape), 'mixture_data_points': mixture.size, 'signature_gene_overlap': overlap, 'signature_gene_overlap_fraction': overlap / len(signature), 'mixture_identical_to_frozen_original': True, 'truth_recount_matches': True, 'training_donor_overlap': 0, 'training_cell_overlap': 0, 'reserved_final_test_donors_used': False, 'unique_test_cells': test_cells.cell_id.nunique(), 'test_donors': sorted(test_cells.Patient.unique().tolist()), 'Job7_actual_web_settings_independently_verified': False, 'Job7_settings_note': 'No runtime log or settings screenshot received; new settings follow original prescribed protocol', 'recommended_web_settings': settings}
(audit / 'preparation_summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)
