# 为120例MSCCR整理独立HiRes上传目录，保留冻结输入与样本顺序，并记录Job17、Job20比例。
# 本机cibersortx_review环境运行；复制冻结文件，保存全基因请求与既有比例核对记录，不启动计算。
from pathlib import Path
import hashlib
import json
import shutil
import numpy as np
import pandas as pd

# Step1. 读取冻结样本、已检查的输入清单和Job17比例
project = Path('D:/IBD_ecotype/IBD_ecotype_analysis')
source = project / 'CIBERSORTx_ready_20261008/02_MSCCR120_resubmit'
outdir = project / 'CIBERSORTx_ready_20261008/03_MSCCR120_HiRes_next'
records = outdir / '核对记录_勿上传'
records.mkdir(parents=True, exist_ok=True)
annotation = pd.read_csv(source / '核对记录_勿上传/MSCCR120_annotation.tsv', sep='\t')
previous_checks = pd.read_csv(source / '核对记录_勿上传/上传文件核对.tsv', sep='\t')
fractions = pd.read_csv(source / '核对记录_勿上传/Job17_既有比例结果_不能当作G.tsv', sep='\t', index_col=0)
assert len(annotation) == 120 and annotation['PatientID'].nunique() == 120
assert fractions.index.tolist() == annotation['ID'].tolist()
print('样本、供者：', len(annotation), annotation['PatientID'].nunique())
print(annotation.groupby(['Region3', 'Disease']).size().to_string())

# Step2. 原样复制输入到HiRes目录，核对与上次完整数值检查的哈希
# 三个矩阵同版复制是必要循环；不再次抽样、筛基因或对矩阵重新归一化。
checks = []
for row in previous_checks.itertuples(index=False):
    input_path = source / row.file
    output_path = outdir / row.file
    source_hash = hashlib.sha256(input_path.read_bytes()).hexdigest()
    assert source_hash == row.sha256, f'冻结输入已变化：{input_path}'
    if output_path.exists():
        assert hashlib.sha256(output_path.read_bytes()).hexdigest() == source_hash, f'同名文件内容不同：{output_path}'
    else:
        shutil.copy2(input_path, output_path)
    output_hash = hashlib.sha256(output_path.read_bytes()).hexdigest()
    assert output_hash == source_hash
    checks.append({'file': row.file, 'genes': int(row.genes), 'columns': int(row.columns), 'data_points': int(row.data_points), 'bytes': output_path.stat().st_size, 'sha256': output_hash, 'unchanged_from_Job17_input': True})
    print(row.file, int(row.genes), int(row.columns), output_path.stat().st_size)

# Step3. 从完整mixture逐行提取全基因请求，不改变表达背景
with (outdir / 'MSCCR_Platform120_CPM.txt').open(encoding='utf-8-sig') as handle:
    mixture_header = handle.readline().rstrip('\r\n').split('\t')
    genes = [line.split('\t', 1)[0] for line in handle]
assert mixture_header[1:] == annotation['ID'].tolist()
assert len(genes) == 17040 and len(set(genes)) == 17040
assert all(gene.strip() == gene and gene != '' for gene in genes)
gene_path = outdir / 'MSCCR120_All17040_GeneSubset.txt'
gene_text = '\n'.join(genes) + '\n'
if gene_path.exists():
    assert gene_path.read_text(encoding='utf-8') == gene_text, '已有全基因清单与mixture不一致'
else:
    gene_path.write_text(gene_text, encoding='utf-8')
checks.append({'file': gene_path.name, 'genes': len(genes), 'columns': 1, 'data_points': None, 'bytes': gene_path.stat().st_size, 'sha256': hashlib.sha256(gene_path.read_bytes()).hexdigest(), 'unchanged_from_Job17_input': None})

# Step4. 留存纯比例与注释，避免将质量指标或比例表错当成G上传
signature_types = pd.read_csv(outdir / 'GJ20v2_Signature20.txt', sep='\t', nrows=0).columns[1:].tolist()
assert len(signature_types) == 20 and set(signature_types).issubset(fractions.columns)
fraction20 = fractions[signature_types].copy()
assert np.isfinite(fraction20.to_numpy()).all() and (fraction20.to_numpy() >= 0).all()
assert np.allclose(fraction20.sum(axis=1).to_numpy(), 1, atol=1e-8)
fraction20.to_csv(records / 'Job17_Fractions20_仅核对网页是否支持导入.txt', sep='\t', index_label='Mixture')
annotation.to_csv(records / 'MSCCR120_annotation.tsv', sep='\t', index=False)
pd.DataFrame({'CellType': signature_types}).to_csv(records / '预期20类.tsv', sep='\t', index=False)
pd.DataFrame(checks).to_csv(records / '上传文件维度与哈希.tsv', sep='\t', index=False)
summary = {'status': 'local_hires_files_prepared_not_uploaded_not_submitted', 'samples': 120, 'distinct_patients': 120, 'types': 20, 'requested_genes': len(genes), 'sample_to_type_ratio': 6, 'upload_files': checks, 'Job17_fractions_reused_as_audit_record_only': True, 'actual_web_fraction_import_field_verified': False, 'actual_web_HiRes_S_mode_fields_verified': False, 'request_scope': 'all 17040 input genes; no arbitrary 1000-marker reduction', 'estimated_runtime': None, 'Job19_error_reported_limit_hours': 0.5, 'Job19_user_observed_elapsed_hours_approximate': 4, 'whole_job_wall_clock_limit_verified': False, 'SourceSubset1060_pairing': 'source-matched diagnostic subset, not the literal 212-column signature-construction input'}
(records / 'HiRes准备状态.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
print('HiRes独立目录：', outdir)
print('本地已准备，尚未上传或提交。')

# Step5. 固定现有120例，补充Job20核对记录；不按炎症状态重新抽样
# 用户于2026-10-09决定保留当前队列组成。两次比例只作审阅记录，不当作HiRes表达输入。
job20_path = project / 'CIBERSORTx_signature_audit/MSCCR_Joint20v2_Job20_review_20261009/Job20_results_with_frozen_annotation.tsv'
job20 = pd.read_csv(job20_path, sep='\t', index_col='sample')
assert job20.index.tolist() == annotation['ID'].tolist()
assert job20['PatientID'].tolist() == annotation['PatientID'].tolist()
for column in ['Disease', 'Inflammation', 'Region3']:
    assert job20[column].tolist() == annotation[column].tolist(), f'冻结注释不一致：{column}'
job20_fraction20 = job20[signature_types].copy()
assert np.isfinite(job20_fraction20.to_numpy()).all() and (job20_fraction20.to_numpy() >= 0).all()
assert np.allclose(job20_fraction20.sum(axis=1).to_numpy(), 1, atol=1e-8)
job20_fraction20.to_csv(records / 'Job20_Fractions20_仅核对网页是否支持导入.txt', sep='\t', index_label='Mixture')
composition = pd.crosstab(annotation['Region3'], annotation['Inflammation']).reindex(['Ileum', 'Colon', 'Rectum'])
assert composition.to_numpy().sum() == 120 and set(annotation['Disease']) == {'CD', 'UC'}
composition.to_csv(records / '固定120例_部位与炎症组成.tsv', sep='\t')
summary['sample_set_fixed_by_user_on'] = '2026-10-09'
summary['disease_counts'] = annotation['Disease'].value_counts().to_dict()
summary['inflammation_counts'] = annotation['Inflammation'].value_counts().to_dict()
summary['control_samples'] = 0
summary['Job20_fractions_reused_as_audit_record_only'] = True
summary['Job20_fraction_file'] = 'Job20_Fractions20_仅核对网页是否支持导入.txt'
summary['web_gene_request_limit_confirmed_20261009'] = 1250
summary['all17040_request_directly_runnable_on_web'] = False
summary['current_submission_guide'] = '04_MSCCR120_HiRes_gene_batches/现在先读_两批HiRes提交.md'
(records / 'HiRes准备状态.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
print('固定样本组成：')
print(composition.to_string())
