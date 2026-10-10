# 整理 MSCCR 的 2,025 个 CD/UC discovery 样本与 Joint20 v2，供网页进行真实数据的细胞比例计算。
# 在本机 cibersortx_review 环境运行；读现成 CPM、注释和 Job13 signature，输出松散 TXT、提交设置及核对记录。
# 首轮沿用 Job15 的未校正设置；不混用旧版 S-mode reference，不改变样本、归一化或 signature。
from pathlib import Path
import hashlib
import json
import shutil
import numpy as np
import pandas as pd

# Step1. 定位现成文件，单独保存 v2 的正式数据提交目录
root = Path('D:/IBD_ecotype/IBD_ecotype_analysis')
ready = root / 'CIBERSORTx_ready_20261007'
out = ready / '09_MSCCR_Joint20v2_Fractions'
audit = out / '核对记录_勿上传'
audit.mkdir(parents=True, exist_ok=True)
(out / '网页结果').mkdir(exist_ok=True)
signature_source = ready / '08_Joint20v2_Natural28_validation/GJ20v2_Signature20.txt'
mixture_source = ready / '03_MSCCR_Joint20_exploratory/MSCCR_CD_UC_2025_GJ20_CPM.txt'
annotation_source = ready / '核对记录_勿上传/MSCCR_CD_UC_2025_annotation.tsv'
signature = pd.read_csv(signature_source, sep='\t', index_col=0)
annotation = pd.read_csv(annotation_source, sep='\t', dtype=str)
header = pd.read_csv(mixture_source, sep='\t', nrows=0).columns[1:]
assert signature.shape == (4474, 20) and signature.index.is_unique and signature.columns.is_unique
assert len(header) == 2025 and header.is_unique and annotation.ID.is_unique
assert header.tolist() == annotation.ID.tolist()
assert set(annotation.Disease) == {'CD', 'UC'} and not annotation[['ID', 'Disease', 'Region3']].isna().any().any()
print('Samples:', len(header), '; diseases:', annotation.Disease.value_counts().to_dict(), '; regions:', annotation.Region3.value_counts().to_dict(), flush=True)

# Step2. 分块回读全部数值，保持完整库 CPM 分母；循环是为避免整张矩阵占用过多内存
genes = []
nonzero_genes = set()
column_sums = np.zeros(len(header))
maximum = 0.0
for block in pd.read_csv(mixture_source, sep='\t', index_col=0, chunksize=500):
    values = block.to_numpy(dtype=float)
    assert np.isfinite(values).all() and (values >= 0).all()
    assert block.columns.tolist() == header.tolist()
    genes.extend(block.index.tolist())
    nonzero_genes.update(block.index[(values > 0).any(axis=1)].tolist())
    column_sums += values.sum(axis=0)
    maximum = max(maximum, float(values.max()))
genes = pd.Index(genes)
assert len(genes) == 17040 and genes.is_unique and not genes.isna().any()
assert len(genes) * len(header) == 34_506_000 and maximum > 50 and (column_sums > 0).all()
coverage = pd.DataFrame({'signature_gene': signature.index, 'in_MSCCR': signature.index.isin(genes), 'nonzero_in_MSCCR': signature.index.isin(nonzero_genes)})
coverage.to_csv(audit / 'v2_signature_MSCCR_gene_coverage.tsv', sep='\t', index=False)
coverage.loc[~coverage.in_MSCCR].to_csv(audit / 'v2_signature_MSCCR_missing_genes.tsv', sep='\t', index=False)
overlap = int(coverage.in_MSCCR.sum())
assert overlap / len(signature) > 0.9  # 本项目输入可用性核对，不代表反卷积准确率通过。
signal = pd.DataFrame({'total_signature_expression': signature.sum(axis=0), 'covered_signature_expression': signature.loc[signature.index.isin(genes)].sum(axis=0)})
signal['covered_expression_fraction'] = signal.covered_signature_expression / signal.total_signature_expression
signal.to_csv(audit / 'v2_reference_expression_coverage_by_type.tsv', sep='\t', index_label='cell_type')
pd.DataFrame({'sample': header, 'exported_background_CPM_sum': column_sums}).to_csv(audit / 'CPM_column_sums.tsv', sep='\t', index=False)
shutil.copyfile(annotation_source, audit / annotation_source.name)
print('Mixture genes × samples:', len(genes), '×', len(header), '; signature overlap:', overlap, '/', len(signature), '; max CPM:', maximum, flush=True)

# Step3. 校验来源哈希，原样复制；两文件不打包，不覆盖不同内容
source_summary = json.loads((root / 'CIBERSORTx_signature_audit/MSCCR_web_input_adaptation_20261006/server_outputs/audit/msccr_ibd_joint_web_20261006/summary.json').read_text(encoding='utf-8'))
with mixture_source.open('rb') as handle:
    mixture_hash = hashlib.file_digest(handle, 'sha256').hexdigest()
assert mixture_hash == source_summary['files'][mixture_source.name]['sha256']
checks = []
for source in [signature_source, mixture_source]:
    destination = out / source.name
    with source.open('rb') as handle:
        digest = hashlib.file_digest(handle, 'sha256').hexdigest()
    if destination.exists():
        with destination.open('rb') as handle:
            assert hashlib.file_digest(handle, 'sha256').hexdigest() == digest, '不覆盖已有不同输入'
    else:
        shutil.copyfile(source, destination)
    with destination.open('rb') as handle:
        assert hashlib.file_digest(handle, 'sha256').hexdigest() == digest
    checks.append({'source': str(source), 'file': source.name, 'bytes': destination.stat().st_size, 'sha256': digest, 'identical_to_source': True})
pd.DataFrame(checks).to_csv(audit / 'web_input_checksums.tsv', sep='\t', index=False)

# Step4. 为后续 HiRes/EcoTyper 留清楚入口，暂不构造细胞表达矩阵或启动远程任务
transition = ready / '10_EcoTyper_Joint20v2_transition'
transition.mkdir(exist_ok=True)
shutil.copyfile(annotation_source, transition / annotation_source.name)
shutil.copyfile(ready / '05_EcoTyper_transition/check_cibersortx_hires_for_ecotyper.py', transition / 'check_cibersortx_hires_for_ecotyper.py')
for mode in ['smoke', 'formal']:
    source = ready / f'05_EcoTyper_transition/config_cibersortx_joint20_region_ibd_{mode}.yml'
    text = source.read_text(encoding='utf-8-sig').replace('MSCCR_CIBERSORTX_JOINT20', 'MSCCR_CIBERSORTX_JOINT20V2')
    (transition / f'config_cibersortx_joint20v2_region_ibd_{mode}.yml').write_text(text, encoding='utf-8')
settings = {'module': 'Impute Cell Fractions', 'analysis_mode': 'Custom', 'signature': signature_source.name, 'mixture': mixture_source.name, 'disable_quantile_normalization': True, 'batch_correction': False, 'run_mode': 'relative', 'permutations': 100, 'suggested_job_label': 'MSCCR_IBD_J20v2_NS'}
summary = {'preparation_complete': True, 'web_job_submitted': False, 'server_accessed': False, 'sample_count': len(header), 'Disease_counts': annotation.Disease.value_counts().to_dict(), 'Region3_counts': annotation.Region3.value_counts().to_dict(), 'gene_count': len(genes), 'data_points': len(genes) * len(header), 'signature_genes': len(signature), 'signature_gene_overlap': overlap, 'signature_gene_overlap_fraction': overlap / len(signature), 'covered_signature_genes_nonzero': int(coverage.nonzero_in_MSCCR.sum()), 'maximum_CPM': maximum, 'sample_order_matches_frozen_discovery': True, 'source_hash_matches_audited_export': True, 'background_CPM_not_renormalized': True, 'recommended_web_settings': settings, 'note': 'Real discovery cohort run; retain validation limitations. Old single-cell reference not paired to this v2 signature.', 'HiRes_available': False, 'EcoTyper_started': False}
(audit / 'preparation_summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
(transition / 'transition_status.json').write_text(json.dumps({'MSCCR_v2_fractions_prepared': True, 'MSCCR_web_job_submitted': False, 'HiRes_available': False, 'EcoTyper_started': False, 'server_accessed': False, 'configuration_parameter_changes': 'Only dataset/path identifier JOINT20 -> JOINT20V2; numerical parameters and Region3 scaling unchanged'}, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)
