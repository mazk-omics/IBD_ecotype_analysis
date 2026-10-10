# 整理下一轮 CIBERSORTx 网页输入：回肠验证、SCP259 结直肠对照和 MSCCR 探索性基线。
# 本机项目 Python 环境运行；使用已冻结的参考和留出数据，保存普通目录与 TXT，不压缩、不提交网页。
# 只截取既有矩阵的样本列，原数值文本、细胞组成、真值与供者划分保持不变。
import csv
import hashlib
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

# Step1. 明确已有输入与输出目录，元数据单独存放，避免误上传真值
project = Path('D:/IBD_ecotype/IBD_ecotype_analysis')
audit = project / 'CIBERSORTx_signature_audit'
validation = audit / 'main_candidate_validation_20261006'
truth_dir = validation / 'server_outputs/validation_joint20_20261006'
candidates = audit / 'three_candidate_comparison_20261006'
pilots = project / 'IBDverse_reference/three_reference_pilots_20261006/run_20261006_183211'
msccr = audit / 'MSCCR_web_input_adaptation_20261006'
out = project / 'CIBERSORTx_ready_20261007'
metadata_dir = out / '核对记录_勿上传'
metadata_dir.mkdir(parents=True, exist_ok=True)

# Step2. 读入冻结的设计与真值；没有重新抽样，也没有接触最终检验供者表达
design = pd.read_csv(truth_dir / 'mixture_design.tsv', sep='\t').set_index('sample')
pairs = pd.read_csv(truth_dir / 'paired_equalRNA_samples.tsv', sep='\t').set_index('sample')
cell_truth = pd.read_csv(truth_dir / 'ground_truth_cell_fractions.tsv', sep='\t', index_col=0)
RNA_truth = pd.read_csv(truth_dir / 'ground_truth_export_UMI_fractions.tsv', sep='\t', index_col=0)
equal_truth = pd.read_csv(truth_dir / 'ground_truth_equalRNA_fractions.tsv', sep='\t', index_col=0)
membership = pd.read_csv(truth_dir / 'mixture_cell_membership.tsv', sep='\t')
natural = pd.read_csv(truth_dir / 'GJ20_Holdout28_Mixture_CPM.txt', sep='\t', index_col=0)
equal = pd.read_csv(truth_dir / 'GJ20_Holdout8_EqualRNA_Mixture_CPM.txt', sep='\t', index_col=0)
export_summary = json.loads((truth_dir / 'summary.json').read_text(encoding='utf-8'))
job6 = pd.read_csv(validation / 'Job6_review/original_outputs/CIBERSORTx_Job6_Results.txt', sep='\t', index_col=0)
job7 = pd.read_csv(validation / 'Job7_review/original_outputs/CIBERSORTx_Job7_Results.txt', sep='\t', index_col=0)
assert natural.shape == (23633, 28) and equal.shape == (23633, 8)
assert natural.index.equals(equal.index) and natural.index.is_unique
assert not membership.duplicated(['sample', 'cell_id']).any()
for source in [truth_dir / 'GJ20_Holdout28_Mixture_CPM.txt', truth_dir / 'GJ20_Holdout8_EqualRNA_Mixture_CPM.txt']:
    with source.open('rb') as handle:
        digest = hashlib.file_digest(handle, 'sha256').hexdigest()
    assert digest == export_summary['files'][source.name]['sha256']
print('Frozen mixtures verified:', natural.shape, equal.shape, flush=True)

# Step3. 复制两套已经审计的 signature，验证数值与原始网页输出完全相同
# 重复处理两个候选是同一套文件交接动作，不额外引入包装函数。
candidate_specs = [
    ('GI19', 'gse_ileum', 'GSE282122_Ileum19_signature_gene_symbols_restored.txt', '01_Ileum19_validation', 'Ileum', 19),
    ('SC18', 'scp_rectum', 'SCP259_ColonRectum18_signature_gene_symbols_restored.txt', '02_SCP259_ColonRectum_validation', 'ColonRectum', 18)
]
upload_records = []
job_records = []
coverage_records = []
for code, candidate, source_name, folder, region, n_types in candidate_specs:
    target_dir = out / folder
    target_dir.mkdir(exist_ok=True)
    signature_source = candidates / candidate / source_name
    candidate_audit = json.loads((candidates / candidate / 'audit_summary.json').read_text(encoding='utf-8'))
    original_signature = next(Path(item['path']) for item in candidate_audit['source_files'] if item['path'].endswith('.bm.K999.txt'))
    with original_signature.open('rb') as handle:
        original_hash = hashlib.file_digest(handle, 'sha256').hexdigest()
    original_record = next(item for item in candidate_audit['source_files'] if item['path'] == str(original_signature))
    assert original_hash == original_record['sha256']
    signature = pd.read_csv(signature_source, sep='\t', index_col=0)
    raw_signature = pd.read_csv(original_signature, sep='\t', index_col=0)
    assert signature.shape == raw_signature.shape and signature.shape[1] == n_types
    assert signature.columns.equals(raw_signature.columns) and np.array_equal(signature.to_numpy(), raw_signature.to_numpy())
    assert signature.index.is_unique and signature.columns.is_unique
    assert np.isfinite(signature.to_numpy()).all() and signature.ge(0).all().all() and signature.to_numpy().max() > 50
    signature_target = target_dir / f'{code}_Signature{n_types}.txt'
    shutil.copy2(signature_source, signature_target)
    with signature_source.open('rb') as handle:
        source_hash = hashlib.file_digest(handle, 'sha256').hexdigest()
    with signature_target.open('rb') as handle:
        copied_hash = hashlib.file_digest(handle, 'sha256').hexdigest()
    assert copied_hash == source_hash
    upload_records.append({'relative_path': str(signature_target.relative_to(out)), 'role': 'signature', 'genes': len(signature), 'columns': n_types, 'bytes': signature_target.stat().st_size, 'sha256': copied_hash, 'source': str(signature_source)})
    covered = signature.index.isin(natural.index)
    coverage_records.append({'candidate': code, 'signature_genes': len(signature), 'mixture_overlap_genes': int(covered.sum()), 'overlap_percent': float(covered.mean() * 100), 'unmatched_genes': int((~covered).sum())})
    pd.Series(signature.index[~covered], name='unmatched_signature_gene').to_csv(metadata_dir / f'{code}_unmatched_signature_genes.tsv', sep='\t', index=False)
    types = signature.columns

    # 按既有部位取列；缺少的类型必须在真值中为零，不能偷偷把缺失类型比例重新归一化。
    main_samples = design.index[design.region.eq(region)]
    equal_samples = pairs.index[pairs.paired_main_sample.isin(main_samples)]
    assert len(main_samples) == 12 and len(equal_samples) == 4
    assert set(types).issubset(cell_truth.columns)
    omitted_types = cell_truth.columns.difference(types)
    assert cell_truth.loc[main_samples, omitted_types].eq(0).all().all()
    assert RNA_truth.loc[main_samples, omitted_types].eq(0).all().all()
    assert equal_truth.loc[equal_samples, omitted_types].eq(0).all().all()

    # 直接复制原 TXT 对应列的数字字符串，不重算 CPM、不删基因、不修改 signature 标度。
    variants = [('Natural12', truth_dir / 'GJ20_Holdout28_Mixture_CPM.txt', main_samples, natural, cell_truth, RNA_truth),
                ('EqualRNA4', truth_dir / 'GJ20_Holdout8_EqualRNA_Mixture_CPM.txt', equal_samples, equal, equal_truth, equal_truth)]
    for variant, source, samples, expression, true_cell, true_RNA in variants:
        target = target_dir / f'{code}_Holdout{variant}_Mixture_CPM.txt'
        with source.open(encoding='utf-8', newline='') as input_handle, target.open('w', encoding='utf-8', newline='') as output_handle:
            rows = csv.reader(input_handle, delimiter='\t')
            header = next(rows)
            positions = [0] + [header.index(sample) for sample in samples]
            output_handle.write('\t'.join(header[position] for position in positions) + '\n')
            for row in rows:
                output_handle.write('\t'.join(row[position] for position in positions) + '\n')
        exported = pd.read_csv(target, sep='\t', index_col=0)
        assert exported.index.equals(expression.index) and exported.columns.equals(samples)
        assert np.array_equal(exported.to_numpy(), expression.loc[:, samples].to_numpy())
        assert np.isfinite(exported.to_numpy()).all() and exported.ge(0).all().all() and exported.to_numpy().max() > 50
        assert np.allclose(exported.sum(axis=0), 1e6, rtol=1e-6) and exported.size < 50000000
        with target.open('rb') as handle:
            digest = hashlib.file_digest(handle, 'sha256').hexdigest()
        upload_records.append({'relative_path': str(target.relative_to(out)), 'role': 'mixture', 'genes': len(exported), 'columns': len(samples), 'bytes': target.stat().st_size, 'sha256': digest, 'source': str(source)})
        job_records.append({'record_label': f'{code}_{variant}', 'priority': 1 if code == 'GI19' else 2, 'module': '2 Impute Cell Fractions', 'signature_file': str(signature_target.relative_to(out)), 'mixture_file': str(target.relative_to(out)), 'samples': len(samples), 'batch_correction': 'disabled', 'quantile_normalization': 'disabled', 'run_mode': 'relative', 'permutations': 100, 'reference_file': '', 'purpose': 'regional heldout reference comparison', 'submitted': False})
        selected_cell_truth = true_cell.loc[samples, types]
        selected_RNA_truth = true_RNA.loc[samples, types]
        assert np.allclose(selected_cell_truth.sum(axis=1), 1) and np.allclose(selected_RNA_truth.sum(axis=1), 1)
        selected_cell_truth.to_csv(metadata_dir / f'{code}_{variant}_ground_truth_cell.tsv', sep='\t')
        selected_RNA_truth.to_csv(metadata_dir / f'{code}_{variant}_ground_truth_RNA.tsv', sep='\t')

        # 保存旧 Joint20 输出和同轴误差，未来比较不能直接混用 20/19/18 类的 MAE 分母。
        old_predictions = job6.loc[samples, cell_truth.columns] if variant == 'Natural12' else job7.loc[samples, cell_truth.columns]
        old_predictions.to_csv(metadata_dir / f'{code}_{variant}_Joint20_baseline_predictions20.tsv', sep='\t')
        baseline = pd.DataFrame(index=samples)
        baseline['cell_MAE_common_types_pp'] = (old_predictions[types] - selected_cell_truth).abs().mean(axis=1) * 100
        baseline['RNA_MAE_common_types_pp'] = (old_predictions[types] - selected_RNA_truth).abs().mean(axis=1) * 100
        baseline['prediction_in_omitted_types'] = old_predictions[omitted_types].sum(axis=1)
        baseline.to_csv(metadata_dir / f'{code}_{variant}_Joint20_same_axis_baseline.tsv', sep='\t')
        print('Exported:', target.name, exported.shape, 'exact numerical equality', flush=True)

    # 供者和细胞泄漏检查针对实际参考细胞；SCP259 是外部来源，供者 ID 不作跨项目同名判断。
    heldout_cells = membership.loc[membership['sample'].isin(main_samples)]
    if code == 'GI19':
        training = pd.read_csv(pilots / 'gse_ileum/selected_cell_manifest.tsv', sep='\t')
        assert not set(training.cell_id) & set(heldout_cells.cell_id)
        assert not set(training.Patient) & set(heldout_cells.Patient)
    assert not set(heldout_cells.Patient) & {'CD12', 'UC17'}
    heldout_cells.to_csv(metadata_dir / f'{code}_natural12_existing_cell_membership.tsv', sep='\t', index=False)
    design.loc[main_samples].to_csv(metadata_dir / f'{code}_natural12_sample_design.tsv', sep='\t')
    pairs.loc[equal_samples].to_csv(metadata_dir / f'{code}_equalRNA4_pairing.tsv', sep='\t')

# Step4. 将 MSCCR 现成输入和 S-mode 所需训练 reference 放在独立目录，保留原文件字节
bulk_dir = out / '03_MSCCR_Joint20_exploratory'
bulk_dir.mkdir(exist_ok=True)
bulk_verification = json.loads((msccr / 'local_delivery_verification.json').read_text(encoding='utf-8'))
training_summary = json.loads((pilots / 'gse_joint/summary.json').read_text(encoding='utf-8'))
joint_input = validation / 'web_upload/GJ20_Signature20.txt'
joint_signature = pd.read_csv(joint_input, sep='\t', index_col=0)
joint_audited = pd.read_csv(candidates / 'gse_joint/GSE282122_Joint20_signature_gene_symbols_restored.txt', sep='\t', index_col=0)
assert joint_signature.equals(joint_audited)
bulk_sources = [(msccr / 'web_upload/MSCCR_CD_UC_2025_GJ20_CPM.txt', 'mixture', 17040, 2025, bulk_verification['txt_sha256']),
                (joint_input, 'signature', 5816, 20, None),
                (pilots / 'gse_joint/GSE282122_Joint20_scRNA_reference.txt', 'scRNA_reference_for_S_mode', 20612, 1000, training_summary['web_input_sha256'])]
for source, role, genes, columns, expected_hash in bulk_sources:
    with source.open('rb') as handle:
        digest = hashlib.file_digest(handle, 'sha256').hexdigest()
    if expected_hash is not None:
        assert digest == expected_hash
    target = bulk_dir / source.name
    shutil.copy2(source, target)
    with target.open('rb') as handle:
        target_digest = hashlib.file_digest(handle, 'sha256').hexdigest()
    assert target_digest == digest and genes * columns < 50000000
    upload_records.append({'relative_path': str(target.relative_to(out)), 'role': role, 'genes': genes, 'columns': columns, 'bytes': target.stat().st_size, 'sha256': digest, 'source': str(source)})
annotation_source = msccr / 'server_outputs/audit/msccr_ibd_joint_web_20261006/MSCCR_CD_UC_2025_annotation.tsv'
annotation = pd.read_csv(annotation_source, sep='\t')
with (bulk_dir / 'MSCCR_CD_UC_2025_GJ20_CPM.txt').open(encoding='utf-8') as handle:
    bulk_samples = handle.readline().rstrip('\r\n').split('\t')[1:]
assert len(bulk_samples) == 2025 and len(set(bulk_samples)) == 2025 and annotation.ID.tolist() == bulk_samples
assert annotation.Disease.value_counts().to_dict() == {'CD': 1154, 'UC': 871}
shutil.copy2(annotation_source, metadata_dir / 'MSCCR_CD_UC_2025_annotation.tsv')
bulk_genes = pd.read_csv(msccr / 'server_outputs/audit/msccr_ibd_joint_web_20261006/exported_gene_identity.tsv', sep='\t')['gene_symbol']
assert bulk_genes.is_unique and len(bulk_genes) == 17040
assert int(joint_signature.index.isin(bulk_genes).sum()) == 5719
job_records.append({'record_label': 'MSCCR_IBD_GJ20_Smode_exploratory', 'priority': 3, 'module': '2 Impute Cell Fractions',
    'signature_file': str((bulk_dir / 'GJ20_Signature20.txt').relative_to(out)), 'mixture_file': str((bulk_dir / 'MSCCR_CD_UC_2025_GJ20_CPM.txt').relative_to(out)),
    'samples': 2025, 'batch_correction': 'S-mode', 'quantile_normalization': 'disabled', 'run_mode': 'relative', 'permutations': 100,
    'reference_file': str((bulk_dir / 'GSE282122_Joint20_scRNA_reference.txt').relative_to(out)), 'purpose': 'exploratory real-bulk baseline; signature not validated', 'submitted': False})

# Step5. 保存上传清单、任务参数和来源 hash，提交时不需要解压任何文件
pd.DataFrame(upload_records).to_csv(metadata_dir / 'upload_file_integrity.tsv', sep='\t', index=False)
pd.DataFrame(job_records).to_csv(out / '提交任务清单.tsv', sep='\t', index=False)
pd.DataFrame(coverage_records).to_csv(metadata_dir / 'heldout_signature_gene_coverage.tsv', sep='\t', index=False)
shutil.copy2(truth_dir / 'donor_roles.tsv', metadata_dir / 'original_donor_roles.tsv')
summary = {'date': '2026-10-07', 'status': 'ready_local_files_only', 'compression_used': False, 'web_jobs_submitted': 0, 'server_access_performed': False,
    'new_sampling_performed': False, 'reserved_final_test_donors_used': False, 'upload_files': len(upload_records), 'prepared_jobs': len(job_records),
    'total_upload_bytes': sum(record['bytes'] for record in upload_records), 'heldout_gene_coverage': coverage_records,
    'MSCCR_samples': 2025, 'MSCCR_signature_overlap_genes': 5719, 'MSCCR_signature_total_genes': 5816,
    'main_reference_accuracy_accepted': False, 'fractions_are_not_cell_type_specific_expression_G': True,
    'directories': [item[3] for item in candidate_specs] + [bulk_dir.name]}
(out / '准备完成记录.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')

# Step6. 写清提交順序和固定参数，网页自动 Job 编号不等于这里用于追溯的任务代号
guide = '''# 2026-10-07 CIBERSORTx 待提交文件

所有上传文件都是普通 TXT，按三个任务目录分开放，不需要解压。文件已经准备、核对并保存；没有登录服务器，也没有代替你提交网页任务。

## 先提交回肠验证

目录：`01_Ileum19_validation`。上传三个文件后，先后建立两个模块 2 任务。

| 任务代号 | Signature | Mixture | 样本数 |
|---|---|---|---:|
| GI19_Natural12 | GI19_Signature19.txt | GI19_HoldoutNatural12_Mixture_CPM.txt | 12 |
| GI19_EqualRNA4 | GI19_Signature19.txt | GI19_HoldoutEqualRNA4_Mixture_CPM.txt | 4 |

两项统一设置：**Impute Cell Fractions → Custom；relative；勾选 Disable quantile normalization；batch correction disabled；100 permutations。** 本轮不选 S-mode，也不需要上传 scRNA reference 或 phenotype classes。沿用 Job6/Job7 的设定，仅改变候选 signature；实际参数请保存截图和 runtime log。

12 个自然样本为 balanced、epithelial_high、immune_high 各 4 个，每个 200 个细胞。4 个等 RNA 样本是 balanced 1/2、epithelial_high 1、immune_high 1 的原配对对照。它们来自昨天的同一批细胞，没有重新抽样；自然 counts 求和后 CPM 与逐细胞 CPM 平均是不同 RNA 贡献机制。

Ileum19 有 5,706 个基因、19 类，缺 Colonic_absorptive；本组回肠混合物该类真值全部为零。测试输入 23,633 个基因保留原数值，5,706 个特征基因全部覆盖。没有预先只保留 signature 基因，也没有按变化后的参考重算 CPM 或真值。

## 第二组：SCP259 结直肠参考对照

目录：`02_SCP259_ColonRectum_validation`。文件已备好，可接着逐个提交两个模块 2 任务。

| 任务代号 | Signature | Mixture | 样本数 |
|---|---|---|---:|
| SC18_Natural12 | SC18_Signature18.txt | SC18_HoldoutNatural12_Mixture_CPM.txt | 12 |
| SC18_EqualRNA4 | SC18_Signature18.txt | SC18_HoldoutEqualRNA4_Mixture_CPM.txt | 4 |

参数与回肠两项相同：**relative、QN disabled、batch correction disabled、100 permutations。** 参考是 SCP259；测试细胞仍来自 GSE282122 的固定结直肠留出混合物，属于跨数据集参考测试，不能称为 SCP259 自身留出验证。

SCP259 有 4,612 个基因、18 类，缺 Paneth 和 Ileal_absorptive；本组结直肠混合物这两类真值均为零，无需修改组成或归一化真值。它不能直接作为全 MSCCR 的完整 20 类参考。

测试输入精确匹配 4,381/4,612 个 signature 基因（94.9913%），231 个基因名未匹配；缺失名单保存在核对记录。此次没有猜测旧名/别名，也没有填零或修改 signature，保留包括 4 个线粒体基因在内的原矩阵。解读需同时考虑跨数据集、类型定义与基因标识差异。

## 第三组：MSCCR 全 CD+UC 探索性基线

目录：`03_MSCCR_Joint20_exploratory`。此组独立保存，优先完成上述小规模验证后再考虑提交；当前 Joint20 仍有回肠偏差，真实 bulk 的结果只能先作为探索性基线。

| 网页字段 | 对应文件 / 设置 |
|---|---|
| 模块 | 2. Impute Cell Fractions / Custom |
| Signature | GJ20_Signature20.txt |
| Mixture | MSCCR_CD_UC_2025_GJ20_CPM.txt |
| Batch correction | S-mode |
| S-mode 单细胞 reference | GSE282122_Joint20_scRNA_reference.txt |
| Run mode | relative |
| Disable quantile normalization | 勾选 |
| Permutations | 100 |
| 记录用任务代号 | MSCCR_IBD_GJ20_Smode_exploratory |

MSCCR 输入是 17,040 基因 × 2,025 样本，即 34,506,000 个数据点；CD 1,154、UC 871，与原 discovery 列顺序一致。未重新做 region 标准化，单位仍为非 log CPM。覆盖 Joint20 的 5,719/5,816 个基因（98.33%）。TXT 约 336 MB，不需要解压；网页上传和计算耗时无法预先保证。

S-mode 所需 reference 是原来训练用的 20,612 基因 × 1,000 个单细胞，原始文件直接复制，约 41.5 MB；不能用 inferred_refsample、signature 或留出 mixture 替代。网页如果已上传同一文件，可直接选已有条目。S-mode 与留出测试的无校正设定不同，需单独解释。此目录不包含回肠/SCP 的训练 reference，因为它们的当前测试不使用 S-mode。

## 上传与结果保存

1. 每个目录只放网页要上传的 TXT；`核对记录_勿上传` 内是真值、来源和核对记录，不需要传到 CIBERSORTx。
2. 每次只运行一个任务；上一轮 GSE 任务并发时出现过结果文件丢失，串行提交已成功。按任务清单逐项记录网页实际 Job 编号，不预先假定会是 Job8 等。
3. 保存 Results.txt、完整 runtime_log.txt、Job Parameters 截图，以及可下载的两份 PDF。提交清单中的任务代号用于记录；如网页没有自定义任务名字段，不必寻找或填写。
4. 输入名称和 type 顺序已核对；真值在本地直接由细胞名单计数，未用 signature 或预测 fraction 生成。共同类型的比较指标已另存，不能直接把 20 类、19 类、18 类 MAE 混用。
5. CD12、UC17 的最终留出表达仍未使用。这一轮是开发验证，不是最终独立性能评估，也不进行临床关联分析。
6. 此轮提交得到的是 fractions；EcoTyper discovery 所需的细胞类型特异表达 G 是后续计算。本包没有把 fraction 当作 G，也没有静态生成高分辨率表达输入。HiRes 接续任务需结合比例结果再准备。

建议回来后先提交 GI19_Natural12，再依次跑 GI19_EqualRNA4、SC18_Natural12、SC18_EqualRNA4。MSCCR 已备齐，作为后续探索项单独保留。无需重复 Job6/Job7 的联合参考测试，也无需重建当前 signature。
'''
(out / '先读_提交顺序和参数.md').write_text(guide, encoding='utf-8')
print(json.dumps(summary, ensure_ascii=True, indent=2), flush=True)
