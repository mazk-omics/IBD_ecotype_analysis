# 整理同一批120例MSCCR的网页重跑文件；不重新抽样、归一化或修改矩阵。
# 在本机cibersortx_review环境运行，读取已冻结的Job17输入，输出独立提交目录和核对记录。
from pathlib import Path
import hashlib
import json
import shutil
import numpy as np
import pandas as pd

# Step1. 读取已冻结的样本和Job17结果
project = Path('D:/IBD_ecotype/IBD_ecotype_analysis')
source = project / 'CIBERSORTx_ready_20261008/01_MSCCR_platform120_pilot'
outdir = project / 'CIBERSORTx_ready_20261008/02_MSCCR120_resubmit'
records = outdir / '核对记录_勿上传'
hires_dir = project / 'CIBERSORTx_ready_20261008/03_MSCCR120_HiRes_next'
outdir.mkdir(exist_ok=True)
records.mkdir(exist_ok=True)
hires_dir.mkdir(exist_ok=True)
annotation = pd.read_csv(source / '核对记录_勿上传/selected_annotation.tsv', sep='\t')
job17 = pd.read_csv(project / 'CIBERSORTx_signature_audit/MSCCR_Joint20v2_Job17_review_20261008/original_outputs/CIBERSORTx_Job17_Adjusted.txt', sep='\t', index_col=0)
assert len(annotation) == 120 and annotation['PatientID'].nunique() == 120
assert job17.index.tolist() == annotation['ID'].tolist(), 'Job17与冻结样本顺序不同'
print(annotation.groupby(['Region3', 'Disease']).size().to_string())
annotation.to_csv(records / 'MSCCR120_annotation.tsv', sep='\t', index=False)
job17.to_csv(records / 'Job17_既有比例结果_不能当作G.tsv', sep='\t')

# Step2. 原样复制三个输入，并检查矩阵方向、数值和来源哈希
# 同版多文件用一个循环；重复细胞类型表头是单细胞参考的有意设计。
file_names = ['GJ20v2_Signature20.txt', 'MSCCR_Platform120_CPM.txt', 'GJ20v2_SourceSubset1060_scRNA_CPM.txt']
expected_shapes = [(4474, 20), (17040, 120), (25703, 1060)]
checks = []
for file_name, expected_shape in zip(file_names, expected_shapes):
    input_path = source / file_name
    output_path = outdir / file_name
    source_hash = hashlib.sha256(input_path.read_bytes()).hexdigest()
    if output_path.exists():
        assert hashlib.sha256(output_path.read_bytes()).hexdigest() == source_hash, f'已有同名文件内容不同：{output_path}'
    else:
        shutil.copy2(input_path, output_path)
    output_hash = hashlib.sha256(output_path.read_bytes()).hexdigest()
    assert source_hash == output_hash
    with output_path.open(encoding='utf-8-sig') as handle:
        header = handle.readline().rstrip('\r\n').split('\t')
        columns = header[1:]
        gene_names = []
        maximum = 0.0
        for line in handle:
            fields = line.rstrip('\r\n').split('\t')
            assert len(fields) == len(header) and all(field != '' for field in fields)
            assert '"' not in line
            values = np.array(fields[1:], dtype=float)
            assert np.isfinite(values).all() and (values >= 0).all()
            gene_names.append(fields[0])
            maximum = max(maximum, float(values.max()))
    shape = (len(gene_names), len(columns))
    assert shape == expected_shape and len(set(gene_names)) == len(gene_names)
    assert maximum >= 50, '非log输入最大值过小，网页可能误判为log'
    if file_name == 'MSCCR_Platform120_CPM.txt':
        assert columns == annotation['ID'].tolist()
        mixture_genes = gene_names
    if file_name == 'GJ20v2_Signature20.txt':
        signature_types = columns
    if file_name == 'GJ20v2_SourceSubset1060_scRNA_CPM.txt':
        assert set(columns) == set(signature_types)
    checks.append({'file': file_name, 'genes': shape[0], 'columns': shape[1], 'data_points': shape[0] * shape[1], 'bytes': output_path.stat().st_size, 'max_expression': maximum, 'sha256': output_hash, 'identical_to_Job17_frozen_input': True})
    print(file_name, shape, output_path.stat().st_size)
pd.DataFrame(checks).to_csv(records / '上传文件核对.tsv', sep='\t', index=False)

# Step3. 保存HiRes全基因请求清单，不按比例结果或预期状态筛基因
# 这是17040基因的请求清单，不是1000基因试运行；单独保存，避免上传到错误文件类型。
gene_request = hires_dir / 'MSCCR120_All17040_GeneSubset.txt'
gene_request.write_text('\n'.join(mixture_genes) + '\n', encoding='utf-8')
summary = {'status': 'local_files_prepared_not_uploaded_not_submitted', 'samples': 120, 'distinct_patients': 120, 'region_counts': annotation['Region3'].value_counts().to_dict(), 'files': checks, 'Job17_fractions_already_available': True, 'HiRes_expression_available': False, 'HiRes_full_gene_request_count': len(mixture_genes), 'HiRes_30_minute_completion_guaranteed': False, 'S_mode_pairing': 'same source-matched 1060-cell diagnostic subset as Job17; exact original pseudobulk construction input equivalence not established'}
(records / '准备状态.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
print('提交文件：', outdir)
print('HiRes后续：', hires_dir)
