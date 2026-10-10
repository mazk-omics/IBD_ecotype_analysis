# 整理120个真实MSCCR样本，用于下一次匹配参考的平台校正对照；本文件只准备测试输入，不提交任务。
# 在本机cibersortx_review环境运行；从现成2025样本CPM中保留原文本数值，输出mixture、Job16同样本基线和选择记录。
# 回肠/结肠/直肠各40例，120个不同供者；不依据预测比例、相关性或marker表达选择样本。
from pathlib import Path
import hashlib
import json
import shutil
import pandas as pd

# Step1. 按固定种子随机选样，部位相同的供者只取一份活检
root = Path('D:/IBD_ecotype/IBD_ecotype_analysis')
source = root / 'CIBERSORTx_ready_20261007/09_MSCCR_Joint20v2_Fractions'
out = root / 'CIBERSORTx_ready_20261008/01_MSCCR_platform120_pilot'
audit = out / '核对记录_勿上传'
audit.mkdir(parents=True, exist_ok=True)
annotation = pd.read_csv(source / '核对记录_勿上传/MSCCR_CD_UC_2025_annotation.tsv', sep='\t', dtype=str)
shuffled = annotation.sample(frac=1, random_state=20261008)
selected_blocks = []
used_patients = set()
for region in ['Ileum', 'Colon', 'Rectum']:
    pool = shuffled.loc[shuffled.Region3.eq(region) & ~shuffled.PatientID.isin(used_patients)]
    pool = pool.drop_duplicates('PatientID')
    assert len(pool) >= 40
    selected = pool.head(40).copy()
    selected_blocks.append(selected)
    used_patients.update(selected.PatientID)
selected_annotation = pd.concat(selected_blocks)
selected_annotation = annotation.loc[annotation.ID.isin(selected_annotation.ID)].copy()
assert len(selected_annotation) == 120 and selected_annotation.PatientID.nunique() == 120
assert selected_annotation.Region3.value_counts().eq(40).all()
selected_annotation.to_csv(audit / 'selected_annotation.tsv', sep='\t', index=False)
print(pd.crosstab(selected_annotation.Region3, selected_annotation.Disease).to_string(), flush=True)

# Step2. 从原文本只选列，不裁剪基因、不重算CPM、不改变已有数值文本
mixture_source = source / 'MSCCR_CD_UC_2025_GJ20_CPM.txt'
destination = out / 'MSCCR_Platform120_CPM.txt'
with mixture_source.open(encoding='utf-8') as handle:
    header = handle.readline().rstrip('\r\n').split('\t')
sample_columns = [i for i, name in enumerate(header) if i == 0 or name in set(selected_annotation.ID)]
assert [header[i] for i in sample_columns[1:]] == selected_annotation.ID.tolist()
line_count = 0
with mixture_source.open(encoding='utf-8') as handle, destination.open('w', encoding='utf-8', newline='') as output:
    for line in handle:
        fields = line.rstrip('\r\n').split('\t')
        assert len(fields) == 2026
        output.write('\t'.join(fields[i] for i in sample_columns) + '\n')
        line_count += 1
assert line_count == 17041
mixture = pd.read_csv(destination, sep='\t', index_col=0)
assert mixture.shape == (17040, 120) and mixture.columns.tolist() == selected_annotation.ID.tolist()
shutil.copyfile(source / 'GJ20v2_Signature20.txt', out / 'GJ20v2_Signature20.txt')

# Step3. 保留同样本Job16未校正基线；配套校正输入尚未确认，暂不让用户上传运行
baseline = pd.read_csv(root / 'CIBERSORTx_signature_audit/MSCCR_Joint20v2_Job16_review_20261008/original_outputs/CIBERSORTx_Job16_Results.txt', sep='\t', index_col=0)
baseline.loc[selected_annotation.ID].to_csv(audit / 'Job16_NS_same120_results.tsv', sep='\t', index_label='sample')
summary = {'preparation_complete': True, 'ready_to_submit_correction': False, 'seed': 20261008, 'samples': 120, 'different_patients': 120, 'genes': 17040, 'data_points': 2044800, 'regions': selected_annotation.Region3.value_counts().to_dict(), 'selection_basis': 'Region-balanced random patient selection; independent of predictions, fit or marker expression', 'original_numeric_strings_preserved': True, 'mixture_bytes': destination.stat().st_size, 'mixture_sha256': hashlib.sha256(destination.read_bytes()).hexdigest(), 'reference': 'Unchanged Joint20v2 Job13', 'correction_reference_prepared': False, 'web_job_submitted': False, 'server_accessed': False}
(audit / 'preparation_summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)
