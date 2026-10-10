# 将联合参考唯一的 S-mode 对照与 EcoTyper 接口文件分开放好，不压缩、不提交网页、不登录服务器。
# 在本地 cibersortx_review 环境运行；复制冻结输入，生成哈希清单、2025样本注释和两级计算配置。
# 配置只有在得到真实 HiRes 表达、完成输入与 Region3 源码检查后才使用。
from pathlib import Path
import hashlib
import json
import shutil
import numpy as np
import pandas as pd

# Step1. 沿用已审计的开发集、参考和 MSCCR 注释，不引入新的供者或基因筛选
root = Path(__file__).resolve().parents[2]
ready = root / 'CIBERSORTx_ready_20261007'
diagnostic = ready / '04_Joint20_Smode_diagnostic'
transition = ready / '05_EcoTyper_transition'
diagnostic.mkdir(parents=True, exist_ok=True)
transition.mkdir(parents=True, exist_ok=True)
validation = root / 'CIBERSORTx_signature_audit/main_candidate_validation_20261006/server_outputs/validation_joint20_20261006'
tasks = [(ready / '03_MSCCR_Joint20_exploratory/GJ20_Signature20.txt', diagnostic / 'GJ20_Signature20.txt'),
    (ready / '03_MSCCR_Joint20_exploratory/GSE282122_Joint20_scRNA_reference.txt', diagnostic / 'GSE282122_Joint20_scRNA_reference.txt'),
    (validation / 'GJ20_Holdout28_Mixture_CPM.txt', diagnostic / 'GJ20_Holdout28_Mixture_CPM.txt'),
    (ready / '核对记录_勿上传/MSCCR_CD_UC_2025_annotation.tsv', transition / 'MSCCR_CD_UC_2025_annotation.tsv'),
    (root / 'scripts/python/check_cibersortx_hires_for_ecotyper.py', transition / 'check_cibersortx_hires_for_ecotyper.py')]
integrity = []
for source, destination in tasks:
    with source.open('rb') as handle:
        source_sha = hashlib.file_digest(handle, 'sha256').hexdigest()
    if destination.exists():
        with destination.open('rb') as handle:
            existing_sha = hashlib.file_digest(handle, 'sha256').hexdigest()
        assert existing_sha == source_sha, f'已有目标与源不同，拒绝覆盖：{destination}'
    else:
        shutil.copy2(source, destination)
    with destination.open('rb') as handle:
        destination_sha = hashlib.file_digest(handle, 'sha256').hexdigest()
    assert source_sha == destination_sha
    integrity.append({'source': str(source), 'destination': str(destination), 'sha256': source_sha, 'identical_copy': True})

# Step2. 核对三个网页文件轴和尺寸；新任务只改变 batch correction 为 S-mode
signature = pd.read_csv(diagnostic / 'GJ20_Signature20.txt', sep='\t', index_col=0)
mixture = pd.read_csv(diagnostic / 'GJ20_Holdout28_Mixture_CPM.txt', sep='\t', index_col=0)
reference = pd.read_csv(diagnostic / 'GSE282122_Joint20_scRNA_reference.txt', sep='\t', index_col=0)
assert signature.index.is_unique and mixture.index.is_unique and reference.index.is_unique
assert signature.index.isin(mixture.index).all() and signature.index.isin(reference.index).all()
assert np.isfinite(mixture.to_numpy()).all() and (mixture.to_numpy() >= 0).all()
assert reference.shape == (20612, 1000) and reference.to_numpy().max() >= 50
assert mixture.shape[1] == 28 and mixture.size <= 50_000_000
truth = pd.read_csv(validation / 'ground_truth_export_UMI_fractions.tsv', sep='\t', index_col=0)
assert mixture.columns.tolist() == truth.index.tolist()
assert set(signature.columns) == set(truth.columns)
pd.DataFrame(integrity).to_csv(transition / 'prepared_file_integrity.tsv', sep='\t', index=False)
pd.DataFrame([{'file': 'GJ20_Signature20.txt', 'role': 'Signature Matrix', 'rows': len(signature), 'columns': signature.shape[1]},
    {'file': 'GJ20_Holdout28_Mixture_CPM.txt', 'role': 'Mixture', 'rows': len(mixture), 'columns': mixture.shape[1]},
    {'file': 'GSE282122_Joint20_scRNA_reference.txt', 'role': 'Single cell reference for S-mode', 'rows': len(reference), 'columns': reference.shape[1]}]).to_csv(transition / 'diagnostic_upload_manifest.tsv', sep='\t', index=False)

# Step3. 先短试运行再完整计算；模型名称隔离，避免覆盖现有 BayesPrism 基线
remote_root = '/home/mazekai/IBD_EcoTyper'
input_root = remote_root + '/05_EcoTyper/input/MSCCR_CIBERSORTX_JOINT20'
for label, restarts, ranks in [('SMOKE', 2, 4), ('FORMAL', 50, 20)]:
    # 2/4 是接口试运行约定；50/20 沿用旧 discovery 方法以保持可比性，不承诺耗时减少。
    config = f'''# 待真实 HiRes 输入检查和服务器 Region3 源码核对通过后使用；不是已启动任务。
default :
  Input :
    Discovery dataset name : "MSCCR_CIBERSORTX_JOINT20_REGION_IBD_{label}"
    Expression matrices : "{input_root}/expression_matrices"
    Annotation file : "{input_root}/annotation/MSCCR_CD_UC_2025_annotation.tsv"
    Annotation file column to scale by : "Region3"
    Annotation file column(s) to plot : ["Disease", "Inflammation", "Region3"]
  Output :
    Output folder : "{remote_root}/05_EcoTyper/output/MSCCR_CIBERSORTX_JOINT20/region_ibd_{label.lower()}"
  Pipeline settings :
    Pipeline steps to skip : [6]
    Filter genes : "cell type specific"
    Number of threads : 18
    Number of NMF restarts : {restarts}
    Maximum number of states per cell type : {ranks}
    Cophenetic coefficient cutoff : 0.95
    Minimum number of states in ecotypes : 3
'''
    config_path = transition / f'config_cibersortx_joint20_region_ibd_{label.lower()}.yml'
    assert not config_path.exists(), '配置已存在，避免覆盖用户之后修改的参数'
    config_path.write_text(config, encoding='utf-8')
annotation = pd.read_csv(transition / 'MSCCR_CD_UC_2025_annotation.tsv', sep='\t')
assert len(annotation) == 2025 and annotation.ID.is_unique
assert annotation.Disease.value_counts().to_dict() == {'CD': 1154, 'UC': 871}
assert annotation.Region3.value_counts().to_dict() == {'Colon': 799, 'Rectum': 676, 'Ileum': 550}
(transition / 'transition_status.json').write_text(json.dumps({'status': 'files_prepared_waiting_for_web_deconvolution_and_HiRes',
    'diagnostic_web_job_submitted': False, 'MSCCR_web_job_submitted': False, 'MSCCR_status_source': 'user_report_not_submitted',
    'HiRes_expression_available': False, 'EcoTyper_started': False, 'server_accessed': False,
    'reference_not_declared_validated': True, 'final_holdout_used': False}, ensure_ascii=False, indent=2), encoding='utf-8')
print('Prepared diagnostic:', diagnostic)
print('Prepared EcoTyper interface:', transition)
