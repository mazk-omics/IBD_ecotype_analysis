# 在服务器为官方Recovery准备本轮独立工作副本，冻结Discovery仅通过只读链接读取。
# 小型官方代码复制到结果目录；输入数值不变，所有新Recovery写在本轮目录。
from pathlib import Path
import shutil, hashlib
import pandas as pd

# Step1. 复制代码并连接冻结模型，不复制或重训大型参考。
root = Path('/home/mazekai/IBD_EcoTyper')
out = root / '06_validation/external_clinical_outcomes_20261015'
src = root / 'software/ecotyper'
runtime = out / 'recovery_runtime'
runtime.mkdir(exist_ok=True)
shutil.copytree(src / 'pipeline', runtime / 'pipeline', dirs_exist_ok=True)
shutil.copy2(src / 'EcoTyper_recovery_bulk.R', runtime / 'EcoTyper_recovery_bulk.R')
discovery = 'MSCCR_CORE12_CTSPEC_REGION_IBD'
model = runtime / 'EcoTyper' / discovery / 'Cell_type_specific_genes'
model.mkdir(parents=True, exist_ok=True)
links = [(runtime / 'utils', src / 'utils'), (model / 'Analysis', src / 'EcoTyper' / discovery / 'Cell_type_specific_genes/Analysis'), (runtime / 'EcoTyper' / discovery / 'config_used.yml', src / 'EcoTyper' / discovery / 'config_used.yml')]
for name in ['Cell_States', 'Ecotypes']:
    (model / name / 'recovery').mkdir(parents=True, exist_ok=True)
    links.append((model / name / 'discovery', src / 'EcoTyper' / discovery / 'Cell_type_specific_genes' / name / 'discovery'))
for link, target in links:
    if not link.exists():
        link.symlink_to(target, target_is_directory=target.is_dir())
    assert link.resolve() == target.resolve()
(runtime / 'datasets/bulk').mkdir(parents=True, exist_ok=True)

# Step2. 唯一文件名必须用真实副本；官方normalizePath会解开符号链接的文件名。
source = root / '06_validation/01_bulk_progression/GSE93624/prepared/06_ecotyper_recovery_input/GSE93624_EcoTyper_recovery_expression.txt'
target = out / '02_RISK_disease_progression/GSE93624_processed245_outcomes_20261015.txt'
if not target.exists():
    shutil.copyfile(source, target)
assert hashlib.sha256(source.read_bytes()).digest() == hashlib.sha256(target.read_bytes()).digest()
header = pd.read_csv(target, sep='\t', nrows=0).columns[1:]
pd.DataFrame({'ID': header, 'Cohort': 'GSE93624'}).to_csv(out / '02_RISK_disease_progression/recovery_annotation.txt', sep='\t', index=False)
(out / 'scripts/seed_recovery.Rprofile').write_text('set.seed(20261015)\n')
print('ISOLATED_RUNTIME_READY')
