#!/usr/bin/env bash
# 在 ecotyper 环境按官方入口恢复既有 Region3 CD+UC 模型；不修改 EcoTyper 源码。
# 先用每类型 100 细胞检查可运行性，再分析全部已抽样细胞并启用官方显著性检验。
# 输出分开保存；正式结果才用于科学图件，烟雾测试不能混入正式统计。
set -euo pipefail

# Step1. 使用现有环境并记录源码校验值
source /home/mazekai/miniconda3/etc/profile.d/conda.sh
conda activate ecotyper
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
root=/home/mazekai/IBD_EcoTyper/06_validation/GSE214695_scRNA_recovery_20261009
input="$root/input"
cd /home/mazekai/IBD_EcoTyper/software/ecotyper
sha256sum EcoTyper_recovery_scRNA.R pipeline/state_recovery_scRNA.R pipeline/ecotypes_recovery_scRNA.R pipeline/lib/ecotyper.R > "$root/source_before.sha256"
date -Is > "$root/run_started.txt"

# Step2. 用独立数据集名执行烟雾测试，防止与正式中间结果互相覆盖
cp "$input/GSE214695_Pilot18_20261009_CPM.txt" "$input/GSE214695_Pilot18_20261009_Smoke_CPM.txt"  # 官方 normalizePath 会解析符号链接，须实际复制以保留独立名称
Rscript EcoTyper_recovery_scRNA.R -d MSCCR_CORE12_CTSPEC_REGION_IBD -m "$input/GSE214695_Pilot18_20261009_Smoke_CPM.txt" -a "$input/GSE214695_Pilot18_20261009_annotation.tsv" -s 100 -z FALSE -t 4 -o "$root/smoke" > "$root/smoke.log" 2>&1
printf '0\n' > "$root/smoke_exit_code.txt"

# Step3. 正式小规模恢复；-s -1 表示保留全部已抽样细胞，-z TRUE 使用当前官方的 30 次打乱
Rscript EcoTyper_recovery_scRNA.R -d MSCCR_CORE12_CTSPEC_REGION_IBD -m "$input/GSE214695_Pilot18_20261009_CPM.txt" -a "$input/GSE214695_Pilot18_20261009_annotation.tsv" -s -1 -z TRUE -t 4 -o "$root/formal" > "$root/formal.log" 2>&1
sha256sum EcoTyper_recovery_scRNA.R pipeline/state_recovery_scRNA.R pipeline/ecotypes_recovery_scRNA.R pipeline/lib/ecotyper.R > "$root/source_after.sha256"
cmp "$root/source_before.sha256" "$root/source_after.sha256"
sha256sum --check "$root/frozen_model_before.sha256" > "$root/model_integrity.log"  # 验证冻结发现模型未被 recovery 改写
date -Is > "$root/run_finished.txt"
printf '0\n' > "$root/formal_exit_code.txt"
