#!/usr/bin/env bash
# 服务器恢复入口：优先复用本轮完整Recovery；全部统计从保存评分及临床映射复算。
# 只写本轮目录。运行前用environment.yml建立ibd_outcomes_20261015，另需既有ecotyper环境。
set -euo pipefail
root=/home/mazekai/IBD_EcoTyper
out="$root/06_validation/external_clinical_outcomes_20261015"
python=/home/mazekai/miniconda3/envs/ibd_outcomes_20261015/bin/python
export PATH=/home/mazekai/miniconda3/envs/ecotyper/bin:$PATH
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$root"

# Step1. 核对原始结局，准备原值输入及隔离的官方运行副本。
"$python" "$out/scripts/01_audit_inputs.py" > "$out/logs/01_audit.log" 2>&1
Rscript "$out/scripts/02_prepare_GSE73661.R" > "$out/logs/02_prepare_GSE73661.log" 2>&1
"$python" "$out/scripts/00_prepare_runtime.py" > "$out/logs/00_prepare_runtime.log" 2>&1
export R_PROFILE_USER="$out/scripts/seed_recovery.Rprofile"
cd "$out/recovery_runtime"
if [ ! -s "$out/03_GSE73661_mucosal_healing/RecoveryOutput/GSE73661_RMA_W0_67_outcomes_20261015/Ecotypes/ecotype_abundance.txt" ]; then
  Rscript EcoTyper_recovery_bulk.R -d MSCCR_CORE12_CTSPEC_REGION_IBD -m "$out/03_GSE73661_mucosal_healing/GSE73661_RMA_W0_67_outcomes_20261015.txt" -a "$out/03_GSE73661_mucosal_healing/recovery_annotation.txt" -t 4 -o "$out/03_GSE73661_mucosal_healing/RecoveryOutput" > "$out/logs/03_GSE73661_recovery.log" 2>&1
fi
if [ ! -s "$out/02_RISK_disease_progression/RecoveryOutput/GSE93624_processed245_outcomes_20261015/Ecotypes/ecotype_abundance.txt" ]; then
  Rscript EcoTyper_recovery_bulk.R -d MSCCR_CORE12_CTSPEC_REGION_IBD -m "$out/02_RISK_disease_progression/GSE93624_processed245_outcomes_20261015.txt" -a "$out/02_RISK_disease_progression/recovery_annotation.txt" -t 4 -o "$out/02_RISK_disease_progression/RecoveryOutput" > "$out/logs/04_RISK_recovery.log" 2>&1
fi

# Step2. 固定检验、数值核验、图件与正式报告依次保存。
"$python" "$out/scripts/03_fit_outcomes.py" > "$out/logs/statistics_all.log" 2>&1
"$python" "$out/scripts/04_recovery_and_source_QC.py" > "$out/logs/07_recovery_source_QC.log" 2>&1
"$python" "$out/scripts/05_independent_Firth_validation.py" > "$out/logs/08_independent_Firth_validation.log" 2>&1
"$python" "$out/scripts/06_plot_outcomes.py" > "$out/logs/09_plots.log" 2>&1
"$python" "$out/scripts/07_write_reports.py" > "$out/logs/10_reports.log" 2>&1
"$python" "$out/scripts/08_delivery_QC.py" > "$out/logs/11_delivery_QC.log" 2>&1
date -Is > "$out/logs/COMPLETED"
