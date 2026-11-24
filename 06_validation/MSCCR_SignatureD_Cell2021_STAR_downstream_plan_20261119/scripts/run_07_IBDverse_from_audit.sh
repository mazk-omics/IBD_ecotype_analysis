#!/usr/bin/env bash
# 本脚本复用已完成的 IBDverse 输入，从输入审计继续官方 recovery 和 100 次 shuffle。
set -eu -o pipefail

RUN=/home/mazekai/IBD_EcoTyper/06_validation/MSCCR_SignatureD_Cell2021_STAR_downstream_plan_20261119
R_PREP=/home/mazekai/miniconda3/envs/ibd_refbuild/bin/Rscript
R_ECOTYPER=/home/mazekai/miniconda3/envs/ecotyper/bin/Rscript
STATUS="$RUN/output/07_IBDverse_full_status.txt"
LOG_DIR="$RUN/output/07_IBDverse_full_logs"
mkdir -p "$LOG_DIR"

CURRENT_STAGE=AUDIT_INPUTS
trap 'code=$?; if [ "$code" -ne 0 ]; then printf "STATUS=FAILED\nSTAGE=%s\nEXIT=%s\nTIME=%s\n" "$CURRENT_STAGE" "$code" "$(date -Is)" >> "$STATUS"; fi' EXIT
printf "STATUS=RUNNING\nDATASET=IBDverse_E_MTAB_16999\nWORKERS=4\nRESUME_FROM=AUDIT_INPUTS\nSTART=%s\n" "$(date -Is)" > "$STATUS"

# Step1. 审计现有冻结输入
printf "STAGE=%s\nTIME=%s\n" "$CURRENT_STAGE" "$(date -Is)" >> "$STATUS"
"$R_PREP" "$RUN/scripts/07b_audit_IBDverse_recovery_inputs.R" "$RUN" > "$LOG_DIR/07b_audit.console.txt" 2>&1

# Step2. 运行官方 recovery 并审计结果
CURRENT_STAGE=OFFICIAL_RECOVERY
printf "STAGE=%s\nTIME=%s\n" "$CURRENT_STAGE" "$(date -Is)" >> "$STATUS"
"$R_ECOTYPER" "$RUN/scripts/07c_run_IBDverse_official_recovery.R" "$RUN" > "$LOG_DIR/07c_recovery.console.txt" 2>&1

CURRENT_STAGE=AUDIT_RECOVERY
printf "STAGE=%s\nTIME=%s\n" "$CURRENT_STAGE" "$(date -Is)" >> "$STATUS"
"$R_PREP" "$RUN/scripts/07d_audit_IBDverse_official_recovery.R" "$RUN" > "$LOG_DIR/07d_audit.console.txt" 2>&1

# Step3. 完成论文标准的 100 次 shuffle
CURRENT_STAGE=SIGNIFICANCE_100
printf "STAGE=%s\nTIME=%s\n" "$CURRENT_STAGE" "$(date -Is)" >> "$STATUS"
"$R_ECOTYPER" "$RUN/scripts/07e_run_IBDverse_recovery_significance_100.R" "$RUN" > "$LOG_DIR/07e_significance.console.txt" 2>&1

CURRENT_STAGE=ACCEPTED
printf "STATUS=ACCEPTED\nSTAGE=%s\nEND=%s\n" "$CURRENT_STAGE" "$(date -Is)" >> "$STATUS"
