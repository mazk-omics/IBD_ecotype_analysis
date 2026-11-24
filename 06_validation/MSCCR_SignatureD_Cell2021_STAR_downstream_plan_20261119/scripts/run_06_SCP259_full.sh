#!/usr/bin/env bash
# 本脚本顺序完成 SCP259 输入构建、官方 recovery、审计和 100 次 shuffle。
set -eu -o pipefail

RUN=/home/mazekai/IBD_EcoTyper/06_validation/MSCCR_SignatureD_Cell2021_STAR_downstream_plan_20261119
R_PREP=/home/mazekai/miniconda3/envs/ibd_refbuild/bin/Rscript
R_ECOTYPER=/home/mazekai/miniconda3/envs/ecotyper/bin/Rscript
STATUS="$RUN/output/06_SCP259_full_status.txt"
LOG_DIR="$RUN/output/06_SCP259_full_logs"
mkdir -p "$LOG_DIR"

CURRENT_STAGE=INIT
trap 'code=$?; if [ "$code" -ne 0 ]; then printf "STATUS=FAILED\nSTAGE=%s\nEXIT=%s\nTIME=%s\n" "$CURRENT_STAGE" "$code" "$(date -Is)" >> "$STATUS"; fi' EXIT
printf "STATUS=RUNNING\nDATASET=SCP259\nWORKERS=4\nSTART=%s\n" "$(date -Is)" > "$STATUS"

# Step1. 构建并审计冻结输入
CURRENT_STAGE=PREPARE_INPUTS
printf "STAGE=%s\nTIME=%s\n" "$CURRENT_STAGE" "$(date -Is)" >> "$STATUS"
"$R_PREP" "$RUN/scripts/06a_prepare_SCP259_recovery_inputs.R" "$RUN" > "$LOG_DIR/06a_prepare.console.txt" 2>&1

CURRENT_STAGE=AUDIT_INPUTS
printf "STAGE=%s\nTIME=%s\n" "$CURRENT_STAGE" "$(date -Is)" >> "$STATUS"
"$R_PREP" "$RUN/scripts/06b_audit_SCP259_recovery_inputs.R" "$RUN" > "$LOG_DIR/06b_audit.console.txt" 2>&1

# Step2. 运行官方 recovery 并审计结果
CURRENT_STAGE=OFFICIAL_RECOVERY
printf "STAGE=%s\nTIME=%s\n" "$CURRENT_STAGE" "$(date -Is)" >> "$STATUS"
"$R_ECOTYPER" "$RUN/scripts/06c_run_SCP259_official_recovery.R" "$RUN" > "$LOG_DIR/06c_recovery.console.txt" 2>&1

CURRENT_STAGE=AUDIT_RECOVERY
printf "STAGE=%s\nTIME=%s\n" "$CURRENT_STAGE" "$(date -Is)" >> "$STATUS"
"$R_PREP" "$RUN/scripts/06d_audit_SCP259_official_recovery.R" "$RUN" > "$LOG_DIR/06d_audit.console.txt" 2>&1

# Step3. 完成论文标准的 100 次 shuffle
CURRENT_STAGE=SIGNIFICANCE_100
printf "STAGE=%s\nTIME=%s\n" "$CURRENT_STAGE" "$(date -Is)" >> "$STATUS"
"$R_ECOTYPER" "$RUN/scripts/06e_run_SCP259_recovery_significance_100.R" "$RUN" > "$LOG_DIR/06e_significance.console.txt" 2>&1

CURRENT_STAGE=ACCEPTED
printf "STATUS=ACCEPTED\nSTAGE=%s\nEND=%s\n" "$CURRENT_STAGE" "$(date -Is)" >> "$STATUS"
