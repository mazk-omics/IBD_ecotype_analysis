#!/usr/bin/env bash
# 分层DE和第一版排名成功后，补齐53状态的缺失留出折，并生成最终证据表。
# 由当前持久SSH会话启动的tmux运行；前一步失败时停止，不自动重连或改模型。
set -euo pipefail
CHAR_OUT=/home/mazekai/IBD_EcoTyper/06_validation/paper_aligned_recovery_characterization_20261012

# Step1. 确认已运行计算成功，保留最初报告和未补齐的逐折表以供追溯
while [ ! -f "$CHAR_OUT/logs/marker_parallel4_runner.exit" ]; do sleep 10; done
test "$(cat "$CHAR_OUT/logs/marker_parallel4_runner.exit")" = 0
cp "$CHAR_OUT/02_external_markers/marker_LOOCV_per_dataset.tsv" "$CHAR_OUT/02_external_markers/marker_LOOCV_per_dataset_initial.tsv"
cp "$CHAR_OUT/09_reports/METHODS_PARITY_AUDIT.md" "$CHAR_OUT/09_reports/METHODS_PARITY_AUDIT_initial.md"

# Step2. 完整逐状态覆盖表；算法和已计算的marker不变，缺失折记录不可评价
/home/mazekai/miniconda3/envs/ibd_refaudit/bin/python -u "$CHAR_OUT/scripts/rank_markers_n4_and_LOOCV_parallel4_v3.py" > "$CHAR_OUT/logs/marker_ranking_LOOCV_final.log" 2>&1
/home/mazekai/miniconda3/envs/ibd_refaudit/bin/python -u "$CHAR_OUT/scripts/assemble_characterization_reports.py" > "$CHAR_OUT/logs/assemble_reports_final.log" 2>&1
/home/mazekai/miniconda3/envs/ibd_refaudit/bin/python -u "$CHAR_OUT/scripts/check_final_marker_outputs.py" > "$CHAR_OUT/logs/final_marker_QC.log" 2>&1
printf 'PASS\n' > "$CHAR_OUT/logs/final_marker_reports_complete.txt"
