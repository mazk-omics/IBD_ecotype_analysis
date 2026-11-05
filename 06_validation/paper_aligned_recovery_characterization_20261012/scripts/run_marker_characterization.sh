#!/usr/bin/env bash
# 全计数导出通过后，运行已验收的Seurat3.1.3 DEG和两版n4留出比较。
# 在同一SSH会话启动的tmux中运行；失败停在该步骤，保留已完成中间结果。
set -euo pipefail
CHAR_OUT=/home/mazekai/IBD_EcoTyper/06_validation/paper_aligned_recovery_characterization_20261012

# Step1. 核对上游检查点；只等待已启动的提取任务，不重启SSH或重新抽样
while [ ! -f "$CHAR_OUT/logs/preflight_full_counts.exit" ]; do sleep 10; done
test "$(cat "$CHAR_OUT/logs/preflight_full_counts.exit")" = 0
test -f "$CHAR_OUT/02_external_markers/full_counts_ready.txt"

# Step2. 正式差异表达；参数和日志由R脚本保存
/home/mazekai/miniconda3/envs/ecotyper_interpretation/bin/Rscript "$CHAR_OUT/scripts/external_markers_seurat313.R" > "$CHAR_OUT/logs/external_DE.log" 2>&1
printf '0' > "$CHAR_OUT/logs/external_DE.exit"

# Step3. 两版n4及逐数据集留出；未完成DE不进入排名
/home/mazekai/miniconda3/envs/ibd_refaudit/bin/python -u "$CHAR_OUT/scripts/rank_markers_n4_and_LOOCV.py" > "$CHAR_OUT/logs/marker_ranking_LOOCV.log" 2>&1
printf '0' > "$CHAR_OUT/logs/marker_ranking_LOOCV.exit"
printf 'PASS\n' > "$CHAR_OUT/logs/marker_workflow_complete.txt"
