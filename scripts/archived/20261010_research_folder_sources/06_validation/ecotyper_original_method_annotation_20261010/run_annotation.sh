#!/usr/bin/env bash
# 等待独立环境安装完成后做数值烟雾测试，再计算文献对应和生态型共关联。
# 使用同一个新环境；失败即记录退出，不改官方EcoTyper或原始输入。
set -euo pipefail
OUT=/home/mazekai/IBD_EcoTyper/06_validation/ecotyper_original_method_annotation_20261010
CONDA=/home/mazekai/miniconda3/bin/conda
while tmux has-session -t ibd_annotation_env 2>/dev/null; do sleep 10; done
"$CONDA" run --no-capture-output -n ecotyper_interpretation Rscript /home/mazekai/IBD_EcoTyper/scripts/R/interpretation_20261010/smoke_paper_annotation_packages_20261010.R > "$OUT/environment/numerical_smoke.log" 2>&1
"$CONDA" run --no-capture-output -n ecotyper_interpretation Rscript /home/mazekai/IBD_EcoTyper/scripts/R/interpretation_20261010/interpret_states_ecotypes_paper_20261010.R > "$OUT/analysis.log" 2>&1
