#!/usr/bin/env bash
# 固定fraction之后的研究流程：独立Group、HiRes入口、worker、收缩及Group辅助回填。
# R+nnls环境运行；尚非完整官方替代，服务器入口目前仅完成Bash语法检查。
set -euo pipefail

# Step1. 保持输入来源可见，不读取官方Group或HiRes输出作为计算输入
PROJECT_DIR="${1:?项目目录}"
INPUT_DIR="${2:?目标输入目录，含fraction和gene_selection}"
FULL_MIXTURE="${3:?原始完整mixture}"
OUTPUT_DIR="${4:?不存在的新输出目录}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
mkdir "$OUTPUT_DIR"

# Step2. Group过滤阈值与HiRes入口分别在完整基因范围计算
Rscript "$PROJECT_DIR/scripts/R/reconstruct_group_expression_20261011.R" "$FULL_MIXTURE" "$INPUT_DIR/fractions_fixed.tsv" "$OUTPUT_DIR/group"
Rscript "$PROJECT_DIR/scripts/R/reconstruct_hires_entry_statistics_20261011.R" "$FULL_MIXTURE" "$INPUT_DIR/fractions_fixed.tsv" "$OUTPUT_DIR/entry"

# Step3. 目标清单的行号来自原始完整mixture，不能使用Adjusted表的遗留行号
Rscript "$PROJECT_DIR/scripts/R/reconstruct_hires_worker_statistics_20261011.R" "$INPUT_DIR" "$INPUT_DIR/mixture_original_target_genes.tsv" "$OUTPUT_DIR/entry/entry_statistics.rds" "$OUTPUT_DIR/worker"
Rscript "$PROJECT_DIR/scripts/R/reconstruct_hires_refinement_20261011.R" "$INPUT_DIR" "$INPUT_DIR/mixture_original_target_genes.tsv" "$OUTPUT_DIR/worker" "$OUTPUT_DIR/refinement"
Rscript "$PROJECT_DIR/scripts/R/reconstruct_hires_auxiliary_refill_20261011.R" "$OUTPUT_DIR/refinement" "$OUTPUT_DIR/group/GEPs_Filtered.tsv" "$OUTPUT_DIR/auxiliary"

# 输出仍依赖上游固定fraction；版本差异需审计，不直接输入正式EcoTyper研究
printf '%s\n' 'Research calculation completed. Review Group/HiRes differences before downstream use.'
