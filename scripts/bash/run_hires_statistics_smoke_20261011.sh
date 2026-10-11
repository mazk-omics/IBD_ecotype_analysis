#!/usr/bin/env bash
# 在独立R环境重建HiRes数值阶段：完整基因入口、冻结24基因统计、收缩和校准。
# 只输出研究原型矩阵，最终NA规则未齐；不连接SSH，不运行官方容器或EcoTyper。
set -euo pipefail

# Step1. 显式指定目录和原始表达；不要把Adjusted导出自动替代原始输入
PROJECT_DIR="${1:?项目目录}"
INPUT_DIR="${2:?24基因固定基准目录}"
FULL_MIXTURE="${3:?原始完整基因mixture文件}"
MIXTURE_PANEL="${4:?同一mixture的冻结24基因子集}"
OUTPUT_DIR="${5:?不存在的新输出目录}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
mkdir "$OUTPUT_DIR"

# Step2. 全基因范围BH先完成，逐基因worker使用同一份全局统计
Rscript "$PROJECT_DIR/scripts/R/reconstruct_hires_entry_statistics_20261011.R" "$FULL_MIXTURE" "$INPUT_DIR/fractions_fixed.tsv" "$OUTPUT_DIR/entry"
Rscript "$PROJECT_DIR/scripts/R/benchmark_hires_numeric_core_20261011.R" "$INPUT_DIR" "$OUTPUT_DIR/core" 1 24 "$MIXTURE_PANEL"
Rscript "$PROJECT_DIR/scripts/R/reconstruct_hires_worker_statistics_20261011.R" "$INPUT_DIR" "$MIXTURE_PANEL" "$OUTPUT_DIR/entry/entry_statistics.rds" "$OUTPUT_DIR/worker"

# Step3. 复用同一份入口与worker统计，接续收缩、弱类型重估和OLS
Rscript "$PROJECT_DIR/scripts/R/reconstruct_hires_refinement_20261011.R" "$INPUT_DIR" "$MIXTURE_PANEL" "$OUTPUT_DIR/worker" "$OUTPUT_DIR/refinement"
