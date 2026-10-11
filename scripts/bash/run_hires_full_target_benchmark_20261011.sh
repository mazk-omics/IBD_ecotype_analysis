#!/usr/bin/env bash
# 独立HiRes研究流程：完整基因入口、固定目标统计、收缩、可选辅助回填。
# 在R+nnls环境运行；比例及可选Group表仍是外部输入，不连接SSH、不调用官方可执行文件。
set -euo pipefail

# Step1. 输入目录须包含fractions_fixed.tsv、gene_selection.tsv及原始目标表达
PROJECT_DIR="${1:?项目目录}"
INPUT_DIR="${2:?完整目标基因基准目录}"
FULL_MIXTURE="${3:?原始完整基因mixture文件}"
OUTPUT_DIR="${4:?不存在的新输出目录}"
AUXILIARY_GROUP="${5:-}"  # 不提供时跳过回填，避免隐式借用官方Group输出
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
mkdir "$OUTPUT_DIR"

# Step2. BH范围为完整输入；worker与收缩依次使用同一份固定统计
Rscript "$PROJECT_DIR/scripts/R/reconstruct_hires_entry_statistics_20261011.R" "$FULL_MIXTURE" "$INPUT_DIR/fractions_fixed.tsv" "$OUTPUT_DIR/entry"
Rscript "$PROJECT_DIR/scripts/R/reconstruct_hires_worker_statistics_20261011.R" "$INPUT_DIR" "$INPUT_DIR/mixture_original_target_genes.tsv" "$OUTPUT_DIR/entry/entry_statistics.rds" "$OUTPUT_DIR/worker"
Rscript "$PROJECT_DIR/scripts/R/reconstruct_hires_refinement_20261011.R" "$INPUT_DIR" "$INPUT_DIR/mixture_original_target_genes.tsv" "$OUTPUT_DIR/worker" "$OUTPUT_DIR/refinement"

# Step3. 整类型全1保护在完整目标集合上判断；提供的Group可靠性表是显式外部依赖
if [[ -n "$AUXILIARY_GROUP" ]]; then
    Rscript "$PROJECT_DIR/scripts/R/reconstruct_hires_auxiliary_refill_20261011.R" "$OUTPUT_DIR/refinement" "$AUXILIARY_GROUP" "$OUTPUT_DIR/auxiliary"
fi
