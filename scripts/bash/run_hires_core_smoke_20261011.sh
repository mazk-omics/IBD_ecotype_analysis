#!/usr/bin/env bash
# 独立 R 数值核心的服务器烟雾入口；不运行官方容器，不连接 SSH，不调用 EcoTyper。
# 用已准备的24基因固定基准，读取全部120样本；输出目录必须是新目录，原件不覆盖。
set -euo pipefail

# Step1. 路径由调用者显式给出，便于沿用服务器的实际项目位置
PROJECT_DIR="${1:?请给出服务器项目目录}"
INPUT_DIR="${2:?请给出固定基准输入目录}"
OUTPUT_DIR="${3:?请给出新的烟雾测试输出目录}"
MIXTURE_PANEL="${4:?请显式给出24基因表达子集路径；本轮推荐原始mixture，不自动使用Adjusted导出}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

# Step2. 当前激活的独立 R 环境需含 nnls；记录由分析脚本保存
Rscript "$PROJECT_DIR/scripts/R/benchmark_hires_numeric_core_20261011.R" "$INPUT_DIR" "$OUTPUT_DIR" 1 24 "$MIXTURE_PANEL"
Rscript "$PROJECT_DIR/scripts/R/check_hires_column_randomization_20261011.R" "$INPUT_DIR" "$OUTPUT_DIR"
