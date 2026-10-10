#!/usr/bin/env bash
# 服务器 ecotyper 环境：检查输入，运行24例测试，通过后恢复完整73例。
# 输出写入本轮独立目录；仅用已有官方 recovery，不重做反卷积或 discovery。
set -euo pipefail
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
run_dir="$1"
conda_bin=/home/mazekai/miniconda3/bin/conda
software_dir=/home/mazekai/IBD_EcoTyper/software/ecotyper
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
stage=preflight
trap 'printf "FAILED stage=%s exit=%s\n" "$stage" "$?" > "$run_dir/FAILED"' ERR

# Step1. 冻结输入和模型，保留中间审计。
"$conda_bin" run --no-capture-output -n ecotyper Rscript "$script_dir/prepare_gse72819_recovery.R" "$run_dir" > "$run_dir/preflight.log" 2>&1
run_id="$(basename -- "$run_dir")"
cd "$software_dir"

# Step2. 2个并行类型运行烟雾测试；通过后才继续，遇到失败不修改分配阈值。
stage=smoke
"$conda_bin" run --no-capture-output -n ecotyper Rscript EcoTyper_recovery_bulk.R -d MSCCR_CORE12_CTSPEC_REGION_IBD -m "$run_dir/GSE72819_BPsmoke24_${run_id}.txt" -a "$run_dir/smoke_annotation.txt" -t 2 -o "$run_dir/smoke_results" > "$run_dir/smoke.log" 2>&1
"$conda_bin" run --no-capture-output -n ecotyper Rscript "$script_dir/validate_gse72819_recovery.R" "$run_dir" smoke > "$run_dir/smoke_validation.log" 2>&1

# Step3. 使用原73例完整数据。4个并行类型是本次资源设置，不改变 NMF 或生态型算法。
stage=full
"$conda_bin" run --no-capture-output -n ecotyper Rscript EcoTyper_recovery_bulk.R -d MSCCR_CORE12_CTSPEC_REGION_IBD -m "$run_dir/GSE72819_BP73_${run_id}.txt" -a "$run_dir/full_annotation.txt" -t 4 -o "$run_dir/full_results" > "$run_dir/full.log" 2>&1
"$conda_bin" run --no-capture-output -n ecotyper Rscript "$script_dir/validate_gse72819_recovery.R" "$run_dir" full > "$run_dir/full_validation.log" 2>&1
date -Is > "$run_dir/COMPLETED"
