#!/usr/bin/env bash
# 在 ecotyper 环境直接调用既有官方 recovery；只新增结果目录，不改源码或冻结模型。
# 本轮健康/nonIBD 输入由独立脚本准备；六个队列顺序运行，避免抢占服务器资源。
set -euo pipefail
root="$1"
source /home/mazekai/miniconda3/etc/profile.d/conda.sh
conda activate ecotyper
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /home/mazekai/IBD_EcoTyper/software/ecotyper
date -Is > "$root/run_started.txt"

# Step1. MSCCR 按 Region3 分开使用官方行标准化，外部回肠 bulk 独立恢复
for name in MSCCR_nonIBD_Ileum_20261010 MSCCR_nonIBD_Colon_20261010 MSCCR_nonIBD_Rectum_20261010 GSE57945_nonIBD_Ileum_20261010; do
    Rscript EcoTyper_recovery_bulk.R -d MSCCR_CORE12_CTSPEC_REGION_IBD -m "$root/input/$name.txt" -a "$root/input/${name}_annotation.tsv" -t 4 -o "$root/formal" > "$root/${name}.log" 2>&1
    printf '0\n' > "$root/${name}_exit_code.txt"
done

# Step2. 沿用已质控的健康细胞；-s -1 保留子集全部细胞，-z TRUE 为官方当前置换设置
for name in IBDverse_Healthy_Ileum_20261010 GSE214695_HC_Colon_20261010; do
    Rscript EcoTyper_recovery_scRNA.R -d MSCCR_CORE12_CTSPEC_REGION_IBD -m "$root/input/$name.txt" -a "$root/input/${name}_annotation.tsv" -s -1 -z TRUE -t 4 -o "$root/formal" > "$root/${name}.log" 2>&1
    printf '0\n' > "$root/${name}_exit_code.txt"
done
date -Is > "$root/run_finished.txt"
