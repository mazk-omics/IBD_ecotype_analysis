#!/usr/bin/env bash
# 外部smoke通过后顺序执行正式生成、BayesPrism和准确性评价，失败即停止。
# 服务器同一SSH会话启动；计算日志和大counts对象都写入独立数据集目录。
set -euo pipefail
parent="$1"
dataset="$2"
folder="$parent/08_external_pseudobulk/$dataset"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1

# Step1. 真正原始counts混合，100份×1000细胞
test -f "$folder/01_smoke/smoke_complete.txt"
date -u > "$folder/logs/formal_start_UTC.txt"
/usr/bin/time -v /home/mazekai/miniconda3/envs/ibd_refaudit/bin/python "$parent/scripts/generate_external.py" "$parent" "$dataset" formal > "$folder/logs/generate_formal.log" 2>&1
printf '0\n' > "$folder/logs/generate_formal.exit_code"

# Step2. 完整20类BayesPrism，不改变参考或默认Gibbs迭代
/usr/bin/time -v /home/mazekai/miniconda3/envs/bayesprism/bin/Rscript "$parent/scripts/run_external.R" "$parent" "$dataset" formal > "$folder/logs/bayesprism_formal.log" 2>&1
printf '0\n' > "$folder/logs/bayesprism_formal.exit_code"

# Step3. 使用实际保留基因真值，三尺度分别评价
/usr/bin/time -v /home/mazekai/miniconda3/envs/ibd_refaudit/bin/python "$parent/scripts/evaluate_external.py" "$parent" "$dataset" > "$folder/logs/accuracy.log" 2>&1
printf '0\n' > "$folder/logs/accuracy.exit_code"
date -u > "$folder/logs/formal_finish_UTC.txt"
