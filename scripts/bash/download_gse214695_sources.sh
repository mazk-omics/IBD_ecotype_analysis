#!/usr/bin/env bash
# 在现有服务器会话中下载公开原始文件，最多 3 路 HTTP；不是 3 个 SSH 登录。
# 下载清单来自已核对的 GEO 元数据；完整文件复用，临时文件支持断点续传。
set -euo pipefail
root=/home/mazekai/IBD_EcoTyper/06_validation/GSE214695_scRNA_recovery_20261009
mkdir -p "$root/source"

# Step1. 文件名和 URL 不含空格，按两列传给现有 curl；三路仅用于独立网络传输
tail -n +2 "$root/download_urls.tsv" | cut -f2,3 | xargs -n 2 -P 3 bash -c '
set -euo pipefail
root=/home/mazekai/IBD_EcoTyper/06_validation/GSE214695_scRNA_recovery_20261009
file="$root/source/$1"
if [ -f "$file" ]; then exit 0; fi
printf "Downloading: %s\n" "$1"
curl --fail --location --retry 2 --retry-delay 5 --connect-timeout 30 --max-time 900 --continue-at - --silent --show-error --output "$file.part" "$2"
gzip -t "$file.part"
mv "$file.part" "$file"
printf "Completed: %s\n" "$1"
' download

# Step2. 下载完成后，仍由原研究脚本逐样本核对和处理
bash "$root/run_all.sh"
