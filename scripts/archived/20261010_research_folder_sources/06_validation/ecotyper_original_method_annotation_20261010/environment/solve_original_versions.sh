#!/usr/bin/env bash
# 先解析原文版本的依赖，不修改已有环境；日志保留安装可行性证据。
set -e
/home/mazekai/miniconda3/bin/conda create -n ecotyper_interpretation -c conda-forge -c bioconda 'r-base=4.0.3' 'r-seurat=3.1.3' 'bioconductor-fgsea=1.14.0' r-data.table r-igraph r-metap --dry-run --json
