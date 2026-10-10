#!/usr/bin/env bash
# 原文Seurat 3.1.3的conda二进制仅适配R3.5/3.6；采用R4.0依赖并源码安装同版Seurat。
# fgsea仍固定原文1.14.0；所有安装操作限于新环境，日志和源码校验随结果保存。
set -euo pipefail
OUT=/home/mazekai/IBD_EcoTyper/06_validation/ecotyper_original_method_annotation_20261010/environment
CONDA=/home/mazekai/miniconda3/bin/conda
CONDA_PLUGINS_USE_SHARDED_REPODATA=false "$CONDA" create --solver libmamba -y -n ecotyper_interpretation --override-channels -c conda-forge -c bioconda 'r-base=4.0.3' 'r-seurat=3.2.3' 'bioconductor-fgsea=1.14.0' r-sdmtools r-rcppprogress r-ape r-tsne r-data.table r-igraph r-metap r-jsonlite compilers
curl --fail --location --max-time 120 https://cran.r-project.org/src/contrib/Archive/Seurat/Seurat_3.1.3.tar.gz -o "$OUT/Seurat_3.1.3.tar.gz"
sha256sum "$OUT/Seurat_3.1.3.tar.gz" > "$OUT/Seurat_source.sha256"
"$CONDA" run --no-capture-output -n ecotyper_interpretation R CMD INSTALL "$OUT/Seurat_3.1.3.tar.gz" > "$OUT/Seurat_source_install.log" 2>&1
"$CONDA" run --no-capture-output -n ecotyper_interpretation Rscript "$OUT/save_environment_20261010.R" > "$OUT/R_sessionInfo_and_smoke.txt" 2>&1
"$CONDA" env export -n ecotyper_interpretation > "$OUT/environment.yml"
"$CONDA" env export -n ecotyper_interpretation --from-history > "$OUT/environment_from_history.yml"
"$CONDA" list -n ecotyper_interpretation --explicit > "$OUT/conda_explicit_linux64.txt"
"$CONDA" list -n ecotyper_interpretation --json > "$OUT/conda_packages.json"
