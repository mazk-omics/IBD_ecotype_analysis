#!/usr/bin/env bash
# 独立解释环境，沿用R分析接口；版本由用户允许使用兼容新版，不改ecotyper环境。
set -euo pipefail
OUT=/home/mazekai/IBD_EcoTyper/06_validation/ecotyper_original_method_annotation_20261010/environment
CONDA=/home/mazekai/miniconda3/bin/conda
"$CONDA" create -y -n ecotyper_interpretation --override-channels -c conda-forge -c bioconda 'r-base=4.3' 'r-seurat=4.4.0' 'bioconductor-fgsea=1.28.0' r-data.table r-igraph r-metap r-jsonlite
"$CONDA" env export -n ecotyper_interpretation > "$OUT/environment.yml"
"$CONDA" env export -n ecotyper_interpretation --from-history > "$OUT/environment_from_history.yml"
"$CONDA" list -n ecotyper_interpretation --explicit > "$OUT/conda_explicit_linux64.txt"
"$CONDA" list -n ecotyper_interpretation --json > "$OUT/conda_packages.json"
"$CONDA" run --no-capture-output -n ecotyper_interpretation Rscript -e 'library(Seurat); library(fgsea); library(data.table); library(igraph); library(metap); library(jsonlite); p <- installed.packages(); write.table(p[,c("Package","Version","Built","LibPath")], file=commandArgs(TRUE)[1], sep="\t", quote=FALSE, row.names=FALSE); print(sessionInfo()); stopifnot(is.function(FindMarkers), is.function(fgseaSimple)); print("PACKAGE_SMOKE_OK")' "$OUT/R_packages.tsv" > "$OUT/R_sessionInfo_and_smoke.txt" 2>&1
