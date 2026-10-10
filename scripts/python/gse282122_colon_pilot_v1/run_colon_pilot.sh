#!/usr/bin/env bash
set -euo pipefail
GSE_PILOT_ROOT=/home/mazekai/IBD_EcoTyper
GSE_PILOT_SCRIPTS="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
GSE_PILOT_OUT="${GSE_PILOT_ROOT}/03_reference/GSE282122/web_pilot/output/colon_pilot_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$(dirname "${GSE_PILOT_OUT}")"
/home/mazekai/miniconda3/bin/conda run --no-capture-output -n ibd_refaudit \
 python -u "${GSE_PILOT_SCRIPTS}/export_gse282122_colon_web_pilot.py" \
 --input /home/mazekai/enteric_glia/data_scRNA/GSE282122/TAURUS_raw_counts_annotated_final.h5ad \
 --candidate "${GSE_PILOT_ROOT}/03_reference/GSE282122/audit_v2/output/04_identity_coverage_v1/mapped_candidate_cells_v1.parquet" \
 --outdir "${GSE_PILOT_OUT}" --cells-per-type 100 --seed 20261006 \
 2>&1 | tee "${GSE_PILOT_OUT}.log"
echo "Colon candidate output: ${GSE_PILOT_OUT}"
