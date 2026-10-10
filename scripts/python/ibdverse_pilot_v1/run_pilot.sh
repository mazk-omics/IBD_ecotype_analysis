#!/usr/bin/env bash
set -euo pipefail

IBD_PILOT_PROJECT="${IBD_PILOT_PROJECT:-/home/mazekai/IBD_EcoTyper}"
IBD_PILOT_CONDA="${IBD_PILOT_CONDA:-/home/mazekai/miniconda3/bin/conda}"
IBD_PILOT_SCRIPTS="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
IBD_PILOT_INPUT="/home/mazekai/enteric_glia/data_scRNA/scIBD_external/IBDverse_E-MTAB-16998/E-MTAB-16999/Full_cohort.h5ad"
IBD_PILOT_COVERAGE="${1:-${IBD_PILOT_PROJECT}/logs/ibdverse_reference_inputs_20261006_161527/label_category_sample_donor_coverage.tsv}"
IBD_PILOT_OUTPUT_ROOT="${IBD_PILOT_PROJECT}/03_reference/ibdverse_cibersortx_pilot/output"
IBD_PILOT_OUTPUT="${IBD_PILOT_OUTPUT_ROOT}/pilot_$(date +%Y%m%d_%H%M%S)"

test -f "${IBD_PILOT_INPUT}" || { echo "Missing source: ${IBD_PILOT_INPUT}" >&2; exit 1; }
test -f "${IBD_PILOT_COVERAGE}" || { echo "Missing audit coverage: ${IBD_PILOT_COVERAGE}" >&2; exit 1; }
mkdir -p "${IBD_PILOT_OUTPUT_ROOT}"

"${IBD_PILOT_CONDA}" run --no-capture-output -n ibd_refaudit \
  python "${IBD_PILOT_SCRIPTS}/export_ibdverse_cibersortx_pilot.py" --version

"${IBD_PILOT_CONDA}" run --no-capture-output -n ibd_refaudit \
  python -u "${IBD_PILOT_SCRIPTS}/export_ibdverse_cibersortx_pilot.py" \
  --input "${IBD_PILOT_INPUT}" \
  --coverage "${IBD_PILOT_COVERAGE}" \
  --outdir "${IBD_PILOT_OUTPUT}" \
  --cells-per-type 200 \
  --holdout-fraction 0.2 \
  --max-per-donor 5 \
  --review-cells-per-pending-label 60 \
  --seed 20261006 \
  2>&1 | tee "${IBD_PILOT_OUTPUT}.log"

echo "Pilot output: ${IBD_PILOT_OUTPUT}"
echo "Download small reports first: ${IBD_PILOT_OUTPUT}/small_audit_reports.tar.gz"
echo "Web input gzip: ${IBD_PILOT_OUTPUT}/IBDverse_TI_18types_scRNA_reference.txt.gz"
