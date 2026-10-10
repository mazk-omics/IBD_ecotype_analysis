#!/usr/bin/env bash

set -Eeuo pipefail

# This launcher reads existing BayesPrism/CORE12 results and writes only:
#   - a new timestamped QC directory under combined_reference/qc/
#   - a new timestamped log under combined_reference/logs/
# It does not modify formal analysis outputs or start/stop EcoTyper.

readonly ROOT="/home/mazekai/IBD_EcoTyper"
readonly R_SCRIPT="${ROOT}/scripts/R/bayesprism_epithelial/audit_bayesprism_qnt_consistency_20260925.R"
readonly LOG_DIR="${ROOT}/03_reference/combined_reference/logs/bayesprism_qnt_audit"
readonly RUN_ID="$(date '+%Y%m%d_%H%M%S')"
readonly LOG_FILE="${LOG_DIR}/qnt_audit_${RUN_ID}.log"

mkdir -p "${LOG_DIR}"

if [[ ! -r "${R_SCRIPT}" ]]; then
  printf '[ERROR] R audit script is missing or unreadable: %s\n' "${R_SCRIPT}" >&2
  exit 1
fi

if ! command -v conda >/dev/null 2>&1; then
  for candidate in \
    "${HOME}/miniconda3/bin/conda" \
    "${HOME}/anaconda3/bin/conda" \
    "/opt/conda/bin/conda"; do
    if [[ -x "${candidate}" ]]; then
      export PATH="$(dirname "${candidate}"):${PATH}"
      break
    fi
  done
fi

if ! command -v conda >/dev/null 2>&1; then
  printf '[ERROR] conda was not found.\n' >&2
  exit 1
fi

CONDA_BASE="$(conda info --base)"
# shellcheck source=/dev/null
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate bayesprism

# Mirror stdout/stderr to the terminal and append it to one timestamped log.
exec > >(tee -a "${LOG_FILE}") 2>&1

printf '============================================================\n'
printf 'BayesPrism q_nt consistency audit\n'
printf 'Started: %s\n' "$(date --iso-8601=seconds)"
printf 'Host: %s\n' "$(hostname)"
printf 'Conda environment: %s\n' "${CONDA_DEFAULT_ENV:-unknown}"
printf 'R script: %s\n' "${R_SCRIPT}"
printf 'Log: %s\n' "${LOG_FILE}"
printf '============================================================\n'

if Rscript --vanilla "${R_SCRIPT}"; then
  status=0
  printf '[PASS] Audit completed at %s\n' "$(date --iso-8601=seconds)"
else
  status=$?
  printf '[FAIL] Audit exited with status %s at %s\n' \
    "${status}" "$(date --iso-8601=seconds)" >&2
fi

exit "${status}"
