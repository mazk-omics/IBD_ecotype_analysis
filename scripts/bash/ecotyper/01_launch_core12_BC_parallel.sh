#!/usr/bin/env bash

set -euo pipefail

ROOT="/home/mazekai/IBD_EcoTyper"
ECOTYPER="${ROOT}/software/ecotyper"
CONDA_SH="/home/mazekai/miniconda3/etc/profile.d/conda.sh"

LOG_DIR="${ROOT}/05_EcoTyper/logs/MSCCR_CORE12"

CFG_B="${ROOT}/05_EcoTyper/config/MSCCR_CORE12/ctspec_region_all/config_discovery_presorted.yml"
CFG_C="${ROOT}/05_EcoTyper/config/MSCCR_CORE12/ctspec_region_ibd/config_discovery_presorted.yml"

SESSION_B="eco_B_region_all"
SESSION_C="eco_C_region_ibd"

DATASET_B="MSCCR_CORE12_CTSPEC_REGION_ALL"
DATASET_C="MSCCR_CORE12_CTSPEC_REGION_IBD"

RUNTIME_B="${ECOTYPER}/EcoTyper/${DATASET_B}"
RUNTIME_C="${ECOTYPER}/EcoTyper/${DATASET_C}"

LOG_B="${LOG_DIR}/ctspec_region_all.log"
LOG_C="${LOG_DIR}/ctspec_region_ibd.log"

mkdir -p "${LOG_DIR}"

section () {
    echo
    echo "============================================================"
    echo "$1"
    echo "============================================================"
}

# ============================================================
# 1. PRE-FLIGHT
# ============================================================

section "1. PRE-FLIGHT"

command -v tmux >/dev/null 2>&1 || {
    echo "[FAIL] tmux not found"
    exit 1
}

[[ -f "${CFG_B}" ]] || {
    echo "[FAIL] Missing B config: ${CFG_B}"
    exit 1
}

[[ -f "${CFG_C}" ]] || {
    echo "[FAIL] Missing C config: ${CFG_C}"
    exit 1
}

[[ -f "${CONDA_SH}" ]] || {
    echo "[FAIL] Missing conda shell initialization: ${CONDA_SH}"
    exit 1
}

# Require the intended thread count.
grep -Eq 'Number of threads[[:space:]]*:[[:space:]]*18' "${CFG_B}" || {
    echo "[FAIL] B config is not set to 18 threads"
    exit 1
}

grep -Eq 'Number of threads[[:space:]]*:[[:space:]]*18' "${CFG_C}" || {
    echo "[FAIL] C config is not set to 18 threads"
    exit 1
}

echo "[PASS] Both configs use 18 threads"

# Avoid accidentally launching duplicate tmux jobs.
if tmux has-session -t "${SESSION_B}" 2>/dev/null; then
    echo "[FAIL] tmux session already exists: ${SESSION_B}"
    exit 1
fi

if tmux has-session -t "${SESSION_C}" 2>/dev/null; then
    echo "[FAIL] tmux session already exists: ${SESSION_C}"
    exit 1
fi

# Avoid silently resuming/contaminating a previous formal run.
if [[ -d "${RUNTIME_B}" ]] && [[ -n "$(find "${RUNTIME_B}" -mindepth 1 -print -quit 2>/dev/null)" ]]; then
    echo "[FAIL] Non-empty B runtime directory already exists:"
    echo "       ${RUNTIME_B}"
    exit 1
fi

if [[ -d "${RUNTIME_C}" ]] && [[ -n "$(find "${RUNTIME_C}" -mindepth 1 -print -quit 2>/dev/null)" ]]; then
    echo "[FAIL] Non-empty C runtime directory already exists:"
    echo "       ${RUNTIME_C}"
    exit 1
fi

echo "[PASS] No pre-existing B/C formal runtime detected"

# System snapshot.
echo
echo "[CPU]"
nproc

echo
echo "[memory]"
free -h

echo
echo "[disk]"
df -h /home

echo
echo "[current EcoTyper/NMF processes]"
ps -ef | grep -E 'EcoTyper_discovery_presorted|state_discovery_NMF' \
    | grep -v grep || true


# ============================================================
# 2. COMMAND BUILDER
# ============================================================

make_command () {

    local config="$1"
    local logfile="$2"
    local label="$3"

    cat <<EOF
source "${CONDA_SH}"
conda activate ecotyper

cd "${ECOTYPER}"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1

echo "============================================================"
echo "${label}"
echo "START: \$(date --iso-8601=seconds)"
echo "HOST: \$(hostname)"
echo "CONFIG: ${config}"
echo "============================================================"

set -o pipefail

Rscript EcoTyper_discovery_presorted.R \
  -c "${config}" \
  2>&1 | tee "${logfile}"

rc=\${PIPESTATUS[0]}

echo
echo "============================================================"
echo "${label}"
echo "END: \$(date --iso-8601=seconds)"
echo "EXIT CODE: \$rc"
echo "============================================================"

if [[ \$rc -eq 0 ]]; then
    echo "FINAL STATUS: SUCCESS"
else
    echo "FINAL STATUS: FAILED"
fi

exec bash
EOF
}


# ============================================================
# 3. START B
# ============================================================

section "2. START B: ALL + REGION3"

CMD_B="$(make_command \
    "${CFG_B}" \
    "${LOG_B}" \
    "B: MSCCR CORE12 ALL + REGION3")"

tmux new-session \
    -d \
    -s "${SESSION_B}" \
    "bash -lc $(printf '%q' "${CMD_B}")"

echo "[PASS] Started: ${SESSION_B}"
echo "       log: ${LOG_B}"


# ============================================================
# 4. STAGGER
# ============================================================

section "3. STAGGER"

echo "Waiting 10 minutes before starting C..."
sleep 600


# ============================================================
# 5. CHECK B STILL EXISTS
# ============================================================

if tmux has-session -t "${SESSION_B}" 2>/dev/null; then
    echo "[PASS] B tmux session still exists"
else
    echo "[WARNING] B tmux session no longer exists."
    echo "Review log before interpreting this as success or failure:"
    echo "${LOG_B}"
fi


# ============================================================
# 6. START C
# ============================================================

section "4. START C: IBD ONLY + REGION3"

CMD_C="$(make_command \
    "${CFG_C}" \
    "${LOG_C}" \
    "C: MSCCR CORE12 IBD ONLY + REGION3")"

tmux new-session \
    -d \
    -s "${SESSION_C}" \
    "bash -lc $(printf '%q' "${CMD_C}")"

echo "[PASS] Started: ${SESSION_C}"
echo "       log: ${LOG_C}"


# ============================================================
# 7. FINAL STATUS
# ============================================================

section "5. ACTIVE TMUX SESSIONS"

tmux ls || true

echo
echo "Expected:"
echo "  ${SESSION_B}"
echo "  ${SESSION_C}"

echo
echo "Monitor with:"
echo "  tmux attach -t ${SESSION_B}"
echo "  tmux attach -t ${SESSION_C}"
echo
echo "or:"
echo "  tail -f ${LOG_B}"
echo "  tail -f ${LOG_C}"

echo
echo "LAUNCH STATUS: B_AND_C_SUBMITTED"
