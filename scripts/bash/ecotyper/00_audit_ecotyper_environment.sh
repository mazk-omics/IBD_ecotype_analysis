#!/usr/bin/env bash

set -u

ROOT="/home/mazekai/IBD_EcoTyper"
REPO="${ROOT}/software/ecotyper"
CONDA="/home/mazekai/miniconda3/bin/conda"

section () {
    echo
    echo "================================================================"
    echo "$1"
    echo "================================================================"
}

section "1. BASIC SYSTEM"

echo "[date]"
date

echo
echo "[hostname]"
hostname

echo
echo "[kernel]"
uname -a

echo
echo "[OS]"
if [[ -f /etc/os-release ]]; then
    cat /etc/os-release
fi

echo
echo "[glibc]"
ldd --version 2>/dev/null | head -n 1 || true

echo
echo "[CPU]"
echo "nproc=$(nproc)"

echo
echo "[memory]"
free -h || true

echo
echo "[disk]"
df -h /home / 2>/dev/null || true


section "2. CONDA"

echo "[conda executable]"
if [[ -x "$CONDA" ]]; then
    echo "$CONDA"
    "$CONDA" --version
else
    echo "ERROR: conda not found at $CONDA"
fi

echo
echo "[conda base]"
"$CONDA" info --base 2>/dev/null || true

echo
echo "[conda environments]"
"$CONDA" env list 2>/dev/null || true

echo
echo "[current shell environment]"
echo "CONDA_DEFAULT_ENV=${CONDA_DEFAULT_ENV:-<none>}"
echo "CONDA_PREFIX=${CONDA_PREFIX:-<none>}"

echo
echo "[current R if any]"
command -v R || true
R --version 2>/dev/null | head -n 2 || true

echo
echo "[current Rscript if any]"
command -v Rscript || true
Rscript --version 2>/dev/null || true


section "3. ECOTYPER REPOSITORY"

echo "[repo path]"
echo "$REPO"

if [[ ! -d "$REPO" ]]; then
    echo "ERROR: repository directory does not exist."
    exit 1
fi

if [[ ! -d "$REPO/.git" ]]; then
    echo "ERROR: $REPO is not a git repository."
    exit 1
fi

echo
echo "[remote]"
git -C "$REPO" remote -v

echo
echo "[branch]"
git -C "$REPO" branch --show-current

echo
echo "[HEAD]"
git -C "$REPO" rev-parse HEAD

echo
echo "[latest local commit]"
git -C "$REPO" log -1 \
    --pretty=format:'commit=%H%ncommit_date=%cI%nsubject=%s'

echo
echo
echo "[working tree]"
git -C "$REPO" status --short

echo
echo "[branch status]"
git -C "$REPO" status -sb

echo
echo "[official master HEAD, network read-only]"
git ls-remote \
    https://github.com/digitalcytometry/ecotyper.git \
    refs/heads/master \
    2>/dev/null || echo "WARNING: could not query GitHub"


section "4. REQUIRED ECOTYPER FILES"

FILES=(
    "EcoTyper_discovery_presorted.R"
    "config_discovery_presorted.yml"
    "pipeline/lib/config.R"
    "pipeline/lib/misc.R"
    "pipeline/lib/multithreading.R"
    "pipeline/state_discovery_presorted_filter_genes.R"
    "pipeline/state_discovery_extract_features.R"
    "pipeline/state_discovery_NMF.R"
    "pipeline/state_discovery_combine_NMF_restarts.R"
    "pipeline/state_discovery_rank_selection.R"
    "pipeline/state_discovery_initial_plots.R"
    "pipeline/state_discovery_first_filter.R"
    "pipeline/ecotypes.R"
    "pipeline/ecotypes_assign_samples.R"
)

for f in "${FILES[@]}"; do
    if [[ -f "$REPO/$f" ]]; then
        printf "PASS  %s\n" "$f"
    else
        printf "FAIL  %s\n" "$f"
    fi
done


section "5. OFFICIAL TUTORIAL 6 INPUT"

TUTORIAL="${REPO}/example_data/Tutorial_6/PresortedDiscovery"

if [[ -d "$TUTORIAL" ]]; then
    echo "PASS: $TUTORIAL"
    echo
    echo "[files]"
    find "$TUTORIAL" -maxdepth 1 -type f \
        -printf '%f\n' | sort
else
    echo "FAIL: Tutorial 6 input directory missing"
fi

echo
echo "[annotation]"
ls -lh \
    "${REPO}/example_data/Tutorial_6/PresortedDiscovery_annotation.txt" \
    2>/dev/null || true


section "6. R PACKAGE REFERENCES IN PRESORTED PIPELINE"

{
    grep -hE '^[[:space:]]*library\(' \
        "$REPO/EcoTyper_discovery_presorted.R" \
        "$REPO/pipeline/"*.R \
        "$REPO/pipeline/lib/"*.R \
        2>/dev/null || true
} \
    | sed -E 's/.*library\(([^,)]+).*/\1/' \
    | tr -d "\"'" \
    | sort -u


section "7. PRESORTED ENTRYPOINT HEADER"

sed -n '1,100p' \
    "$REPO/EcoTyper_discovery_presorted.R"


section "8. PRESORTED CONFIG"

cat "$REPO/config_discovery_presorted.yml"


section "9. CHECK PRESORTED CONFIG FUNCTION"

grep -n \
    -A 90 \
    'check_discovery_configuration_presorted' \
    "$REPO/pipeline/lib/config.R" \
    2>/dev/null || true


section "10. NMF DEPENDENCY CHECK IN SOURCE"

grep -nE \
    'library\(NMF\)|library\(doParallel\)|nmf\(' \
    "$REPO/pipeline/state_discovery_NMF.R" \
    2>/dev/null || true


section "11. IMPORTANT SOURCE ORDER CHECK"

echo "Checking state_discovery_extract_features.R around:"
echo "  expression_top_genes_log2"
echo "  get_variable_genes"

grep -nE \
    'expression_top_genes_log2|get_variable_genes|expression_full_matrix_log2' \
    "$REPO/pipeline/state_discovery_extract_features.R" \
    2>/dev/null || true


section "12. SHELL COMMAND DEPENDENCIES"

for cmd in bash git ln cp mv mkdir sed grep awk sort find; do
    printf "%-10s : " "$cmd"
    command -v "$cmd" || echo "MISSING"
done


section "13. AUDIT SUMMARY"

echo "Repository : $REPO"
echo "Current HEAD:"
git -C "$REPO" rev-parse HEAD

echo
echo "Audit completed."
echo "NO environments were created or modified."
echo "NO repository files were modified."
