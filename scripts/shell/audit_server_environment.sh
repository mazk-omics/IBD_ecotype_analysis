#!/usr/bin/env bash

set -uo pipefail

# ============================================================
# IBD EcoTyper Server Environment & Data Audit
#
# Scan scope:
#   ONLY /home/mazekai/
#
# Purpose:
#   1. Audit current computational environment
#   2. Discover existing Python / Conda / R environments
#   3. Locate IBD EcoTyper-related datasets and software
#   4. Assess readiness for GSE282122 Reference Audit v2
#
# This script DOES NOT:
#   - install software
#   - modify Conda/Mamba environments
#   - modify shell configuration
#   - modify raw/project data
#   - use sudo
#   - download anything
#   - load large h5ad/RDS matrices into memory
#
# The only files created are audit reports under:
#   ~/IBD_EcoTyper/logs/server_audit/
#
# Usage:
#   bash audit_server_environment.sh
#
# ============================================================


# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------

SCAN_ROOT="/home/mazekai"
PROJECT_ROOT="/home/mazekai/IBD_EcoTyper"

TIMESTAMP="$(date '+%Y%m%d_%H%M%S')"

AUDIT_ROOT="${PROJECT_ROOT}/logs/server_audit/${TIMESTAMP}"

REPORT="${AUDIT_ROOT}/server_environment_report.txt"
DATA_REPORT="${AUDIT_ROOT}/project_data_inventory.tsv"
CONDA_REPORT="${AUDIT_ROOT}/conda_environments.tsv"
R_REPORT="${AUDIT_ROOT}/r_installations.tsv"
PYTHON_REPORT="${AUDIT_ROOT}/python_installations.tsv"
SOFTWARE_REPORT="${AUDIT_ROOT}/relevant_software_paths.tsv"
OMICS_REPORT="${AUDIT_ROOT}/relevant_omics_files.tsv"

mkdir -p "${AUDIT_ROOT}"


# ============================================================
# Helper functions
# ============================================================

section() {
    {
        echo
        echo "============================================================"
        echo "$1"
        echo "============================================================"
    } >> "${REPORT}"
}


subsection() {
    {
        echo
        echo "------------------------------------------------------------"
        echo "$1"
        echo "------------------------------------------------------------"
    } >> "${REPORT}"
}


cmd_exists() {
    command -v "$1" >/dev/null 2>&1
}


print_kv() {
    printf "%-30s : %s\n" "$1" "$2" >> "${REPORT}"
}


file_size_bytes() {
    stat -c '%s' "$1" 2>/dev/null || echo "NA"
}


file_size_human() {
    du -h "$1" 2>/dev/null | awk '{print $1}' || echo "NA"
}


# ============================================================
# Header
# ============================================================

{
    echo "IBD EcoTyper Server Environment & Data Audit"
    echo
    echo "Audit time      : $(date -Is 2>/dev/null || date)"
    echo "User scan root  : ${SCAN_ROOT}"
    echo "Project root    : ${PROJECT_ROOT}"
    echo "Audit directory : ${AUDIT_ROOT}"
    echo
    echo "IMPORTANT:"
    echo "Recursive file discovery is restricted to /home/mazekai/"
    echo "No other users' home directories or server-wide filesystems are scanned."
    echo
    echo "READ-ONLY AUDIT"
} > "${REPORT}"


# ============================================================
# A. Basic system
# ============================================================

section "A. BASIC SYSTEM"

print_kv "Hostname" "$(hostname 2>/dev/null || echo NA)"
print_kv "User" "$(whoami 2>/dev/null || echo NA)"
print_kv "UID" "$(id -u 2>/dev/null || echo NA)"
print_kv "Home" "${HOME}"
print_kv "Current directory" "$(pwd)"
print_kv "Shell" "${SHELL:-NA}"
print_kv "Architecture" "$(uname -m 2>/dev/null || echo NA)"
print_kv "Kernel" "$(uname -srmo 2>/dev/null || echo NA)"


if [[ -r /etc/os-release ]]; then

    subsection "Operating system"

    cat /etc/os-release >> "${REPORT}"

fi


subsection "CPU"

if cmd_exists lscpu; then

    lscpu >> "${REPORT}" 2>&1

else

    grep -E 'model name|processor|cpu cores' \
        /proc/cpuinfo \
        | head -100 \
        >> "${REPORT}" 2>&1 || true

fi


subsection "Memory"

if cmd_exists free; then

    free -h >> "${REPORT}" 2>&1

else

    cat /proc/meminfo >> "${REPORT}" 2>&1 || true

fi


subsection "Process limits"

ulimit -a >> "${REPORT}" 2>&1 || true


# ============================================================
# B. Storage
# ============================================================

section "B. STORAGE"

subsection "Filesystem containing /home/mazekai"

df -hT "${SCAN_ROOT}" >> "${REPORT}" 2>&1 || true


subsection "Filesystem containing project"

df -hT "${PROJECT_ROOT}" >> "${REPORT}" 2>&1 || true


subsection "User home usage"

du -sh "${SCAN_ROOT}" >> "${REPORT}" 2>&1 || true


subsection "Top-level usage under /home/mazekai"

du -sh "${SCAN_ROOT}"/* \
    2>/dev/null \
    | sort -h \
    >> "${REPORT}" || true


# ============================================================
# C. Project directory
# ============================================================

section "C. PROJECT DIRECTORY"

if [[ -d "${PROJECT_ROOT}" ]]; then

    print_kv "Project exists" "YES"

    print_kv \
        "Readable" \
        "$([[ -r "${PROJECT_ROOT}" ]] && echo YES || echo NO)"

    print_kv \
        "Writable" \
        "$([[ -w "${PROJECT_ROOT}" ]] && echo YES || echo NO)"

    print_kv \
        "Searchable" \
        "$([[ -x "${PROJECT_ROOT}" ]] && echo YES || echo NO)"


    subsection "Project directory structure"

    find "${PROJECT_ROOT}" \
        -mindepth 1 \
        -maxdepth 3 \
        -type d \
        -print \
        2>/dev/null \
        | sort \
        >> "${REPORT}"


    subsection "Git status"

    if [[ -d "${PROJECT_ROOT}/.git" ]]; then

        (
            cd "${PROJECT_ROOT}" || exit

            echo "Git repository: YES"
            echo

            git status --short --branch 2>&1 || true

            echo
            echo "Remote:"
            git remote -v 2>&1 || true

            echo
            echo "Recent commits:"
            git log -5 --oneline 2>&1 || true

        ) >> "${REPORT}"

    else

        echo "Git repository: NO" >> "${REPORT}"

    fi

else

    print_kv "Project exists" "NO"

fi


# ============================================================
# D. Shell / PATH
# ============================================================

section "D. SHELL AND PATH"

print_kv "SHELL" "${SHELL:-NA}"
print_kv "PATH" "${PATH:-NA}"


subsection "Relevant shell initialization entries"

for f in \
    "${HOME}/.bashrc" \
    "${HOME}/.bash_profile" \
    "${HOME}/.profile" \
    "${HOME}/.zshrc"
do

    if [[ -r "$f" ]]; then

        echo
        echo "### ${f}" >> "${REPORT}"

        grep -nEi \
            'conda|mamba|micromamba|module|R_HOME|R_LIBS|cuda|singularity|apptainer' \
            "$f" \
            >> "${REPORT}" 2>/dev/null || true

    fi

done


# ============================================================
# E. Current Python
# ============================================================

section "E. CURRENT PYTHON ENVIRONMENT"

for exe in python python3 pip pip3; do

    if cmd_exists "$exe"; then

        print_kv "${exe} path" "$(command -v "$exe")"

        "$exe" --version >> "${REPORT}" 2>&1 || true

    else

        print_kv "$exe" "NOT FOUND"

    fi

done


PYTHON_CMD=""

if cmd_exists python; then

    PYTHON_CMD="python"

elif cmd_exists python3; then

    PYTHON_CMD="python3"

fi


if [[ -n "${PYTHON_CMD}" ]]; then

    subsection "Python runtime"

    "${PYTHON_CMD}" - <<'PY' >> "${REPORT}" 2>&1
import os
import sys
import platform
import importlib.metadata as md

print("Python executable :", sys.executable)
print("Python version    :", sys.version.replace("\n", " "))
print("sys.prefix        :", sys.prefix)
print("base_prefix       :", getattr(sys, "base_prefix", "NA"))
print("Platform          :", platform.platform())
print("CONDA_PREFIX      :", os.environ.get("CONDA_PREFIX", ""))
print("CONDA_DEFAULT_ENV :", os.environ.get("CONDA_DEFAULT_ENV", ""))
print("VIRTUAL_ENV       :", os.environ.get("VIRTUAL_ENV", ""))

packages = [
    "numpy",
    "pandas",
    "scipy",
    "pyarrow",
    "anndata",
    "scanpy",
    "matplotlib",
    "h5py",
    "tables",
]

print()
print("Selected package versions:")

for pkg in packages:
    try:
        print(f"{pkg:15s} {md.version(pkg)}")
    except Exception:
        print(f"{pkg:15s} NOT INSTALLED")
PY

fi


# ============================================================
# F. Python installations under /home/mazekai
# ============================================================

section "F. PYTHON INSTALLATIONS UNDER /home/mazekai"

printf "path\tversion\n" > "${PYTHON_REPORT}"


while IFS= read -r -d '' python_path; do

    [[ -x "${python_path}" ]] || continue

    version="$(
        "${python_path}" --version 2>&1 \
        | head -1 \
        | tr '\t' ' '
    )"

    printf "%s\t%s\n" \
        "${python_path}" \
        "${version}" \
        >> "${PYTHON_REPORT}"

done < <(

    find "${SCAN_ROOT}" \
        \( \
            -path '*/.cache' \
            -o -path '*/.cache/*' \
            -o -path '*/.git' \
            -o -path '*/.git/*' \
            -o -path '*/node_modules' \
            -o -path '*/node_modules/*' \
        \) -prune \
        -o \
        -type f \
        \( \
            -name python \
            -o -name python3 \
        \) \
        -perm -u+x \
        -print0 \
        2>/dev/null

)


sort -u "${PYTHON_REPORT}" -o "${PYTHON_REPORT}"


# ============================================================
# G. Conda / Mamba
# ============================================================

section "G. CONDA / MAMBA"

for exe in conda mamba micromamba; do

    if cmd_exists "$exe"; then

        print_kv "${exe} path" "$(command -v "$exe")"

        "$exe" --version >> "${REPORT}" 2>&1 || true

    else

        print_kv "$exe" "NOT FOUND"

    fi

done


printf \
"environment\tpath\tpython_version\tr_version\n" \
> "${CONDA_REPORT}"


if cmd_exists conda; then

    subsection "conda info"

    conda info >> "${REPORT}" 2>&1 || true


    subsection "conda channels"

    conda config --show channels \
        >> "${REPORT}" 2>&1 || true


    subsection "conda env list"

    conda env list \
        >> "${REPORT}" 2>&1 || true


    while IFS= read -r env_path; do

        [[ -z "${env_path}" ]] && continue
        [[ ! -d "${env_path}" ]] && continue

        env_name="$(basename "${env_path}")"

        python_ver="NA"
        r_ver="NA"

        if [[ -x "${env_path}/bin/python" ]]; then

            python_ver="$(
                "${env_path}/bin/python" --version \
                2>&1 | head -1
            )"

        fi

        if [[ -x "${env_path}/bin/R" ]]; then

            r_ver="$(
                "${env_path}/bin/R" --version \
                2>&1 | head -1
            )"

        fi

        printf "%s\t%s\t%s\t%s\n" \
            "${env_name}" \
            "${env_path}" \
            "${python_ver}" \
            "${r_ver}" \
            >> "${CONDA_REPORT}"

    done < <(

        conda env list 2>/dev/null \
            | awk '
                /^[[:space:]]*#/ {next}
                NF == 0 {next}
                {print $NF}
            ' \
            | grep '^/' \
            | sort -u

    )

fi


# ============================================================
# H. Current R
# ============================================================

section "H. CURRENT R ENVIRONMENT"

for exe in R Rscript; do

    if cmd_exists "$exe"; then

        print_kv "${exe} path" "$(command -v "$exe")"

        "$exe" --version >> "${REPORT}" 2>&1 || true

    else

        print_kv "$exe" "NOT FOUND"

    fi

done


if cmd_exists Rscript; then

    subsection "Current R runtime"

    Rscript - <<'RS' >> "${REPORT}" 2>&1

cat("R.version.string:\n")
cat(R.version.string, "\n\n")

cat(".libPaths():\n")
print(.libPaths())

cat("\nR_HOME:\n")
cat(R.home(), "\n")

cat("\nSelected package versions:\n")

pkgs <- c(
    "BiocManager",
    "NMF",
    "data.table",
    "Matrix",
    "Seurat"
)

for (p in pkgs) {

    if (requireNamespace(p, quietly = TRUE)) {

        cat(
            sprintf(
                "%-15s %s\n",
                p,
                as.character(packageVersion(p))
            )
        )

    } else {

        cat(
            sprintf(
                "%-15s NOT INSTALLED\n",
                p
            )
        )

    }
}

RS

fi


# ============================================================
# I. R installations under /home/mazekai
# ============================================================

section "I. R INSTALLATIONS UNDER /home/mazekai"

printf "path\tversion\n" > "${R_REPORT}"


while IFS= read -r -d '' rpath; do

    [[ -x "${rpath}" ]] || continue

    version="$(
        "${rpath}" --version 2>&1 \
        | head -1 \
        | tr '\t' ' '
    )"

    printf "%s\t%s\n" \
        "${rpath}" \
        "${version}" \
        >> "${R_REPORT}"

done < <(

    find "${SCAN_ROOT}" \
        \( \
            -path '*/.cache' \
            -o -path '*/.cache/*' \
            -o -path '*/.git' \
            -o -path '*/.git/*' \
        \) -prune \
        -o \
        -type f \
        \( \
            -name R \
            -o -name Rscript \
        \) \
        -perm -u+x \
        -print0 \
        2>/dev/null

)


sort -u "${R_REPORT}" -o "${R_REPORT}"


# ============================================================
# J. Environment modules
# ============================================================

section "J. ENVIRONMENT MODULES"

if type module >/dev/null 2>&1; then

    echo "Environment Modules detected." >> "${REPORT}"


    subsection "module list"

    module list >> "${REPORT}" 2>&1 || true


    subsection "Relevant available modules"

    module avail 2>&1 \
        | grep -Ei \
        'python|(^|/)r([/-]|$)|conda|anaconda|cuda|singularity|apptainer' \
        | head -500 \
        >> "${REPORT}" || true

else

    echo "Environment Modules not detected." \
        >> "${REPORT}"

fi


# ============================================================
# K. HPC scheduler
# ============================================================

section "K. HPC SCHEDULER"

if cmd_exists sinfo; then

    echo "SLURM detected." >> "${REPORT}"


    subsection "sinfo"

    sinfo >> "${REPORT}" 2>&1 || true


    subsection "Current user's jobs"

    squeue -u "$(whoami)" \
        >> "${REPORT}" 2>&1 || true


    subsection "SLURM environment variables"

    env \
        | grep '^SLURM_' \
        | sort \
        >> "${REPORT}" 2>/dev/null || true


elif cmd_exists qstat; then

    echo "PBS/Torque-like scheduler detected." \
        >> "${REPORT}"

    qstat >> "${REPORT}" 2>&1 || true


elif cmd_exists bjobs; then

    echo "LSF detected." \
        >> "${REPORT}"

    bjobs >> "${REPORT}" 2>&1 || true


else

    echo "No common HPC scheduler detected." \
        >> "${REPORT}"

fi


# ============================================================
# L. GPU / CUDA
# ============================================================

section "L. GPU / CUDA"

if cmd_exists nvidia-smi; then

    nvidia-smi >> "${REPORT}" 2>&1 || true

else

    echo "nvidia-smi not found." \
        >> "${REPORT}"

fi


if cmd_exists nvcc; then

    subsection "CUDA compiler"

    nvcc --version \
        >> "${REPORT}" 2>&1 || true

fi


# ============================================================
# M. Container runtimes
# ============================================================

section "M. CONTAINER RUNTIMES"

for exe in \
    docker \
    podman \
    singularity \
    apptainer
do

    if cmd_exists "$exe"; then

        print_kv "${exe} path" "$(command -v "$exe")"

        "$exe" --version \
            >> "${REPORT}" 2>&1 || true

    else

        print_kv "$exe" "NOT FOUND"

    fi

done


# ============================================================
# N. EcoTyper / CIBERSORTx search under /home/mazekai
# ============================================================

section "N. ECOTYPER / CIBERSORTx SEARCH"

printf "path\ttype\n" > "${SOFTWARE_REPORT}"


while IFS= read -r -d '' path; do

    if [[ -d "${path}" ]]; then

        object_type="DIRECTORY"

    else

        object_type="FILE"

    fi

    printf "%s\t%s\n" \
        "${path}" \
        "${object_type}" \
        >> "${SOFTWARE_REPORT}"

done < <(

    find "${SCAN_ROOT}" \
        \( \
            -path '*/.cache' \
            -o -path '*/.cache/*' \
            -o -path '*/.git' \
            -o -path '*/.git/*' \
        \) -prune \
        -o \
        \( -type f -o -type d \) \
        \( \
            -iname '*ecotyper*' \
            -o -iname '*cibersortx*' \
            -o -iname '*cibersort*' \
        \) \
        -print0 \
        2>/dev/null

)


sort -u "${SOFTWARE_REPORT}" -o "${SOFTWARE_REPORT}"


# ============================================================
# O. Project-relevant data discovery
# ============================================================

section "O. PROJECT-RELEVANT DATA DISCOVERY"

printf \
"path\tsize_bytes\tsize_human\tobject_type\n" \
> "${DATA_REPORT}"


declare -A SEEN_DATA


while IFS= read -r -d '' path; do

    [[ -z "${path}" ]] && continue

    if [[ -n "${SEEN_DATA["${path}"]+x}" ]]; then
        continue
    fi

    SEEN_DATA["${path}"]=1


    if [[ -f "${path}" ]]; then

        bytes="$(file_size_bytes "${path}")"
        human="$(file_size_human "${path}")"

        printf "%s\t%s\t%s\tFILE\n" \
            "${path}" \
            "${bytes}" \
            "${human}" \
            >> "${DATA_REPORT}"


    elif [[ -d "${path}" ]]; then

        printf "%s\tNA\tNA\tDIRECTORY\n" \
            "${path}" \
            >> "${DATA_REPORT}"

    fi


done < <(

    find "${SCAN_ROOT}" \
        \( \
            -path '*/.cache' \
            -o -path '*/.cache/*' \
            -o -path '*/.git' \
            -o -path '*/.git/*' \
            -o -path '*/node_modules' \
            -o -path '*/node_modules/*' \
        \) -prune \
        -o \
        \( -type f -o -type d \) \
        \( \
            -iname '*GSE282122*' \
            -o -iname '*TAURUS*' \
            -o -iname '*MDV*' \
            -o -iname '*SCP1884*' \
            -o -iname '*SCP259*' \
            -o -iname '*GSE193677*' \
            -o -iname '*MSCCR*' \
            -o -iname '*reference*manifest*' \
            -o -iname '*sample*manifest*' \
            -o -iname '*cell*metadata*' \
            -o -iname '*supplementary*table*' \
        \) \
        -print0 \
        2>/dev/null

)


sort -u "${DATA_REPORT}" -o "${DATA_REPORT}"


# ============================================================
# P. Relevant omics files under /home/mazekai
# ============================================================

section "P. RELEVANT OMICS FILE INVENTORY"

printf \
"path\tsize_bytes\tsize_human\n" \
> "${OMICS_REPORT}"


declare -A SEEN_OMICS


while IFS= read -r -d '' path; do

    [[ -z "${path}" ]] && continue

    if [[ -n "${SEEN_OMICS["${path}"]+x}" ]]; then
        continue
    fi

    SEEN_OMICS["${path}"]=1


    bytes="$(file_size_bytes "${path}")"
    human="$(file_size_human "${path}")"


    printf "%s\t%s\t%s\n" \
        "${path}" \
        "${bytes}" \
        "${human}" \
        >> "${OMICS_REPORT}"


done < <(

    find "${SCAN_ROOT}" \
        \( \
            -path '*/.cache' \
            -o -path '*/.cache/*' \
            -o -path '*/.git' \
            -o -path '*/.git/*' \
        \) -prune \
        -o \
        -type f \
        \( \
            -iname '*.h5ad' \
            -o -iname '*.h5' \
            -o -iname '*.loom' \
            -o -iname '*.rds' \
            -o -iname '*.RDS' \
            -o -iname '*.mtx' \
            -o -iname '*.mtx.gz' \
            -o -iname '*.parquet' \
            -o -iname '*.csv.gz' \
            -o -iname '*.tsv.gz' \
            -o -iname '*.txt.gz' \
        \) \
        -print0 \
        2>/dev/null

)


sort -u "${OMICS_REPORT}" -o "${OMICS_REPORT}"


# ============================================================
# Q. GSE282122-specific readiness
# ============================================================

section "Q. GSE282122 REFERENCE AUDIT READINESS"


GSE_COUNT="$(
    grep -Eic \
        'GSE282122|TAURUS|MDV' \
        "${DATA_REPORT}" \
        2>/dev/null \
        || true
)"


print_kv \
    "GSE282122/TAURUS/MDV matches" \
    "${GSE_COUNT}"


subsection "Potential GSE282122-related paths"

grep -Ei \
    'GSE282122|TAURUS|MDV' \
    "${DATA_REPORT}" \
    >> "${REPORT}" 2>/dev/null || true


subsection "Potential GSE282122 omics files"

grep -Ei \
    'GSE282122|TAURUS|MDV' \
    "${OMICS_REPORT}" \
    >> "${REPORT}" 2>/dev/null || true


# ============================================================
# R. Reference Audit v2 software readiness
# ============================================================

section "R. REFERENCE AUDIT V2 SOFTWARE READINESS"

PASS_COUNT=0
WARNING_COUNT=0


if [[ -n "${PYTHON_CMD}" ]]; then

    echo "[PASS] Python available: $(command -v "${PYTHON_CMD}")" \
        >> "${REPORT}"

    PASS_COUNT=$((PASS_COUNT + 1))

else

    echo "[MISSING] Python not found." \
        >> "${REPORT}"

    WARNING_COUNT=$((WARNING_COUNT + 1))

fi


if [[ -n "${PYTHON_CMD}" ]] && \
   "${PYTHON_CMD}" -c "import pandas" \
   >/dev/null 2>&1; then

    echo "[PASS] pandas available." \
        >> "${REPORT}"

    PASS_COUNT=$((PASS_COUNT + 1))

else

    echo "[MISSING] pandas unavailable in current Python." \
        >> "${REPORT}"

    WARNING_COUNT=$((WARNING_COUNT + 1))

fi


if [[ -n "${PYTHON_CMD}" ]] && \
   "${PYTHON_CMD}" -c "import pyarrow" \
   >/dev/null 2>&1; then

    echo "[PASS] pyarrow available." \
        >> "${REPORT}"

    PASS_COUNT=$((PASS_COUNT + 1))

else

    echo "[WARNING] pyarrow unavailable in current Python." \
        >> "${REPORT}"

    WARNING_COUNT=$((WARNING_COUNT + 1))

fi


if [[ -n "${PYTHON_CMD}" ]] && \
   "${PYTHON_CMD}" -c "import anndata" \
   >/dev/null 2>&1; then

    echo "[PASS] anndata available." \
        >> "${REPORT}"

else

    echo "[INFO] anndata unavailable in current Python." \
        >> "${REPORT}"

    echo \
        "       Required only if metadata must be extracted from .h5ad." \
        >> "${REPORT}"

fi


if [[ -r "${PROJECT_ROOT}" && -w "${PROJECT_ROOT}" ]]; then

    echo "[PASS] Project directory is readable and writable." \
        >> "${REPORT}"

    PASS_COUNT=$((PASS_COUNT + 1))

else

    echo "[WARNING] Project directory permission issue." \
        >> "${REPORT}"

    WARNING_COUNT=$((WARNING_COUNT + 1))

fi


echo >> "${REPORT}"

print_kv "PASS count" "${PASS_COUNT}"

print_kv \
    "Warnings / missing" \
    "${WARNING_COUNT}"


# ============================================================
# S. Output summary
# ============================================================

section "S. AUDIT OUTPUT FILES"

cat >> "${REPORT}" <<EOF

Main environment report:
  ${REPORT}

Project/data inventory:
  ${DATA_REPORT}

Conda environments:
  ${CONDA_REPORT}

Python installations:
  ${PYTHON_REPORT}

R installations:
  ${R_REPORT}

EcoTyper / CIBERSORTx paths:
  ${SOFTWARE_REPORT}

Relevant omics files:
  ${OMICS_REPORT}

EOF


# ============================================================
# Terminal output
# ============================================================

echo
echo "============================================================"
echo "IBD EcoTyper environment audit completed"
echo "============================================================"
echo
echo "Scan scope:"
echo "  ${SCAN_ROOT}"
echo
echo "Audit directory:"
echo "  ${AUDIT_ROOT}"
echo
echo "Main report:"
echo "  ${REPORT}"
echo
echo "Project/data inventory:"
echo "  ${DATA_REPORT}"
echo
echo "Conda environments:"
echo "  ${CONDA_REPORT}"
echo
echo "Python installations:"
echo "  ${PYTHON_REPORT}"
echo
echo "R installations:"
echo "  ${R_REPORT}"
echo
echo "EcoTyper / CIBERSORTx paths:"
echo "  ${SOFTWARE_REPORT}"
echo
echo "Relevant omics files:"
echo "  ${OMICS_REPORT}"
echo
echo "No software, environment, raw data, or analysis files were modified."
echo