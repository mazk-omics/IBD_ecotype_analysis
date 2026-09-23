#!/usr/bin/env bash

set -eu

BASE="/home/mazekai/IBD_EcoTyper/06_validation/01_bulk_progression/GSE93624"

RAW_DIR="${BASE}/raw"
META_DIR="${BASE}/metadata"
OUT_DIR="${BASE}/prepared/00_inventory"

mkdir -p \
  "${RAW_DIR}" \
  "${META_DIR}" \
  "${OUT_DIR}"

LOG="${OUT_DIR}/00_public_provenance_and_file_inventory.log"

exec > >(tee "${LOG}") 2>&1


echo "============================================================"
echo "GSE93624 Stage 0"
echo "Public provenance + file inventory"
echo "============================================================"
echo


# ============================================================
# 1. GEO URLs
# ============================================================

GEO_BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE93nnn/GSE93624"

SOFT_URL="${GEO_BASE}/soft/GSE93624_family.soft.gz"
MATRIX_URL="${GEO_BASE}/matrix/GSE93624_series_matrix.txt.gz"
RAW_URL="${GEO_BASE}/suppl/GSE93624_RAW.tar"

SOFT="${META_DIR}/GSE93624_family.soft.gz"
MATRIX="${RAW_DIR}/GSE93624_series_matrix.txt.gz"
RAW_TAR="${RAW_DIR}/GSE93624_RAW.tar"


# ============================================================
# 2. Download helper
# ============================================================

download_if_missing () {

    url="$1"
    output="$2"

    if [ -s "${output}" ]; then
        echo "[SKIP] Existing file:"
        echo "       ${output}"
    else
        echo "[DOWNLOAD]"
        echo "  URL : ${url}"
        echo "  OUT : ${output}"

        curl \
          -L \
          --fail \
          --retry 5 \
          --retry-delay 3 \
          --connect-timeout 30 \
          -o "${output}" \
          "${url}"
    fi
}


echo "============================================================"
echo "DOWNLOAD / PRESENCE CHECK"
echo "============================================================"

download_if_missing \
  "${SOFT_URL}" \
  "${SOFT}"

download_if_missing \
  "${MATRIX_URL}" \
  "${MATRIX}"

download_if_missing \
  "${RAW_URL}" \
  "${RAW_TAR}"

echo


# ============================================================
# 3. Basic file inventory
# ============================================================

echo "============================================================"
echo "LOCAL FILE INVENTORY"
echo "============================================================"

find "${BASE}" \
  -maxdepth 3 \
  -type f \
  -printf '%p\t%s bytes\n' \
  | sort

echo


# ============================================================
# 4. Integrity checks
# ============================================================

echo "============================================================"
echo "FILE INTEGRITY"
echo "============================================================"

gzip -t "${SOFT}"
echo "[PASS] SOFT gzip integrity"

gzip -t "${MATRIX}"
echo "[PASS] series-matrix gzip integrity"

tar -tf "${RAW_TAR}" >/dev/null
echo "[PASS] RAW tar integrity"

echo


# ============================================================
# 5. SHA256 provenance hashes
# ============================================================

echo "============================================================"
echo "SHA256"
echo "============================================================"

sha256sum \
  "${SOFT}" \
  "${MATRIX}" \
  "${RAW_TAR}" \
  | tee "${OUT_DIR}/GSE93624_public_file_sha256.tsv"

echo


# ============================================================
# 6. GEO SOFT global provenance
# ============================================================

echo "============================================================"
echo "SOFT SERIES PROVENANCE"
echo "============================================================"

gzip -dc "${SOFT}" \
  | grep -E \
    '^\^SERIES|^!Series_title|^!Series_status|^!Series_submission_date|^!Series_last_update_date|^!Series_summary|^!Series_overall_design|^!Series_platform_id|^!Series_sample_id' \
  > "${OUT_DIR}/GSE93624_series_provenance_excerpt.txt" || true

grep -E \
  '^\^SERIES|^!Series_title|^!Series_status|^!Series_submission_date|^!Series_last_update_date|^!Series_platform_id' \
  "${OUT_DIR}/GSE93624_series_provenance_excerpt.txt" || true

echo


# ============================================================
# 7. Count GEO sample records
# ============================================================

echo "============================================================"
echo "SOFT SAMPLE COUNT"
echo "============================================================"

N_SOFT_SAMPLES=$(
  gzip -dc "${SOFT}" \
  | grep -c '^\^SAMPLE = GSM' || true
)

echo "SOFT sample records: ${N_SOFT_SAMPLES}"

if [ "${N_SOFT_SAMPLES}" -eq 245 ]; then
    echo "[PASS] Expected 245 GEO samples"
else
    echo "[WARN] Expected 245 GEO samples but found ${N_SOFT_SAMPLES}"
fi

echo


# ============================================================
# 8. Sample titles
# ============================================================

echo "============================================================"
echo "SAMPLE TITLE AUDIT"
echo "============================================================"

gzip -dc "${SOFT}" \
  | grep '^!Sample_title = ' \
  | sed 's/^!Sample_title = //' \
  > "${OUT_DIR}/GSE93624_sample_titles.txt"

echo "Sample titles:"
echo "  N = $(wc -l < "${OUT_DIR}/GSE93624_sample_titles.txt")"

echo
echo "First 10:"
head -n 10 \
  "${OUT_DIR}/GSE93624_sample_titles.txt"

echo
echo "Last 10:"
tail -n 10 \
  "${OUT_DIR}/GSE93624_sample_titles.txt"

echo


# ============================================================
# 9. Metadata characteristic-field inventory
# ============================================================

echo "============================================================"
echo "CHARACTERISTIC FIELD INVENTORY"
echo "============================================================"

gzip -dc "${SOFT}" \
  | grep '^!Sample_characteristics_ch1 = ' \
  | sed 's/^!Sample_characteristics_ch1 = //' \
  | sed 's/:.*$//' \
  | sort \
  | uniq -c \
  | sort -nr \
  | tee "${OUT_DIR}/GSE93624_characteristic_field_counts.txt"

echo


# ============================================================
# 10. Important phenotype distributions
# ============================================================

echo "============================================================"
echo "DIAGNOSIS RAW DISTRIBUTION"
echo "============================================================"

gzip -dc "${SOFT}" \
  | grep '^!Sample_characteristics_ch1 = diagnosis:' \
  | sed 's/^!Sample_characteristics_ch1 = diagnosis:[[:space:]]*//' \
  | sort \
  | uniq -c \
  | tee "${OUT_DIR}/GSE93624_raw_diagnosis_counts.txt"

echo


echo "============================================================"
echo "PROGRESSION RAW DISTRIBUTION"
echo "============================================================"

gzip -dc "${SOFT}" \
  | grep '^!Sample_characteristics_ch1 = progression to complication:' \
  | sed 's/^!Sample_characteristics_ch1 = progression to complication:[[:space:]]*//' \
  | sort \
  | uniq -c \
  | tee "${OUT_DIR}/GSE93624_raw_progression_counts.txt"

echo


echo "============================================================"
echo "TISSUE RAW DISTRIBUTION"
echo "============================================================"

gzip -dc "${SOFT}" \
  | grep '^!Sample_characteristics_ch1 = tissue:' \
  | sed 's/^!Sample_characteristics_ch1 = tissue:[[:space:]]*//' \
  | sort \
  | uniq -c \
  | tee "${OUT_DIR}/GSE93624_raw_tissue_counts.txt"

echo


echo "============================================================"
echo "SEX RAW DISTRIBUTION"
echo "============================================================"

gzip -dc "${SOFT}" \
  | grep '^!Sample_characteristics_ch1 = gender:' \
  | sed 's/^!Sample_characteristics_ch1 = gender:[[:space:]]*//' \
  | sort \
  | uniq -c \
  | tee "${OUT_DIR}/GSE93624_raw_gender_counts.txt"

echo


# ============================================================
# 11. GEO processing description
# ============================================================

echo "============================================================"
echo "DATA PROCESSING DESCRIPTION"
echo "============================================================"

gzip -dc "${SOFT}" \
  | grep '^!Sample_data_processing = ' \
  | sed 's/^!Sample_data_processing = //' \
  | sort -u \
  | tee "${OUT_DIR}/GSE93624_data_processing_unique.txt"

echo


# ============================================================
# 12. RAW tar member inventory
# ============================================================

echo "============================================================"
echo "RAW TAR INVENTORY"
echo "============================================================"

tar -tf "${RAW_TAR}" \
  > "${OUT_DIR}/GSE93624_RAW_tar_members.txt"

N_RAW_MEMBERS=$(
  wc -l < "${OUT_DIR}/GSE93624_RAW_tar_members.txt"
)

echo "RAW tar members: ${N_RAW_MEMBERS}"

if [ "${N_RAW_MEMBERS}" -eq 245 ]; then
    echo "[PASS] RAW tar contains 245 members"
else
    echo "[WARN] RAW tar member count is not 245"
fi

echo
echo "First 10 members:"
head -n 10 \
  "${OUT_DIR}/GSE93624_RAW_tar_members.txt"

echo
echo "Last 10 members:"
tail -n 10 \
  "${OUT_DIR}/GSE93624_RAW_tar_members.txt"

echo


# ============================================================
# 13. Check GSM uniqueness from RAW filenames
# ============================================================

echo "============================================================"
echo "RAW GSM ID AUDIT"
echo "============================================================"

sed -E \
  's#.*/##' \
  "${OUT_DIR}/GSE93624_RAW_tar_members.txt" \
  | grep -oE 'GSM[0-9]+' \
  > "${OUT_DIR}/GSE93624_RAW_GSM_ids.txt" || true

N_GSM=$(
  wc -l < "${OUT_DIR}/GSE93624_RAW_GSM_ids.txt"
)

N_GSM_UNIQ=$(
  sort -u "${OUT_DIR}/GSE93624_RAW_GSM_ids.txt" \
  | wc -l
)

echo "GSM IDs in RAW filenames: ${N_GSM}"
echo "Unique GSM IDs:           ${N_GSM_UNIQ}"

if [ "${N_GSM}" -eq 245 ] && [ "${N_GSM_UNIQ}" -eq 245 ]; then
    echo "[PASS] RAW contains 245 unique GSM IDs"
else
    echo "[WARN] RAW GSM count requires inspection"
fi

echo


# ============================================================
# 14. Inspect first / middle / last processed sample file
#     WITHOUT extracting all files
# ============================================================

echo "============================================================"
echo "REPRESENTATIVE SAMPLE-FILE SCHEMA"
echo "============================================================"

TMP="${OUT_DIR}/tmp_schema"
rm -rf "${TMP}"
mkdir -p "${TMP}"

FIRST_MEMBER=$(
  sed -n '1p' \
  "${OUT_DIR}/GSE93624_RAW_tar_members.txt"
)

MIDDLE_MEMBER=$(
  sed -n '123p' \
  "${OUT_DIR}/GSE93624_RAW_tar_members.txt"
)

LAST_MEMBER=$(
  sed -n '245p' \
  "${OUT_DIR}/GSE93624_RAW_tar_members.txt"
)

inspect_member () {

    member="$1"

    echo
    echo "------------------------------------------------------------"
    echo "MEMBER: ${member}"
    echo "------------------------------------------------------------"

    safe_name=$(
      basename "${member}" \
      | tr '/' '_'
    )

    archived="${TMP}/${safe_name}"
    decoded="${TMP}/${safe_name%.gz}"

    tar -xOf "${RAW_TAR}" "${member}" \
      > "${archived}"

    if echo "${member}" | grep -q '\.gz$'; then

        gzip -t "${archived}"
        echo "[PASS] member gzip integrity"

        gzip -dc "${archived}" \
          > "${decoded}"

    else

        cp "${archived}" "${decoded}"

    fi

    echo "Decoded lines:"
    wc -l "${decoded}"

    echo
    echo "Column count for first 5 lines:"

    awk -F '\t' '
      NR <= 5 {
        print "line", NR, "NF=" NF
      }
    ' "${decoded}"

    echo
    echo "First 5 lines:"

    head -n 5 "${decoded}"
}


inspect_member "${FIRST_MEMBER}"
inspect_member "${MIDDLE_MEMBER}"
inspect_member "${LAST_MEMBER}"

rm -rf "${TMP}"

echo


# ============================================================
# 15. Series matrix quick schema
# ============================================================

echo "============================================================"
echo "SERIES MATRIX QUICK INSPECTION"
echo "============================================================"

echo "Uncompressed line count:"
gzip -dc "${MATRIX}" \
  | wc -l

echo
echo "First 25 lines:"

gzip -dc "${MATRIX}" \
  | head -n 25 || true

echo


# ============================================================
# 16. Final summary
# ============================================================

echo "============================================================"
echo "STAGE 0 SUMMARY"
echo "============================================================"

echo "SOFT samples:      ${N_SOFT_SAMPLES}"
echo "RAW tar members:   ${N_RAW_MEMBERS}"
echo "RAW unique GSMs:   ${N_GSM_UNIQ}"

echo

if \
  [ "${N_SOFT_SAMPLES}" -eq 245 ] && \
  [ "${N_RAW_MEMBERS}" -eq 245 ] && \
  [ "${N_GSM_UNIQ}" -eq 245 ]
then
    echo "[PASS] Public sample-count consistency"
else
    echo "[WARN] Sample-count inconsistency detected"
fi

echo
echo "IMPORTANT:"
echo "- No authoritative expression matrix has been built yet."
echo "- No normalization has been performed."
echo "- No gene annotation has been altered."
echo "- No samples have been excluded."
echo "- No GSE57945/GSE93624 overlap has been assumed."
echo

echo "FINAL STATUS:"
echo "PASS_GSE93624_STAGE0_PUBLIC_PROVENANCE_INVENTORY_COMPLETED"
echo "============================================================"
