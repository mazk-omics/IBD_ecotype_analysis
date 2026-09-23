#!/usr/bin/env bash

set -u
set -o pipefail

# ============================================================
# GSE93624
# Per-query NCBI Datasets forensic
#
# IMPORTANT:
#   One manifest row -> one independent NCBI request.
#
# This preserves query provenance and does NOT assume that
# batch NCBI output order matches input order.
#
# Resume-safe:
#   already successful QueryIDs are skipped.
# ============================================================

BASE="/home/mazekai/IBD_EcoTyper/06_validation/01_bulk_progression/GSE93624"

QUERY_DIR="${BASE}/prepared/02_gene_harmonization/ncbi_external_queries"

MANIFEST="${QUERY_DIR}/GSE93624_NCBI_external_query_manifest.tsv"

RAW_DIR="${QUERY_DIR}/raw_jsonl"
LOG_DIR="${QUERY_DIR}/logs"

STATUS="${QUERY_DIR}/GSE93624_NCBI_external_query_status.tsv"

mkdir -p \
  "${RAW_DIR}" \
  "${LOG_DIR}"


echo "============================================================"
echo "GSE93624 external NCBI Datasets forensic"
echo "============================================================"
echo


# ============================================================
# 1. CLI audit
# ============================================================

if ! command -v datasets >/dev/null 2>&1; then
    echo "[ERROR] datasets CLI not found"
    exit 1
fi

echo "datasets executable:"
command -v datasets

echo
echo "datasets version:"
datasets version 2>/dev/null || datasets --version

echo


# ============================================================
# 2. Initialize status
# ============================================================

if [ ! -f "${STATUS}" ]; then
    printf \
      "QueryID\tQueryType\tQuery\tSourceSymbol\tReason\tExitCode\tResultLines\tStatus\n" \
      > "${STATUS}"
fi


# ============================================================
# 3. Resume helper
# ============================================================

already_passed () {

    query_id="$1"

    awk -F '\t' \
      -v q="${query_id}" \
      '
      NR > 1 && $1 == q && $8 == "PASS" {
          found = 1
      }
      END {
          exit(found ? 0 : 1)
      }
      ' \
      "${STATUS}"
}


# ============================================================
# 4. Run one query with retry
# ============================================================

run_query () {

    query_id="$1"
    query_type="$2"
    query="$3"
    source_symbol="$4"
    reason="$5"

    outfile="${RAW_DIR}/${query_id}.jsonl"
    logfile="${LOG_DIR}/${query_id}.stderr.log"

    if already_passed "${query_id}" && [ -f "${outfile}" ]; then
        echo "[SKIP PASS] ${query_id} ${query_type} ${query}"
        return 0
    fi

    echo "------------------------------------------------------------"
    echo "${query_id}"
    echo "  type   : ${query_type}"
    echo "  query  : ${query}"
    echo "  source : ${source_symbol}"
    echo "  reason : ${reason}"

    rm -f \
      "${outfile}.tmp" \
      "${logfile}"

    success=0
    exit_code=1

    for attempt in 1 2 3; do

        echo "  attempt: ${attempt}"

        if [ "${query_type}" = "GENE_ID" ]; then

            datasets summary gene gene-id \
              "${query}" \
              --limit all \
              --as-json-lines \
              > "${outfile}.tmp" \
              2> "${logfile}"

            exit_code=$?

        elif [ "${query_type}" = "SYMBOL" ]; then

            datasets summary gene symbol \
              "${query}" \
              --taxon human \
              --limit all \
              --as-json-lines \
              > "${outfile}.tmp" \
              2> "${logfile}"

            exit_code=$?

        else

            echo "[ERROR] Unknown QueryType: ${query_type}"
            exit_code=99
        fi

        if [ "${exit_code}" -eq 0 ]; then
            success=1
            break
        fi

        sleep $((attempt * 2))
    done


    # --------------------------------------------------------
    # Preserve response even if it contains zero records.
    # --------------------------------------------------------

    if [ "${success}" -eq 1 ]; then

        mv \
          "${outfile}.tmp" \
          "${outfile}"

        result_lines=$(
          grep -cve '^[[:space:]]*$' \
            "${outfile}" \
          || true
        )

        status="PASS"

    else

        rm -f "${outfile}.tmp"

        result_lines=0
        status="FAILED"
    fi


    printf \
      "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
      "${query_id}" \
      "${query_type}" \
      "${query}" \
      "${source_symbol}" \
      "${reason}" \
      "${exit_code}" \
      "${result_lines}" \
      "${status}" \
      >> "${STATUS}"

    echo "  exit    : ${exit_code}"
    echo "  records : ${result_lines}"
    echo "  status  : ${status}"

    # modest delay to avoid unnecessarily hammering NCBI
    sleep 0.2
}


# ============================================================
# 5. Execute manifest
#
# Skip header.
# IFS tab preserves symbol punctuation safely.
# ============================================================

tail -n +2 "${MANIFEST}" |
while IFS=$'\t' read -r \
    query_id \
    query_type \
    query \
    source_symbol \
    reason \
    expected_geneid \
    previous_symbol

do

    run_query \
      "${query_id}" \
      "${query_type}" \
      "${query}" \
      "${source_symbol}" \
      "${reason}"

done


# ============================================================
# 6. Final audit
# ============================================================

echo
echo "============================================================"
echo "QUERY STATUS"
echo "============================================================"

PASS_N=$(
  awk -F '\t' '
    NR > 1 && $8 == "PASS" {n++}
    END {print n+0}
  ' "${STATUS}"
)

FAILED_N=$(
  awk -F '\t' '
    NR > 1 && $8 == "FAILED" {n++}
    END {print n+0}
  ' "${STATUS}"
)

EXPECTED_N=$(
  awk '
    END {print NR-1}
  ' "${MANIFEST}"
)

echo "Expected manifest queries : ${EXPECTED_N}"
echo "Successful query records  : ${PASS_N}"
echo "Failed query records      : ${FAILED_N}"

# Count unique successfully completed QueryIDs.
UNIQUE_PASS_N=$(
  awk -F '\t' '
    NR > 1 && $8 == "PASS" {
      seen[$1]=1
    }
    END {
      for (x in seen) n++
      print n+0
    }
  ' "${STATUS}"
)

echo "Unique successful QueryIDs: ${UNIQUE_PASS_N}"


if [ "${UNIQUE_PASS_N}" -eq "${EXPECTED_N}" ]; then

    echo
    echo "FINAL STATUS:"
    echo "PASS_ALL_GSE93624_NCBI_EXTERNAL_QUERIES_COMPLETED"

else

    echo
    echo "FINAL STATUS:"
    echo "INCOMPLETE_GSE93624_NCBI_EXTERNAL_QUERIES"
    echo
    echo "Re-run the same script to resume."

fi

echo "============================================================"
