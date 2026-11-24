#!/usr/bin/env bash
set -euo pipefail

RUN=/home/mazekai/IBD_EcoTyper/06_validation/MSCCR_SignatureD_Cell2021_STAR_downstream_plan_20261119
LOG="$RUN/output/G7_G7A_marker_Z_20261123.console.txt"
STATUS="$RUN/output/G7_G7A_marker_Z_20261123.status.txt"
: > "$LOG"
: > "$STATUS"
exec > >(tee -a "$LOG") 2>&1

finish() {
  rc=$?
  if [ "$rc" -eq 0 ]; then
    printf 'STATUS=PASS_G7_G7A_MARKER_Z\nFINISHED_AT=%s\n' "$(date -Is)" > "$STATUS"
  else
    printf 'STATUS=FAILED_G7_G7A_MARKER_Z\nEXIT_CODE=%s\nFINISHED_AT=%s\n' "$rc" "$(date -Is)" > "$STATUS"
  fi
  exit "$rc"
}
trap finish EXIT

printf 'STATUS=RUNNING_G7_G7A_MARKER_Z\nSTARTED_AT=%s\n' "$(date -Is)" > "$STATUS"
for script in 22_run_known_phenotype_enrichment_exact.R 23_run_known_phenotype_enrichment_patient_aware.R 24_run_cell_state_representation.R 25_plot_marker_validation_and_z.R; do
  printf '\n===== START %s %s =====\n' "$script" "$(date -Is)"
  /usr/bin/Rscript "$RUN/scripts/$script" "$RUN"
  printf '===== PASS %s %s =====\n' "$script" "$(date -Is)"
done
