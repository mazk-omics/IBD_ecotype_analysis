#!/usr/bin/env bash
set -uo pipefail

project_root=/home/mazekai/IBD_EcoTyper
run_root="$project_root/06_validation/MSCCR_SignatureD_Cell2021_STAR_downstream_plan_20261119"
job_root="$run_root/output/31_bulk_recovery/GSE57945_current85_20261124"
ecotyper_root="$project_root/software/ecotyper"
input_root="$project_root/06_validation/01_bulk_progression/GSE57945/prepared/06_ecotyper_recovery_input"
expression_source="$input_root/GSE57945_EcoTyper_recovery_expression.txt"
expression_file="$job_root/input/GSE57945_current85_bulk_v2.txt"
annotation_file="$job_root/input/GSE57945_current85_annotation_v2.tsv"
discovery_name=MSCCR_SIGNATURED_IBD2025_TPM_REGION3_20261101
threads=20

mkdir -p "$job_root"
if [ ! -L "$expression_file" ]; then
  mkdir -p "$(dirname "$expression_file")"
  ln -s "$expression_source" "$expression_file"
fi
recovery_name=$(basename "$(readlink -f "$expression_file")")
recovery_name=${recovery_name%.txt}
recovery_name=${recovery_name%.tsv}
if [ ! -s "$annotation_file" ]; then
  printf 'status=FAILED\nend=%s\ndetail=Prepared multi-column annotation is missing\n' "$(date --iso-8601=seconds)" > "$job_root/status.txt"
  exit 2
fi
printf 'status=RUNNING\nstart=%s\n' "$(date --iso-8601=seconds)" > "$job_root/status.txt"
{
  printf 'host=%s\n' "$(hostname)"
  printf 'discovery=%s\n' "$discovery_name"
  printf 'recovery=%s\n' "$recovery_name"
  printf 'threads=%s\n' "$threads"
  sha256sum \
    "$expression_file" \
    "$annotation_file" \
    "$ecotyper_root/EcoTyper/$discovery_name/Custom/Analysis/rank_selection/rank_data.txt" \
    "$ecotyper_root/EcoTyper_recovery_bulk.R" \
    "$ecotyper_root/pipeline/state_recovery_bulk.R" \
    "$ecotyper_root/pipeline/ecotypes_recovery.R"
} > "$job_root/run_identity.txt"

cd "$ecotyper_root"
/home/mazekai/miniconda3/bin/conda run --no-capture-output -n ecotyper \
  Rscript EcoTyper_recovery_bulk.R \
  -d "$discovery_name" \
  -m "$expression_file" \
  -a "$annotation_file" \
  -c Disease,DiseaseSubtype,Tissue,Histopathology \
  -t "$threads" \
  -o "$job_root/RecoveryOutput" \
  > "$job_root/recovery.console.txt" 2>&1
exit_code=$?
printf '%s\n' "$exit_code" > "$job_root/recovery.exit.txt"

result_root="$job_root/RecoveryOutput/$recovery_name"
state_files=$(find "$result_root" -mindepth 2 -maxdepth 2 -name state_abundances.txt -type f 2>/dev/null | wc -l)
assignment_files=$(find "$result_root" -mindepth 2 -maxdepth 2 -name state_assignment.txt -type f 2>/dev/null | wc -l)
ecotype_abundance=0
ecotype_assignment=0
if [ -s "$result_root/Ecotypes/ecotype_abundance.txt" ]; then
  ecotype_abundance=1
fi
if [ -s "$result_root/Ecotypes/ecotype_assignment.txt" ]; then
  ecotype_assignment=1
fi

{
  printf 'exit_code\t%s\n' "$exit_code"
  printf 'state_abundance_files\t%s\n' "$state_files"
  printf 'state_assignment_files\t%s\n' "$assignment_files"
  printf 'ecotype_abundance_present\t%s\n' "$ecotype_abundance"
  printf 'ecotype_assignment_present\t%s\n' "$ecotype_assignment"
} > "$job_root/output_audit.tsv"

if [ "$exit_code" -eq 0 ] && [ "$state_files" -eq 20 ] && [ "$assignment_files" -eq 20 ] && [ "$ecotype_abundance" -eq 1 ] && [ "$ecotype_assignment" -eq 1 ]; then
  printf 'status=ACCEPTED\nend=%s\ndetail=Official bulk recovery completed with 20 cell types and ecotype outputs\n' "$(date --iso-8601=seconds)" > "$job_root/status.txt"
else
  printf 'status=FAILED\nend=%s\ndetail=Recovery exit or output audit failed; inspect logs without automatic rerun\n' "$(date --iso-8601=seconds)" > "$job_root/status.txt"
fi

exit "$exit_code"
