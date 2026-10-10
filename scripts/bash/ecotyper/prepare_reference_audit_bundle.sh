#!/usr/bin/env bash
# Package listed small files only; never run BayesPrism or modify source data.
set -euo pipefail
project_root=/home/mazekai/IBD_EcoTyper
files=(
  "03_reference/combined_reference/output/02_reference_availability/final_reference_ontology.tsv"
  "03_reference/combined_reference/output/05_final_cell_selection/selected_cell_manifest.tsv.gz"
  "03_reference/combined_reference/output/07_unified_reference/unified_reference_identity_counts.tsv"
  "03_reference/combined_reference/output/09_bayesprism_input_harmonization/unified_reference_to_MSCCR_Ensembl75_gene_map.tsv"
  "03_reference/combined_reference/output/09_bayesprism_input_harmonization/MSCCR_reference_gene_overlap_audit.json"
  "03_reference/combined_reference/output/13_bayesprism_broad_input/01_build_broad14_reference.log"
  "03_reference/combined_reference/output/17_bayesprism_epithelial_input/stage2_input_manifest.tsv"
  "03_reference/combined_reference/output/17_bayesprism_epithelial_input/stage2_reference_cell_counts.tsv"
  "03_reference/combined_reference/output/17_bayesprism_epithelial_input/stage2_reference_audit.tsv"
  "03_reference/combined_reference/output/18_bayesprism_epithelial_preprocessing/stage2_preprocessing_manifest.tsv"
  "03_reference/combined_reference/output/18_bayesprism_epithelial_preprocessing/stage2_reference_mixture_common_genes.txt"
  "03_reference/combined_reference/output/19_bayesprism_epithelial_prism/stage2_prism_manifest.tsv"
  "03_reference/combined_reference/output/20_bayesprism_epithelial_full_run/stage2_full_run_manifest.tsv"
  "03_reference/combined_reference/output/20_bayesprism_epithelial_full_run/bayesprism_epithelial7_reference_update.rds"
  "03_reference/combined_reference/output/16_bayesprism_broad14_neutrophil_allocation_audit/primary_neutrophil_marker_aggregate_allocation.tsv"
  "03_reference/combined_reference/output/16_bayesprism_broad14_neutrophil_allocation_audit/primary_neutrophil_gene_allocation_summary.tsv"
  "03_reference/combined_reference/output/16_bayesprism_broad14_neutrophil_allocation_audit/04_neutrophil_posterior_allocation_audit.log"
  "scripts/R/bayesprism_broad/01_build_broad14_reference.R"
  "scripts/R/bayesprism_epithelial/02_cleanup_stage2_epithelial_reference.R"
)
for file in "${files[@]}"; do
  if [[ ! -s "$project_root/$file" ]]; then
    printf 'Missing or empty: %s\n' "$project_root/$file" >&2
    exit 1
  fi
done
bundle_dir=$(mktemp -d "$project_root/03_reference/combined_reference/qc/ecotyper_reference_audit_bundle_XXXXXXXX")
printf '%s\n' "${files[@]}" > "$bundle_dir/files.txt"
tar -czf "$bundle_dir/bayesprism_reference_audit_small.tar.gz" -C "$project_root" --files-from="$bundle_dir/files.txt"
printf 'Download this archive: %s\n' "$bundle_dir/bayesprism_reference_audit_small.tar.gz"
ls -lh "$bundle_dir/bayesprism_reference_audit_small.tar.gz"
