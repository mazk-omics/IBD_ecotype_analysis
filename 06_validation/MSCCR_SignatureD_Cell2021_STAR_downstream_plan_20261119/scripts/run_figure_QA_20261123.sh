#!/usr/bin/env bash
set -u

RUN=/home/mazekai/IBD_EcoTyper/06_validation/MSCCR_SignatureD_Cell2021_STAR_downstream_plan_20261119
PYTHON=/home/mazekai/miniconda3/envs/ibd_refaudit/bin/python
QA="$RUN/output/17_figure_QA"
mkdir -p "$QA"
: > "$QA/figure_QA_results.tsv"
printf 'Figure\tTextAuditExit\tCollisionAuditExit\n' >> "$QA/figure_QA_results.tsv"

while IFS= read -r pdf; do
  relative=${pdf#"$RUN/output/"}
  stem=$(printf '%s' "$relative" | sed 's#[ /]#_#g; s#\.pdf$##')
  "$PYTHON" "$RUN/scripts/figure_qa/audit_pdf_text.py" "$pdf" --min-pt 5 --json > "$QA/$stem.text.json"
  text_exit=$?
  "$PYTHON" "$RUN/scripts/figure_qa/audit_figure_collisions.py" "$pdf" --json-out "$QA/$stem.collision.json" --overlay-pdf "$QA/$stem.collision-overlay.pdf" > "$QA/$stem.collision.console.txt" 2>&1
  collision_exit=$?
  printf '%s\t%s\t%s\n' "$relative" "$text_exit" "$collision_exit" >> "$QA/figure_QA_results.tsv"
done < <(find "$RUN/output/14_validation_figures" "$RUN/output/15_cell_state_representation_figures" -type f -name '*.pdf' | sort)

Rscript - "$QA" <<'RS'
library(data.table)
args <- commandArgs(trailingOnly = TRUE)
qa <- fread(file.path(args[[1L]], "figure_QA_results.tsv"))
summary <- qa[, .(
  PDF = .N,
  TextAuditPass = sum(TextAuditExit == 0L),
  TextAuditFail = sum(TextAuditExit != 0L),
  CollisionPass = sum(CollisionAuditExit == 0L),
  CollisionFail = sum(CollisionAuditExit != 0L)
)]
fwrite(summary, file.path(args[[1L]], "figure_QA_summary.tsv"), sep = "\t")
print(summary)
RS
