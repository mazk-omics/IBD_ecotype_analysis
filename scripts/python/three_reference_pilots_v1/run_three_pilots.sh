#!/usr/bin/env bash
set -euo pipefail
THREE_PILOT_ROOT=/home/mazekai/IBD_EcoTyper
THREE_PILOT_SCRIPTS="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
THREE_PILOT_RUN="${THREE_PILOT_ROOT}/03_reference/web_reference_three_pilots/output/run_$(date +%Y%m%d_%H%M%S)"
THREE_PILOT_CONDA=/home/mazekai/miniconda3/bin/conda
mkdir -p "${THREE_PILOT_RUN}"
printf '%s\n' "${THREE_PILOT_RUN}" > "${THREE_PILOT_SCRIPTS}/last_run_path.txt"
"${THREE_PILOT_CONDA}" run --no-capture-output -n ibd_refaudit python "${THREE_PILOT_SCRIPTS}/verify_three_reference_web_pilots.py" > "${THREE_PILOT_RUN}/selftest.log" 2>&1
tail -n 1 "${THREE_PILOT_RUN}/selftest.log"
for THREE_PILOT_MODE in gse_joint gse_ileum; do
 "${THREE_PILOT_CONDA}" run --no-capture-output -n ibd_refaudit \
 python -u "${THREE_PILOT_SCRIPTS}/export_three_reference_web_pilots.py" \
 --mode "${THREE_PILOT_MODE}" \
 --input /home/mazekai/enteric_glia/data_scRNA/GSE282122/TAURUS_raw_counts_annotated_final.h5ad \
 --candidate "${THREE_PILOT_ROOT}/03_reference/GSE282122/audit_v2/output/04_identity_coverage_v1/mapped_candidate_cells_v1.parquet" \
 --donor-split "${THREE_PILOT_ROOT}/03_reference/GSE282122/web_pilot/output/colon_pilot_20261006_174155/donor_split.tsv" \
 --outdir "${THREE_PILOT_RUN}/${THREE_PILOT_MODE}" --cells-per-type 50 --seed 20261006 \
 2>&1 | tee "${THREE_PILOT_RUN}/${THREE_PILOT_MODE}.log"
done
"${THREE_PILOT_CONDA}" run --no-capture-output -n ibd_refaudit \
 python -u "${THREE_PILOT_SCRIPTS}/export_three_reference_web_pilots.py" \
 --mode scp259 \
 --raw "${THREE_PILOT_ROOT}/01_raw_data/scRNA/SCP259/raw/SCP259" \
 --candidate "${THREE_PILOT_ROOT}/03_reference/SCP259/audit_v1/output/01_candidate_taxonomy_coverage/scp259_candidate_cells_provisional_taxonomy.parquet" \
 --outdir "${THREE_PILOT_RUN}/scp259" --cells-per-type 50 --seed 20261006 \
 2>&1 | tee "${THREE_PILOT_RUN}/scp259.log"
printf 'ALL_THREE_COMPLETED %s\n' "${THREE_PILOT_RUN}"
