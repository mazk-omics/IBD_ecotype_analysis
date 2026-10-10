#!/usr/bin/env bash
set -u
root=/home/mazekai/IBD_EcoTyper/06_validation/GSE214695_scRNA_recovery_20261009
conda=/home/mazekai/miniconda3/bin/conda
if [ -f "$root/download_only_until_tomorrow.txt" ]; then
    date -Is > "$root/download_finished.txt"
    printf 'Download only; recovery deferred at user request.\n'
    exit 0
fi
"$conda" run --no-capture-output -n ibd_refaudit python -u /home/mazekai/IBD_EcoTyper/scripts/python/prepare_gse214695_recovery_pilot.py > "$root/preparation.log" 2>&1
status=$?
printf '%s\n' "$status" > "$root/preparation_exit_code.txt"
if [ "$status" -ne 0 ]; then exit "$status"; fi
bash /home/mazekai/IBD_EcoTyper/scripts/bash/run_gse214695_recovery_pilot_20261009.sh
status=$?
printf '%s\n' "$status" > "$root/recovery_runner_exit_code.txt"
if [ "$status" -ne 0 ]; then exit "$status"; fi
"$conda" run --no-capture-output -n ibd_refaudit python -u /home/mazekai/IBD_EcoTyper/scripts/python/summarize_gse214695_recovery_pilot.py > "$root/audit.log" 2>&1
status=$?
printf '%s\n' "$status" > "$root/audit_exit_code.txt"
exit "$status"
