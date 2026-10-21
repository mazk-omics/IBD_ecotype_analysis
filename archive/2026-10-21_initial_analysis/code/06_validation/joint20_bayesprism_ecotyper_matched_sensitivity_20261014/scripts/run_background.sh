#!/usr/bin/env bash
# 在独立tmux内运行正式推断和候选表达QC；保留退出码，不自动重启。
set -uo pipefail
out=/home/mazekai/IBD_EcoTyper/06_validation/joint20_bayesprism_ecotyper_matched_sensitivity_20261014
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
date -Is > "$out/logs/started_at.txt"
echo "$$" > "$out/logs/launcher.pid"
trap 'rc=$?; printf "%s\n" "$rc" > "$out/logs/exit_status.txt"; date -Is > "$out/logs/ended_at.txt"' EXIT
test -f "$out/00_audit/input_audit.json" || exit 10
test ! -e "$out/02_bayesprism/BayesPrism_result.rds" || exit 11
command_line="/home/mazekai/miniconda3/envs/bayesprism/bin/Rscript --vanilla $out/scripts/run_full_joint20.R $out"
printf '%s\n' "$command_line" > "$out/logs/formal_command.txt"
/usr/bin/time -v -o "$out/logs/formal_resource_usage.txt" /home/mazekai/miniconda3/envs/bayesprism/bin/Rscript --vanilla "$out/scripts/run_full_joint20.R" "$out" > "$out/logs/formal.stdout.log" 2> "$out/logs/formal.stderr.log" &
timing_pid=$!
echo "$timing_pid" > "$out/logs/time.pid"
wait "$timing_pid"
formal_rc=$?
echo "$formal_rc" > "$out/logs/formal_exit_status.txt"
test "$formal_rc" -eq 0 || exit "$formal_rc"
/home/mazekai/miniconda3/envs/bayesprism/bin/Rscript --vanilla "$out/scripts/prepare_glike_qc.R" "$out" > "$out/logs/postprocessing.stdout.log" 2> "$out/logs/postprocessing.stderr.log" &
post_pid=$!
echo "$post_pid" > "$out/logs/postprocessing.pid"
wait "$post_pid"
post_rc=$?
echo "$post_rc" > "$out/logs/postprocessing_exit_status.txt"
test "$post_rc" -eq 0 || exit "$post_rc"
date -Is > "$out/logs/completed_at.txt"
