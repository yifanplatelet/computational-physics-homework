#!/usr/bin/env bash
set -euo pipefail
DIGITS="${1:-10000000}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-8}"
export OMP_PROC_BIND="${OMP_PROC_BIND:-close}"
export OMP_PLACES="${OMP_PLACES:-cores}"

best_time=""
best_cmd=""

for task_depth in 4 5 6; do
  for leaf_size in 96 128 160 192 256; do
    for task_min_size in 2048 4096 8192; do
      for merge_depth in 0 1 2 3; do
        for merge_min_size in 32768 65536 131072; do
          cmd="./pi_bs_yc_inspired $DIGITS $task_depth $leaf_size $task_min_size $merge_depth $merge_min_size"
          echo "RUN: $cmd"
          out="$($cmd)"
          echo "$out"
          t="$(printf '%s\n' "$out" | awk -F= '/^total_time=/{print $2}')"
          if [ -z "$best_time" ] || awk "BEGIN{exit !($t < $best_time)}"; then
            best_time="$t"
            best_cmd="$cmd"
          fi
          echo "----"
        done
      done
    done
  done
done

echo "BEST: $best_cmd"
echo "TIME: $best_time s"
