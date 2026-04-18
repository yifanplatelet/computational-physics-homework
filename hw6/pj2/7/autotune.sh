#!/usr/bin/env bash
set -euo pipefail
DIGITS="${1:-10000000}"
THREADS="${OMP_NUM_THREADS:-8}"
export OMP_NUM_THREADS="$THREADS"
export OMP_PROC_BIND=close
export OMP_PLACES=cores

best_time=""
best_cmd=""

for depth in 4 5 6; do
  for leaf in 64 128 256; do
    for minsize in 2048 4096 8192; do
      cmd="./pi_bs_openmp_rewrite $DIGITS $depth $leaf $minsize"
      echo "RUN: $cmd"
      out=$(eval "$cmd")
      echo "$out"
      t=$(printf '%s\n' "$out" | awk -F= '/^total_time=/{gsub(/ s/,"",$2); print $2}')
      if [[ -z "$best_time" ]] || awk "BEGIN{exit !($t < $best_time)}"; then
        best_time="$t"
        best_cmd="$cmd"
      fi
      echo
    done
  done
done

echo "BEST: $best_cmd"
echo "TIME: $best_time s"
