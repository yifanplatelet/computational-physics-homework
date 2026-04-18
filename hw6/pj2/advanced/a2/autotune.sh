#!/usr/bin/env bash
set -euo pipefail

DIGITS=${1:-10000000}
export OMP_PROC_BIND=${OMP_PROC_BIND:-close}
export OMP_PLACES=${OMP_PLACES:-cores}

best=""
best_time=""

for depth in 8 10 12; do
  for leaf in 128 192 256; do
    for tmin in 2048 4096 8192; do
      for mdepth in 1 2 3; do
        for mmin in 32768 65536 131072; do
          cmd=(./pi_bs_optimized_yc_learned "$DIGITS" "$depth" "$leaf" "$tmin" "$mdepth" "$mmin")
          echo "RUN: ${cmd[*]}"
          out="$(${cmd[@]})"
          time=$(printf '%s
' "$out" | awk -F= '/^total_time=/{print $2}')
          echo "TIME: $time"
          if [[ -z "$best_time" ]] || awk "BEGIN{exit !($time < $best_time)}"; then
            best_time="$time"
            best="${cmd[*]}"
          fi
        done
      done
    done
  done
done

echo
echo "BEST: $best"
echo "TIME: $best_time s"
