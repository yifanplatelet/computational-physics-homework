#!/usr/bin/env bash
set -euo pipefail

DIGITS=${1:-10000000}
export OMP_PROC_BIND=${OMP_PROC_BIND:-close}
export OMP_PLACES=${OMP_PLACES:-cores}

best=""
best_time=""

run_case() {
  local -a cmd=("$@")
  local out time
  echo "RUN: ${cmd[*]}"
  out="$(${cmd[@]})"
  time=$(printf '%s\n' "$out" | awk -F= '/^total_time=/{print $2}')
  echo "TIME: $time"
  if [[ -z "$best_time" ]] || awk "BEGIN{exit !($time < $best_time)}"; then
    best_time="$time"
    best="${cmd[*]}"
  fi
}

for depth in 6 7 8; do
  for leaf in 192 256 320; do
    for tmin in 4096 8192 16384; do
      for mdepth in 0 1 2; do
        for mmin in 65536 131072; do
          run_case ./pi_bs_fixedpoint_yc_adaptive "$DIGITS" "$depth" "$leaf" "$tmin" "$mdepth" "$mmin" 32 0 0 1 0 0
        done
      done
    done
  done
done

for depth in 6 7; do
  for leaf in 256 320; do
    for tmin in 4096 8192; do
      for mdepth in 0 1; do
        for mmin in 65536 131072; do
          run_case ./pi_bs_fixedpoint_yc_adaptive "$DIGITS" "$depth" "$leaf" "$tmin" "$mdepth" "$mmin" 32 0 2 3 1 131072
          run_case ./pi_bs_fixedpoint_yc_adaptive "$DIGITS" "$depth" "$leaf" "$tmin" "$mdepth" "$mmin" 32 0 2 3 2 131072
        done
      done
    done
  done
done

echo
echo "BEST: $best"
echo "TIME: $best_time s"
