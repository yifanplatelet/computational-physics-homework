#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <digits-for-tuning> [threads]" >&2
  exit 1
fi

DIGITS="$1"
THREADS="${2:-$(nproc)}"
export OMP_NUM_THREADS="$THREADS"

BIN="./pi_bs_unified_exact"
if [[ ! -x "$BIN" ]]; then
  echo "Build $BIN first." >&2
  exit 1
fi

best_time=""
best_cmd=""

run_case() {
  local cmd="$1"
  local out t
  out=$(eval "$cmd")
  t=$(printf '%s\n' "$out" | awk -F= '/^total_time=/{print $2}' | awk '{print $1}')
  printf '%s => %s s\n' "$cmd" "$t"
  if [[ -z "$best_time" ]] || awk "BEGIN{exit !($t < $best_time)}"; then
    best_time="$t"
    best_cmd="$cmd"
  fi
}

for depth in 5 6 7; do
  for leaf in 128 192 256 320; do
    for minsize in 2048 4096 8192; do
      run_case "$BIN $DIGITS $depth $leaf $minsize 0 1 0 0 0"
    done
  done
done

# optional light top-skew around the best-known family
for leaf in 192 256; do
  for minsize in 4096 8192; do
    run_case "$BIN $DIGITS 6 $leaf $minsize 2 3 1 131072 0"
    run_case "$BIN $DIGITS 6 $leaf $minsize 2 3 2 131072 0"
  done
done

echo
echo "BEST: $best_cmd"
echo "TIME: $best_time s"
