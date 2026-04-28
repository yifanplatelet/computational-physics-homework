#!/usr/bin/env bash
set -euo pipefail

DIGITS=${1:-100000000}
DEPTH=${2:-12}
LEAF=${3:-256}
TASK_MIN=${4:-4096}
MERGE_DEPTH=${5:-2}
MERGE_MIN=${6:-65536}
GUARD=${7:-32}

export OMP_PROC_BIND=${OMP_PROC_BIND:-close}
export OMP_PLACES=${OMP_PLACES:-cores}

echo "== Baseline =="
./pi_bs_baseline_from_user "$DIGITS"

echo
echo "== Fixed-Point Finish =="
./pi_bs_fixedpoint_yc_adaptive "$DIGITS" "$DEPTH" "$LEAF" "$TASK_MIN" "$MERGE_DEPTH" "$MERGE_MIN" "$GUARD"
