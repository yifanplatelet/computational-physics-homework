#!/usr/bin/env bash
set -euo pipefail
DIGITS="${1:-10000000}"
DEPTH="${2:-5}"
LEAF="${3:-128}"
TASKMIN="${4:-4096}"
PRINTPI="${5:-1}"

export OMP_NUM_THREADS="${OMP_NUM_THREADS:-8}"
export OMP_PROC_BIND="${OMP_PROC_BIND:-close}"
export OMP_PLACES="${OMP_PLACES:-cores}"

exec ./pi_ycruncher_style "$DIGITS" "$DEPTH" "$LEAF" "$TASKMIN" "$PRINTPI"
