#!/usr/bin/env bash
set -euo pipefail

VALUES_CSV="${1:-values.csv}"
RELERR_CSV="${2:-relerr.csv}"
OUTDIR="${3:-plots}"

mkdir -p "$OUTDIR"

if [[ ! -f "$VALUES_CSV" ]]; then
    echo "错误: 找不到 $VALUES_CSV"
    exit 1
fi

if [[ ! -f "$RELERR_CSV" ]]; then
    echo "错误: 找不到 $RELERR_CSV"
    exit 1
fi

gnuplot <<EOF
set datafile separator ","
set datafile missing "nan"
set terminal pngcairo size 1800,1200 enhanced font "Arial,14"

# -------------------------------
# 图1: 函数值比较
# values.csv 列:
# 1:x
# 2:direct_f  3:expand_f  4:horner_f
# 5:direct_d  6:expand_d  7:horner_d
# 8:direct_q  9:expand_q 10:horner_q
# -------------------------------
set output "${OUTDIR}/values_all.png"
set multiplot layout 3,1 title "(x-1)^10 : direct vs expansion vs Horner   (zoom: x in [0.7,1.3])"

set grid
set xrange [0.7:1.3]
set xlabel "x"
set ylabel "y"
set key outside
set format y "%.2e"

set title "float"
plot "${VALUES_CSV}" using 1:2 with lines lw 2 title "direct_f", \
     "${VALUES_CSV}" using 1:3 with lines lw 2 title "expand_f", \
     "${VALUES_CSV}" using 1:4 with lines lw 2 title "horner_f"

set title "double"
plot "${VALUES_CSV}" using 1:5 with lines lw 2 title "direct_d", \
     "${VALUES_CSV}" using 1:6 with lines lw 2 title "expand_d", \
     "${VALUES_CSV}" using 1:7 with lines lw 2 title "horner_d"

set title "quad (__float128)"
plot "${VALUES_CSV}" using 1:8 with lines lw 2 title "direct_q", \
     "${VALUES_CSV}" using 1:9 with lines lw 2 title "expand_q", \
     "${VALUES_CSV}" using 1:10 with lines lw 2 title "horner_q"

unset multiplot


# -------------------------------
# 图2: 相对误差比较
# relerr.csv 列:
# 1:x
# 2:rel-direct_f  3:rel-expand_f  4:rel-horner_f
# 5:rel-direct_d  6:rel-expand_d  7:rel-horner_d
# 8:rel-direct_q  9:rel-expand_q 10:rel-horner_q
# -------------------------------
set output "${OUTDIR}/relerr_all.png"
set multiplot layout 3,1 title "Relative error distribution (log scale)"

set grid
set xrange [0.7:1.3]
set xlabel "x"
set ylabel "relative error"
set key outside
set logscale y
set format y "10^{%L}"

set title "float"
plot "${RELERR_CSV}" using 1:2 with lines lw 2 title "rel-direct_f", \
     "${RELERR_CSV}" using 1:3 with lines lw 2 title "rel-expand_f", \
     "${RELERR_CSV}" using 1:4 with lines lw 2 title "rel-horner_f"

set title "double"
plot "${RELERR_CSV}" using 1:5 with lines lw 2 title "rel-direct_d", \
     "${RELERR_CSV}" using 1:6 with lines lw 2 title "rel-expand_d", \
     "${RELERR_CSV}" using 1:7 with lines lw 2 title "rel-horner_d"

set title "quad (float128)"
plot "${RELERR_CSV}" using 1:8 with lines lw 2 title "rel-direct_q", \
     "${RELERR_CSV}" using 1:9 with lines lw 2 title "rel-expand_q", \
     "${RELERR_CSV}" using 1:10 with lines lw 2 title "rel-horner_q"

unset logscale y
unset multiplot
EOF

echo "绘图完成:"
echo "  ${OUTDIR}/values_all.png"
echo "  ${OUTDIR}/relerr_all.png"
