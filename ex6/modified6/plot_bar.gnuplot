set terminal pngcairo size 1400,800 enhanced font 'Arial,14'
set output 'scheme_comparison.png'
set title 'Comparison of 2026 Predictions by Different Schemes'
set xlabel 'Start year / fitting scheme'
set ylabel 'Predicted value for 2026'
set grid ytics
set yrange [-1100:2000]
set style data histograms
set style fill solid 0.8 border -1
set boxwidth 0.7
set xtics rotate by -45
plot 'scheme_compare.dat' using 2:xtic(1) title 'Prediction'
