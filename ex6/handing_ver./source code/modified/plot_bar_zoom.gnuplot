set terminal pngcairo size 1200,800 enhanced font 'Arial,14'
set output 'scheme_comparison_zoom.png'
set title 'Zoom: 2017-2024 and fit'
set xlabel 'Start year / fitting scheme'
set ylabel 'Predicted value for 2026'
set grid ytics
set boxwidth 0.72
set style fill solid 0.85 border -1
set xrange [6.5:15.5]
set yrange [-1133.000000:2000.000000]
set xtics ('2017' 7, '2018' 8, '2019' 9, '2020' 10, '2021' 11, '2022' 12, '2023' 13, '2024' 14, 'fit' 15) rotate by -45
plot \
'scheme_compare_zoom_pos.dat' using 1:2 with boxes fc rgb 'red' title 'Positive prediction', \
'scheme_compare_zoom_neg.dat' using 1:2 with boxes fc rgb 'blue' title 'Negative prediction'
