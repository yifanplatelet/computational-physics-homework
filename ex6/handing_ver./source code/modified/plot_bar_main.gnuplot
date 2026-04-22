set terminal pngcairo size 1600,950 enhanced font 'Arial,14'
set output 'scheme_comparison_main.png'
set title 'Comparison of 2026 Predictions by Different Schemes'
set xlabel 'Start year / fitting scheme'
set ylabel 'Signed log scale of prediction value'
set grid ytics
set boxwidth 0.72
set style fill solid 0.85 border -1
set xtics rotate by -45
set xrange [-1:16]
set yrange [-6.316305:6.316305]
set xtics ('2010' 0, '2011' 1, '2012' 2, '2013' 3, '2014' 4, '2015' 5, '2016' 6, '2017' 7, '2018' 8, '2019' 9, '2020' 10, '2021' 11, '2022' 12, '2023' 13, '2024' 14, 'fit' 15)
set ytics ('-10^7' -7, '-10^6' -6, '-10^5' -5, '-10^4' -4, '-10^3' -3, '-10^2' -2, '-10^1' -1, '0' 0, '10^1' 1, '10^2' 2, '10^3' 3, '10^4' 4, '10^5' 5, '10^6' 6, '10^7' 7)
plot \
'scheme_compare_all_pos.dat' using 1:2 with boxes fc rgb 'red' title 'Positive prediction', \
'scheme_compare_all_neg.dat' using 1:2 with boxes fc rgb 'blue' title 'Negative prediction'
