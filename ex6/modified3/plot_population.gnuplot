set terminal pngcairo size 1500,900 enhanced font 'Arial,14'
set output 'interpolation_extrapolation.png'
set title 'Interpolation Curves and 2026 Extrapolated Points'
set xlabel 'Year'
set ylabel 'Value'
set grid
set key outside
set xrange [2010:2026.8]
set style data lines
plot \
'curve_2010_2025.dat' using 1:2 with lines lw 1 title '2010-2025 fit', \
'curve_2011_2025.dat' using 1:2 with lines lw 1 title '2011-2025 fit', \
'curve_2012_2025.dat' using 1:2 with lines lw 1 title '2012-2025 fit', \
'curve_2013_2025.dat' using 1:2 with lines lw 1 title '2013-2025 fit', \
'curve_2014_2025.dat' using 1:2 with lines lw 1 title '2014-2025 fit', \
'curve_2015_2025.dat' using 1:2 with lines lw 1 title '2015-2025 fit', \
'curve_2016_2025.dat' using 1:2 with lines lw 1 title '2016-2025 fit', \
'curve_2017_2025.dat' using 1:2 with lines lw 1 title '2017-2025 fit', \
'curve_2018_2025.dat' using 1:2 with lines lw 1 title '2018-2025 fit', \
'curve_2019_2025.dat' using 1:2 with lines lw 1 title '2019-2025 fit', \
'curve_2020_2025.dat' using 1:2 with lines lw 1 title '2020-2025 fit', \
'curve_2021_2025.dat' using 1:2 with lines lw 1 title '2021-2025 fit', \
'curve_2022_2025.dat' using 1:2 with lines lw 1 title '2022-2025 fit', \
'curve_2023_2025.dat' using 1:2 with lines lw 1 title '2023-2025 fit', \
'curve_2024_2025.dat' using 1:2 with lines lw 1 title '2024-2025 fit', \
'actual_data.dat' using 1:2 with points pt 7 ps 1.6 title 'Actual data', \
'all_predictions.dat' using 1:2 with points pt 5 ps 1.6 title 'All 2026 predictions', \
'all_predictions.dat' using 1:2:3 with labels offset 1,0.5 notitle, \
'best_predictions.dat' using 1:2 with points pt 9 ps 2.2 title 'Less extreme predictions', \
'best_predictions.dat' using 1:2:3 with labels offset 1,1.2 notitle
