set terminal pngcairo size 1500,900 enhanced font 'Arial,14'
set output 'interpolation_extrapolation.png'
set title 'Interpolation vs Fitting for 2026 Prediction'
set xlabel 'Year'
set ylabel 'Value'
set grid
set key outside
set xrange [2009.5:2026.8]
set yrange [0:2000]
set style data lines
plot \
'curve_2010_2025.dat' using 1:2 with lines lw 1 title '2010-2025 interp', \
'curve_2011_2025.dat' using 1:2 with lines lw 1 title '2011-2025 interp', \
'curve_2012_2025.dat' using 1:2 with lines lw 1 title '2012-2025 interp', \
'curve_2013_2025.dat' using 1:2 with lines lw 1 title '2013-2025 interp', \
'curve_2014_2025.dat' using 1:2 with lines lw 1 title '2014-2025 interp', \
'curve_2015_2025.dat' using 1:2 with lines lw 1 title '2015-2025 interp', \
'curve_2016_2025.dat' using 1:2 with lines lw 1 title '2016-2025 interp', \
'curve_2017_2025.dat' using 1:2 with lines lw 1 title '2017-2025 interp', \
'curve_2018_2025.dat' using 1:2 with lines lw 1 title '2018-2025 interp', \
'curve_2019_2025.dat' using 1:2 with lines lw 1 title '2019-2025 interp', \
'curve_2020_2025.dat' using 1:2 with lines lw 1 title '2020-2025 interp', \
'curve_2021_2025.dat' using 1:2 with lines lw 1 title '2021-2025 interp', \
'curve_2022_2025.dat' using 1:2 with lines lw 1 title '2022-2025 interp', \
'curve_2023_2025.dat' using 1:2 with lines lw 1 title '2023-2025 interp', \
'curve_2024_2025.dat' using 1:2 with lines lw 1 title '2024-2025 interp', \
'quadratic_fit_curve.dat' using 1:2 with lines lw 3 title 'Quadratic fit curve', \
'actual_data.dat' using 1:2 with points pt 7 ps 1.8 title 'Actual data points', \
'all_predictions.dat' using 1:2 with points pt 5 ps 1.4 title 'All interpolation predictions', \
'all_predictions.dat' using 1:2:3 with labels offset 0.6,0.2 notitle, \
'best_predictions.dat' using 1:2 with points pt 9 ps 2.0 title 'Less extreme interpolation predictions', \
'best_predictions.dat' using 1:2:3 with labels offset 0.6,0.8 notitle, \
'fit_prediction.dat' using 1:2 with points pt 11 ps 2.5 title 'Quadratic-fit prediction', \
'fit_prediction.dat' using 1:2:3 with labels offset 0.6,1.0 notitle
