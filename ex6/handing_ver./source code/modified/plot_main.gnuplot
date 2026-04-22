set terminal pngcairo size 1600,950 enhanced font 'Arial,14'
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
'quadratic_fit_curve.dat' using 1:2 with lines lw 3 lc rgb 'black' title 'Quadratic fit curve', \
'actual_data.dat' using 1:2 with points pt 7 ps 1.8 lc rgb 'purple' title 'Actual data points', \
'positive_predictions.dat' using 1:2 with points pt 5 ps 1.6 lc rgb 'red' title 'Positive interp predictions', \
'positive_predictions.dat' using 1:2:3 with labels offset 0.5,0.2 tc rgb 'red' notitle, \
'negative_predictions.dat' using 1:2 with points pt 5 ps 1.6 lc rgb 'blue' title 'Negative interp predictions', \
'negative_predictions.dat' using 1:2:3 with labels offset 0.5,-0.4 tc rgb 'blue' notitle, \
'best_predictions.dat' using 1:2 with points pt 9 ps 2.0 lc rgb 'orange' title 'Less extreme interp predictions', \
'fit_prediction.dat' using 1:2 with points pt 11 ps 2.5 lc rgb 'dark-green' title 'Quadratic-fit prediction', \
'fit_prediction.dat' using 1:2:3 with labels offset 0.5,0.8 tc rgb 'dark-green' notitle
