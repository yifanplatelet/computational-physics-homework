set terminal pngcairo size 1200,800
set output 'result.png'
set title 'Zener diode V-I: fitting/interpolation comparison'
set xlabel 'Voltage'
set ylabel 'Current'
set grid
plot \
'compare_data.txt' using 1:2 with lines lw 2 title 'Polynomial degree 6', \
'compare_data.txt' using 1:3 with lines lw 2 title 'Natural cubic spline', \
'compare_data.txt' using 1:4 with lines lw 2 title 'Piecewise linear', \
'compare_data.txt' using 1:5 with lines lw 2 title 'PCHIP', \
'compare_data.txt' using 1:6 with lines lw 2 title 'Physics-inspired', \
'raw_points.txt' using 1:2 with points pt 7 ps 1.5 title 'Data points'
