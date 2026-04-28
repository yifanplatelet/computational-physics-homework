set terminal pngcairo size 1400,1800
set output 'runge_comparison.png'
set multiplot layout 3,2 title 'Runge Function: Exact vs Polynomial Interpolation'
set grid
set xlabel 'x'
set ylabel 'y'
set key outside
set title 'n = 5 nodes, degree = 4'
plot 'runge_5_curve.txt' using 1:2 with lines lw 2 title 'Exact f(x)', 'runge_5_curve.txt' using 1:3 with lines lw 2 title 'Interpolation', 'runge_5_nodes.txt' using 1:2 with points pt 7 ps 1.2 title 'Nodes'
set title 'n = 7 nodes, degree = 6'
plot 'runge_7_curve.txt' using 1:2 with lines lw 2 title 'Exact f(x)', 'runge_7_curve.txt' using 1:3 with lines lw 2 title 'Interpolation', 'runge_7_nodes.txt' using 1:2 with points pt 7 ps 1.2 title 'Nodes'
set title 'n = 9 nodes, degree = 8'
plot 'runge_9_curve.txt' using 1:2 with lines lw 2 title 'Exact f(x)', 'runge_9_curve.txt' using 1:3 with lines lw 2 title 'Interpolation', 'runge_9_nodes.txt' using 1:2 with points pt 7 ps 1.2 title 'Nodes'
set title 'n = 17 nodes, degree = 16'
plot 'runge_17_curve.txt' using 1:2 with lines lw 2 title 'Exact f(x)', 'runge_17_curve.txt' using 1:3 with lines lw 2 title 'Interpolation', 'runge_17_nodes.txt' using 1:2 with points pt 7 ps 1.2 title 'Nodes'
set title 'n = 19 nodes, degree = 18'
plot 'runge_19_curve.txt' using 1:2 with lines lw 2 title 'Exact f(x)', 'runge_19_curve.txt' using 1:3 with lines lw 2 title 'Interpolation', 'runge_19_nodes.txt' using 1:2 with points pt 7 ps 1.2 title 'Nodes'
set title 'n = 21 nodes, degree = 20'
plot 'runge_21_curve.txt' using 1:2 with lines lw 2 title 'Exact f(x)', 'runge_21_curve.txt' using 1:3 with lines lw 2 title 'Interpolation', 'runge_21_nodes.txt' using 1:2 with points pt 7 ps 1.2 title 'Nodes'
unset multiplot
