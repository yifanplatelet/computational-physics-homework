set terminal pngcairo size 1400,1600
set output 'problem3_comparison.png'
set multiplot layout 3,2 title 'f(x)=x sin(2*pi*x+1): exact vs polynomial interpolation'
set grid
set xlabel 'x'
set ylabel 'y'
set key outside
set title 'n = 7 nodes, degree = 6'
plot 'problem3_7_curve.txt' using 1:2 with lines lw 2 title 'Exact f(x)', 'problem3_7_curve.txt' using 1:3 with lines lw 2 title 'Interpolation', 'problem3_7_nodes.txt' using 1:2 with points pt 7 ps 1.2 title 'Nodes'
set title 'n = 9 nodes, degree = 8'
plot 'problem3_9_curve.txt' using 1:2 with lines lw 2 title 'Exact f(x)', 'problem3_9_curve.txt' using 1:3 with lines lw 2 title 'Interpolation', 'problem3_9_nodes.txt' using 1:2 with points pt 7 ps 1.2 title 'Nodes'
set title 'n = 17 nodes, degree = 16'
plot 'problem3_17_curve.txt' using 1:2 with lines lw 2 title 'Exact f(x)', 'problem3_17_curve.txt' using 1:3 with lines lw 2 title 'Interpolation', 'problem3_17_nodes.txt' using 1:2 with points pt 7 ps 1.2 title 'Nodes'
set title 'n = 19 nodes, degree = 18'
plot 'problem3_19_curve.txt' using 1:2 with lines lw 2 title 'Exact f(x)', 'problem3_19_curve.txt' using 1:3 with lines lw 2 title 'Interpolation', 'problem3_19_nodes.txt' using 1:2 with points pt 7 ps 1.2 title 'Nodes'
set title 'n = 21 nodes, degree = 20'
plot 'problem3_21_curve.txt' using 1:2 with lines lw 2 title 'Exact f(x)', 'problem3_21_curve.txt' using 1:3 with lines lw 2 title 'Interpolation', 'problem3_21_nodes.txt' using 1:2 with points pt 7 ps 1.2 title 'Nodes'
unset multiplot
