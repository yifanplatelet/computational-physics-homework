set datafile separator ','
set term pngcairo size 1400,900
set grid
set key outside

set output 'zoom1_float.png'
set title 'Zoom 1: x in [0.99, 1.01] (float)'
plot 'poly_zoom1.csv' using 1:2 with lines title '(x-1)^7', \
     'poly_zoom1.csv' using 1:3 with lines title 'expanded'

set output 'zoom1_double.png'
set title 'Zoom 1: x in [0.99, 1.01] (double)'
plot 'poly_zoom1.csv' using 1:4 with lines title '(x-1)^7', \
     'poly_zoom1.csv' using 1:5 with lines title 'expanded'

set output 'zoom1_quad.png'
set title 'Zoom 1: x in [0.99, 1.01] (quad)'
plot 'poly_zoom1.csv' using 1:6 with lines title '(x-1)^7', \
     'poly_zoom1.csv' using 1:7 with lines title 'expanded'

set output 'zoom2_float.png'
set title 'Zoom 2: x in [0.9999, 1.0001] (float)'
plot 'poly_zoom2.csv' using 1:2 with lines title '(x-1)^7', \
     'poly_zoom2.csv' using 1:3 with lines title 'expanded'

set output 'zoom2_double.png'
set title 'Zoom 2: x in [0.9999, 1.0001] (double)'
plot 'poly_zoom2.csv' using 1:4 with lines title '(x-1)^7', \
     'poly_zoom2.csv' using 1:5 with lines title 'expanded'

set output 'zoom2_quad.png'
set title 'Zoom 2: x in [0.9999, 1.0001] (quad)'
plot 'poly_zoom2.csv' using 1:6 with lines title '(x-1)^7', \
     'poly_zoom2.csv' using 1:7 with lines title 'expanded'
