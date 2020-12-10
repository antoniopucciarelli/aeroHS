unset   key
set     xlabel '% x'
set     ylabel 'V [m/s]'
set     title  'velocity'
plot    "VEL_gnuplot.dat" using 1:2 with lines, "VEL_gnuplot.dat" using 3:4 with lines 