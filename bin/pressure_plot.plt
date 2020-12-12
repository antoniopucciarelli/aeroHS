unset   key
set     xlabel '% x'
set     ylabel 'P [Pa]'
set     title  'pressure'
plot    "PRESSURE_gnuplot.dat" using 1:2 with lines, "PRESSURE_gnuplot.dat" using 3:4 with lines 