unset   key
set     xlabel 'alpha'
set     ylabel 'cl [NO UOD]'
set     title  'cl'
set     xrange [-10:10]
set     yrange [-3:3]
plot    'cl_gnuplot.dat' using 1:2 with lines