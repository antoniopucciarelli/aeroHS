unset   key
set     xlabel '% x'
set     ylabel 'cp [NO UOD]'
set     title  'cp'
plot    'cp_gnuplot.dat' using 1:2 with lines, 'cp_gnuplot.dat' using 3:4 with lines