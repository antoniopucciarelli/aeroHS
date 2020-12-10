unset   key
set     xlabel 'x'
set     ylabel 'y'
set     title  'airfoil elements'
set     size    ratio -1
unset   xtics
unset   ytics
plot    'GNUplot_th_norm.dat' using 1:2 with lines, 'GNUplot_tg_norm.dat' using 1:2:5:6 with vectors, 'GNUplot_tg_norm.dat' using 1:2:3:4 with vectors, 