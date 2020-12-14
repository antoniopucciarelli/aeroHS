set contour
set xlabel 'x'
set ylabel 'y'

plot "FLOWfieldGE.dat" u 1:2:6 w image, 'GNUplot_tg_norm.dat' using 1:2 with lines
