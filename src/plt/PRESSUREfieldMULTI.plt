set contour
set xlabel 'x'
set ylabel 'y'

plot "FLOWfieldMULTI.dat" u 1:2:6 w image, 'GNUplot_tg_norm_1.dat' using 1:2 with lines, 'GNUplot_tg_norm_2.dat' using 1:2 with lines
