set contour
set xlabel 'x'
set ylabel 'y'

plot "VELfield.dat" u 1:2:5 w image, 'GNUplot_tg_norm.dat' using 1:2 with lines
