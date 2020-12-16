unset  key 
set    xlabel 'x'
set    ylabel 'y'

plot 'GNUplot_tg_norm1.dat' using 1:2 with lines, 'GNUplot_tg_norm2.dat' using 1:2 with lines
