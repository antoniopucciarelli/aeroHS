set xlabel 'x'
set ylabel 'Cp'

plot 'cp_dataMULTI1.dat' using 1:2 with lines, 'cp_dataMULTI1.dat' using 3:4 with lines, 'cp_dataMULTI2.dat' using 1:2 with lines, 'cp_dataMULTI2.dat' using 3:4 with lines 
