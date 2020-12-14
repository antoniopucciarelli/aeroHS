unset   key
set     xlabel '% x'
set     ylabel 'cp [NO UOD]'
set     title  'cp'
plot    'cp_dataGE.dat' using 1:2 with lines, 'cp_dataGE.dat' using 3:4 with lines
