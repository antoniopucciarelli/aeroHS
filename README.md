# aeroHS_FORTRAN
This program computes the flow properties around NACA4digits airfoil. 

It uses the Hess-Smith method to solve the flow field. 

It can solve 2 airfoils problems and ground effect problems(with the image method and also with a ground panelization method). 

In order to compile the program: GNUplot, LAPACK and BLAS libraries are required. 

> cmake .

> make 

> cd bin/

> ./aeroHS
