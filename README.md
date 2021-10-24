# aeroHS_FORTRAN

This program computes the flow properties around **NACA 4 digits** airfoils. It uses the **Hess-Smith** method to solve the flow field.

It can solve 2 airfoils problems and ground effect problems (with both the image method and with a ground panelization method).

## Code compilation steps

In order to compile the code: **CMake**, **GNUplot**, **LAPACK** and **BLAS** libraries are *required*.

```bash
cmake .

make 

cd bin/

./aeroHS
```

## How to compute flow properties

### Physics

The Hess-Smith method treats **incompressible inviscid potential** flows:
<img src="https://render.githubusercontent.com/render/math?math=\color{red} \boldsymbol{\nabla \cdot u}  = 0")>

<img src="https://render.githubusercontent.com/render/math?math=\color{red}\boldsymbol{\nabla \times u = 0}">

This brings solving the Laplace problem for the scalar function $\phi_{(x, y)}$:

$$
    \text{If } \boldsymbol{u = \nabla} \phi \rightarrow \boldsymbol{\nabla \times ( \nabla} \phi \boldsymbol{) \equiv 0} \text{, so the irrotationality of the flow is satisfied by } \phi \\
    \text{In order to satisfy incompressibility } \boldsymbol{\nabla \cdot (\nabla} \phi \boldsymbol{)} = \boldsymbol{\Delta} \phi = 0
$$

In order to solve the problem, it is used a first order discretization for the potential $\phi$ made by source/sink and vortex distribution on airfoil's panels. The Hess-Smith method uses a constant distribution of vorticity over the panels and a variable distribution of sources/sinks. The Runge-Kutta condition is applied at the trailing edge of the airfoil(s) in order to solve the problem.

### Program

Flow simulation and $C_P$ study of **NACA0012**.

Once in ``` bin/ ```:

```bash
# launching program
./aeroHS
# selecting analysis -> 1 airfoil analysis
1 
# selecting Cp analysis
1
# selecting AOA
0
# selecting airfoil name
naca0012
# selecting chord discretization
120
# leading edge position setup
# selecting leading edge x position
0
# selecting leading edge y position 
0
# grid setup coords for flow visualization
# selecting grid initial x coord
-1 
# selecting grid final x coord
3
# selecting grid initial y coord
-0.5
# selecting grid final y coord
0.5
# flow plot with GNUplot
```
