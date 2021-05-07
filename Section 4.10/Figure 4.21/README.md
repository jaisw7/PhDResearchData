This folder contains data for Fig. 4.21 of the thesis.

### Figures:  

The figures have been named in a format Kn=A/B\_C.pdf where 
- A is the Knudsen number of the flow
- B is the property of interest, for example rho denotes mass-density, rhou denotes momentum-flux.
- C=boltz is the Boltzmann scattering kernel. 

The individal figure files can be found at  
- [Fig. 4.21(a)](Kn=1e-1/rho_boltz.pdf)
- [Fig. 4.21(b)](Kn=1e-1/rhou_boltz.pdf)
- [Fig. 4.21(c)](Kn=1e-1/energy_boltz.pdf)
- [Fig. 4.21(d)](Kn=1e-2/rho_boltz.pdf)
- [Fig. 4.21(e)](Kn=1e-2/rhou_boltz.pdf)
- [Fig. 4.21(f)](Kn=1e-2/energy_boltz.pdf)
- [Fig. 4.21(g)](Kn=1e-3/rho_boltz.pdf)
- [Fig. 4.21(h)](Kn=1e-3/rhou_boltz.pdf)
- [Fig. 4.21(i)](Kn=1e-3/energy_boltz.pdf)

### Data files:  

The NavierStokes and Euler folders contain output from OpenFOAM-2.3.1 post-processing tool *sample*. The results can be reproduced by rerunning their shock-tube case with rhoCentralFoam solver.
The data files are contained in the folders in the format A/B/sets/C/combined.xy, where 
- A is the Knudsen number of the flow
- B is the model, either Euler or NavierStokes
- C is the dimensional-time. In the [rhoCentralFoam](https://doi.org/10.1002/fld.2069) paper, the authors have plotted the results for time=0.007 for example. Therefore C=0.007

In each of the data files, there are eight columns. 
- x x-coordinate of the physical space mesh
- T Temperature
- p Pressure
- rho Mass-density
- x x-coordinate of the physical space mesh
- U:x x-component of velocity
- U:y y-component of velocity
- U:z z-component of velocity

The DGFS data files are contained in their own separate folders. Each of these folders follow the following format: Kn=A/bdf111-sBkCvDm12, where
- A is the Knudsen number of the flow.
- B is the number of elements in the physical space.
- C is the order of the discontinuous Galerkin scheme.
- D is the number of points in the velocity mesh.

The data files have been named in the format temp-x.txt where 
- x is the non-dimensional time. 

In each of the data files, there are eleven columns. 
x rho, U:x, U:y, T, Q:x, Q:y, P:xx, P:yy, P:xy, p
- x x-coordinate of the physical space mesh
- rho Mass-density 
- U:x x-component of velocity 
- U:y y-component of velocity 
- T Temperature 
- Q:x x-component of heat-flux 
- Q:y y-component of heat-flux 
- P:xx xx-component of stress tensor 
- P:yy yy-component of stress tensor
- P:xy xy-component of stress tensor
- p pressure 

All these properties are non-dimensional as defined in the thesis.

### Plot generation 

- [Kn=1e-1/plot.tmd](Kn=1e-1/plot.tmd) is the script for generating figures 4.20 (top row)
- [Kn=1e-2/plot.tmd](Kn=1e-2/plot.tmd) is the script for generating figures 4.20 (middle row)
- [Kn=1e-3/plot.tmd](Kn=1e-3/plot.tmd) is the script for generating figures 4.20 (bottom row)

See [tecCmd](https://github.com/jaisw7/tecCmd) for details on how to run this file.
