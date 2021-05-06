This folder contains data for Fig. 4.4 of the thesis.

### Figures:  

The figures have been named in a format couette_x_s4k3v32m12_normalized_.pdf where 
- x is property of interest, for example: U refers to y-component of velocity, and T refers to temperature.
- s4k3v32m12 refers to s=4 elements in *s*pace, k=3rd order DG, v=32 points in each direction of *v*elocity domain, and m=12 points on sphere.

The individal figure files can be found at  
- [Fig. 4.4(a)](couette_U_s4k3v32m12_normalized_.pdf)
- [Fig. 4.4(b)](couette_T_s4k3v32m12_normalized_.pdf)

### Data files:  

The data file has been named in the format bulksol_dgfs_1d_couette-x.00.txt where 
- x is the non-dimensional time.

In each of the data files, there are fourteen columns. 
- rho:0 Mass-density for the first species
- U:x :0 x-component of velocity for the first species
- U:y:0 y-component of velocity for the first species
- T:0 Temperature for the first species
- Q:x :0 x-component of heat-flux for the first species
- p:0 pressure for the first species
- nden:0 Number density for the first species
- rho:1 Mass-density for the second species
- U:x :1 x-component of velocity for the second species
- U:y:1 y-component of velocity for the second species
- T:1 Temperature for the second species
- Q:x :1 x-component of heat-flux for the second species
- p:1 pressure for the second species
- nden:1 Number density for the second species

All these properties are non-dimensional as per the convention defined in the thesis.

### Plot generation 

[plot_couette.tmd](plot_couette.tmd) is the script for generating 1-D tecplot output. See [tecCmd](https://github.com/jaisw7/tecCmd) for details on how to run this file.
