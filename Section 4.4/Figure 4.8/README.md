This folder contains data for Fig. 4.7 of the thesis.

### Figures:  

The figures have been named in a format U\_x\_Uy\_dt=2e-11_normalized where 
- x={Ar,Kr} is the specie name.
- y={50,500} is is the y-velocity difference between the two wall.

The individal figure files can be found at  
- [Fig. 4.7(a)](U_Ar_U50_dt=2e-11_normalized.pdf)
- [Fig. 4.7(b)](U_Kr_U50_dt=2e-11_normalized.pdf)
- [Fig. 4.7(c)](U_Ar_U500_dt=2e-11_normalized.pdf)
- [Fig. 4.7(d)](U_Kr_U500_dt=2e-11_normalized.pdf)

### Data files:  

The data files are contained in their separate directories named in the format Kn=1_ArKr_Ux where 
- x={50,500} is is the y-velocity difference between the two wall.

The data file has been named in the format oscCouette1DReconstruct_x_speciey_z\_sol\_.txt where 
- x={Ar,Kr} is the species name
- y={1,2} is the species number
- z is the non-dimensional time. 

In each of the data files, there are seven columns. 
- x x-coordinate of the physical space mesh
- rho:0 Mass-density for the first species
- U:x :0 x-component of velocity for the first species
- U:y:0 y-component of velocity for the first species
- T:0 Temperature for the first species
- p:0 pressure for the first species
- Q:x :0 x-component of heat-flux for the first species

All these properties are non-dimensional as defined in the thesis.

### Plot generation 

[plot_oscCouette.tmd](plot_oscCouette.tmd) is the script for generating figures. See [tecCmd](https://github.com/jaisw7/tecCmd) for details on how to run this file.
