This folder contains data for Fig. 4.7 of the thesis.

### Figures:  

The figures have been named in a format T\_x\_s4k3v64\_Ty_deviation where 
- x={Ar,Kr} is the specie name.
- y={20,100} is is the temperature difference between the two ends of the wall.

The individal figure files can be found at  
- [Fig. 4.7(a)](T_Ar_s4k3v64_T20_deviation.pdf)
- [Fig. 4.7(b)](T_Kr_s4k3v64_T20_deviation.pdf)
- [Fig. 4.7(c)](T_Ar_s4k3v64_T100_deviation.pdf)
- [Fig. 4.7(d)](T_Kr_s4k3v64_T100_deviation.pdf)

### Data files:  

The data files are contained in their separate directories named in the format Kn=x\_ArKr\_s4k3v64m12\_Ty where 
- x={0.5,1,5} is the Knudsen number.
- y={20,100} is the temperature difference between the two ends of the wall.

The data file has been named in the format fourier1D\_x\_y\_sol\_.txt where 
- x={Ar,Kr} is the species name
- y is the non-dimensional time. 

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

[plot_fourier20.tmd](plot_fourier20.tmd) is the script for generating figures 4.7(a) and 4.7(b). [plot_fourier100.tmd](plot_fourier100.tmd) is the script for generating figures 4.7(c) and 4.7(d). See [tecCmd](https://github.com/jaisw7/tecCmd) for details on how to run this file.
