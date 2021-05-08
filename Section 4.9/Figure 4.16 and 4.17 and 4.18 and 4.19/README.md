This folder contains data for Fig. 4.16, 4.17, 4.18, and 4.19 of the thesis.

### Figures:  

The figures have been named in a format IMEX\_A\_B.pdf where 
- A is the name of the kinetic model, for example, Boltzmann or Relaxation kinetic model
- B is the property of interest, for example, T: Temperature, U: y-component of velocity

The individal figure files can be found at  
- [Fig. 4.16](IMEX_Relaxation_U.pdf)
- [Fig. 4.17](IMEX_Boltzmann_U.pdf)
- [Fig. 4.18](IMEX_Relaxation_T.pdf)
- [Fig. 4.19](IMEX_Boltzmann_T.pdf)

### Data files:  

The data files have been named in the format Kn=A/B\_C.txt where 
- A is the Knudsen number of the flow
- B is the mode of simulation: dgfs, dsmc, or Navier-Stokes
- C is the scattering model: bgk, esbgk, shakov

In each of the data files for dgfs, there are eleven columns. 
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

### Plot generation 

- [plot_Relaxation_U.tmd](plot_Relaxation_U.tmd) is the script for generating figures 4.16. 
- [plot_Boltzmann_U.tmd](plot_Boltzmann_U.tmd) is the script for generating figures 4.17. 
- [plot_Relaxation_T.tmd](plot_Relaxation_T.tmd) is the script for generating figures 4.18. 
- [plot_Boltzmann_T.tmd](plot_Boltzmann_T.tmd) is the script for generating figures 4.19. 

See [tecCmd](https://github.com/jaisw7/tecCmd) for details on how to run this file.
