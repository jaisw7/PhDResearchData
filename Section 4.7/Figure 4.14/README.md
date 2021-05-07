This folder contains data for Fig. 4.14 of the thesis.

### Figures:  

The figures have been named in a format massDiffusion\_A.pdf where 
- A is the property of interest, for example T denotes temperature.

The individal figure files can be found at  
- [Fig. 4.14(a)](massDiffusion_nfrac.pdf)
- [Fig. 4.14(b)](selfDiffusion_vD.pdf)
- [Fig. 4.14(c)](selfDiffusion_T.pdf)


### Data files:  

The dsmc data files are contained in the directory dsmc. For this case, we vary the Knudsen number. Each of the sub-directory has been named in the format Kn=A\_fracAr0.68\_fracKr0.32\_alpha=1.4 where 
- A is the scattering index.

Here fracAr0.68 contains Ar with 68% concentration and therefore Kr with 32% concentration. 

The DGFS data files are contained in their own separate folders. Each of these folders follow the following format: Kn=A\_ArKr\_sBkCvDm12, where
- A is the Knudsen number of the flow.
- B is the number of elements in the physical space.
- C is the order of the discontinuous Galerkin scheme.
- D is the number of points in the velocity mesh.

The data files have been named in the format massDiffusion1D_x_speciey_z\_sol\_.txt where 
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

[plot_massDiffusion.tmd](plot_massDiffusion.tmd) is the script for generating figures 4.11, 4.12, and 4.13. See [tecCmd](https://github.com/jaisw7/tecCmd) for details on how to run this file.
