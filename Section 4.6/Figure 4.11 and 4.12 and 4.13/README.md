This folder contains data for Fig. 4.11, 4.12 and 4.13 of the thesis.

### Figures:  

The figures have been named in a format selfDiffusion\_A\_sBkCvDm12.pdf where 
- A is the property of interest, for example T denotes temperature; whereas D11_Case01 denotes diffusion coefficient for case 01.
- B is the number of elements in the physical space.
- C is the order of the discontinuous Galerkin scheme.
- D is the number of points in the velocity mesh.

The individal figure files can be found at  
- [Fig. 4.11(a)](selfDiffusion_nfrac_s4k3v32m12.pdf)
- [Fig. 4.11(b)](selfDiffusion_vD_s4k3v32m12.pdf)
- [Fig. 4.12(a)](selfDiffusion_T_s4k3v32m12.pdf)
- [Fig. 4.12(b)](selfDiffusion_T_s8k4v32m12.pdf)
- [Fig. 4.13(a)](selfDiffusion_D11_Case01_s4k3v32m12.pdf)
- [Fig. 4.13(b)](selfDiffusion_D11_Case02_s4k3v32m12.pdf)
- [Fig. 4.13(c)](selfDiffusion_D11_Case01_s4k4v32m12.pdf)
- [Fig. 4.13(d)](selfDiffusion_D11_Case02_s4k4v32m12.pdf)
- [Fig. 4.13(e)](selfDiffusion_D11_Case01_s8k4v32m12.pdf)
- [Fig. 4.13(f)](selfDiffusion_D11_Case02_s8k4v32m12.pdf)


### Data files:  

The dsmc data files are contained in the directory dsmc. For this case, we vary the Knudsen number and scattering index. Each of the sub-directory has been named in the format ppc2000_alpha=A_Kn=B where 
- A is the scattering index.
- B is the Knudsen number of the flow.

Here ppc2000 denotes 2000 particles per cell. This is done to reduce statistical noise. 

The DGFS data files are contained in their own separate folders. Each of these folders follow the following format: Kn=A\_ArAr\_sBkCvDm12\_alpha=E, where
- A is the Knudsen number of the flow.
- B is the number of elements in the physical space.
- C is the order of the discontinuous Galerkin scheme.
- D is the number of points in the velocity mesh.
- E is the scattering index.

The data files have been named in the format selfDiffusion1D_x_speciey_z\_sol\_.txt where 
- x={Ar,Ar} is the species name
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

[plot_selfDiffusion.tmd](plot_selfDiffusion.tmd) is the script for generating figures 4.11, 4.12, and 4.13. See [tecCmd](https://github.com/jaisw7/tecCmd) for details on how to run this file.
