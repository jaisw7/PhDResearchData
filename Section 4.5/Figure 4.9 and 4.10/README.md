This folder contains data for Fig. 4.9 and 4.10 of the thesis.

### Figures:  

The figures have been named in a format rhoTU\_Ma=A\_sBkCvD\_cE\_frac=F\_mr=G\_specieH.pdf where 
- A is the flow Mach number.
- B is the number of elements in the physical space.
- C is the order of the discontinuous Galerkin scheme.
- D is the number of points in the velocity mesh.
- E is the width of the velocity mesh.
- F is the ratio of concentration of species two and one.
- G is the mass ratio of the two species.
- H={1,2} is the specie number.

The individal figure files can be found at  
- [Fig. 4.9(a)](rhoTU_Ma=1.5_s8k3v32_c9_frac=0.5_mr=0.5_specie1.pdf)
- [Fig. 4.9(b)](rhoTU_Ma=1.5_s8k3v32_c9_frac=0.5_mr=0.5_specie2.pdf)
- [Fig. 4.9(c)](rhoTU_Ma=1.5_s16k3v32_c9_frac=0.5_mr=0.5_specie1.pdf)
- [Fig. 4.9(d)](rhoTU_Ma=1.5_s16k3v32_c9_frac=0.5_mr=0.5_specie2.pdf)
- [Fig. 4.9(e)](rhoTU_Ma=1.5_s16k3v48_c15_frac=0.5_mr=0.5_specie1.pdf)
- [Fig. 4.9(f)](rhoTU_Ma=1.5_s16k3v48_c15_frac=0.5_mr=0.5_specie2.pdf)
- [Fig. 4.9(g)](rhoTU_Ma=1.5_s16k3v32_c9_frac=0.1_mr=0.5_specie1.pdf)
- [Fig. 4.9(h)](rhoTU_Ma=1.5_s16k3v32_c9_frac=0.1_mr=0.5_specie2.pdf)
- [Fig. 4.9(g)](rhoTU_Ma=3_s16k3v48_c15_frac=0.1_mr=0.5_specie1.pdf)
- [Fig. 4.9(h)](rhoTU_Ma=3_s16k3v48_c15_frac=0.1_mr=0.5_specie2.pdf)


### Data files:  

The folder Kosuge2001 contains data extracted from [Kosuge 2001](https://doi.org/10.1016/S0997-7546(00)00133-3)

The data files are contained in their separate directories named in the format sAkBvCm6_cD_MaE_He_30mm_fracBF_mr=G where 
- A is the number of elements in the physical space.
- B is the order of the discontinuous Galerkin scheme.
- C is the number of points in the velocity mesh.
- D is the width of the velocity mesh.
- E is the flow Mach number.
- F is the ratio of concentration of species two and one.
- G is the mass ratio of the two species.

Here 30mm denotes the width of the physical domain. 

The data file has been named in the format normalShockKosuge1D_x_speciey_z\_sol\_.txt where 
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

- [plot_Ma=1.5_He_HS_30mm_shifted_nfracB=0.5_mr=0.5.tmd](plot_Ma=1.5_He_HS_30mm_shifted_nfracB=0.5_mr=0.5.tmd) is the script for generating figures 4.9(a,b,c,d). 
- [plot_Ma=1.5_He_HS_30mm_shifted_nfracB=0.5_mr=0.25.tmd](plot_Ma=1.5_He_HS_30mm_shifted_nfracB=0.5_mr=0.25.tmd) is the script for generating figures 4.9(e,f). 
- [plot_Ma=1.5_He_HS_30mm_shifted_nfracB=0.1_mr=0.5.tmd](plot_Ma=1.5_He_HS_30mm_shifted_nfracB=0.1_mr=0.5.tmd) is the script for generating figures 4.10(a,b). 
- [plot_Ma=3_He_HS_30mm_shifted_nfracB=0.1_mr=0.5.tmd](plot_Ma=3_He_HS_30mm_shifted_nfracB=0.1_mr=0.5.tmd) is the script for generating figures 4.10(c,d). 

See [tecCmd](https://github.com/jaisw7/tecCmd) for details on how to run this file.
