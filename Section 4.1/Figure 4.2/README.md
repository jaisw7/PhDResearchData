This folder contains data for Fig. 4.2 of the thesis.

### Figures: The figures have been named in a format variationF\_mr=x\_s=y_L=12.pdf where 
- x is the mass ratio of the two species
- y is the species (first species: s=1, second species: s=2)
- L=12 is the length of the velocity domain.

- [Fig. 4.2(a)](variationF_mr=2_s=1_L=12.pdf)
- [Fig. 4.2(b)](variationF_mr=2_s=2_L=12.pdf)
- [Fig. 4.2(c)](variationF_mr=4_s=1_L=12.pdf)
- [Fig. 4.2(d)](variationF_mr=4_s=2_L=12.pdf)
- [Fig. 4.2(e)](variationF_mr=8_s=1_L=12.pdf)
- [Fig. 4.2(f)](variationF_mr=8_s=2_L=12.pdf)

### Data files:  

The data files for *numerical* solution have been named in a format dist_Nv=64_Nsph=6_Nrho=64_mr=x.000000_t=y.txt where 
- x is the mass ratio of the two species.
- y is the time

In each of the data files, there are three columns. 
- Column 1: The velocity grid point. We use N=64. Therefore the first column varies from 1 to 64.
- Column 2: Value of distribution function for the first species.
- Column 3: Value of distribution function for the second species.


The data files for *analytical* solution have been named in a format dist_Nv=64_Nsph=6_Nrho=64_mr=x.000000_t=y.txt_analytical where 
- x is the mass ratio of the two species.
- y is the time

In each of the data files, there are three columns. 
- Column 1: The velocity grid point. We use N=64. Therefore the first column varies from 1 to 64.
- Column 2: Value of distribution function for the first species.
- Column 3: Value of distribution function for the second species.


### Plot generation 

[plot_variationF.tmd](plot_variationF.tmd) is the script for generating 1-D tecplot output. See [tecCmd](https://github.com/jaisw7/tecCmd) for details on how to run this file.

