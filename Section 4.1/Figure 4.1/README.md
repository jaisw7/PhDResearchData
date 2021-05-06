This folder contains data for Fig. 4.1 of the thesis.

### Figures:  

The figures have been named in a format errorVsMassRatio\_mr=x_L=12.pdf where 
- x is the mass ratio of the two species
- L=12 is the length of the velocity domain.

The individal figure files can be found at  
- [Fig. 4.1(a)](errorVsMassRatio_mr=1_L=12.pdf)
- [Fig. 4.1(b)](errorVsMassRatio_mr=2_L=12.pdf)
- [Fig. 4.1(c)](errorVsMassRatio_mr=4_L=12.pdf)
- [Fig. 4.1(d)](errorVsMassRatio_mr=8_L=12.pdf)

### Data files:  

The data files have been named in a format Linf\_mr=x.0\_M=y\_Nrho=Nv.txt where 
- x is the mass ratio of the two species
- y in number of points on the half sphere (Page 38 of the thesis).

In each of the data files, there are three columns. 
- Column 1: Number of points in velocity space (N)
- Column 2: Linf error between the numerical and analytical solution for the first species (equation 4.5 of thesis with i=1).
- Column 3: Linf error between the numerical and analytical solution for the second species (equation 4.5 of thesis with i=2).

### Plot generation 

[plot_errorVsMassRatio.tmd](plot_errorVsMassRatio.tmd) is the script for generating 1-D tecplot output. See [tecCmd](https://github.com/jaisw7/tecCmd) for details on how to run this file.
