This folder contains data for Fig. 4.15 of the thesis.

### Figures:  

The figures have been named in a format accuracy\_A\_B.pdf where 
- A is the time integration scheme.
- B is the name of the kinetic model

The individal figure files can be found at  
- [Fig. 4.15(a)](accuracy_bdf_bgk.pdf)
- [Fig. 4.15(b)](accuracy_bdf_esbgk.pdf)
- [Fig. 4.15(c)](accuracy_bdf_shakov.pdf)

### Data files:  

The data files have been named in the format data\_accuracy_A.txt where 
- A is the name of the kinetic model

In each of the data files, there is one column.
- The first six rows contain data from BDF-2 model 
- The last six rows contain data from BDF-3 model 

### Plot generation 

[plot_accuracy.tmd](plot_accuracy.tmd) is the script for generating figures 4.15. See [tecCmd](https://github.com/jaisw7/tecCmd) for details on how to run this file.
