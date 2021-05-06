This folder contains data for Fig. 4.6 of the thesis.

### Figures:  

The figures have been named in a format bar\_x\_k\_y.pdf where 
- x=72 is the number of elements in the physical space.
- y={20,32} is the number of points in the velocity space.

The individal figure files can be found at  
- [Fig. 4.6(a)](bar_72_k_20.pdf)
- [Fig. 4.6(b)](bar_72_k_32.pdf)

### Data files:  

The data file has been named in the format abc\_x.k.y.txt where 
- x=72 is the number of elements in the physical space.
- y={20,32} is the number of points in the velocity space.

In each of the data files, there are six columns. 
- K: the order of DG scheme
- FFT: the time taken by FFT calls 
- cosSinMul: time taken by Step 4 of the algorithm
- prod: time taken by Step 6 of the algorithm
- computeQG: time taken by Step 7 of the algorithm
- Others/Remain: time taken by other steps of the algorithm.

All these steps have been defined in Paper [DGFS Multi-species GPU](http://doi.acm.org/10.1145/3324989.3325714)

### Plot generation 

[plot_bar.tmd](plot_bar.tmd) is the script for generating 1-D tecplot output. See [tecCmd](https://github.com/jaisw7/tecCmd) for details on how to run this file.
