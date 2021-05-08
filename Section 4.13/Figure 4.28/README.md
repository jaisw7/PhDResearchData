This folder contains data for Fig. 4.28 of the thesis.

### Figures:  

The figures have been named in a format Kn=A/residual.pdf where 
- A is the flow Knudsen number.

The individal figure files can be found at  
- [Fig. 4.28(a)](Kn=0.1/residual.pdf)
- [Fig. 4.28(b)](Kn=1/residual.pdf)


### Data files:  

The data files have been named in a format Kn=A\_sBvC/residual.txt where 
- A is the flow Knudsen number.
- B is the number of elements in the physical space.
- C is the number of points in the velocity mesh.

The data file is essentially the output from my solver [frfs](https://github.com/jaisw7/frfs)
Observe the lines in the output file that are in the format dgfsresidualstd :  A, B
- A is the non-dimensional time
- B is the distribution norm

We are essentially plotting A on x-axis and B on y-axis. Don't worry. The plotting software will do the work for you.

### Plot generation 

- [Kn=0.1/plot_residual.tmd](Kn=0.1/plot_residual.tmd) is the script for generating figure 4.28(a). 
- [Kn=1/plot_residual.tmd](Kn=1/plot_residual.tmd) is the script for generating figure 4.28(b). 

See [tecCmd](https://github.com/jaisw7/tecCmd) for details on how to run this file.
