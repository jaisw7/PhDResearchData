This folder contains data for Fig. 4.30 and 4.32 of the thesis.

### Figures:  

The figures have been named in a format Kn=A/B_C.png where 
- A is the flow Knudsen number.
- B is the bulk property of interest.
- C is the simulation mode, for example, DGFS or DSMC.

The individal figure files can be found at  
- [Fig. 4.30(a)](Kn=0.1/Rho_DSMC.png)
- [Fig. 4.30(b)](Kn=0.1/Rho_DGFS.png)
- [Fig. 4.30(c)](Kn=0.1/Pxy_DSMC.png)
- [Fig. 4.30(d)](Kn=0.1/Pxy_DGFS.png)
- [Fig. 4.30(e)](Kn=0.1/Pyz_DSMC.png)
- [Fig. 4.30(f)](Kn=0.1/Pyz_DGFS.png)
- [Fig. 4.32(a)](Kn=1/Rho_DSMC.png)
- [Fig. 4.32(b)](Kn=1/Rho_DGFS.png)
- [Fig. 4.32(c)](Kn=1/Pxy_DSMC.png)
- [Fig. 4.32(d)](Kn=1/Pxy_DGFS.png)
- [Fig. 4.32(e)](Kn=1/Pyz_DSMC.png)
- [Fig. 4.32(f)](Kn=1/Pyz_DGFS.png)

### Data files:  

The DGFS data files have been named in a format Kn=A/sBvC/sol.vtu where 
- A is the flow Knudsen number.
- B is the number of elements in the physical space in each coordinate direction.
- C is the number of points in the velocity mesh.
These files can be opened in Tecplot.

The DSMC data files have been named in the format Kn=A/dsmc/sol.plt where 
- A is the flow Knudsen number.
These files can be opened in Tecplot.

### Plot generation 

- [Kn=1/Rho_DGFS.lay](Kn=1/Rho_DGFS.lay) is the Tecplot Layout file for generating figure 4.32(a). This has been generated from the data contained in Kn=1/dgfs/sol.vtu
- [Kn=1/Rho_DSMC.lay](Kn=1/Rho_DSMC.lay) is the Tecplot Layout file for generating figure 4.32(b). This has been generated from the data contained in Kn=1/dsmc/sol.plt
- [Kn=1/Pxy_DGFS.lay](Kn=1/Pxy_DGFS.lay) is the Tecplot Layout file for generating figure 4.32(c). This has been generated from the data contained in Kn=1/dgfs/sol.vtu
- [Kn=1/Pxy_DSMC.lay](Kn=1/Pxy_DSMC.lay) is the Tecplot Layout file for generating figure 4.32(d). This has been generated from the data contained in Kn=1/dsmc/sol.plt
- [Kn=1/Pyz_DGFS.lay](Kn=1/Pyz_DGFS.lay) is the Tecplot Layout file for generating figure 4.32(e). This has been generated from the data contained in Kn=1/dgfs/sol.vtu
- [Kn=1/Pyz_DSMC.lay](Kn=1/Pyz_DSMC.lay) is the Tecplot Layout file for generating figure 4.32(f). This has been generated from the data contained in Kn=1/dsmc/sol.plt

You may similarly go about Fig. 4.30. The files sol.vtu and sol.plt can be opened in Tecplot.