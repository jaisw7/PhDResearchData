This folder contains data for Fig. 4.29 and 4.31 of the thesis.

### Figures:  

The figures have been named in a format Kn=A/residual.pdf where 
- A is the flow Knudsen number.

The individal figure files can be found at  
- [Fig. 4.29(a)](Kn=0.1/Pxx_diag.pdf)
- [Fig. 4.29(b)](Kn=0.1/Pxy_diag.pdf)
- [Fig. 4.29(c)](Kn=0.1/Pyy_diag.pdf)
- [Fig. 4.29(d)](Kn=0.1/Pzz_diag.pdf)
- [Fig. 4.29(e)](Kn=0.1/Pxz_diag.pdf)
- [Fig. 4.29(f)](Kn=0.1/Pyz_diag.pdf)
- [Fig. 4.31(a)](Kn=1/Pxx_diag.pdf)
- [Fig. 4.31(b)](Kn=1/Pxy_diag.pdf)
- [Fig. 4.31(c)](Kn=1/Pyy_diag.pdf)
- [Fig. 4.31(d)](Kn=1/Pzz_diag.pdf)
- [Fig. 4.31(e)](Kn=1/Pxz_diag.pdf)
- [Fig. 4.31(f)](Kn=1/Pyz_diag.pdf) 

### Data files:  

The DGFS data files have been named in a format Kn=A/sBvC/data\_line\_1,1,0-0,0,1.txt where 
- A is the flow Knudsen number.
- B is the number of elements in the physical space.
- C is the number of points in the velocity mesh.

The data files have twenty-two columns. 
- x x-coordinate of the physical space mesh
- y y-coordinate of the physical space mesh
- z z-coordinate of the physical space mesh
- rho Mass-density 
- U:x x-component of velocity 
- U:y y-component of velocity 
- U:z z-component of velocity 
- T Temperature 
- Q:x x-component of heat-flux 
- Q:y y-component of heat-flux 
- Q:z z-component of heat-flux 
- P:xx xx-component of stress tensor 
- P:xy xy-component of stress tensor
- P:xz xz-component of stress tensor
- P:yy yy-component of stress tensor
- P:yz yz-component of stress tensor
- P:zz zz-component of stress tensor
- p pressure 
- x x-coordinate of the physical space mesh in milli-meters
- y y-coordinate of the physical space mesh in milli-meters
- z z-coordinate of the physical space mesh in milli-meters
- arc length

All these properties are non-dimensional as defined in the thesis.


The DSMC data files have been named in the format Kn=A/dsmc/data\_line\_1,1,0-0,0,1.txt where 
- A is the flow Knudsen number.

These data files have twenty-eight columns. 
- x x-coordinate of the physical space mesh
- y y-coordinate of the physical space mesh
- z z-coordinate of the physical space mesh
- id mesh cell id
- xlo x-coordinate of the left edge of a square cell in physical space mesh (see SPARTA manual)
- ylo y-coordinate of the left edge of a square cell in physical space mesh (see SPARTA manual)
- xhi x-coordinate of the right edge of a square cell in physical space mesh (see SPARTA manual)
- yhi y-coordinate of the right edge of a square cell in physical space mesh (see SPARTA manual)
- n number of particles
- nden number density 
- rho Mass-density 
- U:x x-component of velocity 
- U:y y-component of velocity 
- T Temperature 
- p pressure 
- P:xy xy-component of stress tensor
- P:yz yz-component of stress tensor
- P:xz xz-component of stress tensor
- P:xx xx-component of stress tensor 
- P:yy yy-component of stress tensor
- P:zz zz-component of stress tensor
- Q:x x-component of heat-flux 
- Q:y y-component of heat-flux 
- Q:z z-component of heat-flux 
- x x-coordinate of the physical space mesh in milli-meters
- y y-coordinate of the physical space mesh in milli-meters
- z z-coordinate of the physical space mesh in milli-meters
- arc length

### Plot generation 

- [Kn=0.1/plot.tmd](Kn=0.1/plot.tmd) is the script for generating figure 4.29. 
- [Kn=1/plot.tmd](Kn=1/plot.tmd) is the script for generating figure 4.31. 

See [tecCmd](https://github.com/jaisw7/tecCmd) for details on how to run this file.
