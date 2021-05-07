This folder contains data for Fig. 4.25 of the thesis.

### Figures:  

The figures have been named in a format A\_B\_hor\_centerline.pdf where 
- A is the property plotted on the left y-axis of the figure
- B is the property plotted on the right y-axis of the figure
hor\_centerline refers to horizontal centerline of the domain.

The individal figure files can be found at  
- [Fig. 4.25(a)](nden_T_hor_centerline.pdf)
- [Fig. 4.25(b)](Speed_Pxy_hor_centerline.pdf)

### Data files:  

The DGFS data file has thirteen columns. 
- x x-coordinate of the physical space mesh
- y y-coordinate of the physical space mesh
- z z-coordinate of the physical space mesh
- rho Mass-density 
- U:x x-component of velocity 
- U:y y-component of velocity 
- T Temperature 
- Q:x x-component of heat-flux 
- Q:y y-component of heat-flux 
- P:xx xx-component of stress tensor 
- P:xy xy-component of stress tensor
- P:yy yy-component of stress tensor
- p pressure 

All these properties are non-dimensional as defined in the thesis.


The DSMC data file has nineteen columns. 
- x x-coordinate of the physical space mesh
- y y-coordinate of the physical space mesh
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
- P:xx xx-component of stress tensor 
- P:yy yy-component of stress tensor
- Q:x x-component of heat-flux 
- Q:y y-component of heat-flux 

### Plot generation 

[plot.tmd](plot.tmd) is the script for generating figures 4.25. See [tecCmd](https://github.com/jaisw7/tecCmd) for details on how to run this file.
