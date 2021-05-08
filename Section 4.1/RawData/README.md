This folder contains a self contained Python script to run multi-species BKW solution on GPUs.

This is a Python-3 script that you can run by typing 
> python bkwbi_vss_cufft 1.0 16 6 4 gj 4.0 64

The script needs seven inputs:
- mr mass ratio of the two species
- N number of points in each direction of the velocity mesh
- M number of points on sphere
- Nrho Number of points in the radial direction
- quad Quadrature type used for integration
- t0 Initial time for Krook-Wu solution
- blkSize CUDA block size

Refer to the PhDThesis to understand the meaning of these parameters. 

The inputs are purely meant for benchmarking. You can write a simple bash loop to test the effect of input on the accuracy. To run the script, you will need to install the following Python modules installed on your system:
- mako
- pycuda
- numpy