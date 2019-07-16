# The Game of Life

This code implements the Game of Life <https://en.wikipedia.org/wiki/Conway%27s_Game_of_Life> in C++ and MPI using a custom Multi-Dimensional Array class (Array.h), and outputs results to file in netCDF format using the parallel-netcdf library <https://trac.mcs.anl.gov/projects/parallel-netcdf> 

To build the code, edit the `Makefile` to have the correct path to the parallel-netcdf library on your machine, and then type `make`. To run the code, type

`mpirun -n [# tasks] ./gameOfLife`

This will run the code in parallel. To easily view the results, you can use the `ncview` tool <http://meteora.ucsd.edu/~pierce/ncview_home_page.html>.

