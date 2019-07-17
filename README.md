# The Game of Life

This code implements the Game of Life <https://en.wikipedia.org/wiki/Conway%27s_Game_of_Life> in C++ and MPI using a custom Multi-Dimensional Array class (Array.h), and outputs results to file in netCDF format using the parallel-netcdf library <https://trac.mcs.anl.gov/projects/parallel-netcdf> 

To build the code, edit the `Makefile` to have the correct path to the parallel-netcdf library on your machine, and then type `make`. To run the code, type

`mpirun -n [# tasks] ./gameOfLife`

This will run the code in parallel. To easily view the results, you can use the `ncview` tool <http://meteora.ucsd.edu/~pierce/ncview_home_page.html>.

## Data Structure Indexing

I made the choice of using so-called "halo" cells to make neighbor calculations easier. Different MPI tasks will send their edge data to their neighbors to fill the neighboring halo cells. Halo cells are also used to implement the periodic boundary conditions. I have the left / bottom halo cell as index 0 in each dimension. Then the current MPI task's domain is in the indices `1,...,nx` and `1,...,ny`. Finally, the right / top halo cell is at index `nx+1` and `ny+1`. The left halo belongs to the left neighbor's right-most cell. The right halo belongs to the right neighbor's left-most cell. The same goes for the bottom and top. In 2-D, the grid looks like:

```
0000000000
0********0
0********0
0********0
0********0
0********0
0********0
0000000000
```

where a `0` is a halo cell and a `*` is a domain cell for this MPI task.
