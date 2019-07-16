
#include <iostream>
#include "stdlib.h"
#include "time.h"
#include "pnetcdf.h"
#include "Array.h"


// User editable parameters
int    const nxglob   = 1000;    // Global number of cells in the x-direction
int    const nyglob   = 1000;    // Global number of cells in the y-direction
double const sparsity = 0.2;     // Initial sparsity of living cells
int    const nsteps   = 100;     // How many steps to run
int    const outFreq  = 1;       // How frequently to dump output to file in terms of steps


// Other simulation parameters
Array<unsigned char> grid;  // This process's local board of automata plus a halo of 1 cell in each direction for convenience
Array<unsigned char> next;  // Store the next step's data in a separate array so we don't overwrite the current grid while we're still computing the next step
Array<int> outData;         // For dumping out the data
int numOut;                 // Number of outputs so far
int nx;                     // The local number of x cells
int const ny = nyglob;      // The local number of y cells (I don't decompose in the y-dimension, so this is the global number)
int nranks;                 // The number of MPI ranks
int myrank;                 // My MPI Rank ID
int left_rank;              // The MPI rank ID of my left  neighbor (periodic domain)
int right_rank;             // The MPI rank ID of my right neighbor (periodic domain)
int i_beg, i_end;           // begining and ending global cell indices in the x-direction for my MPI rank
int const hs = 1;           // Halo size is always 1 cell
unsigned char sendbuf_l[ny+2];  // buffer for sending   data to   my left  neighbor
unsigned char sendbuf_r[ny+2];  // buffer for sending   data to   my right neighbor
unsigned char recvbuf_l[ny+2];  // buffer for receiving data from my left  neighbor
unsigned char recvbuf_r[ny+2];  // buffer for receiving data from my right neighbor


// Function definitions
unsigned char random_uniform(double const ratio);
void initialize();
void exchangeHalos();
void advance();
void ncwrap(int ierr, int line);
void output(int nstep);


// Main Driver
int main(int argc, char** argv) {
  int istep = 0;
  MPI_Init( &argc , &argv );   // Initialize MPI
  initialize();                // Initialize the model
  output(0);                   // Output the initial state, and create the file

  while (istep < nsteps) {
    exchangeHalos();           // Exchange halos with my left and right neighboring MPI tasks
    advance();                 // Advance one step based on the Game of Life rules
    istep++;

    // If it's time to output, the do output
    if (istep%outFreq == 0) {
      output(istep);
    }
  }

  MPI_Finalize();              // Finalize the MPI
}



// Generate a random alive / dead state based on a sparsity parameter
// If ratio is 0.1, then approximately 10% of the board will be alive at the beginning
unsigned char random_uniform(double const ratio) {
  float num = (double) rand() / (double) RAND_MAX;  // Random number \in [0,1]
  if (num < ratio) {
    return 1;
  } else {
    return 0;
  }
}



// Initialize the MPI and the grid of living and dead cells
void initialize() {
  srand(time(NULL));       // Set the initial random seed based on the time

  MPI_Comm_size(MPI_COMM_WORLD,&nranks);  // Get the number of MPI Ranks
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);  // Get my MPI Rank ID

  // Determine my rank's beginning and ending index in the x-direction
  double nper = ( (double) nxglob ) / nranks;   // Get the number of cells per rank
  i_beg = round( nper* (myrank)    );
  i_end = round( nper*((myrank)+1) )-1;

  // Determine my rank's number of x-direction cells to work on
  nx = i_end - i_beg + 1;

  // Determine my left and right neighbor rank IDs
  left_rank  = myrank - 1;
  if (left_rank == -1) left_rank = nranks-1;
  right_rank = myrank + 1;
  if (right_rank == nranks) right_rank = 0;

  // Allocate the local grid
  // Index zero is the left / bottom halo cell in the x / y direction
  // Indices [1,nx] and [1,ny] are the local MPI Rank's domain
  // Index nx+1 and ny+1 are the right / top halo cell in the x / y direction
  grid.setup(ny+2,nx+2);
  next.setup(ny+2,nx+2);

  // Randomly initialize this MPI rank's domain using the sparsity parameter
  for (int j=0; j<ny; j++) {
    for (int i=0; i<nx; i++) {
      grid(hs+j,hs+i) = random_uniform(sparsity);
    }
  }

  numOut = 0;
}



// Exchange halos between the different ranks in the x-direction
void exchangeHalos() {
  int ierr;
  MPI_Request req_r[2], req_s[2];

  //Prepost receives
  ierr = MPI_Irecv(recvbuf_l,ny+2,MPI_UNSIGNED_CHAR, left_rank,0,MPI_COMM_WORLD,&req_r[0]);
  ierr = MPI_Irecv(recvbuf_r,ny+2,MPI_UNSIGNED_CHAR,right_rank,1,MPI_COMM_WORLD,&req_r[1]);

  //Pack the send buffers
  for (int j=0; j<ny+2; j++) {
    sendbuf_l[j] = grid(j,1 );
    sendbuf_r[j] = grid(j,nx);
  }

  //Fire off the sends
  ierr = MPI_Isend(sendbuf_l,ny+2,MPI_UNSIGNED_CHAR, left_rank,1,MPI_COMM_WORLD,&req_s[0]);
  ierr = MPI_Isend(sendbuf_r,ny+2,MPI_UNSIGNED_CHAR,right_rank,0,MPI_COMM_WORLD,&req_s[1]);

  //Wait for receives to finish
  ierr = MPI_Waitall(2,req_r,MPI_STATUSES_IGNORE);

  //Unpack the receive buffers
  for (int j=0; j<ny+2; j++) {
    grid(j,0   ) = recvbuf_l[j];
    grid(j,nx+1) = recvbuf_r[j];
  }

  //Wait for sends to finish
  ierr = MPI_Waitall(2,req_s,MPI_STATUSES_IGNORE);

  // Exchange periodic halos in y-direction
  for (int i=0; i<nx+2; i++) {
    grid(0   ,i) = grid(ny,i);
    grid(ny+1,i) = grid(1 ,i);
  }
}



// Advance the game by one step based on the following rules
// 1. Any live cell with fewer than two live neighbours dies, as if by underpopulation.
// 2. Any live cell with two or three live neighbours lives on to the next generation.
// 3. Any live cell with more than three live neighbours dies, as if by overpopulation.
// 4. Any dead cell with three live neighbours becomes a live cell, as if by reproduction.
void advance() {

  // Calculate the next generation grid
  for (int j=0; j<ny; j++) {
    for (int i=0; i<nx; i++) {
      // Calculate the number of living neighbors among my 8 neighbors
      int numLivingNeighbors = 0;
      for (int jj=0; jj<3; jj++) {
        for (int ii=0; ii<3; ii++) {
          if (grid(j+jj,i+ii)) { numLivingNeighbors++; }
        }
      }
      // Take out my cell's contribution
      if (grid(hs+j,hs+i)) { numLivingNeighbors--; }

      if (grid(hs+j,hs+i)) {
        // My cell is currently alive
        if (numLivingNeighbors < 2) {
          next(hs+j,hs+i) = 0;    // Die by underpopulation
        } else if (numLivingNeighbors == 2 || numLivingNeighbors == 3) {
          next(hs+j,hs+i) = 1;    // Continue living in sustainable population
        } else {
          next(hs+j,hs+i) = 0;    // Die by overpopulation
        }
      } else {
        if (numLivingNeighbors == 3) {
          next(hs+j,hs+i) = 1;    // Alive by reproduction from neighbors
        } else {
          next(hs+j,hs+i) = 0;    // Continue to be dead
        }
      }
      
    }
  }

  // Store the next generation grid into the current grid
  for (int j=0; j<ny; j++) {
    for (int i=0; i<nx; i++) {
      grid(hs+j,hs+i) = next(hs+j,hs+i);
    }
  }

}



// Write the current grid out to a netCDF file in parallel
void output(int const nstep) {
  int dimids[3];
  MPI_Offset st[3], ct[3];
  int ncid, xDim, yDim, sDim, gVar, sVar;

  if (nstep == 0) {
    // Create the file
    ncwrap( ncmpi_create( MPI_COMM_WORLD , "output.nc" , NC_CLOBBER , MPI_INFO_NULL , &ncid ) , __LINE__ );

    // Define the dimensions
    ncwrap( ncmpi_def_dim( ncid , "step" , (MPI_Offset) NC_UNLIMITED , &sDim ) , __LINE__ );
    ncwrap( ncmpi_def_dim( ncid , "x"    , (MPI_Offset) nxglob       , &xDim ) , __LINE__ );
    ncwrap( ncmpi_def_dim( ncid , "y"    , (MPI_Offset) nyglob       , &yDim ) , __LINE__ );

    // Define the variables (step and cells)
    dimids[0] = sDim;
    ncwrap( ncmpi_def_var( ncid , "step"  , NC_INT , 1 , dimids , &sVar ) , __LINE__ );
    dimids[0] = sDim; dimids[1] = yDim; dimids[2] = xDim;
    ncwrap( ncmpi_def_var( ncid , "cells" , NC_INT , 3 , dimids , &gVar ) , __LINE__ );

    // End "define" mode
    ncwrap( ncmpi_enddef( ncid ) , __LINE__ );

    // Allocate the array for dumping out data
    outData.setup(ny,nx);
  } else {
    // Open the file
    ncwrap( ncmpi_open( MPI_COMM_WORLD , "output.nc" , NC_WRITE , MPI_INFO_NULL , &ncid ) , __LINE__ );

    // Get the variable IDs
    ncwrap( ncmpi_inq_varid( ncid , "step"  , &sVar  ) , __LINE__ );
    ncwrap( ncmpi_inq_varid( ncid , "cells" , &gVar  ) , __LINE__ );
  }

  for (int j=0; j<ny; j++) {
    for (int i=0; i<nx; i++) {
      outData(j,i) = grid(hs+j,hs+i);
    }
  }
  st[0] = numOut; st[1] = 0 ; st[2] = i_beg;
  ct[0] = 1     ; ct[1] = ny; ct[2] = nx   ;
  ncwrap( ncmpi_put_vara_int_all( ncid , gVar , st , ct , outData.get_data() ) , __LINE__ );

  ncwrap( ncmpi_begin_indep_data(ncid) , __LINE__ );
  st[0] = numOut;
  ncwrap( ncmpi_put_var1_int( ncid , sVar , st , &nstep ) , __LINE__ );
  ncwrap( ncmpi_end_indep_data(ncid) , __LINE__ );

  numOut++;
}



//Error reporting routine for the PNetCDF I/O
void ncwrap( int ierr , int line ) {
  if (ierr != NC_NOERR) {
    printf("NetCDF Error at line: %d\n", line);
    printf("%s\n",ncmpi_strerror(ierr));
    exit(-1);
  }
}

