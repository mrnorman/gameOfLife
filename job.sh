#!/bin/bash
#SBATCH --time=10   #minutes
#SBATCH --nodes=1
#SBATCH -C haswell
#SBATCH --account=m1266
#SBATCH --reservation=csgfhpc_thu

srun -n 32 ./gameOfLife

