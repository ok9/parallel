#!/bin/bash
#SBATCH -p main
#SBATCH -n2
#SBATCH --time=0-5
module load openmpi/4.0

mpirun  ./optimize
