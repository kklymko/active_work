#!/bin/bash

#job name:

#SBATCH --job-name=lam-1.00

#

# Account:

#SBATCH --account=co_noneq

#

# Partition:

#SBATCH --partition=savio3

#

# Wall clock limit:

#SBATCH --time=848:00:00

#

# Number of MPI tasks needed for use case (example):

#SBATCH --ntasks=32
#

# Processors per task:

#SBATCH --cpus-per-task=1

mpirun -np 32 ./single -1.0 10 3200000
