#!/bin/bash

#job name:

#SBATCH --job-name=lam1

#

# Account:

#SBATCH --account=co_noneq

#

# Partition:

#SBATCH --partition=savio3

#

# Wall clock limit:

#SBATCH --time=496:00:00

#

# Number of MPI tasks needed for use case (example):

#SBATCH --ntasks=32
#

# Processors per task:

#SBATCH --cpus-per-task=1



mpirun -np 32 ./single 1.0 120 3200000 

