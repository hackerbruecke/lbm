#!/bin/bash
# Job file for Stampede supercomputer
#SBATCH -J LBM_MPI # job name
#SBATCH -o LBM_MPI_OUT_%j # output file name (%j expands to jobID)
#SBATCH -e LBM_MPI_ERR_%j # error file name (%j expands to jobID)
#SBATCH -N 1                # total number of nodes requested (16 cores/node)
#SBATCH -n 4                # 1 task
#SBATCH -p development # queue (partition) -- normal, development, etc.
#SBATCH -t 00:30:00         # run time (hh:mm:ss)

ibrun tacc_affinity ./isotropy_mpi config.cfg 
