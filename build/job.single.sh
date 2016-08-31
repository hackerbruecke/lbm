#!/bin/bash
# Job file for Stampede supercomputer
#SBATCH -J LBM_SINGLE # job name
#SBATCH -o LBM_SINGLE_OUT_%j # output file name (%j expands to jobID)
#SBATCH -e LBM_SINGLE_ERR_%j # error file name (%j expands to jobID)
#SBATCH -N 1                # total number of nodes requested (16 cores/node)
#SBATCH -n 1                # 1 task
#SBATCH -p development # queue (partition) -- normal, development, etc.
#SBATCH -t 00:30:00         # run time (hh:mm:ss)

module load boost
./lbm config.cfg
