#!/bin/bash
#SBATCH -c 1                               # Core Number 
#SBATCH -t 0-12:00                         # Runtime in D-HH:MM format
#SBATCH -p defq                           # Partition to run in
#SBATCH --mem=8G                           # Memory

ml load Anaconda3/5.0.1
ml load GCC/8.2.0-2.31.1
source activate pysb

python run_MARM1_simulation_reduced.py "$@"
