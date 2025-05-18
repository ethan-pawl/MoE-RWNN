#!/bin/bash

# Example of a SLURM script which can be 
# used to run the simulations

# You will need to specify more SBATCH directives
# to run on your HPC instance

#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=500MB
#SBATCH --time=02:00:00
#SBATCH --output=./logs/%x_%A_%a.log
#SBATCH --error=./logs/%x_%A_%a.err
#SBATCH --array=1-140
#SBATCH --job-name=simulation

# Load necessary modules

mkdir successes
mkdir failures

# Loop through all simulation scenarios while keeping track of 
# which jobs fail and which run successfully
Rscript sim_04.R $SLURM_ARRAY_TASK_ID && echo "success" > successes/${SLURM_ARRAY_TASK_ID}.txt || echo "failure" > failures/${SLURM_ARRAY_TASK_ID}.txt