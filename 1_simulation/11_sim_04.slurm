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

###

# Load necessary modules

###

# Create some directories to keep track of which jobs fail
# and which run successfully
# Feel free to adjust file paths as per your convenience
success_dir="1_simulation/job_tracking/${cv_step}/successes"
failure_dir="1_simulation/job_tracking/${cv_step}/failures"

mkdir -p "${success_dir}"
mkdir -p "${failure_dir}"

# Loop through all simulation scenarios while keeping track of 
# which jobs fail and which run successfully
Rscript sim_04.R $SLURM_ARRAY_TASK_ID && echo "success" > ${success_dir}/${SLURM_ARRAY_TASK_ID}.txt || echo "failure" > ${failure_dir}/${SLURM_ARRAY_TASK_ID}.txt