#!/bin/bash

# Example of a SLURM script which can be 
# used to run the simulations

# The file path separators may need to be adjusted if you're not using a 
# Unix-like operating system. (For example, replace with a \)

# You will need to specify more SBATCH directives
# to run on your HPC instance

#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=512x64-ib
#SBATCH --qos=pi-shyun
#SBATCH --account=pi-shyun
#SBATCH --time=07-00:00:00
#SBATCH --output=./logs/%x_%A_%a.log
#SBATCH --error=./logs/%x_%A_%a.err
#SBATCH --array=1-144
#SBATCH --job-name=sim-all
#SBATCH --mail-user=epawl@ucsc.edu
#SBATCH --mail-type=ALL

###

# Load necessary modules
module load miniconda3/3.12
conda activate R-4.4.3
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
Rscript 1_simulation/sim_04.R $SLURM_ARRAY_TASK_ID TRUE && echo "success" > ${success_dir}/${SLURM_ARRAY_TASK_ID}.txt || echo "failure" > ${failure_dir}/${SLURM_ARRAY_TASK_ID}.txt
