#!/bin/bash

# You will need to specify more SBATCH directives
# to run on your HPC instance

#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4GB
#SBATCH --time=07-00:00:00
#SBATCH --output=./logs/%x_%A_%a.log
#SBATCH --error=./logs/%x_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --job-name=gam_cv
#SBATCH --array=1-50

# Load modules here

# Loop through all cross-validation jobs and save success/fail state
Rscript 5_competitors/50_sim/501_gam/gam_sim_cv_hpc.R ${SLURM_ARRAY_TASK_ID}
