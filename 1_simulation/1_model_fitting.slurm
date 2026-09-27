#!/bin/bash

# Example of a SLURM script which can be 
# used to run the application

# The file path separators may need to be adjusted if you're not using a 
# Unix-like operating system. (For example, replace with a \)

# You will need to specify more SBATCH directives
# to run on your HPC instance

#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=07-00:00:00
#SBATCH --output=./logs/%x_%A_%a.log
#SBATCH --error=./logs/%x_%A_%a.err
#SBATCH --array=1-140
#SBATCH --job-name=model_fitting

Rscript 1_simulation/model_fitting.R $SLURM_ARRAY_TASK_ID