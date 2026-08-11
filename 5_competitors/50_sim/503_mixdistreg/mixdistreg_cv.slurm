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
#SBATCH --job-name=mixdistreg_cv_redo
#SBATCH --array=1-500

# Load modules

Rscript 5_competitors/50_sim/503_mixdistreg/mixdistreg_cv.R ${SLURM_ARRAY_TASK_ID}
