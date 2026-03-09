#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G
#SBATCH --time=07-00:00:00
#SBATCH --output=./logs/%x_%A_%a.log
#SBATCH --error=./logs/%x_%A_%a.err
#SBATCH --array=1-12
#SBATCH --job-name=n_PCs

# Load modules here

Rscript 2_application/23_PCs/pcs.R $SLURM_ARRAY_TASK_ID