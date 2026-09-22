#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=07-00:00:00
#SBATCH --output=./logs/%x_%A_%a.log
#SBATCH --error=./logs/%x_%A_%a.err
#SBATCH --array=1-140
#SBATCH --job-name=model_fitting
#SBATCH --mail-type=ALL

Rscript 1_simulation/model_fitting.R $SLURM_ARRAY_TASK_ID