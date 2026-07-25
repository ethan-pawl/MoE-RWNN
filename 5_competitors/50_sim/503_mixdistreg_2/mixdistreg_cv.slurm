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
#SBATCH --partition=windfall
#SBATCH --qos=windfall
#SBATCH --account=pi-shyun
#SBATCH --mail-user=epawl@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=mixdistreg_cv_redo
#SBATCH --array=1-500

# Hard-coding the path to Rscript is more robust on Hummingbird 
# than loading the conda environment and calling Rscript

# Activate your conda environment and 
# put the output of `whereis Rscript` below
R_BIN='/mnt/beegfs/home/epawl/.conda/envs/R-4.4.3-2026-2/bin/Rscript'

# Loop through all cross-validation jobs and save success/fail state
"${R_BIN}" 5_competitors/50_sim/503_mixdistreg/mixdistreg_cv.R ${SLURM_ARRAY_TASK_ID}