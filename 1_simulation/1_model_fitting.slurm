#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --partition=windfall
#SBATCH --qos=windfall
#SBATCH --account=pi-shyun
#SBATCH --time=07-00:00:00
#SBATCH --output=./logs/%x_%A_%a.log
#SBATCH --error=./logs/%x_%A_%a.err
#SBATCH --array=1-140
#SBATCH --job-name=model_fitting
#SBATCH --mail-user=epawl@ucsc.edu
#SBATCH --mail-type=ALL

R_BIN="/mnt/beegfs/home/epawl/.conda/envs/R-4.4.3-2026-2/bin/Rscript"

"$R_BIN" 1_simulation/model_fitting.R $SLURM_ARRAY_TASK_ID