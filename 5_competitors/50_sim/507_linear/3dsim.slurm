#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --mem-per-cpu=4GB
#SBATCH --time=07-00:00:00
#SBATCH --output=./logs/%x_%A.log
#SBATCH --error=./logs/%x_%A.err
#SBATCH --partition=windfall
#SBATCH --qos=windfall
#SBATCH --account=pi-shyun
#SBATCH --mail-user=epawl@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=3dsim_linear

# Hard-coding the path to Rscript is more robust on Hummingbird 
# than loading the conda environment and calling Rscript

# Activate your conda environment and 
# put the output of `whereis Rscript` below
R_BIN='/mnt/beegfs/home/epawl/.conda/envs/R-4.4.3-2026/bin/Rscript'

"${R_BIN}" 5_competitors/50_sim/507_linear/3dsim.R
