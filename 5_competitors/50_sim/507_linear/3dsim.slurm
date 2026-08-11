#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --mem-per-cpu=4GB
#SBATCH --time=07-00:00:00
#SBATCH --output=./logs/%x_%A.log
#SBATCH --error=./logs/%x_%A.err
#SBATCH --job-name=3dsim_linear

# Load necessary modules

Rscript 5_competitors/50_sim/507_linear/3dsim.R
