#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=63
#SBATCH --mem-per-cpu=3GB
#SBATCH --time=07-00:00:00
#SBATCH --output=./logs/%x_%A.log
#SBATCH --error=./logs/%x_%A.err
#SBATCH --job-name=3dsim

# Load modules here

Rscript 4_3dsim/3dsim.R
