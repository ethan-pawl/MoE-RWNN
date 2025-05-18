#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=500MB
#SBATCH --partition=512x64-ib
#SBATCH --qos=pi-shyun
#SBATCH --account=pi-shyun
#SBATCH --time=07-00:00:00
#SBATCH --mail-user=epawl@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH --output=./logs/%x_%A_%a.log
#SBATCH --error=./logs/%x_%A_%a.err
#SBATCH --array=1-140
#SBATCH --job-name=sim_04

module load miniconda3/3.12
conda activate main-env

mkdir successes
mkdir failures

Rscript sim_04.R $SLURM_ARRAY_TASK_ID && echo "success" > successes/${SLURM_ARRAY_TASK_ID}.txt || echo "failure" > failures/${SLURM_ARRAY_TASK_ID}.txt