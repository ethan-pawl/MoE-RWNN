#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4GB
#SBATCH --partition=512x64-ib               # IB partition
#SBATCH --qos=pi-shyun                        # Both of these lines need to be present in order to use the 
#SBATCH --account=pi-shyun                  #   IB network and parititons
#SBATCH --time=2:00:00
#SBATCH --mail-user=epawl@ucsc.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=./logs/%x_%A_%a.log
#SBATCH --error=./logs/%x_%A_%a.err

# Set the session up to use R
module load miniconda3/3.12
conda activate main-env

mkdir -p ${cv_step}/successes
mkdir -p ${cv_step}/failures

Rscript pcr_oos.R $cv_step $SLURM_ARRAY_TASK_ID && echo "success" > ${cv_step}/successes/${SLURM_ARRAY_TASK_ID}.txt || echo "failure" > ${cv_step}/failures/${SLURM_ARRAY_TASK_ID}.txt