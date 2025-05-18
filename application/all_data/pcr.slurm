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
#SBATCH --mail-type=ALL
#SBATCH --output=./logs/%x_%a.log
#SBATCH --error=./logs/%x_%a.err

# Set the session up to use R
module load miniconda3/3.12
conda activate main-env

mkdir -p ${cv_step}/successes
mkdir -p ${cv_step}/failures

Rscript pcr.R $cv_step $SLURM_ARRAY_TASK_ID && echo "Job $cv_step $SLURM_ARRAY_TASK_ID was a success" > ${cv_step}/successes/${SLURM_ARRAY_TASK_ID}.txt || echo "Job $cv_step $SLURM_ARRAY_TASK_ID was a failure" > ${cv_step}/failures/${SLURM_ARRAY_TASK_ID}.txt
