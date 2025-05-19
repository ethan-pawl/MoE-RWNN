#!/bin/bash

# Example of a SLURM script which can be 
# used to run the application

# The file path separators may need to be adjusted if you're not using a 
# Unix-like operating system. (For example, replace with a \)

# You will need to specify more SBATCH directives
# to run on your HPC instance

#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4GB
#SBATCH --time=2:00:00
#SBATCH --output=./logs/%x_%A_%a.log
#SBATCH --error=./logs/%x_%A_%a.err

###

# Load any necessary modules here 
# (may not be necessary depending on your HPC instance)

###

# Create some directories to keep track of which jobs fail
# and which run successfully
# Feel free to adjust file paths as per your convenience
success_dir="2_application/20_nested_cv/job_tracking/${cv_step}/successes"
failure_dir="2_application/20_nested_cv/job_tracking/${cv_step}/failures"

mkdir -p "${success_dir}"
mkdir -p "${failure_dir}"

Rscript 2_application/20_nested_cv/pcr_oos.R $cv_step $SLURM_ARRAY_TASK_ID && echo "success" > ${success_dir}/${SLURM_ARRAY_TASK_ID}.txt || echo "failure" > ${failure_dir}/${SLURM_ARRAY_TASK_ID}.txt