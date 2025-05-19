#!/bin/bash

# Example of a bash script which calls the SLURM script 
# to sequence the cross-validation steps

# The max array indices would need to be adjusted if 
# the cross-validation settings are altered

# The file path separators may need to be adjusted if you're not using a 
# Unix-like operating system. (For example, replace with a \)

script_path="2_application/20_nested_cv/pcr_oos.slurm"

jid1=$(sbatch --parsable --export=cv_step='maxres' --job-name=pcr_oos_maxres --array=1-10 "$script_path")
jid2=$(sbatch --parsable --dependency=afterok:$jid1 --export=cv_step='cv' --job-name=pcr_oos_cv --array=1-5000 "$script_path")
jid3=$(sbatch --parsable --dependency=afterok:$jid2 --export=cv_step='refit' --job-name=pcr_oos_refit --array=1-1000 "$script_path")
sbatch --dependency=afterok:$jid3 --export=cv_step='summary' --job-name=pcr_oos_summary --array=1-10 "$script_path"
