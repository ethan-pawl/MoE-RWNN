#!/bin/bash

script_path="3_seedsim/seedsim.slurm"

jid1=$(sbatch --parsable --export=cv_step='maxres' --job-name=seedsim_maxres --array=1-5 "$script_path")
jid2=$(sbatch --parsable --dependency=afterok:$jid1 --export=cv_step='cv' --job-name=seedsim_cv --array=1-5000%225 "$script_path")
jid3=$(sbatch --parsable --dependency=afterok:$jid2 --export=cv_step='refit' --job-name=seedsim_refit --array=1-500%225 "$script_path")
sbatch --dependency=afterok:$jid3 --export=cv_step='summary' --job-name=seedsim_summary --array=1-5 "$script_path"