#!/bin/bash

script_path="2_application/21_all_data/pcr.slurm"

jid1=$(sbatch --parsable --export=cv_step='maxres' --job-name=app_maxres --array=1-10 "$script_path")
jid2=$(sbatch --parsable --dependency=afterok:$jid1 --export=cv_step='cv' --job-name=app_cv --array=1-5000 "$script_path")
jid3=$(sbatch --parsable --dependency=afterok:$jid2 --export=cv_step='refit' --job-name=app_refit --array=1-1000 "$script_path")
sbatch --dependency=afterok:$jid3 --export=cv_step='summary' --job-name=app_summary --array=1-10 "$script_path"
