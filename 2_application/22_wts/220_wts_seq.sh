#!/bin/bash

script_path="2_application/22_wts/wts.slurm"

jid1=$(sbatch --parsable --export=cv_step='maxres' --job-name=no_wts_maxres --array=1 "$script_path")
jid2=$(sbatch --parsable --dependency=afterok:$jid1 --export=cv_step='cv' --job-name=no_wts_cv --array=1-5000%300 "$script_path")
jid3=$(sbatch --parsable --dependency=afterok:$jid2 --export=cv_step='refit' --job-name=no_wts_refit --array=1-100 "$script_path")
sbatch --dependency=afterok:$jid3 --export=cv_step='summary' --job-name=no_wts_summary --array=1 "$script_path"