#!/bin/bash
jid1=$(sbatch --parsable --export=cv_step='maxres' --job-name=pcr_maxres --array=1-2 pcr.slurm)
jid2=$(sbatch --parsable --dependency=afterok:$jid1 --export=cv_step='cv' --job-name=pcr_cv --array=1-5000 pcr.slurm)
jid3=$(sbatch --parsable --dependency=afterok:$jid2 --export=cv_step='refit' --job-name=pcr_refit --array=1-200 pcr.slurm)
sbatch --dependency=afterok:$jid3 --export=cv_step='summary' --job-name=pcr_summary--array=1-2 pcr.slurm
