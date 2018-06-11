#!/bin/bash
#SBATCH --job-name=update_m2
#SBATCH --output=logs/update_m2_%a_%A.out
#SBATCH --error=logs/update_m2_%a_%A.err
#SBATCH --time=1:00:00
#SBATCH --partition=rbaltman
#SBATCH --array=8,9
#SBATCH --mem=8000

Rscript update_snp_estimates.R $SLURM_ARRAY_TASK_ID
