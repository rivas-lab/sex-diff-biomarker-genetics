#!/bin/bash
#SBATCH --job-name=h_snp
#SBATCH --output=logs/h_snp_%A_%a.out
#SBATCH --error=logs/h_snp_%A_%a.err
#SBATCH --time=2:00:00
#SBATCH --array=1-5
#SBATCH -p rbaltman
#SBATCH --mem=10000

ml r-rstan
Rscript extract_snp_tabs.R $SLURM_ARRAY_TASK_ID

