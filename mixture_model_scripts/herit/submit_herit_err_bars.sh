#!/bin/bash
#SBATCH --job-name=h_err
#SBATCH --output=logs/h_err_%A_%a.out
#SBATCH --error=logs/h_err_%A_%a.err
#SBATCH --time=1:30:00
#SBATCH --array=4
#SBATCH -p rbaltman
#SBATCH --mem=8000

ml r-rstan
Rscript heritErrBars.R ${SLURM_ARRAY_TASK_ID}

