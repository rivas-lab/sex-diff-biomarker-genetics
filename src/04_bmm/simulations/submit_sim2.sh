#!/bin/bash
#SBATCH --job-name=simM1
#SBATCH --output=logs/sim_%A_%a.out
#SBATCH --error=logs/sim_%A_%a.err
#SBATCH --time=3:00:00 
#SBATCH --partition=rbaltman,owners
#SBATCH --mem=8000


ml r-rstan
Rscript simM2.R $SLURM_ARRAY_TASK_ID
