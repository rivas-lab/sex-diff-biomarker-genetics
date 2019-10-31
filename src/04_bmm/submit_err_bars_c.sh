#!/bin/bash
#SBATCH --job-name=err_%A
#SBATCH --output=logs/err_%A_%a.out
#SBATCH --error=logs/err_%A_%a.err
#SBATCH --time=40:00 
#SBATCH --partition=rbaltman,owners
#SBATCH --mem=8000

trait=$1
dir=$2

ml r-rstan
Rscript err_bars_complete.R $trait $dir $SLURM_ARRAY_TASK_ID
