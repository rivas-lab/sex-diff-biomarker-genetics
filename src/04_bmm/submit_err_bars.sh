#!/bin/bash
#SBATCH --job-name=err_%A
#SBATCH --output=logs/err_%A_%a.out
#SBATCH --error=logs/err_%A_%a.err
#SBATCH --time=2:00:00 
#SBATCH --mem=8000

trait=$1
dir=$2

ml purge
ml load R/3.5.1
ml load gcc/8.1.0
ml load python/3.6

Rscript src/04_bmm/calc_err_bars.R $trait $dir $SLURM_ARRAY_TASK_ID
