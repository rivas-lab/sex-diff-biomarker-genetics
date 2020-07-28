#!/bin/bash
#SBATCH --job-name=ll%A
#SBATCH --output=logs/ll_%A_%a.out
#SBATCH --error=logs/ll_%A_%a.err
#SBATCH --nodes=4
#SBATCH --time=2:00:00 
#SBATCH --mem=20000

ml purge
ml load R/3.5.1
ml load gcc/8.1.0

Rscript src/05_loo/extract_loo.R $SLURM_ARRAY_TASK_ID
