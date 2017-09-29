#!/bin/bash
#SBATCH --job-name=models_bin_arr
#SBATCH --output=models_bin_%A_%a.out
#SBATCH --error=models_bin_%A_%a.err
#SBATCH --array=1-3
#SBATCH --time=48:00:00 
#SBATCH --partition=normal
#SBATCH --mem=50000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=erflynn@stanford.edu

echo $SLURM_ARRAY_TASK_ID

Rscript bin_model.R $SLURM_ARRAY_TASK_ID

