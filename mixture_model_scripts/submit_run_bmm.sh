#!/bin/bash
#SBATCH --job-name=bmm_arr
#SBATCH --output=bmm_%A_%a.out
#SBATCH --error=bmm_%A_%a.err
#SBATCH --array=1-3
#SBATCH --time=48:00:00 
#SBATCH --partition=normal
#SBATCH --mem=50000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=erflynn@stanford.edu

echo $SLURM_ARRAY_TASK_ID


Rscript run_bmm.R $SLURM_ARRAY_TASK_ID '21001' 'quant' 'TRUE' 
#Rscript run_bmm.R $SLURM_ARRAY_TASK_ID 'RH107' 'binary'  