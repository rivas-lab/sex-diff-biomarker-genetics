#!/bin/bash
#SBATCH --job-name=bmm_meno_%A
#SBATCH --output=logs/bmm_meno_%A_%a.out
#SBATCH --error=logs/bmm_meno_%A_%a.err
#SBATCH --array=1
#SBATCH --time=48:00:00 
#SBATCH --partition=rbaltman
#SBATCH --mem=50000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=erflynn@stanford.edu

echo $SLURM_ARRAY_TASK_ID


#Rscript run_bmm.R $SLURM_ARRAY_TASK_ID '21001' 'quant' 'TRUE' 'opt'
#Rscript run_bmm.R $SLURM_ARRAY_TASK_ID 'RH107' 'binary' 'TRUE' 'opt'
#Rscript run_bmm.R $SLURM_ARRAY_TASK_ID 'RH107' 'binary' 'TRUE' 'vb'
Rscript run_bmm_m1.R $SLURM_ARRAY_TASK_ID $1 'quant' $2 $3
