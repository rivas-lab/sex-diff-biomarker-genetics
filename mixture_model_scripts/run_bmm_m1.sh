#!/bin/bash
#SBATCH --job-name=bmm_arr
#SBATCH --output=logs/bmm_%A_%a.out
#SBATCH --error=logs/bmm_%A_%a.err
#SBATCH --array=1
#SBATCH --time=48:00:00 
#SBATCH --partition=rbaltman
#SBATCH --mem=50000

echo $SLURM_ARRAY_TASK_ID

trait=$1
trait_type=$2
Rscript run_bmm_m1.R $SLURM_ARRAY_TASK_ID ${trait} ${trait_type}
#Rscript run_bmm.R $SLURM_ARRAY_TASK_ID ${trait} 'quant' #'FALSE' 'opt' $SLURM_ARRAY_TASK_ID
#Rscript run_bmm.R $SLURM_ARRAY_TASK_ID 'RH107' 'binary' 'TRUE' 'opt'
#Rscript run_bmm.R $SLURM_ARRAY_TASK_ID 'RH107' 'binary' 'TRUE' 'vb'
