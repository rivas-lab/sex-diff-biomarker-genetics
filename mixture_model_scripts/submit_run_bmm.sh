#!/bin/bash
#SBATCH --job-name=bmm_bio_%A
#SBATCH --output=logs/bmm_bio_%A_%a.out
#SBATCH --error=logs/bmm_bio_%A_%a.err
#SBATCH --array=1
#SBATCH --time=48:00:00 
#SBATCH --partition=rbaltman
#SBATCH --mem=10000
# usage:
#    sbatch submit_run_bmm.sh <trait> <ndim> <downsampled>

echo $SLURM_ARRAY_TASK_ID

ml r-rstan

#Rscript run_bmm.R $SLURM_ARRAY_TASK_ID '21001' 'quant' 'TRUE' 'opt'
#Rscript run_bmm.R $SLURM_ARRAY_TASK_ID 'RH107' 'binary' 'TRUE' 'opt'
#Rscript run_bmm.R $SLURM_ARRAY_TASK_ID 'RH107' 'binary' 'TRUE' 'vb'
Rscript run_bmm_m1.R $SLURM_ARRAY_TASK_ID $1 'quant' $2 $3 $4
