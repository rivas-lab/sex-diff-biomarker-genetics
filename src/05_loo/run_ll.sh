#!/bin/bash
#SBATCH --job-name=ll%A
#SBATCH --output=logs/ll_%A_%a.out
#SBATCH --error=logs/ll_%A_%a.err
#SBATCH --nodes=4
#SBATCH --time=2:00:00 
#SBATCH --mem=10000
# usage:
#    sbatch submit_run_bmm.sh <model> <trait> <ndim>  <outdir> <indir>



#ml r-rstan
ml purge
ml load R/3.5.1
ml load gcc/8.1.0

model_id=$1
params_id=$2

Rscript src/05_loo/compute_LL.R ${model_id} ${params_id} $SLURM_ARRAY_TASK_ID
