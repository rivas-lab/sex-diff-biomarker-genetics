#!/bin/bash
#SBATCH --job-name=run_bmm_%A
#SBATCH --output=logs/run_bmm_%A_%a.out
#SBATCH --error=logs/run_bmm_%A_%a.err
#SBATCH --nodes=4
#SBATCH --time=24:00:00 
#SBATCH --partition=rbaltman,owners
#SBATCH --mem=30000
# usage:
#    sbatch submit_run_bmm.sh <model> <trait> <ndim>  <outdir> <indir>



ml r-rstan

model=$SLURM_ARRAY_TASK_ID
trait=$1
ndim=$2
outdir=$3
indir=$4
suffix=$5
Rscript run_bmm.R ${model} ${trait} ${ndim} ${outdir} ${indir} ${suffix}
