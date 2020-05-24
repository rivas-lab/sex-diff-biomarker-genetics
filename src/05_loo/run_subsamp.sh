#!/bin/bash
#SBATCH --job-name=ss%A
#SBATCH --output=logs/ss_%A_%a.out
#SBATCH --error=logs/ss_%A_%a.err
#SBATCH --nodes=4
#SBATCH --time=2:00:00 
#SBATCH --partition=rbaltman,owners
#SBATCH --mem=20000
# usage:
#    sbatch submit_run_bmm.sh <model> <trait> <ndim>  <outdir> <indir>



#ml r-rstan
ml purge
ml load R/3.5.1
ml load gcc/8.1.0

Rscript src/05_loo/run_subsampled_loo.R
