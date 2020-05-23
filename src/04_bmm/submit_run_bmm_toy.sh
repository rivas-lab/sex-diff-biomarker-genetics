#!/bin/bash
#SBATCH --job-name=run_bmm_%A
#SBATCH --output=logs/run_bmm_%A_%a.out
#SBATCH --error=logs/run_bmm_%A_%a.err
#SBATCH --nodes=4
#SBATCH --time=4:00:00 
#SBATCH --partition=rbaltman,owners
#SBATCH --mem=10000
# usage:
#    sbatch submit_run_bmm.sh <model> <trait> <ndim>  <outdir> <indir>



#ml r-rstan
ml purge
ml load R/3.5.1
ml load gcc/8.1.0

model=$SLURM_ARRAY_TASK_ID
trait=Test_chr22
ndim=2
outdir=data/toy_run/
indir=data/gwas/
Rscript src/04_bmm/run_bmm_vary_params.R ${model} ${trait} ${ndim} ${outdir} ${indir} 
