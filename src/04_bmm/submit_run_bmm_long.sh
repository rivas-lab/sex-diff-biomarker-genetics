#!/bin/bash
#SBATCH --job-name=lr_bmm_%A
#SBATCH --output=logs/lr_bmm_%A_%a.out
#SBATCH --error=logs/lr_bmm_%A_%a.err
#SBATCH --nodes=4
#SBATCH --time=48:00:00 
#SBATCH --partition=rbaltman,owners
#SBATCH --mem=30000
# usage:
#    sbatch submit_run_bmm.sh <model> <trait> <ndim>  <outdir> <indir>



ml purge
ml load R/3.5.1
ml load gcc/8.1.0

model=2
trait=$1
ndim=2
outdir=data2/
indir=data/gwas_0522/
#suffix=$5
Rscript src/04_bmm/run_bmm_long.R ${model} ${trait} ${ndim} ${outdir} ${indir} #${suffix}
